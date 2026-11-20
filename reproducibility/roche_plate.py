"""Replay the Roche Plate measurements from the two extracted Zenodo archives."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import capacities, solve_lp, frontier, zero_frontier

PV_KW = 0.75
WIND_KW = 0.25
PV_DERATING = 0.80
WIND_CUT_IN, WIND_RATED, WIND_CUT_OUT = 2.0, 8.0, 25.0
CHARGE_KW, DISCHARGE_KW, LOW_DISCHARGE_KW = 0.5, 1.2, 0.05
EFFICIENCY = 0.95


def read_measurements(path, columns):
    data = pd.read_csv(path, skiprows=1, usecols=["Date", *columns], parse_dates=["Date"])
    grouped = data.dropna(subset=["Date"]).set_index("Date").sort_index().groupby(level=0)
    conflict = grouped[columns].nunique().gt(1).any(axis=1)
    return grouped.first().loc[~conflict]


def weather_inputs(source):
    weather = read_measurements(source / "Meteo_dataset.csv", ["Gincl", "Ws"])
    weather = weather.where(weather.ge(0))
    irradiance = weather.Gincl
    speed = weather.Ws
    weather["PV_kW"] = PV_KW * (PV_DERATING * irradiance).clip(upper=1)
    wind_cf = ((speed**3 - WIND_CUT_IN**3) / (WIND_RATED**3 - WIND_CUT_IN**3)).clip(0, 1)
    weather["Wind_kW"] = WIND_KW * wind_cf.where(speed.lt(WIND_CUT_OUT) | speed.isna(), 0)
    return weather


def get_chronology(source):
    weather = weather_inputs(source)
    demand = read_measurements(source / "Demand_house_1_dataset.csv", ["Input4"])
    demand = demand.where(demand.ge(0)).rename(columns={"Input4": "Load_kW"}) / 1000
    wh_hours = weather.resample("h")
    load_hours = demand.resample("h")
    means = wh_hours.mean().join(load_hours.mean(), how="inner")
    counts = wh_hours.count().join(load_hours.count(), how="inner")
    valid = (counts.Gincl >= 171) & (counts.Ws >= 171) & (counts.Load_kW >= 57)
    valid = valid.loc["2020-08-01":]
    valid_groups = (~valid).cumsum()
    runs = valid_groups[valid].value_counts(sort=False)
    best_run = runs.idxmax()
    index = valid.index[valid & valid_groups.eq(best_run)]
    hourly = means.loc[index].copy()
    hourly["Gincl_count"] = counts.loc[index, "Gincl"]
    hourly["Ws_count"] = counts.loc[index, "Ws"]
    hourly["Load_count"] = counts.loc[index, "Load_kW"]
    return hourly.rename_axis("Recorded_time")


def solve_case(generation, load, discharge, allowance):
    cut_size, history = frontier(
        generation, load, CHARGE_KW, discharge, allowance, EFFICIENCY, EFFICIENCY
    )
    lp_size, result = solve_lp(
        generation, load, CHARGE_KW, discharge, allowance, EFFICIENCY, EFFICIENCY
    )
    charge, required, delivered = capacities(
        generation, load, CHARGE_KW, discharge, EFFICIENCY, EFFICIENCY
    )
    target = max(0, required.sum() - allowance * load.sum() / EFFICIENCY)
    if np.isinf(cut_size):
        assert np.isinf(lp_size) and delivered.sum() < target <= charge.sum()
        return {"allowance": allowance, "discharge_kW": discharge, "cut_kWh": None,
                "lp_kWh": None, "storage_delivery_target_kWh": float(target),
                "discharge_energy_bound_kWh": float(delivered.sum()),
                "unserved_floor_kWh": float(EFFICIENCY * (required.sum() - delivered.sum()))}, None

    assert abs(cut_size - lp_size) < 2e-7
    if allowance == 0:
        interval_result = zero_frontier(
            generation, load, CHARGE_KW, discharge, EFFICIENCY, EFFICIENCY
        )
        assert abs(cut_size - interval_result) < 2e-7
    low_size, _ = solve_lp(
        generation, load, CHARGE_KW, discharge, allowance, EFFICIENCY, EFFICIENCY,
        Efixed=0.99*cut_size
    )
    assert np.isinf(low_size)
    n = len(load)
    inventory = result.x[:n]
    shortfall = result.x[3*n:4*n]
    assert inventory.min() >= -1e-8 and inventory.max() <= cut_size + 2e-7
    assert shortfall.sum() <= allowance * load.sum() + 2e-7
    binding = history[-2]
    assert abs(binding["C"] + binding["k"]*cut_size - target) < 2e-7
    return {"allowance": allowance, "discharge_kW": discharge,
            "cut_kWh": cut_size, "lp_kWh": lp_size,
            "ens_kWh": float(shortfall.sum()), "limiting_cut_blocks": binding["k"],
            "charge_cap_kWh": float(charge.sum())}, inventory


def find_critical_interval(generation, load):
    n = len(load)
    a, b, _ = capacities(generation, load, CHARGE_KW, DISCHARGE_KW, EFFICIENCY, EFFICIENCY)
    cumulative = np.r_[0.0, np.cumsum(np.tile(b - a, 2))]
    draws = np.array([cumulative[i+1:i+n] - cumulative[i] for i in range(n)])
    start, length_index = np.unravel_index(draws.argmax(), draws.shape)
    length = int(length_index + 1)
    steps = (start + np.arange(length)) % n
    return {"start_index": int(start), "length_h": length,
            "draw_kWh": float(draws[start, length_index]), "recharge_kWh": float(a[steps].sum())}


def plot_inventory(hourly, critical, target):
    x = np.arange(len(hourly) + 1)
    figure, axis = plt.subplots(figsize=(7.2, 2.5), layout="constrained")
    axis.step(x, np.r_[hourly.Zero_shortfall_SOC_kWh, hourly.Zero_shortfall_SOC_kWh.iloc[0]],
              where="post", linewidth=1.15, color="#12667c", label="Zero shortfall")
    axis.step(x, np.r_[hourly.One_percent_SOC_kWh, hourly.One_percent_SOC_kWh.iloc[0]],
              where="post", linewidth=0.9, alpha=0.8, color="#a04e31", label="1% allowance")
    start, length = critical["start_index"], critical["length_h"]
    axis.axvspan(start, min(start+length, len(hourly)), color="#8499a4", alpha=0.20,
                 label="Zero-shortfall critical interval")
    if start + length > len(hourly):
        axis.axvspan(0, start+length-len(hourly), color="#8499a4", alpha=0.20)
    axis.set(xlim=(0, len(hourly)), xlabel="Elapsed hour", ylabel="Stored energy (kWh)")
    axis.grid(alpha=0.17)
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, 1.19), ncol=3, fontsize=8, frameon=False)
    figure.savefig(target)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("figure", type=Path)
    args = parser.parse_args()
    hourly = get_chronology(args.raw_dir)
    generation = (hourly.PV_kW + hourly.Wind_kW).to_numpy()
    load = hourly.Load_kW.to_numpy()
    zero, hourly["Zero_shortfall_SOC_kWh"] = solve_case(generation, load, DISCHARGE_KW, 0)
    positive, hourly["One_percent_SOC_kWh"] = solve_case(generation, load, DISCHARGE_KW, 0.01)
    limited, _ = solve_case(generation, load, LOW_DISCHARGE_KW, 0.01)
    critical = find_critical_interval(generation, load)
    assert abs(critical["draw_kWh"]-zero["cut_kWh"]) < 2e-7
    critical["recorded_start"] = str(hourly.index[critical["start_index"]])
    summary = dict(start=str(hourly.index[0]), last=str(hourly.index[-1]), hours=len(hourly),
                   load_kWh=float(load.sum()), peak_load_kW=float(load.max()), pv_kWh=float(hourly.PV_kW.sum()),
                   wind_kWh=float(hourly.Wind_kW.sum()), mean_wind_ms=float(hourly.Ws.mean()),
                   critical=critical, cases=[zero, positive, limited])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    hourly.to_csv(args.output_dir / "roche_plate_hourly.csv", float_format="%.12g")
    (args.output_dir / "roche_plate_results.json").write_text(json.dumps(summary, indent=2)+"\n")
    plot_inventory(hourly, critical, args.figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
