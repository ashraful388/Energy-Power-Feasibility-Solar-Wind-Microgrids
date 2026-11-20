import json
import sys
from pathlib import Path
import numpy as np

source = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(source))
from annual_study import generate, interval_frontier
from analysis import capacities

rows = []
for rho in (0.95, 0.99):
    pv, wind, load, _ = generate(20261017, rho=rho)
    generation = 300 * (pv + wind)
    result = interval_frontier(generation, load, 200, 160)
    a,b,d = capacities(generation, load, 200, 160, 0.95, 0.95)
    i = (np.arange(result["length"]) + result["start"]) % len(load)
    assert abs((b[i]-a[i]).sum()-result["E"]) < 1e-7
    rows.append(dict(rho=rho, **result, annual_generation=generation.sum(), annual_pv=(pv*300).sum(), annual_wind=(wind*300).sum(), annual_load=load.sum(), interval_residual=(b[i]*0.95).sum(), interval_withdrawal=b[i].sum(), interval_recharge=a[i].sum(), interval_charge_hours=int((a[i]>0).sum()), interval_deficit_hours=int((b[i]>0).sum())))
print(json.dumps(rows, indent=2))