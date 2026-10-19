# Energy–Power Feasibility Boundaries for Battery Storage in Stand-Alone Solar–Wind Microgrids

Simulation code and synthetic data supporting the manuscript of this title.
No measured site data are included. These profiles represent constructed scenarios, not observations or forecasts for a specific location.

**Authors:** Rajae Gaamouche $^{1,2,*}$, Md. Ashraful Babu $^{3,4}$, and Mohammed Lahby $^{5}$  
$^{1}$ Department of Electronics and Electromechanics, National School of Mineral Industry (ENIM), Rabat 10000, Morocco  
$^{2}$ SmartiLab Laboratory, Moroccan School of Sciences of Engineering (EMSI), Rabat 10000, Morocco  
$^{3}$ Department of Physical Sciences, Independent University, Bangladesh, Dhaka 1245, Bangladesh  
$^{4}$ Explainable Artificial Intelligence and Optimization Research (XAiOR), Independent University, Bangladesh, Dhaka 1245, Bangladesh  
$^{5}$ Laboratory of Mathematics, Artificial Intelligence, and Digital Learning, Higher Normal School, Hassan II University, Casablanca 50069, Morocco  

*Correspondence: `r.gaamouche@enim.ac.ma`

---

## Contents and provenance

- `reproducibility/`: original analysis scripts, dependency specifications and archived numerical outputs supplied with the study. The scripts and archived CSV/JSON files are preserved unchanged.
- `reproducibility/results/synthetic_profiles.csv`: original 336-hour synthetic input trace.
- `reproducibility/annual_results/annual_base_profiles.csv`: original 8,760-hour base profiles.
- `reproducibility/annual_results/*dispatch*.csv`: original saved dispatch trajectories, including the six-step example and shortfall allowances.
- `reproducibility/annual_results/*.json` and `reproducibility/results/results.json`: archived experiment results and summaries. Python JSON permits `Infinity` for infeasible cases; strict JSON readers must handle this explicitly.
- `datasets/annual_seed_profiles.zip`: regenerated CSV inputs for all 30 original annual seeds (20261017–20261046).
- `datasets/scenario_profiles.zip`: regenerated generation/demand inputs for the 14 original scenarios.
- `datasets/quarter_hour_profiles.csv.gz`: regenerated 35,040-step, 15-minute input profiles.
- `figures/`: original manuscript figure assets.
- `export_datasets.py`: additional export and verification utility; it calls the original generator without modifying it.
- `verification.json`: checks performed when preparing this repository. The original base CSV is compared against the regenerated base; all 330 annual mixture frontiers, 14 scenarios and the quarter-hour frontier are compared with archived results. A base-case LP is solved independently of the interval calculation. Archived timings are not new benchmarks.
- `SHA256SUMS.json`: checksums for the release files (excluding this checksum file and Git metadata).

The exported inputs were regenerated for this release from the original code and seeds; they are not presented as previously archived files. Results for other power ratings, efficiencies, reliability allowances and mixtures use the same inputs with the transformations in `annual_study.py`. Small factorial/random tests and timing traces are generated within `analysis.py` and `extended_checks.py`.

## Reproduction

Python 3.13 was used for release verification. Install dependencies in an isolated environment:

```sh
python -m venv .venv
# Activate the environment using the command appropriate for your operating system.
python -m pip install -r reproducibility/requirements.txt
python export_datasets.py
```

To rerun the complete original study (these commands overwrite archived numerical outputs locally):

```sh
python reproducibility/analysis.py
python reproducibility/annual_study.py
python reproducibility/extended_checks.py
python reproducibility/make_figures.py
python reproducibility/report_annual.py
```

The full study includes many optimisation problems and can take substantial time. Numerical solver results can differ slightly between platforms; hardware-dependent runtimes are not expected to match. The release verification does not rerun every archived experiment.

## Data dictionary

CSV headers identify all columns. `hour` is elapsed time from the beginning of the synthetic horizon, not a calendar timestamp. `pv_pu` and `wind_pu` are dimensionless capacity factors; `wind_speed_ms` is m/s; generation, demand, charging, discharging, unserved load and curtailment are kW; stored energy is kWh. Annual exported profiles use a 600 kW combined nameplate rating and a 0.5 photovoltaic fraction for `generation_kw`. Other fractions are calculated as `600 * (alpha * pv_pu + (1-alpha) * wind_pu)`. The hourly step is 1 h; the quarter-hour step is 0.25 h. Dispatch state values refer to the beginning of each interval and use a cyclic boundary. Scenario files contain `hour,generation_kw,load_kw`.

The model uses idealised renewable profiles and storage operation with specified efficiencies and power limits. The results are conditional on these assumptions and do not establish site-specific engineering feasibility.
