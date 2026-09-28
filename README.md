# hiv_kenya

Agent-based HIV transmission model for Kenya, built on
[STIsim](https://github.com/starsimhub/stisim) and
[Starsim](https://github.com/starsimhub/starsim). Sibling project to
[`hiv_zambia`](https://github.com/starsimhub/hiv_zambia).

Structured sexual network with FSW segmentation, three testing arms
(FSW-targeted, general population, low-CD4 opportunistic), ART with
future coverage targets, and PrEP.

See `CLAUDE.md` for scope and current state.

## Install

```bash
pip install -e .
```

or `bash install_python.sh`.

## Run a single sim

```bash
python hiv_model.py
```

## Calibration

Optuna calibration against UNAIDS/national HIV data in
`data/kenya_hiv_calib.csv`. Calibrates 7 parameters (dot-notation keys
routed by stisim's `default_build_fn`):

| Parameter                    | Range        |
|------------------------------|--------------|
| `hiv.beta_m2f`               | 0.008 – 0.02 |
| `hiv.eff_condom`             | 0.5 – 0.95   |
| `structuredsexual.prop_f0`   | 0.55 – 0.9   |
| `structuredsexual.prop_m0`   | 0.50 – 0.9   |
| `structuredsexual.f1_conc`   | 0.01 – 0.2   |
| `structuredsexual.m1_conc`   | 0.01 – 0.2   |
| `structuredsexual.p_pair_form` | 0.4 – 0.9  |

```bash
python run_hiv_calibration.py     # 1000 trials, 50 workers (edit at top)
python plot_calibrations.py       # figures/hiv_calib_*.png
```

Set `debug=True` at the top of `run_hiv_calibration.py` for a 2-trial
smoke run.

Outputs (`results/`): `kenya_hiv_calib.obj` (shrunk to top 500 draws),
`kenya_hiv_calib_stats.df` (percentiles by year),
`kenya_hiv_par_stats.df` (posterior parameter distributions).

## Tests

```bash
pytest -v
```

## Repository layout

```
hiv_kenya/
  hiv_model.py              # sim builder
  interventions.py          # testing + ART + PrEP
  run_hiv_calibration.py    # Optuna calibration entry point
  plot_sims.py              # single-sim / ensemble figure
  plot_calibrations.py      # calibration figure
  utils.py                  # plotting helpers
  install_python.sh         # pip install -e .
  data/
    init_prev_hiv.csv       # initial HIV prevalence by risk/sex/SW
    condom_use.csv          # condom use by partnership type
    n_art.csv               # ART coverage counts by year
    n_vmmc.csv              # VMMC coverage by year
    kenya_hiv_calib.csv     # UNAIDS/national targets 1990–2024
    kenya_age_1985.csv      # initial age distribution
    kenya_asfr.csv          # age-specific fertility
    kenya_deaths.csv        # age/sex mortality
    kenya_migration.csv     # net migration
  results/                  # calibration + sim outputs (gitignored bulk)
  assets/                   # plotting fonts
```
