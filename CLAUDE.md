# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run the model (single sim)
python hiv_model.py

# Run tests (Python)
cd tests && python test_model.py

# Run tests (R) -- from repo root
Rscript tests/test_model.R

# Run calibration (slow; use debug=True in the file for quick local runs)
python run_hiv_calibration.py

# Plot calibration results
python plot_calibrations.py
```

**Install dependencies:**
```bash
pip install starsim stisim sciris
```

## Architecture

This is an agent-based HIV transmission model for Kenya, built on [STIsim](https://github.com/starsimhub/stisim)/[Starsim](https://github.com/starsimhub/starsim).

### Core model (`hiv_model.py`)

- **`make_sim(**kwargs)`** — entry point; creates a `sti.Sim` configured for Kenya. Auto-loads `init_prev` and `condom_use` from `data/` via `DataLoader`. Custom interventions (testing, ART, PrEP) are always added, and user-provided `interventions`/`analyzers` kwargs are merged in.
- **`make_custom_interventions()`** — builds FSW-targeted testing, general-population testing, low-CD4 opportunistic testing, ART (with future coverage), and PrEP. Testing coverage scales linearly from 1990 to 2020, then continues to 2050.
- **`make_sim_pars(sim, calib_pars)`** — applies calibrated parameters to a sim; parameters prefixed `hiv_` route to `sim.diseases.hiv.pars`, those prefixed `nw_` route to `sim.networks.structuredsexual.pars`.
- **`run_msim(use_calib, n_pars)`** — runs an ensemble via `ss.parallel()`, optionally applying rows from the calibration posterior.
- **`save_stats(sims)`** — extracts age/sex stratified prevalence and incidence, plus SW stats, saving to `results/epi_df.df` and `results/sw_df.df`.

### Calibration (`run_hiv_calibration.py`)

Uses `sti.Calibration` (Optuna-based) to fit 7 parameters against UNAIDS/national data in `data/kenya_hiv_calib.csv`. Set `debug = True` at the top for a quick 2-trial local run. Outputs to `results/kenya_hiv_calib.obj`, `results/kenya_hiv_calib_stats.df`, and `results/kenya_hiv_par_stats.df`.

### R interface (`hiv_model.R`)

Wraps the Python model via `reticulate`/`rstarsim`. All key Python functions (`make_sim`, `make_sim_pars`, `run_msim`, `save_stats`) have R equivalents. Results (sciris `.df` objects) are pandas DataFrames accessible in R via reticulate.

### Parameter naming convention

Calibration parameters use prefixes that map to model components:
- `hiv_*` → `sim.diseases.hiv.pars` (e.g., `hiv_beta_m2f`, `hiv_eff_condom`)
- `nw_*` → `sim.networks.structuredsexual.pars` (e.g., `nw_prop_f0`, `nw_p_pair_form`)

### Results persistence

Simulation outputs are saved with `sc.saveobj()` as sciris binary objects (`.obj` or `.df` extension). Load with `sc.loadobj()`. DataFrames are resampled to yearly frequency via `sim.to_df(resample='year', use_years=True, sep='.')`.
