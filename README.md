# HIV Kenya

Agent-based model of HIV transmission and intervention delivery in Kenya, built on [STIsim](https://github.com/starsimhub/stisim) and [Starsim](https://github.com/starsimhub/starsim). The model includes structured sexual networks with risk groups, HIV testing (FSW-targeted, general population, and opportunistic), ART, and PrEP.

Both Python and R interfaces are provided. The R interface uses [rstarsim](https://github.com/starsimhub/rstarsim), which calls the Python engine via reticulate.


## Prerequisites

### Python environment

The Python packages are required regardless of whether you use the R or Python interface:

```bash
pip install starsim stisim sciris
```

### R packages (for R users)

```r
install.packages(c("reticulate", "devtools"))
devtools::install_github("starsimhub/rstarsim")
```

On first use, `rstarsim` will set up a conda environment automatically if needed. To use an existing environment instead:

```r
library(starsim)
load_starsim("my_env_name")
```


## Quick start (R)

Set your working directory to the repo root, then:

```r
library(starsim)
load_starsim()

source("hiv_model.R")

# Create and run a single simulation
sim <- make_sim(verbose = 1/12)
sim$run()

# View results
df <- sim$to_df(resample = "year", use_years = TRUE, sep = ".")
head(df)
```

### Overriding defaults

`make_sim()` accepts any Sim parameter as a named argument:

```r
sim <- make_sim(
  rand_seed = 42L,
  stop      = 2025L,
  verbose   = -1L
)
```


## Running with calibration (R)

Calibration is run in Python (see below) and produces `results/kenya_hiv_calib.obj`. To apply calibrated parameters in R:

```r
source("hiv_model.R")

sim <- make_sim(verbose = 1/12)

# Load calibration and apply best-fit parameters
calib      <- sc$loadobj("results/kenya_hiv_calib.obj")
calib_pars <- calib$df$iloc[0L]$to_dict()
sim$init()
sim <- make_sim_pars(sim, calib_pars)

sim$run()
```

### Multi-run ensemble

Run multiple parameter sets from the calibration posterior:

```r
source("hiv_model.R")

sims <- run_msim(use_calib = TRUE, n_pars = 50L, do_save = TRUE)
# Saves results/msim.df

save_stats(sims)
# Saves results/epi_df.df and results/sw_df.df
```


## Calibration (Python)

Calibration uses [Optuna](https://optuna.org/) via `sti.Calibration` and should be run in Python:

```bash
python run_hiv_calibration.py
```

This fits 7 parameters against UNAIDS/national data in `data/kenya_hiv_calib.csv`:

| Parameter | Description | Range |
|-----------|-------------|-------|
| `hiv_beta_m2f` | Male-to-female transmission rate | 0.008 -- 0.02 |
| `hiv_eff_condom` | Condom effectiveness | 0.5 -- 0.95 |
| `nw_prop_f0` | Proportion of females in low-risk group | 0.55 -- 0.9 |
| `nw_prop_m0` | Proportion of males in low-risk group | 0.50 -- 0.9 |
| `nw_f1_conc` | Female mid-risk concurrency | 0.01 -- 0.2 |
| `nw_m1_conc` | Male mid-risk concurrency | 0.01 -- 0.2 |
| `nw_p_pair_form` | Partnership formation probability | 0.4 -- 0.9 |

Settings at the top of `run_hiv_calibration.py`:
- `n_trials`: number of Optuna trials (default 1000)
- `n_workers`: parallel workers (default 50)
- `debug = True`: quick local run with 2 trials / 1 worker

Outputs saved to `results/`:
- `kenya_hiv_calib.obj` -- calibration object (used by `make_sim_pars`)
- `kenya_hiv_calib_stats.df` -- result percentiles by year
- `kenya_hiv_par_stats.df` -- posterior parameter distributions


## Plotting results

### From Python

```bash
python plot_calibrations.py
```

### From R

Results are pandas DataFrames accessible via reticulate. To plot with ggplot2:

```r
library(ggplot2)

source("hiv_model.R")

sim <- make_sim(verbose = 1/12)
sim$run()
df <- sim$to_df(resample = "year", use_years = TRUE, sep = ".")
df <- as.data.frame(df)

# HIV prevalence over time
ggplot(df, aes(x = timevec)) +
  geom_line(aes(y = hiv.prevalence_15_49 * 100)) +
  labs(x = "Year", y = "HIV prevalence (%)", title = "Kenya HIV prevalence (15-49)")
```


## Repository structure

```
hiv_kenya/
  hiv_model.py              # Model definition (Python)
  hiv_model.R               # Model definition (R, calls Python via reticulate)
  run_hiv_calibration.py    # Calibration script (Python)
  plot_sims.py              # Plotting functions
  plot_calibrations.py      # Plot calibration results
  utils.py                  # Plotting utilities
  data/
    init_prev_hiv.csv       # Initial HIV prevalence by risk group/sex/SW status
    condom_use.csv           # Condom use by partnership type over time
    n_art.csv                # ART coverage (absolute numbers) by year
    n_vmmc.csv               # VMMC coverage by year
    kenya_hiv_calib.csv      # Calibration targets (UNAIDS data 1990-2024)
    kenya_age_1985.csv       # Initial age distribution
    kenya_asfr.csv           # Age-specific fertility rates
    kenya_deaths.csv         # Age/sex-specific mortality rates
    kenya_migration.csv      # Net migration data
  results/                   # Saved calibration and simulation outputs
  assets/                    # Fonts for plotting
```


## Key functions

| Function | File | Description |
|----------|------|-------------|
| `make_sim(...)` | `hiv_model.py/R` | Create a configured Kenya HIV simulation |
| `make_sim_pars(sim, calib_pars)` | `hiv_model.py/R` | Apply calibration parameters to a sim |
| `make_custom_interventions()` | `hiv_model.py/R` | Build testing, ART, and PrEP interventions |
| `run_msim(use_calib, n_pars)` | `hiv_model.py/R` | Run an ensemble of simulations |
| `save_stats(sims)` | `hiv_model.py/R` | Save age/sex stratified results |
| `run_calibration(n_trials, n_workers)` | `run_hiv_calibration.py` | Run Optuna calibration (Python only) |


## Python-to-R translation reference

| Python | R |
|--------|---|
| `import stisim as sti` | `sti <- import("stisim")` |
| `sti.Sim(...)` | `sti$Sim(...)` or `do.call(sti$Sim, args)` |
| `dict(a=1, b=2)` | `list(a = 1, b = 2)` |
| `sim.diseases.hiv` | `sim$diseases$hiv` |
| `~array` (bitwise NOT) | `!array` |
| `array1 & array2` | `array1 & array2` |
| `sc.thispath() / 'data'` | `file.path(getwd(), "data")` |
| `pd.read_csv(f).set_index('y')` | `pd$read_csv(f)$set_index("y")` |
