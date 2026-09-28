# CLAUDE.md

Agent-based HIV transmission model for Kenya. Sibling project to
`hiv_zambia`; same stisim base, same modernization pattern.

See `README.md` for install and run instructions.

## Repo layout

Root-level Python for a small sprint-scale project.

- `hiv_model.py` — model builder
- `interventions.py` — HIV testing (FSW / general population / low-CD4) + ART + PrEP
- `run_hiv_calibration.py` — Optuna calibration entry point
- `plot_sims.py` / `plot_calibrations.py` — figure scripts
- `utils.py` — plotting helpers
- `data/` — Kenya demography, condom use, ART counts, national HIV surveillance targets
- `results/` — calibration + sim outputs (gitignored bulk; committable summaries only)

## State of play

**Modernization complete + first calibration (branch `modernize`, 2026-09-28).**
- R interface removed (`hiv_model.R`, `plot_sims.R`, `install_R.sh`, `test_model.R`, `sync-r-py` sub-agent, R CI workflow).
- Dot-notation calibration parameters routed via stisim's `default_build_fn` (no custom `make_sim_pars`).
- Interventions split into `interventions.py` (FSW / general / low-CD4 HIV testing + ANC testing + ART with 2024+ 0.97 projection + PrEP).
- 5-year `age_bins` on `sti.HIV`.
- `plot_sims.py` rewritten to hiv_zambia's cleaner `_load_data` pattern (previous version read a non-existent `data/kenya_hiv_data.csv`).

**First calibration.** 1000-trial Optuna TPE (50 workers, ~5 min wall time), shrunk to top 500 draws. Mismatch of best trial: 21.77. Ensemble brackets UNAIDS at every panel (population, PLHIV, prevalence 15-49, new infections, HIV-related deaths, on ART) — see `figures/hiv_calib.png`. Posterior parameter summary:

| Parameter                        | Mean  | 5%–95%       |
|----------------------------------|-------|--------------|
| `hiv.beta_m2f`                   | 0.012 | 0.011–0.013  |
| `hiv.eff_condom`                 | 0.931 | 0.904–0.949  |
| `structuredsexual.prop_f0`       | 0.573 | 0.550–0.603  |
| `structuredsexual.prop_m0`       | 0.627 | 0.516–0.677  |
| `structuredsexual.f1_conc`       | 0.103 | 0.031–0.151  |
| `structuredsexual.m1_conc`       | 0.128 | 0.022–0.196  |
| `structuredsexual.p_pair_form`   | 0.747 | 0.437–0.864  |

`eff_condom` sits near the prior's upper edge — worth revisiting the range (currently 0.5–0.95) once a research question crystallises.

**Research question: TBD.** Calibrated baseline is ready; the research question and downstream analysis will be scoped in a follow-up session.

## Intake

**Model.** `sti.HIV` + `sti.StructuredSexual` (FSW-segmented) +
`MaternalNet`, `demographics='kenya'`, 10k agents, 1985 start.
Interventions: FSW / general-population / CD4 < 200 HIV testing arms
with historical scale-up curves; ART with `n_art.csv` historical
counts and a 0.97 projected proportion from 2024; PrEP scaling to 80%
by 2025.

**Question.** TBD. First deliverable is a calibrated baseline model
ready to build a research question on top of.

**Data.** `data/kenya_hiv_calib.csv` (UNAIDS/national surveillance
1990–2024) is the primary calibration target. Age × sex validation
data (KENPHIA / KDHS) not yet incorporated — flagged as a next step
once the research question requires it.

**Constraints.** Solo (Robyn). End-of-day 2026-09-28 for the modernized
model + first calibration.

## Environment

- stisim 1.7.0 (editable at `/home/robyn/stisim/`)
- starsim 3.6.1
- Python at `/home/robyn/miniconda/bin/python`; no conda env activation required

## Conventions

- Any modification to the editable `stisim` install is committed, pushed, and PR'd immediately per `stisim:editable-dep-hygiene`.
- Downstream vs upstream decisions run through `stisim:extending-stisim` — real bugs go upstream; opt-in project knobs stay downstream.
- Comment discipline in shared library code per `stisim:comment-hygiene` — no project-history in stisim source.

## Related projects

- [`hiv_zambia`](https://github.com/starsimhub/hiv_zambia) — sibling calibration + partner-notification analysis.
- [`hiv_fp_kenya`](https://github.com/starsimhub/hiv_fp_kenya) — FPsim + STIsim postpartum "one-stop shop" demo. Distinct research track; no shared code with this repo.
