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

**Modernization in progress (branch `modernize`).** Bringing hiv_kenya
up to the hiv_zambia pattern: dot-notation calibration parameters via
stisim's `default_build_fn` (no custom `make_sim_pars`), interventions
split into `interventions.py`, ZAMPHIA-aligned 5-year `age_bins` on
`sti.HIV`. R interface and `sync-r-py` sync agent removed.

**Research question: TBD.** The current scope is calibration
modernization + a first Optuna fit against `data/kenya_hiv_calib.csv`.
Research question will be defined once a decent calibrated baseline is
in place.

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
