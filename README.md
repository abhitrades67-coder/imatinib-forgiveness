# Imatinib metabolite-aware forgiveness analysis (CML), version 1.1.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21899205.svg)](https://doi.org/10.5281/zenodo.21899205)

**Authors:** Abhishek Anil, Aswini Saravanan, Archana A. S., Debasish Hota, Anand Srinivasan
**Version:** 1.1.0 (revision in response to peer review, October 2026) · **License:** MIT ·
**Contact:** anandsrinivasan@aiimsbhubaneswar.edu.in
**Archive:** Zenodo concept DOI [10.5281/zenodo.21899205](https://doi.org/10.5281/zenodo.21899205)
(resolves to the newest version). Version 1.0.0
([10.5281/zenodo.21899206](https://doi.org/10.5281/zenodo.21899206)) reproduces the originally
submitted analysis.

This archive contains the analysis code and generated result datasets for the article
*"The active metabolite N-desmethyl-imatinib adds little to imatinib forgiveness of missed
doses: a pharmacokinetic simulation in chronic myeloid leukemia"* (Cancer
Chemotherapy and Pharmacology, revised manuscript).

## What the analysis does
A literature-assembled simulation framework (two-compartment parent imatinib with first-order
absorption, coupled to the active metabolite N-desmethyl-imatinib, CGP74588) is applied to 1000
virtual chronic-phase CML patients at 400 and 600 mg once daily. After a 21-day run-in, 0-10
consecutive doses are omitted. The metabolite is scored on an active-moiety criterion whose target
is scaled consistently with the parent target, so that it can contribute only through its slower
decay (no double counting). Endpoints: maintenance of the population-derived 1000 ng/mL trough
target (nadir criterion), excess time and area below target, exposure deficit (48 h to 28 days),
and the exposure headroom needed to bridge missed doses. The primary analysis is restricted to
patients at target with uninterrupted dosing.

## How to run
```r
# from this folder (the pipeline root)
Rscript R/run_all.R        # full run, N = 1000 virtual patients (about 20 minutes on 8 cores)
Rscript R/run_all.R 30     # smoke test with 30 patients (about 2 minutes)
```
R 4.5.2 with rxode2 5.0.1, dplyr, tidyr, ggplot2 and data.table (versions in
`results/package_versions.csv` and `results/sessionInfo.txt`; `renv.lock` pins the environment:
install renv and run `renv::restore()`). The analysis is deterministic: every source of randomness
is seeded (population `20260606`; delayed-dose days `+7`; random adherence patterns `+12...+15`;
residual error `+31`; inter-occasion variability `+41`; uncertainty `+99` and `+199`).

## Code (`R/`)
| File | Purpose |
|---|---|
| `00_config.R` | Every parameter with its source and status (adopted / calibrated / assumed); uncertainty table with the source of each RSE; sensitivity settings. |
| `01_model.R` | rxode2 models: base case (linear formation), Michaelis-Menten formation, Menon-Andersen two-compartment metabolite. |
| `02_population.R` | Virtual population (base case: Schmidli 2005 weight exponents and CL-Vc covariance; conventional allometric variant; IIV on Q/Vp variant); deterministic and random (independent / Markov-clustered) adherence patterns after the run-in. |
| `03_forgiveness.R` | Missed-dose experiments; nadir and trajectory metrics; conditional and whole-cohort summaries; threshold sweep; Michaelis-Menten normalisation; inter-occasion variability. |
| `04_figures.R` | Figures for the article (main and supplementary). |
| `05_qualification.R` | Consistency checks against published exposure (with residual error), run-in adequacy, calibration anchors, parameter-uncertainty propagation. |
| `06_sensitivity.R` | Structural, kinetic, variability, covariate, saturable-formation and exploratory CYP3A4 (parent-only) sensitivity analyses; omission behaviour (late intake, doubling). |
| `07_individualized_threshold_check.R` | Per-patient versus population-ratio scaling of the moiety target (post hoc). |
| `08_external_validation.R` | External consistency figure (IRIS, TOPS, real-world window). |
| `09_benchmark.R` | Typical-patient benchmark against the published parent-metabolite models (Delbaldo 2006, Petain 2008, Menon-Andersen 2009). |
| `10_tables_numbers.R` | Parameter-provenance table, population table, pattern summary and `results/manuscript_numbers.csv` (every number quoted in the article). |
| `run_all.R` | Runs everything in order. |

## Outputs and article items
| Article item | Output |
|---|---|
| Table 1 | `results/parameter_provenance.csv` |
| Table 2 | `results/adherence_pattern_summary.csv` |
| Table 3 | `results/base_summary.csv`, `results/trajectory_by_k.csv` |
| Table 4 | `results/sensitivity_summary.csv`, `threshold_sweep.csv`, `calibration_anchors.csv`, `uncertainty_summary.csv` |
| Fig. 1 | `figures/Fig1_headroom.png` (`headroom_summary.csv`) |
| Fig. 2 | `figures/Fig2_benchmark.png` (`benchmark_profiles.csv`, `benchmark_published_models.csv`) |
| Tables S1-S11 | `benchmark_published_models.csv`, `population_characteristics.csv`, (S3 is descriptive), `runin_adequacy.csv`, `uncertainty_specification.csv`, `model_qualification.csv`, `external_consistency.csv`, `sensitivity_summary.csv`, `threshold_sweep.csv`, `calibration_anchors.csv`, `omission_behaviour.csv` |
| Figs S1-S10 | `figures/FigS1_*.png` ... `figures/FigS10_*.png` |
| Journal artwork | `figures/submission/Fig1.tiff` ... `FigS10.tiff` (1000 dpi, LZW) and matching PNGs |

Per-patient results: `base_individual_long.csv` (every patient, dose and number of omitted doses),
`individual_tolerance.csv`, `adherence_pattern_patient.csv`, `qualification_per_subject.csv`,
`sensitivity/*.rds` (per-patient results of each sensitivity configuration).

## Changes from version 1.0.0
- Base-case covariate model: Schmidli 2005 weight exponents (0.301, 0.405), no sex effect, CL-Vc
  covariance 0.071 (r = 0.707). The v1.0.0 model is the `submitted` sensitivity configuration and
  reproduces the v1.0.0 results exactly (e.g. 55.7% at target and 78.1% losing the target after one
  missed dose at 400 mg).
- Conditional (at-target) endpoints made primary; trajectory metrics (excess time and area below
  target; exposure deficit) and exposure headroom added.
- Adherence patterns now start after the run-in; random (independent / clustered) patterns added.
- New sensitivity analyses: published metabolite structures, metabolite half-life 20/90 h,
  Michaelis-Menten formation, IIV on Q/Vp, inter-occasion variability, metabolite IIV x0/x2,
  late intake and doubling, threshold sweep 500-1500 ng/mL; calibration anchors with conditional
  metrics; uncertainty propagation with sourced RSEs plus a stress test.
- CYP3A4 interaction analysis is parent-only (the v1.0.0 metabolite-formation multipliers were
  removed; the reported v1.0.0 interaction results were parent-only and are unaffected).

## Note
Simulation only: no individual patient data are fitted. The framework is assembled from published
estimates, with calibrated disposition parameters and assumed metabolite variability; comparisons
with published exposure are consistency checks, not validation of the missed-dose trajectories.

## License
MIT (see `LICENSE`).
