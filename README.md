# Imatinib metabolite-aware forgiveness analysis (CML)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21899206.svg)](https://doi.org/10.5281/zenodo.21899206)

**Authors:** Abhishek Anil, Aswini S, Archana A. S., Debasish Hota, Anand Srinivasan
**Version:** 1.0.0 · **License:** MIT · **Contact:** anandsrinivasan@aiimsbhubaneswar.edu.in
**Archived release:** [10.5281/zenodo.21899206](https://doi.org/10.5281/zenodo.21899206) (v1.0.0)
· concept DOI [10.5281/zenodo.21899205](https://doi.org/10.5281/zenodo.21899205) always resolves
to the newest version.

This archive is the analysis code and generated result datasets underlying the article
*"Does the long-lived active metabolite make imatinib more forgiving? A metabolite-aware
pharmacokinetic simulation of missed doses in chronic myeloid leukemia."* Please cite the
deposit via its Zenodo DOI, 10.5281/zenodo.21899206 (see `CITATION.cff`).

A from-scratch, simulation-based pharmacokinetic **forgiveness** analysis of imatinib in chronic
myeloid leukemia that explicitly models the **long-lived active metabolite N-desmethyl-imatinib
(CGP74588, t½ ≈ 40 h)** alongside the parent drug (t½ ≈ 18 h), and asks whether that metabolite
makes imatinib more tolerant of missed doses.

## Headline finding
When the metabolite's contribution is scored **consistently** (active-moiety threshold scaled so
parent and moiety criteria coincide at steady state — no double-counting), the long-lived metabolite
adds only a **modest, potency-dependent** amount of forgiveness. At the standard **400 mg** dose
imatinib is intrinsically unforgiving (≈88% of patients cannot miss even one dose without the nadir
dropping below the 1000 ng/mL target at the base-case calibration, ranging from ~60% to ~96% across
plausible population calibrations), and the metabolite does not change that. The benefit is
larger at **600 mg**, where the parent already provides buffer.

## How to run
Everything is reproduced from a single command (R 4.5.2; exact package versions pinned in `renv.lock`):

```r
# from the project root
Rscript R/run_all.R          # full run, N = 1000 virtual patients (~50 min)
Rscript R/run_all.R 50       # quick smoke test with 50 patients (a few minutes)
```

## Reproducibility
The analysis is **fully deterministic**: every source of randomness is seeded (the virtual
population, the delayed-dose scenario, and the 200 parameter-uncertainty replicates), so a re-run
regenerates every result **byte-for-byte identical**. This was verified end to end on
9 August 2026: a clean `Rscript R/run_all.R` reproduced all 20 CSV/text outputs in `results/`,
all 14 figures in `figures/` and all seven journal artwork files byte-for-byte. To recreate the exact software environment,
install [renv](https://rstudio.github.io/renv/) and run `renv::restore()` (reads `renv.lock`: R 4.5.2,
rxode2 5.0.1, dplyr, tidyr, ggplot2, knitr, data.table, rmarkdown and their dependencies).
`results/sessionInfo.txt` and `results/package_versions.csv` record the environment used to produce
the archived outputs.

## Code (`R/`)
| File | Purpose |
|------|---------|
| `00_config.R` | All evidence-based parameters, **each cited inline** (parent PK, metabolite, efficacy target, covariates, DDI states, uncertainty settings). |
| `01_model.R` | rxode2 structural model: 2-compartment parent + concentration-driven metabolite (identifiability-safe). |
| `02_population.R` | Virtual CML population (with metabolite IIV) + five adherence scenarios at 400/600 mg; per-day-minimum troughs; attainment across thresholds. |
| `03_forgiveness.R` | Individual-level forgiveness (parent vs active moiety, RPF range), threshold sensitivity, and exploratory CYP3A4 DDI. |
| `04_figures.R` | Publication + supplementary figures, model schematic, graphical abstract. |
| `05_qualification.R` | Model qualification (VPC-style) and parameter-uncertainty propagation. |
| `06_metabolite_ratio_sensitivity.R` | Sensitivity of the metabolite forgiveness benefit to the steady-state ratio R (AUC-based 0.15 vs trough-based 0.25). |
| `07_individualized_threshold_check.R` | Robustness check (post-hoc, reads saved results): per-patient vs population-R moiety threshold — benefit moves ≤1 pp, no conclusion changes. |
| `08_external_validation.R` | External validation (reads saved qualification output): simulated trough distribution and target attainment vs independent published CML cohorts (IRIS/Larson 2008; Picard 2007; Huang 2023 real-world TDM). |
| `run_all.R` | Orchestrates the whole pipeline and prints headline numbers. |
| `regen_dose_parity.R` | Maintenance-only utility (not part of the pipeline, not sourced by `run_all.R`): re-runs qualification and overwrites `calibration_sensitivity.csv` and three figures. Not needed to reproduce any reported result. |

### Manuscript figure ↔ file mapping
The manuscript numbers figures by order of first appearance, which does **not** match the
generation order in `04_figures.R`:

| Manuscript | File in `figures/` | Produced by |
|---|---|---|
| Fig. 1 (model schematic) | `fig0_model_schematic.png` | `04_figures.R` |
| Fig. 2 (external validation) | `fig13_external_validation.png` / `.tiff` | `08_external_validation.R` |
| Fig. 3 (target attainment) | `fig2_target_attainment.png` | `04_figures.R` |
| Fig. 4 (nadir vs missed doses) | `fig5_cohort_nadir.png` | `04_figures.R` |
| Fig. 5 (dose-holiday time course) | `fig1_concentration_profile.png` | `04_figures.R` |
| Fig. 6 (potency sensitivity) | `fig6_potency_sensitivity.png` | `04_figures.R` |

`fig7`–`fig12` and `graphical_abstract.png` are additional analyses retained in this archive;
they are not numbered figures in the article.

## Outputs
- `results/` — virtual population, scenario metrics, individual tolerance, cohort summaries, model
  qualification, threshold sensitivity, exploratory DDI, uncertainty replicates/summary, package
  versions, `sessionInfo.txt`, and a saved `analysis_objects.rds`.
- `figures/` — the six numbered manuscript figures plus additional analyses (model schematic,
  concentration profile, target attainment, nadir curves, potency sensitivity, threshold sensitivity,
  model qualification, exploratory DDI, parameter uncertainty, calibration sensitivity,
  time-below-threshold, external validation) plus `graphical_abstract.png`.

The article itself is published separately by the journal; this archive contains the analysis code
and the outputs needed to reproduce every reported result. The article carries no supplementary
files, so the additional analyses above are available only here.

## Key methodological choices (evidence base)
- **Efficacy anchored to PARENT trough ≥ 1000 ng/mL** (Picard 2007, threshold 1002 ng/mL).
- **Metabolite parameterised** so its 40 h half-life and ~0.15 steady-state metabolite:parent ratio
  match literature exactly, folding the unidentifiable formation-fraction × bioavailability product
  into the calibrated ratio.
- **Relative potency (RPF) treated as a sensitivity axis (0.1–1.0)** because it is contested:
  "comparable" (FDA label) vs "significantly less active" (Manley et al. 2013).
- **Steady-state clearance (~10.5 L/h)** used to reflect imatinib auto-inhibition (Schmidli 2005);
  disposition parameters (Vc/F, Vp/F, Q/F, kₐ) **calibrated** to reproduce the observed steady-state
  400 mg profile (see manuscript Table 1).

All parameter choices are documented with inline citations in `R/00_config.R`; the full
bibliography accompanies the article rather than this archive.

## Note
This is **simulation only** — no individual patient data are fitted; inter-individual variability is
imposed from published estimates. See the manuscript Discussion for the full limitations.

## License
Released under the MIT License (see `LICENSE`).
