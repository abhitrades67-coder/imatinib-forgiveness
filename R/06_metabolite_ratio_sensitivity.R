# =============================================================================
# 06_metabolite_ratio_sensitivity.R
# Sensitivity on the steady-state metabolite:parent ratio R.
# The base case calibrates R to the AUC ratio (~0.15; Gleevec label). The reported
# IRIS *trough* ratio is higher (~0.25; Larson 2008), because the slowly-eliminated
# metabolite is flatter than the parent. Because the metabolite reservoir scales
# with R, we recompute the metabolite forgiveness benefit (Table 5) at R = 0.25 to
# confirm the qualitative conclusion is not an artefact of the AUC-based ratio.
#
# Reuses the EXACT code path (same seed, population builder, forgiveness engine and
# consistently-scaled moiety threshold), changing only MET$ss_ratio.
# Run standalone:  Rscript R/06_metabolite_ratio_sensitivity.R
# or it is sourced (and executed once) by R/run_all.R.
# =============================================================================

# Load dependencies only if not already in the environment (e.g. standalone run).
if (!exists("build_population")) {
  source("R/00_config.R")
  source("R/01_model.R")
  source("R/02_population.R")
  source("R/03_forgiveness.R")
}

run_ratio_case <- function(R_value) {
  MET$ss_ratio <<- R_value          # affects kf in build_population AND moiety_threshold()
  pop <- build_population()
  rows <- list()
  for (d in SIM$doses_mg) {
    res <- run_individual_forgiveness(pop, d)
    tol <- res$tolerance
    for (rp in MET$rpf_sens) {
      tr <- tol[tol$dose_mg == d & tol$rpf == rp, ]
      rows[[paste(d, rp)]] <- data.frame(
        R = R_value, dose_mg = d, rpf = rp,
        mean_gain = mean(tr$metab_gain),
        pct_gain_ge1 = 100 * mean(tr$metab_gain >= 1))
    }
  }
  dplyr::bind_rows(rows)
}

run_metabolite_ratio_sensitivity <- function() {
  cat("\n=== Metabolite SS-ratio sensitivity (R = 0.15 base vs 0.25 trough-based) ===\n")
  R_keep <- MET$ss_ratio
  out <- dplyr::bind_rows(run_ratio_case(0.15), run_ratio_case(0.25))
  MET$ss_ratio <<- R_keep           # restore base case so downstream steps are unaffected
  print(out, row.names = FALSE)
  write.csv(out, file.path(DIRS$results, "metabolite_ratio_sensitivity.csv"), row.names = FALSE)
  cat(sprintf("Wrote %s\n", file.path(DIRS$results, "metabolite_ratio_sensitivity.csv")))
  invisible(out)
}

run_metabolite_ratio_sensitivity()
