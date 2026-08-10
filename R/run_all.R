# =============================================================================
# run_all.R  -- master script for the imatinib metabolite-aware forgiveness analysis
# Usage:  Rscript R/run_all.R           (full run, N = SIM$n_subjects)
#         Rscript R/run_all.R 30        (smoke test with 30 subjects)
# =============================================================================
t_start <- Sys.time()
source("R/00_config.R")
source("R/01_model.R")
source("R/02_population.R")
source("R/03_forgiveness.R")
source("R/04_figures.R")
source("R/05_qualification.R")
source("R/08_external_validation.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  SIM$n_subjects <- as.integer(args[1])
  SIM$n_uncertainty <- 10           # fast smoke test
  cat(">>> Override N =", SIM$n_subjects, "(smoke test; uncertainty reps = 10)\n")
}

# 1. Virtual population -------------------------------------------------------
pop <- build_population()
write.csv(pop, file.path(DIRS$results, "virtual_population.csv"), row.names = FALSE)
cat(sprintf("Population: N=%d, median CL/F=%.1f L/h, median wt=%.0f kg\n",
            nrow(pop), median(pop$cl), median(pop$wt)))

# 2. Adherence-scenario forgiveness ------------------------------------------
cat("\n[Adherence scenarios]\n")
scn <- run_all_scenarios(pop)
write.csv(scn$troughs, file.path(DIRS$results, "scenario_troughs_long.csv"), row.names = FALSE)
write.csv(scn$metrics, file.path(DIRS$results, "forgiveness_metrics.csv"), row.names = FALSE)
for (d in SIM$doses_mg) {
  write.csv(scn$metrics[scn$metrics$dose_mg == d, ],
            file.path(DIRS$results, sprintf("forgiveness_metrics_%dmg.csv", d)), row.names = FALSE)
}

# 3. Individual-level metabolite-aware forgiveness ---------------------------
cat("\n[Individual-level forgiveness]\n")
ind <- list()
for (d in SIM$doses_mg) ind[[as.character(d)]] <- run_individual_forgiveness(pop, d)

tol_all    <- bind_rows(lapply(ind, `[[`, "tolerance"))
cohort_all <- bind_rows(lapply(ind, `[[`, "cohort"))
long_all   <- bind_rows(lapply(ind, `[[`, "long"))
thr_all    <- bind_rows(lapply(ind, `[[`, "threshold"))
write.csv(tol_all,    file.path(DIRS$results, "individual_tolerance.csv"), row.names = FALSE)
write.csv(cohort_all, file.path(DIRS$results, "cohort_summary_by_missed.csv"), row.names = FALSE)
write.csv(long_all,   file.path(DIRS$results, "individual_nadir_long.csv"), row.names = FALSE)
write.csv(thr_all,    file.path(DIRS$results, "threshold_sensitivity.csv"), row.names = FALSE)

# 3b. Model qualification, DDI (exploratory) and parameter uncertainty ----------
cat("\n[Model qualification]\n")
qual <- qualify_model(pop)
print(qual$table)

cat("\n[External validation vs independent published cohorts]\n")
extval <- external_validation(qual$per)
print(extval$table)

cat("\n[Exploratory CYP3A4 DDI sensitivity]\n")
ddi <- bind_rows(lapply(SIM$doses_mg, function(d) run_ddi_sensitivity(pop, d)))
write.csv(ddi, file.path(DIRS$results, "ddi_exploratory.csv"), row.names = FALSE)
print(ddi)

cat("\n[Parameter-uncertainty propagation]\n")
unc_by_dose <- lapply(SIM$doses_mg, function(d) {
  cat(sprintf("  -- %d mg --\n", d)); u <- run_uncertainty(d); u$summary$dose_mg <- d; u
})
unc <- list(summary = bind_rows(lapply(unc_by_dose, `[[`, "summary")),
            replicates = bind_rows(lapply(seq_along(SIM$doses_mg), function(i) {
              r <- unc_by_dose[[i]]$replicates; r$dose_mg <- SIM$doses_mg[i]; r })))
write.csv(unc$summary,    file.path(DIRS$results, "uncertainty_summary.csv"), row.names = FALSE)
write.csv(unc$replicates, file.path(DIRS$results, "uncertainty_replicates.csv"), row.names = FALSE)
print(unc$summary)

cat("\n[Calibration-target sensitivity]\n")
cal <- bind_rows(lapply(SIM$doses_mg, function(d)
  run_calibration_sensitivity(pop, d, targets = CAL$targets * (d / 400))))
write.csv(cal, file.path(DIRS$results, "calibration_sensitivity.csv"), row.names = FALSE)
print(cal)

cat("\n[Metabolite-ratio sensitivity]\n")
source("R/06_metabolite_ratio_sensitivity.R")   # writes metabolite_ratio_sensitivity.csv; restores MET$ss_ratio

cat("\n[Individualized-threshold robustness check]\n")
# Sourced into its own environment: the script assigns to `pop`, `d`, `long`,
# `tol` and `summ`, which would otherwise clobber the objects used below.
source("R/07_individualized_threshold_check.R", local = new.env())

# 4. Headline numbers ---------------------------------------------------------
cat("\n================ HEADLINE RESULTS ================\n")
for (d in SIM$doses_mg) {
  tb <- tol_all[tol_all$dose_mg == d & tol_all$rpf == MET$rpf_base, ]   # base-case rpf
  cat(sprintf("\n-- %d mg QD  (base-case relative potency = %.1f) --\n", d, MET$rpf_base))
  cat(sprintf("  Max tolerable missed doses (PARENT-only)     : mean %.2f, median %g, range %g-%g\n",
              mean(tb$max_tol_parent), median(tb$max_tol_parent),
              min(tb$max_tol_parent), max(tb$max_tol_parent)))
  cat(sprintf("  Max tolerable missed doses (METABOLITE-aware): mean %.2f, median %g, range %g-%g\n",
              mean(tb$max_tol_moiety), median(tb$max_tol_moiety),
              min(tb$max_tol_moiety), max(tb$max_tol_moiety)))
  cat(sprintf("  Patients tolerating 0 missed doses (parent)  : %.1f%%\n",
              100 * mean(tb$max_tol_parent == 0)))
  cat("  Metabolite forgiveness gain across potency range:\n")
  for (rp in MET$rpf_sens) {
    tr <- tol_all[tol_all$dose_mg == d & tol_all$rpf == rp, ]
    cat(sprintf("     RPF=%.1f : mean gain %.2f doses, %.1f%% gain >=1 extra dose\n",
                rp, mean(tr$metab_gain), 100 * mean(tr$metab_gain >= 1)))
  }
}

# 5. Figures ------------------------------------------------------------------
cat("\n[Figures]\n")
make_all_figures(pop, scn$metrics, tol_all, cohort_all,
                 thr_all = thr_all, qual = qual$per, ddi = ddi, unc = unc$replicates,
                 cal = cal)

# 6. Reproducibility (B8) ------------------------------------------------------
si <- capture.output(sessionInfo())
writeLines(si, file.path(DIRS$results, "sessionInfo.txt"))
pkgs <- c("rxode2","dplyr","tidyr","ggplot2","knitr","data.table")
pv <- data.frame(package = pkgs,
                 version = sapply(pkgs, function(p) as.character(packageVersion(p))))
write.csv(pv, file.path(DIRS$results, "package_versions.csv"), row.names = FALSE)

saveRDS(list(pop = pop, scenarios = scn, individual = ind,
             qualification = qual, ddi = ddi, uncertainty = unc),
        file.path(DIRS$results, "analysis_objects.rds"))
cat(sprintf("\nDone in %.1f min. Results in %s/\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")), DIRS$results))
