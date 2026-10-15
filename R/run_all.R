# =============================================================================
# run_all.R  -- master script, revision 1.1.0 (response to reviewers)
# Usage (from the pipeline root):
#   Rscript R/run_all.R        full run, N = 1000 virtual patients
#   Rscript R/run_all.R 30     smoke test with 30 patients (few uncertainty reps)
# Every source of randomness is seeded; outputs are written to results/ and
# figures/ (journal artwork in figures/submission/).
# =============================================================================
t_start <- Sys.time()
source("R/00_config.R")
source("R/01_model.R")
source("R/02_population.R")
source("R/03_forgiveness.R")
source("R/04_figures.R")
source("R/05_qualification.R")
source("R/06_sensitivity.R")
source("R/07_individualized_threshold_check.R")
source("R/08_external_validation.R")
source("R/09_benchmark.R")
source("R/10_tables_numbers.R")

args <- commandArgs(trailingOnly = TRUE)
SMOKE <- length(args) >= 1
if (SMOKE) {
  SIM$n_subjects <- as.integer(args[1]); SIM$n_uncertainty <- 4; SIM$n_unc_subjects <- 50
  cat(">>> Smoke test: N =", SIM$n_subjects, "\n")
}
step <- function(msg) cat(sprintf("\n[%s | %.1f min] %s\n", format(Sys.time(), "%H:%M:%S"),
                                  as.numeric(difftime(Sys.time(), t_start, units = "mins")), msg))
wcsv <- function(x, f) write.csv(x, file.path(DIRS$results, f), row.names = FALSE)

# 1. Population and provenance ----------------------------------------------
step("Virtual population (base case)")
pop <- build_population("base")
wcsv(pop, "virtual_population.csv")
popt <- make_population_table(pop); print(popt)
prov <- make_provenance()

# 2. Consistency checks, run-in adequacy ----------------------------------------
step("Consistency checks against published exposure")
qual <- qualify_model(pop); print(qual$table)
typ <- data.frame(id = 1, ka = PK$ka, cl = PK$CL, vc = PK$Vc, q = PK$Q, vp = PK$Vp,
                  kf = MET$ss_ratio * log(2) / MET$t_half, kme = log(2) / MET$t_half, rpf = MET$rpf_base)
typ_trough <- min(ss_exposure(typ, 400)$trough)
ext <- external_checks(qual$per, typ_trough); print(ext)

step("Run-in adequacy")
pop90 <- build_population("base", t_half = 90)
runin <- bind_rows(runin_adequacy(pop, SIM$ss_days, label = "base"),
                   runin_adequacy(pop90, SIM$ss_days, label = "thalf90"),
                   runin_adequacy(pop90, SIM$ss_days_long, label = "thalf90"))
wcsv(runin, "runin_adequacy.csv"); print(runin)

# 3. Base-case missed-dose experiments (k = 0..10) ------------------------------
step("Base case: individual missed-dose experiments")
base_long <- rbindlist(lapply(SIM$doses_mg, function(d) run_forgiveness(pop, d, sweep = TRUE)))
fwrite(base_long, file.path(DIRS$results, "base_individual_long.csv"))
base_summary <- bind_rows(lapply(SIM$doses_mg, function(d) summarise_forgiveness(base_long[dose_mg == d])))
wcsv(base_summary, "base_summary.csv")
tol <- rbindlist(lapply(SIM$doses_mg, function(d) tolerance_table(base_long[dose_mg == d])))
fwrite(tol, file.path(DIRS$results, "individual_tolerance.csv"))
traj <- bind_rows(lapply(SIM$doses_mg, function(d) trajectory_by_k(base_long[dose_mg == d])))
wcsv(traj, "trajectory_by_k.csv")
sweep <- bind_rows(lapply(SIM$doses_mg, function(d) threshold_sweep(base_long[dose_mg == d])))
wcsv(sweep, "threshold_sweep.csv")
hd <- headroom_data(base_long[dose_mg == 400], ks = 1:3)
headroom <- hd %>% group_by(criterion, k) %>%
  summarise(median = median(h), q1 = quantile(h, .25), q3 = quantile(h, .75),
            p10 = quantile(h, .10), p90 = quantile(h, .90), .groups = "drop")
wcsv(headroom, "headroom_summary.csv"); print(headroom)
med_headroom <- sapply(SIM$doses_mg, function(d) median(base_long[dose_mg == d & k == 0]$nadir_parent) / EFF$thr_parent)
names(med_headroom) <- paste0(SIM$doses_mg, " mg median")
cohort <- base_long[, .(nadir_parent_mean = mean(nadir_parent), nadir_parent_sd = sd(nadir_parent),
                        nadir_moiety_mean = mean(`nadir_moiety_r0.5`), nadir_moiety_sd = sd(`nadir_moiety_r0.5`)),
                    by = .(dose_mg, k)][, rpf := MET$rpf_base]
wcsv(cohort, "cohort_summary_by_missed.csv")
print(t(base_summary[, c("dose_mg", "pct_at_target", "cond_pct_lose_after1_parent",
                         "cond_pct_gain_ge1_r0.1", "cond_pct_gain_ge1_r0.5", "cond_pct_gain_ge1_r1.0",
                         "all_pct_fail_after1_parent", "cond_mean_excess_tb_parent_168",
                         "cond_mean_excess_tb_moiety_r0.5_168", "cond_median_deficit_parent_168",
                         "cond_median_deficit_moiety_r0.5_168")]))

# 4. Adherence patterns (deterministic and stochastic) --------------------------
step("Adherence patterns")
pat <- run_adherence_patterns(pop)
fwrite(pat, file.path(DIRS$results, "adherence_pattern_patient.csv"))
pat_summary <- summarise_patterns(pat); print(pat_summary)

# 5. Omission behaviours --------------------------------------------------------
step("Omission behaviours (late intake, doubling)")
omit <- run_omission_behaviour(pop, base_long); wcsv(omit, "omission_behaviour.csv"); print(omit)
omit_prof <- omission_profiles(400); wcsv(omit_prof, "omission_profiles_typical.csv")

# 6. Sensitivity analyses ----------------------------------------------------------
step("Sensitivity analyses (structure, kinetics, variability, interactions)")
sens <- run_all_sensitivity()
base_rows <- base_summary %>% mutate(group = "Base case", description = "Base case", R_thr = MET$ss_ratio)
sens_all <- bind_rows(base_rows, sens)
wcsv(sens_all, "sensitivity_summary.csv")

# 7. Calibration anchors (exposure headroom) ---------------------------------------
step("Calibration-anchor sensitivity")
cal <- bind_rows(lapply(SIM$doses_mg, function(d)
  run_calibration_anchors(pop, d, targets = CAL$targets * d / 400)))
cal$anchor_400 <- rep(CAL$targets, length(SIM$doses_mg))
wcsv(cal, "calibration_anchors.csv")
print(cal[, c("dose_mg", "anchor", "median_ss_trough", "pct_at_target", "cond_pct_lose_after1_parent",
              "all_pct_fail_after1_parent", "cond_pct_gain_ge1_r0.5", "cond_pct_gain_ge1_r1.0")])

# 8. Parameter uncertainty ------------------------------------------------------------
step("Parameter-uncertainty propagation (sourced RSEs; 25% stress)")
u1 <- run_uncertainty(stress = FALSE); u2 <- run_uncertainty(stress = TRUE)
unc_rep <- bind_rows(u1$replicates, u2$replicates); unc_sum <- bind_rows(u1$summary, u2$summary)
wcsv(unc_rep, "uncertainty_replicates.csv"); wcsv(unc_sum, "uncertainty_summary.csv")
wcsv(UNC, "uncertainty_specification.csv"); print(as.data.frame(unc_sum))

# 9. Published-model benchmark --------------------------------------------------------
step("Benchmark against published parent-metabolite models")
bench <- run_benchmark(); print(bench$table)

# 10. Post hoc checks --------------------------------------------------------------------
step("Individualized-threshold check")
indiv <- individualized_threshold_check(DIRS$results)

# 11. Figures ---------------------------------------------------------------------------------
step("Figures")
# Main-text figures are numbered by first citation: Fig. 1 headroom, Fig. 2 benchmark
f1 <- fig_headroom(hd %>% filter(k <= 2), med_headroom); gsave(f1, "Fig1_headroom.png", h = 5); tsave(f1, "Fig1", h = 5)
f2 <- fig_benchmark(bench$profiles);            gsave(f2, "Fig2_benchmark.png", h = 5.2);  tsave(f2, "Fig2", h = 5.2)
fs1 <- fig_schematic();                         gsave(fs1, "FigS1_schematic.png", h = 4); tsave(fs1, "FigS1", h = 4)
fs2 <- fig_external(qual$per);                  gsave(fs2, "FigS2_external.png");          tsave(fs2, "FigS2")
fs3 <- fig_patterns(pat_summary);               gsave(fs3, "FigS3_patterns.png", w = 9.5, h = 7.5); tsave(fs3, "FigS3", w = 9.5, h = 7.5)
fs4 <- fig_cohort_nadir(cohort);                gsave(fs4, "FigS4_nadir.png");             tsave(fs4, "FigS4")
fs5 <- fig_trajectory(traj);                    gsave(fs5, "FigS5_trajectory.png", h = 6.5); tsave(fs5, "FigS5", h = 6.5)
fs6 <- fig_sens_forest(sens_all);               gsave(fs6, "FigS6_sensitivity.png", w = 9, h = 6.5); tsave(fs6, "FigS6", w = 9, h = 6.5)
fs7 <- fig_threshold_sweep(sweep);              gsave(fs7, "FigS7_threshold.png");         tsave(fs7, "FigS7")
fs8 <- fig_calibration(cal);                    gsave(fs8, "FigS8_calibration.png");       tsave(fs8, "FigS8")
fs9 <- fig_omission(omit_prof);                 gsave(fs9, "FigS9_omission.png");          tsave(fs9, "FigS9")
fs10 <- fig_uncertainty(unc_rep);               gsave(fs10, "FigS10_uncertainty.png", h = 6); tsave(fs10, "FigS10", h = 6)

# 12. Numbers for the article ------------------------------------------------------------------
step("Numbers file")
v1_ref <- c(pct_at_target_400 = 55.7, cond_lose_400 = 78.1, pct_at_target_600 = 83.1, cond_lose_600 = 63.1)
sub <- sens_all %>% filter(label == "submitted")
extra <- data.frame(key = c("n_subjects", "typ_trough_400", "med_headroom_400", "med_headroom_600",
                            "spearman_trough_auc_400", "pct_below1000_obs_400", "pct_window_obs_400",
                            "median_trough_obs_400", "q1_trough_obs_400", "q3_trough_obs_400",
                            "median_trough_obs_600", "median_trough_400", "median_trough_600",
                            "q1_trough_400", "q3_trough_400",
                            "mean_trough_ratio_400", "sd_trough_ratio_400", "median_ratio_3to8h_400",
                            "median_troughm_400", "median_cmax_400", "median_auc_400", "median_t12",
                            "median_tmax_400", "median_auc_ratio_400", "mean_troughm_400", "cv_troughm_400",
                            "sub_reproduces_v1_400", "sub_reproduces_v1_600", "elapsed_min"),
  value = c(nrow(pop), typ_trough, med_headroom[1], med_headroom[2],
            cor(qual$per[dose_mg == 400]$trough, qual$per[dose_mg == 400]$AUC24, method = "spearman"),
            100 * mean(qual$per[dose_mg == 400]$trough_obs < 1000),
            100 * mean(qual$per[dose_mg == 400]$trough_obs >= 860 & qual$per[dose_mg == 400]$trough_obs <= 1500),
            median(qual$per[dose_mg == 400]$trough_obs), quantile(qual$per[dose_mg == 400]$trough_obs, .25),
            quantile(qual$per[dose_mg == 400]$trough_obs, .75), median(qual$per[dose_mg == 600]$trough_obs),
            median(qual$per[dose_mg == 400]$trough), median(qual$per[dose_mg == 600]$trough),
            quantile(qual$per[dose_mg == 400]$trough, .25), quantile(qual$per[dose_mg == 400]$trough, .75),
            mean(qual$per[dose_mg == 400]$trough_ratio), sd(qual$per[dose_mg == 400]$trough_ratio),
            median(qual$per[dose_mg == 400]$ratio_3to8h), median(qual$per[dose_mg == 400]$troughm),
            median(qual$per[dose_mg == 400]$Cmax), median(qual$per[dose_mg == 400]$AUC24),
            median(qual$per[dose_mg == 400]$t12), median(qual$per[dose_mg == 400]$Tmax),
            median(qual$per[dose_mg == 400]$met_ratio_auc),
            mean(qual$per[dose_mg == 400]$troughm),
            100 * sd(qual$per[dose_mg == 400]$troughm) / mean(qual$per[dose_mg == 400]$troughm),
            as.numeric(abs(sub$pct_at_target[sub$dose_mg == 400] - v1_ref[1]) < 0.05 &
                       abs(sub$cond_pct_lose_after1_parent[sub$dose_mg == 400] - v1_ref[2]) < 0.05),
            as.numeric(abs(sub$pct_at_target[sub$dose_mg == 600] - v1_ref[3]) < 0.05 &
                       abs(sub$cond_pct_lose_after1_parent[sub$dose_mg == 600] - v1_ref[4]) < 0.05),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
kv <- write_numbers(list(base_summary = base_summary, sens = sens_all, patterns = pat_summary, cal = cal,
                         sweep = sweep, unc_summary = unc_sum, traj = traj, omit = omit, bench = bench$table,
                         runin = runin, headroom = headroom, indiv = indiv, extra = extra))
cat(sprintf("  %d numbers written\n", nrow(kv)))

# 13. Reproducibility ------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(DIRS$results, "sessionInfo.txt"))
pkgs <- c("rxode2", "dplyr", "tidyr", "ggplot2", "data.table")
wcsv(data.frame(package = pkgs, version = sapply(pkgs, function(p) as.character(packageVersion(p)))),
     "package_versions.csv")
step(sprintf("Done. Results in %s/, figures in %s/", DIRS$results, DIRS$figures))
