# =============================================================================
# 10_tables_numbers.R
# Parameter-provenance table (Table 1), population characteristics (Table S1),
# adherence-pattern summary, and a key-value file of every number quoted in the
# revised article (results/manuscript_numbers.csv), so that the text can be
# regenerated from the results without manual transcription.
# =============================================================================
suppressMessages({library(dplyr); library(data.table)})

make_provenance <- function() {
  tab <- data.frame(
    Parameter = c("CL/F (L/h)", "Vc/F (L)", "Vp/F (L)", "Q/F (L/h)", "ka (1/h)",
                  "IIV CL/F; Vc/F (omega)", "Correlation eta(CL)-eta(Vc)", "IIV ka (omega)",
                  "IIV Q/F, Vp/F", "Inter-occasion variability", "Residual error",
                  "Weight effect on CL/F, Q/F; on Vc/F, Vp/F", "Sex effect on CL/F",
                  "Metabolite t1/2 (h)", "Metabolite:parent ratio R", "IIV R; IIV metabolite t1/2",
                  "Metabolite formation", "Relative potency (RPF)", "Population-derived trough target (ng/mL)"),
    Value = c("10.5", "160", "110", "50", "0.5",
              "0.32; 0.31", "0.707", "0.40",
              "none (base case)", "none (base case)", "none (true concentrations)",
              "(WT/70)^0.301; (WT/70)^0.405", "none",
              "40", "0.15", "0.30; 0.20",
              "first order, kf = R x kme", "0.1, 0.5, 1.0", "1000"),
    Status = c("Adopted", "Calibrated", "Calibrated", "Calibrated", "Calibrated",
               "Adopted", "Adopted", "Assumed",
               "Sensitivity analysis", "Sensitivity analysis", "Consistency checks only",
               "Adopted", "Adopted (no effect)",
               "Adopted", "Adopted", "Assumed",
               "Assumed; tested", "Assumption (sensitivity axis)", "Adopted"),
    Source = c("Schmidli 2005 (day-29 value 9.6 L/h at 70 kg); Golabchifar 2014 (10.8 L/h)",
               "Calibrated to the 400 mg steady-state profile and terminal half-life; cf. 139 (Golabchifar 2016)",
               "As above; cf. 130 (Golabchifar 2016)",
               "As above; cf. 29.6 (Golabchifar 2016)",
               "As above; cf. 0.94 (Golabchifar 2016) and 1.1 (Schmidli 2005)",
               "Schmidli 2005: omega^2 0.102 and 0.099",
               "Schmidli 2005: covariance 0.071 (applied to Vc/F)",
               "Golabchifar 2016: 31.2% on ka, 47-48% on its zero-order input terms",
               "Golabchifar 2016: 41.8% and 29.7%, correlation 0.921",
               "F per dose, CV 28% (TOPS intra-patient trough CV 27.8%) and 44% (Menon-Andersen 2009)",
               "Schmidli 2005: 26% proportional; added only when compared with measured concentrations",
               "Schmidli 2005 final model (exponents 0.301, SE 0.096; 0.405, SE 0.118)",
               "Sex not retained in the Schmidli 2005 final model",
               "Prescribing information; 20 h (Gschwind 2005) and 90 h (le Coutre 2004) as sensitivity analyses",
               "Prescribing-information AUC ratio; 0.25 (trough-based) as sensitivity analysis",
               "No published between-subject estimate; IRIS trough-ratio CV 33% (Larson 2008)",
               "Linear in TOPS (dose-proportional CGP74588 troughs) and in three published models; Michaelis-Menten sensitivity (Nebot 2010)",
               "Prescribing information ('comparable') versus Manley 2013 ('significantly less active')",
               "Picard 2007 (1002 ng/mL); OPTIM-imatinib; IATDMCT consensus (Clarke 2021); swept 500-1500"),
    Uncertainty_RSE = c("7% (Golabchifar 2016)", "4% (Golabchifar 2016); 25% stress", "9% (Golabchifar 2016); 25% stress",
                        "11% (Golabchifar 2016); 25% stress", "10% (Golabchifar 2016); 25% stress",
                        "14%; 16% (Schmidli 2005)", "15% on covariance (Schmidli 2005)", "-",
                        "-", "-", "-", "SE 0.096; 0.118 (Schmidli 2005)", "-",
                        "10% (assumed)", "15% (assumed)", "-", "-", "-", "-"),
    check.names = FALSE, stringsAsFactors = FALSE)
  write.csv(tab, file.path(DIRS$results, "parameter_provenance.csv"), row.names = FALSE)
  tab
}

make_population_table <- function(pop) {
  q <- function(x, f = "%.1f") sprintf(paste0(f, " [", f, "-", f, "]"), median(x), quantile(x, .25), quantile(x, .75))
  t12 <- thalf_terminal(pop$cl, pop$vc, pop$vp, pop$q)
  tab <- data.frame(
    Characteristic = c("Age (years)", "Body weight (kg)", "Female (%)", "CL/F (L/h)", "Vc/F (L)",
                       "Vp/F (L)", "Q/F (L/h)", "ka (1/h)", "Parent terminal half-life (h)",
                       "Metabolite:parent ratio R", "Metabolite half-life (h)"),
    Value = c(q(pop$age, "%.0f"), q(pop$wt, "%.0f"), sprintf("%.1f", 100 * mean(pop$sex == 1)),
              q(pop$cl), q(pop$vc, "%.0f"), q(pop$vp, "%.0f"), q(pop$q, "%.0f"), q(pop$ka, "%.2f"),
              q(t12), q(pop$met_ratio, "%.3f"), q(pop$met_thalf)),
    stringsAsFactors = FALSE)
  write.csv(tab, file.path(DIRS$results, "population_characteristics.csv"), row.names = FALSE)
  tab
}

summarise_patterns <- function(pat) {
  thr <- EFF$thr_parent
  perf <- pat[scenario == "perfect", .(id, dose_mg, at = mean_trough_parent >= thr,
                                       auc0 = auc_parent, aucm0 = auc_metab)]
  d <- merge(pat, perf, by = c("id", "dose_mg"))
  d[, deficit_parent := 100 * (1 - auc_parent / auc0)]
  out <- d[, .(n = .N, adherence = 100 * mean(adherence),
               median_trough = median(mean_trough_parent),
               pi_lo = quantile(mean_trough_parent, .025), pi_hi = quantile(mean_trough_parent, .975),
               median_trough_metab = median(mean_trough_metab),
               pct_ge_target = 100 * mean(mean_trough_parent >= thr),
               cond_pct_days_below_parent = mean(pct_days_below_parent[at]),
               cond_pct_time_below_parent = mean(pct_time_below_parent[at]),
               `cond_pct_time_below_moiety_r0.5` = mean(`pct_time_below_moiety_r0.5`[at]),
               `cond_pct_time_below_moiety_r1.0` = mean(`pct_time_below_moiety_r1.0`[at]),
               cond_median_deficit_parent = median(deficit_parent[at]),
               all_pct_time_below_parent = mean(pct_time_below_parent)),
           by = .(dose_mg, scenario)]
  write.csv(out, file.path(DIRS$results, "adherence_pattern_summary.csv"), row.names = FALSE)
  out
}

# ----------------------------------------------------------------------------
# Key-value numbers file. Keys are built systematically:
#   base_<column>_<dose>          base-case summary (summarise_forgiveness)
#   sens_<config>_<column>_<dose> sensitivity configurations
#   pat_<scenario>_<column>_<dose>, cal_<anchor>_<column>_<dose>,
#   sweep_<threshold>_<column>_<dose>, unc_<run>_<column>_<dose>,
#   traj_k<k>_<column>_<dose>, omit_<case>_<column>_<dose>, bench_<model>_<column>
# ----------------------------------------------------------------------------
kv_from <- function(df, prefix, id_cols, by_dose = TRUE) {
  df <- as.data.frame(df)
  num <- names(df)[sapply(df, is.numeric)]
  num <- setdiff(num, c(id_cols, "dose_mg"))
  out <- list()
  for (i in seq_len(nrow(df))) {
    idp <- paste(sapply(id_cols, function(c) gsub("[^A-Za-z0-9.]+", "", as.character(df[i, c]))), collapse = "_")
    for (c in num) {
      key <- paste0(prefix, if (nchar(idp)) paste0(idp, "_") else "", c,
                    if (by_dose && "dose_mg" %in% names(df)) paste0("_", df$dose_mg[i]) else "")
      out[[length(out) + 1]] <- data.frame(key = key, value = df[i, c])
    }
  }
  bind_rows(out)
}

write_numbers <- function(objs) {
  kv <- bind_rows(
    kv_from(objs$base_summary, "base_", character(0)),
    kv_from(objs$sens, "sens_", "label"),
    kv_from(objs$patterns, "pat_", "scenario"),
    kv_from(objs$cal, "cal_", "anchor"),
    kv_from(objs$sweep, "sweep_", "threshold"),
    kv_from(objs$unc_summary, "unc_", "run"),
    kv_from(objs$traj, "traj_", "k"),
    kv_from(objs$omit, "omit_", "behaviour"),
    kv_from(objs$bench, "bench_", "model", by_dose = FALSE),
    kv_from(objs$runin, "runin_", c("label", "run_in_days"), by_dose = FALSE),
    kv_from(objs$headroom, "headroom_", c("criterion", "k"), by_dose = FALSE),
    kv_from(objs$indiv, "indiv_", "rpf"),
    objs$extra)
  kv <- kv[!duplicated(kv$key), ]
  write.csv(kv, file.path(DIRS$results, "manuscript_numbers.csv"), row.names = FALSE)
  kv
}
