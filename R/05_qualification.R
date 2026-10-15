# =============================================================================
# 05_qualification.R
# Consistency checks of the assembled framework against published exposure
# (described as CONSISTENCY CHECKS, not validation of the missed-dose
# trajectories), run-in adequacy, parameter-uncertainty propagation with sourced
# RSEs, and the calibration-anchor (exposure-headroom) sensitivity analysis.
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(data.table)})

# ----------------------------------------------------------------------------
# Per-subject steady-state exposure (day 30 interval of uninterrupted dosing).
# ----------------------------------------------------------------------------
ss_exposure <- function(pop, dose_mg, model_key = "base") {
  tau <- SIM$tau
  ev <- shared_events((0:29) * tau, dose_mg, seq(29 * tau, 30 * tau, by = 0.25))
  s <- solve_pop(model_key, pop[, MODEL_PARAMS[[model_key]]], ev)[, .(id, time, Cp, Cm)]
  setkey(s, id, time)
  per <- s[, .(trough = min(Cp), troughm = Cm[which.min(Cp)], troughm_min = min(Cm),
               Cmax = max(Cp), Tmax = time[which.max(Cp)] - 29 * tau,
               AUC24 = sum(diff(time) * (head(Cp, -1) + tail(Cp, -1)) / 2) / 1000,
               AUCm  = sum(diff(time) * (head(Cm, -1) + tail(Cm, -1)) / 2) / 1000,
               ratio_3to8h = mean(Cm[time - 29 * tau >= 3 & time - 29 * tau <= 8] /
                                  Cp[time - 29 * tau >= 3 & time - 29 * tau <= 8])),
           by = id]
  per[, `:=`(met_ratio_auc = AUCm / AUC24, trough_ratio = troughm / trough, dose_mg = dose_mg)]
  per <- merge(per, as.data.table(pop[, c("id", "cl", "vc", "vp", "q")]), by = "id")
  per[, t12 := thalf_terminal(cl, vc, vp, q)]
  per[]
}

# Measured-type concentrations: add Schmidli 2005 proportional residual error
# (26%), log-normal form, with a fixed seed.
add_residual <- function(x, seed) { set.seed(seed); x * exp(rnorm(length(x), 0, PK$resid_prop)) }

qualify_model <- function(pop) {
  per <- rbindlist(lapply(SIM$doses_mg, function(d) ss_exposure(pop, d)))
  per[, trough_obs := add_residual(trough, SIM$seed + 31)]
  mcv <- function(x) sprintf("%.0f (%.0f%%)", mean(x), 100 * sd(x) / mean(x))
  q3  <- function(x, f = "%.0f") sprintf(paste0(f, " [", f, "-", f, "]"), median(x), quantile(x, .25), quantile(x, .75))
  col <- function(p) c(q3(p$trough), q3(p$trough_obs), mcv(p$trough), mcv(p$troughm),
                       sprintf("%.2f (SD %.2f)", mean(p$trough_ratio), sd(p$trough_ratio)),
                       q3(p$ratio_3to8h, "%.2f"), q3(p$Cmax), q3(p$Tmax, "%.1f"),
                       q3(p$AUC24, "%.1f"), q3(p$t12, "%.1f"), q3(p$met_ratio_auc, "%.2f"))
  tab <- data.frame(
    Metric = c("Parent trough, median [IQR] (ng/mL)",
               "Parent trough with 26% residual error, median [IQR] (ng/mL)",
               "Parent trough, mean (CV%) (ng/mL)",
               "CGP74588 trough, mean (CV%) (ng/mL)",
               "Metabolite:parent trough ratio, mean (SD)",
               "Metabolite:parent ratio 3-8 h post-dose, median [IQR]",
               "Cmax, median [IQR] (ng/mL)", "Tmax, median [IQR] (h)",
               "AUC0-24, median [IQR] (ug.h/mL)",
               "Parent terminal half-life, median [IQR] (h)",
               "Metabolite:parent AUC ratio, median [IQR]"),
    Simulated_400mg = col(per[dose_mg == 400]),
    Published_400mg = c("IRIS: median 879 [647-1170], n = 351 (Larson 2008); TOPS: median 1200 (Guilhot 2012)",
                        "IRIS: median 879 [647-1170] (measured concentrations)",
                        "IRIS: mean 979 (CV 54%) (Larson 2008)",
                        "IRIS: 242 (CV 44%) (Larson 2008)",
                        "IRIS: 0.27 (SD 0.09) (Larson 2008); TOPS ~0.20-0.25 (Guilhot 2012)",
                        "0.17 (CV 74%), single samples 3-8 h (le Coutre 2004)",
                        "1920-2040 (le Coutre 2004; Egorin 2009)",
                        "2-4 (prescribing information)",
                        "34-40 (Egorin 2009; le Coutre 2004)",
                        "18.2 (le Coutre 2004)",
                        "~0.15 (prescribing information)"),
    Simulated_600mg = col(per[dose_mg == 600]),
    Published_600mg = c("TOPS: median 1935 at an actual 600 mg/day (Guilhot 2012)",
                        "TOPS: median 1935 (measured concentrations)",
                        "-", "-", "TOPS: ~0.20-0.25 at all doses (Guilhot 2012)",
                        "-", "-", "-", "-", "-", "-"),
    check.names = FALSE)
  write.csv(tab, file.path(DIRS$results, "model_qualification.csv"), row.names = FALSE)
  write.csv(per, file.path(DIRS$results, "qualification_per_subject.csv"), row.names = FALSE)
  list(table = tab, per = per)
}

# External consistency checks against independent published values.
external_checks <- function(per, typ_trough) {
  p4 <- per[dose_mg == 400]; p6 <- per[dose_mg == 600]
  w <- c(860, 1500)
  tab <- data.frame(
    Check = c("Median [IQR] parent trough, 400 mg, with residual error (ng/mL)",
              "Median parent trough, 600 mg, with residual error (ng/mL)",
              "Patients below 1000 ng/mL on 400 mg (%)",
              "Patients within 860-1500 ng/mL on 400 mg (%)",
              "Typical-patient parent trough, 400 mg (ng/mL)",
              "Metabolite:parent trough ratio, 400 mg, mean (SD)",
              "Metabolite:parent ratio 3-8 h post-dose, 400 mg, median",
              "Spearman correlation, steady-state trough vs AUC0-24 (400 mg)"),
    Model = c(sprintf("%.0f [%.0f-%.0f]", median(p4$trough_obs), quantile(p4$trough_obs, .25), quantile(p4$trough_obs, .75)),
              sprintf("%.0f", median(p6$trough_obs)),
              sprintf("%.1f", 100 * mean(p4$trough_obs < 1000)),
              sprintf("%.1f", 100 * mean(p4$trough_obs >= w[1] & p4$trough_obs <= w[2])),
              sprintf("%.0f", typ_trough),
              sprintf("%.2f (%.2f)", mean(p4$trough_ratio), sd(p4$trough_ratio)),
              sprintf("%.2f", median(p4$ratio_3to8h)),
              sprintf("%.3f", cor(p4$trough, p4$AUC24, method = "spearman"))),
    Published = c("879 [647-1170] (IRIS, Larson 2008, n = 351)",
                  "1935 (TOPS, Guilhot 2012)",
                  "64.6 (OPTIM-imatinib, Johnson-Ansah 2022; ~2 weeks of treatment)",
                  "41.4 (real-world TDM, Huang 2023, n = 102)",
                  "555-1388; grand median 870 across nine population models (Gotta 2013)",
                  "0.27 (0.09) (IRIS, Larson 2008)",
                  "0.17 (CV 74%) (le Coutre 2004)",
                  "trough representative of overall exposure (IATDMCT, Clarke 2021)"),
    check.names = FALSE)
  write.csv(tab, file.path(DIRS$results, "external_consistency.csv"), row.names = FALSE)
  tab
}

# ----------------------------------------------------------------------------
# Run-in adequacy: trough just before the first perturbed dose, as % of the
# true steady-state trough (taken after a very long run).
# ----------------------------------------------------------------------------
runin_adequacy <- function(pop, ss_days, dose_mg = 400, long_days = 100, label = "base") {
  tau <- SIM$tau
  ev <- shared_events((0:(long_days - 1)) * tau, dose_mg, c(ss_days * tau, long_days * tau) - 1e-6)
  s <- solve_pop("base", pop[, MODEL_PARAMS$base], ev)[, .(id, time, Cp, Cm)]
  a <- s[time < ss_days * tau]; b <- s[time > ss_days * tau]
  m <- merge(a, b, by = "id", suffixes = c("_a", "_b"))
  data.frame(label = label, run_in_days = ss_days,
             parent_min_pct = 100 * min(m$Cp_a / m$Cp_b),
             parent_median_pct = 100 * median(m$Cp_a / m$Cp_b),
             metab_min_pct = 100 * min(m$Cm_a / m$Cm_b),
             metab_p01_pct = 100 * unname(quantile(m$Cm_a / m$Cm_b, 0.01)),
             metab_median_pct = 100 * median(m$Cm_a / m$Cm_b),
             metab_thalf_max = max(pop$met_thalf))
}

# ----------------------------------------------------------------------------
# Calibration-anchor (exposure-headroom) sensitivity. CL/F is scaled so that
# the simulated MEDIAN steady-state trough equals each anchor (root finding),
# then the forgiveness analysis is repeated (k = 0..6).
# ----------------------------------------------------------------------------
median_ss_trough <- function(pop, dose_mg) {
  tau <- SIM$tau; ss <- SIM$ss_days
  ev <- shared_events((0:(ss - 1)) * tau, dose_mg, seq((ss - 1) * tau, ss * tau, by = 0.25))
  s <- solve_pop("base", pop[, MODEL_PARAMS$base], ev)[, .(id, Cp)]
  median(s[, .(m = min(Cp)), by = id]$m)
}
run_calibration_anchors <- function(pop, dose_mg, targets) {
  bind_rows(lapply(targets, function(tg) {
    f <- function(lg) { pp <- pop; pp$cl <- pop$cl * exp(lg); log(median_ss_trough(pp, dose_mg) / tg) }
    lg <- uniroot(f, c(-1.5, 1.5), tol = 1e-3)$root
    pp <- pop; pp$cl <- pop$cl * exp(lg)
    long <- run_forgiveness(pp, dose_mg, K = SIM$max_consec_sens, label = sprintf("anchor %g", tg), verbose = FALSE)
    s <- summarise_forgiveness(long, label = sprintf("anchor_%g", tg))
    s$anchor <- tg; s$cl_scale <- exp(lg)
    s
  }))
}

# ----------------------------------------------------------------------------
# Parameter-uncertainty propagation (sourced RSEs; UNC table in 00_config.R).
# Each replicate draws a typical-value vector, variances, covariance and weight
# exponents, builds a population of SIM$n_unc_subjects patients and runs a
# single-omission experiment (pre-omission trough = steady-state trough).
# ----------------------------------------------------------------------------
draw_unc <- function(stress = FALSE) {
  u <- UNC; rse <- if (stress) u$stress_rse else u$rse
  x <- setNames(numeric(nrow(u)), u$parameter)
  for (i in seq_len(nrow(u))) {
    x[i] <- if (u$distribution[i] == "log-normal") u$estimate[i] * exp(rnorm(1, 0, rse[i])) else
      { v <- -1; while (v <= 0) v <- rnorm(1, u$estimate[i], rse[i] * u$estimate[i]); v }
  }
  x
}
build_population_unc <- function(N, th) {
  wt <- pmin(130, pmax(45, rnorm(N, COV$wt_mean, COV$wt_sd)))
  z1 <- rnorm(N); z2 <- rnorm(N); z3 <- rnorm(N); z4 <- rnorm(N); z5 <- rnorm(N)
  om_cl <- sqrt(th["omega2_CL"]); om_v <- sqrt(th["omega2_V"])
  r <- max(-0.99, min(0.99, th["cov_CL_V"] / (om_cl * om_v)))
  wr <- wt / COV$wt_ref
  cl <- th["CL"] * wr^th["allo_CL"] * exp(om_cl * z1)
  vc <- th["Vc"] * wr^th["allo_V"] * exp(om_v * (r * z1 + sqrt(1 - r^2) * z2))
  vp <- th["Vp"] * wr^th["allo_V"]
  q  <- th["Q"]  * wr^th["allo_CL"]
  ka <- th["ka"] * exp(PK$iiv_ka * z3)
  ratio <- th["ss_ratio"] * exp(MET$iiv_ratio * z4)
  thalf <- th["t_half"] * exp(MET$iiv_thalf * z5)
  kme <- log(2) / thalf
  data.frame(id = 1:N, ka, cl, vc, q, vp, kf = ratio * kme, kme, rpf = MET$rpf_base)
}
one_omission <- function(pp, dose_mg) {
  tau <- SIM$tau; ss <- SIM$ss_days
  days <- 1:(ss + 3); amt <- rep(dose_mg, length(days)); amt[ss + 1] <- 0
  ev <- shared_events((days[amt > 0] - 1) * tau, dose_mg, seq((ss - 1) * tau, (ss + 2) * tau, by = 0.25))
  s <- solve_pop("base", pp[, MODEL_PARAMS$base], ev)[, .(id, time, Cp, Cm)]
  s[, .(ss_trough = min(Cp[time <= ss * tau]), nadir = min(Cp),
        nadir_m05 = min(Cp + 0.5 * Cm), nadir_m10 = min(Cp + 1.0 * Cm),
        ss_m05 = min((Cp + 0.5 * Cm)[time <= ss * tau])), by = id]
}
run_uncertainty <- function(stress = FALSE, reps = SIM$n_uncertainty, N = SIM$n_unc_subjects) {
  set.seed(SIM$seed + if (stress) 199 else 99)
  thr <- EFF$thr_parent
  out <- vector("list", reps)
  for (i in seq_len(reps)) {
    th <- draw_unc(stress)
    pp <- build_population_unc(N, th)
    rows <- lapply(SIM$doses_mg, function(d) {
      x <- one_omission(pp, d); at <- x$ss_trough >= thr
      data.frame(rep = i, dose_mg = d, median_trough = median(x$ss_trough),
                 pct_at_target = 100 * mean(at),
                 cond_pct_lose_after1_parent = 100 * mean(x$nadir[at] < thr),
                 cond_pct_lose_after1_moiety_r0.5 = 100 * mean(x$nadir_m05[at] < moiety_threshold(0.5, thr, th["ss_ratio"])),
                 cond_pct_lose_after1_moiety_r1.0 = 100 * mean(x$nadir_m10[at] < moiety_threshold(1.0, thr, th["ss_ratio"])),
                 all_pct_fail_after1_parent = 100 * mean(x$nadir < thr))
    })
    out[[i]] <- bind_rows(rows)
    if (i %% 50 == 0) cat(sprintf("    uncertainty (%s) rep %d/%d\n", if (stress) "stress" else "sourced", i, reps))
  }
  rep_tab <- bind_rows(out)
  summ <- rep_tab %>% group_by(dose_mg) %>%
    summarise(across(-rep, list(median = median, p2.5 = ~quantile(.x, .025), p97.5 = ~quantile(.x, .975))),
              .groups = "drop") %>% mutate(run = if (stress) "stress" else "sourced")
  list(replicates = rep_tab %>% mutate(run = if (stress) "stress" else "sourced"), summary = summ)
}
