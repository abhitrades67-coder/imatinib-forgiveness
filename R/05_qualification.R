# =============================================================================
# 05_qualification.R
# (A5) Model qualification against published exposure metrics + VPC-style data.
# (B7) Parameter-uncertainty propagation of the headline result.
# =============================================================================
suppressMessages({library(rxode2); library(dplyr)})

# ----------------------------------------------------------------------------
# A5. Qualify the model: simulate perfect 400 mg dosing to steady state and
# compare simulated exposure (median [IQR]) against published values.
# ----------------------------------------------------------------------------
# Analytic terminal half-life of the two-compartment parent disposition (per
# subject), so the model's terminal t1/2 can be qualified against the measured
# value (le Coutre 2004, 18.2 h) rather than left unverified.
thalf_terminal <- function(cl, vc, vp, q) {
  k10 <- cl / vc; k12 <- q / vc; k21 <- q / vp
  b <- k10 + k12 + k21; cc <- k10 * k21
  beta <- (b - sqrt(b^2 - 4 * cc)) / 2           # smaller (terminal) eigenvalue
  log(2) / beta
}

# Per-subject steady-state exposure summary at a given dose.
ss_exposure <- function(pop, dose_mg) {
  ev <- et(time = (0:29) * SIM$tau, amt = dose_mg, cmt = "depot")
  ev <- add.sampling(ev, seq(0, 30 * SIM$tau, by = 0.25))
  sim <- rxSolve(imatinib_model, params = pop[, PARAM_COLS], events = ev,
                 returnType = "data.frame", cores = 4)
  names(sim)[names(sim) == "sim.id"] <- "id"
  ss <- sim[sim$time >= 29 * SIM$tau & sim$time <= 30 * SIM$tau, ]
  per <- ss %>% group_by(id) %>%
    summarise(
      trough  = min(Cp),
      troughm = min(Cm),
      Cmax    = max(Cp),
      Tmax    = time[which.max(Cp)] - 29 * SIM$tau,
      AUC24   = sum(diff(time) * (head(Cp, -1) + tail(Cp, -1)) / 2) / 1000,
      AUCm    = sum(diff(time) * (head(Cm, -1) + tail(Cm, -1)) / 2) / 1000,
      .groups = "drop") %>%
    mutate(met_ratio = AUCm / AUC24, dose_mg = dose_mg)
  disp <- pop[, c("id", "cl", "vc", "vp", "q")]
  per <- merge(per, disp, by = "id")
  per$t12 <- thalf_terminal(per$cl, per$vc, per$vp, per$q)
  per
}

# Qualify the model against published 400 mg exposure; additionally report the
# model's dose-proportional prediction at the other simulated dose(s). Published
# reference ranges exist only at 400 mg, so 600 mg is labelled as a prediction.
qualify_model <- function(pop, ref_dose = 400, pred_doses = setdiff(SIM$doses_mg, 400)) {
  per_all <- bind_rows(lapply(SIM$doses_mg, function(d) ss_exposure(pop, d)))
  per <- per_all[per_all$dose_mg == ref_dose, ]

  mcv <- function(x) sprintf("%.0f (%.0f%%)", mean(x), 100*sd(x)/mean(x))   # mean (CV%)
  q3  <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, .25), quantile(x, .75))
  q3r <- function(x) sprintf("%.2f [%.2f-%.2f]", median(x), quantile(x, .25), quantile(x, .75))
  fmt_col <- function(p) c(mcv(p$trough), mcv(p$troughm), q3(p$Cmax),
      sprintf("%.1f [%.1f-%.1f]", median(p$Tmax), quantile(p$Tmax,.25), quantile(p$Tmax,.75)),
      sprintf("%.1f [%.1f-%.1f]", median(p$AUC24), quantile(p$AUC24,.25), quantile(p$AUC24,.75)),
      sprintf("%.1f [%.1f-%.1f]", median(p$t12), quantile(p$t12,.25), quantile(p$t12,.75)),
      q3r(p$met_ratio))
  tab <- data.frame(
    Metric = c("SS trough imatinib, mean (CV%) (ng/mL)",
               "SS trough CGP74588, mean (CV%) (ng/mL)",
               "Cmax (ng/mL), median [IQR]", "Tmax (h), median [IQR]",
               "AUC0-24 parent (ug.h/mL), median [IQR]",
               "Terminal half-life parent (h), median [IQR]",
               "Metabolite:parent AUC ratio, median [IQR]"),
    Simulated = fmt_col(per),
    Published_reference = c("published population means 979-1452 (Larson 2008 IRIS; Picard 2007)",
                            "242 (44%) (Larson 2008)",
                            "1920-2040 (le Coutre 2004; Egorin 2009)",
                            "2-4 (Gleevec label)",
                            "34-40 (Egorin 2009; le Coutre 2004)",
                            "18.2 (le Coutre 2004; ~18 h consensus)",
                            "~0.15 (Gleevec label; Larson 2008 trough ratio ~0.25)"),
    check.names = FALSE)
  for (d in pred_doses)
    tab[[sprintf("Predicted_%dmg", d)]] <- fmt_col(per_all[per_all$dose_mg == d, ])

  write.csv(tab, file.path(DIRS$results, "model_qualification.csv"), row.names = FALSE)
  write.csv(per_all, file.path(DIRS$results, "qualification_per_subject.csv"), row.names = FALSE)
  list(table = tab, per = per_all)
}

# ----------------------------------------------------------------------------
# B7. Parameter-uncertainty propagation. Resample the population-typical values
# (CL, Vc, ka, metabolite ratio, metabolite t1/2) from their reported relative
# standard errors, rebuild a population, and recompute the headline metrics:
#   - median steady-state parent trough at 400 mg
#   - % of patients above the 1000 ng/mL target
#   - % of patients tolerating 0 consecutive missed doses (parent criterion)
# Efficient: only k = 0 (SS) and k = 1 simulations are needed for these metrics.
# ----------------------------------------------------------------------------
build_population_typ <- function(N, typ) {
  wt  <- pmin(130, pmax(45, rnorm(N, COV$wt_mean, COV$wt_sd)))
  sex_cl <- ifelse(rbinom(N, 1, COV$female_frac) == 1, COV$female_cl, 1)
  cl <- typ["CL"] * (wt/COV$wt_ref)^COV$allo_CL * sex_cl * exp(rnorm(N, 0, PK$iiv_CL))
  vc <- typ["Vc"] * (wt/COV$wt_ref)^COV$allo_V  * exp(rnorm(N, 0, PK$iiv_Vc))
  ka <- typ["ka"] * exp(rnorm(N, 0, PK$iiv_ka))
  vp <- PK$Vp * (wt/COV$wt_ref)^COV$allo_V
  q  <- PK$Q  * (wt/COV$wt_ref)^COV$allo_CL
  thalf <- typ["t_half"] * exp(rnorm(N, 0, MET$iiv_thalf))
  ratio <- typ["ss_ratio"] * exp(rnorm(N, 0, MET$iiv_ratio))
  kme <- log(2)/thalf; kf <- ratio * kme
  data.frame(id = 1:N, ka, cl, vc, q, vp, kf, kme, rpf = MET$rpf_base)
}

# ----------------------------------------------------------------------------
# Calibration-target sensitivity. The population central tendency of imatinib
# trough is itself uncertain (IRIS population mean ~980 ng/mL, CV 54%, Larson
# 2008; responder-subgroup means ~1450). We scale CL/F so the simulated MEDIAN
# trough hits a range of targets and recompute the headline (% tolerating 0
# missed doses), showing how sensitive the conclusion is to the calibration
# anchor rather than only to parameter noise around a single anchor.
# ----------------------------------------------------------------------------
run_calibration_sensitivity <- function(pop, dose_mg = 400, targets = CAL$targets) {
  thr <- EFF$thr_parent
  base_med <- median(simulate_missed(pop, 0, dose_mg) %>%
                       filter(rpf == MET$rpf_base) %>% pull(nadir_parent))
  out <- lapply(targets, function(tg) {
    pp <- pop; pp$cl <- pop$cl * (base_med / tg)         # scale CL to hit target median trough
    s0 <- simulate_missed(pp, 0, dose_mg) %>% filter(rpf == MET$rpf_base)
    s1 <- simulate_missed(pp, 1, dose_mg) %>% filter(rpf == MET$rpf_base)
    data.frame(target_trough = tg, dose_mg = dose_mg,
               median_trough = median(s0$nadir_parent),
               pct_above = 100 * mean(s0$nadir_parent >= thr),
               pct_tol0  = 100 * mean(s1$nadir_parent < thr))
  })
  res <- bind_rows(out)
  write.csv(res, file.path(DIRS$results, "calibration_sensitivity.csv"), row.names = FALSE)
  res
}

run_uncertainty <- function(dose_mg = 400, N = 300, reps = SIM$n_uncertainty) {
  set.seed(SIM$seed + 99)
  base <- c(CL = PK$CL, Vc = PK$Vc, ka = PK$ka, ss_ratio = MET$ss_ratio, t_half = MET$t_half)
  rse  <- SIM$unc_rse[names(base)]
  thr  <- EFF$thr_parent
  res <- vector("list", reps)
  for (i in seq_len(reps)) {
    typ <- base * exp(rnorm(length(base), 0, rse))
    pp  <- build_population_typ(N, typ)
    s0  <- simulate_missed(pp, 0, dose_mg)   # SS
    s1  <- simulate_missed(pp, 1, dose_mg)   # one missed dose
    s0b <- s0[s0$rpf == MET$rpf_base, ]; s1b <- s1[s1$rpf == MET$rpf_base, ]
    res[[i]] <- data.frame(
      rep = i,
      median_trough = median(s0b$nadir_parent),
      pct_above = 100 * mean(s0b$nadir_parent >= thr),
      pct_tol0  = 100 * mean(s1b$nadir_parent < thr)
    )
    if (i %% 50 == 0) cat(sprintf("    uncertainty rep %d/%d\n", i, reps))
  }
  out <- bind_rows(res)
  write.csv(out, file.path(DIRS$results, "uncertainty_replicates.csv"), row.names = FALSE)
  summ <- data.frame(
    Metric = c("Median SS trough (ng/mL)", "% above 1000 ng/mL target",
               "% tolerating 0 missed doses"),
    Median = c(median(out$median_trough), median(out$pct_above), median(out$pct_tol0)),
    P2.5   = c(quantile(out$median_trough, .025), quantile(out$pct_above, .025), quantile(out$pct_tol0, .025)),
    P97.5  = c(quantile(out$median_trough, .975), quantile(out$pct_above, .975), quantile(out$pct_tol0, .975))
  )
  write.csv(summ, file.path(DIRS$results, "uncertainty_summary.csv"), row.names = FALSE)
  list(replicates = out, summary = summ)
}
