# =============================================================================
# 03_forgiveness.R
# Individual-level forgiveness: bring each patient to steady state, then withhold
# k = 0..K consecutive doses. For each k, capture the nadir of PARENT and of the
# ACTIVE MOIETY during the drug-free interval, and the time spent below threshold.
#
# Headline metric: maximum tolerable consecutive missed doses
#   - parent-only criterion   : nadir(Cp)     >= thr_parent
#   - metabolite-aware (moiety): nadir(moiety) >= thr_moiety (consistently scaled)
# The DIFFERENCE between the two is the forgiveness bought by the long-lived
# active metabolite -- evaluated without double-counting, ACROSS the contested
# relative-potency range (MET$rpf_sens).
# =============================================================================
suppressMessages({library(rxode2); library(dplyr)})

# Consistently-scaled active-moiety threshold for a given relative potency factor.
moiety_threshold <- function(rpf) EFF$thr_parent * (1 + rpf * MET$ss_ratio)

# Simulate the whole population for k consecutive missed doses at a dose level.
# Returns per-subject nadirs: parent once, moiety for each rpf in MET$rpf_sens.
simulate_missed <- function(pop, k, dose_mg) {
  ss <- SIM$ss_days
  recover <- 12
  total_days <- ss + k + recover

  day_idx <- 1:total_days
  taken <- rep(TRUE, total_days)
  if (k > 0) taken[(ss + 1):(ss + k)] <- FALSE
  dose_times <- (day_idx[taken] - 1) * SIM$tau

  ev <- et(time = dose_times, amt = dose_mg, cmt = "depot")
  t0 <- (ss - 1) * SIM$tau
  t_end <- total_days * SIM$tau
  samp <- sort(unique(c(seq(0, t_end, by = 1), seq(t0, t_end, by = 0.25))))
  ev <- add.sampling(ev, samp)

  sim <- rxSolve(imatinib_model, params = pop[, PARAM_COLS], events = ev,
                 returnType = "data.frame", cores = 4)
  names(sim)[names(sim) == "sim.id"] <- "id"

  win_start <- ss * SIM$tau
  win_end   <- (ss + k) * SIM$tau + SIM$tau
  w <- sim[sim$time >= win_start - SIM$tau & sim$time <= win_end,
           c("id", "time", "Cp", "Cm")]
  dt <- if (nrow(w) > 1) median(diff(sort(unique(w$time)))) else 0.25

  # parent nadir / time-below (rpf-independent)
  base <- w %>% group_by(id) %>%
    summarise(nadir_parent = min(Cp),
              nadir_metab  = min(Cm),
              time_below_parent = sum(Cp < EFF$thr_parent) * dt,
              .groups = "drop")

  # moiety nadir / time-below for each candidate relative potency factor
  out <- lapply(MET$rpf_sens, function(rpf) {
    thr_m <- moiety_threshold(rpf)
    m <- w %>% mutate(moiety = Cp + rpf * Cm) %>% group_by(id) %>%
      summarise(nadir_moiety = min(moiety),
                time_below_moiety = sum(moiety < thr_m) * dt,
                .groups = "drop")
    base %>% left_join(m, by = "id") %>%
      mutate(rpf = rpf, thr_moiety = thr_m, k = k, dose_mg = dose_mg)
  })
  bind_rows(out)
}

run_individual_forgiveness <- function(pop, dose_mg) {
  K <- SIM$max_consec_miss
  per_k <- lapply(0:K, function(k) {
    cat(sprintf("  %d mg: withholding %d consecutive dose(s)...\n", dose_mg, k))
    simulate_missed(pop, k, dose_mg)
  })
  long <- bind_rows(per_k)

  # max tolerable consecutive misses per subject x rpf
  tol <- long %>% group_by(id, rpf) %>%
    summarise(
      max_tol_parent = { ok <- k[nadir_parent >= EFF$thr_parent]; if (length(ok)) max(ok) else 0L },
      max_tol_moiety = { ok <- k[nadir_moiety >= thr_moiety[1]];   if (length(ok)) max(ok) else 0L },
      .groups = "drop") %>%
    mutate(dose_mg = dose_mg, metab_gain = max_tol_moiety - max_tol_parent)

  # threshold sensitivity (A1): max tolerable missed doses on the PARENT criterion
  # evaluated at each efficacy threshold in EFF$thr_sens
  thr_tol <- lapply(EFF$thr_sens, function(th) {
    long %>% filter(rpf == MET$rpf_base) %>% group_by(id) %>%
      summarise(max_tol = { ok <- k[nadir_parent >= th]; if (length(ok)) max(ok) else 0L },
                .groups = "drop") %>%
      summarise(threshold = th, dose_mg = dose_mg,
                mean_max_tol = mean(max_tol), median_max_tol = median(max_tol),
                pct_tol0 = 100 * mean(max_tol == 0),
                pct_tol_ge3 = 100 * mean(max_tol >= 3))
  }) %>% bind_rows()

  # cohort summary by number of missed doses x rpf
  cohort <- long %>% group_by(k, rpf) %>%
    summarise(
      nadir_parent_mean = mean(nadir_parent), nadir_parent_sd = sd(nadir_parent),
      nadir_moiety_mean = mean(nadir_moiety), nadir_moiety_sd = sd(nadir_moiety),
      time_below_parent_mean = mean(time_below_parent),
      time_below_moiety_mean = mean(time_below_moiety),
      pct_above_parent = 100 * mean(nadir_parent >= EFF$thr_parent),
      pct_above_moiety = 100 * mean(nadir_moiety >= thr_moiety[1]),
      .groups = "drop") %>%
    mutate(dose_mg = dose_mg)

  list(long = long, tolerance = tol, cohort = cohort, threshold = thr_tol)
}

# ----------------------------------------------------------------------------
# A6 (EXPLORATORY): individual forgiveness under CYP3A4 DDI states.
# Applies a CL/F multiplier and a metabolite-formation-ratio multiplier to the
# population, then recomputes parent-criterion max tolerable missed doses.
# Clearly labelled exploratory; run at 400 mg only.
# ----------------------------------------------------------------------------
run_ddi_sensitivity <- function(pop, dose_mg = 400) {
  K <- SIM$max_consec_miss
  out <- lapply(DDI$states, function(st) {
    pp <- pop
    pp$cl <- pop$cl * DDI$cl_mult[[st]]
    pp$kf <- pop$kf * DDI$ratio_mult[[st]]      # scale formation (metabolite) with CYP3A4
    per_k <- lapply(0:K, function(k) simulate_missed(pp, k, dose_mg))
    long <- bind_rows(per_k)
    long %>% filter(rpf == MET$rpf_base) %>% group_by(id) %>%
      summarise(
        max_tol_parent = { ok <- k[nadir_parent >= EFF$thr_parent]; if (length(ok)) max(ok) else 0L },
        ss_trough = nadir_parent[k == 0],
        .groups = "drop") %>%
      summarise(state = st, dose_mg = dose_mg,
                median_ss_trough = median(ss_trough),
                mean_max_tol = mean(max_tol_parent),
                pct_tol0 = 100 * mean(max_tol_parent == 0),
                pct_above_target = 100 * mean(ss_trough >= EFF$thr_parent))
  })
  bind_rows(out)
}
