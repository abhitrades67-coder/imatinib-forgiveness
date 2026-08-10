# =============================================================================
# 02_population.R
# Generate the virtual CML population and simulate adherence scenarios at
# 400 and 600 mg QD. Records parent, metabolite and active-moiety troughs.
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(tidyr)})

# ----------------------------------------------------------------------------
# 1. Virtual population with evidence-based IIV and allometric covariates
# ----------------------------------------------------------------------------
build_population <- function() {
  set.seed(SIM$seed)
  N <- SIM$n_subjects
  age <- pmin(85, pmax(18, rnorm(N, COV$age_mean, COV$age_sd)))  # descriptor only
  wt  <- pmin(130, pmax(45, rnorm(N, COV$wt_mean, COV$wt_sd)))
  sex <- rbinom(N, 1, COV$female_frac)               # 1 = female

  eta_cl <- rnorm(N, 0, PK$iiv_CL)
  eta_vc <- rnorm(N, 0, PK$iiv_Vc)
  eta_ka <- rnorm(N, 0, PK$iiv_ka)

  sex_cl <- ifelse(sex == 1, COV$female_cl, 1)        # modest lower CL/F in females
  cl <- PK$CL * (wt / COV$wt_ref)^COV$allo_CL * sex_cl * exp(eta_cl)
  vc <- PK$Vc * (wt / COV$wt_ref)^COV$allo_V  * exp(eta_vc)
  vp <- PK$Vp * (wt / COV$wt_ref)^COV$allo_V
  q  <- PK$Q  * (wt / COV$wt_ref)^COV$allo_CL
  ka <- PK$ka * exp(eta_ka)

  # metabolite IIV: each patient has their own formation ratio and half-life
  ratio_i <- MET$ss_ratio  * exp(rnorm(N, 0, MET$iiv_ratio))
  thalf_i <- MET$t_half    * exp(rnorm(N, 0, MET$iiv_thalf))
  kme <- log(2) / thalf_i                            # metabolite elim rate (median t1/2 = 40 h)
  kf  <- ratio_i * kme                               # formation rate (median SS ratio = R)
  rpf <- rep(MET$rpf_base, N)

  data.frame(id = 1:N, age, wt, sex,
             ka, cl, vc, q, vp, kf, kme, rpf,
             met_ratio = ratio_i, met_thalf = thalf_i)
}

# ----------------------------------------------------------------------------
# 2. Adherence-scenario dosing schedules (nominal dose indices 1..ndose)
#    Returns a data.frame of administered doses: time (h), amt (mg)
# ----------------------------------------------------------------------------
make_schedule <- function(scenario, dose_mg, ndays = SIM$followup_days) {
  tau <- SIM$tau
  idx <- 1:ndays                       # dose index = day number
  nominal_time <- (idx - 1) * tau      # 0, 24, 48, ...
  taken <- rep(TRUE, ndays)
  delay <- rep(0, ndays)

  if (scenario == "perfect") {
    # all taken on time
  } else if (scenario == "moderate") {        # 80%: drop every 5th dose
    taken[idx %% 5 == 0] <- FALSE
  } else if (scenario == "poor") {            # 40%: 2 taken, 3 missed, repeating
    taken[((idx - 1) %% 5) >= 2] <- FALSE
  } else if (scenario == "weekend") {         # drop Sat+Sun (day1 = Monday)
    dow <- ((idx - 1) %% 7) + 1               # 1=Mon ... 6=Sat,7=Sun
    taken[dow %in% c(6, 7)] <- FALSE
  } else if (scenario == "delayed") {         # 8 random doses delayed 6 h
    set.seed(SIM$seed + 7)
    dl <- sample(idx, 8)
    delay[dl] <- 6
  } else stop("unknown scenario: ", scenario)

  data.frame(time = nominal_time[taken] + delay[taken],
             amt  = dose_mg,
             nominal_day = idx[taken])
}

# ----------------------------------------------------------------------------
# 3. Simulate one scenario for the whole population
# ----------------------------------------------------------------------------
PARAM_COLS <- c("id", "ka", "cl", "vc", "q", "vp", "kf", "kme", "rpf")

simulate_scenario <- function(pop, scenario, dose_mg) {
  sched <- make_schedule(scenario, dose_mg)
  ev <- et(time = sched$time, amt = sched$amt, cmt = "depot")
  # fine grid over the whole horizon so the per-day trough (daily minimum) and the
  # SS-day AUC are both well resolved regardless of dose timing.
  samp <- sort(unique(c(seq(0, 24 * SIM$followup_days, by = 0.25))))
  ev <- add.sampling(ev, samp)

  sim <- rxSolve(imatinib_model, params = pop[, PARAM_COLS], events = ev,
                 returnType = "data.frame", cores = 4)
  names(sim)[names(sim) == "sim.id"] <- "id"

  # --- TROUGH = daily minimum concentration within each nominal 24 h day, over the
  #     evaluation window (days 15-30). This pre-dose-equivalent metric is robust to
  #     shifted dose timing (eliminates the delayed-dosing sampling artifact). ---
  eval_lo <- 14 * SIM$tau
  w <- sim[sim$time >= eval_lo & sim$time < SIM$followup_days * SIM$tau, ]
  w$day <- floor(w$time / SIM$tau) + 1
  tro <- w %>% group_by(id, day) %>%
    summarise(Cp = min(Cp), Cm = min(Cm), moiety = min(moiety), .groups = "drop") %>%
    mutate(scenario = scenario, dose_mg = dose_mg)

  # --- AUC0-24 and Cmax on representative SS day (day 25) per subject ---
  d <- sim[sim$time >= 24 * 24 & sim$time <= 24 * 25, ]
  auc <- d %>% group_by(id) %>%
    summarise(AUC_24h = sum(diff(time) * (head(Cp, -1) + tail(Cp, -1)) / 2) / 1000,
              Cmax = max(Cp), .groups = "drop")

  list(troughs = tro, auc = auc)
}

# ----------------------------------------------------------------------------
# 4. Run all scenarios x doses; summarise forgiveness metrics
# ----------------------------------------------------------------------------
SCENARIOS <- c("perfect", "moderate", "poor", "weekend", "delayed")
SCEN_LABEL <- c(perfect = "Perfect (100%)", moderate = "Moderate (80%)",
                poor = "Poor (40%)", weekend = "Weekend holiday",
                delayed = "Delayed (6 h)")

run_all_scenarios <- function(pop) {
  all_tro <- list(); all_metrics <- list()
  for (dose in SIM$doses_mg) {
    for (scn in SCENARIOS) {
      cat(sprintf("  simulating %d mg / %s ...\n", dose, scn))
      res <- simulate_scenario(pop, scn, dose)
      tro <- merge(res$troughs, res$auc, by = "id")
      all_tro[[paste(dose, scn)]] <- tro

      # per-subject mean daily trough over evaluation window, then summarise cohort
      subj <- res$troughs %>% group_by(id) %>%
        summarise(Cp = mean(Cp), Cm = mean(Cm), moiety = mean(moiety), .groups = "drop") %>%
        left_join(res$auc, by = "id")

      thr_p <- EFF$thr_parent
      thr_m <- EFF$thr_parent * (1 + MET$rpf_base * MET$ss_ratio)   # consistent moiety threshold

      # target attainment across the threshold-sensitivity set (A1)
      att <- as.list(setNames(
        sapply(EFF$thr_sens, function(th) 100 * mean(subj$Cp >= th)),
        paste0("Pct_Above_", EFF$thr_sens)))

      all_metrics[[paste(dose, scn)]] <- data.frame(
        dose_mg = dose, scenario = scn, label = SCEN_LABEL[[scn]],
        Mean_Parent_Trough   = mean(subj$Cp),
        SD_Parent_Trough     = sd(subj$Cp),
        Median_Parent_Trough = median(subj$Cp),
        # 95% population PREDICTION interval (percentiles), not a CI
        PI_Lower = quantile(subj$Cp, 0.025), PI_Upper = quantile(subj$Cp, 0.975),
        Mean_Metab_Trough    = mean(subj$Cm),
        Median_Metab_Trough  = median(subj$Cm),
        Mean_Moiety_Trough   = mean(subj$moiety),
        Median_Moiety_Trough = median(subj$moiety),
        Mean_AUC_24h = mean(subj$AUC_24h), Median_AUC_24h = median(subj$AUC_24h),
        Mean_Cmax = mean(subj$Cmax), Median_Cmax = median(subj$Cmax),
        Pct_Above_Parent = 100 * mean(subj$Cp >= thr_p),
        Pct_Above_Moiety = 100 * mean(subj$moiety >= thr_m),
        as.data.frame(att)
      )
    }
  }
  list(troughs = bind_rows(all_tro),
       metrics = bind_rows(all_metrics))
}
