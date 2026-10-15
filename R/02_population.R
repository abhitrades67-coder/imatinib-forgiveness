# =============================================================================
# 02_population.R
# Virtual CML population (with variability-model variants) and adherence
# patterns (deterministic illustrative patterns and stochastic, per-patient
# patterns), evaluated AFTER a steady-state run-in.
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(tidyr); library(data.table)})

N_CORES <- max(1L, min(8L, parallel::detectCores() - 1L))

# ----------------------------------------------------------------------------
# 1. Virtual population
# ----------------------------------------------------------------------------
# Random numbers are drawn in the same order as v1.0.0 (age, weight, sex,
# eta_CL, eta_Vc, eta_ka, metabolite ratio, metabolite half-life), so the
# "submitted" variant reproduces the v1.0.0 population exactly; the draws for
# IIV on Q and Vp are appended after them.
#
# variant:
#   "base"       Schmidli 2005 weight exponents (0.301 / 0.405), no sex effect,
#                CL-Vc correlation 0.707                    (revised base case)
#   "submitted"  allometric 0.75 / 1.0, female CL x0.92, no correlation (v1.0.0)
#   "base_qvp"   base + IIV on Q (41.8%) and Vp (29.7%), correlation 0.921
#                (Golabchifar 2016)
build_population <- function(variant = "base", N = SIM$n_subjects, seed = SIM$seed,
                             t_half = MET$t_half, ss_ratio = MET$ss_ratio,
                             met_iiv_mult = 1, rpf = MET$rpf_base) {
  set.seed(seed)
  age <- pmin(85, pmax(18, rnorm(N, COV$age_mean, COV$age_sd)))   # descriptor only
  wt  <- pmin(130, pmax(45, rnorm(N, COV$wt_mean, COV$wt_sd)))
  sex <- rbinom(N, 1, COV$female_frac)                              # 1 = female
  z_cl <- rnorm(N); z_vc <- rnorm(N); z_ka <- rnorm(N)
  z_ratio <- rnorm(N); z_thalf <- rnorm(N)
  z_q <- rnorm(N); z_vp <- rnorm(N)

  if (variant == "submitted") {
    aCL <- COV_SUBMITTED$allo_CL; aV <- COV_SUBMITTED$allo_V
    sex_cl <- ifelse(sex == 1, COV_SUBMITTED$female_cl, 1); r <- COV_SUBMITTED$corr_CL_Vc
  } else {
    aCL <- COV$allo_CL; aV <- COV$allo_V
    sex_cl <- ifelse(sex == 1, COV$female_cl, 1); r <- PK$corr_CL_Vc
  }
  eta_cl <- PK$iiv_CL * z_cl
  eta_vc <- PK$iiv_Vc * (r * z_cl + sqrt(1 - r^2) * z_vc)
  if (variant == "base_qvp") {
    rq <- PK$corr_Q_Vp
    eta_q  <- PK$iiv_Q  * z_q
    eta_vp <- PK$iiv_Vp * (rq * z_q + sqrt(1 - rq^2) * z_vp)
  } else { eta_q <- 0; eta_vp <- 0 }

  wr <- wt / COV$wt_ref
  cl <- PK$CL * wr^aCL * sex_cl * exp(eta_cl)
  vc <- PK$Vc * wr^aV * exp(eta_vc)
  vp <- PK$Vp * wr^aV * exp(eta_vp)
  q  <- PK$Q  * wr^aCL * exp(eta_q)
  ka <- PK$ka * exp(PK$iiv_ka * z_ka)

  ratio_i <- ss_ratio * exp(met_iiv_mult * MET$iiv_ratio * z_ratio)
  thalf_i <- t_half   * exp(met_iiv_mult * MET$iiv_thalf * z_thalf)
  kme <- log(2) / thalf_i
  kf  <- ratio_i * kme

  data.frame(id = 1:N, age, wt, sex, ka, cl, vc, q, vp, kf, kme, rpf = rpf,
             met_ratio = ratio_i, met_thalf = thalf_i)
}

# Analytic terminal half-life of the two-compartment parent disposition.
thalf_terminal <- function(cl, vc, vp, q) {
  k10 <- cl / vc; k12 <- q / vc; k21 <- q / vp
  b <- k10 + k12 + k21
  log(2) / ((b - sqrt(b^2 - 4 * k10 * k21)) / 2)
}

# ----------------------------------------------------------------------------
# 2. Adherence patterns over the pattern window (days 1..ndays AFTER run-in)
# ----------------------------------------------------------------------------
# Deterministic, ILLUSTRATIVE patterns, identical for every patient:
#   moderate  every fifth dose missed                   (80%)
#   poor      two taken, three missed, repeating       (40%)
#   weekend   Saturday and Sunday missed; the window starts on a Monday (22/30 = 73%)
#   delayed   eight of thirty doses taken 6 h late (same days for every patient)
make_pattern <- function(scenario, ndays = SIM$pattern_days) {
  idx <- 1:ndays
  taken <- rep(TRUE, ndays); delay <- rep(0, ndays)
  if (scenario == "perfect") {
  } else if (scenario == "moderate") {
    taken[idx %% 5 == 0] <- FALSE
  } else if (scenario == "poor") {
    taken[((idx - 1) %% 5) >= 2] <- FALSE
  } else if (scenario == "weekend") {
    dow <- ((idx - 1) %% 7) + 1                # 1 = Monday ... 7 = Sunday
    taken[dow %in% c(6, 7)] <- FALSE
  } else if (scenario == "delayed") {
    set.seed(SIM$seed + 7)
    delay[sample(idx, 8)] <- 6
  } else stop("unknown scenario: ", scenario)
  list(taken = taken, delay = delay)
}

# Stochastic, per-patient patterns: two-state Markov chain on the omission
# indicator with stationary omission probability pi = 1 - adherence and lag-1
# autocorrelation rho (rho = 0: independent omissions). Returns an N x ndays
# logical matrix of doses TAKEN.
stochastic_taken <- function(N, ndays, adherence, rho, seed) {
  set.seed(seed)
  p <- 1 - adherence
  p11 <- p + rho * (1 - p)        # P(miss today | missed yesterday)
  p01 <- p * (1 - rho)            # P(miss today | taken yesterday)
  miss <- matrix(FALSE, N, ndays)
  miss[, 1] <- runif(N) < p
  for (d in 2:ndays) {
    u <- runif(N)
    miss[, d] <- ifelse(miss[, d - 1], u < p11, u < p01)
  }
  !miss
}

SCENARIOS <- c("perfect", "moderate", "poor", "weekend", "delayed")
SCEN_LABEL <- c(perfect = "Perfect (100%)", moderate = "Every fifth dose missed (80%)",
                poor = "Two taken, three missed (40%)",
                weekend = "Weekend omissions (73%)", delayed = "Eight doses 6 h late (100%)")

# ----------------------------------------------------------------------------
# 3. Event-table helpers
# ----------------------------------------------------------------------------
# Shared schedule (same for every patient): et() event table.
shared_events <- function(dose_time, dose_amt, obs) {
  ev <- et(time = dose_time, amt = dose_amt, cmt = "depot")
  add.sampling(ev, obs)
}
# Per-patient schedule: data frame with one block of dose records per id.
patient_events <- function(dose_df, obs, ids) {
  # dose_df: data.frame(id, time, amt)
  d <- data.frame(id = dose_df$id, time = dose_df$time, amt = dose_df$amt,
                  evid = 1L, cmt = "depot")
  o <- data.frame(id = rep(ids, each = length(obs)), time = rep(obs, times = length(ids)),
                  amt = 0, evid = 0L, cmt = "depot")
  ev <- rbind(d, o)
  ev[order(ev$id, ev$time, -ev$evid), ]
}
solve_pop <- function(model_key, params, ev) {
  s <- rxSolve(MODELS[[model_key]], params = params, events = ev,
               returnType = "data.frame", cores = N_CORES)
  if ("sim.id" %in% names(s)) names(s)[names(s) == "sim.id"] <- "id"
  s <- as.data.table(s)
  # a single-subject solve returns no id column
  if (!"id" %in% names(s)) s[, id := if (is.data.frame(params) && "id" %in% names(params)) params$id[1] else 1L]
  s
}

# ----------------------------------------------------------------------------
# 4. Adherence-pattern simulation (run-in, then a 30-day pattern window)
# ----------------------------------------------------------------------------
# taken_mat: optional N x ndays matrix (stochastic patterns); otherwise the
# deterministic pattern 'scenario' is used for every patient.
# Per patient, over the 30 nominal days of the pattern window:
#   mean daily trough   mean of the 30 daily minimum concentrations
#   pct_days_below      % of days whose minimum is below the target
#   pct_time_below      % of time below the target (linear interpolation)
#   auc                 AUC over the window (for the exposure deficit)
simulate_pattern <- function(pop, dose_mg, scenario = "perfect", taken_mat = NULL,
                             model_key = "base", rpf_set = MET$rpf_sens) {
  ss <- SIM$ss_days; tau <- SIM$tau; nd <- SIM$pattern_days
  t0 <- ss * tau; t1 <- (ss + nd) * tau
  obs <- seq(t0, t1, by = 0.25)
  if (is.null(taken_mat)) {
    pat <- make_pattern(scenario, nd)
    dt_run <- (0:(ss - 1)) * tau
    dt_pat <- t0 + (which(pat$taken) - 1) * tau + pat$delay[pat$taken]
    ev <- shared_events(c(dt_run, dt_pat), dose_mg, obs)
    adh <- rep(mean(pat$taken), nrow(pop))
  } else {
    N <- nrow(pop)
    run <- data.frame(id = rep(pop$id, each = ss), time = rep((0:(ss - 1)) * tau, N), amt = dose_mg)
    w <- which(taken_mat, arr.ind = TRUE)
    pat <- data.frame(id = pop$id[w[, 1]], time = t0 + (w[, 2] - 1) * tau, amt = dose_mg)
    ev <- patient_events(rbind(run, pat), obs, pop$id)
    adh <- rowMeans(taken_mat)
  }
  s <- solve_pop(model_key, pop[, MODEL_PARAMS[[model_key]]], ev)
  s <- s[, .(id, time, Cp, Cm)]
  thr <- EFF$thr_parent
  s[, day := pmin(floor((time - t0) / tau) + 1, nd)]
  daily <- s[, .(Cp = min(Cp), Cm = min(Cm)), by = .(id, day)]
  out <- daily[, .(mean_trough_parent = mean(Cp), mean_trough_metab = mean(Cm),
                   pct_days_below_parent = 100 * mean(Cp < thr)), by = id]
  for (rp in rpf_set) {
    thm <- moiety_threshold(rp)
    tmp <- daily[, .(v = 100 * mean(Cp + rp * Cm < thm)), by = id]
    out[tmp, on = "id", (sprintf("pct_days_below_moiety_r%.1f", rp)) := i.v]
  }
  # time below target and AUC with linear interpolation between samples
  setkey(s, id, time)
  s[, `:=`(t2 = shift(time, -1), p2 = shift(Cp, -1), m2 = shift(Cm, -1)), by = id]
  s <- s[!is.na(t2)]
  s[, dtt := t2 - time]
  tb <- s[, .(pct_time_below_parent = 100 * sum(dtt * frac_below(Cp, p2, thr)) / (t1 - t0),
              auc_parent = sum(dtt * (Cp + p2) / 2), auc_metab = sum(dtt * (Cm + m2) / 2)), by = id]
  out <- out[tb, on = "id"]
  for (rp in rpf_set) {
    thm <- moiety_threshold(rp)
    tmp <- s[, .(v = 100 * sum(dtt * frac_below(Cp + rp * Cm, p2 + rp * m2, thm)) / (t1 - t0)), by = id]
    out[tmp, on = "id", (sprintf("pct_time_below_moiety_r%.1f", rp)) := i.v]
  }
  out[, `:=`(adherence = adh[match(id, pop$id)], dose_mg = dose_mg, scenario = scenario)]
  out[]
}

# Fraction of the interval [t_i, t_i+1] spent below threshold c, assuming the
# concentration varies linearly between the two samples (vectorised).
frac_below <- function(y1, y2, c) {
  s1 <- y1 - c; s2 <- y2 - c
  f <- numeric(length(s1))
  both <- s1 < 0 & s2 < 0
  f[both] <- 1
  up <- s1 < 0 & s2 >= 0                       # crosses upward
  f[up] <- (-s1[up]) / (s2[up] - s1[up])
  dn <- s1 >= 0 & s2 < 0                       # crosses downward
  f[dn] <- (-s2[dn]) / (s1[dn] - s2[dn])
  f
}
# Area below threshold c over the interval (same linear interpolation), per unit time.
area_below <- function(y1, y2, c) {
  s1 <- y1 - c; s2 <- y2 - c
  a <- numeric(length(s1))
  both <- s1 < 0 & s2 < 0
  a[both] <- (-(s1[both] + s2[both])) / 2
  up <- s1 < 0 & s2 >= 0
  fu <- (-s1[up]) / (s2[up] - s1[up]); a[up] <- 0.5 * fu * (-s1[up])
  dn <- s1 >= 0 & s2 < 0
  fd <- (-s2[dn]) / (s1[dn] - s2[dn]); a[dn] <- 0.5 * fd * (-s2[dn])
  a
}

# Run deterministic and stochastic patterns at both doses.
run_adherence_patterns <- function(pop) {
  res <- list()
  for (dose in SIM$doses_mg) {
    for (scn in SCENARIOS) {
      cat(sprintf("  pattern %d mg / %s ...\n", dose, scn))
      res[[paste(dose, scn)]] <- simulate_pattern(pop, dose, scn)
    }
    i <- 0
    for (a in STOCH$adherence) for (rho in STOCH$rho) {
      i <- i + 1
      lab <- sprintf("stoch_%02d_%s", round(100 * a), if (rho == 0) "independent" else "clustered")
      cat(sprintf("  pattern %d mg / %s ...\n", dose, lab))
      tm <- stochastic_taken(nrow(pop), SIM$pattern_days, a, rho, SIM$seed + 11 + i)
      res[[paste(dose, lab)]] <- simulate_pattern(pop, dose, lab, taken_mat = tm)
    }
  }
  rbindlist(res, fill = TRUE)
}
