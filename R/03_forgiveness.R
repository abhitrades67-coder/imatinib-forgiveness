# =============================================================================
# 03_forgiveness.R
# Individual-level missed-dose experiments. Each patient is brought to steady
# state by a run-in of perfect dosing; then k consecutive scheduled doses are
# omitted (or a dose is taken late, or omitted and the next one doubled) and
# dosing resumes on schedule.
#
# Endpoints (per patient, per k):
#   * nadir of the parent and of the active moiety over the observation window
#     (last dose before the omission -> one interval after dosing resumes,
#     i.e. k + 2 dosing intervals), and time below target in that window
#     (v1.0.0 definitions, retained for continuity);
#   * trajectory metrics over fixed horizons starting at the first omitted dose
#     (48 h, 7, 14 and 28 days): time below target and area below target
#     (linear interpolation), and AUC, from which the EXPOSURE DEFICIT relative
#     to uninterrupted dosing (k = 0) is derived.
# Summaries are reported CONDITIONAL on being at target under uninterrupted
# dosing (primary) and for the whole cohort (secondary).
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(data.table)})

moiety_threshold <- function(rpf, thr = EFF$thr_parent, R = MET$ss_ratio) thr * (1 + rpf * R)

# Dose schedule for a missed-dose experiment (days are 1-indexed; the dose of
# day d is scheduled at (d - 1) * tau; the first perturbed dose is day ss + 1).
missed_schedule <- function(k, dose_mg, ss, end_day, behaviour = "omit", late_h = 0) {
  days <- 1:end_day
  time <- (days - 1) * SIM$tau
  amt  <- rep(dose_mg, end_day)
  if (behaviour == "omit" && k > 0) amt[(ss + 1):(ss + k)] <- 0
  if (behaviour == "late")   time[ss + 1] <- time[ss + 1] + late_h
  if (behaviour == "double") { amt[ss + 1] <- 0; amt[ss + 2] <- 2 * dose_mg }
  keep <- amt > 0
  data.frame(day = days[keep], time = time[keep], amt = amt[keep])
}

# ----------------------------------------------------------------------------
# Core experiment: one dosing perturbation, whole population.
#   model_key  "base" | "mm" | "ma"           (see 01_model.R)
#   ss         run-in length (days)
#   behaviour  "omit" (k doses) | "late" (one dose late_h late) | "double"
#   fmat       optional N x days matrix of relative bioavailability (IOV on F)
#   thr        parent target; R_thr = ratio used to scale the moiety target
#   sweep      also return 7-day time below target across EFF$thr_sweep
# ----------------------------------------------------------------------------
simulate_missed <- function(pop, k, dose_mg, model_key = "base", ss = SIM$ss_days,
                            behaviour = "omit", late_h = 0, fmat = NULL,
                            thr = EFF$thr_parent, R_thr = MET$ss_ratio,
                            rpf_set = MET$rpf_sens, sweep = FALSE) {
  tau <- SIM$tau
  H <- SIM$horizons_h
  end_day <- ss + max(H) / tau + 1
  kk <- if (behaviour == "omit") k else 1          # width of the nadir window
  sched <- missed_schedule(k, dose_mg, ss, end_day, behaviour, late_h)

  w0 <- (ss - 1) * tau                             # last dose before the perturbation
  w1 <- (ss + kk + 1) * tau                        # one interval after dosing resumes
  t_start <- ss * tau                              # first perturbed (omitted) dose
  t_end <- t_start + max(H)
  obs <- sort(unique(c(seq(w0, w1, by = 0.25), seq(w1, t_end, by = 0.5))))

  if (is.null(fmat)) {
    ev <- shared_events(sched$time, sched$amt, obs)
  } else {
    N <- nrow(pop)
    dd <- data.frame(id = rep(pop$id, each = nrow(sched)),
                     time = rep(sched$time, N),
                     amt = as.vector(t(fmat[, sched$day, drop = FALSE])) * rep(sched$amt, N))
    ev <- patient_events(dd, obs, pop$id)
  }
  s <- solve_pop(model_key, pop[, MODEL_PARAMS[[model_key]]], ev)
  s <- s[, .(id, time, Cp, Cm)]
  setkey(s, id, time)

  # --- nadir window (v1.0.0 definitions) ---
  win <- s[time >= w0 & time <= w1]
  out <- win[, .(nadir_parent = min(Cp), nadir_metab = min(Cm), cmax_parent = max(Cp),
                 tb_win_parent = sum(Cp < thr) * 0.25), by = id]
  for (rp in rpf_set) {
    thm <- moiety_threshold(rp, thr, R_thr)
    tmp <- win[, .(n = min(Cp + rp * Cm), t = sum(Cp + rp * Cm < thm) * 0.25), by = id]
    out[tmp, on = "id", `:=`(nm = i.n, tm = i.t)]
    setnames(out, c("nm", "tm"), sprintf(c("nadir_moiety_r%.1f", "tb_win_moiety_r%.1f"), rp))
  }

  # --- trajectory metrics over fixed horizons from the first perturbed dose ---
  tr <- s[time >= t_start & time <= t_end]
  tr[, `:=`(t2 = shift(time, -1), p2 = shift(Cp, -1), m2 = shift(Cm, -1)), by = id]
  tr <- tr[!is.na(t2)]
  tr[, dtt := t2 - time]
  tr[, `:=`(tbp = dtt * frac_below(Cp, p2, thr), abp = dtt * area_below(Cp, p2, thr),
            aucp = dtt * (Cp + p2) / 2, aucm = dtt * (Cm + m2) / 2)]
  for (rp in rpf_set) {
    thm <- moiety_threshold(rp, thr, R_thr)
    y1 <- tr$Cp + rp * tr$Cm; y2 <- tr$p2 + rp * tr$m2
    set(tr, j = sprintf("tbm%.1f", rp), value = tr$dtt * frac_below(y1, y2, thm))
    set(tr, j = sprintf("abm%.1f", rp), value = tr$dtt * area_below(y1, y2, thm))
  }
  mcols <- c("tbp", "abp", "aucp", "aucm", sprintf("tbm%.1f", rpf_set), sprintf("abm%.1f", rpf_set))
  for (h in H) {
    agg <- tr[time < t_start + h, lapply(.SD, sum), by = id, .SDcols = mcols]
    setnames(agg, mcols, paste0(mcols, "_", h))
    out <- out[agg, on = "id"]
  }
  if (sweep) {
    h7 <- tr[time < t_start + 168]
    for (th in EFF$thr_sweep) {
      v <- h7[, .(p = sum(dtt * frac_below(Cp, p2, th)),
                  m05 = sum(dtt * frac_below(Cp + 0.5 * Cm, p2 + 0.5 * m2, moiety_threshold(0.5, th, R_thr))),
                  m10 = sum(dtt * frac_below(Cp + 1.0 * Cm, p2 + 1.0 * m2, moiety_threshold(1.0, th, R_thr)))),
              by = id]
      setnames(v, c("p", "m05", "m10"), sprintf(c("sw_tbp_%d", "sw_tbm0.5_%d", "sw_tbm1.0_%d"), th))
      out <- out[v, on = "id"]
    }
  }
  out[, `:=`(k = k, dose_mg = dose_mg, behaviour = behaviour, late_h = late_h)]
  out[]
}

# ----------------------------------------------------------------------------
# Run k = 0..K for one dose and summarise.
# ----------------------------------------------------------------------------
run_forgiveness <- function(pop, dose_mg, K = SIM$max_consec_miss, label = "base",
                            verbose = TRUE, ...) {
  per_k <- lapply(0:K, function(k) {
    if (verbose) cat(sprintf("    [%s] %d mg: %d consecutive dose(s) omitted\n", label, dose_mg, k))
    simulate_missed(pop, k, dose_mg, ...)
  })
  rbindlist(per_k, fill = TRUE)
}

# Max tolerable consecutive omissions: largest k whose nadir stays >= target
# (0 if the target is not met even after one omission). The nadir falls
# monotonically with k, so this is also the first failing k minus one.
max_tol <- function(k, nadir, th) { ok <- k[nadir >= th]; if (length(ok)) max(ok) else 0L }

# Per-patient tolerance and conditional status from a run_forgiveness() table.
tolerance_table <- function(long, thr = EFF$thr_parent, R_thr = MET$ss_ratio,
                            rpf_set = MET$rpf_sens) {
  base <- long[k == 0, .(id, ss_trough = nadir_parent, ss_trough_metab = nadir_metab)]
  tol <- long[, .(max_tol_parent = max_tol(k, nadir_parent, thr)), by = .(id, dose_mg)]
  for (rp in rpf_set) {
    col <- sprintf("nadir_moiety_r%.1f", rp); thm <- moiety_threshold(rp, thr, R_thr)
    tmp <- long[, .(v = max_tol(k, get(col), thm)), by = id]
    tol[tmp, on = "id", (sprintf("max_tol_moiety_r%.1f", rp)) := i.v]
    tol[, (sprintf("gain_r%.1f", rp)) := get(sprintf("max_tol_moiety_r%.1f", rp)) - max_tol_parent]
  }
  tol <- tol[base, on = "id"]
  tol[, at_target := ss_trough >= thr]
  tol[]
}

# Headline summary of one configuration at one dose (conditional + whole cohort).
summarise_forgiveness <- function(long, thr = EFF$thr_parent, R_thr = MET$ss_ratio,
                                  rpf_set = MET$rpf_sens, label = "base") {
  tol <- tolerance_table(long, thr, R_thr, rpf_set)
  dose <- tol$dose_mg[1]
  at <- tol$at_target
  k1 <- long[k == 1][match(tol$id, id)]
  k0 <- long[k == 0][match(tol$id, id)]
  K <- max(long$k)
  res <- data.frame(label = label, dose_mg = dose, n = nrow(tol),
    pct_at_target = 100 * mean(at),
    median_ss_trough = median(tol$ss_trough),
    # conditional (primary): among patients at target under uninterrupted dosing
    cond_pct_lose_after1_parent = 100 * mean(k1$nadir_parent[at] < thr),
    cond_median_max_tol_parent = median(tol$max_tol_parent[at]),
    cond_q1_max_tol_parent = unname(quantile(tol$max_tol_parent[at], .25)),
    cond_q3_max_tol_parent = unname(quantile(tol$max_tol_parent[at], .75)),
    cond_pct_tol_ge1_parent = 100 * mean(tol$max_tol_parent[at] >= 1),
    cond_pct_tol_ge3_parent = 100 * mean(tol$max_tol_parent[at] >= 3),
    cond_mean_tb_win_after1_parent = mean(k1$tb_win_parent[at]),
    # whole cohort (secondary composite of target attainment and forgiveness)
    all_pct_fail_after1_parent = 100 * mean(k1$nadir_parent < thr),
    all_median_max_tol_parent = median(tol$max_tol_parent),
    all_pct_tol_ge3_parent = 100 * mean(tol$max_tol_parent >= 3),
    all_mean_tb_win_after1_parent = mean(k1$tb_win_parent),
    n_at_cap = sum(tol$max_tol_parent >= K),
    stringsAsFactors = FALSE)
  for (rp in rpf_set) {
    s <- sprintf("%.1f", rp); thm <- moiety_threshold(rp, thr, R_thr)
    nm <- k1[[paste0("nadir_moiety_r", s)]]
    res[[paste0("cond_pct_lose_after1_moiety_r", s)]] <- 100 * mean(nm[at] < thm)
    res[[paste0("cond_pct_gain_ge1_r", s)]] <- 100 * mean(tol[[paste0("gain_r", s)]][at] >= 1)
    res[[paste0("cond_mean_gain_r", s)]] <- mean(tol[[paste0("gain_r", s)]][at])
    res[[paste0("all_pct_gain_ge1_r", s)]] <- 100 * mean(tol[[paste0("gain_r", s)]] >= 1)
    res[[paste0("cond_median_max_tol_moiety_r", s)]] <- median(tol[[paste0("max_tol_moiety_r", s)]][at])
    res[[paste0("cond_mean_tb_win_after1_moiety_r", s)]] <- mean(k1[[paste0("tb_win_moiety_r", s)]][at])
    res[[paste0("all_pct_fail_after1_moiety_r", s)]] <- 100 * mean(nm < thm)
    res[[paste0("pct_at_target_moiety_r", s)]] <- 100 * mean(k0[[paste0("nadir_moiety_r", s)]] >= thm)
    res[[paste0("all_mean_tb_win_after1_moiety_r", s)]] <- mean(k1[[paste0("tb_win_moiety_r", s)]])
  }
  # trajectory metrics after 1 missed dose, conditional: excess time and area
  # below target over 7 days, and exposure deficit over each horizon
  for (h in SIM$horizons_h) {
    dp <- 1 - k1[[paste0("aucp_", h)]] / k0[[paste0("aucp_", h)]]
    dm <- 1 - k1[[paste0("aucm_", h)]] / k0[[paste0("aucm_", h)]]
    res[[paste0("cond_median_deficit_parent_", h)]] <- 100 * median(dp[at])
    res[[paste0("cond_median_deficit_metab_", h)]] <- 100 * median(dm[at])
    for (rp in rpf_set) {
      s <- sprintf("%.1f", rp)
      mo1 <- k1[[paste0("aucp_", h)]] + rp * k1[[paste0("aucm_", h)]]
      mo0 <- k0[[paste0("aucp_", h)]] + rp * k0[[paste0("aucm_", h)]]
      res[[paste0("cond_median_deficit_moiety_r", s, "_", h)]] <- 100 * median((1 - mo1 / mo0)[at])
    }
    xp <- k1[[paste0("tbp_", h)]] - k0[[paste0("tbp_", h)]]
    res[[paste0("cond_mean_excess_tb_parent_", h)]] <- mean(xp[at])
    res[[paste0("cond_mean_excess_ab_parent_", h)]] <- mean((k1[[paste0("abp_", h)]] - k0[[paste0("abp_", h)]])[at])
    for (rp in rpf_set) {
      s <- sprintf("%.1f", rp)
      xm <- k1[[paste0("tbm", s, "_", h)]] - k0[[paste0("tbm", s, "_", h)]]
      res[[paste0("cond_mean_excess_tb_moiety_r", s, "_", h)]] <- mean(xm[at])
      res[[paste0("cond_mean_excess_ab_moiety_r", s, "_", h)]] <-
        mean((k1[[paste0("abm", s, "_", h)]] - k0[[paste0("abm", s, "_", h)]])[at])
    }
  }
  res
}

# Conditional trajectory metrics for every k (for figures and Table 3).
trajectory_by_k <- function(long, thr = EFF$thr_parent, rpf_set = MET$rpf_sens) {
  k0 <- long[k == 0]
  at_ids <- k0$id[k0$nadir_parent >= thr]
  out <- list()
  for (kk in sort(unique(long$k))) {
    x <- long[k == kk][match(k0$id, id)]
    at <- x$id %in% at_ids
    row <- data.frame(dose_mg = x$dose_mg[1], k = kk,
      cond_pct_maintain_parent = 100 * mean(x$nadir_parent[at] >= thr),
      all_pct_maintain_parent = 100 * mean(x$nadir_parent >= thr),
      cond_mean_nadir_parent = mean(x$nadir_parent[at]),
      all_mean_nadir_parent = mean(x$nadir_parent), all_sd_nadir_parent = sd(x$nadir_parent))
    for (rp in rpf_set) {
      s <- sprintf("%.1f", rp); thm <- moiety_threshold(rp, thr)
      nm <- x[[paste0("nadir_moiety_r", s)]]
      row[[paste0("cond_pct_maintain_moiety_r", s)]] <- 100 * mean(nm[at] >= thm)
      row[[paste0("all_mean_nadir_moiety_r", s)]] <- mean(nm)
      row[[paste0("all_sd_nadir_moiety_r", s)]] <- sd(nm)
    }
    for (h in SIM$horizons_h) {
      row[[paste0("cond_mean_excess_tb_parent_", h)]] <- mean((x[[paste0("tbp_", h)]] - k0[[paste0("tbp_", h)]])[at])
      row[[paste0("cond_mean_excess_ab_parent_", h)]] <- mean((x[[paste0("abp_", h)]] - k0[[paste0("abp_", h)]])[at])
      row[[paste0("cond_median_deficit_parent_", h)]] <- 100 * median((1 - x[[paste0("aucp_", h)]] / k0[[paste0("aucp_", h)]])[at])
      row[[paste0("cond_median_deficit_metab_", h)]] <- 100 * median((1 - x[[paste0("aucm_", h)]] / k0[[paste0("aucm_", h)]])[at])
      for (rp in rpf_set) {
        s <- sprintf("%.1f", rp)
        row[[paste0("cond_mean_excess_tb_moiety_r", s, "_", h)]] <-
          mean((x[[paste0("tbm", s, "_", h)]] - k0[[paste0("tbm", s, "_", h)]])[at])
        mo1 <- x[[paste0("aucp_", h)]] + rp * x[[paste0("aucm_", h)]]
        mo0 <- k0[[paste0("aucp_", h)]] + rp * k0[[paste0("aucm_", h)]]
        row[[paste0("cond_median_deficit_moiety_r", s, "_", h)]] <- 100 * median((1 - mo1 / mo0)[at])
      }
    }
    out[[length(out) + 1]] <- row
  }
  bind_rows(out)
}

# Threshold sweep from stored nadirs (and the in-loop 7-day time-below sweep).
threshold_sweep <- function(long, rpf_set = c(0.5, 1.0)) {
  bind_rows(lapply(EFF$thr_sweep, function(th) {
    tol <- tolerance_table(long, thr = th, rpf_set = rpf_set)
    at <- tol$at_target
    k1 <- long[k == 1][match(tol$id, id)]; k0 <- long[k == 0][match(tol$id, id)]
    r <- data.frame(dose_mg = tol$dose_mg[1], threshold = th,
                    pct_at_target = 100 * mean(at),
                    cond_pct_lose_after1_parent = if (any(at)) 100 * mean(k1$nadir_parent[at] < th) else NA,
                    all_pct_fail_after1_parent = 100 * mean(k1$nadir_parent < th),
                    cond_median_max_tol_parent = if (any(at)) median(tol$max_tol_parent[at]) else NA)
    for (rp in rpf_set) {
      s <- sprintf("%.1f", rp)
      r[[paste0("cond_pct_gain_ge1_r", s)]] <- if (any(at)) 100 * mean(tol[[paste0("gain_r", s)]][at] >= 1) else NA
    }
    cp <- sprintf("sw_tbp_%d", th)
    if (cp %in% names(long)) {
      r$cond_mean_excess_tb7d_parent <- mean((k1[[cp]] - k0[[cp]])[at])
      for (s in c("0.5", "1.0")) {
        cm <- sprintf("sw_tbm%s_%d", s, th)
        r[[paste0("cond_mean_excess_tb7d_moiety_r", s)]] <- mean((k1[[cm]] - k0[[cm]])[at])
      }
    }
    r
  }))
}

# ----------------------------------------------------------------------------
# Michaelis-Menten formation: per-patient normalisation of kfmax so that the
# steady-state 400 mg metabolite:parent AUC ratio equals the patient's R.
# ----------------------------------------------------------------------------
mm_normalise <- function(pop, fdrive, km = MET$km_cyp2c8, dose_ref = 400) {
  ev <- shared_events((0:(SIM$ss_days - 1)) * SIM$tau, dose_ref,
                      seq((SIM$ss_days - 1) * SIM$tau, SIM$ss_days * SIM$tau, by = 0.1))
  s <- solve_pop("base", pop[, MODEL_PARAMS$base], ev)[, .(id, time, Cp)]
  setkey(s, id, time)
  g <- s[, {
    sat <- Cp / (1 + fdrive * Cp / MET$mw_imatinib / km)
    .(g = sum(diff(time) * (head(sat, -1) + tail(sat, -1)) / 2) /
          sum(diff(time) * (head(Cp, -1) + tail(Cp, -1)) / 2))
  }, by = id]
  p <- pop
  p$kfmax <- pop$kf / g$g[match(pop$id, g$id)]
  p$fdrive <- fdrive; p$mw <- MET$mw_imatinib; p$km <- km
  p
}

# Menon-Andersen metabolite parameters (allometric to 70 kg as published).
ma_population <- function(pop) {
  p <- pop; w <- pop$wt / 70
  p$clm <- MET$ma$clm * w^0.75; p$v1m <- MET$ma$v1m * w
  p$qm  <- MET$ma$qm  * w^0.75; p$v2m <- MET$ma$v2m * w
  p$fm  <- pop$met_ratio * p$clm / pop$cl     # AUC ratio = fm * CL / CLm = R_i
  p
}

# Inter-occasion variability in relative bioavailability: mean-preserving
# log-normal per dose, F = exp(eta - omega^2/2), omega^2 = log(1 + CV^2).
iov_matrix <- function(N, ndays, cv, seed) {
  set.seed(seed)
  om <- sqrt(log(1 + cv^2))
  matrix(exp(rnorm(N * ndays, 0, om) - om^2 / 2), N, ndays)
}
