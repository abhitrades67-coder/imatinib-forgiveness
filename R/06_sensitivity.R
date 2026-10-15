# =============================================================================
# 06_sensitivity.R
# Structural, kinetic, variability, covariate, formation-kinetics, omission-
# behaviour and drug-interaction sensitivity analyses (reviewer 1 comments 1, 5,
# 8; reviewer 2 comments 4, 5, 7, 8, 10). Each configuration repeats the
# missed-dose experiment (k = 0..6) in the same 1000 virtual patients and
# reports the same conditional and whole-cohort summaries as the base case.
# =============================================================================
suppressMessages({library(dplyr); library(data.table)})

SENS_CONFIGS <- list(
  # --- metabolite structure and kinetics ---
  list(name = "thalf20",  group = "Metabolite structure and kinetics", label = "Metabolite half-life 20 h (Gschwind 2005)", t_half = 20),
  list(name = "thalf90",  group = "Metabolite structure and kinetics", label = "Metabolite half-life 90 h (le Coutre 2004); 35-day run-in", t_half = 90, ss = 35),
  list(name = "ratio025", group = "Metabolite structure and kinetics", label = "Metabolite:parent ratio R = 0.25", ss_ratio = 0.25),
  list(name = "formlim",  group = "Metabolite structure and kinetics", label = "Formation-rate-limited metabolite (Delbaldo 2006; Petain 2008)", t_half = MET$thalf_formlim, ss_ratio = MET$ratio_formlim),
  list(name = "ma015",    group = "Metabolite structure and kinetics", label = "Two-compartment metabolite (Menon-Andersen 2009), R = 0.15", model = "ma"),
  list(name = "ma025",    group = "Metabolite structure and kinetics", label = "Two-compartment metabolite (Menon-Andersen 2009), R = 0.25", model = "ma", ss_ratio = 0.25),
  list(name = "mm_unbound", group = "Metabolite structure and kinetics", label = "Michaelis-Menten formation, K_m 1.4 uM, unbound concentration", model = "mm", fdrive = MET$fu_plasma),
  list(name = "mm_total", group = "Metabolite structure and kinetics", label = "Michaelis-Menten formation, K_m 1.4 uM, total concentration (extreme bound)", model = "mm", fdrive = 1),
  # --- variability and covariates ---
  list(name = "submitted", group = "Variability and covariate model", label = "Submitted (v1.0.0) model: allometry 0.75/1, female CL x0.92, no CL-Vc correlation", variant = "submitted"),
  list(name = "base_qvp", group = "Variability and covariate model", label = "Added IIV on Q (41.8%) and Vp (29.7%) (Golabchifar 2016)", variant = "base_qvp"),
  list(name = "iov28",    group = "Variability and covariate model", label = "Added inter-occasion variability in F, CV 28%", iov = IOV$cv_F[1]),
  list(name = "iov44",    group = "Variability and covariate model", label = "Added inter-occasion variability in F, CV 44%", iov = IOV$cv_F[2]),
  list(name = "metiiv0",  group = "Variability and covariate model", label = "Metabolite IIV removed", met_iiv_mult = 0),
  list(name = "metiiv2",  group = "Variability and covariate model", label = "Metabolite IIV doubled", met_iiv_mult = 2),
  # --- exploratory CYP3A4 interaction (parent clearance only) ---
  list(name = "ddi_inducer",   group = "Exploratory CYP3A4 interaction", label = "Strong inducer: CL/F x3.85 (Bolton 2004)", cl_mult = DDI$cl_mult[["inducer"]]),
  list(name = "ddi_inhibitor", group = "Exploratory CYP3A4 interaction", label = "Strong inhibitor: CL/F x0.71 (Dutreix 2004; upper bound at steady state)", cl_mult = DDI$cl_mult[["inhibitor"]])
)

run_sensitivity_config <- function(cfg, K = SIM$max_consec_sens) {
  g <- function(x, d) if (is.null(cfg[[x]])) d else cfg[[x]]
  model_key <- g("model", "base"); ss <- g("ss", SIM$ss_days)
  R <- g("ss_ratio", MET$ss_ratio)
  pop <- build_population(variant = g("variant", "base"), t_half = g("t_half", MET$t_half),
                          ss_ratio = R, met_iiv_mult = g("met_iiv_mult", 1))
  if (!is.null(cfg$cl_mult)) pop$cl <- pop$cl * cfg$cl_mult
  if (model_key == "ma") pop <- ma_population(pop)
  if (model_key == "mm") pop <- mm_normalise(pop, cfg$fdrive)
  fmat <- if (!is.null(cfg$iov)) iov_matrix(nrow(pop), 80, cfg$iov, SIM$seed + 41) else NULL
  rows <- list(); longs <- list()
  for (d in SIM$doses_mg) {
    # moiety target scaled with the configuration's own steady-state ratio; for
    # saturable formation the realised ratio is dose dependent, so it is taken
    # from the simulated steady state at that dose
    R_thr <- if (model_key == "mm") median(ss_exposure(pop, d, "mm")$met_ratio_auc) else R
    long <- run_forgiveness(pop, d, K = K, label = cfg$name, verbose = FALSE,
                            model_key = model_key, ss = ss, fmat = fmat, R_thr = R_thr)
    s <- summarise_forgiveness(long, R_thr = R_thr, label = cfg$name)
    s$group <- cfg$group; s$description <- cfg$label; s$R_thr <- R_thr
    rows[[as.character(d)]] <- s; longs[[as.character(d)]] <- long
    cat(sprintf("  [%s] %d mg: at target %.1f%%, conditional loss after 1 %.1f%%, gain>=1 (RPF 0.5/1.0) %.1f%% / %.1f%%\n",
                cfg$name, d, s$pct_at_target, s$cond_pct_lose_after1_parent,
                s$cond_pct_gain_ge1_r0.5, s$cond_pct_gain_ge1_r1.0))
  }
  saveRDS(rbindlist(longs), file.path(DIRS$results, "sensitivity", paste0(cfg$name, ".rds")))
  bind_rows(rows)
}

run_all_sensitivity <- function() {
  dir.create(file.path(DIRS$results, "sensitivity"), showWarnings = FALSE)
  bind_rows(lapply(SENS_CONFIGS, run_sensitivity_config))
}

# ----------------------------------------------------------------------------
# Omission behaviour (reviewer 1 comments 4-5): the first perturbed dose is
# taken 6, 12 or 18 h late (next dose on schedule), or omitted with the next
# dose doubled. Compared with complete omission (k = 1) and uninterrupted
# dosing (k = 0) in the same patients.
# ----------------------------------------------------------------------------
run_omission_behaviour <- function(pop, base_long) {
  thr <- EFF$thr_parent
  out <- list()
  for (d in SIM$doses_mg) {
    k0 <- base_long[dose_mg == d & k == 0]; setkey(k0, id)
    at_ids <- k0$id[k0$nadir_parent >= thr]
    cases <- c(list(list(b = "omit", h = 0, lab = "Omitted (base case)")),
               lapply(OMIT$late_h, function(h) list(b = "late", h = h, lab = sprintf("Taken %d h late", h))),
               list(list(b = "double", h = 0, lab = "Omitted, next dose doubled")))
    for (cs in cases) {
      x <- if (cs$b == "omit") base_long[dose_mg == d & k == 1] else
        simulate_missed(pop, 1, d, behaviour = cs$b, late_h = cs$h)
      x <- x[match(k0$id, id)]; at <- x$id %in% at_ids
      r <- data.frame(dose_mg = d, behaviour = cs$lab,
                      cond_pct_lose_target_parent = 100 * mean(x$nadir_parent[at] < thr),
                      cond_median_nadir_pct_of_ss = 100 * median((x$nadir_parent / k0$nadir_parent)[at]),
                      cond_mean_excess_tb_parent_168 = mean((x$tbp_168 - k0$tbp_168)[at]),
                      cond_median_deficit_parent_168 = 100 * median((1 - x$aucp_168 / k0$aucp_168)[at]),
                      median_peak_pct_of_ss_cmax = 100 * median(x$cmax_parent / k0$cmax_parent),
                      p95_peak_ng_ml = unname(quantile(x$cmax_parent, .95)))
      for (rp in c(0.5, 1.0)) {
        s <- sprintf("%.1f", rp)
        r[[paste0("cond_pct_lose_target_moiety_r", s)]] <-
          100 * mean(x[[paste0("nadir_moiety_r", s)]][at] < moiety_threshold(rp))
      }
      out[[length(out) + 1]] <- r
    }
  }
  bind_rows(out)
}

# Typical-patient concentration profiles for the omission behaviours (figure).
omission_profiles <- function(dose_mg = 400) {
  typ <- data.frame(id = 1, ka = PK$ka, cl = PK$CL, vc = PK$Vc, q = PK$Q, vp = PK$Vp,
                    kf = MET$ss_ratio * log(2) / MET$t_half, kme = log(2) / MET$t_half, rpf = MET$rpf_base)
  ss <- SIM$ss_days; tau <- SIM$tau; end_day <- ss + 6
  obs <- seq((ss - 1) * tau, (ss + 4) * tau, by = 0.1)
  cases <- list(c("uninterrupted", "omit", 0, 0), c("omitted", "omit", 1, 0),
                c("6 h late", "late", 1, 6), c("12 h late", "late", 1, 12), c("18 h late", "late", 1, 18),
                c("omitted, next doubled", "double", 1, 0))
  bind_rows(lapply(cases, function(cs) {
    sc <- missed_schedule(as.integer(cs[3]), dose_mg, ss, end_day, cs[2], as.numeric(cs[4]))
    s <- solve_pop("base", typ, shared_events(sc$time, sc$amt, obs))
    data.frame(case = cs[1], time = s$time - ss * tau, Cp = s$Cp, Cm = s$Cm)
  }))
}
