# =============================================================================
# 09_benchmark.R
# Typical-patient benchmark of the simulated parent and metabolite trajectories
# against the three published parent-metabolite population models (reviewer 2,
# comment 5), at 400 mg once daily, during and after two consecutive omissions.
#
#   This study  2-cpt parent (CL/F 10.5, Vc/F 160, Q/F 50, Vp/F 110, ka 0.5);
#               metabolite t1/2 40 h, R 0.15
#   Menon-Andersen 2009 (CCAP 63:229-38; children and young adults, 70 kg):
#               1-cpt parent CL/F 10.8 L/h, V/F 284 L, zero-order input 1.67 h;
#               2-cpt metabolite CL/F 9.65 L/h, V1/F 11.6 L, Q 2.9 L/h,
#               V2/F 256 L; formation fraction fixed at 0.2 (as published)
#   Delbaldo 2006 (Clin Cancer Res 12:6073-8; adults with GIST, AAG 1.15 g/L):
#               1-cpt parent CL/F 7.97 L/h, V/F 168 L, zero-order input 2.53 h;
#               metabolite CLm/fm 58.6 L/h x 0.55 at steady state, Vm/fm 15.4 L
#   Petain 2008 (Clin Cancer Res 14:7102-9; children and adults; 70 kg, typical
#               AGP and albumin): 1-cpt parent ka 1.03 /h,
#               CL/F 7.29 x (70/54)^0.56 L/h, V/F 202 L; metabolite CLm/fm
#               52.2 L/h x 0.70 at steady state, Vm/fm 32.2 L
# Concentrations are compared as a percentage of each model's own steady-state
# trough, so that the shapes of the declines (which determine forgiveness) are
# compared independently of the differing absolute exposures.
# =============================================================================
suppressMessages({library(rxode2); library(dplyr)})

bm_ours <- rxode2({
  d/dt(depot) = -ka * depot
  d/dt(centr) = ka * depot - (cl / vc) * centr - (q / vc) * centr + (q / vp) * peri
  d/dt(peri)  = (q / vc) * centr - (q / vp) * peri
  Cp = 1000 * centr / vc
  d/dt(Cm) = kf * Cp - kme * Cm
})
# one-compartment parent; apparent metabolite amount Am' = Am/fm in a single
# compartment of volume Vm/fm and clearance CLm/fm (Delbaldo, Petain)
bm_1c1c <- rxode2({
  d/dt(depot) = -ka * depot
  d/dt(centr) = ka * depot - (cl / v) * centr
  Cp = 1000 * centr / v
  d/dt(met) = (cl / v) * centr - (clmf / vmf) * met
  Cm = 1000 * met / vmf
})
bm_ma <- rxode2({
  d/dt(centr) = -(cl / v) * centr
  Cp = 1000 * centr / v
  d/dt(m1c) = fm * (cl / v) * centr - (clm / v1m) * m1c - (qm / v1m) * m1c + (qm / v2m) * m2p
  d/dt(m2p) = (qm / v1m) * m1c - (qm / v2m) * m2p
  Cm = 1000 * m1c / v1m
})

BENCH <- list(
  list(name = "This study (base case)", mod = bm_ours, cmt = "depot", dur = NA,
       p = c(ka = PK$ka, cl = PK$CL, vc = PK$Vc, q = PK$Q, vp = PK$Vp,
             kf = MET$ss_ratio * log(2) / MET$t_half, kme = log(2) / MET$t_half)),
  list(name = "Menon-Andersen 2009", mod = bm_ma, cmt = "centr", dur = 1.67,
       p = c(cl = 10.8, v = 284, fm = 0.2, clm = 9.65, v1m = 11.6, qm = 2.9, v2m = 256)),
  list(name = "Delbaldo 2006", mod = bm_1c1c, cmt = "centr", dur = 2.53,
       p = c(ka = 1, cl = 7.97, v = 168, clmf = 58.6 * 0.55, vmf = 15.4)),
  list(name = "Petain 2008", mod = bm_1c1c, cmt = "depot", dur = NA,
       p = c(ka = 1.03, cl = 7.29 * (70 / 54)^0.56, v = 202, clmf = 52.2 * 0.70, vmf = 32.2))
)

run_benchmark <- function(dose_mg = 400, n_miss = 2) {
  tau <- SIM$tau; ss <- SIM$ss_days; end_day <- ss + n_miss + 6
  prof <- list(); tab <- list()
  for (b in BENCH) {
    sim_k <- function(k) {
      days <- 1:end_day; tk <- rep(TRUE, end_day); if (k > 0) tk[(ss + 1):(ss + k)] <- FALSE
      tt <- (days[tk] - 1) * tau
      ev <- if (is.na(b$dur)) et(time = tt, amt = dose_mg, cmt = b$cmt) else
        et(time = tt, amt = dose_mg, dur = b$dur, cmt = b$cmt)
      ev <- add.sampling(ev, seq((ss - 1) * tau, end_day * tau, by = 0.1))
      as.data.frame(rxSolve(b$mod, b$p, ev))
    }
    s0 <- sim_k(0); s1 <- sim_k(1); s2 <- sim_k(n_miss)
    w0 <- s0[s0$time >= (ss - 1) * tau & s0$time <= ss * tau, ]
    ctp <- min(w0$Cp); ctm <- w0$Cm[which.min(w0$Cp)]
    aucp <- sum(diff(w0$time) * (head(w0$Cp, -1) + tail(w0$Cp, -1)) / 2)
    aucm <- sum(diff(w0$time) * (head(w0$Cm, -1) + tail(w0$Cm, -1)) / 2)
    nad <- function(s, k, y) min(s[[y]][s$time >= (ss - 1) * tau & s$time <= (ss + k + 1) * tau])
    mo <- function(s, rp) s$Cp + rp * s$Cm
    rel_gain <- function(s, k, rp) {
      w <- s$time >= (ss - 1) * tau & s$time <= (ss + k + 1) * tau
      m0 <- min(mo(w0, rp))
      (min(mo(s, rp)[w]) / m0) / (min(s$Cp[w]) / ctp) - 1
    }
    # apparent metabolite half-life over the second drug-free day of the 2-dose omission
    g <- function(h) approx(s2$time, s2$Cm, xout = (ss - 1) * tau + h)$y
    app_m <- 24 * log(2) / log(g(24) / g(48))
    gp <- function(h) approx(s2$time, s2$Cp, xout = (ss - 1) * tau + h)$y
    app_p <- 24 * log(2) / log(gp(24) / gp(48))
    tab[[b$name]] <- data.frame(model = b$name,
      ss_trough_parent = ctp, ss_trough_metab = ctm, trough_ratio = ctm / ctp, auc_ratio = aucm / aucp,
      parent_nadir_after1_pct = 100 * nad(s1, 1, "Cp") / ctp,
      metab_nadir_after1_pct = 100 * nad(s1, 1, "Cm") / ctm,
      parent_nadir_after2_pct = 100 * nad(s2, n_miss, "Cp") / ctp,
      metab_nadir_after2_pct = 100 * nad(s2, n_miss, "Cm") / ctm,
      apparent_t12_parent_day2 = app_p, apparent_t12_metab_day2 = app_m,
      moiety_rel_gain_after1_rpf0.5 = 100 * rel_gain(s1, 1, 0.5),
      moiety_rel_gain_after1_rpf1.0 = 100 * rel_gain(s1, 1, 1.0),
      moiety_rel_gain_after2_rpf1.0 = 100 * rel_gain(s2, n_miss, 1.0))
    prof[[b$name]] <- data.frame(model = b$name, time = s2$time - ss * tau,
                                 parent_pct = 100 * s2$Cp / ctp, metab_pct = 100 * s2$Cm / ctm,
                                 Cp = s2$Cp, Cm = s2$Cm)
  }
  tab <- bind_rows(tab); prof <- bind_rows(prof)
  write.csv(tab, file.path(DIRS$results, "benchmark_published_models.csv"), row.names = FALSE)
  write.csv(prof, file.path(DIRS$results, "benchmark_profiles.csv"), row.names = FALSE)
  list(table = tab, profiles = prof)
}
