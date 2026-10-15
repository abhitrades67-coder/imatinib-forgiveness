# =============================================================================
# 00_config.R  (revision 1.1.0, October 2026 -- response to reviewers)
# Imatinib forgiveness analysis in CML -- parent + active metabolite (CGP74588)
# Evidence-based parameter configuration (every value cited below)
# =============================================================================
#
# DESIGN PRINCIPLE
# ----------------
# This is a model-based simulation study of PHARMACOKINETIC forgiveness: the
# dosing history is known by design and the question is how long exposure stays
# at or above a population-derived trough target after omitted doses. It is not
# an adherence-detection (inverse) problem. Imatinib has a major active
# metabolite, N-desmethyl-imatinib (CGP74588). Its terminal half-life is
# reported between ~20 h (single dose; Gschwind 2005) and ~85-95 h (after
# stopping long-term dosing; le Coutre 2004); the prescribing information gives
# ~40 h, which is the base case here and 20 h / 90 h are sensitivity analyses.
#
# No individual concentrations are fitted. The framework is ASSEMBLED from
# published estimates; four disposition parameters are CALIBRATED; several
# variability terms are ASSUMED. The status of every parameter (adopted /
# calibrated / assumed) and its source are written to
# results/parameter_provenance.csv, which is Table 1 of the revised article.
#
# CHANGES FROM v1.0.0 (made in response to peer review)
#   * Covariate model now uses the empirical estimates of the source population
#     analysis (Schmidli 2005, Table 4): weight exponents 0.301 (CL) and 0.405
#     (V); no sex effect (not retained in that final model). The conventional
#     allometric model of v1.0.0 is retained as a sensitivity analysis.
#   * The CL-V covariance of Schmidli 2005 (0.071; r = 0.707) is now applied.
#   * Uncertainty RSEs are now each traced to a source (UNC table below).
#   * Michaelis-Menten metabolite formation, published metabolite structures,
#     metabolite half-life, variability, omission-behaviour, stochastic-adherence
#     and threshold-sweep sensitivity analyses added (R/06_sensitivity.R).
#   * The CYP3A4 DDI analysis is now parent-only (the v1.0.0 metabolite-formation
#     multipliers did not reproduce the cited studies; the reported v1.0.0 DDI
#     results were parent-only and are unaffected).
# =============================================================================

# ----------------------------------------------------------------------------
# PARENT IMATINIB PHARMACOKINETICS (apparent, oral; steady state in CML)
# ----------------------------------------------------------------------------
# CL/F: clearance falls over the first month of dosing. Schmidli 2005 (phase III
#   IRIS popPK, n = 371 CP-CML; Br J Clin Pharmacol 60:35-44,
#   doi:10.1111/j.1365-2125.2005.02372.x), Table 4:
#     CL = (13.8 - 3.81*OCC) * (WT/80)^0.301 * (Hb/13)^0.897 * (WBC/16)^-0.105
#     V  = (252  - 7.82*OCC) * (WT/80)^0.405 * (Hb/13)^0.676 * (WBC/16)^-0.070
#     omega^2 CL 0.102 (SE 0.0142), omega^2 V 0.099 (SE 0.0157),
#     cov(CL,V) 0.071 (SE 0.0104) -> r = 0.707; residual 26% proportional.
#   Day-29 (steady-state) typical CL = 9.99 L/h at 80 kg (9.6 L/h at 70 kg).
#   Golabchifar 2014 (SS CML, n = 61): CL/F 10.8 L/h (CCAP 74:85-93).
#   --> base-case steady-state CL/F = 10.5 L/h at 70 kg (within 9.6-10.8).
# Two-compartment disposition: Golabchifar 2016 (healthy volunteers, NLMEM;
#   Eur J Drug Metab Pharmacokinet 41:527-39): ka 0.94/h (IIV 31.2%), Vc/F 139 L
#   (21.5%), Vp/F 130 L (29.7%), Q/F 29.6 L/h (41.8%); RSE CL/F 7%, ka 10%,
#   Vc/F 4%, Q/F 11%, Vp/F 9%; eta(Vp)-eta(Q) correlation 0.921.
#   Vc, Vp, Q and ka were CALIBRATED (not adopted) so that the combined set
#   reproduces the observed 400 mg steady-state Cmax, AUC and trough and the
#   terminal half-life (le Coutre 2004: 18.2 h).
PK <- list(
  CL    = 10.5,   # L/h  adopted: Schmidli 2005 (9.6 at 70 kg) / Golabchifar 2014 (10.8)
  Vc    = 160.0,  # L    CALIBRATED (cf. Golabchifar 2016, 139)
  Vp    = 110.0,  # L    CALIBRATED (cf. Golabchifar 2016, 130)
  Q     = 50.0,   # L/h  CALIBRATED (cf. Golabchifar 2016, 29.6)
  ka    = 0.5,    # 1/h  CALIBRATED; below published (Schmidli 1.1; Golabchifar 0.94)
  # Inter-individual variability: SD of the log-normal random effect (omega)
  iiv_CL = 0.32,  # Schmidli 2005 omega^2 0.102 -> omega 0.319
  iiv_Vc = 0.31,  # Schmidli 2005 omega^2 0.099 -> omega 0.315 (applied to Vc)
  iiv_ka = 0.40,  # ASSUMED; Golabchifar 2016: 31.2% on ka, 47-48% on its zero-order input terms
  corr_CL_Vc = 0.707, # Schmidli 2005 cov 0.071 / sqrt(0.102 * 0.099)
  # Sensitivity only (variant "base_qvp"): IIV on Q and Vp (Golabchifar 2016)
  iiv_Q  = sqrt(log(1 + 0.418^2)),  # CV 41.8%
  iiv_Vp = sqrt(log(1 + 0.297^2)),  # CV 29.7%
  corr_Q_Vp = 0.921,                # Golabchifar 2016 eta(Vp)-eta(Q)
  # Residual error: NOT part of the forgiveness simulations (true concentrations);
  # added only when simulated concentrations are compared with MEASURED ones.
  resid_prop = 0.26   # Schmidli 2005, 26% proportional
)

# ----------------------------------------------------------------------------
# ACTIVE METABOLITE  N-desmethyl-imatinib (CGP74588)
# ----------------------------------------------------------------------------
# Formation: N-demethylation by CYP2C8 and CYP3A4 (Nebot 2010; Filppula 2013).
# Half-life: ~40 h (Gleevec US prescribing information, rev. 7/2026);
#   20.6 +/- 1.7 h after a single 239 mg dose in 4 healthy men (Gschwind 2005);
#   84.4 and 94.5 h in two patients after stopping imatinib (le Coutre 2004).
# Ratio R (steady-state metabolite:parent): AUC ratio ~15% (prescribing
#   information); IRIS trough ratio 0.27 +/- 0.09 (Larson 2008; CV 33%, which
#   supports the 30% IIV assumed on R); TOPS trough ratio ~20-25% at all doses
#   (Guilhot 2012).
# Relative potency (RPF) is CONTESTED (label "comparable"; Manley 2013
#   "significantly less active") -> sensitivity axis 0.1 / 0.5 / 1.0.
# Identifiability-safe parameterisation: dCm/dt = kf*Cp - kme*Cm with
#   kme = ln2/t_half and kf = R*kme.
MET <- list(
  t_half = 40.0,        # h  prescribing information (base case)
  thalf_sens = c(20, 90), # h  Gschwind 2005 (single dose); le Coutre 2004 (after stopping)
  ss_ratio = 0.15,      # -  prescribing-information AUC ratio (base case)
  ratio_sens = 0.25,    # -  trough-based ratio (IRIS 0.27; TOPS 0.20-0.25; Delbaldo 2006 AUC 0.25)
  rpf_base = 0.5,       # -  explicitly arbitrary midpoint of the contested range
  rpf_sens = c(0.1, 0.5, 1.0),
  iiv_ratio  = 0.30,    # ASSUMED; consistent with the IRIS trough-ratio CV of 33%
  iiv_thalf  = 0.20,    # ASSUMED (no between-subject estimate published)
  # Michaelis-Menten formation sensitivity (Nebot 2010; recombinant enzymes)
  mw_imatinib = 493.6,  # g/mol, imatinib free base (1 uM = 493.6 ng/mL)
  km_cyp2c8 = 1.4,      # uM  recombinant CYP2C8 (Nebot 2010)
  km_cyp3a4 = 44,       # uM  recombinant CYP3A4 (Nebot 2010)
  fu_plasma = 0.05,     # -   ~95% plasma protein binding (prescribing information)
  # Formation-rate-limited structure (Delbaldo 2006: CLm/fm 58.6 L/h x 0.55 at
  # steady state, Vm/fm 15.4 L -> t1/2 0.33 h; AUC ratio 0.25 at steady state;
  # Petain 2008 similar: t1/2 0.61 h)
  thalf_formlim = 0.693 * 15.4 / (58.6 * 0.55),
  ratio_formlim = 0.25,
  # Menon-Andersen 2009 two-compartment metabolite disposition at 70 kg
  # (CCAP 63:229-38): CL/F 9.65 L/h, V1/F 11.6 L, Q 2.9 L/h, V2/F 256 L;
  # allometric exponents fixed at 0.75 (CL, Q) and 1 (V).
  ma = list(clm = 9.65, v1m = 11.6, qm = 2.9, v2m = 256)
)

# ----------------------------------------------------------------------------
# POPULATION-DERIVED TROUGH TARGET  (metabolite-aware, NON double-counting)
# ----------------------------------------------------------------------------
# 1000 ng/mL: Picard 2007 (n = 68; 50 chronic-phase patients on 400 mg and 18
#   accelerated-phase patients on 600 mg; ROC cut-off 1002 ng/mL for MMR,
#   sensitivity 77%, specificity 71%). Prospectively used as a dosing target in
#   OPTIM-imatinib (Johnson-Ansah 2022) and recommended by the IATDMCT consensus
#   (Clarke 2021; Cmin >= 1000 ng/mL). It is a POPULATION-DERIVED exposure
#   target, not an individual pharmacodynamic boundary.
# Active moiety = Cp + RPF*Cm is scored against the consistently scaled target
#   thr_moiety = thr_parent * (1 + RPF * R), so that parent and moiety criteria
#   coincide at steady state and the metabolite can contribute only through its
#   slower decay (no double counting).
EFF <- list(
  thr_parent = 1000.0,
  thr_sens   = c(750, 1000, 1257),     # v1.0.0 sensitivity set (1257: Golabchifar 2014)
  thr_sweep  = sort(unique(c(seq(500, 1500, by = 100), 750, 1257)))
)

# ----------------------------------------------------------------------------
# CYP3A4 DRUG-DRUG INTERACTION (EXPLORATORY; PARENT ONLY)
# ----------------------------------------------------------------------------
# rifampicin (Bolton 2004): imatinib AUC0-inf -74% -> CL/F x3.85
# ketoconazole (Dutreix 2004; single 200 mg imatinib dose): CL/F -28.6% -> x0.71.
# At steady state, imatinib's own mechanism-based CYP3A4 inhibition (Filppula
# 2012) leaves less CYP3A4 activity to inhibit: ritonavir did not change
# steady-state imatinib exposure (van Erp 2007). The inhibitor multiplier is
# therefore an UPPER bound for steady state. The metabolite is not modified.
DDI <- list(
  states = c("none", "inducer", "inhibitor"),
  cl_mult = c(none = 1.0, inducer = 3.85, inhibitor = 0.71)
)

# ----------------------------------------------------------------------------
# COVARIATES
# ----------------------------------------------------------------------------
# Base case: empirical weight effects of Schmidli 2005 (exponent 0.301 on CL,
# applied also to Q; 0.405 on V, applied to Vc and Vp); no sex effect (sex was
# significant on univariate screening but not retained in the final model).
# Hb and WBC (also in that final model) are not represented: the virtual
# population is characterised by weight only, so their effects are subsumed in
# the random effects. Age and renal function were not significant covariates.
COV <- list(
  wt_ref   = 70,
  allo_CL  = 0.301, # Schmidli 2005, Table 4 (SE 0.096)
  allo_V   = 0.405, # Schmidli 2005, Table 4 (SE 0.118)
  female_cl = 1.0,  # no sex effect (not retained in Schmidli 2005 final model)
  wt_mean  = 70, wt_sd = 13,
  age_mean = 55, age_sd = 13,   # descriptor only
  female_frac = 0.42
)
# The v1.0.0 (submitted) covariate/variability model, kept as a sensitivity analysis.
COV_SUBMITTED <- list(allo_CL = 0.75, allo_V = 1.0, female_cl = 0.92, corr_CL_Vc = 0)

# ----------------------------------------------------------------------------
# CALIBRATION-ANCHOR (EXPOSURE-HEADROOM) SENSITIVITY
# ----------------------------------------------------------------------------
# CL/F is scaled so that the simulated MEDIAN 400 mg trough hits each anchor:
#   860 (Huang 2023 window floor), 980 (IRIS mean 979), 1120, 1300, 1450
#   (Picard responder mean 1452); 600 mg anchors scaled x1.5.
CAL <- list(targets = c(860, 980, 1120, 1300, 1450))

# ----------------------------------------------------------------------------
# SIMULATION SETTINGS
# ----------------------------------------------------------------------------
SIM <- list(
  n_subjects = 1000,
  doses_mg   = c(400, 600),
  tau        = 24,
  seed       = 20260606,
  ss_days     = 21,       # run-in of perfect dosing before any omission (adequacy is checked)
  ss_days_long = 35,      # run-in used for the 90-h metabolite half-life case
  pattern_days = 30,      # adherence-pattern window (after the run-in)
  max_consec_miss = 10,   # base case k = 0..10 (no patient tolerates more than ~6)
  max_consec_sens = 6,    # sensitivity analyses k = 0..6 (truncation reported)
  horizons_h  = c(48, 168, 336, 672),  # trajectory-metric horizons from the first omitted dose
  n_uncertainty = 200,
  n_unc_subjects = 300
)

# Omission behaviours (R1.4/R1.5): late intake by 6/12/18 h, or omission with
# the next dose doubled (the label advises against doubling).
OMIT <- list(late_h = c(6, 12, 18))

# Inter-occasion variability in relative bioavailability (per dose, log-normal)
#   28%: intra-patient trough CV 27.8% at 400 mg in TOPS (Guilhot 2012);
#        intra-individual AUC CV 27.8% (Golabchifar 2016)
#   44%: IOV on F estimated by Menon-Andersen 2009
IOV <- list(cv_F = c(0.278, 0.44))

# Stochastic adherence (R1.3/R2.9): per-patient random omissions at matched mean
# adherence, independent (Bernoulli) or clustered (two-state Markov chain with
# lag-1 autocorrelation rho of the omission indicator; Girard 1998; Wong 2003).
STOCH <- list(adherence = c(0.8, 0.4), rho = c(0, 0.5))

# ----------------------------------------------------------------------------
# PARAMETER UNCERTAINTY (R2.7): every RSE traced to a source
# ----------------------------------------------------------------------------
# Typical values are sampled log-normally (log-scale SD = RSE); variances and
# weight exponents normally (truncated at > 0). The four calibrated disposition
# parameters have no published RSE that applies to their calibrated values, so
# a second run applies a 25% stress RSE to them.
UNC <- data.frame(
  parameter = c("CL", "Vc", "ka", "Q", "Vp", "ss_ratio", "t_half",
                "omega2_CL", "omega2_V", "cov_CL_V", "allo_CL", "allo_V"),
  estimate  = c(10.5, 160, 0.5, 50, 110, 0.15, 40,
                0.102, 0.099, 0.071, 0.301, 0.405),
  rse       = c(0.07, 0.04, 0.10, 0.11, 0.09, 0.15, 0.10,
                0.0142/0.102, 0.0157/0.099, 0.0104/0.071, 0.096/0.301, 0.118/0.405),
  stress_rse = c(0.07, 0.25, 0.25, 0.25, 0.25, 0.15, 0.10,
                 0.0142/0.102, 0.0157/0.099, 0.0104/0.071, 0.096/0.301, 0.118/0.405),
  distribution = c(rep("log-normal", 7), rep("normal, truncated > 0", 5)),
  source = c("Golabchifar 2016 (RSE CL/F 7%); Schmidli 2005 day-29 CL implies ~6.5%",
             "Golabchifar 2016 (RSE Vc/F 4%); value calibrated",
             "Golabchifar 2016 (RSE ka 10%); value calibrated",
             "Golabchifar 2016 (RSE Q/F 11%); value calibrated",
             "Golabchifar 2016 (RSE Vp/F 9%); value calibrated",
             "assumed (between-source spread 0.15-0.27)",
             "assumed (between-source range 20-95 h examined separately)",
             "Schmidli 2005 Table 4 (SE 0.0142)", "Schmidli 2005 Table 4 (SE 0.0157)",
             "Schmidli 2005 Table 4 (SE 0.0104)", "Schmidli 2005 Table 4 (SE 0.096)",
             "Schmidli 2005 Table 4 (SE 0.118)"),
  stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------
# PATHS  (scripts are run from the pipeline root)
# ----------------------------------------------------------------------------
DIRS <- list(results = "results", figures = "figures")
for (d in DIRS) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("Config loaded (v1.1.0): parent CL/F", PK$CL, "L/h; weight exponents",
    COV$allo_CL, "/", COV$allo_V, "; CL-V correlation", PK$corr_CL_Vc,
    "; metabolite t1/2", MET$t_half, "h; R", MET$ss_ratio,
    "; target", EFF$thr_parent, "ng/mL\n")
