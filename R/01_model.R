# =============================================================================
# 01_model.R
# rxode2 structural models: parent imatinib (2-compartment, first-order
# absorption) coupled to the active metabolite CGP74588.
#
#   imatinib_model     base case: concentration-driven, linear (first-order)
#                      formation, mono-exponential elimination
#   imatinib_model_mm  sensitivity: Michaelis-Menten (saturable) formation
#   imatinib_model_ma  sensitivity: Menon-Andersen 2009 two-compartment
#                      metabolite disposition on the same parent
# =============================================================================
suppressMessages(library(rxode2))

# Base case. The metabolite is tracked in CONCENTRATION units (ng/mL):
#   dCm/dt = kf*Cp - kme*Cm
# so that  metabolite t1/2 = ln2 / kme  and  steady-state ratio Cm/Cp = kf/kme = R.
# The unidentifiable formation-fraction/bioavailability product is folded into
# R, and formation does not deplete the parent (the parent is calibrated to its
# own observed exposure, which already reflects all elimination pathways).
imatinib_model <- rxode2({
  d/dt(depot)  = -ka * depot
  d/dt(centr)  =  ka * depot - (cl / vc) * centr - (q / vc) * centr + (q / vp) * peri
  d/dt(peri)   =  (q / vc) * centr - (q / vp) * peri
  Cp = 1000 * centr / vc                 # mg/L -> ng/mL
  d/dt(Cm) = kf * Cp - kme * Cm
  moiety = Cp + rpf * Cm
})

# Michaelis-Menten formation (reviewer 1, comment 1). Formation saturates with
# the driving concentration Cd = fdrive * Cp (fdrive = unbound fraction, or 1
# for the extreme "total concentration" bound), expressed in uM:
#   formation = kfmax * Cp / (1 + Cd / km)
# kfmax is set per patient so that the steady-state metabolite:parent AUC ratio
# at 400 mg equals the patient's ratio in the linear model (see
# 03_forgiveness.R, mm_normalise()); at other doses the ratio is then free to
# fall if formation saturates.
imatinib_model_mm <- rxode2({
  d/dt(depot)  = -ka * depot
  d/dt(centr)  =  ka * depot - (cl / vc) * centr - (q / vc) * centr + (q / vp) * peri
  d/dt(peri)   =  (q / vc) * centr - (q / vp) * peri
  Cp = 1000 * centr / vc
  Cd = fdrive * Cp / mw                  # ng/mL / (g/mol) = uM
  d/dt(Cm) = kfmax * Cp / (1 + Cd / km) - kme * Cm
  moiety = Cp + rpf * Cm
})

# Menon-Andersen 2009 metabolite disposition (two compartments, amounts in mg
# of parent-equivalent) on the same parent model. fm is set per patient so that
# the steady-state AUC ratio equals the patient's R (fm itself is not
# identifiable from parent and metabolite concentrations alone).
imatinib_model_ma <- rxode2({
  d/dt(depot)  = -ka * depot
  d/dt(centr)  =  ka * depot - (cl / vc) * centr - (q / vc) * centr + (q / vp) * peri
  d/dt(peri)   =  (q / vc) * centr - (q / vp) * peri
  Cp = 1000 * centr / vc
  d/dt(m1c) = fm * (cl / vc) * centr - (clm / v1m) * m1c - (qm / v1m) * m1c + (qm / v2m) * m2p
  d/dt(m2p) = (qm / v1m) * m1c - (qm / v2m) * m2p
  Cm = 1000 * m1c / v1m
  moiety = Cp + rpf * Cm
})

MODEL_PARAMS <- list(
  base = c("id", "ka", "cl", "vc", "q", "vp", "kf", "kme", "rpf"),
  mm   = c("id", "ka", "cl", "vc", "q", "vp", "kfmax", "kme", "rpf", "fdrive", "mw", "km"),
  ma   = c("id", "ka", "cl", "vc", "q", "vp", "fm", "clm", "v1m", "qm", "v2m", "rpf"))
MODELS <- list(base = imatinib_model, mm = imatinib_model_mm, ma = imatinib_model_ma)

cat("rxode2 models compiled: base, Michaelis-Menten formation, Menon-Andersen metabolite.\n")
