# =============================================================================
# 01_model.R
# rxode2 structural model: parent imatinib (2-compartment) + active metabolite
# CGP74588 (concentration-driven, identifiability-safe).
# =============================================================================
suppressMessages(library(rxode2))

# Model parameters (per individual) supplied at solve time:
#   ka, cl, vc, q, vp   parent PK
#   kf, kme             metabolite formation / elimination rate constants
#   rpf                 relative potency factor (for active-moiety output)
#
# Metabolite is tracked directly in CONCENTRATION units (ng/mL):
#   dCm/dt = kf*Cp - kme*Cm
# so that:
#   metabolite t1/2 = ln2 / kme           (set kme = ln2 / 40 h)
#   SS ratio Cm/Cp  = kf / kme = R        (set kf  = R * kme)
# This makes both metabolite half-life and SS ratio exactly match literature,
# while folding the unidentifiable formation-fraction/F product into the
# calibrated R (no fm/F identifiability problem).

imatinib_model <- rxode2({
  # --- parent disposition (amounts in mg) ---
  d/dt(depot)  = -ka * depot
  d/dt(centr)  =  ka * depot - (cl / vc) * centr - (q / vc) * centr + (q / vp) * peri
  d/dt(peri)   =  (q / vc) * centr - (q / vp) * peri

  # parent plasma concentration: mg/L -> ng/mL (x1000); 1 mg/L = 1000 ng/mL
  Cp = 1000 * centr / vc

  # --- metabolite (concentration units, ng/mL) ---
  d/dt(Cm) = kf * Cp - kme * Cm

  # --- active moiety (parent + potency-weighted metabolite) ---
  moiety = Cp + rpf * Cm
})

cat("rxode2 parent+metabolite model compiled.\n")
