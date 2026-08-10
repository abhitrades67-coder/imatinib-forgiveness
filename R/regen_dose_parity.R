# =============================================================================
# regen_dose_parity.R -- regenerate the (previously 400 mg-only) outputs at BOTH
# doses, reusing the saved virtual population (no full forgiveness re-run needed).
#   - concentration-profile figure (both doses)
#   - model qualification: 400 mg vs reference + 600 mg model prediction
#   - calibration-target sensitivity (both doses)
# =============================================================================
source("R/00_config.R")
source("R/01_model.R")
source("R/02_population.R")
source("R/03_forgiveness.R")
source("R/05_qualification.R")
source("R/04_figures.R")
suppressMessages(library(dplyr))

pop <- read.csv(file.path(DIRS$results, "virtual_population.csv"))
cat(sprintf("Loaded population N=%d\n", nrow(pop)))

cat("\n[Qualification: 400 mg measured + 600 mg predicted]\n")
qual <- qualify_model(pop)
print(qual$table)

cat("\n[Calibration-target sensitivity: both doses]\n")
cal <- bind_rows(lapply(SIM$doses_mg, function(d)
  run_calibration_sensitivity(pop, d, targets = CAL$targets * (d / 400))))
write.csv(cal, file.path(DIRS$results, "calibration_sensitivity.csv"), row.names = FALSE)
print(cal)

cat("\n[Regenerating figures at both doses]\n")
gsave(fig_concentration_profile(pop), "fig1_concentration_profile.png")
gsave(fig_qualification(qual$per), "fig8_model_qualification.png", w = 9)
gsave(fig_calibration(cal), "fig11_calibration_sensitivity.png")
cat("regen_dose_parity.R complete\n")
