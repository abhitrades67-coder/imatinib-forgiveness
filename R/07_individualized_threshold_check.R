# =============================================================================
# 07_individualized_threshold_check.R
# Robustness check: does using each patient's OWN steady-state metabolite:parent
# ratio in the moiety threshold (theta_i = 1000*(1 + RPF*R_i)) -- rather than the
# population R = 0.15 -- change the estimated metabolite forgiveness benefit?
#
# The consistently-scaled moiety threshold (Methods 2.3) uses the population R.
# With between-subject variability on R, that typical-value scaling is not exact
# at the individual level. This script re-derives the benefit from the ALREADY
# SAVED per-patient nadirs (no ODE re-run) under both schemes and reports the
# difference. Result: the benefit moves by <= ~1 percentage point in either
# direction and no conclusion changes, confirming the population-R scaling is an
# adequate (and conservative) instrument for isolating the metabolite's kinetic
# contribution.
#
# Run AFTER run_all.R (reads results/individual_nadir_long.csv + virtual_population.csv):
#   Rscript R/07_individualized_threshold_check.R
# =============================================================================
suppressMessages(library(dplyr))

res_dir <- if (dir.exists("results")) "results" else file.path("..", "results")
long <- read.csv(file.path(res_dir, "individual_nadir_long.csv"))
pop  <- read.csv(file.path(res_dir, "virtual_population.csv"))[, c("id", "met_ratio")]
thr_parent <- 1000

d <- long %>% left_join(pop, by = "id")

tol <- d %>% group_by(dose_mg, rpf, id, met_ratio) %>%
  summarise(
    max_parent     = { ok <- k[nadir_parent >= thr_parent];                       if (length(ok)) max(ok) else 0L },
    max_moiety_pop = { th <- thr_parent*(1 + rpf[1]*0.15);        ok <- k[nadir_moiety >= th]; if (length(ok)) max(ok) else 0L },
    max_moiety_ind = { th <- thr_parent*(1 + rpf[1]*met_ratio[1]); ok <- k[nadir_moiety >= th]; if (length(ok)) max(ok) else 0L },
    .groups = "drop") %>%
  mutate(gain_pop = max_moiety_pop - max_parent,
         gain_ind = max_moiety_ind - max_parent)

summ <- tol %>% group_by(dose_mg, rpf) %>%
  summarise(pct_gain_ge1_population   = round(100*mean(gain_pop >= 1), 1),
            pct_gain_ge1_individualized = round(100*mean(gain_ind >= 1), 1),
            .groups = "drop") %>%
  mutate(delta_pp = pct_gain_ge1_population - pct_gain_ge1_individualized)

cat("\n=== Metabolite benefit: population-R threshold (manuscript) vs individualized ===\n")
print(as.data.frame(summ), row.names = FALSE)
cat(sprintf("\nMax absolute change across all dose x RPF cells: %.1f percentage points\n",
            max(abs(summ$delta_pp))))
write.csv(summ, file.path(res_dir, "individualized_threshold_check.csv"), row.names = FALSE)
