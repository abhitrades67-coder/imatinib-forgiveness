# =============================================================================
# 07_individualized_threshold_check.R
# Robustness check: does scaling the moiety target with each patient's OWN
# metabolite:parent ratio (theta_i = 1000 * (1 + RPF * R_i)) rather than the
# population R change the estimated metabolite benefit? Post hoc: re-derives
# the benefit from the saved base-case nadirs (no ODE re-run).
# Run after run_all.R (reads results/base_individual_long.csv and
# results/virtual_population.csv), or sourced by run_all.R.
# =============================================================================
suppressMessages({library(dplyr); library(data.table)})

individualized_threshold_check <- function(res_dir = "results", thr = 1000, R_pop = 0.15) {
  long <- fread(file.path(res_dir, "base_individual_long.csv"))
  pop  <- fread(file.path(res_dir, "virtual_population.csv"))[, .(id, met_ratio)]
  d <- merge(long, pop, by = "id")
  out <- list()
  for (rp in c(0.1, 0.5, 1.0)) {
    col <- sprintf("nadir_moiety_r%.1f", rp)
    tol <- d[, .(at = nadir_parent[k == 0] >= thr,
                 mp  = { ok <- k[nadir_parent >= thr]; if (length(ok)) max(ok) else 0L },
                 mpop = { th <- thr * (1 + rp * R_pop); ok <- k[get(col) >= th]; if (length(ok)) max(ok) else 0L },
                 mind = { th <- thr * (1 + rp * met_ratio[1]); ok <- k[get(col) >= th]; if (length(ok)) max(ok) else 0L }),
             by = .(dose_mg, id)]
    out[[length(out) + 1]] <- tol[at == TRUE, .(rpf = rp,
      cond_pct_gain_ge1_population = 100 * mean(mpop - mp >= 1),
      cond_pct_gain_ge1_individualized = 100 * mean(mind - mp >= 1)), by = dose_mg]
  }
  res <- rbindlist(out)
  res[, delta_pp := cond_pct_gain_ge1_population - cond_pct_gain_ge1_individualized]
  write.csv(res, file.path(res_dir, "individualized_threshold_check.csv"), row.names = FALSE)
  print(res)
  invisible(res)
}

if (sys.nframe() == 0) individualized_threshold_check()
