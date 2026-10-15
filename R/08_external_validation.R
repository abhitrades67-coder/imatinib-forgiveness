# =============================================================================
# 08_external_validation.R
# External CONSISTENCY check (not validation of missed-dose trajectories): the
# simulated steady-state trough distribution, with 26% proportional residual
# error added so that it is comparable with measured concentrations, against
# independent published CML cohorts (all values verified in the full texts):
#   IRIS / Larson 2008 (Blood 111:4022-8): day-29 trough median 879 ng/mL
#     [IQR 647-1170], n = 351
#   TOPS / Guilhot 2012 (Haematologica 97:731-8): median Cmin 1200 ng/mL at an
#     actual 400 mg/day and 1935 ng/mL at 600 mg/day (months 1-12 pooled)
#   Huang 2023 (J Clin Pharmacol 63:1438-47): 41.4% of a real-world TDM cohort
#     (n = 102) within 860-1500 ng/mL
# Reads the per-subject qualification output; writes Fig. S2.
# =============================================================================
suppressMessages({library(dplyr); library(ggplot2)})

EXT_REF <- data.frame(
  dose = c(400, 400, 600),
  source = c("IRIS (Larson 2008)", "TOPS (Guilhot 2012)", "TOPS (Guilhot 2012)"),
  median = c(879, 1200, 1935), q1 = c(647, NA, NA), q3 = c(1170, NA, NA))

fig_external <- function(per) {
  d <- per %>% transmute(dose = dose_lab(dose_mg), trough = trough_obs)
  pts <- EXT_REF %>% mutate(dose = dose_lab(dose),
                            xpos = ifelse(grepl("IRIS", source), 0.82, 1.18))
  ggplot(d, aes(x = 1, y = trough)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 860, ymax = 1500, fill = "#e8a33d", alpha = 0.15) +
    geom_hline(yintercept = EFF$thr_parent, linetype = "dashed", colour = "grey40") +
    geom_violin(fill = "#7fb3d5", colour = NA, alpha = 0.5, width = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.3) +
    geom_errorbar(data = pts, aes(x = xpos, ymin = q1, ymax = q3, colour = source),
                  inherit.aes = FALSE, width = 0.06, na.rm = TRUE) +
    geom_point(data = pts, aes(x = xpos, y = median, colour = source, shape = source),
               inherit.aes = FALSE, size = 3) +
    facet_wrap(~dose) +
    scale_colour_manual(values = c("IRIS (Larson 2008)" = "#1b9e77", "TOPS (Guilhot 2012)" = "#d95f02")) +
    scale_shape_manual(values = c("IRIS (Larson 2008)" = 16, "TOPS (Guilhot 2012)" = 17)) +
    coord_cartesian(ylim = c(0, as.numeric(quantile(per$trough_obs[per$dose_mg == 600], .97)))) +
    labs(x = NULL, y = "Steady-state parent trough (ng/mL)", colour = NULL, shape = NULL) +
    THEME + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
