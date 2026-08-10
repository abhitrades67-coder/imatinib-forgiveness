# =============================================================================
# 08_external_validation.R
# External validation of the simulated steady-state trough DISTRIBUTION against
# INDEPENDENT published CML cohorts (beyond the summary-statistic qualification
# in 05_qualification.R / main-text Table 2). Sources (all values verified):
#   - IRIS / Larson 2008 (Blood 2008;111:4022-8; doi:10.1182/blood-2007-10-116475):
#       mean steady-state imatinib trough 979 ng/mL (CV 54%); CGP74588 trough
#       242 ng/mL (CV 44%); trough metabolite:parent ratio ~0.25.
#   - Picard 2007 (Blood 2007;109:3496-9; doi:10.1182/blood-2006-07-036012):
#       responder-subgroup mean trough ~1452 ng/mL.
#   - Huang 2023 (J Clin Pharmacol 2023;63:1438-47; doi:10.1002/jcph.2329):
#       real-world CML TDM cohort (n=102); 41.4% of patients within the
#       recommended 860-1500 ng/mL window.
# Reads results/qualification_per_subject.csv (produced by 05_qualification.R);
# NO re-simulation is required. Writes figures/fig13_external_validation.png and
# results/external_validation.csv.
# =============================================================================
suppressMessages({library(dplyr); library(ggplot2)})

# Published external comparators (verified; see header for full citations).
EXT_REF <- list(
  larson_parent_mean  = 979, larson_parent_cv = 0.54,   # IRIS, Larson 2008
  larson_met_mean     = 242, larson_met_cv    = 0.44,   # IRIS, Larson 2008
  larson_trough_ratio = 0.25,                           # IRIS trough ratio
  picard_responder    = 1452,                           # Picard 2007 responders
  huang_window        = c(860, 1500),                   # Huang 2023 recommended window
  huang_pct_in_window = 41.4, huang_n = 102             # Huang 2023 real-world TDM
)

external_validation <- function(per, dose = 400) {
  d   <- per[per$dose_mg == dose, ]
  thr <- if (exists("EFF")) EFF$thr_parent else 1000
  w   <- EXT_REF$huang_window

  med_iqr <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, .25), quantile(x, .75))
  pct_in  <- 100 * mean(d$trough >= w[1] & d$trough <= w[2])
  pct_ge  <- 100 * mean(d$trough >= thr)
  ratio   <- d$troughm / d$trough

  tab <- data.frame(
    Metric = c("Median [IQR] parent trough (ng/mL)",
               "% within 860-1500 ng/mL window",
               sprintf("%% >= %g ng/mL target", thr),
               "Median [IQR] CGP74588 trough (ng/mL)",
               "Median [IQR] metabolite:parent trough ratio"),
    `Model (400 mg, n=1000)` = c(
      med_iqr(d$trough),
      sprintf("%.0f", pct_in),
      sprintf("%.0f", pct_ge),
      med_iqr(d$troughm),
      sprintf("%.2f [%.2f-%.2f]", median(ratio), quantile(ratio, .25), quantile(ratio, .75))),
    `Independent published value (source)` = c(
      "979 mean, CV 54% (IRIS, Larson 2008); 1452 responders (Picard 2007)",
      sprintf("%.1f%% (real-world TDM; Huang 2023, n=%d)", EXT_REF$huang_pct_in_window, EXT_REF$huang_n),
      "-",
      "242 mean, CV 44% (IRIS, Larson 2008)",
      sprintf("%.2f trough ratio (IRIS, Larson 2008)", EXT_REF$larson_trough_ratio)),
    check.names = FALSE)
  write.csv(tab, file.path(DIRS$results, "external_validation.csv"), row.names = FALSE)

  # ---- figure: simulated 400 mg trough distribution vs published cohorts ----
  src_lev <- c("Model\n(this study)\nn=1000", "IRIS / Larson\n2008", "Picard 2007\n(responders)")
  sim <- data.frame(source = factor(src_lev[1], levels = src_lev), trough = d$trough)
  pts <- data.frame(
    source = factor(c(src_lev[2], src_lev[3]), levels = src_lev),
    mean   = c(EXT_REF$larson_parent_mean, EXT_REF$picard_responder),
    sd     = c(EXT_REF$larson_parent_mean * EXT_REF$larson_parent_cv, NA))

  p <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = w[1], ymax = w[2],
             fill = "#e8a33d", alpha = 0.18) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "grey40") +
    geom_violin(data = sim, aes(source, trough), fill = "#7fb3d5", colour = NA, alpha = 0.5) +
    geom_boxplot(data = sim, aes(source, trough), width = 0.15, outlier.size = 0.4) +
    geom_point(data = pts, aes(source, mean), size = 3, colour = "#1b9e77") +
    geom_errorbar(data = pts, aes(source, ymin = mean - sd, ymax = mean + sd),
                  width = 0.12, colour = "#1b9e77", na.rm = TRUE) +
    coord_cartesian(ylim = c(0, as.numeric(quantile(d$trough, .99)))) +
    # No in-image title/subtitle: this is a numbered main-text figure, and the
    # journal requires the caption title to live in the manuscript, not the image.
    labs(x = NULL, y = "Steady-state parent trough (ng/mL)") +
    (if (exists("THEME")) THEME else theme_bw(base_size = 12))

  gsave2 <- if (exists("gsave")) gsave else
    function(p, f, w = 8, h = 5.5) ggsave(file.path(DIRS$figures, f), p, width = w, height = h, dpi = 300)
  gsave2(p, "fig13_external_validation.png", w = 9, h = 5.5)
  # Submission-ready TIFF (LZW). This is line/halftone artwork, not a photograph,
  # so it is exported at 1000 dpi rather than the 300 dpi photograph minimum.
  ggsave(file.path(DIRS$figures, "fig13_external_validation.tiff"), p,
         width = 8, height = 5.5, dpi = 1000,
         device = grDevices::tiff, compression = "lzw", type = "cairo")
  if (exists("tsave")) tsave(p, "Figure_2.tiff", w = 8, h = 5.5)
  list(table = tab, figure = p)
}
