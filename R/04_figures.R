# =============================================================================
# 04_figures.R  -- publication figures for the metabolite-aware forgiveness study
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(tidyr); library(ggplot2)})

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))
PAL <- c(Parent = "#1b9e77", Metabolite = "#d95f02", `Active moiety` = "#7570b3")
# Redundant (non-colour) encoding so every series stays distinguishable in
# greyscale and under colour-vision deficiency, per the journal's accessibility rule.
LTY <- c(Parent = "solid", Metabolite = "dotted", `Active moiety` = "22")
SHP <- c(Parent = 16, Metabolite = 15, `Active moiety` = 17)
gsave <- function(p, file, w = 8, h = 5.5, dpi = 300)
  ggsave(file.path(DIRS$figures, file), p, width = w, height = h, dpi = dpi)

# Journal-ready artwork. The manuscript figures are bitmapped line drawings /
# line-halftone combinations, so they are exported at 1000 dpi (the strictest
# Elsevier class) rather than the 300 dpi photograph minimum.
SUBMIT_FIG_DIR <- file.path("Biochemical pharmacology", "Figures")
tsave <- function(p, file, w = 8, h = 5.5, dpi = 1000, dir = SUBMIT_FIG_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  # Uses only grDevices (no extra package dependency). Every figure paints an
  # opaque background, so the exported TIFFs carry no alpha channel, which some
  # print workflows render as black.
  ggsave(file.path(dir, file), p, width = w, height = h, dpi = dpi,
         device = grDevices::tiff, compression = "lzw", type = "cairo")
}

SCEN_ORDER <- c("perfect", "moderate", "poor", "weekend", "delayed")
SCEN_LAB   <- c(perfect="Perfect\n(100%)", moderate="Moderate\n(80%)",
                poor="Poor\n(40%)", weekend="Weekend\nholiday", delayed="Delayed\n(6 h)")

# ---------------------------------------------------------------------------
# Figure 1 -- conceptual: parent vs metabolite vs moiety during a dose holiday
#   Representative (median-CL) patient, 400 mg QD, 3 consecutive missed doses.
# ---------------------------------------------------------------------------
fig_concentration_profile <- function(pop) {
  med <- pop[which.min(abs(pop$cl - median(pop$cl))), ]
  ss <- SIM$ss_days; k <- 3; recover <- 8
  total <- ss + k + recover
  taken <- rep(TRUE, total); taken[(ss+1):(ss+k)] <- FALSE
  dose_times <- (which(taken) - 1) * SIM$tau
  rpf <- MET$rpf_base
  prof <- function(dose) {
    ev <- add.sampling(et(time = dose_times, amt = dose, cmt = "depot"),
                       seq(0, total*SIM$tau, by = 0.1))
    s <- as.data.frame(rxSolve(imatinib_model, unlist(med[c("ka","cl","vc","q","vp","kf","kme","rpf")]), ev))
    s$day <- s$time / 24
    s %>% transmute(day, Parent = Cp, Metabolite = Cm, `Active moiety` = Cp + rpf*Cm,
                    dose = paste0(dose, " mg QD")) %>%
      pivot_longer(c(Parent, Metabolite, `Active moiety`), names_to = "Analyte", values_to = "conc")
  }
  df <- bind_rows(lapply(SIM$doses_mg, prof))
  df$Analyte <- factor(df$Analyte, levels = names(PAL))
  df$dose <- factor(df$dose, levels = paste0(SIM$doses_mg, " mg QD"))

  ggplot(df, aes(day, conc, colour = Analyte, linetype = Analyte)) +
    annotate("rect", xmin = ss, xmax = ss + k, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_hline(yintercept = EFF$thr_parent, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.8) +
    facet_wrap(~dose) +
    scale_colour_manual(values = PAL) +
    scale_linetype_manual(values = LTY) +
    coord_cartesian(xlim = c(ss - 4, ss + k + 6)) +
    labs(x = "Time (days)", y = "Concentration (ng/mL)",
         colour = NULL, linetype = NULL) +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 2 -- target attainment by scenario (parent vs moiety), faceted by dose
# ---------------------------------------------------------------------------
fig_target_attainment <- function(metrics) {
  d <- metrics %>%
    select(dose_mg, scenario, Pct_Above_Parent, Pct_Above_Moiety) %>%
    pivot_longer(c(Pct_Above_Parent, Pct_Above_Moiety),
                 names_to = "criterion", values_to = "pct") %>%
    mutate(criterion = recode(criterion, Pct_Above_Parent = "Parent trough",
                              Pct_Above_Moiety = "Active moiety"),
           scenario = factor(scenario, levels = SCEN_ORDER, labels = SCEN_LAB[SCEN_ORDER]),
           dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(scenario, pct, fill = criterion, linetype = criterion)) +
    geom_col(position = position_dodge(0.75), width = 0.7,
             colour = "grey20", linewidth = 0.35) +
    facet_wrap(~dose) +
    scale_fill_manual(values = c("Parent trough" = "#1b9e77", "Active moiety" = "#7570b3")) +
    scale_linetype_manual(values = c("Parent trough" = "solid", "Active moiety" = "22")) +
    labs(x = NULL, y = "Patients above target (%)", fill = NULL, linetype = NULL) +
    ylim(0, 100) + THEME
}

# ---------------------------------------------------------------------------
# Figure 3 -- waterfall of parent-trough reduction vs perfect adherence
# ---------------------------------------------------------------------------
fig_waterfall <- function(metrics) {
  d <- metrics %>% group_by(dose_mg) %>%
    mutate(perfect = Median_Parent_Trough[scenario == "perfect"],
           pct_red = 100 * (Median_Parent_Trough / perfect - 1)) %>%
    ungroup() %>% filter(scenario != "perfect") %>%
    mutate(scenario = factor(scenario, levels = SCEN_ORDER, labels = SCEN_LAB[SCEN_ORDER]),
           dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(scenario, pct_red, fill = dose)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", pct_red)),
              position = position_dodge(0.75), vjust = 1.2, size = 3) +
    scale_fill_manual(values = c("400 mg QD" = "#386cb0", "600 mg QD" = "#f0027f")) +
    labs(x = NULL, y = "Change in median parent trough vs perfect adherence (%)",
         fill = NULL, title = "Impact of non-adherence on parent trough concentration") +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 4 -- distribution of max tolerable missed doses (parent vs moiety)
# ---------------------------------------------------------------------------
fig_max_tolerable <- function(tol_all) {
  d <- tol_all %>% filter(rpf == MET$rpf_base) %>%
    select(dose_mg, max_tol_parent, max_tol_moiety) %>%
    pivot_longer(c(max_tol_parent, max_tol_moiety),
                 names_to = "criterion", values_to = "doses") %>%
    mutate(criterion = recode(criterion, max_tol_parent = "Parent only",
                              max_tol_moiety = "Metabolite-aware"),
           dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(factor(doses), fill = criterion)) +
    geom_bar(position = position_dodge(preserve = "single"), width = 0.8) +
    facet_wrap(~dose) +
    scale_fill_manual(values = c("Parent only" = "#1b9e77", "Metabolite-aware" = "#7570b3")) +
    labs(x = "Max tolerable consecutive missed doses",
         y = "Number of patients", fill = NULL,
         title = "Individual forgiveness distribution",
         subtitle = sprintf("RPF = %.1f; tolerable = nadir stays >= target", MET$rpf_base)) +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 5 -- cohort nadir vs number of missed doses (parent vs moiety)
# ---------------------------------------------------------------------------
fig_cohort_nadir <- function(cohort_all) {
  d <- cohort_all %>% filter(rpf == MET$rpf_base) %>%
    select(dose_mg, k, nadir_parent_mean, nadir_parent_sd,
           nadir_moiety_mean, nadir_moiety_sd) %>%
    pivot_longer(-c(dose_mg, k),
                 names_to = c("analyte", ".value"),
                 names_pattern = "nadir_(parent|moiety)_(mean|sd)") %>%
    mutate(analyte = recode(analyte, parent = "Parent", moiety = "Active moiety"),
           dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(k, mean, colour = analyte, fill = analyte)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = analyte), linewidth = 0.9) +
    geom_point(aes(shape = analyte), size = 1.8) +
    facet_wrap(~dose) +
    geom_hline(yintercept = EFF$thr_parent, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("Parent" = "#1b9e77", "Active moiety" = "#7570b3")) +
    scale_fill_manual(values = c("Parent" = "#1b9e77", "Active moiety" = "#7570b3")) +
    scale_linetype_manual(values = c("Parent" = "solid", "Active moiety" = "22")) +
    scale_shape_manual(values = c("Parent" = 16, "Active moiety" = 17)) +
    scale_x_continuous(breaks = 0:SIM$max_consec_miss) +
    labs(x = "Consecutive missed doses", y = "Mean nadir concentration (ng/mL)",
         colour = NULL, fill = NULL, linetype = NULL, shape = NULL) +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 6 -- metabolite forgiveness gain across the contested potency range
# ---------------------------------------------------------------------------
fig_potency_sensitivity <- function(tol_all) {
  d <- tol_all %>% group_by(dose_mg, rpf) %>%
    summarise(mean_gain = mean(metab_gain),
              pct_gain1 = 100 * mean(metab_gain >= 1), .groups = "drop") %>%
    mutate(dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(factor(rpf), mean_gain, fill = dose, linetype = dose)) +
    geom_col(position = position_dodge(0.7), width = 0.6,
             colour = "grey20", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.1f%%\n≥ 1", pct_gain1)),
              position = position_dodge(0.7), vjust = -0.2, size = 3) +
    scale_fill_manual(values = c("400 mg QD" = "#386cb0", "600 mg QD" = "#f0027f")) +
    scale_linetype_manual(values = c("400 mg QD" = "solid", "600 mg QD" = "22")) +
    # Headroom above the tallest bar so the two-line labels are never clipped.
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = "Assumed relative potency factor of CGP74588\n(contested: 'significantly less active' to 'potency similar to parent')",
         y = "Mean extra tolerable missed doses\n(metabolite-aware vs parent-only)",
         fill = NULL, linetype = NULL) +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 7 -- threshold sensitivity: % tolerating 0 missed doses across targets
# ---------------------------------------------------------------------------
fig_threshold_sensitivity <- function(thr_all) {
  d <- thr_all %>% mutate(dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(factor(threshold), pct_tol0, fill = dose)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%", pct_tol0)),
              position = position_dodge(0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("400 mg QD" = "#386cb0", "600 mg QD" = "#f0027f")) +
    labs(x = "Efficacy threshold (ng/mL parent trough)",
         y = "Patients tolerating 0 missed doses (%)", fill = NULL,
         title = "Forgiveness is robust across plausible efficacy thresholds",
         subtitle = "Thresholds: 750 (lower bound), 1000 (Picard 2007), 1257 (Golabchifar 2014)") +
    expand_limits(y = 0) + THEME
}

# ---------------------------------------------------------------------------
# Figure 8 -- model qualification: simulated SS exposure vs published ranges
# ---------------------------------------------------------------------------
fig_qualification <- function(per) {
  if (!"dose_mg" %in% names(per)) per$dose_mg <- 400
  lev <- c("SS trough (ng/mL)", "Cmax (ng/mL)", "AUC0-24 (ug.h/mL)")
  dose_lev <- paste0(sort(unique(per$dose_mg)), " mg QD")
  d <- data.frame(
    panel = rep(lev, each = nrow(per)),
    value = c(per$trough, per$Cmax, per$AUC24),
    dose  = rep(paste0(per$dose_mg, " mg QD"), times = 3))
  d$panel <- factor(d$panel, levels = lev)
  d$dose  <- factor(d$dose, levels = dose_lev)
  # Published reference ranges exist only for 400 mg (no independent 600 mg
  # population-exposure dataset); the 400 mg band is drawn over the 400 mg violin only.
  ref <- data.frame(
    panel = factor(lev, levels = lev),
    lo = c(979, 1920, 34), hi = c(1452, 2040, 40),
    xpos = which(dose_lev == "400 mg QD"))
  ggplot(d, aes(x = dose, y = value)) +
    geom_rect(data = ref, inherit.aes = FALSE,
              aes(xmin = xpos - 0.45, xmax = xpos + 0.45, ymin = lo, ymax = hi),
              fill = "#e8a33d", alpha = 0.35) +
    geom_violin(aes(fill = dose), alpha = 0.5, colour = NA, show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.size = 0.4) +
    facet_wrap(~panel, scales = "free_y") +
    scale_fill_manual(values = c("400 mg QD" = "#7fb3d5", "600 mg QD" = "#a6d854")) +
    labs(x = NULL, y = NULL,
         title = "Model qualification (400 mg) and dose-proportional prediction (600 mg)",
         subtitle = "Orange band = published 400 mg reference range (le Coutre 2004; Egorin 2009; Picard 2007); 600 mg is a model prediction with no independent reference") +
    THEME
}

# ---------------------------------------------------------------------------
# Figure 9 -- EXPLORATORY CYP3A4 DDI effect on forgiveness (400 mg)
# ---------------------------------------------------------------------------
fig_ddi <- function(ddi) {
  d <- ddi %>% mutate(state = factor(state, levels = c("inhibitor","none","inducer"),
                       labels = c("Strong\ninhibitor", "No\ninteraction", "Strong\ninducer")),
                       dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(state, pct_above_target, fill = state)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.0f%% above\n%.0f%% tol 0", pct_above_target, pct_tol0)),
              vjust = -0.2, size = 2.7) +
    facet_wrap(~dose) +
    scale_fill_manual(values = c("#2c7fb8", "grey60", "#d95f0e")) +
    labs(x = NULL, y = "Patients above 1000 ng/mL target (%)",
         title = "Exploratory: CYP3A4 interactions dominate the forgiveness picture",
         subtitle = "Upper-bound CL/F effects (Bolton 2004; FDA label). Inhibitor e.g. ketoconazole; inducer e.g. rifampicin. Exploratory only.") +
    expand_limits(y = 115) + THEME
}

# ---------------------------------------------------------------------------
# Figure 10 -- parameter-uncertainty propagation of the headline metric
# ---------------------------------------------------------------------------
fig_uncertainty <- function(unc) {
  d <- unc; if (!"dose_mg" %in% names(d)) d$dose_mg <- 400
  d$dose <- paste0(d$dose_mg, " mg QD")
  ggplot(d, aes(pct_tol0)) +
    geom_histogram(bins = 30, fill = "#7570b3", colour = "white") +
    facet_wrap(~dose) +
    labs(x = "Patients tolerating 0 missed doses (%)",
         y = "Parameter-uncertainty replicates",
         title = "Headline result is robust to parameter uncertainty",
         subtitle = "Distribution across resampled parameter sets (each dose)") +
    THEME
}

# ---------------------------------------------------------------------------
# Model schematic (compartment diagram), drawn with ggplot primitives
# ---------------------------------------------------------------------------
fig_schematic <- function() {
  bw <- 1.5; bh <- 0.85
  box <- function(x, y, lab, fill) data.frame(x, y, w = bw, h = bh, lab, fill)
  boxes <- rbind(
    box(1.0, 3.0, "Depot\n(gut)",        "#eaf3fb"),
    box(3.5, 3.0, "Central\nparent",     "#cfe3f5"),
    box(3.5, 1.3, "Peripheral\nparent",  "#eaf3fb"),
    box(6.2, 3.0, "Metabolite\nCGP74588","#fde3c8"))
  seg <- function(x, y, xe, ye) data.frame(x, y, xe, ye)
  arrows_df <- rbind(
    seg(1.0 + bw/2, 3.0, 3.5 - bw/2, 3.0),    # depot -> central (ka)
    seg(3.5 + bw/2, 3.0, 6.2 - bw/2, 3.0),    # central -> metabolite (kf)
    seg(3.35, 3.0 - bh/2, 3.35, 1.3 + bh/2),  # central -> peripheral (Q)
    seg(3.65, 1.3 + bh/2, 3.65, 3.0 - bh/2))  # peripheral -> central (Q)
  ggplot() +
    geom_tile(data = boxes, aes(x, y, width = w, height = h, fill = I(fill)), colour = "grey30") +
    geom_text(data = boxes, aes(x, y, label = lab), size = 3.8, colour = "grey15") +
    geom_segment(data = arrows_df, aes(x = x, y = y, xend = xe, yend = ye),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    # elimination arrows pointing to explicit "Eliminated" sinks
    annotate("segment", x = 3.5, y = 3.0 + bh/2, xend = 3.5, yend = 3.74,
             arrow = arrow(length = unit(0.18,"cm")), linewidth = 0.5) +
    annotate("text", x = 3.5, y = 3.92, label = "Eliminated", size = 3, fontface = "italic", colour = "grey35") +
    annotate("segment", x = 6.2, y = 3.0 - bh/2, xend = 6.2, yend = 1.78,
             arrow = arrow(length = unit(0.18,"cm")), linewidth = 0.5) +
    annotate("text", x = 6.2, y = 1.6, label = "Eliminated", size = 3, fontface = "italic", colour = "grey35") +
    # rate-constant labels
    annotate("text", x = 2.25, y = 3.18, label = "k[a]", parse = TRUE, size = 3.3) +
    annotate("text", x = 4.85, y = 3.20, label = "k[f]~'(CYP3A4)'", parse = TRUE, size = 3.3) +
    annotate("text", x = 3.0,  y = 2.15, label = "Q", size = 3.3) +
    annotate("text", x = 3.85, y = 3.72, label = "CL/F", size = 3.0) +
    annotate("text", x = 6.55, y = 2.4,  label = "k[me]", parse = TRUE, size = 3.0) +
    annotate("text", x = 2.6,  y = 3.62, label = "parent t1/2 ~ 18 h", size = 3.1, colour = "grey40") +
    annotate("text", x = 6.2,  y = 3.62, label = "metabolite t1/2 ~ 40 h", size = 3.1, colour = "grey40") +
    # xlim trimmed on the right so the schematic sits centred in the canvas.
    coord_cartesian(xlim = c(0.1, 7.12), ylim = c(0.6, 4.1)) +
    theme_void() +
    # theme_void() leaves the canvas transparent; paint it so the exported TIFF
    # carries no alpha channel.
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

# ---------------------------------------------------------------------------
# Graphical abstract (single-panel key message)
# ---------------------------------------------------------------------------
fig_graphical_abstract <- function(cohort_all, tol_all) {
  d <- cohort_all %>% filter(rpf == MET$rpf_base, dose_mg == 400) %>%
    select(k, Parent = nadir_parent_mean, `Active moiety` = nadir_moiety_mean) %>%
    tidyr::pivot_longer(-k, names_to = "Analyte", values_to = "nadir")
  p1 <- ggplot(d, aes(k, nadir, colour = Analyte)) +
    geom_hline(yintercept = 1000, linetype = "dashed", colour = "grey40") +
    annotate("text", x = 10, y = 1000, hjust = 1, vjust = -0.6, size = 5,
             colour = "grey40", label = "Efficacy target (1000 ng/mL)") +
    geom_line(aes(linetype = Analyte), linewidth = 1.1) +
    geom_point(aes(shape = Analyte), size = 2.6) +
    scale_colour_manual(values = c("Parent" = "#1b9e77", "Active moiety" = "#7570b3")) +
    scale_linetype_manual(values = c("Parent" = "solid", "Active moiety" = "22")) +
    scale_shape_manual(values = c("Parent" = 16, "Active moiety" = 17)) +
    scale_x_continuous(breaks = 0:10) +
    labs(x = "Consecutive missed doses (400 mg QD)", y = "Mean nadir (ng/mL)",
         colour = NULL, linetype = NULL, shape = NULL,
         title = "Does imatinib's long-lived metabolite make it more forgiving?",
         subtitle = "Largely no: the active moiety crosses the target only marginally later than parent") +
    # Type is sized so every element still clears 7 pt when the graphical
    # abstract is displayed at the journal's stated 5 x 13 cm.
    theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          strip.background = element_rect(fill = "grey92"),
          plot.title = element_text(size = 22, face = "bold"),
          plot.subtitle = element_text(size = 17),
          # keep the legend tight against the axis so the panel fills the canvas
          legend.box.spacing = unit(4, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(6, 10, 4, 6))
  p1
}

# ---------------------------------------------------------------------------
# Time below threshold vs consecutive missed doses (clinical co-metric)
# ---------------------------------------------------------------------------
fig_time_below <- function(cohort_all) {
  d <- cohort_all %>% filter(rpf == MET$rpf_base) %>%
    select(dose_mg, k, Parent = time_below_parent_mean, `Active moiety` = time_below_moiety_mean) %>%
    tidyr::pivot_longer(c(Parent, `Active moiety`), names_to = "Analyte", values_to = "hours") %>%
    mutate(dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(k, hours, colour = Analyte)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    facet_wrap(~dose) +
    scale_colour_manual(values = c("Parent" = "#1b9e77", "Active moiety" = "#7570b3")) +
    scale_x_continuous(breaks = 0:SIM$max_consec_miss) +
    labs(x = "Consecutive missed doses", y = "Mean time below 1000 ng/mL (h)",
         colour = NULL,
         title = "Cumulative sub-therapeutic time increases with missed doses",
         subtitle = "A graded, clinically interpretable companion to the binary nadir criterion") +
    THEME
}

# ---------------------------------------------------------------------------
# Calibration-target sensitivity of the headline metric
# ---------------------------------------------------------------------------
fig_calibration <- function(cal) {
  if (!"dose_mg" %in% names(cal)) cal$dose_mg <- 400
  cal$dose <- factor(paste0(cal$dose_mg, " mg QD"),
                     levels = paste0(sort(unique(cal$dose_mg)), " mg QD"))
  # IRIS population-mean anchor (979 ng/mL at 400 mg), scaled dose-proportionally
  vl <- data.frame(dose = factor(paste0(sort(unique(cal$dose_mg)), " mg QD"),
                                 levels = paste0(sort(unique(cal$dose_mg)), " mg QD")),
                   x = 979 * (sort(unique(cal$dose_mg)) / 400))
  ggplot(cal, aes(median_trough, pct_tol0)) +
    geom_vline(data = vl, aes(xintercept = x), linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.9, colour = "#386cb0") +
    geom_point(size = 2.2, colour = "#386cb0") +
    geom_text(aes(label = sprintf("%.0f%%", pct_tol0)), vjust = -0.8, size = 3) +
    facet_wrap(~dose, scales = "free_x") +
    labs(x = "Calibrated median steady-state trough (ng/mL)",
         y = "Patients tolerating 0 missed doses (%)",
         title = "Headline depends on the assumed population trough, but stays high",
         subtitle = "CL/F scaled to hit each target median trough; dashed line = IRIS population mean (979 ng/mL at 400 mg, scaled dose-proportionally)") +
    expand_limits(y = c(0, 100)) + THEME
}

make_all_figures <- function(pop, metrics, tol_all, cohort_all,
                             thr_all = NULL, qual = NULL, ddi = NULL, unc = NULL,
                             cal = NULL) {
  cat("  figure 1: concentration profile...\n");  gsave(fig_concentration_profile(pop), "fig1_concentration_profile.png")
  cat("  figure 2: target attainment...\n");      gsave(fig_target_attainment(metrics), "fig2_target_attainment.png")
  # fig3_trough_waterfall (was Supplementary Fig. S8) and fig4_max_tolerable_distribution
  # (was Supplementary Fig. S7) removed from the manuscript as redundant; generation disabled.
  cat("  figure 5: cohort nadir...\n");           gsave(fig_cohort_nadir(cohort_all), "fig5_cohort_nadir.png")
  cat("  figure 6: potency sensitivity...\n");    gsave(fig_potency_sensitivity(tol_all), "fig6_potency_sensitivity.png")
  if (!is.null(thr_all)) { cat("  figure 7: threshold sensitivity...\n"); gsave(fig_threshold_sensitivity(thr_all), "fig7_threshold_sensitivity.png") }
  if (!is.null(qual))    { cat("  figure 8: model qualification...\n");   gsave(fig_qualification(qual), "fig8_model_qualification.png", w = 9) }
  if (!is.null(ddi))     { cat("  figure 9: DDI exploratory...\n");       gsave(fig_ddi(ddi), "fig9_ddi_exploratory.png") }
  if (!is.null(unc))     { cat("  figure 10: uncertainty...\n");          gsave(fig_uncertainty(unc), "fig10_uncertainty.png") }
  if (!is.null(cal))     { cat("  figure 11: calibration sensitivity...\n"); gsave(fig_calibration(cal), "fig11_calibration_sensitivity.png") }
  cat("  figure 12: time below threshold...\n");                          gsave(fig_time_below(cohort_all), "fig12_time_below_threshold.png")
  cat("  schematic + graphical abstract...\n")
  gsave(fig_schematic(), "fig0_model_schematic.png", w = 8, h = 4)
  # graphical abstract: landscape aspect (~2.5:1, readable at 5 x 13 cm)
  gsave(fig_graphical_abstract(cohort_all, tol_all), "graphical_abstract.png", w = 10, h = 4)
  cat("  figures written to", DIRS$figures, "\n")

  # Journal submission artwork. Manuscript Fig. 2 is written by
  # 08_external_validation.R; the remaining five plus the graphical abstract
  # are exported here at 1000 dpi (line-drawing class) / 600 dpi (abstract).
  cat("  submission artwork (TIFF)...\n")
  tsave(fig_schematic(),                        "Figure_1.tiff", w = 8, h = 4)
  tsave(fig_target_attainment(metrics),         "Figure_3.tiff")
  tsave(fig_cohort_nadir(cohort_all),           "Figure_4.tiff")
  tsave(fig_concentration_profile(pop),         "Figure_5.tiff")
  tsave(fig_potency_sensitivity(tol_all),       "Figure_6.tiff")
  tsave(fig_graphical_abstract(cohort_all, tol_all), "5. Graphical Abstract.tiff",
        w = 10, h = 4, dpi = 600, dir = "Biochemical pharmacology")
  cat("  submission artwork written to", SUBMIT_FIG_DIR, "\n")
}
