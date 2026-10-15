# =============================================================================
# 04_figures.R  -- figures for the revised article (v1.1.0)
#   Main text:  Fig. 1 published-model benchmark;  Fig. 2 exposure headroom
#   Supplement: Fig. S1 schematic; S2 external consistency; S3 adherence
#               patterns; S4 mean nadir vs missed doses; S5 trajectory metrics;
#               S6 metabolite benefit across sensitivity analyses; S7 target
#               sweep; S8 calibration anchors; S9 omission behaviours;
#               S10 parameter uncertainty
# =============================================================================
suppressMessages({library(rxode2); library(dplyr); library(tidyr); library(ggplot2)})

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        plot.background = element_rect(fill = "white", colour = NA))
PAL <- c(Parent = "#1b9e77", Metabolite = "#d95f02", `Active moiety` = "#7570b3")
gsave <- function(p, file, w = 8, h = 5.5, dpi = 300)
  ggsave(file.path(DIRS$figures, file), p, width = w, height = h, dpi = dpi)

# Journal artwork: 1000 dpi LZW TIFF (line/halftone class) plus an 1800-px-wide
# PNG rendering of the same image for embedding in the Word files.
SUBMIT_FIG_DIR <- file.path(DIRS$figures, "submission")
tsave <- function(p, name, w = 8, h = 5.5) {
  if (!dir.exists(SUBMIT_FIG_DIR)) dir.create(SUBMIT_FIG_DIR, recursive = TRUE)
  ggsave(file.path(SUBMIT_FIG_DIR, paste0(name, ".tiff")), p, width = w, height = h, dpi = 1000,
         device = grDevices::tiff, compression = "lzw", type = "cairo")
  ggsave(file.path(SUBMIT_FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 1800 / w)
}

dose_lab <- function(d) factor(paste0(d, " mg QD"), levels = paste0(SIM$doses_mg, " mg QD"))

# ---------------------------------------------------------------------------
# Fig. 1 -- typical-patient benchmark against published parent-metabolite models
# ---------------------------------------------------------------------------
BM_COL <- c("This study (base case)" = "#1b9e77", "Menon-Andersen 2009" = "#d95f02",
            "Delbaldo 2006" = "#7570b3", "Petain 2008" = "#e7298a")
BM_LTY <- c("This study (base case)" = "solid", "Menon-Andersen 2009" = "22",
            "Delbaldo 2006" = "42", "Petain 2008" = "12")
fig_benchmark <- function(prof) {
  d <- prof %>% filter(time >= -24, time <= 96) %>%
    pivot_longer(c(parent_pct, metab_pct), names_to = "analyte", values_to = "pct") %>%
    mutate(analyte = factor(recode(analyte, parent_pct = "Imatinib (parent)",
                                   metab_pct = "CGP74588 (metabolite)"),
                            levels = c("Imatinib (parent)", "CGP74588 (metabolite)")),
           model = factor(model, levels = names(BM_COL)))
  ggplot(d, aes(time, pct, colour = model, linetype = model)) +
    annotate("rect", xmin = 0, xmax = 48, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_hline(yintercept = 100, linetype = "dotted", colour = "grey40") +
    geom_line(linewidth = 0.8) +
    facet_wrap(~analyte) +
    scale_colour_manual(values = BM_COL) + scale_linetype_manual(values = BM_LTY) +
    scale_x_continuous(breaks = seq(-24, 96, 24)) +
    labs(x = "Time from the first omitted dose (h)",
         y = "Concentration (% of steady-state trough)",
         colour = NULL, linetype = NULL) +
    THEME + guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
}

# ---------------------------------------------------------------------------
# Fig. 2 -- exposure headroom required to keep the target after k omissions
# For a linear model the ratio nadir_k / steady-state trough is a patient-
# specific kinetic property, independent of dose. A patient keeps the target
# after k omissions when (steady-state trough / target) >= required headroom,
# where required headroom = trough / nadir_k (parent) or
# (theta_moiety / theta) * trough / moiety nadir_k (active moiety).
# ---------------------------------------------------------------------------
headroom_data <- function(long, ks = 1:2) {
  k0 <- long[k == 0]
  bind_rows(lapply(ks, function(kk) {
    x <- long[k == kk][match(k0$id, id)]
    bind_rows(
      data.frame(k = kk, criterion = "Parent", h = k0$nadir_parent / x$nadir_parent),
      data.frame(k = kk, criterion = "Active moiety, RPF 0.5",
                 h = (moiety_threshold(0.5) / EFF$thr_parent) * k0$nadir_parent / x$`nadir_moiety_r0.5`),
      data.frame(k = kk, criterion = "Active moiety, RPF 1.0",
                 h = (moiety_threshold(1.0) / EFF$thr_parent) * k0$nadir_parent / x$`nadir_moiety_r1.0`))
  }))
}
fig_headroom <- function(hd, med_headroom) {
  hd <- hd %>% mutate(kl = factor(paste0(k, ifelse(k == 1, " missed dose", " consecutive missed doses"))),
                      criterion = factor(criterion, levels = c("Parent", "Active moiety, RPF 0.5", "Active moiety, RPF 1.0")))
  vl <- data.frame(x = med_headroom, lab = names(med_headroom))
  ggplot(hd, aes(h, colour = criterion, linetype = criterion)) +
    geom_vline(data = vl, aes(xintercept = x), colour = "grey45", linetype = "dotted") +
    geom_text(data = vl, aes(x = x, y = 100, label = lab), inherit.aes = FALSE,
              angle = 90, vjust = -0.4, hjust = 1, size = 3, colour = "grey30") +
    stat_ecdf(aes(y = after_stat(ecdf) * 100), linewidth = 0.85, pad = FALSE) +
    facet_wrap(~kl) +
    scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 6, 8, 10)) +
    coord_cartesian(xlim = c(1, 10)) +   # zoom only: the ECDF uses every patient
    scale_colour_manual(values = c("Parent" = "#1b9e77", "Active moiety, RPF 0.5" = "#7570b3",
                                   "Active moiety, RPF 1.0" = "#e7298a")) +
    scale_linetype_manual(values = c("Parent" = "solid", "Active moiety, RPF 0.5" = "22",
                                     "Active moiety, RPF 1.0" = "12")) +
    labs(x = "Steady-state trough as a multiple of the 1000 ng/mL target (log scale)",
         y = "Patients who would keep the target (%)", colour = NULL, linetype = NULL) +
    THEME
}

# ---------------------------------------------------------------------------
# Fig. S1 -- model schematic
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
    seg(1.0 + bw/2, 3.0, 3.5 - bw/2, 3.0),
    seg(3.5 + bw/2, 3.0, 6.2 - bw/2, 3.0),
    seg(3.35, 3.0 - bh/2, 3.35, 1.3 + bh/2),
    seg(3.65, 1.3 + bh/2, 3.65, 3.0 - bh/2))
  ggplot() +
    geom_tile(data = boxes, aes(x, y, width = w, height = h, fill = I(fill)), colour = "grey30") +
    geom_text(data = boxes, aes(x, y, label = lab), size = 3.8, colour = "grey15") +
    geom_segment(data = arrows_df, aes(x = x, y = y, xend = xe, yend = ye),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    annotate("segment", x = 3.5, y = 3.0 + bh/2, xend = 3.5, yend = 3.74,
             arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    annotate("text", x = 3.5, y = 3.92, label = "Eliminated", size = 3, fontface = "italic", colour = "grey35") +
    annotate("segment", x = 6.2, y = 3.0 - bh/2, xend = 6.2, yend = 1.78,
             arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    annotate("text", x = 6.2, y = 1.6, label = "Eliminated", size = 3, fontface = "italic", colour = "grey35") +
    annotate("text", x = 2.25, y = 3.18, label = "k[a]", parse = TRUE, size = 3.3) +
    annotate("text", x = 4.85, y = 3.20, label = "k[f]", parse = TRUE, size = 3.3) +
    annotate("text", x = 4.85, y = 2.80, label = "(CYP2C8, CYP3A4)", size = 2.7) +
    annotate("text", x = 3.0,  y = 2.15, label = "Q", size = 3.3) +
    annotate("text", x = 3.85, y = 3.72, label = "CL/F", size = 3.0) +
    annotate("text", x = 6.55, y = 2.4,  label = "k[me]", parse = TRUE, size = 3.0) +
    annotate("text", x = 2.6,  y = 3.62, label = "parent t1/2 ~ 18 h", size = 3.1, colour = "grey40") +
    annotate("text", x = 6.2,  y = 3.62, label = "metabolite t1/2 40 h (20-90 h)", size = 3.1, colour = "grey40") +
    coord_cartesian(xlim = c(0.1, 7.3), ylim = c(0.6, 4.1)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
}

# ---------------------------------------------------------------------------
# Fig. S3 -- adherence patterns: time below target in at-target patients
# ---------------------------------------------------------------------------
PATTERN_ORDER <- c("perfect", "delayed", "moderate", "stoch_80_independent", "stoch_80_clustered",
                   "weekend", "poor", "stoch_40_independent", "stoch_40_clustered")
PATTERN_LAB <- c(perfect = "Perfect\n(100%)", delayed = "8 doses\n6 h late", moderate = "Every 5th\nmissed (80%)",
                 stoch_80_independent = "Random\n80%, indep.", stoch_80_clustered = "Random\n80%, clustered",
                 weekend = "Weekends\nmissed (73%)", poor = "2 taken,\n3 missed (40%)",
                 stoch_40_independent = "Random\n40%, indep.", stoch_40_clustered = "Random\n40%, clustered")
fig_patterns <- function(pat_summary) {
  d <- pat_summary %>% select(dose_mg, scenario, Parent = cond_pct_time_below_parent,
                              `Active moiety, RPF 0.5` = `cond_pct_time_below_moiety_r0.5`,
                              `Active moiety, RPF 1.0` = `cond_pct_time_below_moiety_r1.0`) %>%
    pivot_longer(-c(dose_mg, scenario), names_to = "criterion", values_to = "pct") %>%
    mutate(scenario = factor(scenario, levels = PATTERN_ORDER, labels = PATTERN_LAB[PATTERN_ORDER]),
           dose = dose_lab(dose_mg),
           criterion = factor(criterion, levels = c("Parent", "Active moiety, RPF 0.5", "Active moiety, RPF 1.0")))
  ggplot(d, aes(scenario, pct, fill = criterion)) +
    geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey20", linewidth = 0.3) +
    facet_wrap(~dose, ncol = 1) +
    scale_fill_manual(values = c("Parent" = "#1b9e77", "Active moiety, RPF 0.5" = "#7570b3",
                                 "Active moiety, RPF 1.0" = "#e7298a")) +
    labs(x = NULL, y = "Time below target over 30 days (%)\n(patients at target with uninterrupted dosing)",
         fill = NULL) +
    THEME + theme(axis.text.x = element_text(size = 8.5))
}

# ---------------------------------------------------------------------------
# Fig. S4 -- cohort nadir vs number of missed doses (parent vs moiety)
# ---------------------------------------------------------------------------
fig_cohort_nadir <- function(cohort_all) {
  # Each series must be scored against ITS OWN criterion: the parent against the
  # 1000 ng/mL target, the active moiety against the consistently scaled
  # threshold thr_parent * (1 + RPF * ss_ratio). Drawing a single line at 1000
  # would score the moiety against the parent threshold and so re-introduce
  # exactly the double-counting that the scaled threshold exists to prevent.
  thr_m <- EFF$thr_parent * (1 + MET$rpf_base * MET$ss_ratio)
  d <- cohort_all %>% filter(rpf == MET$rpf_base) %>%
    select(dose_mg, k, nadir_parent_mean, nadir_parent_sd,
           nadir_moiety_mean, nadir_moiety_sd) %>%
    pivot_longer(-c(dose_mg, k),
                 names_to = c("analyte", ".value"),
                 names_pattern = "nadir_(parent|moiety)_(mean|sd)") %>%
    mutate(analyte = recode(analyte, parent = "Parent", moiety = "Active moiety"),
           dose = paste0(dose_mg, " mg QD"))
  ggplot(d, aes(k, mean, colour = analyte, fill = analyte)) +
    # Mean +/- SD with the lower bound truncated at zero: nadir concentrations
    # are strictly positive and right-skewed, so an unclamped ribbon dips below
    # the axis into physically impossible values.
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = analyte), linewidth = 0.9) +
    geom_point(aes(shape = analyte), size = 1.8) +
    facet_wrap(~dose) +
    geom_hline(yintercept = EFF$thr_parent, linetype = "dashed",
               colour = "#1b9e77", linewidth = 0.5) +
    geom_hline(yintercept = thr_m, linetype = "dotted",
               colour = "#7570b3", linewidth = 0.7) +
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
# Fig. S5 -- trajectory metrics vs number of missed doses (at-target patients)
# ---------------------------------------------------------------------------
fig_trajectory <- function(traj) {
  sel <- function(prefix_p, prefix_m, metric) {
    traj %>% filter(k >= 1) %>%
      transmute(dose_mg, k, Parent = .data[[prefix_p]],
                `Active moiety, RPF 0.1` = .data[[sprintf(prefix_m, "0.1")]],
                `Active moiety, RPF 0.5` = .data[[sprintf(prefix_m, "0.5")]],
                `Active moiety, RPF 1.0` = .data[[sprintf(prefix_m, "1.0")]]) %>%
      pivot_longer(-c(dose_mg, k), names_to = "criterion", values_to = "value") %>%
      mutate(metric = metric)
  }
  d <- bind_rows(
    sel("cond_mean_excess_tb_parent_672", "cond_mean_excess_tb_moiety_r%s_672",
        "Excess time below target
over 28 days (h)"),
    sel("cond_median_deficit_parent_168", "cond_median_deficit_moiety_r%s_168",
        "Exposure deficit
over 7 days (%)")) %>%
    mutate(dose = dose_lab(dose_mg),
           criterion = factor(criterion, levels = c("Parent", "Active moiety, RPF 0.1",
                                                    "Active moiety, RPF 0.5", "Active moiety, RPF 1.0")))
  ggplot(d, aes(k, value, colour = criterion, linetype = criterion)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
    facet_grid(metric ~ dose, scales = "free_y", switch = "y") +
    scale_colour_manual(values = c("Parent" = "#1b9e77", "Active moiety, RPF 0.1" = "#a6cee3",
                                   "Active moiety, RPF 0.5" = "#7570b3", "Active moiety, RPF 1.0" = "#e7298a")) +
    scale_linetype_manual(values = c("Parent" = "solid", "Active moiety, RPF 0.1" = "42",
                                     "Active moiety, RPF 0.5" = "22", "Active moiety, RPF 1.0" = "12")) +
    scale_x_continuous(breaks = 1:SIM$max_consec_miss) +
    labs(x = "Consecutive missed doses", y = NULL, colour = NULL, linetype = NULL) +
    THEME + theme(strip.placement = "outside") + guides(colour = guide_legend(nrow = 2))
}

# ---------------------------------------------------------------------------
# Fig. S6 -- metabolite benefit across all sensitivity analyses
# ---------------------------------------------------------------------------
SENS_SHORT <- c(base = "Base case", thalf20 = "Metabolite t1/2 20 h", thalf90 = "Metabolite t1/2 90 h",
                ratio025 = "Ratio R = 0.25", formlim = "Formation-rate-limited metabolite",
                ma015 = "Two-compartment metabolite, R 0.15", ma025 = "Two-compartment metabolite, R 0.25",
                mm_unbound = "Saturable formation (unbound)", mm_total = "Saturable formation (total; extreme)",
                submitted = "Submitted (v1.0.0) covariate model", base_qvp = "IIV on Q and Vp",
                iov28 = "IOV on F, CV 28%", iov44 = "IOV on F, CV 44%", metiiv0 = "No metabolite IIV",
                metiiv2 = "Metabolite IIV doubled", ddi_inhibitor = "Strong CYP3A4 inhibitor (upper bound)")
fig_sens_forest <- function(sens) {
  sens <- sens %>% filter(!is.na(cond_pct_gain_ge1_r0.5)) %>%
    mutate(short = ifelse(label %in% names(SENS_SHORT), SENS_SHORT[label], description))
  d <- sens %>% select(short, group, dose_mg, `RPF 0.5` = cond_pct_gain_ge1_r0.5,
                       `RPF 1.0` = cond_pct_gain_ge1_r1.0) %>%
    pivot_longer(c(`RPF 0.5`, `RPF 1.0`), names_to = "rpf", values_to = "pct") %>%
    mutate(dose = dose_lab(dose_mg),
           short = factor(short, levels = rev(unique(sens$short))))
  ggplot(d, aes(pct, short, colour = rpf, shape = rpf)) +
    geom_point(size = 2.2) +
    facet_wrap(~dose) +
    scale_colour_manual(values = c("RPF 0.5" = "#7570b3", "RPF 1.0" = "#e7298a")) +
    labs(x = "At-target patients gaining at least one tolerable missed dose (%)",
         y = NULL, colour = NULL, shape = NULL) +
    THEME + theme(axis.text.y = element_text(size = 8))
}

# ---------------------------------------------------------------------------
# Fig. S7 -- target sweep 500-1500 ng/mL;  Fig. S8 -- calibration anchors
# ---------------------------------------------------------------------------
METRIC_COL <- c("At target with uninterrupted dosing" = "grey30",
                "Lose target after one missed dose (at-target patients)" = "#1b9e77",
                "Gain >= 1 tolerable missed dose, RPF 0.5" = "#7570b3",
                "Gain >= 1 tolerable missed dose, RPF 1.0" = "#e7298a")
metric_long <- function(d, x) {
  d %>% transmute(dose = dose_lab(dose_mg), x = .data[[x]],
                  `At target with uninterrupted dosing` = pct_at_target,
                  `Lose target after one missed dose (at-target patients)` = cond_pct_lose_after1_parent,
                  `Gain >= 1 tolerable missed dose, RPF 0.5` = cond_pct_gain_ge1_r0.5,
                  `Gain >= 1 tolerable missed dose, RPF 1.0` = cond_pct_gain_ge1_r1.0) %>%
    pivot_longer(-c(dose, x), names_to = "metric", values_to = "pct") %>%
    mutate(metric = factor(metric, levels = names(METRIC_COL)))
}
fig_threshold_sweep <- function(sw) {
  ggplot(metric_long(sw, "threshold"), aes(x, pct, colour = metric, linetype = metric)) +
    geom_vline(xintercept = EFF$thr_parent, colour = "grey60", linetype = "dotted") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    facet_wrap(~dose) +
    scale_colour_manual(values = METRIC_COL) +
    scale_linetype_manual(values = c("solid", "solid", "22", "12")) +
    labs(x = "Population-derived trough target (ng/mL)", y = "Patients (%)", colour = NULL, linetype = NULL) +
    ylim(0, 100) + THEME + guides(colour = guide_legend(nrow = 2))
}
fig_calibration <- function(cal) {
  ggplot(metric_long(cal, "median_ss_trough"), aes(x, pct, colour = metric, linetype = metric)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    facet_wrap(~dose, scales = "free_x") +
    scale_colour_manual(values = METRIC_COL) +
    scale_linetype_manual(values = c("solid", "solid", "22", "12")) +
    labs(x = "Calibrated median steady-state trough (ng/mL)", y = "Patients (%)",
         colour = NULL, linetype = NULL) +
    ylim(0, 100) + THEME + guides(colour = guide_legend(nrow = 2))
}

# ---------------------------------------------------------------------------
# Fig. S9 -- omission behaviours (typical patient, 400 mg)
# ---------------------------------------------------------------------------
fig_omission <- function(prof) {
  lev <- c("uninterrupted", "6 h late", "12 h late", "18 h late", "omitted", "omitted, next doubled")
  d <- prof %>% filter(time >= -24, time <= 72) %>% mutate(case = factor(case, levels = lev))
  ggplot(d, aes(time, Cp, colour = case, linetype = case)) +
    geom_hline(yintercept = EFF$thr_parent, colour = "grey40", linetype = "dashed") +
    geom_line(linewidth = 0.75) +
    scale_colour_manual(values = c("grey20", "#a6cee3", "#1f78b4", "#08306b", "#d95f02", "#e7298a")) +
    scale_linetype_manual(values = c("solid", "42", "22", "12", "solid", "dotdash")) +
    scale_x_continuous(breaks = seq(-24, 72, 12)) +
    labs(x = "Time from the scheduled time of the perturbed dose (h)",
         y = "Parent concentration (ng/mL)", colour = NULL, linetype = NULL) +
    THEME + guides(colour = guide_legend(nrow = 2))
}

# ---------------------------------------------------------------------------
# Fig. S10 -- parameter uncertainty (sourced RSEs and 25% stress test)
# ---------------------------------------------------------------------------
fig_uncertainty <- function(unc) {
  d <- unc %>% mutate(dose = dose_lab(dose_mg),
                      run = recode(run, sourced = "Sourced RSEs", stress = "25% stress on calibrated parameters"))
  ggplot(d, aes(cond_pct_lose_after1_parent)) +
    geom_histogram(bins = 30, fill = "#1b9e77", colour = "white") +
    facet_grid(run ~ dose) +
    labs(x = "At-target patients losing the target after one missed dose (%)",
         y = "Uncertainty replicates") +
    THEME
}
