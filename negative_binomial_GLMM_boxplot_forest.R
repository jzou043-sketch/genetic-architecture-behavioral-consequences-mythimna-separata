suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(broom.mixed)
  library(performance)
  library(DHARMa)
})

options(contrasts = c("contr.treatment", "contr.poly"))
set.seed(20261015)

input_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Field trapping reveals lower attraction of hybrid-type pheromone signals than parental types/Field attraction data.xlsx"
output_dir <- "C:/Users/15280/OneDrive/文档/New project/outputs/Field_trapping_analysis_code_current_input/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

captures <- read_excel(input_file, sheet = "Captures")
trials <- read_excel(input_file, sheet = "Trials")

trial_key <- trials %>%
  transmute(
    trial_id = as.character(trial_id),
    trial_site = as.character(city_or_county),
    trial_start_month = as.character(trial_start_month),
    year = substr(trial_start_month, 1, 4)
  )

dat <- captures %>%
  mutate(trial_id = as.character(trial_id)) %>%
  filter(record_status == "observed", lure_group %in% c("A", "E", "AE", "EA")) %>%
  left_join(trial_key, by = "trial_id") %>%
  mutate(
    sampling_site = if_else(!is.na(trial_site) & trial_site != "", trial_site, as.character(sampling_site)),
    lure = as.character(lure_group),
    lure = factor(lure, levels = c("A", "E", "AE", "EA")),
    sampling_site = factor(sampling_site),
    year = factor(year),
    trial_id = factor(trial_id),
    trap_id = factor(trap_id),
    check_order_num = as.numeric(check_order),
    check_order_z = as.numeric(scale(check_order_num)),
    males_caught = as.integer(males_caught)
  ) %>%
  filter(!is.na(males_caught), !is.na(sampling_site), !is.na(year), !is.na(trap_id), !is.na(check_order_z))

stopifnot(nrow(dat) == 9593L)
stopifnot(sum(dat$males_caught) == 4815L)

# Trial-level standardized catch: mean catch per observed trap-check record.
trial_level <- dat %>%
  group_by(trial_id, trial_start_month, sampling_site, lure) %>%
  summarise(
    total_males = sum(males_caught),
    n_trap_checks = n(),
    n_traps = n_distinct(trap_id),
    n_check_events = n_distinct(check_id),
    standardized_catch = total_males / n_trap_checks,
    .groups = "drop"
  ) %>%
  arrange(lure, trial_start_month, sampling_site)

write.csv(
  trial_level,
  file.path(output_dir, "Trial_level_standardized_catches.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Negative-binomial GLMM. The 33-level trial factor is a fixed block combining
# sampling site and trial-start time; check order adjusts for within-trial time.
# Repeated observations from the same trap share a random intercept.
model_full <- glmmTMB(
  males_caught ~ lure + trial_id + check_order_z + (1 | trap_id),
  family = nbinom2(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS"),
    optCtrl = list(maxit = 20000)
  )
)

model_no_lure <- glmmTMB(
  males_caught ~ trial_id + check_order_z + (1 | trap_id),
  family = nbinom2(link = "log"),
  data = dat,
  control = glmmTMBControl(
    optimizer = optim,
    optArgs = list(method = "BFGS"),
    optCtrl = list(maxit = 20000)
  )
)

lure_lrt <- anova(model_no_lure, model_full)

fixed_effects <- broom.mixed::tidy(
  model_full,
  effects = "fixed",
  component = "cond",
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  rename(rate_ratio = estimate, lower_95_CI = conf.low, upper_95_CI = conf.high)

write.csv(
  fixed_effects,
  file.path(output_dir, "NB_GLMM_fixed_effects_exponentiated.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

emm_link <- emmeans(model_full, ~ lure)
vs_a_raw <- as.data.frame(summary(
  contrast(emm_link, method = "trt.vs.ctrl", ref = 1),
  type = "response", infer = c(TRUE, TRUE)
))
all_pairwise_raw <- as.data.frame(summary(
  pairs(emm_link, reverse = TRUE), type = "response", infer = c(TRUE, TRUE)
))

standardize_ratio_table <- function(x) {
  ratio_col <- intersect(c("ratio", "response", "odds.ratio"), names(x))[1]
  lower_col <- intersect(c("asymp.LCL", "lower.CL", "LCL"), names(x))[1]
  upper_col <- intersect(c("asymp.UCL", "upper.CL", "UCL"), names(x))[1]
  p_col <- intersect(c("p.value", "p.value."), names(x))[1]
  out <- x
  out$rate_ratio <- out[[ratio_col]]
  out$lower_95_CI <- out[[lower_col]]
  out$upper_95_CI <- out[[upper_col]]
  if (!is.na(p_col)) out$p_value <- out[[p_col]]
  out
}

vs_a <- standardize_ratio_table(vs_a_raw) %>%
  mutate(
    comparison = gsub(" / ", " vs ", contrast, fixed = TRUE),
    comparison = factor(comparison, levels = c("E vs A", "AE vs A", "EA vs A"))
  ) %>%
  arrange(comparison)

all_pairwise <- standardize_ratio_table(all_pairwise_raw) %>%
  mutate(
    comparison = gsub(" / ", " vs ", contrast, fixed = TRUE),
    comparison = factor(
      comparison,
      levels = c("E vs A", "AE vs A", "EA vs A", "AE vs E", "EA vs E", "EA vs AE")
    )
  ) %>%
  arrange(comparison)

write.csv(
  vs_a,
  file.path(output_dir, "NB_GLMM_IRR_vs_A.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_pairwise,
  file.path(output_dir, "NB_GLMM_all_pairwise_IRR.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Diagnostics.
simulation <- simulateResiduals(model_full, n = 1000, plot = FALSE, seed = 20261015)
uniformity_test <- testUniformity(simulation)
dispersion_test <- testDispersion(simulation)
zero_inflation_test <- testZeroInflation(simulation)
overdispersion_check <- performance::check_overdispersion(model_full)

random_sd <- as.data.frame(VarCorr(model_full)$cond)$sdcor[1]
lrt_p <- lure_lrt$`Pr(>Chisq)`[2]

diagnostic_lines <- c(
  "Negative-binomial GLMM diagnostics",
  sprintf("Rows analyzed: %d", nrow(dat)),
  sprintf("Trials: %d", n_distinct(dat$trial_id)),
  sprintf("Sampling sites: %d", n_distinct(dat$sampling_site)),
  sprintf("Trap IDs: %d", n_distinct(dat$trap_id)),
  sprintf("Converged (positive-definite Hessian): %s", model_full$sdr$pdHess),
  sprintf("Optimizer convergence code: %s", model_full$fit$convergence),
  sprintf("Random-intercept SD (trap ID): %.6f", random_sd),
  sprintf("Global lure likelihood-ratio test: Chisq=%.4f, df=%d, p=%.6g", lure_lrt$Chisq[2], lure_lrt$`Chi Df`[2], lrt_p),
  sprintf("DHARMa uniformity test p=%.6g", uniformity_test$p.value),
  sprintf("DHARMa dispersion test p=%.6g", dispersion_test$p.value),
  sprintf("DHARMa zero-inflation test p=%.6g", zero_inflation_test$p.value),
  capture.output(print(overdispersion_check)),
  "",
  "Model formula:",
  deparse(formula(model_full)),
  "",
  "Model summary:",
  capture.output(summary(model_full)),
  "",
  "Likelihood-ratio test for lure:",
  capture.output(print(lure_lrt))
)
writeLines(diagnostic_lines, file.path(output_dir, "NB_GLMM_diagnostics.txt"), useBytes = TRUE)

# Figure.
palette_lure <- c(A = "#D55E00", E = "#0072B2", AE = "#009E73", EA = "#8E44AD")
base_theme <- theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.tag = element_text(face = "bold", size = 15),
    legend.position = "none"
  )

p_box <- ggplot(trial_level, aes(x = lure, y = standardized_catch, color = lure, fill = lure)) +
  geom_boxplot(width = 0.58, alpha = 0.18, outlier.shape = NA, linewidth = 0.75) +
  geom_jitter(width = 0.13, height = 0, size = 2.1, alpha = 0.78) +
  stat_summary(fun = median, geom = "crossbar", width = 0.42, linewidth = 0.7, color = "black") +
  scale_color_manual(values = palette_lure) +
  scale_fill_manual(values = palette_lure) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.03, base = 10),
    breaks = c(0, 0.03, 0.1, 0.3, 1, 3, 10, 30),
    labels = label_number(accuracy = 0.01, trim = TRUE),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    title = "Trial-level standardized catches",
    x = "Lure group",
    y = "Males per observed trap check"
  ) +
  base_theme

forest_data <- all_pairwise %>%
  mutate(
    comparison = factor(
      as.character(comparison),
      levels = rev(c("E vs A", "AE vs A", "EA vs A", "AE vs E", "EA vs E", "EA vs AE"))
    ),
    label = sprintf("%.2f (%.2f, %.2f)", rate_ratio, lower_95_CI, upper_95_CI)
  )

x_max <- max(forest_data$upper_95_CI, na.rm = TRUE)
x_min <- min(forest_data$lower_95_CI, na.rm = TRUE)
label_x <- x_max * 1.16

p_forest <- ggplot(forest_data, aes(y = comparison, x = rate_ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.65, color = "grey40") +
  geom_errorbar(aes(xmin = lower_95_CI, xmax = upper_95_CI), width = 0.16, linewidth = 0.85, color = "#1F4E79", orientation = "y") +
  geom_point(shape = 21, size = 3.8, stroke = 0.9, fill = "#2F75B5", color = "#17365D") +
  geom_text(aes(x = label_x, label = label), hjust = 0, family = "Times New Roman", size = 3.7) +
  scale_x_log10(
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8),
    labels = label_number(accuracy = 0.01),
    limits = c(max(0.03, x_min / 1.35), max(64, x_max * 20)),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "All pairwise relative capture effects",
    x = "Capture ratio (95% CI; log scale)",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(face = "bold", size = 11))

global_note <- if (is.na(lrt_p)) {
  "Global lure effect: likelihood-ratio test unavailable"
} else if (lrt_p < 0.001) {
  "Global lure effect: likelihood-ratio test, P < 0.001"
} else {
  sprintf("Global lure effect: likelihood-ratio test, P = %.3f", lrt_p)
}

combined <- (p_box | p_forest) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      text = element_text(family = "Times New Roman"),
      plot.tag = element_text(face = "bold", size = 15)
    )
  )

ggsave(
  file.path(output_dir, "Figure_Boxplot_Forest_6Pairwise_NB_GLMM.tiff"),
  combined, width = 8.2, height = 4.4, units = "in", dpi = 1000,
  compression = "lzw", bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_Boxplot_Forest_6Pairwise_NB_GLMM.png"),
  combined, width = 8.2, height = 4.4, units = "in", dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_Boxplot_Forest_6Pairwise_NB_GLMM.pdf"),
  combined, width = 8.2, height = 4.4, units = "in", device = cairo_pdf, bg = "white"
)

cat("MODEL_COMPLETE\n")
cat(sprintf("n=%d; trials=%d; sites=%d; traps=%d\n", nrow(dat), n_distinct(dat$trial_id), n_distinct(dat$sampling_site), n_distinct(dat$trap_id)))
cat(sprintf("pdHess=%s; global_lure_p=%.8g\n", model_full$sdr$pdHess, lrt_p))
print(trial_level %>% group_by(lure) %>% summarise(n_trials=n(), median=median(standardized_catch), mean=mean(standardized_catch), .groups="drop"))
print(vs_a %>% select(any_of(c("comparison", "rate_ratio", "lower_95_CI", "upper_95_CI", "p_value"))))
cat(sprintf("DHARMa: uniformity p=%.6g; dispersion p=%.6g; zero inflation p=%.6g\n", uniformity_test$p.value, dispersion_test$p.value, zero_inflation_test$p.value))
# Display all six pairwise comparisons
print(
  all_pairwise %>%
    select(
      comparison,
      rate_ratio,
      lower_95_CI,
      upper_95_CI,
      p_value
    )
)

# Display the boxplot and forest plot in RStudio
print(combined)
