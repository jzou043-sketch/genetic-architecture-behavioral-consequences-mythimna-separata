# ============================================================
# Chemical phenotype comparison, PCA and segregation analysis
# ============================================================

# Install once if needed:
# install.packages(c(
#   "readxl", "dplyr", "tidyr", "ggplot2",
#   "patchwork", "cowplot"
# ))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# ============================================================
# Part 1. File paths
# ============================================================

analysis_dir <- paste0(
  "C:/Users/15280/OneDrive/Desktop/Data availability/",
  "Parental pheromone types differ in chemical composition and gene expression"
)

chemical_file <- file.path(
  analysis_dir,
  "supplement_A_E_only.xlsx"
)

figure_file <- file.path(
  analysis_dir,
  "Figure_2bc_chemical_PCA_1000dpi.tif"
)

chi_square_file <- file.path(
  analysis_dir,
  "Chi_square_segregation_results.csv"
)

# ============================================================
# Part 2. Read chemical phenotype data
# ============================================================

chemical_data <- read_excel(
  chemical_file,
  sheet = "Individual data"
) %>%
  mutate(
    Type = factor(
      Type,
      levels = c("A", "E"),
      labels = c("A-Type", "E-Type")
    )
  )

# Check sample sizes
print(table(chemical_data$Type))

# Expected:
# A-Type = 250
# E-Type = 216

type_colors <- c(
  "A-Type" = "#FF6B0A",
  "E-Type" = "#1389E9"
)

# ============================================================
# Part 3. Panel b: chemical phenotype comparisons
# ============================================================

chemical_metrics <- c(
  "16:Ald",
  "Z11-16:OH",
  "Z11-16:AC",
  "16:AC",
  "Z11-16:Ald"
)

# Four continuously measured variables:
# Welch's two-sample t-test
#
# Z11-16:Ald:
# Fisher's exact test because it is fixed at 100 in A-Type
# and 0 in E-Type.

calculate_p_value <- function(metric) {
  
  test_data <- chemical_data %>%
    select(Type, value = all_of(metric)) %>%
    filter(!is.na(Type), !is.na(value))
  
  if (metric == "Z11-16:Ald") {
    
    presence_table <- table(
      test_data$Type,
      test_data$value > 0
    )
    
    p_value <- fisher.test(presence_table)$p.value
    
  } else {
    
    values_A <- test_data$value[
      test_data$Type == "A-Type"
    ]
    
    values_E <- test_data$value[
      test_data$Type == "E-Type"
    ]
    
    p_value <- t.test(
      values_A,
      values_E,
      alternative = "two.sided",
      var.equal = FALSE
    )$p.value
  }
  
  return(p_value)
}

chemical_p_values <- sapply(
  chemical_metrics,
  calculate_p_value
)

print(chemical_p_values)

format_p_value <- function(p) {
  
  if (p < 0.001) {
    return("italic(P) < 0.001")
  }
  
  paste0(
    "italic(P) == ",
    formatC(p, digits = 3, format = "f")
  )
}

make_chemical_plot <- function(metric) {
  
  plot_data <- chemical_data %>%
    select(Type, value = all_of(metric)) %>%
    filter(!is.na(Type), !is.na(value))
  
  p_value <- chemical_p_values[[metric]]
  
  observed_max <- max(
    plot_data$value,
    na.rm = TRUE
  )
  
  y_upper <- ifelse(
    observed_max == 0,
    1,
    observed_max * 1.18
  )
  
  bracket_y <- y_upper * 0.88
  bracket_tip <- y_upper * 0.035
  label_y <- y_upper * 0.96
  
  ggplot(
    plot_data,
    aes(
      x = Type,
      y = value,
      color = Type,
      fill = Type
    )
  ) +
    geom_boxplot(
      width = 0.55,
      linewidth = 0.75,
      alpha = 0.15,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.13,
      height = 0,
      size = 1.25,
      alpha = 0.88
    ) +
    
    # Significance bracket
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      linewidth = 0.65,
      color = "black"
    ) +
    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = bracket_y,
      yend = bracket_y - bracket_tip,
      linewidth = 0.65,
      color = "black"
    ) +
    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = bracket_y,
      yend = bracket_y - bracket_tip,
      linewidth = 0.65,
      color = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = label_y,
      label = format_p_value(p_value),
      parse = TRUE,
      size = 3.2,
      fontface = "bold"
    ) +
    scale_color_manual(
      values = type_colors
    ) +
    scale_fill_manual(
      values = type_colors
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = metric,
      x = NULL,
      y = "Relative amount"
    ) +
    theme_classic(
      base_family = "Arial",
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        size = 10,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      axis.title.y = element_text(
        size = 9,
        face = "bold",
        margin = margin(r = 5)
      ),
      axis.text.x = element_text(
        size = 9,
        face = "bold",
        color = c("#FF6B0A", "#1389E9")
      ),
      axis.text.y = element_text(
        size = 8,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.6,
        color = "black"
      ),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}

chemical_plot_list <- lapply(
  chemical_metrics,
  make_chemical_plot
)

panel_b <- wrap_plots(
  chemical_plot_list,
  nrow = 1
)

# ============================================================
# Part 4. Panel c: principal component analysis
# ============================================================

# Z11-16:Ald is excluded from PCA because it was assigned
# directly according to phenotype type:
# A-Type = 100 and E-Type = 0.
#
# Including it would create circular separation.

pca_variables <- c(
  "16:Ald",
  "Z11-16:OH",
  "Z11-16:AC",
  "16:AC"
)

pca_input <- chemical_data %>%
  select(
    Type,
    all_of(pca_variables)
  ) %>%
  drop_na()

pca_fit <- prcomp(
  pca_input[, pca_variables],
  center = TRUE,
  scale. = TRUE
)

variance_explained <- (
  pca_fit$sdev^2 /
    sum(pca_fit$sdev^2)
) * 100

print(variance_explained)

# Expected:
# PC1 approximately 58.2%
# PC2 approximately 19.3%

pca_scores <- as.data.frame(
  pca_fit$x[, 1:2]
)

pca_scores$Type <- pca_input$Type

# Set PC1 direction to match the figure:
# A-Type on the left and E-Type on the right.
if (
  mean(
    pca_scores$PC1[
      pca_scores$Type == "A-Type"
    ]
  ) >
  mean(
    pca_scores$PC1[
      pca_scores$Type == "E-Type"
    ]
  )
) {
  pca_scores$PC1 <- -pca_scores$PC1
}

pca_sample_sizes <- table(
  pca_scores$Type
)

pca_legend_labels <- c(
  "A-Type" = paste0(
    "A-Type (n = ",
    pca_sample_sizes["A-Type"],
    ")"
  ),
  "E-Type" = paste0(
    "E-Type (n = ",
    pca_sample_sizes["E-Type"],
    ")"
  )
)

panel_c <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = Type
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#AAB7C4"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#AAB7C4"
  ) +
  stat_ellipse(
    aes(fill = Type),
    geom = "polygon",
    type = "norm",
    level = 0.95,
    alpha = 0.08,
    color = NA,
    show.legend = FALSE
  ) +
  stat_ellipse(
    type = "norm",
    level = 0.95,
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    size = 1.45,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = type_colors,
    labels = pca_legend_labels
  ) +
  scale_fill_manual(
    values = type_colors
  ) +
  labs(
    title = "PCA scores",
    x = paste0(
      "PC1 (",
      sprintf("%.1f", variance_explained[1]),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      sprintf("%.1f", variance_explained[2]),
      "%)"
    ),
    color = NULL
  ) +
  theme_bw(
    base_family = "Arial",
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 9.5,
      face = "bold"
    ),
    axis.text = element_text(
      size = 8.5,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    panel.grid.major = element_line(
      linewidth = 0.35,
      color = "#E3EBF2"
    ),
    panel.grid.minor = element_blank(),
    legend.position = c(0.18, 0.88),
    legend.background = element_rect(
      fill = "white",
      color = "#C9D3DD",
      linewidth = 0.4
    ),
    legend.text = element_text(
      size = 8.5,
      face = "bold"
    ),
    plot.margin = margin(6, 8, 5, 6)
  )

# ============================================================
# Part 5. Combine panels b and c
# ============================================================

figure_bc <- plot_grid(
  panel_b,
  panel_c,
  ncol = 1,
  labels = c("b", "c"),
  label_fontfamily = "Arial",
  label_fontface = "bold",
  label_size = 15,
  rel_heights = c(1, 1.15),
  align = "v",
  axis = "lr"
)

print(figure_bc)

ggsave(
  filename = figure_file,
  plot = figure_bc,
  device = "tiff",
  width = 18,
  height = 12.5,
  units = "cm",
  dpi = 1000,
  compression = "lzw",
  bg = "white"
)

cat(
  "\nFigure saved to:\n",
  figure_file,
  "\n"
)

# ============================================================
# Part 6. Chi-square goodness-of-fit test
# ============================================================
# ============================================================
# Chi-square goodness-of-fit tests
# Mythimna separata crossing data
# Column order: AE, A, E, N
# ============================================================


# -------------------------
# 1. F2: AE × AE
# Observed = 274, 91, 101, 35
# Expected ratio = 9:3:3:1
# -------------------------

F2_AE <- c(
  AE = 274,
  A  = 91,
  E  = 101,
  N  = 35
)

test_F2_AE <- chisq.test(
  x = F2_AE,
  p = c(9, 3, 3, 1) / 16
)

test_F2_AE
test_F2_AE$expected


# -------------------------
# 2. F2: EA × EA
# Observed = 253, 87, 105, 40
# Expected ratio = 9:3:3:1
# -------------------------

F2_EA <- c(
  AE = 253,
  A  = 87,
  E  = 105,
  N  = 40
)

test_F2_EA <- chisq.test(
  x = F2_EA,
  p = c(9, 3, 3, 1) / 16
)

test_F2_EA
test_F2_EA$expected


# ============================================================
# Backcrosses
# For all backcrosses, test only the two expected phenotype classes
# Expected ratio = 1:1
# correct = FALSE gives Pearson chi-square without Yates correction
# ============================================================


# -------------------------
# 3. EE female × AE male
# AE:E = 54:47
# -------------------------

BC1 <- c(
  AE = 54,
  E  = 47
)

test_BC1 <- chisq.test(
  BC1,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC1


# -------------------------
# 4. EE female × EA male
# AE:E = 50:46
# -------------------------

BC2 <- c(
  AE = 50,
  E  = 46
)

test_BC2 <- chisq.test(
  BC2,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC2


# -------------------------
# 5. EA female × AA male
# AE:A = 41:59
# -------------------------

BC3 <- c(
  AE = 41,
  A  = 59
)

test_BC3 <- chisq.test(
  BC3,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC3


# -------------------------
# 6. EA female × EE male
# AE:E = 58:42
# -------------------------

BC4 <- c(
  AE = 58,
  E  = 42
)

test_BC4 <- chisq.test(
  BC4,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC4


# -------------------------
# 7. AA female × AE male
# AE:A = 28:48
# -------------------------

BC5 <- c(
  AE = 28,
  A  = 48
)

test_BC5 <- chisq.test(
  BC5,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC5


# -------------------------
# 8. AA female × EA male
# AE:A = 58:52
# -------------------------

BC6 <- c(
  AE = 58,
  A  = 52
)

test_BC6 <- chisq.test(
  BC6,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC6


# -------------------------
# 9. AE female × AA male
# AE:A = 56:45
# -------------------------

BC7 <- c(
  AE = 56,
  A  = 45
)

test_BC7 <- chisq.test(
  BC7,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC7


# -------------------------
# 10. AE female × EE male
# AE:E = 40:58
# -------------------------

BC8 <- c(
  AE = 40,
  E  = 58
)

test_BC8 <- chisq.test(
  BC8,
  p = c(0.5, 0.5),
  correct = FALSE
)

test_BC8
# ============================================================
# Summary table
# ============================================================

get_result <- function(cross, obs, test, expected_ratio) {
  
  data.frame(
    Cross = cross,
    N = sum(obs),
    Observed = paste(
      names(obs),
      obs,
      sep = "=",
      collapse = "; "
    ),
    Expected_ratio = expected_ratio,
    Chi_square = round(as.numeric(test$statistic), 3),
    df = as.numeric(test$parameter),
    P_value = signif(test$p.value, 4)
  )
}


results <- rbind(
  
  get_result(
    "F2: AE x AE",
    F2_AE,
    test_F2_AE,
    "9:3:3:1"
  ),
  
  get_result(
    "F2: EA x EA",
    F2_EA,
    test_F2_EA,
    "9:3:3:1"
  ),
  
  get_result(
    "Backcross: EE x AE",
    BC1,
    test_BC1,
    "1:1"
  ),
  
  get_result(
    "Backcross: EE x EA",
    BC2,
    test_BC2,
    "1:1"
  ),
  
  get_result(
    "Backcross: EA x AA",
    BC3,
    test_BC3,
    "1:1"
  ),
  
  get_result(
    "Backcross: EA x EE",
    BC4,
    test_BC4,
    "1:1"
  ),
  
  get_result(
    "Backcross: AA x AE",
    BC5,
    test_BC5,
    "1:1"
  ),
  
  get_result(
    "Backcross: AA x EA",
    BC6,
    test_BC6,
    "1:1"
  ),
  
  get_result(
    "Backcross: AE x AA",
    BC7,
    test_BC7,
    "1:1"
  ),
  
  get_result(
    "Backcross: AE x EE",
    BC8,
    test_BC8,
    "1:1"
  )
)

print(results)

