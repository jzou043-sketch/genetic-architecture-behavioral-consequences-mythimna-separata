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

segregation_data <- data.frame(
  Population = c(
    rep("Parent", 2),
    rep("F2", 2),
    rep("Backcross", 8)
  ),
  
  Female = c(
    "AA", "EE",
    "AE", "EA",
    "EE", "EE",
    "EA", "EA",
    "AA", "AA",
    "AE", "AE"
  ),
  
  Male = c(
    "AA", "EE",
    "AE", "EA",
    "AE", "EA",
    "AA", "EE",
    "AE", "EA",
    "AA", "EE"
  ),
  
  AE = c(
    0, 0,
    274, 253,
    53, 50,
    41, 58,
    28, 58,
    57, 40
  ),
  
  A = c(
    45, 0,
    91, 87,
    0, 0,
    59, 0,
    48, 52,
    45, 0
  ),
  
  E = c(
    0, 51,
    101, 105,
    48, 46,
    0, 39,
    0, 0,
    0, 58
  ),
  
  Null = c(
    0, 0,
    35, 40,
    0, 0,
    0, 0,
    0, 0,
    0, 0
  ),
  
  Expected_ratio = c(
    "0:1:0:0",
    "0:0:1:0",
    "9:3:3:1",
    "9:3:3:1",
    "1:0:1:0",
    "1:0:1:0",
    "1:1:0:0",
    "1:0:1:0",
    "1:1:0:0",
    "1:1:0:0",
    "1:1:0:0",
    "1:0:1:0"
  ),
  
  stringsAsFactors = FALSE
)

# The observed-category order is:
# AE : A : E : Null

chi_square_test <- function(observed, ratio_text) {
  
  expected_ratio <- as.numeric(
    strsplit(
      ratio_text,
      ":",
      fixed = TRUE
    )[[1]]
  )
  
  # Retain only categories with non-zero expected probabilities
  keep <- expected_ratio > 0
  
  observed <- observed[keep]
  expected_ratio <- expected_ratio[keep]
  
  # No valid chi-square test when only one category is expected
  if (length(observed) < 2) {
    
    return(
      c(
        Chi_square = NA_real_,
        df = 0,
        P_value = NA_real_
      )
    )
  }
  
  test_result <- chisq.test(
    x = observed,
    p = expected_ratio / sum(expected_ratio)
  )
  
  return(
    c(
      Chi_square = unname(
        test_result$statistic
      ),
      df = unname(
        test_result$parameter
      ),
      P_value = test_result$p.value
    )
  )
}

chi_square_results <- t(
  sapply(
    seq_len(nrow(segregation_data)),
    function(i) {
      
      observed <- as.numeric(
        segregation_data[
          i,
          c("AE", "A", "E", "Null")
        ]
      )
      
      chi_square_test(
        observed,
        segregation_data$Expected_ratio[i]
      )
    }
  )
)

final_chi_square <- cbind(
  segregation_data,
  as.data.frame(chi_square_results)
)

final_chi_square$Chi_square <- round(
  final_chi_square$Chi_square,
  digits = 4
)

final_chi_square$P_value <- round(
  final_chi_square$P_value,
  digits = 4
)

print(
  final_chi_square,
  row.names = FALSE
)

write.csv(
  final_chi_square,
  chi_square_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nChi-square results saved to:\n",
  chi_square_file,
  "\n"
)