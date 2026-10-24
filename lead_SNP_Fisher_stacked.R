# =====================================================
# Fisher's exact tests and visualization
# No additional R packages required
# =====================================================

# Contingency tables
ald_table <- matrix(
  c(46, 7,
    24, 23),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Genotype = c("TT", "GG+GT"),
    Phenotype = c("Positive", "Negative")
  )
)

ester_table <- matrix(
  c(77, 10,
    4, 9),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Genotype = c("AC+AA", "CC"),
    Phenotype = c("Positive", "Negative")
  )
)

# Fisher's exact tests
ald_test <- fisher.test(ald_table)
ester_test <- fisher.test(ester_table)

print(ald_test)
print(ester_test)

# Plotting function
draw_figure <- function() {
  
  par(
    mfrow = c(1, 2),
    mar = c(5.5, 5, 4.5, 1.5),
    family = "sans",
    las = 1
  )
  
  # ---------- SNP12656 ----------
  ald_counts <- rbind(
    "Absent (0)"  = c(7, 23),
    "Present (1)" = c(46, 24)
  )
  
  ald_prop <- sweep(
    ald_counts,
    2,
    colSums(ald_counts),
    "/"
  )
  
  bp1 <- barplot(
    ald_prop,
    col = c("#D9EAF7", "#087FC3"),
    border = "white",
    names.arg = c("TT", "GG+GT"),
    ylim = c(0, 1.16),
    width = 0.65,
    space = 0.45,
    axes = FALSE,
    main = "Aldehyde production | SNP12656",
    xlab = "Genotype group",
    ylab = "Proportion of individuals (%)",
    font.main = 2,
    font.lab = 2,
    cex.main = 1.05,
    cex.lab = 1.00
  )
  
  axis(
    2,
    at = seq(0, 1, 0.2),
    labels = paste0(seq(0, 100, 20)),
    font = 2,
    lwd = 1.2
  )
  
  axis(1, at = bp1, labels = c("TT", "GG+GT"), font = 2)
  box(bty = "l", lwd = 1.2)
  
  # Percentages
  for (j in seq_len(ncol(ald_prop))) {
    cumulative <- cumsum(ald_prop[, j])
    midpoint <- cumulative - ald_prop[, j] / 2
    
    for (i in seq_len(nrow(ald_prop))) {
      text(
        bp1[j],
        midpoint[i],
        sprintf("%.1f%%", ald_prop[i, j] * 100),
        col = ifelse(i == 2, "white", "black"),
        font = 2,
        cex = 0.90
      )
    }
    
    text(
      bp1[j],
      1.06,
      paste0("n = ", sum(ald_counts[, j])),
      font = 2,
      cex = 0.95
    )
  }
  
  mtext(
    paste0(
      "Two-sided Fisher's exact test, ",
      format.pval(ald_test$p.value, digits = 3, eps = 0.001)
    ),
    side = 3,
    line = 0.3,
    cex = 0.85
  )
  
  legend(
    "bottom",
    inset = c(0, -0.27),
    legend = rownames(ald_counts),
    fill = c("#D9EAF7", "#087FC3"),
    border = NA,
    horiz = TRUE,
    xpd = TRUE,
    bty = "n",
    cex = 0.85
  )
  
  # ---------- SNP39462 ----------
  ester_counts <- rbind(
    "Absent (0)"  = c(10, 9),
    "Present (1)" = c(77, 4)
  )
  
  ester_prop <- sweep(
    ester_counts,
    2,
    colSums(ester_counts),
    "/"
  )
  
  bp2 <- barplot(
    ester_prop,
    col = c("#FCE5C3", "#F28E00"),
    border = "white",
    names.arg = c("AC+AA", "CC"),
    ylim = c(0, 1.16),
    width = 0.65,
    space = 0.45,
    axes = FALSE,
    main = "Ester production | SNP39462",
    xlab = "Genotype group",
    ylab = "Proportion of individuals (%)",
    font.main = 2,
    font.lab = 2,
    cex.main = 1.05,
    cex.lab = 1.00
  )
  
  axis(
    2,
    at = seq(0, 1, 0.2),
    labels = paste0(seq(0, 100, 20)),
    font = 2,
    lwd = 1.2
  )
  
  axis(1, at = bp2, labels = c("AC+AA", "CC"), font = 2)
  box(bty = "l", lwd = 1.2)
  
  for (j in seq_len(ncol(ester_prop))) {
    cumulative <- cumsum(ester_prop[, j])
    midpoint <- cumulative - ester_prop[, j] / 2
    
    for (i in seq_len(nrow(ester_prop))) {
      text(
        bp2[j],
        midpoint[i],
        sprintf("%.1f%%", ester_prop[i, j] * 100),
        col = ifelse(i == 2, "white", "black"),
        font = 2,
        cex = 0.90
      )
    }
    
    text(
      bp2[j],
      1.06,
      paste0("n = ", sum(ester_counts[, j])),
      font = 2,
      cex = 0.95
    )
  }
  
  mtext(
    paste0(
      "Two-sided Fisher's exact test, ",
      format.pval(ester_test$p.value, digits = 3, eps = 0.001)
    ),
    side = 3,
    line = 0.3,
    cex = 0.85
  )
  
  legend(
    "bottom",
    inset = c(0, -0.27),
    legend = rownames(ester_counts),
    fill = c("#FCE5C3", "#F28E00"),
    border = NA,
    horiz = TRUE,
    xpd = TRUE,
    bty = "n",
    cex = 0.85
  )
}

# 1. Directly display in R
windows(width = 9, height = 4.8)
draw_figure()

# 2. Save as 1000-DPI TIFF
tiff(
  filename = "Binary_Fisher_exact_test.tif",
  width = 9,
  height = 4.8,
  units = "in",
  res = 1000,
  compression = "lzw",
  bg = "white"
)

draw_figure()
dev.off()
