############################################################
# Four-class multinomial SNP association analysis
#
# Phenotypes:
#   AE   = reference
#   A
#   E
#   Null
#
# SNP coding:
#   0 = REF/REF
#   1 = REF/ALT
#   2 = ALT/ALT
#
# Model:
#   Multinomial logistic regression
#
# Global SNP test:
#   Likelihood-ratio Chi-square
#
# Class-specific tests:
#   Wald Chi-square
############################################################


############################################################
# 1. Load packages
############################################################

library(readxl)
library(dplyr)
library(writexl)
library(nnet)

cat("============================================\n")
cat("Four-class multinomial association analysis\n")
cat("============================================\n")
dir.create("results", recursive = TRUE, showWarnings = FALSE)


############################################################
# 2. Input files
############################################################

pheno_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/1.Phenotype.xlsx"
geno_file  <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/2.SNP_GENOTYPE1224_012.xlsx"


############################################################
# 3. Read phenotype and genotype data
############################################################

pheno <- read_excel(pheno_file)
geno  <- read_excel(geno_file)

cat("\nPhenotype:", nrow(pheno), "individuals\n")
cat("Genotype:", nrow(geno), "SNPs\n")


############################################################
# 4. Construct four phenotype classes
#
# A- + ee  = A
# A- + E-  = AE
# aa + E-  = E
# aa + ee  = Null
############################################################

pheno <- pheno %>%
  mutate(

    Phenotype4 = case_when(

      `A-Type` == "A-" &
        `E-Type` == "ee" ~ "A",

      `A-Type` == "A-" &
        `E-Type` == "E-" ~ "AE",

      `A-Type` == "aa" &
        `E-Type` == "E-" ~ "E",

      `A-Type` == "aa" &
        `E-Type` == "ee" ~ "Null",

      TRUE ~ NA_character_
    )
  )


############################################################
# 5. Set phenotype factor
#
# AE is reference phenotype
############################################################

pheno$Phenotype4 <- factor(
  pheno$Phenotype4,
  levels = c(
    "AE",
    "A",
    "E",
    "Null"
  )
)

cat("\nFour phenotype counts:\n")

print(
  table(
    pheno$Phenotype4,
    useNA = "ifany"
  )
)


############################################################
# 6. Check duplicated phenotype IDs
############################################################

if(anyDuplicated(pheno$ID)) {

  stop(
    "Duplicated IDs found in phenotype file."
  )
}


############################################################
# 7. Match genotype and phenotype IDs
############################################################

sample_cols <- colnames(geno)[
  colnames(geno) %in% pheno$ID
]

cat(
  "\nMatched individuals:",
  length(sample_cols),
  "\n"
)

if(length(sample_cols) == 0) {

  stop(
    "No matched IDs between phenotype and genotype files."
  )
}


cat("\nPhenotype IDs missing from genotype:\n")

print(
  setdiff(
    pheno$ID,
    sample_cols
  )
)


############################################################
# 8. Reorder phenotype according to genotype sample order
############################################################

pheno2 <- pheno[
  match(
    sample_cols,
    pheno$ID
  ),
]

stopifnot(
  all(
    pheno2$ID == sample_cols
  )
)

cat("\nID matching successful.\n")


############################################################
# 9. Construct genotype matrix
############################################################

gmat <- as.matrix(
  geno[, sample_cols]
)

storage.mode(gmat) <- "numeric"

cat(
  "\nGenotype matrix:",
  nrow(gmat),
  "SNP x",
  ncol(gmat),
  "individuals\n"
)

cat("\nOverall genotype counts:\n")

print(
  table(
    gmat,
    useNA = "ifany"
  )
)


############################################################
# 10. Check genotype coding
############################################################

observed_values <- sort(
  unique(
    as.vector(gmat)
  )
)

cat("\nObserved genotype codes:\n")

print(
  observed_values
)

bad_values <- observed_values[
  !is.na(observed_values) &
    !(observed_values %in% c(0, 1, 2))
]

if(length(bad_values) > 0) {

  stop(
    paste(
      "Unexpected genotype codes:",
      paste(
        bad_values,
        collapse = ", "
      )
    )
  )
}


############################################################
# 11. Phenotype vector
############################################################

y_all <- pheno2$Phenotype4

if(any(is.na(y_all))) {

  warning(
    "Missing phenotype values detected."
  )
}


############################################################
# 12. Fit common null model
#
# Null:
#   Phenotype ~ 1
#
# If all genotype values are observed,
# this same null likelihood can be reused for every SNP.
############################################################

null_dat <- data.frame(
  y = y_all
)

null_fit <- multinom(

  y ~ 1,

  data = null_dat,

  trace = FALSE,

  maxit = 500,

  MaxNWts = 10000
)

null_logLik <- as.numeric(
  logLik(null_fit)
)

cat(
  "\nCommon null log-likelihood:",
  null_logLik,
  "\n"
)


############################################################
# 13. Define exact output column names
############################################################

expected_names <- c(

  "N",
  "N0",
  "N1",
  "N2",

  "LRChiSq",
  "df_global",
  "P_global",

  "Beta_A",
  "SE_A",
  "OR_A",
  "WaldChiSq_A",
  "P_A",

  "Beta_E",
  "SE_E",
  "OR_E",
  "WaldChiSq_E",
  "P_E",

  "Beta_Null",
  "SE_Null",
  "OR_Null",
  "WaldChiSq_Null",
  "P_Null",

  "Converged"
)


############################################################
# 14. Function for failed / unavailable SNP
############################################################

empty_result <- function(
  N,
  N0,
  N1,
  N2,
  convergence = NA
) {

  unname(
    c(

      N,
      N0,
      N1,
      N2,

      NA,
      3,
      NA,

      NA,
      NA,
      NA,
      NA,
      NA,

      NA,
      NA,
      NA,
      NA,
      NA,

      NA,
      NA,
      NA,
      NA,
      NA,

      convergence
    )
  )
}


############################################################
# 15. Single-SNP multinomial analysis
############################################################

multinom_snp <- function(g) {

  ##########################################################
  # Remove missing genotype / phenotype
  ##########################################################

  ok <- !is.na(g) &
        !is.na(y_all)

  g2 <- g[ok]
  y2 <- y_all[ok]


  ##########################################################
  # Sample size
  ##########################################################

  N <- length(y2)


  ##########################################################
  # Genotype counts
  ##########################################################

  N0 <- sum(g2 == 0)
  N1 <- sum(g2 == 1)
  N2 <- sum(g2 == 2)


  ##########################################################
  # Monomorphic SNP cannot be tested
  ##########################################################

  if(length(unique(g2)) < 2) {

    return(
      empty_result(
        N,
        N0,
        N1,
        N2
      )
    )
  }


  ##########################################################
  # All four phenotype classes must exist
  ##########################################################

  if(length(unique(y2)) < 4) {

    return(
      empty_result(
        N,
        N0,
        N1,
        N2
      )
    )
  }


  ##########################################################
  # Analysis dataframe
  ##########################################################

  dat <- data.frame(

    y = factor(
      y2,
      levels = c(
        "AE",
        "A",
        "E",
        "Null"
      )
    ),

    g = g2
  )


  ##########################################################
  # Full multinomial logistic model
  #
  # log[P(A)/P(AE)]       = aA + bA * G
  # log[P(E)/P(AE)]       = aE + bE * G
  # log[P(Null)/P(AE)]    = aN + bN * G
  ##########################################################

  fit <- try(

    suppressWarnings(

      multinom(

        y ~ g,

        data = dat,

        trace = FALSE,

        maxit = 500,

        MaxNWts = 10000
      )
    ),

    silent = TRUE
  )


  ##########################################################
  # Model fitting failed
  ##########################################################

  if(inherits(fit, "try-error")) {

    return(
      empty_result(
        N,
        N0,
        N1,
        N2,
        convergence = 0
      )
    )
  }


  ##########################################################
  # Convergence
  #
  # 1 = converged
  # 0 = not converged
  ##########################################################

  converged <- as.numeric(
    fit$convergence == 0
  )


  ##########################################################
  # Full model log-likelihood
  ##########################################################

  full_logLik <- as.numeric(
    logLik(fit)
  )


  ##########################################################
  # Null log-likelihood
  #
  # If all individuals retained:
  # use common null model
  #
  # Otherwise:
  # refit null model using same individuals
  ##########################################################

  if(all(ok)) {

    ll0 <- null_logLik

  } else {

    null_local <- try(

      suppressWarnings(

        multinom(

          y ~ 1,

          data = dat,

          trace = FALSE,

          maxit = 500,

          MaxNWts = 10000
        )
      ),

      silent = TRUE
    )


    if(inherits(null_local, "try-error")) {

      return(
        empty_result(
          N,
          N0,
          N1,
          N2,
          convergence = converged
        )
      )
    }


    ll0 <- as.numeric(
      logLik(
        null_local
      )
    )
  }


  ##########################################################
  # Global likelihood-ratio test
  #
  # LR = 2 * (LLfull - LLnull)
  ##########################################################

  LRChiSq <- 2 * (
    full_logLik -
      ll0
  )


  ##########################################################
  # Protect against tiny negative numerical errors
  ##########################################################

  if(
    is.finite(LRChiSq) &&
      LRChiSq < 0 &&
      abs(LRChiSq) < 1e-8
  ) {

    LRChiSq <- 0
  }


  ##########################################################
  # Global degrees of freedom
  #
  # Three additional SNP slopes:
  #
  # beta_A
  # beta_E
  # beta_Null
  ##########################################################

  df_global <- 3


  ##########################################################
  # Global P-value
  ##########################################################

  P_global <- pchisq(

    LRChiSq,

    df = df_global,

    lower.tail = FALSE
  )


  ##########################################################
  # Extract coefficient and standard error matrices
  ##########################################################

  sm <- summary(fit)

  beta_mat <- sm$coefficients
  se_mat   <- sm$standard.errors


  ##########################################################
  # Check expected phenotype rows
  ##########################################################

  required_rows <- c(
    "A",
    "E",
    "Null"
  )

  if(
    !all(
      required_rows %in%
        rownames(beta_mat)
    )
  ) {

    return(
      empty_result(
        N,
        N0,
        N1,
        N2,
        convergence = converged
      )
    )
  }


  ##########################################################
  # Function for class-specific SNP effect
  ##########################################################

  extract_class <- function(class_name) {

    beta <- beta_mat[
      class_name,
      "g"
    ]

    se <- se_mat[
      class_name,
      "g"
    ]


    ########################################################
    # Check numerical validity
    ########################################################

    if(
      !is.finite(beta) ||
        !is.finite(se) ||
        se <= 0
    ) {

      return(
        c(
          NA,
          NA,
          NA,
          NA,
          NA
        )
      )
    }


    ########################################################
    # Odds ratio
    ########################################################

    OR <- exp(beta)


    ########################################################
    # Wald Chi-square
    ########################################################

    WaldChiSq <- (
      beta /
        se
    )^2


    ########################################################
    # Wald P-value
    ########################################################

    P <- pchisq(

      WaldChiSq,

      df = 1,

      lower.tail = FALSE
    )


    return(
      unname(
        c(
          beta,
          se,
          OR,
          WaldChiSq,
          P
        )
      )
    )
  }


  ##########################################################
  # A vs AE
  ##########################################################

  A_res <- extract_class(
    "A"
  )


  ##########################################################
  # E vs AE
  ##########################################################

  E_res <- extract_class(
    "E"
  )


  ##########################################################
  # Null vs AE
  ##########################################################

  Null_res <- extract_class(
    "Null"
  )


  ##########################################################
  # Final fixed-length result vector
  ##########################################################

  result_vector <- unname(

    c(

      N,
      N0,
      N1,
      N2,

      LRChiSq,
      df_global,
      P_global,

      A_res,

      E_res,

      Null_res,

      converged
    )
  )


  ##########################################################
  # Safety check
  ##########################################################

  if(
    length(result_vector) !=
      length(expected_names)
  ) {

    stop(
      paste0(
        "Result length error: expected ",
        length(expected_names),
        ", obtained ",
        length(result_vector)
      )
    )
  }


  return(
    result_vector
  )
}


############################################################
# 16. Test SNP1 first
############################################################

cat("\n============================================\n")
cat("Testing SNP1 before genome-wide scan\n")
cat("============================================\n")

test1 <- multinom_snp(
  gmat[1, ]
)

names(test1) <- expected_names

print(
  test1
)


############################################################
# 17. Genome-wide scan
############################################################

cat("\n============================================\n")
cat("Running genome-wide multinomial scan\n")
cat("============================================\n")

start_time <- Sys.time()

multi_scan <- t(

  apply(

    gmat,

    1,

    multinom_snp
  )
)

end_time <- Sys.time()

cat("\nGenome-wide scan finished.\n")

cat("Elapsed time:\n")

print(
  end_time -
    start_time
)


############################################################
# 18. Convert scan result to dataframe
############################################################

multi_scan <- as.data.frame(

  multi_scan,

  check.names = FALSE
)


############################################################
# 19. Force exact column names
############################################################

if(
  ncol(multi_scan) !=
    length(expected_names)
) {

  cat(
    "\nActual number of result columns:",
    ncol(multi_scan),
    "\n"
  )

  stop(
    paste0(
      "Unexpected multinomial output. Expected ",
      length(expected_names),
      " columns but obtained ",
      ncol(multi_scan)
    )
  )
}


colnames(
  multi_scan
) <- expected_names


cat("\nMultinomial scan column names:\n")

print(
  colnames(
    multi_scan
  )
)


############################################################
# 20. Add SNP marker information
############################################################

result <- data.frame(

  SNP = paste0(
    "SNP",
    seq_len(
      nrow(geno)
    )
  ),

  CHROM =
    geno$`#CHROM`,

  POS =
    geno$POS,

  REF =
    geno$REF,

  ALT =
    geno$ALT,

  multi_scan,

  check.names = FALSE
)


############################################################
# 21. Confirm P_global exists
############################################################

if(
  !"P_global" %in%
    names(result)
) {

  stop(
    "P_global column is missing."
  )
}


############################################################
# 22. Global multiple-testing correction
############################################################

result <- result %>%
  mutate(

    FDR_global =
      p.adjust(
        P_global,
        method = "BH"
      ),

    Bonferroni_global =
      p.adjust(
        P_global,
        method = "bonferroni"
      ),

    minus_log10P_global =
      -log10(
        P_global
      )
  )


############################################################
# 23. A vs AE multiple-testing correction
############################################################

result <- result %>%
  mutate(

    FDR_A =
      p.adjust(
        P_A,
        method = "BH"
      ),

    Bonferroni_A =
      p.adjust(
        P_A,
        method = "bonferroni"
      ),

    minus_log10P_A =
      -log10(
        P_A
      )
  )


############################################################
# 24. E vs AE multiple-testing correction
############################################################

result <- result %>%
  mutate(

    FDR_E =
      p.adjust(
        P_E,
        method = "BH"
      ),

    Bonferroni_E =
      p.adjust(
        P_E,
        method = "bonferroni"
      ),

    minus_log10P_E =
      -log10(
        P_E
      )
  )


############################################################
# 25. Null vs AE multiple-testing correction
############################################################

result <- result %>%
  mutate(

    FDR_Null =
      p.adjust(
        P_Null,
        method = "BH"
      ),

    Bonferroni_Null =
      p.adjust(
        P_Null,
        method = "bonferroni"
      ),

    minus_log10P_Null =
      -log10(
        P_Null
      )
  )


############################################################
# 26. Sort by global P-value
############################################################

result_sorted <- result %>%
  arrange(
    P_global
  )


############################################################
# 27. Top 20 global SNPs
############################################################

top20_global <- result_sorted %>%

  filter(
    !is.na(
      P_global
    )
  ) %>%

  slice_head(
    n = 20
  )


cat("\n============================================\n")
cat("Top 20 global SNPs\n")
cat("============================================\n")

print(

  top20_global %>%

    select(

      SNP,
      CHROM,
      POS,
      REF,
      ALT,

      N,
      N0,
      N1,
      N2,

      LRChiSq,
      df_global,
      P_global,
      FDR_global,
      Bonferroni_global
    )
)


############################################################
# 28. Top 20 A vs AE
############################################################

top20_A <- result %>%

  filter(
    !is.na(
      P_A
    )
  ) %>%

  arrange(
    P_A
  ) %>%

  slice_head(
    n = 20
  )


############################################################
# 29. Top 20 E vs AE
############################################################

top20_E <- result %>%

  filter(
    !is.na(
      P_E
    )
  ) %>%

  arrange(
    P_E
  ) %>%

  slice_head(
    n = 20
  )


############################################################
# 30. Top 20 Null vs AE
############################################################

top20_Null <- result %>%

  filter(
    !is.na(
      P_Null
    )
  ) %>%

  arrange(
    P_Null
  ) %>%

  slice_head(
    n = 20
  )


############################################################
# 31. Genome-wide summary
############################################################

n_valid <- sum(
  !is.na(
    result$P_global
  )
)

bonf_threshold <- 0.05 /
  n_valid


cat("\n============================================\n")
cat("Genome-wide summary\n")
cat("============================================\n")


cat(
  "Valid SNP tests:",
  n_valid,
  "\n"
)


cat(
  "Bonferroni threshold:",
  bonf_threshold,
  "\n"
)


cat(
  "Global P < 0.05:",
  sum(
    result$P_global < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Global P < 0.001:",
  sum(
    result$P_global < 0.001,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Global FDR < 0.05:",
  sum(
    result$FDR_global < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Global Bonferroni < 0.05:",
  sum(
    result$Bonferroni_global < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Non-converged SNPs:",
  sum(
    result$Converged == 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "SNPs without valid global test:",
  sum(
    is.na(
      result$P_global
    )
  ),
  "\n"
)


############################################################
# 32. Save Excel
############################################################

write_xlsx(

  list(

    Multinomial_Association =
      result_sorted,

    Top20_Global =
      top20_global,

    Top20_A_vs_AE =
      top20_A,

    Top20_E_vs_AE =
      top20_E,

    Top20_Null_vs_AE =
      top20_Null
  ),

  "results/Pheromone_4Class_Multinomial_Association_v2.xlsx"
)


############################################################
# 33. Save CSV
############################################################

write.csv(

  result_sorted,

  "results/Pheromone_4Class_Multinomial_Association_v2.csv",

  row.names = FALSE
)


############################################################
# 34. Finished
############################################################

cat("\n============================================\n")
cat("Multinomial association completed successfully\n")
cat("============================================\n")

cat("\nOutput files:\n")

cat(
  "results/Pheromone_4Class_Multinomial_Association_v2.xlsx\n"
)

cat(
  "results/Pheromone_4Class_Multinomial_Association_v2.csv\n"
)

