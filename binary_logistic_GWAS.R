library(readxl)
library(dplyr)
library(writexl)

cat("Starting association analysis...\n")
dir.create("results", recursive = TRUE, showWarnings = FALSE)

############################################################
# 1. Read data
############################################################

pheno_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/1.Phenotype.xlsx"
geno_file  <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/2.SNP_GENOTYPE1224_012.xlsx"

pheno <- read_excel(pheno_file, sheet = "Sheet1")
geno  <- read_excel(geno_file, sheet = "Sheet2")

cat("Phenotype:", nrow(pheno), "individuals\n")
cat("Genotype:", nrow(geno), "SNPs\n")

############################################################
# 2. Define binary phenotypes
############################################################

pheno <- pheno %>%
  mutate(
    Ald = ifelse(`A-Type` == "A-", 1, 0),
    Ester = ifelse(`E-Type` == "E-", 1, 0)
  )

cat("\nAld phenotype:\n")
print(table(pheno$Ald))

cat("\nEster phenotype:\n")
print(table(pheno$Ester))

############################################################
# 3. Match IDs
############################################################

sample_cols <- colnames(geno)[colnames(geno) %in% pheno$ID]

cat("\nMatched individuals:", length(sample_cols), "\n")

if(length(sample_cols) == 0){
  stop("No matching IDs between phenotype and genotype files.")
}

cat("\nPhenotype IDs missing from genotype:\n")
print(setdiff(pheno$ID, sample_cols))

cat("\nGenotype IDs missing from phenotype:\n")
print(setdiff(sample_cols, pheno$ID))

if(anyDuplicated(pheno$ID)){
  stop("Duplicated phenotype IDs.")
}

if(anyDuplicated(sample_cols)){
  stop("Duplicated genotype IDs.")
}

pheno2 <- pheno[
  match(sample_cols, pheno$ID),
]

stopifnot(all(pheno2$ID == sample_cols))

cat("ID matching successful.\n")

############################################################
# 4. Genotype matrix
############################################################

gmat <- as.matrix(geno[, sample_cols])

storage.mode(gmat) <- "numeric"

cat(
  "Genotype matrix:",
  nrow(gmat), "SNP x",
  ncol(gmat), "individuals\n"
)

cat("\nGenotype counts:\n")
print(table(gmat, useNA = "ifany"))

############################################################
# 5. Single-SNP additive logistic regression
#
# logit(P) = beta0 + beta1 * G
# G = 0,1,2
############################################################

logistic_snp <- function(g, y) {

  ok <- !is.na(g) & !is.na(y)

  g <- g[ok]
  y <- y[ok]

  N  <- length(y)
  N0 <- sum(g == 0)
  N1 <- sum(g == 1)
  N2 <- sum(g == 2)

  # SNP must contain variation
  if(length(unique(g)) < 2) {
    return(c(
      N = N,
      N0 = N0,
      N1 = N1,
      N2 = N2,
      Beta = NA,
      SE = NA,
      OR = NA,
      WaldChiSq = NA,
      P = NA
    ))
  }

  fit <- try(
    suppressWarnings(
      glm(
        y ~ g,
        family = binomial(link = "logit")
      )
    ),
    silent = TRUE
  )

  if(inherits(fit, "try-error")) {
    return(c(
      N = N,
      N0 = N0,
      N1 = N1,
      N2 = N2,
      Beta = NA,
      SE = NA,
      OR = NA,
      WaldChiSq = NA,
      P = NA
    ))
  }

  sm <- summary(fit)$coefficients

  if(!"g" %in% rownames(sm)) {
    return(c(
      N = N,
      N0 = N0,
      N1 = N1,
      N2 = N2,
      Beta = NA,
      SE = NA,
      OR = NA,
      WaldChiSq = NA,
      P = NA
    ))
  }

  beta <- sm["g", "Estimate"]
  se   <- sm["g", "Std. Error"]

  if(!is.finite(beta) || !is.finite(se) || se <= 0) {
    return(c(
      N = N,
      N0 = N0,
      N1 = N1,
      N2 = N2,
      Beta = NA,
      SE = NA,
      OR = NA,
      WaldChiSq = NA,
      P = NA
    ))
  }

  OR <- exp(beta)

  WaldChiSq <- (beta / se)^2

  P <- pchisq(
    WaldChiSq,
    df = 1,
    lower.tail = FALSE
  )

  c(
    N = N,
    N0 = N0,
    N1 = N1,
    N2 = N2,
    Beta = beta,
    SE = se,
    OR = OR,
    WaldChiSq = WaldChiSq,
    P = P
  )
}

############################################################
# 6. Ald scan
############################################################

cat("\nRunning Ald scan...\n")
t1 <- Sys.time()

ald_scan <- t(
  apply(
    gmat,
    1,
    logistic_snp,
    y = pheno2$Ald
  )
)

ald_scan <- as.data.frame(ald_scan)

cat("Ald finished:", as.character(Sys.time() - t1), "\n")

############################################################
# 7. Ester scan
############################################################

cat("\nRunning Ester scan...\n")
t1 <- Sys.time()

ester_scan <- t(
  apply(
    gmat,
    1,
    logistic_snp,
    y = pheno2$Ester
  )
)

ester_scan <- as.data.frame(ester_scan)

cat("Ester finished:", as.character(Sys.time() - t1), "\n")

############################################################
# 8. Add marker information
############################################################

snp_id <- paste0("SNP", seq_len(nrow(geno)))

ald_result <- data.frame(
  SNP = snp_id,
  CHROM = geno$`#CHROM`,
  POS = geno$POS,
  REF = geno$REF,
  ALT = geno$ALT,
  ald_scan
)

ester_result <- data.frame(
  SNP = snp_id,
  CHROM = geno$`#CHROM`,
  POS = geno$POS,
  REF = geno$REF,
  ALT = geno$ALT,
  ester_scan
)

############################################################
# 9. Multiple-testing correction
############################################################

ald_result <- ald_result %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Bonferroni = p.adjust(P, method = "bonferroni"),
    minus_log10P = -log10(P)
  ) %>%
  arrange(P)

ester_result <- ester_result %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Bonferroni = p.adjust(P, method = "bonferroni"),
    minus_log10P = -log10(P)
  ) %>%
  arrange(P)

############################################################
# 10. Output
############################################################

write_xlsx(
  list(
    Ald_Association = ald_result,
    Ester_Association = ester_result,
    Top20_Ald = head(ald_result, 20),
    Top20_Ester = head(ester_result, 20)
  ),
  "results/Pheromone_Logistic_Wald_Association.xlsx"
)

write.csv(
  ald_result,
  "results/Ald_Association.csv",
  row.names = FALSE
)

write.csv(
  ester_result,
  "results/Ester_Association.csv",
  row.names = FALSE
)

cat("\nTop Ald SNPs:\n")
print(head(ald_result, 10))

cat("\nTop Ester SNPs:\n")
print(head(ester_result, 10))

cat("\nAnalysis completed successfully.\n")

