# Genetic Architecture and Behavioral Consequences of Sex Pheromone Polymorphism in *Mythimna separata*

This repository contains the R scripts used for chemical phenotype analysis, segregation analysis, genome-wide association analyses, lead-SNP genotype–phenotype summaries, and field-trapping analyses in *Mythimna separata*.

## Scripts

- `Chemical phenotype comparison, PCA and segregation analysis.R`: chemical phenotype comparisons, principal component analysis, and chi-square segregation tests.
- `binary_logistic_GWAS.R`: additive logistic-regression GWAS for binary pheromone-production traits.
- `four_class_multinomial_GWAS.R`: multinomial GWAS for the four-category pheromone phenotype.
- `lead_SNP_Fisher_stacked.R`: lead-SNP genotype–phenotype association tests and proportional stacked bar plots.
- `trial_total_heatmap.R`: trial-level total-catch heatmap for field trapping.
- `negative_binomial_GLMM_boxplot_forest.R`: trial-level standardized catches and negative-binomial GLMM analyses of field-trapping data.

## Software

The scripts use R and the following packages where applicable: `readxl`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `cowplot`, `writexl`, `nnet`, `glmmTMB`, `emmeans`, `scales`, `broom.mixed`, `performance`, and `DHARMa`.

## Data paths

Input and output paths are defined near the beginning of each script. Users should update these paths to match their local data directory before running the analyses.

## Citation

Please cite the associated article, **“Genetic Architecture and Behavioral Consequences of Sex Pheromone Polymorphism in *Mythimna separata*”**, when using this code.
