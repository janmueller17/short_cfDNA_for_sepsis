#!/usr/bin/env Rscript

# R script to run ROC analysis included in manuscript of identified top hits from differential enrichment analyses versus clinical markers.

# Provided data files in references have to be extracted, before execution of this script.

# Example usage: Rscript ROC_LF_top_hits_vs_clinical_data.R

# Clear workspace
rm(list = ls())

# Load required libraries
library(readr)
library(dplyr)
library(pROC)

# Define paths and constants
base_path <- "path/to/differential_enrichment/"
out_path <- file.path(base_path, "/ROC/")
top_to_select <- 10
conditions_in_early_death <- factor(c(rep("Recovery", 21), rep("Early death", 22)))
conditions_in_infection_type <- factor(c(rep("Viral", 8), rep("Bacterial", 10)))

# Create out_path if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Load data
cohort_early_death_metadata <- read_csv(file.path(base_path, "/references/cohort_metadata_early_death.csv"), na = "NA")
cohort_infection_type_metadata <- read_csv(file.path(base_path, "/references/cohort_metadata_infection_type.csv"), na = "NA")
cohort_early_death_LF <- read_csv(file.path(base_path, "/references/early_death_edgeR_norm_reads_test_results.csv"))
cohort_infection_type_LF <- read_csv(file.path(base_path, "/references/infection_type_edgeR_norm_reads_test_results.csv"))

# Process early death metadata
cohort_early_death_metadata <- cohort_early_death_metadata %>%
  mutate(Survival = factor(Survival, levels = c("recovery", "early_death"), labels = c("Recovery", "Early death")),
         Gender = factor(Gender, levels = c("m", "f"), labels = c("Male", "Female")),
         Septic_shock = factor(Septic_shock, levels = c("no", "yes"), labels = c("No", "Yes")),
         Vasopressors_inotropes = factor(Vasopressors_inotropes, levels = c("no", "yes"), labels = c("No", "Yes")),
         RRT = factor(RRT, levels = c("no", "yes"), labels = c("No", "Yes")))

# Process infection type metadata
cohort_infection_type_metadata <- cohort_infection_type_metadata %>%
  mutate(Survival = factor(Survival, levels = c("recovery", "early_death"), labels = c("Recovery", "Early death")),
         Infection_type = factor(Infection_type, levels = c("bacterial", "viral"), labels = c("Bacterial", "Viral")),
         Gender = factor(Gender, levels = c("m", "f"), labels = c("Male", "Female")),
         Septic_shock = factor(Septic_shock, levels = c("no", "yes"), labels = c("No", "Yes")),
         Vasopressors_inotropes = factor(Vasopressors_inotropes, levels = c("no", "yes"), labels = c("No", "Yes")),
         RRT = factor(RRT, levels = c("no", "yes"), labels = c("No", "Yes")))

# Select top genes for early death
select_top_genes <- function(df, top_n) {
  df_FCpos <- df %>% filter(logFC > 0) %>% arrange(PValue) %>% head(top_n / 2)
  df_FCneg <- df %>% filter(logFC < 0) %>% arrange(PValue) %>% head(top_n / 2)
  rbind(df_FCpos, df_FCneg)
}

cohort_early_death_LF_top <- select_top_genes(cohort_early_death_LF, top_to_select)
cohort_infection_type_LF_top <- select_top_genes(cohort_infection_type_LF, top_to_select)

# Perform PCA and create data frames
perform_pca <- function(df, conditions) {
  prcomp_result <- prcomp(t(df), center = TRUE, scale. = TRUE)
  data.frame(PC1 = prcomp_result$x[, 1],
             PC2 = prcomp_result$x[, 2],
             PC3 = prcomp_result$x[, 3],
             classification = conditions)
}

cohort_early_death_pca_df <- perform_pca(cohort_early_death_LF_top[, 2:44], conditions_in_early_death)
cohort_infection_type_pca_df <- perform_pca(cohort_infection_type_LF_top[, 2:19], conditions_in_infection_type)

# Plot ROC curves
pdf(file = file.path(out_path, "ROC_LF_top_hits_vs_clinical_data.pdf"), width = 4, height = 4)

roc(response = cohort_early_death_pca_df$classification, predictor = cohort_early_death_pca_df$PC2,
    ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2,
    main = "Early death - PC2, SOFA",
    print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
plot.roc(cohort_early_death_metadata$Survival, cohort_early_death_metadata$SOFA, add = TRUE, ci = TRUE,
         print.auc = TRUE, print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f", lty = 2)

roc(response = cohort_infection_type_pca_df$classification, predictor = cohort_infection_type_pca_df$PC1,
    ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2,
    main = "Infection type - PC1, PCT",
    print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
plot.roc(cohort_infection_type_metadata$Infection_type, cohort_infection_type_metadata$PCT, add = TRUE, ci = TRUE,
         print.auc = TRUE, print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f", lty = 2)

dev.off()