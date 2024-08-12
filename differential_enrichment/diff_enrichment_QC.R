#!/usr/bin/env Rscript

# Run QC for differential enrichment analysis of short cfDNA data (before actual analysis). Used to define thresholds for filtering data -> staged execution.

# Example usage ./diff_enrichment_QC.R -b "/path/to/base" -o "/path/to/output" -n infection_type -m 0.1 -z 8 -c "Viral,Viral,Viral,Viral,Viral,Viral,Viral,Viral,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial"


# Load required libraries
library(EDASeq)
library(readr)
library(corrplot)
library(ggplot2)
library(pals)
library(MASS)
library(optparse)

# Command line argument parsing
option_list <- list(
  make_option(c("-b", "--base_path"), type = "character", default = NULL,
              help = "Base path for the project", metavar = "character"),
  make_option(c("-o", "--out_path"), type = "character", default = NULL,
              help = "Output path for the results", metavar = "character"),
  make_option(c("-n", "--name"), type = "character", default = NULL,
              help = "Name for the analysis", metavar = "character"),
  make_option(c("-c", "--conditions"), type = "character", default = NULL,
              help = "Conditions for the samples (comma separated)", metavar = "character"),
  make_option(c("-m", "--mean_thresh"), type = "numeric", default = 0.1,
              help = "Mean reads per region threshold", metavar = "numeric"),
  make_option(c("-z", "--zeros_thresh"), type = "integer", default = NULL,
              help = "Threshold for number of zeros per condition", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$base_path) || is.null(opt$out_path) || is.null(opt$name) || is.null(opt$conditions)) {
  print_help(opt_parser)
  stop("All arguments must be supplied (base_path, out_path, name, conditions).", call. = FALSE)
}

# Parse conditions
conditions <- unlist(strsplit(opt$conditions, ","))
conditions <- factor(conditions)

# Define paths
base_path <- opt$base_path
out_path <- opt$out_path
name <- opt$name

mean_reads_per_region_thresh <- opt$mean_thresh
zeros_per_cond <- opt$zeros_thresh

# Function to check if a file exists
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
}

# Function to read and prepare data
read_and_prepare_data <- function(base_path, name) {
  gc_content_path <- file.path(base_path, paste0(name, "_GC.bed"))
  count_matrix_path <- file.path(base_path, paste0(name, ".txt"))

  check_file_exists(gc_content_path)
  check_file_exists(count_matrix_path)

  gc_content_df <- read_delim(gc_content_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  count_matrix_df <- read_table(count_matrix_path)

  gc_content_df$peak <- paste0(gc_content_df$`#1_usercol`, ":", gc_content_df$`2_usercol`, "-", gc_content_df$`3_usercol`)
  gc_content_df$length <- abs(gc_content_df$`2_usercol` - gc_content_df$`3_usercol`)

  count_matrix_df <- count_matrix_df[order(count_matrix_df$peak), ]
  gc_content_df <- gc_content_df[order(gc_content_df$peak), ]

  list(gc_content_df = gc_content_df, count_matrix_df = count_matrix_df)
}

# Function to create Annotated DataFrame
create_annotated_dataframe <- function(count_matrix_df, conditions) {
  pheno_df <- AnnotatedDataFrame(data.frame(conditions = conditions))
  rownames(pheno_df) <- colnames(count_matrix_df[, -1])
  pheno_df
}

# Function to filter data
filter_data <- function(count_matrix_arr, gc_content_df, mean_reads_per_region_thresh, zeros_per_cond, conditions) {
  mean_filter <- rowMeans(count_matrix_arr) > mean_reads_per_region_thresh
  count_matrix_arr_filt_temp <- count_matrix_arr[mean_filter, ]
  gc_content_df_filt_temp <- gc_content_df[mean_filter, ]

  design_matrix <- model.matrix(~conditions, data = pData(pheno_df))
  zero_filter <- (rowSums(count_matrix_arr_filt_temp[, design_matrix[, 2] == 0] == 0) <= zeros_per_cond) &
                (rowSums(count_matrix_arr_filt_temp[, design_matrix[, 2] == 1] == 0) <= zeros_per_cond)
  count_matrix_arr_filt <- count_matrix_arr_filt_temp[zero_filter, ]
  gc_content_df_filt <- gc_content_df_filt_temp[zero_filter, ]

  list(count_matrix_arr_filt = count_matrix_arr_filt, gc_content_df_filt = gc_content_df_filt)
}

# Function to create SeqExpressionSet
create_seqexpression_set <- function(count_matrix_arr_filt, gc_content_df_filt, pheno_df) {
  feature_df <- AnnotatedDataFrame(data.frame(gc = gc_content_df_filt$`5_pct_gc`, length = gc_content_df_filt$length, row.names = gc_content_df_filt$peak))
  newSeqExpressionSet(counts = count_matrix_arr_filt, featureData = feature_df, phenoData = pheno_df)
}

# Function for QC plots
qc_plots <- function(data_seq, data_normalized, out_path, name) {
  pdf(file = file.path(out_path, paste0(name, "_avg_read_distro_after_norm.pdf")), width = 5, height = 4)
  hist(log2(rowMeans(normCounts(data_normalized))), freq = TRUE, breaks = 50, xlim = c(0, 8))
  dev.off()

  color_palette <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

  pdf(file.path(out_path, paste0(name, "_read_distros.pdf")), width = 5, height = 4)
  boxplot(data_seq)
  biasPlot(data_seq, "length", log = TRUE, ylim = c(0, log(max(counts(data_normalized)))))
  biasPlot(data_seq, "gc", log = TRUE, ylim = c(0, log(max(counts(data_normalized)))))
  correlation_matrix <- cor(counts(data_seq), method = "pearson")
  corrplot(correlation_matrix, method = "color", col = color_palette(200), type = "upper", tl.col = "black", diag = FALSE)
  dev.off()

  pdf(file.path(out_path, paste0(name, "_norm_read_distros.pdf")), width = 5, height = 4)
  boxplot(data_within_lane)
  boxplot(data_normalized)
  biasPlot(data_normalized, "length", log = TRUE, ylim = c(0, log(max(counts(data_normalized)))))
  biasPlot(data_normalized, "gc", log = TRUE, ylim = c(0, log(max(counts(data_normalized)))))
  correlation_matrix <- cor(normCounts(data_normalized), method = "pearson")
  corrplot(correlation_matrix, method = "color", col = color_palette(200), type = "upper", tl.col = "black", diag = FALSE)
  dev.off()
}

# Function for PCA analysis
pca_analysis <- function(data_normalized, out_path, name, conditions) {
  pca_input <- as.data.frame(t(normCounts(data_normalized)))
  constant_columns <- names(pca_input)[apply(pca_input, 2, function(column) var(column) == 0)]
  pca_input_filt <- pca_input[, !(names(pca_input) %in% constant_columns)]

  prcomp_result <- prcomp(pca_input_filt, center = TRUE, scale. = TRUE)
  pca_df <- data.frame(
    PC1 = prcomp_result$x[, 1],
    PC2 = prcomp_result$x[, 2],
    classification = conditions
  )
  variance_explained <- (prcomp_result$sdev)^2 / sum((prcomp_result$sdev)^2)

  pdf(file = file.path(out_path, paste0(name, "_norm_read_PCA.pdf")), width = 5, height = 3)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = classification)) +
    geom_point() +
    scale_color_manual(values = glasbey(n = length(unique(pca_df$classification))), name = "Experiment_ID") +
    xlab(paste0("PC1(", round(variance_explained[1] * 100, 2), "%)")) +
    ylab(paste0("PC2(", round(variance_explained[2] * 100, 2), "%)")) +
    theme_bw()
  print(p)
  dev.off()

  prcomp_result <- prcomp(as.data.frame(t(counts(data_normalized))), center = TRUE, scale. = TRUE)
  pca_df <- data.frame(
    PC1 = prcomp_result$x[, 1],
    PC2 = prcomp_result$x[, 2],
    classification = conditions
  )
  variance_explained <- (prcomp_result$sdev)^2 / sum((prcomp_result$sdev)^2)

  pdf(file = file.path(out_path, paste0(name, "_PCA.pdf")), width = 5, height = 3)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = classification)) +
    geom_point() +
    scale_color_manual(values = glasbey(n = length(unique(pca_df$classification))), name = "Experiment_ID") +
    xlab(paste0("PC1(", round(variance_explained[1] * 100, 2), "%)")) +
    ylab(paste0("PC2(", round(variance_explained[2] * 100, 2), "%)")) +
    theme_bw()
  print(p)
  dev.off()
}

# Main script
data <- read_and_prepare_data(base_path, name)
pheno_df <- create_annotated_dataframe(data$count_matrix_df, conditions)
count_matrix_arr <- as.matrix(data$count_matrix_df[, -1])
rownames(count_matrix_arr) <- data$count_matrix_df$peak

filtered_data <- filter_data(count_matrix_arr, data$gc_content_df, mean_reads_per_region_thresh, zeros_per_cond, conditions)
seq_set <- create_seqexpression_set(filtered_data$count_matrix_arr_filt, filtered_data$gc_content_df_filt, pheno_df)

data_within_lane <- withinLaneNormalization(seq_set, "gc", which = "full")
data_normalized <- betweenLaneNormalization(data_within_lane, which = "full")

qc_plots(seq_set, data_normalized, out_path, name)
pca_analysis(data_normalized, out_path, name, conditions)
