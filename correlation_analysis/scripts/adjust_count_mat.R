#!/usr/bin/env Rscript

# Adjust count matrix for GC content of regions and sequencing depth of samples.

# Load necessary libraries
library(EDASeq)
library(readr)
library(data.table)
library(tools)

# Function to compute TPM normalization
tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  coef <- sum(rpk) / 1e6
  rpk / coef
}

# Main script execution
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Insufficient arguments provided. Usage: adjust_count_mat.R <GC_content_file> <count_matrix_file> <output_path>")
  }

  gc_content_file <- args[1]
  count_matrix_file <- args[2]
  out_path <- args[3]

  # Extract the base name of the count matrix file to use in output file names
  count_matrix_base <- file_path_sans_ext(basename(count_matrix_file))

  # Read input files
  GC_content_df <- read_delim(gc_content_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  count_mat_df <- read_table(count_matrix_file)

  # Check for required columns in GC content file
  required_columns <- c("#1_usercol", "2_usercol", "3_usercol", "5_pct_gc", "12_seq_len")
  missing_columns <- setdiff(required_columns, colnames(GC_content_df))
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the GC content file:", paste(missing_columns, collapse = ", ")))
  }

  # Rename first three columns of count matrix for consistency
  colnames(count_mat_df)[1:3] <- c("chr", "start", "end")
  colnames(count_mat_df) <- gsub(x = colnames(count_mat_df), pattern = "'", replacement = "")

  # Create region and name columns in GC content data frame
  GC_content_df$region <- paste0(GC_content_df$`#1_usercol`, ":", GC_content_df$`2_usercol`, "-", GC_content_df$`3_usercol`)
  GC_content_df$name <- paste0(GC_content_df$region, "_", GC_content_df$`12_seq_len`)

  # Create region column in count matrix data frame
  count_mat_df$region <- paste0(count_mat_df$chr, ":", count_mat_df$start, "-", count_mat_df$end)

  # Filter and sort count matrix to include only regions present in GC content data frame
  count_mat_df <- count_mat_df[count_mat_df$region %in% GC_content_df$region, ]
  count_mat_df <- count_mat_df[order(count_mat_df$region), ]
  GC_content_df <- GC_content_df[order(GC_content_df$region), ]

  # Prepare count matrix and feature data for EDASeq
  count_mat_arr <- as.matrix(subset(count_mat_df, select = -c(1:3, ncol(count_mat_df))))
  rownames(count_mat_arr) <- GC_content_df$name

  gc_col_name <- "5_pct_gc"
  feature <- AnnotatedDataFrame(data.frame(gc = GC_content_df[[gc_col_name]], length = GC_content_df$`12_seq_len`, row.names = GC_content_df$name))
  data <- newSeqExpressionSet(counts = count_mat_arr, featureData = feature)

  # Add a single level factor variable to the pData for biasPlot
  pData(data)$single_level_factor <- factor("all")

  # Perform normalization
  dataWithin <- withinLaneNormalization(data, "gc", which = "full", round = FALSE, offset = FALSE) #"loess","median","upper","full"
  dataNorm <- dataWithin # Copy object
  normCounts(dataNorm) <- apply(normCounts(dataWithin), 2, function(x) tpm(x, featureData(dataWithin)$length)) # TPM normalization

  # Generate QC plots
  generate_qc_plots(data, dataWithin, dataNorm, out_path, count_matrix_base)

  # Save normalized counts to CSV files
  save_normalized_counts(data, dataWithin, dataNorm, out_path, count_matrix_base)
}

# Function to generate QC plots
generate_qc_plots <- function(data, dataWithin, dataNorm, out_path, count_matrix_base) {
  pdf(paste0(out_path, "/", count_matrix_base, "_read_QC_plots.pdf"), width = 5, height = 4)
  boxplot(counts(data), main = "Raw counts")
  biasPlot(data, "gc", "single_level_factor", log = FALSE, ylim = c(0, max(counts(data))), main = "Raw counts")
  dev.off()

  pdf(paste0(out_path, "/", count_matrix_base, "_adjusted_read_QC_plots.pdf"), width = 5, height = 4)
  boxplot(normCounts(dataWithin), main = "GC adjusted counts")
  boxplot(normCounts(dataNorm), main = "GC and seq depth adjusted counts")
  biasPlot(dataWithin, "gc", "single_level_factor", log = FALSE, ylim = c(0, max(normCounts(dataWithin))), main = "GC adjusted counts")
  biasPlot(dataNorm, "gc", "single_level_factor", log = FALSE, ylim = c(0, max(normCounts(dataNorm))), main = "GC and Seq-depth adjusted counts")
  dev.off()
}

# Function to save normalized counts to CSV files
save_normalized_counts <- function(data, dataWithin, dataNorm, out_path, count_matrix_base) {
  write.csv(x = counts(data), file = paste0(out_path, "/", count_matrix_base, "_raw_counts.csv"), quote = FALSE, row.names = TRUE)
  write.csv(x = normCounts(dataWithin), file = paste0(out_path, "/", count_matrix_base, "_GCadj_counts.csv"), quote = FALSE, row.names = TRUE)
  write.csv(x = normCounts(dataNorm), file = paste0(out_path, "/", count_matrix_base, "_GCadj_DEPTHadj_counts.csv"), quote = FALSE, row.names = TRUE)
}

# Run the main script
main()