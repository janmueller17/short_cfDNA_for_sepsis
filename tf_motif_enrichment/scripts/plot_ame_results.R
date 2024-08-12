#!/usr/bin/env Rscript

# Create a bar plot of top motifs from AME analysis

# Load necessary libraries
library(ggpubr)
library(readr)
library(dplyr)
library(tidyverse)
library(scales)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-d", "--data"), type = "character", default = NULL, help = "Input AME data file path", metavar = "character"),
  make_option(c("-c", "--condition"), type = "character", default = "Condition", help = "Condition name", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output PDF file path", metavar = "character"),
  make_option(c("-t", "--top"), type = "integer", default = 10, help = "Number of top motifs to display", metavar = "integer"),
  make_option(c("-s", "--significance"), type = "numeric", default = 0.05, help = "Significance threshold", metavar = "numeric"),
  make_option(c("-l", "--color"), type = "character", default = "#000000", help = "Condition color", metavar = "character")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Read input data
process_tsv <- function(file_path) {
  data <- read_tsv(file_path, quote = "#") %>%
    head(-1) %>%
    mutate(
      motif_alt_ID = str_split_fixed(motif_ID, "_", 2)[, 1],
      logFC = log2(((TP + 0.1) / pos) / ((FP + 0.1) / neg)),
      log_adj_p = -log10(`adj_p-value`)
    )
  return(data)
}

# Create plot
create_plot <- function(data, condition_name, condition_color, sign_thresh, top_select, out_path) {
  data_top <- data %>%
    arrange(desc(log_adj_p)) %>%
    slice_head(n = top_select) %>%
    mutate(color = ifelse(log_adj_p > -log10(sign_thresh), "Significant", "Not Significant")) %>%
    mutate(color = factor(color, levels = c("Significant", "Not Significant")))

  p <- ggbarplot(
    data = data_top,
    x = "motif_alt_ID",
    y = "log_adj_p",
    orientation = "horiz",
    order = data_top$motif_alt_ID[order(data_top$log_adj_p, decreasing = FALSE)],
    fill = "color",
    xlab = "Transcription factor",
    title = condition_name) +
    ylab(expression(-log[10](adj.p))) +
    scale_fill_manual(
      values = c("Significant" = condition_color, "Not Significant" = alpha(condition_color, 0.5)),
      name = "Significance")

  ggsave(out_path, plot = p, width = 3, height = 5)
}

# Main execution
data_in <- process_tsv(opt$data)
create_plot(data_in, opt$condition, opt$color, opt$significance, opt$top, opt$output)