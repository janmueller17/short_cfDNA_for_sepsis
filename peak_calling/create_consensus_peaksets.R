#!/usr/bin/env Rscript

# Example usage
# ./create_consensus_peaksets.R --input <sample_001_peaks.narrowPeak>,<sample_002_peaks.narrowPeak>,<sample_003_peaks.narrowPeak> --output <consensus_peaks_narrow.bed> --merging_gap_size 31 --peak_percentage 0.5

# Load necessary libraries
library(rtracklayer)
library(GenomicRanges)
library(BRGenomics)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Comma-separated list of input peak files (.narrowPeak)"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output file path (.bed)"),
  make_option(c("-m", "--merging_gap_size"), type = "integer", default = 31, help = "Merging gap size [default: 31]"),
  make_option(c("-p", "--peak_percentage"), type = "numeric", default = 0.5, help = "Percentage of input files required to have a peak [default: 0.5]")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check mandatory arguments
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input files and output file must be provided.", call. = FALSE)
}

# Convert input from comma-separated string to a vector
peak_input_files <- unlist(strsplit(opt$input, ","))

# Check if input files exist
missing_files <- peak_input_files[!file.exists(peak_input_files)]
if (length(missing_files) > 0) {
  stop(paste("The following input files do not exist:", paste(missing_files, collapse = ", ")), call. = FALSE)
}

# Calculate min_samples_with_peak as a percentage of the number of input files
min_samples_with_peak <- round(length(peak_input_files) * opt$peak_percentage)

# Import peak files as GRanges objects
peak_granges <- lapply(peak_input_files, import)
peak_grangeslist <- GRangesList(peak_granges)

# Find genome regions covered by at least min_samples_with_peak sets of peaks
peak_coverage <- coverage(peak_grangeslist)
covered_ranges <- slice(peak_coverage, lower = min_samples_with_peak, rangesOnly = TRUE)
covered_granges <- GRanges(covered_ranges)
covered_granges <- reduce(covered_granges, min.gapwidth = opt$merging_gap_size)

# Keep only default chromosomes
covered_granges <- tidyChromosomes(covered_granges, keep.X = FALSE, keep.Y = FALSE)

# Add names to the regions
names(covered_granges) <- paste0("region_", seq_along(covered_granges))

# Export the regions as a .bed file
export(covered_granges, opt$output)