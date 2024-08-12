#!/usr/bin/env bash

# Shell script to run AME analysis on two BED files and generate plots.

# Non-standard tools required (not including respective dependencies):
# BEDTools - https://bedtools.readthedocs.io/
# AME - https://meme-suite.org/meme/doc/ame.html
# R with ggpubr, readr, dplyr, tidyverse, scales, optparse
# Python with sys

# Example usage: ./run_ame_motif_enrichment.sh /path/to/bed1.bed /path/to/bed2.bed -r /path/to/reference.fa -m /path/to/motif_db.meme -e 24 -o /output/dir -s 0.05 -t 10 -n1 condition1 -n2 condition2 -cmp comparison_name -c1 "#CC79A7" -c2 "#E69F00" -a /path/to/ame -b /path/to/bedtools]

# Function to display usage message
usage() {
    echo "Usage: $0 <consensus_peaks_bed_1> <consensus_peaks_bed_2> -r <reference_fasta> -m <motif_db> -e <extension_length> -o <output_dir_base> [-a <ame_executable>] [-b <bedtools_executable>] [-s <significance_threshold>] [-t <top_to_plot>] -n1 <condition_name1> -n2 <condition_name2> -c1 <condition_color1> -c2 <condition_color2> -cmp <comparison_name>"
    echo "  -r, --reference           Path to the reference FASTA file (required)"
    echo "  -m, --motif-database      Path to the motif database file (required)"
    echo "  -e, --extension-length    Length to extend BED files (required)"
    echo "  -o, --output-dir          Base output directory for results (required)"
    echo "  -a, --ame-executable      Path to the AME executable (default: ame; assumes AME is installed and in PATH)"
    echo "  -b, --bedtools-executable Path to the bedtools executable (default: bedtools; assumes BEDTools is installed and in PATH)"
    echo "  -s, --significance        Significance threshold for plotting (default: 0.05)"
    echo "  -t, --top                 Number of top motifs to plot (default: 10)"
    echo "  -n1, --condition-name1    Name of condition 1 (required)"
    echo "  -n2, --condition-name2    Name of condition 2 (required)"
    echo "  -c1, --condition-color1   Color for condition 1 in plots (required)"
    echo "  -c2, --condition-color2   Color for condition 2 in plots (required)"
    echo "  -cmp, --comparison-name   Name of the comparison (required)"
    exit 1
}

# Default values
ame_executable="ame"
bedtools_executable="bedtools"
significance_threshold=0.05
top_to_plot=10

# Parse CLI arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -r|--reference) reference_fasta="$2"; shift ;;
        -m|--motif-database) motif_db="$2"; shift ;;
        -e|--extension-length) extension_length="$2"; shift ;;
        -o|--output-dir) output_dir_base="$2"; shift ;;
        -a|--ame-executable) ame_executable="$2"; shift ;;
        -b|--bedtools-executable) bedtools_executable="$2"; shift ;;
        -s|--significance) significance_threshold="$2"; shift ;;
        -t|--top) top_to_plot="$2"; shift ;;
        -n1|--condition-name1) condition_name1="$2"; shift ;;
        -n2|--condition-name2) condition_name2="$2"; shift ;;
        -c1|--condition-color1) condition_color1="$2"; shift ;;
        -c2|--condition-color2) condition_color2="$2"; shift ;;
        -cmp|--comparison-name) comparison_name="$2"; shift ;;
        -h|--help) usage ;;
        *) [ -z "$consensus_peaks_bed_1" ] && consensus_peaks_bed_1="$1" || consensus_peaks_bed_2="$1" ;;
    esac
    shift
done

# Check for required arguments
if [ -z "${consensus_peaks_bed_1}" ] || [ -z "${consensus_peaks_bed_2}" ] || [ -z "${reference_fasta}" ] || [ -z "${motif_db}" ] || [ -z "${extension_length}" ] || [ -z "${output_dir_base}" ] || [ -z "${condition_name1}" ] || [ -z "${condition_name2}" ] || [ -z "${condition_color1}" ] || [ -z "${condition_color2}" ] || [ -z "${comparison_name}" ]; then
    usage
fi

# Check if files exist
if [ ! -f "${consensus_peaks_bed_1}" ]; then echo "Error: ${consensus_peaks_bed_1} does not exist."; exit 1; fi
if [ ! -f "${consensus_peaks_bed_2}" ]; then echo "Error: ${consensus_peaks_bed_2} does not exist."; exit 1; fi
if [ ! -f "${reference_fasta}" ]; then echo "Error: ${reference_fasta} does not exist."; exit 1; fi
if [ ! -f "${motif_db}" ]; then echo "Error: ${motif_db} does not exist."; exit 1; fi

# Define common paths
script_dir="$(dirname "$0")"
extend_bed_script="${script_dir}/scripts/extend_bed.py"
plot_script="${script_dir}/scripts/plot_ame_results.R"

# Extend BED files
echo "Extending BED files..."
python "${extend_bed_script}" "${consensus_peaks_bed_1}" "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.bed" "${extension_length}"
python "${extend_bed_script}" "${consensus_peaks_bed_2}" "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.bed" "${extension_length}"

# Check if the extended BED files were created
if [ ! -f "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.bed" ]; then echo "Error: Failed to create extended BED file for ${consensus_peaks_bed_1}."; exit 1; fi
if [ ! -f "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.bed" ]; then echo "Error: Failed to create extended BED file for ${consensus_peaks_bed_2}."; exit 1; fi

# Generate FASTA files
echo "Generating FASTA files..."
${bedtools_executable} getfasta -fi "$reference_fasta" -bed "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.bed" -fo "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.fasta"
${bedtools_executable} getfasta -fi "$reference_fasta" -bed "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.bed" -fo "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.fasta"

# Check if the FASTA files were created
if [ ! -f "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.fasta" ]; then echo "Error: Failed to create FASTA file for ${consensus_peaks_bed_1}."; exit 1; fi
if [ ! -f "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.fasta" ]; then echo "Error: Failed to create FASTA file for ${consensus_peaks_bed_2}."; exit 1; fi

# Function to create output directory
create_output_dir() {
    local output_subdir="$1"
    if [ ! -d "${output_dir_base}/${output_subdir}" ]; then
        mkdir -p "${output_dir_base}/${output_subdir}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create output directory ${output_dir_base}/${output_subdir}."
            exit 1
        fi
    fi
}

# Create output directories
create_output_dir "${comparison_name}_${condition_name1}"
create_output_dir "${comparison_name}_${condition_name2}"

# Run AME analysis
run_ame() {
    local control_file="$1"
    local test_file="$2"
    local output_subdir="$3"

    echo "Running AME analysis for ${control_file} vs ${test_file}..."
    "${ame_executable}" --control "${control_file}" --scoring avg --method fisher --noseq --evalue-report-threshold 100000 --oc "${output_dir_base}/${output_subdir}" "${test_file}" "${motif_db}"
}

run_ame "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.fasta" "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.fasta" "${comparison_name}_${condition_name1}"
run_ame "${consensus_peaks_bed_1%.bed}_extended_${extension_length}.fasta" "${consensus_peaks_bed_2%.bed}_extended_${extension_length}.fasta" "${comparison_name}_${condition_name2}"

# Check if the AME output files were created
if [ ! -f "${output_dir_base}/${comparison_name}_${condition_name1}/ame.tsv" ]; then echo "Error: AME output file for ${condition_name1} not found."; exit 1; fi
if [ ! -f "${output_dir_base}/${comparison_name}_${condition_name2}/ame.tsv" ]; then echo "Error: AME output file for ${condition_name2} not found."; exit 1; fi

# Plot results
echo "Generating plots..."
Rscript "$plot_script" \
  --data "${output_dir_base}/${comparison_name}_${condition_name1}/ame.tsv" \
  --condition "$condition_name1" \
  --output "${output_dir_base}/${comparison_name}_${condition_name1}/barchart.pdf" \
  --top "$top_to_plot" \
  --significance "$significance_threshold" \
  --color "$condition_color1"

Rscript "$plot_script" \
  --data "${output_dir_base}/${comparison_name}_${condition_name2}/ame.tsv" \
  --condition "$condition_name2" \
  --output "${output_dir_base}/${comparison_name}_${condition_name2}/barchart.pdf" \
  --top "$top_to_plot" \
  --significance "$significance_threshold" \
  --color "$condition_color2"

# Check if the plot files were created
if [ ! -f "${output_dir_base}/${comparison_name}_${condition_name1}/barchart.pdf" ]; then echo "Error: Plot for ${condition_name1} not created."; exit 1; fi
if [ ! -f "${output_dir_base}/${comparison_name}_${condition_name2}/barchart.pdf" ]; then echo "Error: Plot for ${condition_name2} not created."; exit 1; fi

echo "Analysis and plotting completed successfully."