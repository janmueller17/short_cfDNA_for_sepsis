#!/usr/bin/env bash

# Annotate, count, and plot proportions for entries in a BED file. Annotations from references directory have to be decompressed before running this script.

# Non-standard tools required (not including respective dependencies):
# BEDTools - https://bedtools.readthedocs.io/
# Python with matplotlib, pandas, os, sys

# Usage:
# ./annotate_bed.sh -i /path/to/consensus_set.bed -o /path/to/output_directory [-a /path/to/bedtools]

# Default values
bedtools_executable="bedtools"
script_dir="$(dirname "$0")"
plot_script_path="${script_dir}/scripts/plot_piecharts_annotation.py"

annotations=(
    "hg19-cCREs.PLS.sort.bed"
    "hg19-cCREs.pELS.sort.bed"
    "hg19-cCREs.dELS.sort.bed"
    "hg19-cCREs.CTCF-only.sort.bed"
    "hg19-cCREs.DNase-H3K4me3.sort.bed"
    "cpgIslandExt.hg19.chrnr.bed"
    "Txn_Factr_ChIP_E3.bed"
    "gene_bodies.bed"
)

# Function to display usage message
usage() {
    echo "Usage: $0 -i <input_bed> -o <output_dir> [-a <bedtools_executable>]"
    echo "  -i, --input               Path to the input BED file (required)"
    echo "  -o, --output              Path to the output directory (required)"
    echo "  -a, --bedtools-executable Path to the bedtools executable (default: bedtools; assumes BEDTools is installed and in PATH)"
    exit 1
}

# Parse CLI arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input_bed="$2"; shift ;;
        -o|--output) output_dir="$2"; shift ;;
        -a|--bedtools-executable) bedtools_executable="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check for required arguments
if [ -z "${input_bed}" ] || [ -z "${output_dir}" ]; then
    usage
fi

# Check if input BED file exists
if [ ! -f "${input_bed}" ]; then
    echo "Error: Input BED file ${input_bed} does not exist."
    exit 1
fi

# Check if bedtools executable exists
if [ ! -x "$(command -v ${bedtools_executable})" ]; then
    echo "Error: bedtools executable ${bedtools_executable} not found or not executable."
    exit 1
fi

# Check if plot script exists
if [ ! -f "${plot_script_path}" ]; then
    echo "Error: Plot script ${plot_script_path} does not exist."
    exit 1
fi

# Check if annotation files exist
for annotation in "${annotations[@]}"; do
    if [ ! -f "${script_dir}/references/${annotation}" ]; then
        echo "Error: Annotation file ${script_dir}/references/${annotation} does not exist."
        exit 1
    fi
done

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Generate output filenames based on input filename
base_name=$(basename "${input_bed}" .bed)
annotated_bed_input="${output_dir}/${base_name}_annotated.bed"
output_file="${output_dir}/${base_name}_annotation_counts.txt"

# Annotate BED file
echo "Annotating BED file..."
if ! ${bedtools_executable} intersect -a "${input_bed}" \
    -b "${annotations[@]/#/${script_dir}\/references\/}" \
    -wao > "${annotated_bed_input}"; then
    echo "Error: Failed to annotate BED file."
    exit 1
fi

# Check if the annotated BED file was created
if [ ! -f "${annotated_bed_input}" ]; then
    echo "Error: Failed to create annotated BED file ${annotated_bed_input}."
    exit 1
fi

# Count total lines in the consensus bed file
consensus_total=$(wc -l < "${input_bed}")
if [ $? -ne 0 ]; then
    echo "Error: Failed to count lines in input BED file."
    exit 1
fi

# Function to count unique entries based on the first three columns
count_unique_entries() {
    awk '!seen[$1, $2, $3]++' "$1" | wc -l
}

# Count lines matching pattern for non-annotated entries
consensus_no_anno=$(grep "\.\s*\.\s*-1\s*-1\s*\.\s*0" "${annotated_bed_input}" | count_unique_entries /dev/stdin)
if [ $? -ne 0 ]; then
    echo "Error: Failed to count non-annotated entries."
    exit 1
fi

# Use function to count unique CpG and CRE entries
consensus_gene=$(grep "ENSG" "${annotated_bed_input}" | count_unique_entries /dev/stdin)
consensus_CpG=$(grep "CpG_" "${annotated_bed_input}" | count_unique_entries /dev/stdin)
consensus_CRE=$(grep "EH38E" "${annotated_bed_input}" | count_unique_entries /dev/stdin)
consensus_TFBS=$(grep -v -e "\.\s*\.\s*-1\s*-1\s*\.\s*0" -e "CpG_" -e "EH38E" "${annotated_bed_input}" | count_unique_entries /dev/stdin)

consensus_anno=$((consensus_total - consensus_no_anno))
consensus_anno_non_gene=$((consensus_anno - consensus_gene))
consensus_anno_non_CpG=$((consensus_anno - consensus_CpG))
consensus_anno_non_CRE=$((consensus_anno - consensus_CRE))
consensus_anno_non_TFBS=$((consensus_anno - consensus_TFBS))

# Output the result
{
    echo -e "Total\t$consensus_total"
    echo -e "Annotated\t$consensus_anno"
    echo -e "Non-annotated\t$consensus_no_anno"
    echo -e "Genes\t$consensus_gene"
    echo -e "Non-genes\t$consensus_anno_non_gene"
    echo -e "CpG\t$consensus_CpG"
    echo -e "Non-CpG\t$consensus_anno_non_CpG"
    echo -e "CRE\t$consensus_CRE"
    echo -e "Non-CRE\t$consensus_anno_non_CRE"
    echo -e "TFBS\t$consensus_TFBS"
    echo -e "Non-TFBS\t$consensus_anno_non_TFBS"
} > "${output_file}"

if [ $? -ne 0 ]; then
    echo "Error: Failed to write annotation counts to file."
    exit 1
fi

# Plot the results
echo "Generating plots..."
if ! python "${plot_script_path}" "${output_file}"; then
    echo "Error: Failed to generate plots."
    exit 1
fi

echo "Annotation and plotting completed successfully."