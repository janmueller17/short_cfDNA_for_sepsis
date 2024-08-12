#!/usr/bin/env bash
# Calculate GC content for BED files using BEDTools and output results

# Non-standard tools required (not including respective dependencies): 
# BEDTools - https://bedtools.readthedocs.io/

# Usage:
# ./calculate_gc_content_from_bed.sh <in_path> <out_path> <reference_fasta> <bed_files> <output_suffix> [bedtools_path]

# Check the number of arguments
if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <in_path> <out_path> <reference_fasta> <bed_files> <output_suffix> [bedtools_path]"
    echo "Example: $0 /path/to/in /path/to/out /path/to/reference.fa 'file1.bed,file2.bed' _GC.bed [bedtools_path]"
    exit 1
fi

# Assign input arguments to variables
in_path="$1"
out_path="$2"
reference_fasta="$3" # Takes a genome .fasta as the reference, e.g. hg19.fa
IFS=',' read -r -a bed_files <<< "$4"
output_suffix="$5"
bedtools_path="${6:-bedtools}"  # Default to "bedtools" if no path is provided

# Check if BEDTools is accessible
if ! command -v "$bedtools_path" &> /dev/null; then
    echo "Error: BEDTools is not installed or not found in the specified path '${bedtools_path}'."
    exit 1
fi

# Check if input and reference files exist
if [ ! -d "$in_path" ]; then
    echo "Error: Input directory ${in_path} does not exist."
    exit 1
fi

if [ ! -d "$out_path" ]; then
    echo "Error: Output directory ${out_path} does not exist."
    exit 1
fi

if [ ! -f "${reference_fasta}" ]; then
    echo "Error: Reference fasta file ${reference_fasta} does not exist."
    exit 1
fi

# Process each BED file
for bed_file in "${bed_files[@]}"; do
    input_file="${in_path}/${bed_file}"
    output_file="${out_path}/${bed_file%.*}${output_suffix}"

    if [ ! -f "${input_file}" ]; then
        echo "Error: Input BED file ${input_file} does not exist."
        continue
    fi

    echo "Processing $input_file ..."
    "${bedtools_path}" nuc -fi "${reference_fasta}" -bed "${input_file}" | awk -v OFS='\t' '{print $1,$2,$3,$5}' > "${output_file}"
    if [ $? -ne 0 ]; then
        echo "Error: BEDTools nuc command failed for ${input_file}."
        continue
    fi

    echo "Output written to ${output_file}"
done

echo "GC content calculation completed for all valid BED files."
