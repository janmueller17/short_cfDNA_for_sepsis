#!/usr/bin/env bash
# Merge two consensus peak sets (or any two BED files) using BEDTools

# Non-standard tools required (not including respective dependencies): 
# BEDTools - https://bedtools.readthedocs.io/

# Usage:
# ./combine_consensus_peaksets.sh <input_file1> <input_file2> <output_file> [bedtools_path]

# Check the number of arguments
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <input_file1> <input_file2> <output_file> [bedtools_path]"
    exit 1
fi

# Assign input arguments to variables
input_file1="$1"
input_file2="$2"
output_file="$3"
bedtools_path="${4:-bedtools}"  # Default to "bedtools" if no path is provided

# Check if BEDTools is accessible
if ! command -v "$bedtools_path" &> /dev/null; then
    echo "Error: BEDTools is not installed or not found in the provided path."
    exit 1
fi

# Merge the BED files
cat "$input_file1" "$input_file2" | sort -k1,1 -k2,2n | "$bedtools_path" merge -i stdin -d 31 > "$output_file"

# Check if the output file was created successfully
if [ $? -eq 0 ]; then
    echo "Merged BED file created successfully: $output_file"
else
    echo "Error: Failed to create merged BED file."
    exit 1
fi