#!/usr/bin/env bash
# Call peaks for short cfDNA Seq data with MACS2.

# Non-standard tools required (not including respective dependencies): 
# MACS2 - https://pypi.org/project/MACS2/

# Usage:
# ./call_peaks.sh <sample_path> <out_path_base> <sample_names> [macs2_path]

# Check the number of arguments
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <sample_path> <out_path_base> <sample_names> [macs2_path]"
    echo "Example: $0 /path/to/in /path/to/out sample_x,sample_y,sample_z [macs2_path]"
    exit 1
fi

# Assign input arguments to variables
sample_path="$1"
out_path_base="$2"
IFS=',' read -r -a sample_names <<< "$3"
macs2_path="${4:-macs2}"  # Default to "macs2" if no path is provided

# Check if MACS2 is accessible
if ! command -v "$macs2_path" &> /dev/null; then
    echo "Error: MACS2 is not installed or not found in the provided path."
    exit 1
fi

# Loop through each sample name and process
for sample_name in "${sample_names[@]}"; do
    out_path="${out_path_base}/${sample_name}"
    
    # Check if the output directory already exists
    if [ -d "${out_path}" ]; then
        echo "Output directory '${out_path}' exists. Rename folder to prevent overwriting."
        continue
    fi

    echo "Processing ${sample_name} ..."
    
    # Create the output directory
    mkdir -p "${out_path}" || { echo "Failed to create output directory '${out_path}'"; exit 1; }

    # Call peaks using MACS2
    "$macs2_path" callpeak -t "${sample_path}/${sample_name}.bam" \
        --outdir "${out_path}" \
        2> "${out_path}/${sample_name}_narrow_peakcall.log" \
        -f BAM -g 2864785220 -n "${sample_name}_macs2_narrow" \
        --nomodel --extsize 32 --call-summits --min-length 30 -q 0.05 || { echo "MACS2 narrow peak calling failed for ${sample_name}"; exit 1; }
done

echo "Peak calling completed for all samples."
