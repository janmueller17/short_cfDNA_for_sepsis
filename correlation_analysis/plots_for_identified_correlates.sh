#!/usr/bin/env bash

# This script executes the plotting for identified correlates between clinical metrics and loci. 

# Non-standard tools required (not including respective dependencies):
# Python with pandas, numpy, scipy, matplotlib, argparse, os

# Usage:
# ./plots_for_identified_correlates.sh \
#     --cfoi <clinical metric> \
#     --locus_interest <selected locus> \
#     --count_data_csv /path/to/readcount_data.csv \
#     --samples_metadata ./references/ident_patients_IDs.csv \
#     --corr_val_df /path/to/combined_corr_val_df.csv \
#     --corr_p_adj_df /path/to/combined_corr_p_adj_df.csv \
#     --output_base_path /path/to/output/ \
#     --clinical_metrics ./references/patient_metadata.csv \
#     --patients S10 S33 S39 S42 S43 S49 \
#     --python_path /path/to/python

# Basic usage statement
usage() {
    echo "Usage: $0 [options] --patients patient1 patient2 ..."
    echo "Options:"
    echo "  --cfoi VALUE"
    echo "  --locus_interest VALUE"
    echo "  --count_data_csv PATH"
    echo "  --samples_metadata PATH"
    echo "  --corr_val_df PATH"
    echo "  --corr_p_adj_df PATH"
    echo "  --output_base_path PATH"
    echo "  --clinical_metrics PATH"
    echo "  --python_path PATH (optional, default: /fungen/funhome/Software/miniconda3/envs/r_env_jam/bin/python)"
    exit 1
}

# Check if at least one argument is passed
if [ $# -lt 1 ]; then
    usage
fi

# Set defaults 
PYTHON_PATH="python"
script_dir="$(dirname "$0")"

# Parse command line arguments
while [[ "$1" != "" ]]; do
    case $1 in
        --cfoi ) shift; CFOI=$1 ;;
        --locus_interest ) shift; LOCUS_INTEREST=$1 ;;
        --count_data_csv ) shift; COUNT_DATA_CSV=$1 ;;
        --samples_metadata ) shift; SAMPLES_METADATA=$1 ;;
        --corr_val_df ) shift; CORR_VAL_DF=$1 ;;
        --corr_p_adj_df ) shift; CORR_P_ADJ_DF=$1 ;;
        --output_base_path ) shift; OUTPUT_BASE_PATH=$1 ;;
        --clinical_metrics ) shift; CLINICAL_METRICS=$1 ;;
        --patients ) shift; PATIENTS=(); while [[ "$1" != "" && "$1" != --* ]]; do PATIENTS+=("$1"); shift; done ;;
        --python_path ) shift; PYTHON_PATH=$1 ;;
        * ) usage ;;
    esac
    if [[ "$1" != --* ]]; then
        shift
    fi
done

# Check if all required arguments are set
if [ -z "$CFOI" ] || [ -z "$LOCUS_INTEREST" ] || [ -z "$COUNT_DATA_CSV" ] || [ -z "$SAMPLES_METADATA" ] || [ -z "$CORR_VAL_DF" ] || [ -z "$CORR_P_ADJ_DF" ] || [ -z "$OUTPUT_BASE_PATH" ] || [ -z "$CLINICAL_METRICS" ]; then
    usage
fi

# Define common arguments as an array
common_args=(--cfoi "$CFOI" --locus_interest "$LOCUS_INTEREST" --count_data_csv "$COUNT_DATA_CSV" --samples_metadata "$SAMPLES_METADATA" --corr_val_df "$CORR_VAL_DF" --corr_p_adj_df "$CORR_P_ADJ_DF" --output_base_path "$OUTPUT_BASE_PATH" --clinical_metrics "$CLINICAL_METRICS")

# Run script without patient argument
"${PYTHON_PATH}" "${script_dir}/scripts/plot_clinical_feature_vs_readcount.py" "${common_args[@]}"

# Run script for each patient
for patient in "${PATIENTS[@]}"; do
    echo "Running for patient: ${patient}"
    "${PYTHON_PATH}" "${script_dir}/scripts/plot_clinical_feature_vs_readcount.py" "${common_args[@]}" --patient "${patient}"
done