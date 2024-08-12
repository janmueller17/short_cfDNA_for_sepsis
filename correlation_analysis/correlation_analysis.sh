#!/usr/bin/env bash

# This script performs a comprehensive correlation analysis to find genomic locations whose sequencing read counts correlate to clinical metrics. 

# Non-standard tools required (not including respective dependencies):
# BEDTools - https://bedtools.readthedocs.io/
# deepTools - https://github.com/deeptools/deepTools
# Python with pandas, numpy, scipy, statsmodels, matplotlib, argparse, multiprocessing, os, sklearn
# R with EDASeq, data.table, tools, readr

# Usage: ./correlation_analysis.sh -i ./references/ident_patients_bam_paths.txt -t ./references/test_patients_bam_paths.txt -m ./references/patient_metadata.csv -I ./references/ident_patients_IDs.csv -T ./references/test_patients_IDs.csv -r </path/to/genome.fasta> -o </path/to/output_dir> -n <num_threads> -b <bin_size> -c "Leukocytes Erythrocytes Creatinine Bilirubin Albumine Urea ALT AST ALP Hemoglobine Thrombocytes Quick aPTT CRP PCT" [-B <bedtools_path>] [-M <multiBamSummary_path>] [-P <python_path>] [-R <Rscript_path>]

# Detailed usage statement
usage() {
  echo "Usage: $0 -i <ident_samples_file> -t <test_samples_file> -o <output_dir> -r <reference_fa> -m <metadata_csv> -I <ident_ids_csv> -T <test_ids_csv> -n <num_threads> -b <bin_size> -c <clinical_metrics> [-B <bedtools_path>] [-M <multiBamSummary_path>] [-P <python_path>] [-R <Rscript_path>]"
  echo
  echo "Arguments:"
  echo "  -i|--ident-samples          File containing labels and BAM file paths for ident samples"
  echo "  -t|--test-samples           File containing labels and BAM file paths for test samples"
  echo "  -o|--output-dir             Directory to store output files"
  echo "  -r|--reference-fa           Reference fasta file"
  echo "  -m|--clinical-metrics-csv   Metadata CSV file"
  echo "  -I|--ident-ids-csv          Ident IDs CSV file"
  echo "  -T|--test-ids-csv           Test IDs CSV file"
  echo "  -n|--num-threads            Number of threads to use"
  echo "  -b|--bin-size               Bin size for multiBamSummary"
  echo "  -c|--clinical-metrics       Quoted whitespace-separated list of clinical metrics to analyze"
  echo "  -B|--bedtools-path          (optional) Path to bedtools executable (default: bedtools)"
  echo "  -M|--multiBamSummary-path   (optional) Path to multiBamSummary executable (default: multiBamSummary)"
  echo "  -P|--python-path            (optional) Path to python executable (default: python)"
  echo "  -R|--Rscript-path           (optional) Path to Rscript executable (default: Rscript)"
  echo "  -h|--help                   Show this help message and exit"
  echo
  exit 1
}

# Parse CLI arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--ident-samples) ident_samples_file="$2"; shift ;;
        -t|--test-samples)  test_samples_file="$2"; shift ;;
        -o|--output-dir) output_dir="$2"; shift ;;
        -r|--reference-fa) reference_fa="$2"; shift ;;
        -m|--clinical-metrics-csv) clinical_metrics_csv="$2"; shift ;;
        -I|--ident-ids-csv) ident_ids_csv="$2"; shift ;;
        -T|--test-ids-csv) test_ids_csv="$2"; shift ;;
        -n|--num-threads) num_threads="$2"; shift ;;
        -b|--bin-size) bin_size="$2"; shift ;;
        -c|--clinical-metrics) clinical_metrics=("$2"); shift ;;
        -B|--bedtools-path) bedtools_path="$2"; shift ;;
        -M|--multiBamSummary-path) multiBamSummary_path="$2"; shift ;;
        -P|--python-path) python_path="$2"; shift ;;
        -R|--Rscript-path) Rscript_path="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "${ident_samples_file}" ] || [ -z "${test_samples_file}" ] || [ -z "${output_dir}" ] || [ -z "${reference_fa}" ] || [ -z "${clinical_metrics_csv}" ] || [ -z "${ident_ids_csv}" ] || [ -z "${test_ids_csv}" ] || [ -z "${num_threads}" ] || [ -z "${bin_size}" ] || [ -z "${clinical_metrics[*]}" ]; then
  usage
fi

# Set default paths if not provided
bedtools_path=${bedtools_path:-bedtools}
multiBamSummary_path=${multiBamSummary_path:-multiBamSummary}
python_path=${python_path:-python}
Rscript_path=${Rscript_path:-Rscript}

script_dir="$(dirname "$0")"


# Check if input files exist
if [ ! -f "${ident_samples_file}" ]; then
  echo "Error: ident_samples_file '${ident_samples_file}' not found."
  exit 1
fi

if [ ! -f "${test_samples_file}" ]; then
  echo "Error: test_samples_file '${test_samples_file}' not found."
  exit 1
fi

if [ ! -f "${reference_fa}" ]; then
  echo "Error: reference_fa '${reference_fa}' not found."
  exit 1
fi

if [ ! -f "${clinical_metrics_csv}" ]; then
  echo "Error: clinical_metrics_csv '${clinical_metrics_csv}' not found."
  exit 1
fi

if [ ! -f "${ident_ids_csv}" ]; then
  echo "Error: ident_ids_csv '${ident_ids_csv}' not found."
  exit 1
fi

if [ ! -f "${test_ids_csv}" ]; then
  echo "Error: test_ids_csv '${test_ids_csv}' not found."
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Function to read labels and BAM files
read_labels_and_bamfiles() {
  local file=$1
  local -n labels_ref=$2
  local -n bamfiles_ref=$3
  while IFS= read -r line; do
    labels_ref+=("$(echo "${line}" | awk '{print $1}')")
    bamfiles_ref+=("$(echo "${line}" | awk '{print $2}')")
  done < "${file}"
}

# Arrays to hold labels and BAM files
labels_ident=()
bamfiles_ident=()

# Read labels and BAM files for ident samples
read_labels_and_bamfiles "${ident_samples_file}" labels_ident bamfiles_ident

# Define output CSV file paths
ident_counts_tsv="${output_dir}/readcounts_ident_samples_${bin_size}_bin.tsv"
ident_gc_bed="${output_dir}/${bin_size}_bin_GC.bed"
ident_counts_adj_csv="${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts.csv"

test_counts_tsv="${output_dir}/readcounts_test_samples_${bin_size}_bin.tsv"
test_counts_adj_csv="${output_dir}/readcounts_test_samples_${bin_size}_bin_GCadj_DEPTHadj_counts.csv"

# Verbose output
echo "Running multiBamSummary for ident samples..."
echo "BAM files:" 
echo "${bamfiles_ident[@]}"
echo "Labels: "
echo "${labels_ident[@]}"

# Run multiBamSummary for ident samples
"${multiBamSummary_path}" bins \
  --numberOfProcessors "${num_threads}" \
  --binSize "${bin_size}" \
  --labels "${labels_ident[@]}" \
  --bamfiles "${bamfiles_ident[@]}" \
  --outRawCounts "$ident_counts_tsv"

# Check if multiBamSummary for ident samples completed successfully
if [ $? -ne 0 ] || [ ! -f "$ident_counts_tsv" ]; then
  echo "Error: multiBamSummary for ident samples failed."
  exit 1
fi

# Extract bins
echo "Extracting bins..."
awk 'NR>1 {print $1"\t"$2"\t"$3}' "${ident_counts_tsv}" > "${output_dir}/${bin_size}_bin.bed"

# Run bedtools nuc
echo "Running bedtools nuc..."
"${bedtools_path}" nuc \
  -fi "${reference_fa}" -bed "${output_dir}/${bin_size}_bin.bed" | \
  awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' > "${ident_gc_bed}"

# Check if bedtools nuc completed successfully
if [ $? -ne 0 ] || [ ! -f "$ident_gc_bed" ]; then
  echo "Error: bedtools nuc failed."
  exit 1
fi

# Adjust count matrix with Rscript
echo "Adjusting count matrix with Rscript..."
"${Rscript_path}" "${script_dir}/scripts/adjust_count_mat.R" \
  "${ident_gc_bed}" "${ident_counts_tsv}" "${output_dir}"

# Check if Rscript for adjusting count matrix completed successfully
if [ $? -ne 0 ] || [ ! -f "$ident_counts_adj_csv" ]; then
  echo "Error: Rscript for adjusting count matrix failed."
  exit 1
fi

# Arrays to hold labels and BAM files
labels_test=()
bamfiles_test=()

# Read labels and BAM files for test samples
read_labels_and_bamfiles "${test_samples_file}" labels_test bamfiles_test

# Run multiBamSummary for test samples
echo "Running multiBamSummary for test samples..."
echo "BAM files:" 
echo "${bamfiles_test[@]}"
echo "Labels: "
echo "${labels_test[@]}"

"${multiBamSummary_path}" BED-file \
  --BED "${ident_gc_bed}" --numberOfProcessors "${num_threads}" --labels "${labels_test[@]}" \
  --bamfiles "${bamfiles_test[@]}" \
  --outRawCounts "${test_counts_tsv}"

# Check if multiBamSummary for test samples completed successfully
if [ $? -ne 0 ] || [ ! -f "$test_counts_tsv" ]; then
  echo "Error: multiBamSummary for test samples failed."
  exit 1
fi

# Adjust test count matrix with Rscript
echo "Adjusting test count matrix with Rscript..."
"${Rscript_path}" "${script_dir}/scripts/adjust_count_mat.R" \
  "${ident_gc_bed}" "${test_counts_tsv}" "${output_dir}"

# Check if Rscript for adjusting test count matrix completed successfully
if [ $? -ne 0 ] || [ ! -f "$test_counts_adj_csv" ]; then
  echo "Error: Rscript for adjusting test count matrix failed."
  exit 1
fi

# Run correlation analysis
echo "Running correlation analysis..."

"${python_path}" "${script_dir}/scripts/correlation_readcounts_clinical.py" \
  -t "${num_threads}" -m "${clinical_metrics_csv}" -p "${ident_ids_csv}" -c pearson \
  -i "${ident_counts_adj_csv}" \
  -v "${clinical_metrics}" \
  -o "${output_dir}"

# Check if correlation analysis completed successfully
if [ $? -ne 0 ]; then
  echo "Error: correlation analysis failed."
  exit 1
fi

# Process correlation results
echo "Processing correlation results..."
"${python_path}" "${script_dir}/scripts/process_correlation_results.py" \
  -a 0.05 \
  -v "${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts_combined_pearson_corr_vals.csv" \
  -p "${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts_combined_pearson_corr_p.csv"

# Check if correlation results processing completed successfully
if [ $? -ne 0 ]; then
  echo "Error: processing correlation results failed."
  exit 1
fi

# Run linear model analysis
echo "Running linear model analysis..."
"${python_path}" "${script_dir}/scripts/linear_model_clinical_predictor.py" \
  --clinical_metrics "${clinical_metrics}" \
  --base_path "${output_dir}" \
  --test_samples_csv "${test_ids_csv}" \
  --ident_samples_csv "${ident_ids_csv}" \
  --ident_avg_correlation_csv "${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts_pearson_individuals_avg_corr_vals.csv" \
  --ident_adj_p_csv "${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts_combined_pearson_corr_p_adj.csv" \
  --ident_comb_corr_csv "${output_dir}/readcounts_ident_samples_${bin_size}_bin_GCadj_DEPTHadj_counts_combined_pearson_corr_vals.csv" \
  --clinical_metrics_csv "${clinical_metrics_csv}" \
  --alpha 0.05 \
  --test_samples_count_mat_csv "${test_counts_adj_csv}" \
  --ident_samples_count_mat_csv "${ident_counts_adj_csv}"

# Check if linear model analysis completed successfully
if [ $? -ne 0 ]; then
  echo "Error: linear model analysis failed."
  exit 1
fi

echo "Analysis complete. Output files are in ${output_dir}"