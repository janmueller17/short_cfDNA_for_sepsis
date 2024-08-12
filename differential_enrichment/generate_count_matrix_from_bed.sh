#!/usr/bin/env bash
# Generate a count matrix from BAM files using multiBamSummary and process it

# .bed file is expected to have region names (not bed3 format).

# Non-standard tools required (not including respective dependencies): 
# multiBamSummary (part of deepTools) - https://deeptools.readthedocs.io/

# Usage:
# ./generate_count_matrix_from_bed.sh <count_matrix_path> <bed_file> <labels> <bam_files> [multiBamSummary_path] [processors]

# Check the number of arguments
if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <count_matrix_path> <bed_file> <labels> <bam_files> [multiBamSummary_path] [processors]"
    echo "Example: $0 /path/to/count_matrix.txt /path/to/reference.bed label1,label2 /path/to/bam1.bam,/path/to/bam2.bam [multiBamSummary_path] [processors]"
    exit 1
fi

# Assign input arguments to variables
count_matrix_path="$1"
bed_file="$2"
IFS=',' read -r -a labels <<< "$3"
IFS=',' read -r -a bam_files <<< "$4"
multiBamSummary_path="${5:-multiBamSummary}"  # Default to "multiBamSummary" if no path is provided
processors="${6:-8}"  # Default to 8 processors if not provided

# Check if multiBamSummary is accessible
if ! command -v "$multiBamSummary_path" &> /dev/null; then
    echo "Error: multiBamSummary is not installed or not found in the provided path."
    exit 1
fi

# Check if BED file exists
if [ ! -f "$bed_file" ]; then
    echo "Error: BED file $bed_file does not exist."
    exit 1
fi

# Check if BAM files exist
for bam_file in "${bam_files[@]}"; do
    if [ ! -f "$bam_file" ]; then
        echo "Error: BAM file $bam_file does not exist."
        exit 1
    fi
done

# Construct labels and bam files strings
labels_str=$(IFS=' '; echo "${labels[*]}")
bam_files_str=$(IFS=' '; echo "${bam_files[*]}")

echo ${labels_str}
echo ${bam_files_str}
echo ${bed_file}
echo "${multiBamSummary_path} BED-file --BED \"$bed_file\" --labels \"$labels_str\" --bamfiles \"$bam_files_str\" --outRawCounts \"$count_matrix_path\" --numberOfProcessors \"$processors\""

# Run multiBamSummary
"${multiBamSummary_path}" BED-file --BED "${bed_file}" --labels ${labels_str} --bamfiles ${bam_files_str} --outRawCounts "${count_matrix_path}" --numberOfProcessors "${processors}"
if [ $? -ne 0 ]; then
    echo "Error: multiBamSummary command failed."
    exit 1
fi

# Process the consensus matrix
awk -i inplace -F '\t' 'NR==1 {print $0} NR>1 {output = $1 ":" $2 "-" $3; for (i=4; i<=NF; i++) output = output "\t" $i; print output}' "${count_matrix_path}"
if [ $? -ne 0 ]; then
    echo "Error: awk command 1 failed."
    exit 1
fi

sed -i "1s/#//; s/'//g" "${count_matrix_path}"
if [ $? -ne 0 ]; then
    echo "Error: sed command failed."
    exit 1
fi

awk -i inplace 'NR==1 {$1="peak"; $2=$3=""; print} NR!=1' "${count_matrix_path}"
if [ $? -ne 0 ]; then
    echo "Error: awk command 2 failed."
    exit 1
fi

echo "Consensus matrix generated and processed successfully: ${count_matrix_path}"
