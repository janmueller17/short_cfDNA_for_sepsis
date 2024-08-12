#!/usr/bin/env bash
# Process and map short cfDNA seq data
# This version keeps all temp files (greppable with '_temp'), these can be removed if not needed.

# Non-standard tools required (not including respective dependencies): 
# FastQC - https://github.com/s-andrews/FastQC
# bbduk - https://github.com/BioInfoTools/BBMap
# prinseq-lite - https://github.com/uwb-linux/prinseq
# NextGenMap (ngm) - https://github.com/Cibiv/NextGenMap
# samtools - https://github.com/samtools/samtools
# deepTools - https://github.com/deeptools/deepTools

# Default paths and parameters
default_bbduk_path="bbduk.sh"
default_fastqc_path="fastqc"
default_prinseq_path="prinseq-lite.pl"
default_ngm_path="ngm"
default_samtools_path="samtools"
default_deeptools_path="bamCoverage"
default_genome_reference="/path/to/genome_reference/genome.fasta"
default_genome_reference_eff_size=2864785220
default_seq_adapter_reference="/path/to/adapter_reference/adapter.fasta"
default_blacklist_bed="/path/to/ENCODE_consensusBlacklist.bed"
default_threads=4

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --sample-path) sample_path="$2"; shift ;;
    --out-path-base) out_path_base="$2"; shift ;;
    --threads) threads="$2"; shift ;;
    --bbduk-path) bbduk_path="$2"; shift ;;
    --fastqc-path) fastqc_path="$2"; shift ;;
    --prinseq-path) prinseq_path="$2"; shift ;;
    --ngm-path) ngm_path="$2"; shift ;;
    --samtools-path) samtools_path="$2"; shift ;;
    --deeptools-path) deeptools_path="$2"; shift ;;
    --genome-reference) genome_reference="$2"; shift ;;
    --genome-reference-eff-size) genome_reference_eff_size="$2"; shift ;;
    --seq-adapter-reference) seq_adapter_reference="$2"; shift ;;
    --blacklist-bed) blacklist_bed="$2"; shift ;;
    --sample-names) sample_names="$2"; shift ;;
    --sample-name) sample_name="$2"; single_sample=true; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

# Set defaults if not provided
bbduk_path="${bbduk_path:-$default_bbduk_path}"
fastqc_path="${fastqc_path:-$default_fastqc_path}"
prinseq_path="${prinseq_path:-$default_prinseq_path}"
ngm_path="${ngm_path:-$default_ngm_path}"
samtools_path="${samtools_path:-$default_samtools_path}"
deeptools_path="${deeptools_path:-$default_deeptools_path}"
genome_reference="${genome_reference:-$default_genome_reference}"
genome_reference_eff_size="${genome_reference_eff_size:-$default_genome_reference_eff_size}"
seq_adapter_reference="${seq_adapter_reference:-$default_seq_adapter_reference}"
blacklist_bed="${blacklist_bed:-$default_blacklist_bed}"
threads="${threads:-$default_threads}"

# Check required arguments
if [[ -z "$sample_path" || -z "$out_path_base" ]]; then
  echo "Missing required arguments: --sample-path and --out-path-base are required."
  exit 1
fi

if [[ -n "$single_sample" && -z "$sample_name" ]]; then
  echo "Missing required argument: --sample-name is required when --single-sample is used."
  exit 1
fi

if [[ -z "$single_sample" && -z "$sample_names" ]]; then
  echo "Missing required argument: --sample-names is required when not using --single-sample."
  exit 1
fi

# Convert sample_names string to array if provided
if [[ -z "$single_sample" ]]; then
  IFS=',' read -r -a sample_names <<< "$sample_names"
else
  sample_names=("$sample_name")
fi

# Function to check if a file exists
check_file() {
  if [[ ! -f "$1" ]]; then
    echo "File not found: $1"
    exit 1
  fi
}

# Check if essential files exist
check_file "$genome_reference"
check_file "$seq_adapter_reference"
check_file "$blacklist_bed"

# Process each sample
for sample_name in "${sample_names[@]}"; do

  out_path="${out_path_base}/${sample_name}"

  if [ -d "${out_path}" ]; then
    echo " "
    echo "Output directory '${out_path}' exists. Rename folder to prevent overwriting."
    continue

  else
    echo "Processing ${sample_name} ..."
    mkdir -p "${out_path}" || { echo "Failed to create output directory '${out_path}'"; exit 1; }

    # Initial QC of raw data
    "${fastqc_path}" -o "${out_path}" --noextract -t "${threads}" "${sample_path}/${sample_name}_R1.fastq.gz" || { echo "Initial QC failed"; exit 1; }
    echo "Initial QC finished"

    # First read processing: artifact filtering (terminal A, polyG, etc.)
    "${bbduk_path}" trimpolyg=10 hdist=1 mink=12 threads=${threads} maxns=10 minlen=20 maxlength=60 trimq=20 qtrim=t ktrim=r k=28 ref="${seq_adapter_reference}" in="${sample_path}/${sample_name}_R1.fastq.gz" out="${out_path}/${sample_name}_R1_processed_temp1.fastq.gz" || { echo "BBDuk step 1 failed"; exit 1; } # polyG trimming, size filter, adapter trimming
    "${bbduk_path}" forcetrimright2=1 threads=${threads} ref="${seq_adapter_reference}" in="${out_path}/${sample_name}_R1_processed_temp1.fastq.gz" out="${out_path}/${sample_name}_R1_processed_temp2.fastq.gz" || { echo "BBDuk step 2 failed"; exit 1; } # 1 base right side trim to remove single A from library prep
    
    gunzip "${out_path}/${sample_name}_R1_processed_temp2.fastq.gz" || { echo "gunzip failed"; exit 1; }
    "${prinseq_path}" -lc_method dust -lc_threshold 7 -fastq "${out_path}/${sample_name}_R1_processed_temp2.fastq" -out_good "${out_path}/${sample_name}_R1_processed" -out_bad "${out_path}/${sample_name}_R1_processed_bad" || { echo "PRINSEQ-lite processing failed"; exit 1; } # low complexity filtering
    gzip "${out_path}/${sample_name}_R1_processed.fastq" || { echo "gzip failed"; exit 1; }
    
    "${fastqc_path}" -o "${out_path}" --noextract -t "${threads}" "${out_path}/${sample_name}_R1_processed.fastq.gz" || { echo "QC after processing failed"; exit 1; } # get QC stats before processing
    echo "First read processing finished"

    # Read mapping to genome and basic processing steps (remove unmapped, sort, and index)
    "${ngm_path}" -t "${threads}" --no-progress -b -o "${out_path}/${sample_name}.bam" -r "${genome_reference}" -q "${out_path}/${sample_name}_R1_processed.fastq.gz" 2> "${out_path}/${sample_name}_ngm_mapping_report.txt" || { echo "NGM mapping failed"; exit 1; }
    "${samtools_path}" view -@ "${threads}" -b -f 4 "${out_path}/${sample_name}.bam" > "${out_path}/${sample_name}_unmapped.bam" || { echo "Samtools view failed"; exit 1; }
    "${samtools_path}" sort -@ "${threads}" -o "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}.bam" || { echo "Samtools sort failed"; exit 1; }
    "${samtools_path}" index -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" || { echo "Samtools index failed"; exit 1; }
    "${samtools_path}" flagstat -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" > "${out_path}/${sample_name}_stats_before_processing.txt" || { echo "Samtools flagstat before processing failed"; exit 1; }
    echo "Mapping finished"

    # Further processing
    "${samtools_path}" markdup -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}_dup_flagged.bam" || { echo "Samtools markdup flag failed"; exit 1; } # flag duplicates 
    "${samtools_path}" markdup -r -s -@ "${threads}" "${out_path}/${sample_name}_sorted.bam" "${out_path}/${sample_name}_dup_rem_temp.bam" || { echo "Samtools markdup remove failed"; exit 1; } # remove duplicates 
    "${samtools_path}" view -H "${out_path}/${sample_name}_dup_rem_temp.bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | "${samtools_path}" reheader - "${out_path}/${sample_name}_dup_rem_temp.bam" > "${out_path}/${sample_name}_chr_fix_temp.bam" || { echo "Reheadering failed"; exit 1; } # filter for autosomal and sex chromosomes only (chr 1 to 22 X and Y) and change chromosome names to UCSC style (1 -> chr1)
    "${samtools_path}" view "${out_path}/${sample_name}_chr_fix_temp.bam" -@ "${threads}" -o "${out_path}/${sample_name}_reads_in_bl_temp.bam" -U "${out_path}/${sample_name}_bl_rem_temp.bam" -L "${blacklist_bed}" || { echo "Blacklist removal failed"; exit 1; } # remove reads mapped to the ENCODE blacklist
    "${samtools_path}" view "${out_path}/${sample_name}_bl_rem_temp.bam" -@ "${threads}" -q 1 -o "${out_path}/${sample_name}_processed.bam" || { echo "Samtools view with mapQ filter failed"; exit 1; } # keep only reads with mapQ > 1
    "${samtools_path}" index -@ "${threads}" "${out_path}/${sample_name}_processed.bam" || { echo "Samtools index after processing failed"; exit 1; }
    "${samtools_path}" flagstat -@ "${threads}" "${out_path}/${sample_name}_processed.bam" > "${out_path}/${sample_name}_stats_after_processing.txt" || { echo "Samtools flagstat after processing failed"; exit 1; } # get QC stats after processing
    "${deeptools_path}" --bam "${out_path}/${sample_name}_processed.bam" -o "${out_path}/${sample_name}_processed.bw" --outFileFormat bigwig --normalizeUsing CPM --binSize 10 --numberOfProcessors "${threads}" --effectiveGenomeSize "${genome_reference_eff_size}" || { echo "bamCoverage failed"; exit 1; } # generate bigwig from processed reads
    echo "Processing finished"
  fi

done