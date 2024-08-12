# About 
---
# short_cfDNA_for_sepsis: Host response Liquid Footprint for improved diagnostics in septic patients

## Introduction
This repository features custom scripts and workflows for the study described in the manuscript: [Host response Liquid Footprint for improved diagnostics in septic patients](https://www.placeholder.org/). It provides code for analyzing short double-stranded cfDNA fragments enriched from human plasma and reproduction of results in the manuscript.

## Abstract
We illustrate the usage of short double-stranded cfDNA, approximately 40 base pairs in length, from cfDNA for high-throughput DNA sequencing. These cfDNA fragments, enriched at regulatory genomic loci and transcription factor binding sites, enable genome-wide DNA footprinting from liquid biopsies. We demonstrate the utility and  its diagnostic potential of short cfDNA for clinical sepsis patients.

## Contents
The repository is organized as follows:
1. **Data Processing:** Scripts for preprocessing, mapping, and postprocessing of short cfDNA sequencing data. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/map_and_process) - "Sequencing data processing"
2. **Peak Calling:** Methods for identifying significant peaks. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/peak_calling) - "Peak calling"
3. **Correlation analysis:** Identification of correlations between clinical metrics and short cfDNA. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/correlation_analysis) - "Correlation analysis"
4. **Differential enrichmen analysis:** Comparison of two experimental groups to determine genome regions with different enrichment of short cfDNA. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/differential_enrichment) - "Differential enrichment analysis"
5. **Annoatation:** Annotation of genomic regions to references of genomic elements. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/annotation) - "Annoatation to genomic functional elements"
6. **TF Motif Enrichment:** Analysis of transcription factor motifs within short cfDNA. [View](https://github.com/janmueller17/short_cfDNA_for_sepsis/tree/main/tf_motif_enrichment) - "Transcription factor motif enrichment analysis"

## Usage
Ensure all dependencies are installed as specified in setup.md (or by installing individual requirements from each script's header).
In the Peak Calling, Correlation Analysis, and TF Motif Enrichment directories, there is a main executable bash script for running these analyses, supported by additional helper scripts located in their respective scripts subfolders.
Reference data are provided where necessary. 
All analyses are designed for execution on Linux systems.

## License
For more information, see the [LICENSE](./LICENSE) file.

## Contact
For further inquiries, please contact us via github or our [website](https://www.igb.fraunhofer.de/en/research/in-vitro-diagnostics.html).

