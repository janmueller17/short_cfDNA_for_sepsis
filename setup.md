# Setup Instructions

## Conda Environment

### Create the Conda Environment

To create the Conda environment with all the required dependencies, run the following command:

```bash
conda env create -f conda_env.yml
```

### Activate the Environment

After creating the environment, activate it using:

```bash
conda activate short_cfDNA_for_sepsis
```

### Deactivating the Environment

When you are done working in the environment, you can deactivate it with:

```bash
conda deactivate
```

### Updating the Environment

If you need to update the environment with new dependencies, you can modify the `conda_env.yml` file and then run:

```bash
conda env update -f conda_env.yml --prune
```

This will update the existing environment with any new dependencies while removing any dependencies that are no longer required.

### Removing the Environment

If you need to remove the environment, you can do so with:

```bash
conda remove --name short_cfDNA_for_sepsis --all
```

This will delete the environment and all its packages.

## Verifying Installation

### Python Packages

To check if all Python packages are installed correctly, run the following command:

```bash
python -c "import pandas, numpy, scipy, matplotlib, argparse, multiprocessing, os, sys, sklearn, statsmodels"
```

### R Packages

To check if all R packages are installed correctly, run the following command:
```bash
Rscript check_packages.R
```

### Non-standard Tools

To check if non-standard tools are installed correctly, use their respective commands:

```bash
bedtools --version
deeptools --version
fastqc --version
bbmap.sh --version
prinseq-lite.pl --version
ngm --version
samtools --version
macs2 --version
ame --version
```

This should ensure that your environment is properly set up and ready for use.