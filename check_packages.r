# check_packages.R

required_packages <- c(
    "ggplot2", "ggpubr", "dplyr", "tidyverse", "scales", "optparse",
    "EDASeq", "data.table", "readr", "gplots", "corrplot", "Rfast",
    "pals", "MASS", "edgeR", "RColorBrewer", "rtracklayer",
    "GenomicRanges", "BRGenomics", "pROC"
)

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste("Package not found:", pkg))
    } else {
        cat(paste("Package loaded successfully:", pkg, "\n"))
    }
}