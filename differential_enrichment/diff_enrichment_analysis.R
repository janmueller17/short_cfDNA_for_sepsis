#!/usr/bin/env Rscript

# This script performs differential enrichment analysis using EDASeq on a short cfDNA seq data count matrix. It reads the data, preprocesses it, normalizes it, performs differential analysis using EDASeq (and edgeR), and generates plots for visualization.

# Example usage:
# Rscript ./diff_enrichment_analysis.R -b /path/to/input -o /path/to/output -n infection_type -c "Viral,Viral,Viral,Viral,Viral,Viral,Viral,Viral,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial" -z 8 -t 10 -m 0.1 -s 0.05 -f 1.4 -x "#0D9258,#000000"


# Load necessary libraries
library(EDASeq)
library(readr)
library(edgeR)
library(ggplot2)
library(gplots)
library(data.table)
library(RColorBrewer)
library(Rfast)
library(pROC)
library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-b", "--base_path"), type = "character", default = "path/to/project/",
              help = "Base path [default = %default]", metavar = "character"),
  make_option(c("-o", "--out_path"), type = "character", default = NULL,
              help = "Output path [default = <base_path>/out/]", metavar = "character"),
  make_option(c("-n", "--name"), type = "character", default = "infection_type",
              help = "Name of the analysis [default = %default]", metavar = "character"),
  make_option(c("-c", "--conditions"), type = "character", default = "Viral,Viral,Viral,Viral,Viral,Viral,Viral,Viral,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial,Bacterial",
              help = "Conditions formatted as a comma-separated list [default = %default]", metavar = "character"),
  make_option(c("-z", "--zeros_per_cond"), type = "integer", default = 8,
              help = "Zeros per condition threshold [default = %default]", metavar = "integer"),
  make_option(c("-t", "--top_to_select"), type = "integer", default = 10,
              help = "Top loci to select [default = %default]", metavar = "integer"),
  make_option(c("-m", "--mean_reads_per_region_thresh"), type = "numeric", default = 0.1,
              help = "Mean reads per region threshold [default = %default]", metavar = "numeric"),
  make_option(c("-s", "--sig_thresh"), type = "numeric", default = 0.05,
              help = "Significance threshold [default = %default]", metavar = "numeric"),
  make_option(c("-f", "--fc_thresh"), type = "numeric", default = 1.4,
              help = "Fold change threshold [default = %default]", metavar = "numeric"),
  make_option(c("-x", "--colors"), type = "character", default = "#CC79A7,#E69F00",
              help = "Colors for conditions separated by comma [default = %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set paths and parameters
base_path <- opt$base_path
out_path <- ifelse(is.null(opt$out_path), paste0(base_path, "Differential_analysis/"), opt$out_path)
if (!file.exists(out_path)) dir.create(out_path, recursive = TRUE)

mean_reads_per_region_thresh <- opt$mean_reads_per_region_thresh
sig_thresh <- opt$sig_thresh
fc_thresh <- opt$fc_thresh
top_to_select <- opt$top_to_select
name <- opt$name

# Parse conditions
conditions_in <- factor(strsplit(opt$conditions, ",")[[1]])

zeros_per_cond <- opt$zeros_per_cond

# Set colors for conditions
color_list <- strsplit(opt$colors, ",")[[1]]
if (length(color_list) != length(levels(conditions_in))) {
  stop("Number of colors provided does not match number of conditions.")
}
col_1 <- color_list[1]
col_2 <- color_list[2]

# Function to read and preprocess data
read_and_preprocess_data <- function(base_path, name, conditions_in, mean_reads_per_region_thresh, zeros_per_cond) {
  GC_content_df <- read_delim(paste0(base_path, name, "_GC.bed"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  count_mat_df <- read_table(paste0(base_path, name, ".txt"))

  GC_content_df$peak <- paste0(GC_content_df$`#1_usercol`, ":", GC_content_df$`2_usercol`, "-", GC_content_df$`3_usercol`)
  GC_content_df$length <- abs(GC_content_df$`2_usercol` - GC_content_df$`3_usercol`)

  count_mat_df <- count_mat_df[order(count_mat_df$peak), ]
  GC_content_df <- GC_content_df[order(GC_content_df$peak), ]

  pheno_df <- AnnotatedDataFrame(data.frame(conditions = conditions_in))
  rownames(pheno_df) <- colnames(count_mat_df[, -1])

  count_mat_arr <- as.matrix(count_mat_df[, -1])
  rownames(count_mat_arr) <- count_mat_df$peak

  filter1 <- apply(count_mat_arr, 1, function(x) mean(x) > mean_reads_per_region_thresh)
  count_mat_arr_filt_temp <- count_mat_arr[filter1, ]
  GC_content_df_filt_temp <- GC_content_df[filter1, ]

  design <- model.matrix(~conditions, data = pData(pheno_df))
  filter2 <- (rowSums(count_mat_arr_filt_temp[, design[, 2] == 0] == 0) <= zeros_per_cond) & (rowSums(count_mat_arr_filt_temp[, design[, 2] == 1] == 0) <= zeros_per_cond)
  print(sum(filter2))
  print(length(filter2))
  count_mat_arr_filt <- count_mat_arr_filt_temp[filter2, ]
  GC_content_df_filt <- GC_content_df_filt_temp[filter2, ]

  feature <- AnnotatedDataFrame(data.frame(gc = GC_content_df_filt$`5_pct_gc`, length = GC_content_df_filt$length, row.names = GC_content_df_filt$peak))
  data <- newSeqExpressionSet(counts = count_mat_arr_filt, featureData = feature, phenoData = pheno_df)

  list(data = data, pheno_df = pheno_df)
}

# Function for normalization
normalize_data <- function(data) {
  data_offset <- withinLaneNormalization(data, "gc", which = "full", offset = TRUE)
  data_offset <- betweenLaneNormalization(data_offset, which = "full", offset = TRUE)
  data_offset
}

# Function for differential analysis
differential_analysis <- function(data_offset, conditions_in, sig_thresh, fc_thresh, top_to_select, out_path, name) {
  design <- model.matrix(~conditions, data = pData(data_offset))
  y <- DGEList(counts = counts(data_offset), group = pData(data_offset)$conditions)
  y$offset <- -offst(data_offset)
  y <- estimateDisp(y, design)

  pdf(file = paste0(out_path, "/", name, "_edgeR_BCV.pdf"), width = 5, height = 5)
  plotBCV(y)
  dev.off()

  fit <- glmFit(y, design)
  test <- glmLRT(fit, coef = 2)
  res_eR <- topTags(test, n = dim(y)[1])

  log2fc_eR <- res_eR$table$logFC
  pval_eR <- res_eR$table$FDR

  res_eR_FDR <- res_eR[res_eR$table$FDR <= sig_thresh, ]
  res_eR_FDR_FC <- res_eR_FDR[which(res_eR_FDR$table$logFC > log2(fc_thresh) | res_eR_FDR$table$logFC < -log2(fc_thresh)), ]
  diff_loci_id_list <- row.names(res_eR_FDR_FC$table)
  diff_loci_adj_counts <- test$fitted.values[diff_loci_id_list, ]

  cond1_cv <- rowcvs(x = test$fitted.values[, design[, 2] == 0])
  cond2_cv <- rowcvs(x = test$fitted.values[, design[, 2] == 1])
  var_cond_df <- data.frame(cbind(row.names(test$fitted.value), cond1_cv, cond2_cv))
  colnames(var_cond_df)[1] <- "rn"
  var_cond_df$cond1_cv <- as.numeric(var_cond_df$cond1_cv)
  var_cond_df$cond2_cv <- as.numeric(var_cond_df$cond2_cv)

  export_df <- data.frame(merge(x = as.data.table(test$fitted.values, keep.rownames = TRUE), y = as.data.table(res_eR$table, keep.rownames = TRUE), sort = FALSE, by = "rn"))
  export_df <- data.frame(merge(x = export_df, y = var_cond_df, sort = FALSE, by = "rn"))

  write.csv(x = export_df, file = paste0(out_path, name, "_edgeR_norm_reads_test_results.csv"), quote = FALSE, row.names = FALSE)


  res_eR_FCpos <- res_eR[which(res_eR$table$logFC > 0), ]
  res_eR_FCneg <- res_eR[which(res_eR$table$logFC < 0), ]
  res_eR_FCpos <- res_eR_FCpos[order(res_eR_FCpos$table$PValue), ]
  res_eR_FCneg <- res_eR_FCneg[order(res_eR_FCneg$table$PValue), ]

  res_eR_FCpos_top <- head(res_eR_FCpos, top_to_select %/% 2)
  res_eR_FCneg_top <- head(res_eR_FCneg, top_to_select %/% 2)

  top_loci_id_list <- c(row.names(res_eR_FCneg_top$table), row.names(res_eR_FCpos_top$table))
  top_loci_export_df <- export_df[export_df$rn %in% top_loci_id_list, ]
  
  write.csv(x = top_loci_export_df, file = paste0(out_path, name, "_edgeR_norm_reads_test_results_top.csv"), quote = FALSE, row.names = FALSE)


  list(log2fc_eR = log2fc_eR, pval_eR = pval_eR, diff_loci_adj_counts = diff_loci_adj_counts, res_eR = res_eR, design = design, test = test)
}

# Function to generate plots
generate_plots <- function(log2fc_eR, pval_eR, diff_loci_adj_counts, conditions_in, col_1, col_2, sig_thresh, fc_thresh, design, out_path, name, top_to_select, res_eR, test) {
  colors <- ifelse(log2fc_eR > log2(fc_thresh) & pval_eR <= sig_thresh, col_1, ifelse(log2fc_eR < -log2(fc_thresh) & pval_eR <= sig_thresh, col_2, "#979797"))

  pdf(file = paste0(out_path, name, "_edgeR_volcano.pdf"), width = 4, height = 5)
  plot(log2fc_eR, -log10(pval_eR), pch = ifelse(pval_eR <= sig_thresh, 16, 1), cex = 0.5,
       col = colors, xlab = "log2 Fold Change", ylab = "-log10(Adjusted p-value)", main = name)
  dev.off()

  if (length(diff_loci_adj_counts) > dim(design)[1]) {
    color_vector <- ifelse(conditions_in == as.character(conditions_in[1]), col_1, col_2)

    hc_row <- hclust(dist(log2(diff_loci_adj_counts + 1), method = "euclidean"), method = "ward.D2")
    hc_col <- hclust(dist(t(log2(diff_loci_adj_counts + 1)), method = "euclidean"), method = "ward.D2")

    pdf(file = paste0(out_path, name, "_edgeR_heatmap_diffs.pdf"), width = 4, height = 6)
    heatmap.2(x = log2(diff_loci_adj_counts + 1),
              Rowv = as.dendrogram(hc_row),
              Colv = as.dendrogram(hc_col),
              dendrogram = "both",
              scale = "row",
              labRow = "",
              col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25)),
              main = name,
              srtCol = 30,
              trace = "none",
              ColSideColors = color_vector)

    heatmap.2(x = log2(diff_loci_adj_counts + 1),
              Colv = FALSE,
              dendrogram = "none",
              scale = "row",
              labRow = "",
              col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25)),
              main = name,
              srtCol = 30,
              trace = "none",
              ColSideColors = color_vector)
    dev.off()

    prcomp_df <- prcomp(as.data.frame(t(diff_loci_adj_counts)), center = TRUE, scale. = TRUE)
    pca_df <- data.frame(
      PC1 = prcomp_df$x[, 1],
      PC2 = prcomp_df$x[, 2],
      PC3 = prcomp_df$x[, 3],
      classification = conditions_in
    )
    var_explained <- (prcomp_df$sdev)^2 / sum((prcomp_df$sdev)^2)

    pdf(file = paste0(out_path, name, "_diff_hits_PCA.pdf"), width = 3, height = 3)
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = classification)) +
      geom_point() +
      scale_color_manual(values = c(col_2, col_1), name = "Condition") +
      xlab(paste0("PC1(", round(var_explained[1] * 100, 2), "%)")) +
      ylab(paste0("PC2(", round(var_explained[2] * 100, 2), "%)")) +
      theme_bw() +
      theme(legend.position = "top")
    plot(p)
    dev.off()

    pdf(file = paste0(out_path, name, "_diff_hits_ROC_based_on_PC.pdf"), width = 4, height = 4)
    roc(response = pca_df$classification, predictor = pca_df$PC1,
        ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2, main = "PC1",
        print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
    roc(response = pca_df$classification, predictor = pca_df$PC2,
        ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2, main = "PC2",
        print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
    dev.off()
  } else {
    print("Less than 2 hits identified. No Clustering or PCA.")
  }

  res_eR_FCpos <- res_eR[which(res_eR$table$logFC > 0), ]
  res_eR_FCneg <- res_eR[which(res_eR$table$logFC < 0), ]
  res_eR_FCpos <- res_eR_FCpos[order(res_eR_FCpos$table$PValue), ]
  res_eR_FCneg <- res_eR_FCneg[order(res_eR_FCneg$table$PValue), ]

  res_eR_FCpos_top <- head(res_eR_FCpos, top_to_select %/% 2)
  res_eR_FCneg_top <- head(res_eR_FCneg, top_to_select %/% 2)

  top_loci_id_list <- c(row.names(res_eR_FCneg_top$table), row.names(res_eR_FCpos_top$table))
  top_loci_adj_counts <- test$fitted.values[top_loci_id_list, ]

  color_vector <- ifelse(conditions_in == as.character(conditions_in[1]), col_1, col_2)

  hc_row <- hclust(dist(log2(top_loci_adj_counts + 1), method = "euclidean"), method = "ward.D2")
  hc_col <- hclust(dist(t(log2(top_loci_adj_counts + 1)), method = "euclidean"), method = "ward.D2")

  pdf(file = paste0(out_path, name, "_edgeR_heatmap_top.pdf"), width = 4, height = 6)
  heatmap.2(x = log2(top_loci_adj_counts + 1),
            Rowv = as.dendrogram(hc_row),
            Colv = as.dendrogram(hc_col),
            dendrogram = "both",
            scale = "row",
            labRow = "",
            col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25)),
            main = name,
            srtCol = 30,
            trace = "none",
            ColSideColors = color_vector)

  heatmap.2(x = log2(top_loci_adj_counts + 1),
            Colv = FALSE,
            dendrogram = "none",
            scale = "row",
            labRow = "",
            col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(25)),
            main = name,
            srtCol = 30,
            trace = "none",
            ColSideColors = color_vector)
  dev.off()

  prcomp_df <- prcomp(as.data.frame(t(top_loci_adj_counts)), center = TRUE, scale. = TRUE)
  pca_df <- data.frame(
    PC1 = prcomp_df$x[, 1],
    PC2 = prcomp_df$x[, 2],
    PC3 = prcomp_df$x[, 3],
    classification = conditions_in
  )
  var_explained <- (prcomp_df$sdev)^2 / sum((prcomp_df$sdev)^2)

  pdf(file = paste0(out_path, name, "_top_hits_PCA.pdf"), width = 3, height = 3)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = classification)) +
    geom_point() +
    scale_color_manual(values = c(col_2, col_1), name = "Condition") +
    xlab(paste0("PC1(", round(var_explained[1] * 100, 2), "%)")) +
    ylab(paste0("PC2(", round(var_explained[2] * 100, 2), "%)")) +
    theme_bw() +
    theme(legend.position = "top")
  plot(p)
  dev.off()

  pdf(file = paste0(out_path, name, "_top_hits_ROC_based_on_PC.pdf"), width = 4, height = 4)
  roc(response = pca_df$classification, predictor = pca_df$PC1,
      ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2, main = "PC1",
      print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
  roc(response = pca_df$classification, predictor = pca_df$PC2,
      ci = TRUE, plot = TRUE, print.auc = TRUE, print.auc.x = 0.6, print.auc.y = 0.2, main = "PC2",
      print.auc.pattern = "AUC: %.3f\nCI: %.3f-%.3f")
  dev.off()
}

# Main function to run the analysis
run_analysis <- function(base_path, out_path, name, conditions_in, mean_reads_per_region_thresh, zeros_per_cond, sig_thresh, fc_thresh, top_to_select, col_1, col_2) {
  data_list <- read_and_preprocess_data(base_path, name, conditions_in, mean_reads_per_region_thresh, zeros_per_cond)
  data_offset <- normalize_data(data_list$data)
  analysis_results <- differential_analysis(data_offset, conditions_in, sig_thresh, fc_thresh, top_to_select, out_path, name)
  generate_plots(analysis_results$log2fc_eR, analysis_results$pval_eR, analysis_results$diff_loci_adj_counts, conditions_in, col_1, col_2, sig_thresh, fc_thresh, analysis_results$design, out_path, name, top_to_select, analysis_results$res_eR, analysis_results$test)
}

# Run the main function
run_analysis(base_path, out_path, name, conditions_in, mean_reads_per_region_thresh, zeros_per_cond, sig_thresh, fc_thresh, top_to_select, col_1, col_2)