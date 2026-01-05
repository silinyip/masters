# Single-cell RNA-seq - differential gene expression analysis - INTEGRATED VERSION
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer
# Based off HBC tutorial https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md


# Clean environment
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)
set.seed(42)


# Load libraries
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(ashr)
library(apeglm)
library(png)
library(ggrepel)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(future)


# Setup parallelization
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)


# Load annotated Seurat object
seurat_annotated <- readRDS("results/Clustering/data/annotation/seurat_annotated.rds")

Idents(seurat_annotated) <- seurat_annotated$celltype_major
seurat_annotated$cluster_id <- Idents(seurat_annotated)

seurat_annotated$isolate <- NA
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196636_S1"))] <- "sensitive"
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196637_S2"))] <- "sensitive"
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196638_S3"))] <- "sensitive"
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196639_R1"))] <- "resistant"
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196640_R2"))] <- "resistant"
seurat_annotated$isolate[which(str_detect(seurat_annotated$sample, "SAMN33196641_R3"))] <- "resistant"

seurat_annotated$sex <- NA
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196636_S1"))] <- "F"
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196637_S2"))] <- "M"
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196638_S3"))] <- "F"
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196639_R1"))] <- "M"
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196640_R2"))] <- "M"
seurat_annotated$sex[which(str_detect(seurat_annotated$sample, "SAMN33196641_R3"))] <- "F"

seurat_annotated$tissue <- NA
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196636_S1"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196637_S2"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196638_S3"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196639_R1"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196640_R2"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN33196641_R3"))] <- "CRC"

seurat_annotated$single_cell_chemistry <- NA
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196636_S1"))] <- "10X_v2"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196637_S2"))] <- "10X_v2"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196638_S3"))] <- "10X_v2"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196639_R1"))] <- "10X_v2"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196640_R2"))] <- "10X_v2"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN33196641_R3"))] <- "10X_v2"

seurat_annotated$isolate <- as.factor(seurat_annotated$isolate)
seurat_annotated$sex <- as.factor(seurat_annotated$sex)
seurat_annotated$tissue <- as.factor(seurat_annotated$tissue)
seurat_annotated$single_cell_chemistry <- as.factor(seurat_annotated$single_cell_chemistry)


# Set column as sample_id
seurat_annotated$sample_id <- factor(seurat_annotated$sample)


# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(seurat_annotated, assay = "RNA", layer = "counts")
metadata <- seurat_annotated@meta.data


# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat_annotated@active.ident)


# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)


# Extract unique names of clusters and samples
cluster_names <- levels(colData(sce)$cluster_id)
sample_names <- levels(colData(sce)$sample_id)

print(paste("Total clusters:", length(cluster_names)))
print(paste("Total samples:", length(sample_names)))


# Aggregate counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 
aggr_counts <- t(aggr_counts)


# Create counts list for each cluster
counts_ls <- list()
for (i in 1:length(cluster_names)) {
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}


# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(sample_id, isolate, tissue, sex, single_cell_chemistry)

metadata <- metadata[!duplicated(metadata), ]
rownames(metadata) <- metadata$sample_id


# Number of cells per sample and cluster
t <- table(colData(sce)$sample_id, colData(sce)$cluster_id)


# Creating metadata list
metadata_ls <- list()
for (i in 1:length(counts_ls)) {
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id <- sub("^[^_]*_", "", df$cluster_sample_id)
  
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  cell_counts <- cell_counts[cell_counts > 0]
  
  sample_order <- match(df$sample_id, names(cell_counts))
  df$cell_count <- cell_counts[sample_order]
  
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  rownames(df) <- df$cluster_sample_id
  
  # Add group_id column for compatibility with the function
  df$group_id <- df$isolate
  
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
}


# Function to run DESeq2 Wald Test and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function
## A is the sample group to compare (e.g. stimulated condition)
## B is the sample group to compare against (base/control level)
## padj_cutoff defines the adjusted p-value cutoff for significance (set to 0.05 by default)

get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
  
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
    return(NULL)
  }
  
  # Check if we have enough samples
  if (nrow(cluster_metadata) < 2) {
    print(paste("Skipping", clustx, "- not enough samples"))
    return(NULL)
  }
  
  # Check if we have both conditions
  group_table <- table(cluster_metadata$group_id)
  if (length(group_table) < 2 || !A %in% names(group_table) || !B %in% names(group_table)) {
    print(paste("Skipping", clustx, "- need both", A, "and", B, "samples"))
    return(NULL)
  }
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "group_id")
  if (!dir.exists("results/DGEA/figures")) { dir.create("results/DGEA/figures", recursive = TRUE) }
  ggsave(paste0("results/DGEA/figures/", clustx, "_specific_PCAplot.png"))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  png(paste0("results/DGEA/figures/", clustx, "_specific_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop = FALSE])
  dev.off()
  
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  png(paste0("results/DGEA/figures/", clustx, "_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  # Check available contrasts
  print("Available contrasts:")
  print(resultsNames(dds))
  
  # Set reference level to B (control/baseline)
  dds$group_id <- relevel(dds$group_id, ref = B)
  
  # Re-run DESeq2 with proper reference level
  dds <- DESeq(dds)
  
  # Get the contrast name - it should be the last one in resultsNames
  contrast_name <- resultsNames(dds)[length(resultsNames(dds))]
  print(paste("Using contrast:", contrast_name))
  
  # Get results using the proper contrast name
  res <- results(dds, name = contrast_name, alpha = padj_cutoff)
  
  # Use apeglm for shrinkage (more robust than ashr for this case)
  # If apeglm is not available, comment out the line below
  tryCatch({
    res <- lfcShrink(dds, coef = contrast_name, res = res, type = "apeglm")
  }, error = function(e) {
    print("apeglm not available, using ashr instead")
    res <<- lfcShrink(dds, coef = contrast_name, res = res, type = "ashr")
  })
  
  # Create a readable contrast name for file naming
  contrast <- paste(c("group_id", A, "vs", B), collapse = "_")
  
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0("results/DGEA/data/", clustx, "_", contrast, "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  ## Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0("results/DGEA/data/", clustx, "_", contrast, "_signif_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  
  # Generate results visualization plots
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  
  ## Determine how many significant genes we actually have
  n_sig_genes <- nrow(sig_res)
  n_genes_to_plot <- min(n_sig_genes, 20)  # Use available genes, max 20
  
  ## Extract top N DEG from sig_res (make sure to order by padj)
  top_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = n_genes_to_plot)
  
  if (length(top_sig_genes) > 0) {
    ## Extract matching normalized count values from matrix
    top_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top_sig_genes, ]
    
    ## Convert wide matrix to long data frame for ggplot2
    top_sig_df <- data.frame(top_sig_counts)
    top_sig_df$gene <- rownames(top_sig_counts)
    
    # Fix column names to be valid R names before melting
    colnames(top_sig_df) <- make.names(colnames(top_sig_df))
    
    # Check if gene column exists before melting
    if (!"gene" %in% colnames(top_sig_df)) {
      print("Warning: gene column missing, skipping visualization")
      return(list(
        cluster = clustx,
        n_genes = nrow(cluster_counts),
        n_samples = ncol(cluster_counts),
        n_significant = nrow(sig_res),
        contrast = contrast
      ))
    }
    
    top_sig_df <- melt(setDT(top_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "cluster_sample_id") %>% 
      data.frame()
    
    ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
    top_sig_df$cluster_sample_id <- gsub("\\.", " ", top_sig_df$cluster_sample_id)
    top_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top_sig_df$cluster_sample_id)
    
    ## Join counts data frame with metadata
    top_sig_df <- plyr::join(top_sig_df, as.data.frame(colData(dds)),
                             by = "cluster_sample_id")
    
    ## Generate plot with dynamic title
    plot_title <- ifelse(n_genes_to_plot == 20, 
                         "Top 20 Significant DE Genes",
                         paste("Top", n_genes_to_plot, "Significant DE Genes"))
    
    ggplot(top_sig_df, aes(y = value, x = group_id, col = group_id)) +
      geom_jitter(height = 0, width = 0.15) +
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("") +  # Remove x-axis title
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),      # Remove x-axis text
            axis.ticks.x = element_blank()) +   # Remove x-axis ticks
      facet_wrap(~ gene)
    
    # Dynamic filename based on actual number of genes plotted
    filename <- paste0("results/DGEA/figures/", clustx, "_", contrast, "_top", n_genes_to_plot, "_DE_genes.png")
    ggsave(filename)
    
    print(paste("Generated plot for", n_genes_to_plot, "genes in cluster", clustx))
    
  } else {
    print(paste("No significant genes found for cluster", clustx, "- skipping visualization"))
  }
  
  # Additional comprehensive plots
  # Create safe filename
  safe_name <- gsub("[^A-Za-z0-9_]", "_", clustx)
  
  # Additional PCA plots with different groupings
  p2 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "group_id")
  ggsave(paste0("results/DGEA/figures/", safe_name, "_isolate_PCAplot.png"), plot = p2)
  
  # If other variables exist in metadata, create PCA plots for them too
  if ("tissue" %in% colnames(cluster_metadata)) {
    p3 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "tissue")
    ggsave(paste0("results/DGEA/figures/", safe_name, "_tissue_PCAplot.png"), plot = p3)
  }
  
  if ("sex" %in% colnames(cluster_metadata)) {
    p4 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "sex")
    ggsave(paste0("results/DGEA/figures/", safe_name, "_sex_PCAplot.png"), plot = p4)
  }
  
  if ("single_cell_chemistry" %in% colnames(cluster_metadata)) {
    p5 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "single_cell_chemistry")
    ggsave(paste0("results/DGEA/figures/", safe_name, "_single_cell_chemistry_PCAplot.png"), plot = p5)
  }
  
  if ("cell_count" %in% colnames(cluster_metadata)) {
    p6 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
    ggsave(paste0("results/DGEA/figures/", safe_name, "_cell_count_PCAplot.png"), plot = p6)
  }
  
  # Enhanced correlation heatmap
  png(paste0("results/DGEA/figures/", safe_name, "_correlation_heatmap.png"),
      width = 800, height = 600)
  pheatmap(rld_cor, 
           annotation = cluster_metadata[, c("group_id"), drop = FALSE], 
           main = paste("Correlation:", clustx))
  dev.off()
  
  # Heatmap of significant genes (select top genes if too many)
  if (nrow(sig_res) > 0) {
    # If more than 50 genes, select top 50 most significant
    if (nrow(sig_res) > 50) {
      top_genes_for_heatmap <- sig_res %>%
        arrange(padj) %>%
        head(50) %>%
        pull(gene)
      print(paste("Using top 50 out of", nrow(sig_res), "significant genes for heatmap"))
    } else {
      top_genes_for_heatmap <- sig_res$gene
      print(paste("Using all", nrow(sig_res), "significant genes for heatmap"))
    }
    
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% top_genes_for_heatmap, ]
    heat_colors <- rev(brewer.pal(11, "PuOr"))
    
    png(paste0("results/DGEA/figures/", safe_name, "_significant_genes_heatmap.png"), 
        width = 800, height = 600)
    pheatmap(sig_counts, 
             color = heat_colors, 
             cluster_rows = TRUE, 
             show_rownames = TRUE,
             annotation = cluster_metadata[, c("group_id"), drop = FALSE], 
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 8,
             main = paste("Top Significant Genes:", clustx))
    dev.off()
  }
  
  # Volcano plot
  padj_cutoff_volcano <- 0.05
  log2fc_cutoff <- 0.58
  
  # Label only the most significant genes
  res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
    mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
  
  # Get top genes to label (most significant + highest fold change)
  top_genes_to_label <- res_table_thres %>%
    filter(threshold == TRUE) %>%  # Only significant genes
    arrange(padj) %>%              # Order by significance
    head(10)                       # Top 10 most significant
  
  p_volcano <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    geom_text_repel(data = top_genes_to_label,
                    aes(x = log2FoldChange, y = -log10(padj), label = gene),
                    color = "black",
                    fontsize = 3,
                    max.overlaps = 20,
                    box.padding = 0.3,
                    point.padding = 0.3) +
    ggtitle(paste("Volcano plot:", clustx)) +
    xlab("log2 fold change") +
    xlim(-4.5, 12) +
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0, 250)) +
    scale_color_manual(values = c("grey60", "red3")) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.3), hjust = 0.5),
          axis.title = element_text(size = rel(1.15)))
  
  ggsave(paste0("results/DGEA/figures/", safe_name, "_volcano.png"), 
         plot = p_volcano, width = 10, height = 8)
  
  # Count up and down regulated genes for summary
  n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>% nrow()
  n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>% nrow()
  
  print(paste("Significant genes:", nrow(sig_res), 
              "Up:", n_sig_up, "Down:", n_sig_dn))
  
  # Return summary information
  return(list(
    cluster = clustx,
    n_genes = nrow(cluster_counts),
    n_samples = ncol(cluster_counts),
    n_significant = nrow(sig_res),
    n_upregulated = n_sig_up,
    n_downregulated = n_sig_dn,
    contrast = contrast
  ))
}


# Run the script on all clusters comparing resistant condition relative to sensitive condition
results_summary <- map(cluster_names, get_dds_resultsAvsB, A = "resistant", B = "sensitive")

# Convert results to data frame
summary_df <- do.call(rbind, lapply(results_summary[!sapply(results_summary, is.null)], data.frame))

# Save summary
write.csv(summary_df, "results/DGEA/data/function_analysis_summary.csv", 
          quote = FALSE, row.names = FALSE)

