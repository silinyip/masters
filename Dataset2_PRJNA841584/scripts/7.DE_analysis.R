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

seurat_annotated$tissue <- NA
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600722_LiverMet_44"))] <- "LiverMet"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600723_Liver_44"))] <- "Liver"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600724_Tumour_44"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600725_Tumoroid_28"))] <- "Tumoroid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600726_Colonoid_28"))] <- "Colonoid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600727_Tumour_28"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600728_Colonoid_26"))] <- "Colonoid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600729_Tumoroid_24"))] <- "Tumoroid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600730_Colonoid_24"))] <- "Colonoid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600731_LiverMet_27"))] <- "LiverMet"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600732_Tumour_27"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600733_Tumour_26"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600734_Tumour_24"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600735_Tumoroid_08"))] <- "Tumoroid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600736_Colonoid_08"))] <- "Colonoid"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600737_LiverMet_09"))] <- "LiverMet"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600738_Liver_09"))] <- "Liver"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600739_Tumour_08"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600740_Tumour_09"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600741_Tumour_07"))] <- "CRC"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600742_Colon_08"))] <- "Colon"
seurat_annotated$tissue[which(str_detect(seurat_annotated$sample, "SAMN28600743_Colon_09"))] <- "Colon"
seurat_annotated$tissue[which(str_detect(seurat_annotated$cells, "SAMN28600744_Colon_07"))] <- "Colon"

seurat_annotated$single_cell_chemistry <- NA
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600722_LiverMet_44"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600723_Liver_44"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600724_Tumour_44"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600725_Tumoroid_28"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600726_Colonoid_28"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600727_Tumour_28"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600728_Colonoid_26"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600729_Tumoroid_24"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600730_Colonoid_24"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600731_LiverMet_27"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600732_Tumour_27"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600733_Tumour_26"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600734_Tumour_24"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600735_Tumoroid_08"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600736_Colonoid_08"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600737_LiverMet_09"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600738_Liver_09"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600739_Tumour_08"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600740_Tumour_09"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600741_Tumour_07"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600742_Colon_08"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$sample, "SAMN28600743_Colon_09"))] <- "10X_v3"
seurat_annotated$single_cell_chemistry[which(str_detect(seurat_annotated$cells, "SAMN28600744_Colon_07"))] <- "10X_v3"

# Create sample_type column - broader categorization
seurat_annotated$sample_type <- NA

# Primary tumours (in vivo samples)
seurat_annotated$sample_type[seurat_annotated$tissue %in% c("CRC", "LiverMet", "Liver", "Colon")] <- "primary_tumours"

# Organoids (in vitro samples)  
seurat_annotated$sample_type[seurat_annotated$tissue %in% c("Tumoroid", "Colonoid")] <- "organoids"

# Convert to factor
seurat_annotated$sample_type <- as.factor(seurat_annotated$sample_type)

# Check the distribution
print("Sample type distribution:")
print(table(seurat_annotated$sample_type, seurat_annotated$tissue))

# Verify no missing values
print("Missing values in sample_type:")
print(sum(is.na(seurat_annotated$sample_type)))

# Convert existing columns to factors
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


# Extract sample-level variables - UPDATED TO INCLUDE SAMPLE_TYPE
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(sample_id, tissue, sample_type, single_cell_chemistry)

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
  df$group_id <- df$tissue
  
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
  
  # DEBUG: Add this line to see available tissue types
  print(paste("Available tissue types in", clustx, ":", paste(names(group_table), collapse = ", ")))
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "group_id")
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
  
  # Debug: print available levels before releveling
  print(paste("Available tissue types in", clustx, ":", paste(levels(dds$group_id), collapse = ", ")))
  print(paste("Trying to set reference level to:", B))
  
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
    
    ## Generate plot with dynamic title and increased font sizes
    plot_title <- ifelse(n_genes_to_plot == 20, 
                         "Top 20 Significant DE Genes",
                         paste("Top", n_genes_to_plot, "Significant DE Genes"))
    
    ggplot(top_sig_df, aes(y = value, x = group_id, col = group_id)) +
      geom_jitter(height = 0, width = 0.15) +
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("") +  # Remove x-axis title
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),     # Increased title size
            axis.text.x = element_blank(),                         # Remove x-axis text
            axis.ticks.x = element_blank(),                        # Remove x-axis ticks
            axis.title.y = element_text(size = 14),                # Increased y-axis title size
            axis.text.y = element_text(size = 12),                 # Increased y-axis text size
            strip.text = element_text(size = 12),                  # Increased facet label size
            legend.title = element_text(size = 14),                # Increased legend title size
            legend.text = element_text(size = 12)) +               # Increased legend text size
      facet_wrap(~ gene)
    
    # Dynamic filename based on actual number of genes plotted
    filename <- paste0("results/DGEA/figures/", clustx, "_", contrast, "_top", n_genes_to_plot, "_DE_genes.png")
    ggsave(filename, width = 14, height = 10)  # Increased plot size
    
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
  
  if ("sample_type" %in% colnames(cluster_metadata)) {
    p4 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_type")
    ggsave(paste0("results/DGEA/figures/", safe_name, "_sample_type_PCAplot.png"), plot = p4)
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
  
  # Heatmap of significant genes (only if not too many genes) - IMPROVED LAYOUT
  if (nrow(sig_res) > 0 && nrow(sig_res) <= 50) {
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
    heat_colors <- rev(brewer.pal(11, "PuOr"))
    
    # Calculate dynamic height based on number of genes to prevent squashing
    n_genes_heatmap <- nrow(sig_counts)
    heatmap_height <- max(600, n_genes_heatmap * 25)  # Minimum 600px, 25px per gene
    heatmap_width <- max(800, ncol(sig_counts) * 50)   # Minimum 800px, 50px per sample
    
    png(paste0("results/DGEA/figures/", safe_name, "_significant_genes_heatmap.png"), 
        width = heatmap_width, height = heatmap_height)
    pheatmap(sig_counts, 
             color = heat_colors, 
             cluster_rows = TRUE, 
             show_rownames = TRUE,
             annotation = cluster_metadata[, c("group_id"), drop = FALSE], 
             border_color = NA, 
             fontsize = 12,                    # Increased base font size
             scale = "row", 
             fontsize_row = 10,                # Increased row font size (gene names)
             fontsize_col = 10,                # Increased column font size (sample names)
             main = paste("Significant genes:", clustx),
             cellwidth = 20,                   # Fixed cell width to prevent squashing
             cellheight = 15)                  # Fixed cell height to prevent squashing
    dev.off()
  }
  
  # Volcano plot with gene labels - IMPROVED FONT SIZES
  padj_cutoff_volcano <- 0.05
  log2fc_cutoff <- 0.58
  
  res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
    mutate(threshold = padj < padj_cutoff_volcano & abs(log2FoldChange) >= log2fc_cutoff)
  
  # Identify top genes for labeling (most significant and high fold change)
  top_genes_to_label <- res_table_thres %>%
    filter(threshold == TRUE) %>%
    arrange(padj) %>%
    head(15) %>%  # Label top 15 significant genes
    pull(gene)
  
  # Create a column for gene labels (only label top genes)
  res_table_thres <- res_table_thres %>%
    mutate(gene_label = ifelse(gene %in% top_genes_to_label, gene, ""))
  
  p_volcano <- ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = threshold), alpha = 0.7) +
    geom_text_repel(aes(label = gene_label), 
                    size = 4,                    # Increased label size from 3 to 4
                    max.overlaps = 20,
                    box.padding = 0.3,
                    point.padding = 0.3,
                    segment.color = "grey50",
                    segment.size = 0.2) +
    ggtitle(paste("Volcano plot:", clustx)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    scale_color_manual(values = c("grey60", "red3"), 
                       name = "Significant",
                       labels = c("No", "Yes")) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 18, hjust = 0.5),      # Increased title size
          axis.title = element_text(size = 16),                   # Increased axis title size
          axis.text = element_text(size = 14),                    # Increased axis text size
          legend.title = element_text(size = 14),                 # Increased legend title size
          legend.text = element_text(size = 12),                  # Increased legend text size
          panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA)) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff_volcano), 
               linetype = "dashed", alpha = 0.5)
  
  ggsave(paste0("results/DGEA/figures/", safe_name, "_volcano.png"), 
         plot = p_volcano, width = 14, height = 12)              # Increased plot size
  
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


# Run the script on all clusters comparing different tissue types
# Check tissue distribution across clusters first
tissue_cluster_summary <- colData(sce) %>%
  as.data.frame() %>%
  group_by(cluster_id, tissue) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = n_cells, values_fill = 0)

print("Tissue distribution across clusters:")
print(tissue_cluster_summary)

# Check sample_type distribution across clusters
sample_type_cluster_summary <- colData(sce) %>%
  as.data.frame() %>%
  group_by(cluster_id, sample_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample_type, values_from = n_cells, values_fill = 0)

print("Sample type distribution across clusters:")
print(sample_type_cluster_summary)

# Find clusters that have both tissues for each comparison
clusters_CRC_Colon <- tissue_cluster_summary %>%
  filter(CRC > 0 & Colon > 0) %>%
  pull(cluster_id)

clusters_LiverMet_CRC <- tissue_cluster_summary %>%
  filter(LiverMet > 0 & CRC > 0) %>%
  pull(cluster_id)

clusters_Tumoroid_Colonoid <- tissue_cluster_summary %>%
  filter(Tumoroid > 0 & Colonoid > 0) %>%
  pull(cluster_id)

clusters_CRC_Tumoroid <- tissue_cluster_summary %>%
  filter(CRC > 0 & Tumoroid > 0) %>%
  pull(cluster_id)

clusters_Colon_Colonoid <- tissue_cluster_summary %>%
  filter(Colon > 0 & Colonoid > 0) %>%
  pull(cluster_id)

# Find clusters that have both sample types for broad comparison
clusters_primary_organoids <- sample_type_cluster_summary %>%
  filter(primary_tumours > 0 & organoids > 0) %>%
  pull(cluster_id)

print("Clusters available for each comparison:")
print(paste("CRC vs Colon:", length(clusters_CRC_Colon), "clusters"))
print(paste("LiverMet vs CRC:", length(clusters_LiverMet_CRC), "clusters"))
print(paste("Tumoroid vs Colonoid:", length(clusters_Tumoroid_Colonoid), "clusters"))
print(paste("CRC vs Tumoroid:", length(clusters_CRC_Tumoroid), "clusters"))
print(paste("Colon vs Colonoid:", length(clusters_Colon_Colonoid), "clusters"))
print(paste("Primary tumours vs Organoids:", length(clusters_primary_organoids), "clusters"))


# Function to run DESeq2 with sample_type as grouping variable
get_dds_resultsAvsB_sample_type <- function(clustx, A, B, padj_cutoff = 0.05) {
  
  print(paste("Analyzing", clustx, "for sample_type comparison")) # useful for debugging
  
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
  
  # Check if we have both sample types
  sample_type_table <- table(cluster_metadata$sample_type)
  if (length(sample_type_table) < 2 || !A %in% names(sample_type_table) || !B %in% names(sample_type_table)) {
    print(paste("Skipping", clustx, "- need both", A, "and", B, "sample types"))
    return(NULL)
  }
  
  # DEBUG: Add this line to see available sample types
  print(paste("Available sample types in", clustx, ":", paste(names(sample_type_table), collapse = ", ")))
  
  # Use sample_type for DESeq2 design
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sample_type)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "sample_type")
  ggsave(paste0("results/DGEA/figures/", clustx, "_sample_type_PCAplot.png"))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  png(paste0("results/DGEA/figures/", clustx, "_sample_type_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("sample_type"), drop = FALSE])
  dev.off()
  
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  png(paste0("results/DGEA/figures/", clustx, "_sample_type_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  # Check available contrasts
  print("Available contrasts:")
  print(resultsNames(dds))
  
  # Debug: print available levels before releveling
  print(paste("Available sample types in", clustx, ":", paste(levels(dds$sample_type), collapse = ", ")))
  print(paste("Trying to set reference level to:", B))
  
  # Set reference level to B (control/baseline)
  dds$sample_type <- relevel(dds$sample_type, ref = B)
  
  # Re-run DESeq2 with proper reference level
  dds <- DESeq(dds)
  
  # Get the contrast name - it should be the last one in resultsNames
  contrast_name <- resultsNames(dds)[length(resultsNames(dds))]
  print(paste("Using contrast:", contrast_name))
  
  # Get results using the proper contrast name
  res <- results(dds, name = contrast_name, alpha = padj_cutoff)
  
  # Use apeglm for shrinkage (more robust than ashr for this case)
  tryCatch({
    res <- lfcShrink(dds, coef = contrast_name, res = res, type = "apeglm")
  }, error = function(e) {
    print("apeglm not available, using ashr instead")
    res <<- lfcShrink(dds, coef = contrast_name, res = res, type = "ashr")
  })
  
  # Create a readable contrast name for file naming
  contrast <- paste(c("sample_type", A, "vs", B), collapse = "_")
  
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
                         paste("Top 20 Significant DE Genes:", A, "vs", B),
                         paste("Top", n_genes_to_plot, "Significant DE Genes:", A, "vs", B))
    
    ggplot(top_sig_df, aes(y = value, x = sample_type, col = sample_type)) +
      geom_jitter(height = 0, width = 0.15) +
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("") +  # Remove x-axis title
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),     # Increased title size
            axis.text.x = element_blank(),                         # Remove x-axis text
            axis.ticks.x = element_blank(),                        # Remove x-axis ticks
            axis.title.y = element_text(size = 14),                # Increased y-axis title size
            axis.text.y = element_text(size = 12),                 # Increased y-axis text size
            strip.text = element_text(size = 12),                  # Increased facet label size
            legend.title = element_text(size = 14),                # Increased legend title size
            legend.text = element_text(size = 12)) +               # Increased legend text size
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
  
  # Enhanced correlation heatmap
  png(paste0("results/DGEA/figures/", safe_name, "_sample_type_correlation_heatmap.png"),
      width = 800, height = 600)
  pheatmap(rld_cor, 
           annotation = cluster_metadata[, c("sample_type"), drop = FALSE], 
           main = paste("Sample Type Correlation:", clustx))
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
    
    # Calculate dynamic height based on number of genes to prevent squashing
    n_genes_heatmap <- nrow(sig_counts)
    heatmap_height <- max(600, n_genes_heatmap * 25)  # Minimum 600px, 25px per gene
    heatmap_width <- max(800, ncol(sig_counts) * 50)   # Minimum 800px, 50px per sample
    
    png(paste0("results/DGEA/figures/", safe_name, "_sample_type_significant_genes_heatmap.png"), 
        width = heatmap_width, height = heatmap_height)
    pheatmap(sig_counts, 
             color = heat_colors, 
             cluster_rows = TRUE, 
             show_rownames = TRUE,
             annotation = cluster_metadata[, c("sample_type"), drop = FALSE], 
             fontsize = 12,                    # Increased base font size
             scale = "row", 
             fontsize_row = 10,                # Increased row font size (gene names)
             fontsize_col = 10,                # Increased column font size (sample names)
             cellwidth = 20,                   # Fixed cell width to prevent squashing
             cellheight = 15,                  # Fixed cell height to prevent squashing
             main = paste("Top Significant Genes (Sample Type):", clustx))
    dev.off()
  }
  
  # Volcano plot with gene labels
  padj_cutoff_volcano <- 0.05
  log2fc_cutoff <- 0.58
  
  res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
    mutate(threshold = padj < padj_cutoff_volcano & abs(log2FoldChange) >= log2fc_cutoff)
  
  # Identify top genes for labeling (most significant and high fold change)
  top_genes_to_label <- res_table_thres %>%
    filter(threshold == TRUE) %>%
    arrange(padj) %>%
    head(15) %>%  # Label top 15 significant genes
    pull(gene)
  
  # Create a column for gene labels (only label top genes)
  res_table_thres <- res_table_thres %>%
    mutate(gene_label = ifelse(gene %in% top_genes_to_label, gene, ""))
  
  p_volcano <- ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = threshold), alpha = 0.7) +
    geom_text_repel(aes(label = gene_label), 
                    size = 4, 
                    max.overlaps = 20,
                    box.padding = 0.3,
                    point.padding = 0.3,
                    segment.color = "grey50",
                    segment.size = 0.2) +
    ggtitle(paste("Volcano plot (Sample Type):", clustx, "-", A, "vs", B)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    scale_color_manual(values = c("grey60", "red3"), 
                       name = "Significant",
                       labels = c("No", "Yes")) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 18, hjust = 0.5),      # Increased title size
          axis.title = element_text(size = 16),                   # Increased axis title size
          axis.text = element_text(size = 14),                    # Increased axis text size
          legend.title = element_text(size = 14),                 # Increased legend title size
          legend.text = element_text(size = 12),                  # Increased legend text size
          panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA)) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff_volcano), 
               linetype = "dashed", alpha = 0.5)
  
  ggsave(paste0("results/DGEA/figures/", safe_name, "_sample_type_volcano.png"), 
         plot = p_volcano, width = 12, height = 10)
  
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


# Run analyses only on clusters that have both conditions
# Tumour vs Normal tissue
results_CRC_vs_Colon <- map(clusters_CRC_Colon, get_dds_resultsAvsB, A = "CRC", B = "Colon")

# Metastasis vs Primary Tumour  
results_LiverMet_vs_CRC <- map(clusters_LiverMet_CRC, get_dds_resultsAvsB, A = "LiverMet", B = "CRC")

# Cancerous Organoids vs Normal Organoids
results_Tumoroid_vs_Colonoid <- map(clusters_Tumoroid_Colonoid, get_dds_resultsAvsB, A = "Tumoroid", B = "Colonoid")

# In vivo vs in vitro (specific tissues)
results_CRC_vs_Tumoroid <- map(clusters_CRC_Tumoroid, get_dds_resultsAvsB, A = "CRC", B = "Tumoroid")
results_Colon_vs_Colonoid <- map(clusters_Colon_Colonoid, get_dds_resultsAvsB, A = "Colon", B = "Colonoid")

# Primary tumours vs Organoids (broad comparison)
results_primary_vs_organoids <- map(clusters_primary_organoids, get_dds_resultsAvsB_sample_type, A = "primary_tumours", B = "organoids")


# Combine all results
all_results <- c(results_CRC_vs_Colon, results_LiverMet_vs_CRC, results_Tumoroid_vs_Colonoid, 
                 results_CRC_vs_Tumoroid, results_Colon_vs_Colonoid, results_primary_vs_organoids)


# Remove NULL results and convert to data frames
valid_results <- list()
for (i in seq_along(all_results)) {
  if (!is.null(all_results[[i]])) {
    # Convert to data frame and ensure consistent structure
    df <- data.frame(all_results[[i]])
    # Add missing columns with NA if needed
    expected_cols <- c("cluster", "n_genes", "n_samples", "n_significant", 
                       "n_upregulated", "n_downregulated", "contrast")
    missing_cols <- setdiff(expected_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }
    # Reorder columns to match expected order
    df <- df[expected_cols]
    valid_results[[length(valid_results) + 1]] <- df
  }
}


# Convert results to data frame
if (length(valid_results) > 0) {
  summary_df <- do.call(rbind, valid_results)
  rownames(summary_df) <- NULL
} else {
  print("No valid results to combine")
  summary_df <- data.frame(cluster = character(0), 
                           n_genes = numeric(0),
                           n_samples = numeric(0),
                           n_significant = numeric(0),
                           n_upregulated = numeric(0),
                           n_downregulated = numeric(0),
                           contrast = character(0))
}

print(paste("Total valid results:", nrow(summary_df)))

# Save summary
write.csv(summary_df, "results/DGEA/data/function_analysis_summary.csv", 
          quote = FALSE, row.names = FALSE)

# Print summary by contrast type
print("Summary by contrast type:")
contrast_summary <- summary_df %>%
  group_by(contrast) %>%
  summarise(
    n_clusters = n(),
    total_significant_genes = sum(n_significant, na.rm = TRUE),
    avg_significant_per_cluster = mean(n_significant, na.rm = TRUE),
    .groups = "drop"
  )
print(contrast_summary)



# Define marker genes of interest
epithelial_markers <- c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19")
goblet_markers <- c("MUC2", "TFF3", "FCGBP", "SPDEF")
enterocyte_markers <- c("VILLIN1", "VIL1", "ALPI", "SI", "APOA1", "APOA4")
stem_markers <- c("LGR5", "OLFM4", "ASCL2", "SMOC2", "AXIN2")
proliferation_markers <- c("MKI67", "TOP2A", "PCNA", "CDK1")
paneth_markers <- c("LYZ", "DEFA5", "DEFA6", "REG3A")
enteroendocrine_markers <- c("CHGA", "CHGB", "TAC1", "TPH1", "NEUROD1")
housekeeping <- c("ACTB", "GAPDH", "B2M", "HPRT1")

all_markers <- unique(c(epithelial_markers, goblet_markers, enterocyte_markers, 
                        stem_markers, proliferation_markers, paneth_markers, 
                        enteroendocrine_markers, housekeeping))

# Create a list to organize markers by category
marker_categories <- list(
  "Epithelial" = epithelial_markers,
  "Goblet" = goblet_markers,
  "Enterocyte" = enterocyte_markers,
  "Stem" = stem_markers,
  "Proliferation" = proliferation_markers,
  "Paneth" = paneth_markers,
  "Enteroendocrine" = enteroendocrine_markers,
  "Housekeeping" = housekeeping
)

# Filter markers that exist in your dataset
available_markers <- all_markers[all_markers %in% rownames(seurat_annotated)]
missing_markers <- all_markers[!all_markers %in% rownames(seurat_annotated)]

print("=== MARKER AVAILABILITY ===")
print(paste("Available markers:", length(available_markers), "out of", length(all_markers)))
if (length(missing_markers) > 0) {
  print("Missing markers:")
  print(missing_markers)
}

# Create directory for marker analysis results
dir.create("results/DGEA/marker_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/DGEA/marker_analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/DGEA/marker_analysis/data", showWarnings = FALSE, recursive = TRUE)

# 1. EXTRACT NORMALIZED EXPRESSION FOR MARKERS
# Get normalized counts for available markers
marker_expression <- GetAssayData(seurat_annotated, assay = "RNA", layer = "data")[available_markers, ]

# Create metadata frame
marker_metadata <- seurat_annotated@meta.data %>%
  dplyr::select(sample_id, tissue, sample_type, celltype_major)

# Combine expression with metadata
marker_df <- as.data.frame(t(as.matrix(marker_expression)))
marker_df$cell_id <- rownames(marker_df)
marker_df <- merge(marker_df, marker_metadata, by.x = "cell_id", by.y = "row.names")

# 2. CALCULATE AVERAGE EXPRESSION BY TISSUE AND CELL TYPE
avg_expression <- marker_df %>%
  pivot_longer(cols = all_of(available_markers), 
               names_to = "gene", 
               values_to = "expression") %>%
  group_by(tissue, celltype_major, gene) %>%
  summarise(
    mean_expr = mean(expression),
    median_expr = median(expression),
    pct_expressing = sum(expression > 0) / n() * 100,
    n_cells = n(),
    .groups = "drop"
  )

# Save detailed results
write.csv(avg_expression, 
          "results/DGEA/marker_analysis/data/marker_expression_by_tissue_celltype.csv",
          row.names = FALSE)

# 3. COMPARE ORGANOIDS VS PRIMARY TISSUES
comparison_summary <- avg_expression %>%
  mutate(tissue_category = case_when(
    tissue %in% c("Tumoroid", "Colonoid") ~ "Organoid",
    tissue %in% c("CRC", "Colon", "LiverMet", "Liver") ~ "Primary",
    TRUE ~ "Other"
  )) %>%
  group_by(gene, celltype_major, tissue_category) %>%
  summarise(
    avg_mean_expr = mean(mean_expr),
    avg_pct_expressing = mean(pct_expressing),
    total_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = tissue_category,
    values_from = c(avg_mean_expr, avg_pct_expressing, total_cells)
  ) %>%
  mutate(
    expr_ratio = avg_mean_expr_Organoid / avg_mean_expr_Primary,
    pct_diff = avg_pct_expressing_Organoid - avg_pct_expressing_Primary,
    captured = case_when(
      expr_ratio >= 0.5 & expr_ratio <= 2 & abs(pct_diff) <= 30 ~ "Well Captured",
      expr_ratio >= 0.3 & expr_ratio <= 3 & abs(pct_diff) <= 50 ~ "Moderately Captured",
      TRUE ~ "Poorly Captured"
    )
  )

write.csv(comparison_summary, 
          "results/DGEA/marker_analysis/data/organoid_vs_primary_marker_summary.csv",
          row.names = FALSE)

# 4. VISUALIZATION: HEATMAP OF MARKER EXPRESSION
# Create matrix for heatmap
heatmap_data <- avg_expression %>%
  dplyr::select(tissue, celltype_major, gene, mean_expr) %>%
  unite("sample", tissue, celltype_major, sep = "_") %>%
  pivot_wider(names_from = sample, values_from = mean_expr, values_fill = 0) %>%
  column_to_rownames("gene")

# Create annotation for heatmap
annotation_col <- data.frame(
  Tissue = str_extract(colnames(heatmap_data), "^[^_]+"),
  CellType = str_extract(colnames(heatmap_data), "[^_]+$")
)
rownames(annotation_col) <- colnames(heatmap_data)
annotation_col$Type <- ifelse(annotation_col$Tissue %in% c("Tumoroid", "Colonoid"), 
                              "Organoid", "Primary")

# Define colors
ann_colors <- list(
  Type = c("Organoid" = "#E64B35FF", "Primary" = "#4DBBD5FF"),
  Tissue = setNames(
    rainbow(length(unique(annotation_col$Tissue))),
    unique(annotation_col$Tissue)
  )
)

png("results/DGEA/marker_analysis/figures/marker_expression_heatmap.png",
    width = 16, height = 12, units = "in", res = 300)
pheatmap(heatmap_data,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         fontsize_row = 14,
         fontsize_col = 12,
         fontsize = 14,
         main = "Marker Gene Expression: Organoids vs Primary Tissues",
         cellwidth = 15,
         cellheight = 15)
dev.off()

# 5. VIOLIN PLOTS FOR KEY MARKERS
for (marker in c("EPCAM", "MUC2", "LGR5", "MKI67", "CHGA", "LYZ", "ACTB")) {
  if (marker %in% available_markers) {
    p <- ggplot(marker_df, aes(x = tissue, y = .data[[marker]], fill = sample_type)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      scale_y_continuous(trans = 'log1p') +
      labs(title = paste("Expression of", marker),
           y = "Normalized Expression (log1p)",
           x = "Tissue Type") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 14)) +
      facet_wrap(~ celltype_major, scales = "free_y")
    
    ggsave(paste0("results/DGEA/marker_analysis/figures/", marker, "_violin_plot.png"),
           plot = p, width = 14, height = 10)
  }
}

# 6. SUMMARY STATISTICS BY MARKER CATEGORY
category_summary <- list()
for (cat_name in names(marker_categories)) {
  cat_markers <- marker_categories[[cat_name]]
  cat_markers_available <- cat_markers[cat_markers %in% available_markers]
  
  if (length(cat_markers_available) > 0) {
    cat_summary <- comparison_summary %>%
      filter(gene %in% cat_markers_available) %>%
      group_by(captured) %>%
      summarise(
        n_markers = n(),
        markers = paste(gene, collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(category = cat_name)
    
    category_summary[[cat_name]] <- cat_summary
  }
}

category_summary_df <- bind_rows(category_summary)
write.csv(category_summary_df,
          "results/DGEA/marker_analysis/data/marker_capture_by_category.csv",
          row.names = FALSE)

# 7. CREATE SUMMARY REPORT
print("\n=== MARKER CAPTURE SUMMARY ===")
print(table(comparison_summary$captured))

print("\n=== BY MARKER CATEGORY ===")
print(category_summary_df)

print("\n=== POORLY CAPTURED MARKERS ===")
poorly_captured <- comparison_summary %>%
  filter(captured == "Poorly Captured") %>%
  arrange(expr_ratio)
print(poorly_captured)

# 8. CORRELATION ANALYSIS
# Calculate correlation between organoid and primary expression
correlation_data <- comparison_summary %>%
  filter(!is.na(avg_mean_expr_Organoid) & !is.na(avg_mean_expr_Primary)) %>%
  group_by(celltype_major) %>%
  summarise(
    correlation = cor(avg_mean_expr_Organoid, avg_mean_expr_Primary, 
                      use = "complete.obs"),
    n_markers = n(),
    .groups = "drop"
  )

write.csv(correlation_data,
          "results/DGEA/marker_analysis/data/organoid_primary_correlation.csv",
          row.names = FALSE)

print("\n=== CORRELATION BY CELL TYPE ===")
print(correlation_data)

# 9. SCATTER PLOT: ORGANOID VS PRIMARY EXPRESSION
# Load DESeq2 results to identify significant DEGs
deg_files <- list.files("results/DGEA/data/Primary_vs_Organoid/", 
                        pattern = "sample_type_primary_tumours_vs_organoids_signif_genes.csv",
                        full.names = TRUE)

significant_degs <- c()
if (length(deg_files) > 0) {
  for (file in deg_files) {
    deg_data <- read.csv(file)
    significant_degs <- c(significant_degs, deg_data$gene)
  }
  significant_degs <- unique(significant_degs)
  print(paste("Found", length(significant_degs), "significant DEGs from DESeq2 analysis"))
} else {
  print("Warning: No DESeq2 results found. Make sure to run the primary_tumours vs organoids comparison first.")
}

# Filter comparison_summary to only include markers that are also significant DEGs
comparison_summary_with_sig <- comparison_summary %>%
  mutate(is_significant_deg = gene %in% significant_degs)

# Create scatter plot - only label markers that are BOTH poorly captured AND significant DEGs
p_scatter <- ggplot(comparison_summary_with_sig %>% 
                      filter(!is.na(avg_mean_expr_Organoid) & 
                               !is.na(avg_mean_expr_Primary)),
                    aes(x = avg_mean_expr_Primary, 
                        y = avg_mean_expr_Organoid,
                        color = captured)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(data = . %>% filter(is_significant_deg == TRUE),
                  aes(label = gene), 
                  size = 5, 
                  max.overlaps = 25,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.size = 0.4,
                  fontface = "bold",
                  min.segment.length = 0) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Well Captured" = "green3",
                                "Moderately Captured" = "orange",
                                "Poorly Captured" = "red")) +
  labs(title = "Marker Expression: Organoids vs Primary Tissues\n(Labels show DESeq2-significant genes)",
       x = "Average Expression in Primary Tissues (log10)",
       y = "Average Expression in Organoids (log10)",
       color = "Capture Quality") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_wrap(~ celltype_major, scales = "free")

ggsave("results/DGEA/marker_analysis/figures/organoid_vs_primary_scatter.png",
       plot = p_scatter, width = 28, height = 20, dpi = 600)

# Save summary with DEG information
comparison_summary_with_sig_export <- comparison_summary_with_sig %>%
  arrange(desc(is_significant_deg), captured)

write.csv(comparison_summary_with_sig_export,
          "results/DGEA/marker_analysis/data/markers_with_DEG_status.csv",
          row.names = FALSE)

print(paste("\nMarkers that are also significant DEGs:", 
            sum(comparison_summary_with_sig$is_significant_deg, na.rm = TRUE)))

print("\n=== ANALYSIS COMPLETE ===")
print("Results saved to results/DGEA/marker_analysis/")