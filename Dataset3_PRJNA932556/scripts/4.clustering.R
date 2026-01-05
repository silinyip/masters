# Single-cell RNA-seq - clustering
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer


# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(42)


# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(scales)
library(metap)
library(RPresto)
library(future)


# Setup parallelization to speed up computational tasks
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)


# Load seurat_integrated object
seurat_integrated <- readRDS("results/Normalisation/SCT/data/integrated_seurat.rds")


# Explore heatmap of PCs
PC_heatmap <- DimHeatmap(seurat_integrated, 
                         dims = 1:9, 
                         cells = 500, 
                         balanced = TRUE)
png("results/Clustering/figures/integrated/QC/PC_heatmap.png", width = 2000, height = 2000, res = 300)
print(DimHeatmap(seurat_integrated, dims = 1:9, cells = 500, balanced = TRUE))
dev.off()


# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)


# Plot the elbow plot
elbow_plot <- ElbowPlot(object = seurat_integrated, ndims = 40) + theme_bw()
ggsave("results/Clustering/figures/integrated/QC/elbow_plot.png", plot = elbow_plot)


# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)


# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


# Explore resolutions
seurat_integrated@meta.data %>% 
  View()


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

saveRDS(seurat_integrated, file = "results/Clustering/data/clustering/seurat_integrated.rds")


# Plot the UMAP
# Reorder factor levels of identities
seurat_integrated$seurat_clusters <- factor(seurat_integrated$seurat_clusters,
                                            levels = sort(as.numeric(levels(seurat_integrated$seurat_clusters))))
# Now re-plot
UMAP_clusters_seurat <- DimPlot(seurat_integrated,
                                reduction = "umap",
                                label = TRUE,
                                label.size = 6) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/UMAP_clusters_seurat.png", plot = UMAP_clusters_seurat, width = 10, height = 10, limitsize = FALSE)


# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)


# Barplot of number of cells per cluster by sample
cells_per_cluster_barplot <- ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/cells_per_cluster_barplot.png", plot = cells_per_cluster_barplot, width = 55, height = 10, limitsize = FALSE)


# UMAP of cells in each cluster by sample
UMAP_cluster_split_by_sample <- DimPlot(seurat_integrated, label = TRUE, split.by = "sample")  + NoLegend()
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/UMAP_cluster_split_by_sample.png", plot = UMAP_cluster_split_by_sample, width = 40, height = 10, limitsize = FALSE)


# UMAP of cells in each cluster by cell cycle phase
UMAP_cell_cycle_phase <- DimPlot(seurat_integrated, label = TRUE, split.by = "Phase")  + NoLegend()
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/UMAP_cell_cycle_phase.png", plot = UMAP_cell_cycle_phase, width = 40, height = 20, limitsize = FALSE)


# Barplot of proportion of cells in each cluster by sample
cell_proportion_barplot <- ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), position=position_fill())
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/cell_proportion_barplot.png", plot = cell_proportion_barplot, width = 20, height = 10, limitsize = FALSE)


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

feature_plot <- FeaturePlot(seurat_integrated, 
                            reduction = "umap", 
                            features = metrics,
                            pt.size = 0.8, 
                            order = TRUE,
                            min.cutoff = 'q10',
                            label = TRUE) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 20)  # For facet titles if multiple plots
  )
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/feature_plot.png", plot = feature_plot, width = 25, height = 30, limitsize = FALSE)


PCA_group_by_sample <- PCAPlot(seurat_integrated, group.by = "sample")
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/PCA_group_by_sample.png", plot = PCA_group_by_sample, width = 15, height = 15, limitsize = FALSE)


PCA_split_by_sample <- PCAPlot(seurat_integrated, split.by = "sample")
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/PCA_split_by_sample.png", plot = PCA_split_by_sample, width = 40, height = 10, limitsize = FALSE)


# Boxplot of nGene per cluster
nGene_per_cluster <- ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, fill=integrated_snn_res.0.8)) + NoLegend()
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/nGene_per_cluster.png", plot = nGene_per_cluster, width = 40, height = 15, limitsize = FALSE)


# Boxplot of mitoRatio per cluster
mitoRatio_per_cluster <- ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=mitoRatio, fill=integrated_snn_res.0.8)) + NoLegend()
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/mitoRatio_per_cluster.png", plot = mitoRatio_per_cluster, width = 40, height = 15, limitsize = FALSE)


# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "umap_1", "umap_2")


# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)


# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "umap_1", "umap_2"))  %>%
  group_by(ident) %>%
  dplyr::summarise(x=mean(umap_1), y=mean(umap_2))


# Plotting a UMAP plot for each of the PCs
UMAP_of_each_PC <- map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(umap_1, umap_2)) +
    geom_point(aes(color = !!sym(pc)), 
               alpha = 0.7) +
    scale_color_gradient(guide = "none", 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data = umap_label, 
              aes(label = ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/UMAP_of_each_PC.png", plot = UMAP_of_each_PC, width = 30, height = 40, limitsize = FALSE)


# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)


UMAP_no_legend <- DimPlot(object = seurat_integrated, reduction = "umap", label = TRUE) + NoLegend()
ggsave("results/Clustering/figures/integrated/QC/resolution_0.8/UMAP_no_legend.png", plot = UMAP_no_legend, width = 10, height = 10, limitsize = FALSE)
