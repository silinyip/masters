# Single-cell RNA-seq - normalization
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
library(future)


# Setup parallelization to speed up computational tasks
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)


# Load filtered_seurat_merged object
load("results/QC/data/seurat_filtered.RData")


# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)


# Load cell cycle markers
load("data/cycle.rda")


# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)


# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)                                


# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)


# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]


# Plot the average expression and variance of these genes with labels to indicate which genes are in the top 15
top_15_HVG <- VariableFeaturePlot(seurat_phase) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),        # Axis tick labels
    axis.title = element_text(size = 25),       # Axis titles
    legend.text = element_text(size = 20),      # Legend labels
    legend.title = element_text(size = 25),     # Legend title
  )
top_15_HVG <- LabelPoints(plot = top_15_HVG, points = top_genes, repel = TRUE, size = 6)
ggsave("results/Normalisation/SCT/figures/PCA_with_cells/No_regressed_variables/top_15_HVG.png", plot = top_15_HVG, width = 15, height = 15, limitsize = FALSE)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)


# Plot the PCA colored by cell cycle phase
PCA_cell_cycle_phase <- DimPlot(seurat_phase,
                                reduction = "pca",
                                group.by= "Phase",
                                split.by = "Phase")
ggsave("results/Normalisation/SCT/figures/PCA_with_cells/No_regressed_variables/PCA_cell_cycle_phase.png", plot = PCA_cell_cycle_phase, width = 15, height = 15, limitsize = FALSE)


# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))


# Plot the PCA colored by mitoFr
PCA_mitoRatio <- DimPlot(seurat_phase,
                         reduction = "pca",
                         group.by= "mitoFr",
                         split.by = "mitoFr")
ggsave("results/Normalisation/SCT/figures/PCA_with_cells/No_regressed_variables/PCA_mitoRatio.png", plot = PCA_mitoRatio, width = 15, height = 15, limitsize = FALSE)


# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"), vst.flavor = "v2")
}


# Save the split seurat object
saveRDS(split_seurat, "results/Normalisation/SCT/data/split_seurat.rds")


# Run UMAP
UMAP_unintegrated_for_comparison <- RunUMAP(seurat_phase, dims = 1:40,reduction = "pca")


# Plot UMAP
UMAP_unintegrated_plot <- DimPlot(UMAP_unintegrated_for_comparison)
ggsave("results/Normalisation/SCT/figures/unintegrated/UMAP_unintegrated_no_regressed.png", plot = UMAP_unintegrated_plot, width = 10, height = 10, limitsize = FALSE)


split_seurat <- readRDS("results/Normalisation/SCT/data/split_seurat.rds")
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)


# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        reduction = "rpca",
                                        k.anchor = 95)


# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")


# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA with larger text
PCA_plot_integrated <- PCAPlot(seurat_integrated, 
                               split.by = "sample") +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.key.size = unit(1.2, "cm"),
    strip.text = element_text(size = 24),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("results/Normalisation/SCT/figures/integrated/PCA_plot_integrated.png", 
       plot = PCA_plot_integrated, width = 40, height = 10, dpi = 300, limitsize = FALSE)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP with larger text
UMAP_integrated <- DimPlot(seurat_integrated) +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 32, hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("results/Normalisation/SCT/figures/integrated/UMAP_integrated.png", 
       plot = UMAP_integrated, width = 20, height = 20, dpi = 300, limitsize = FALSE)

# Plot UMAP split by sample with larger text
UMAP_integrated_split_by_sample <- DimPlot(seurat_integrated, split.by = "sample") +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.key.size = unit(1.2, "cm"),
    strip.text = element_text(size = 24),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave("results/Normalisation/SCT/figures/integrated/UMAP_integrated_split_by_sample.png", 
       plot = UMAP_integrated_split_by_sample, width = 55, height = 10, dpi = 300, limitsize = FALSE)


# Save integrated seurat object
saveRDS(seurat_integrated, "results/Normalisation/SCT/data/integrated_seurat.rds")