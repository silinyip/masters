# Single-cell RNA-seq - manual annotation
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer
# Markers from Wu et al., 2023


# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(42)


# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(scales)
library(tidyverse)
library(RCurl)
library(cowplot)
library(future)


# Setup parallelization to speed up computational tasks
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)


# Load integrated Seurat object
seurat_integrated <- readRDS("results/Clustering/data/clustering/seurat_integrated.rds")


Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

Idents(seurat_integrated) <- factor(Idents(seurat_integrated), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                                                          "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26"))

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"


# Join layers in the RNA assay
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])


T_Cells <- VlnPlot(seurat_integrated,
                   features = c("CD3D", "CD3E", "CD3G", "CD2", "CD7", "TRAC", "D4", "IL7R", "CD8A", "CD8B", "GZMA", "GZMK", "PRF1", "FOXP3", "IL2RA", "CD25", "CTLA4"),
                   assay = "RNA",
                   layer = "data",
                   stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/T_Cells.pdf", plot = T_Cells)


Epithelial_Cells <- VlnPlot(seurat_integrated,
                            features = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN3", "CLDN4", "MUC1", "MUC2"),
                            assay = "RNA",
                            layer = "data",
                            stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Epithelial_Cells.pdf", plot = Epithelial_Cells)


B_Cells <- VlnPlot(seurat_integrated,
                   features = c("CD19", "CD79A", "CD79B", "MS4A1", "CD20", "MZB1", "JCHAIN", "SDC1", "CD138"),
                   assay = "RNA",
                   layer = "data",
                   stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/B_Cells.pdf", plot = B_Cells)


Dendritic_Cells <- VlnPlot(seurat_integrated,
                           features = c("CD1C", "BDCA1", "CLEC9A", "BATF3", "LILRA4", "IL3RA", "CD123", "HLA-DRA", "ITGAX", "CD11c"),
                           assay = "RNA",
                           layer = "data",
                           stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Dendritic_Cells.pdf", plot = Dendritic_Cells)


Hematopoietic_Stem_Cells <- VlnPlot(seurat_integrated,
                                    features = c("CD34", "PROM1", "CD133", "KIT", "CD117", "FLT3", "GATA2", "THY1", "CD90"),
                                    assay = "RNA",
                                    layer = "data",
                                    stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Hematopoietic_Stem_Cells.pdf", plot = Hematopoietic_Stem_Cells)


Monocytes <- VlnPlot(seurat_integrated,
                     features = c("CD14", "FCGR3A", "CD16", "LYZ", "S100A8", "S100A9"),
                     assay = "RNA",
                     layer = "data",
                     stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Monocytes.pdf", plot = Monocytes)


Fibroblasts <- VlnPlot(seurat_integrated,
                       features = c("DCN", "COL1A1", "COL3A1", "PDGFRA", "PDGFRB", "FAP", "ACTA2", "αSMA"),
                       assay = "RNA",
                       layer = "data",
                       stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Fibroblasts.pdf", plot = Fibroblasts)


Myocytes <- VlnPlot(seurat_integrated,
                    features = c("ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "LMOD1"),
                    assay = "RNA",
                    layer = "data",
                    stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Myocytes.pdf", plot = Myocytes)


Endothelial_Cells <- VlnPlot(seurat_integrated,
                             features = c("PECAM1", "CD31", "VWF", "ENG", "CD105", "CLDN5", "CDH5", "PROX1", "LYVE1", "PDPN"),
                             assay = "RNA",
                             layer = "data",
                             stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Endothelial_Cells.pdf", plot = Endothelial_Cells)


##### ANNOTATION #####
# Annotate the clusters and save the annotated seurat object
seurat_integrated@meta.data <- mutate(
  seurat_integrated@meta.data,
  celltype_major = case_when(
    integrated_snn_res.0.4 %in% c(
      5, 9, 13, 15, 22
    ) ~ "B Cells",
    integrated_snn_res.0.4 %in% c(
      1, 2
    ) ~ "T Cells",
    integrated_snn_res.0.4 %in% c(
      16, 24
    ) ~ "Dendritic Cells",
    integrated_snn_res.0.4 %in% c(
      19
    ) ~ "Endothelial Cells",
    integrated_snn_res.0.4 %in% c(
      3, 4, 6, 8, 12, 14, 17, 18, 21, 26
    ) ~ "Epithelial Cells",
    integrated_snn_res.0.4 %in% c(
      10
    ) ~ "Fibroblasts",
    integrated_snn_res.0.4 %in% c(
      20
    ) ~ "Hematopoietic Cells",
    integrated_snn_res.0.4 %in% c(
      0, 7, 11, 25
    ) ~ "Monocytes",
    integrated_snn_res.0.4 %in% c(
      23
    ) ~ "Myocytes"
  )
)


saveRDS(seurat_integrated, "results/Clustering/data/annotation/seurat_annotated.rds")


annotated_clusters <- DimPlot(seurat_integrated, 
                              group.by = "celltype_major",
                              reduction = "umap",
                              label = FALSE) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave("results/Clustering/figures/integrated/annotation/manual/annotated_clusters.pdf", plot = annotated_clusters)
