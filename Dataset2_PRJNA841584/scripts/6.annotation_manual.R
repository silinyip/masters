# Single-cell RNA-seq - manual annotation
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer
# Markers from Li et al., 2023


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


Epithelial_Cells <- VlnPlot(seurat_integrated,
                            features = c("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4", "MUC1", "CEACAM5", "CEACAM6", "CDH1", "VIL1", "MUC2"),
                            assay = "RNA",
                            layer = "data",
                            stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Epithelial_Cells.pdf", plot = Epithelial_Cells)


Fibroblasts <- VlnPlot(seurat_integrated,
                       features = c("VIM", "COL1A1", "COL1A2", "PDGFRA", "PDGFRB", "DCN", "COL3A1", "FAP", "THY1", "CD90"),
                       assay = "RNA",
                       layer = "data",
                       stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Fibroblasts.pdf", plot = Fibroblasts)


Myofibroblast <- VlnPlot(seurat_integrated,
                         features = c("VIM", "ACTA2","TAGLN", "MYH11", "POSTN", "FAP"),
                         assay = "RNA",
                         layer = "data",
                         stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Myofibroblast.pdf", plot = Myofibroblast)


Endothelial_Cells <- VlnPlot(seurat_integrated,
                             features = c("CLDN5", "CDH5", "VMF", "PECAM1", "CD31", "VWF", "TIE1", "KDR", "VEGFR2", "ENG", "ESAM"),
                             assay = "RNA",
                             layer = "data",
                             stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Endothelial_Cells.pdf", plot = Endothelial_Cells)


T_Cells <- VlnPlot(seurat_integrated,
                   features = c("CD3D", "CD3E", "CD3G", "CD2", "CD7", "TRAC", "D4", "IL7R", "CD8A", "CD8B", "GZMA", "GZMK", "PRF1", "FOXP3", "IL2RA", "CD25", "CTLA4"),
                   assay = "RNA",
                   layer = "data",
                   stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/T_Cells.pdf", plot = T_Cells)


B_Cells <- VlnPlot(seurat_integrated,
                   features = c("CD19", "CD79A", "CD79B", "MS4A1", "CD20", "MZB1", "JCHAIN", "SDC1", "CD138", "CD22", "BANK1"),
                   assay = "RNA",
                   layer = "data",
                   stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/B_Cells.pdf", plot = B_Cells)


Plasma_Cells <- VlnPlot(seurat_integrated,
                        features = c("JCHAIN", "IGHA", "IGHG", "IGHM", "MZB1", "SDC1", "CD138", "PRDM1", "XBP1", "IGHG1", "IGHA1", "CD79A"),
                        assay = "RNA",
                        layer = "data",
                        stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Plasma_Cells.pdf", plot = Plasma_Cells)


Macrophage <- VlnPlot(seurat_integrated,
                      features = c("CD68", "LYZ", "CD163", "MRC1", "CD206", "CSF1R", "C1QA", "MARCO", "MSR1"),
                      assay = "RNA",
                      layer = "data",
                      stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Macrophage.pdf", plot = Macrophage)


Dendritic_Cells <- VlnPlot(seurat_integrated,
                           features = c("CD1C", "BATF3", "ITGAX", "CD11c", "CLEC9A", "XCR1", "IRF8", "LILRA4", "GZMB", "IRF7", "TCF4", "CLEC4C"),
                           assay = "RNA",
                           layer = "data",
                           stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Dendritic_Cells.pdf", plot = Dendritic_Cells)


Mast_Cells <- VlnPlot(seurat_integrated,
                      features = c("TPSAB1", "TPSD1", "CPA3", "KIT", "CD117", "TPSB2", "HDC", "MS4A2", "FCER1A"),
                      assay = "RNA",
                      layer = "data",
                      stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Mast_Cells.pdf", plot = Mast_Cells)


Liver_Cells <- VlnPlot(seurat_integrated,
                       features = c("ALT", "AST", "ALP", "GGT", "AFP", "CK7", "CK19", "CD147", "PPIA", "TMSB10", "ALB", "APOB", "CYP3A4", "TTR", "TF", "HNF4A", "KRT7", "KRT19", "SOX9"),
                       assay = "RNA",
                       layer = "data",
                       stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Liver_Cells.pdf", plot = Liver_Cells)


##### ANNOTATION #####
# Annotate the clusters and save the annotated seurat object
seurat_integrated@meta.data <- mutate(
  seurat_integrated@meta.data,
  celltype_major = case_when(
    integrated_snn_res.0.4 %in% c(
      9
    ) ~ "B Cells",
    integrated_snn_res.0.4 %in% c(
      24
    ) ~ "Dendritic Cells",
    integrated_snn_res.0.4 %in% c(
      10
    ) ~ "Endothelial Cells",
    integrated_snn_res.0.4 %in% c(
      0, 1, 3, 4, 6, 11, 14, 15, 16, 18, 19, 20, 21, 22, 26
    ) ~ "Epithelial Cells",
    integrated_snn_res.0.4 %in% c(
      12, 13
    ) ~ "Fibroblasts",
    integrated_snn_res.0.4 %in% c(
      2, 8
    ) ~ "Macrophages",
    integrated_snn_res.0.4 %in% c(
      25
    ) ~ "Mast Cells",
    integrated_snn_res.0.4 %in% c(
      23
    ) ~ "Myofibroblasts",
    integrated_snn_res.0.4 %in% c(
      7, 17
    ) ~ "Plasma",
    integrated_snn_res.0.4 %in% c(
      5
    ) ~ "T Cells"
  )
)


saveRDS(seurat_integrated, file = "results/Clustering/data/annotation/seurat_annotated.rds")

seurat_integrated <- readRDS("results/Clustering/data/annotation/seurat_annotated.rds")

annotated_clusters <- DimPlot(seurat_integrated, 
                              group.by = "celltype_major",
                              reduction = "umap",
                              label = FALSE) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave("results/Clustering/figures/integrated/annotation/manual/annotated_clusters.pdf", plot = annotated_clusters)
