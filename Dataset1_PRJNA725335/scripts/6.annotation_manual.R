# Single-cell RNA-seq - manual annotation
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer
# Markers from Che et al., 2021


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

Idents(seurat_integrated) <- factor(Idents(seurat_integrated), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                                                                          "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"))

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


B_Cells <- VlnPlot(seurat_integrated,
                   features = c("CD19", "CD79A", "CD79B", "MS4A1", "CD20", "CD22", "BANK1"),
                   assay = "RNA",
                   layer = "data",
                   stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/B_Cells.pdf", plot = B_Cells)


Plasma_Cells <- VlnPlot(seurat_integrated,
                        features = c("MZB1", "SDC1", "CD138", "PRDM1", "XBP1", "JCHAIN", "IGHG1", "IGHA1", "CD79A"),
                        assay = "RNA",
                        layer = "data",
                        stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Plasma_Cells.pdf", plot = Plasma_Cells)


Myeloid_Cells <- VlnPlot(seurat_integrated,
                         features = c("LYZ", "CD14", "CD68", "ITGAM", "CD11b", "FCGR3A", "CD16", "CD14", "CD68", "CD163", "MRC1", "CD206", "MARCO", "ITGAX", "CD11c", "CLEC9A", "FLT3", "IRF8", "LAMP3"),
                         assay = "RNA",
                         layer = "data",
                         stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Myeloid_Cells.pdf", plot = Myeloid_Cells)


NK_Cells <- VlnPlot(seurat_integrated,
                    features = c("NCAM1", "CD56", "KLRD1", "CD94", "NKG7", "GNLY", "GZMB", "PRF1", "FCGR3A", "KLRF1", "FGFBP2"),
                    layer = "data",
                    stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/NK_Cells.pdf", plot = NK_Cells)


CAFs <- VlnPlot(seurat_integrated,
                features = c("ACTA2", "PDGFRA", "PDGFRB", "FAP", "TAGLN", "COL1A1", "COL3A1", "MMP2", "POSTN", "DCN"),
                assay = "RNA",
                layer = "data",
                stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/CAFs.pdf", plot = CAFs)


Endothelial_Cells <- VlnPlot(seurat_integrated,
                             features = c("CLDN5", "CDH5", "VMF", "PECAM1", "CD31", "VWF", "TIE1", "KDR", "VEGFR2", "ENG", "ESAM"),
                             assay = "RNA",
                             layer = "data",
                             stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Endothelial_Cells.pdf", plot = Endothelial_Cells)


Mast_Cells <- VlnPlot(seurat_integrated,
                      features = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC", "FCER1A"),
                      assay = "RNA",
                      layer = "data",
                      stack = TRUE, flip = TRUE, fill.by = "ident") +
  theme(axis.text.x=element_text(size=7, angle = 0), legend.position = "none") +
  xlab("cluster")

ggsave("results/Clustering/figures/integrated/annotation/manual/Mast_Cells.pdf", plot = Mast_Cells)


##### ANNOTATION #####
# Annotate the clusters and save the annotated seurat object
seurat_integrated@meta.data <- mutate(
  seurat_integrated@meta.data,
  celltype_major = case_when(
    integrated_snn_res.0.4 %in% c(
      8
    ) ~ "B Cells",
    integrated_snn_res.0.4 %in% c(
      12
    ) ~ "CAFs",
    integrated_snn_res.0.4 %in% c(
      21
    ) ~ "Endothelial Cells",
    integrated_snn_res.0.4 %in% c(
      15
    ) ~ "Mast Cells",
    integrated_snn_res.0.4 %in% c(
      2, 10
    ) ~ "Myeloid Cells",
    integrated_snn_res.0.4 %in% c(
      3, 14
    ) ~ "NK Cells",
    integrated_snn_res.0.4 %in% c(
      18, 20
    ) ~ "Plasma Cells",
    integrated_snn_res.0.4 %in% c(
      0, 1, 4, 5, 6, 7, 9, 11, 13, 16, 17, 19
    ) ~ "T Cells"
  )
)


saveRDS(seurat_integrated, file = "results/Clustering/data/annotation/seurat_annotated.rds")


annotated_clusters <- DimPlot(seurat_integrated, 
                              group.by = "celltype_major",
                              reduction = "umap",
                              label = FALSE) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave("results/Clustering/figures/integrated/annotation/manual/annotated_clusters.pdf", plot = annotated_clusters)
