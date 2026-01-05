# Single-cell RNA-seq - marker annotation
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer


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


Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"


# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined.markers <- FindAllMarkers(object = seurat_integrated, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25) 

save(combined.markers, file= "results/Clustering/data/marker_identification/resolution_0.8/combined_markers_res_08.RData")

annotations <- read.csv("data/annotation.csv")

top10_combined <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


# Extract the gene descriptions for each gene
gene_descriptions <- unique(annotations[, c("gene_name", "description")])


# Keep only the first occurrence of each gene name
gene_descriptions_unique <- gene_descriptions %>%
  distinct(gene_name, .keep_all = TRUE)


# Merge gene descriptions with markers
annotate_combined_markers <- left_join(x = combined.markers,
                                       y = gene_descriptions_unique,
                                       by = c("gene" = "gene_name"))


save(top10_combined, file= "results/Clustering/data/marker_identification/resolution_0.8/top10_markers_per_cluster_res_08.RData")

write.csv(top10_combined, file= "results/Clustering/data/marker_identification/resolution_0.8/top10_markers_per_cluster_res_08.csv")

write.csv(annotate_combined_markers, "results/Clustering/data/marker_identification/resolution_0.8/markers_SCSA_tool_res_08.csv")
