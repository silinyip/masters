# Single-cell RNA-seq analysis - QC
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
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(future)
library(ggplot2)
library(gridExtra)


# Setup parallelization to speed up computational tasks
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)


# Load the PRJNA932556 dataset and initialise the Seurat object with the raw (non-normalized data)
SAMN33196636_S1.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196636/outs/raw_feature_bc_matrix/")
SAMN33196636_S1 <- CreateSeuratObject(counts = SAMN33196636_S1.data, project = "SAMN33196636_S1", min.features = 100)

SAMN33196637_S2.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196637/outs/raw_feature_bc_matrix/")
SAMN33196637_S2 <- CreateSeuratObject(counts = SAMN33196637_S2.data, project = "SAMN33196637_S2", min.features = 100)

SAMN33196638_S3.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196638/outs/raw_feature_bc_matrix/")
SAMN33196638_S3 <- CreateSeuratObject(counts = SAMN33196638_S3.data, project = "SAMN33196638_S3", min.features = 100)

SAMN33196639_R1.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196639/outs/raw_feature_bc_matrix/")
SAMN33196639_R1 <- CreateSeuratObject(counts = SAMN33196639_R1.data, project = "SAMN33196639_R1", min.features = 100)

SAMN33196640_R2.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196640/outs/raw_feature_bc_matrix/")
SAMN33196640_R2 <- CreateSeuratObject(counts = SAMN33196640_R2.data, project = "SAMN33196640_R2", min.features = 100)

SAMN33196641_R3.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA932556/CellRanger/SAMN33196641/outs/raw_feature_bc_matrix/")
SAMN33196641_R3 <- CreateSeuratObject(counts = SAMN33196641_R3.data, project = "SAMN33196641_R3", min.features = 100)


# Merge the Seurat objects - storing all information in a single object for ease of use
merged_seurat <- merge(x = SAMN33196636_S1, 
                       y = c(SAMN33196637_S2, SAMN33196638_S3, SAMN33196639_R1, SAMN33196640_R2, SAMN33196641_R3), 
                       add.cell.ids = c("SAMN33196636_S1", "SAMN33196637_S2", "SAMN33196638_S3", "SAMN33196639_R1", "SAMN33196640_R2", "SAMN33196641_R3"))


# Concatenate the count matrices of samples together
merged_seurat <- JoinLayers(merged_seurat)


# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)


# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100


# Create metadata dataframe
metadata <- merged_seurat@meta.data


# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196636_S1_"))] <- "SAMN33196636_S1"
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196637_S2_"))] <- "SAMN33196637_S2"
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196638_S3_"))] <- "SAMN33196638_S3"
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196639_R1_"))] <- "SAMN33196639_R1"
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196640_R2_"))] <- "SAMN33196640_R2"
metadata$sample[which(str_detect(metadata$cells, "^SAMN33196641_R3_"))] <- "SAMN33196641_R3"



# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file = "results/QC/data/merged_filtered_seurat.RData")


# Number of cells per sample
table(merged_seurat$sample)


# Number of cells in the merged Seurat object
cat("Number of cells in merged_seurat:", ncol(merged_seurat), "\n")


# Visualize the number of cell counts per sample
unfiltered_cell_counts_per_sample <- metadata %>%
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("results/QC/figures/unfiltered/unfiltered_cell_counts_per_sample.png", plot = unfiltered_cell_counts_per_sample, width = 15, height = 15, limitsize = FALSE)


# Visualize the number UMIs/transcripts per cell
UMI_counts_per_cell <- metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("results/QC/figures/unfiltered/UMI_counts_per_cell.png", plot = UMI_counts_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the distribution of genes detected per cell via histogram
gene_counts_per_cell <- metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("results/QC/figures/unfiltered/gene_counts_per_cell.png", plot = gene_counts_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
gene_expression_complexity <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("results/QC/figures/unfiltered/gene_expression_complexity.png", plot = gene_expression_complexity, width = 15, height = 15, limitsize = FALSE)


# Visualize the distribution of mitochondrial gene expression detected per cell
mito_gene_expression_per_cell <- metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave("results/QC/figures/unfiltered/mito_gene_expression_per_cell.png", plot = mito_gene_expression_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
genes_UMI_correlation <- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
ggsave("results/QC/figures/unfiltered/genes_UMI_correlation.png", plot = genes_UMI_correlation, width = 30, height = 30, dpi = 600, limitsize = FALSE)


# Visualise QC metrics as a violin plot
consistent_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18, margin = margin(t = 15)), # Increase top margin to 15
  axis.text.y = element_text(size = 25),
  axis.title.x = element_text(size = 20, margin = margin(t = 20), vjust = -1), # Increase top margin for x-axis title to 20 and adjust vertical position
  axis.title.y = element_text(size = 20),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
  strip.text = element_text(size = 18),
  plot.margin = unit(c(1, 1, 2.5, 7), "cm"), # Increase left margin
  legend.position = "none" # Remove legends
)

# Create individual plots with consistent theme and updated sample names
plot_nGene <- VlnPlot(merged_seurat, features = "nGene", pt.size = 0.0001, group.by = "sample", raster = FALSE) + consistent_theme
plot_nUMI <- VlnPlot(merged_seurat, features = "nUMI", pt.size = 0.0001, group.by = "sample", raster = FALSE) + consistent_theme
plot_mitoRatio <- VlnPlot(merged_seurat, features = "mitoRatio", pt.size = 0.0001, group.by = "sample", raster = FALSE) + consistent_theme

# Combine plots
violin_plot <- gridExtra::grid.arrange(plot_nGene, plot_nUMI, plot_mitoRatio, ncol = 3)
ggsave("results/QC/figures/unfiltered/violin_plot.png", plot = violin_plot, width = 40, height = 20, limitsize = FALSE)


########## Cell-level filtering ########## 
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset = (nGene >= 200)) 


########## Gene-level filtering ########## 
# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data


# Number of cells per sample
table(filtered_seurat$sample)


# Number of cells in the merged Seurat object
cat("Number of cells in filtered_seurat:", ncol(filtered_seurat), "\n")


# Visualize the number of cell counts per sample
filtered_cell_counts_per_sample <- metadata_clean %>%
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("results/QC/figures/filtered/filtered_cell_counts_per_sample.png", plot = filtered_cell_counts_per_sample, width = 15, height = 15, limitsize = FALSE)


# Visualize the number UMIs/transcripts per cell
UMI_counts_per_cell <- metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("results/QC/figures/filtered/UMI_counts_per_cell.png", plot = UMI_counts_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the distribution of genes detected per cell via histogram
gene_counts_per_cell <- metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave("results/QC/figures/filtered/gene_counts_per_cell.png", plot = gene_counts_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
gene_expression_complexity <- metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("results/QC/figures/filtered/gene_expression_complexity.png", plot = gene_expression_complexity, width = 15, height = 15, limitsize = FALSE)


# Visualize the distribution of mitochondrial gene expression detected per cell
mito_gene_expression_per_cell <- metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave("results/QC/figures/filtered/mito_gene_expression_per_cell.png", plot = mito_gene_expression_per_cell, width = 15, height = 15, limitsize = FALSE)


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
genes_UMI_correlation <- metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
ggsave("results/QC/figures/filtered/genes_UMI_correlation.png", plot = genes_UMI_correlation, width = 30, height = 30, dpi = 600, limitsize = FALSE)


# Visualise QC metrics as a violin plot
violin_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18, 
                             margin = margin(t = 15)),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_text(size = 22, margin = margin(t = 25)),
  axis.title.y = element_text(size = 22, margin = margin(r = 20)),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 24, 
                            margin = margin(b = 15)),
  strip.text = element_text(size = 20),
  plot.margin = unit(c(1.5, 1.5, 2, 2), "cm"),
  legend.position = "none"
)

# Create individual violin plots with enhanced theme
plot_nGene <- VlnPlot(filtered_seurat, features = "nGene", pt.size = 0.0001, 
                      group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("Genes per Cell") +
  ylab("Number of Genes")

plot_nUMI <- VlnPlot(filtered_seurat, features = "nUMI", pt.size = 0.0001, 
                     group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("UMIs per Cell") +
  ylab("Number of UMIs")

plot_mitoRatio <- VlnPlot(filtered_seurat, features = "mitoRatio", pt.size = 0.0001, 
                          group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("Mitochondrial Ratio") +
  ylab("Mitochondrial Ratio")

# Combine plots with extra space
violin_plot <- gridExtra::grid.arrange(plot_nGene, plot_nUMI, plot_mitoRatio, ncol = 3)
ggsave("results/QC/figures/filtered/violin_plot.png", 
       plot = violin_plot, width = 45, height = 16, dpi = 300, limitsize = FALSE)


# Create .RData object to load at any time
save(filtered_seurat, file = "results/QC/data/seurat_filtered.RData")
