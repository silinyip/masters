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
SAMN28600722_LiverMet_44.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600722/outs/raw_feature_bc_matrix/")
SAMN28600722_LiverMet_44 <- CreateSeuratObject(counts = SAMN28600722_LiverMet_44.data, project = "SAMN28600722_LiverMet_44", min.features = 100)

SAMN28600723_Liver_44.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600723/outs/raw_feature_bc_matrix/")
SAMN28600723_Liver_44 <- CreateSeuratObject(counts = SAMN28600723_Liver_44.data, project = "SAMN28600723_Liver_44", min.features = 100)

SAMN28600724_Tumour_44.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600724/outs/raw_feature_bc_matrix/")
SAMN28600724_Tumour_44 <- CreateSeuratObject(counts = SAMN28600724_Tumour_44.data, project = "SAMN28600724_Tumour_44", min.features = 100)

SAMN28600725_Tumoroid_28.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600725/outs/raw_feature_bc_matrix/")
SAMN28600725_Tumoroid_28 <- CreateSeuratObject(counts = SAMN28600725_Tumoroid_28.data, project = "SAMN28600725_Tumoroid_28", min.features = 100)

SAMN28600726_Colonoid_28.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600726/outs/raw_feature_bc_matrix/")
SAMN28600726_Colonoid_28 <- CreateSeuratObject(counts = SAMN28600726_Colonoid_28.data, project = "SAMN28600726_Colonoid_28", min.features = 100)

SAMN28600727_Tumour_28.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600727/outs/raw_feature_bc_matrix/")
SAMN28600727_Tumour_28 <- CreateSeuratObject(counts = SAMN28600727_Tumour_28.data, project = "SAMN28600727_Tumour_28", min.features = 100)

SAMN28600728_Colonoid_26.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600728/outs/raw_feature_bc_matrix/")
SAMN28600728_Colonoid_26 <- CreateSeuratObject(counts = SAMN28600728_Colonoid_26.data, project = "SAMN28600728_Colonoid_26", min.features = 100)

SAMN28600729_Tumoroid_24.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600729/outs/raw_feature_bc_matrix/")
SAMN28600729_Tumoroid_24 <- CreateSeuratObject(counts = SAMN28600729_Tumoroid_24.data, project = "SAMN28600729_Tumoroid_24", min.features = 100)

SAMN28600730_Colonoid_24.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600730/outs/raw_feature_bc_matrix/")
SAMN28600730_Colonoid_24 <- CreateSeuratObject(counts = SAMN28600730_Colonoid_24.data, project = "SAMN28600730_Colonoid_24", min.features = 100)

SAMN28600731_LiverMet_27.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600731/outs/raw_feature_bc_matrix/")
SAMN28600731_LiverMet_27 <- CreateSeuratObject(counts = SAMN28600731_LiverMet_27.data, project = "SAMN28600731_LiverMet_27", min.features = 100)

SAMN28600732_Tumour_27.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600732/outs/raw_feature_bc_matrix/")
SAMN28600732_Tumour_27 <- CreateSeuratObject(counts = SAMN28600732_Tumour_27.data, project = "SAMN28600732_Tumour_27", min.features = 100)

SAMN28600733_Tumour_26.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600733/outs/raw_feature_bc_matrix/")
SAMN28600733_Tumour_26 <- CreateSeuratObject(counts = SAMN28600733_Tumour_26.data, project = "SAMN28600733_Tumour_26", min.features = 100)

SAMN28600734_Tumour_24.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600734/outs/raw_feature_bc_matrix/")
SAMN28600734_Tumour_24 <- CreateSeuratObject(counts = SAMN28600734_Tumour_24.data, project = "SAMN28600734_Tumour_24", min.features = 100)

SAMN28600735_Tumoroid_08.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600735/outs/raw_feature_bc_matrix/")
SAMN28600735_Tumoroid_08 <- CreateSeuratObject(counts = SAMN28600735_Tumoroid_08.data, project = "SAMN28600735_Tumoroid_08", min.features = 100)

SAMN28600736_Colonoid_08.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600736/outs/raw_feature_bc_matrix/")
SAMN28600736_Colonoid_08 <- CreateSeuratObject(counts = SAMN28600736_Colonoid_08.data, project = "SAMN28600736_Colonoid_08", min.features = 100)

SAMN28600737_LiverMet_09.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600737/outs/raw_feature_bc_matrix/")
SAMN28600737_LiverMet_09 <- CreateSeuratObject(counts = SAMN28600737_LiverMet_09.data, project = "SAMN28600737_LiverMet_09", min.features = 100)

SAMN28600738_Liver_09.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600738/outs/raw_feature_bc_matrix/")
SAMN28600738_Liver_09 <- CreateSeuratObject(counts = SAMN28600738_Liver_09.data, project = "SAMN28600738_Liver_09", min.features = 100)

SAMN28600739_Tumour_08.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600739/outs/raw_feature_bc_matrix/")
SAMN28600739_Tumour_08 <- CreateSeuratObject(counts = SAMN28600739_Tumour_08.data, project = "SAMN28600739_Tumour_08", min.features = 100)

SAMN28600740_Tumour_09.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600740/outs/raw_feature_bc_matrix/")
SAMN28600740_Tumour_09 <- CreateSeuratObject(counts = SAMN28600740_Tumour_09.data, project = "SAMN28600740_Tumour_09", min.features = 100)

SAMN28600741_Tumour_07.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600741/outs/raw_feature_bc_matrix/")
SAMN28600741_Tumour_07 <- CreateSeuratObject(counts = SAMN28600741_Tumour_07.data, project = "SAMN28600741_Tumour_07", min.features = 100)

SAMN28600742_Colon_08.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600742/outs/raw_feature_bc_matrix/")
SAMN28600742_Colon_08 <- CreateSeuratObject(counts = SAMN28600742_Colon_08.data, project = "SAMN28600742_Colon_08", min.features = 100)

SAMN28600743_Colon_09.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600743/outs/raw_feature_bc_matrix/")
SAMN28600743_Colon_09 <- CreateSeuratObject(counts = SAMN28600743_Colon_09.data, project = "SAMN28600744_Colon_09", min.features = 100)

SAMN28600744_Colon_07.data <- Read10X(data.dir = "/media/silin/53e5b279-0573-4370-a1b9-a9f22dbfa3be/yard/PRJNA841584/CellRanger/SAMN28600744/outs/raw_feature_bc_matrix/")
SAMN28600744_Colon_07 <- CreateSeuratObject(counts = SAMN28600744_Colon_07.data, project = "SAMN28600744_Colon_07", min.features = 100)


# Merge the Seurat objects - storing all information in a single object for ease of use
merged_seurat <- merge(x = SAMN28600722_LiverMet_44, 
                       y = c(SAMN28600723_Liver_44, SAMN28600724_Tumour_44, SAMN28600725_Tumoroid_28, SAMN28600726_Colonoid_28, SAMN28600727_Tumour_28, SAMN28600728_Colonoid_26, SAMN28600729_Tumoroid_24, SAMN28600730_Colonoid_24, SAMN28600731_LiverMet_27, SAMN28600732_Tumour_27, SAMN28600733_Tumour_26, SAMN28600734_Tumour_24, SAMN28600735_Tumoroid_08, SAMN28600736_Colonoid_08, SAMN28600737_LiverMet_09, SAMN28600738_Liver_09, SAMN28600739_Tumour_08, SAMN28600740_Tumour_09, SAMN28600741_Tumour_07, SAMN28600742_Colon_08, SAMN28600743_Colon_09, SAMN28600744_Colon_07), 
                       add.cell.ids = c("SAMN28600722_LiverMet_44", "SAMN28600723_Liver_44", "SAMN28600724_Tumour_44", "SAMN28600725_Tumoroid_28", "SAMN28600726_Colonoid_28", "SAMN28600727_Tumour_28", "SAMN28600728_Colonoid_26", "SAMN28600729_Tumoroid_24", "SAMN28600730_Colonoid_24", "SAMN28600731_LiverMet_27", "SAMN28600732_Tumour_27", "SAMN28600733_Tumour_26", "SAMN28600734_Tumour_24", "SAMN28600735_Tumoroid_08", "SAMN28600736_Colonoid_08", "SAMN28600737_LiverMet_09", "SAMN28600738_Liver_09", "SAMN28600739_Tumour_08", "SAMN28600740_Tumour_09", "SAMN28600741_Tumour_07", "SAMN28600742_Colon_08", "SAMN28600743_Colon_09", "SAMN28600744_Colon_07"))


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
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600744_Colon_07_"))] <- "SAMN28600744_Colon_07"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600743_Colon_09_"))] <- "SAMN28600743_Colon_09"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600742_Colon_08_"))] <- "SAMN28600742_Colon_08"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600741_Tumour_07_"))] <- "SAMN28600741_Tumour_07"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600740_Tumour_09_"))] <- "SAMN28600740_Tumour_09"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600739_Tumour_08_"))] <- "SAMN28600739_Tumour_08"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600738_Liver_09_"))] <- "SAMN28600738_Liver_09"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600737_LiverMet_09_"))] <- "SAMN28600737_LiverMet_09"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600736_Colonoid_08_"))] <- "SAMN28600736_Colonoid_08"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600735_Tumoroid_08_"))] <- "SAMN28600735_Tumoroid_08"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600734_Tumour_24_"))] <- "SAMN28600734_Tumour_24"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600733_Tumour_26_"))] <- "SAMN28600733_Tumour_26"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600732_Tumour_27_"))] <- "SAMN28600732_Tumour_27"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600731_LiverMet_27_"))] <- "SAMN28600731_LiverMet_27"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600730_Colonoid_24_"))] <- "SAMN28600730_Colonoid_24"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600729_Tumoroid_24_"))] <- "SAMN28600729_Tumoroid_24"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600728_Colonoid_26_"))] <- "SAMN28600728_Colonoid_26"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600727_Tumour_28_"))] <- "SAMN28600727_Tumour_28"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600726_Colonoid_28_"))] <- "SAMN28600726_Colonoid_28"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600725_Tumoroid_28_"))] <- "SAMN28600725_Tumoroid_28"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600724_Tumour_44_"))] <- "SAMN28600724_Tumour_44"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600723_Liver_44_"))] <- "SAMN28600723_Liver_44"
metadata$sample[which(str_detect(metadata$cells, "^SAMN28600722_LiverMet_44_"))] <- "SAMN28600722_LiverMet_44"


# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata


# Create .RData object to load at any time
#save(merged_seurat, file = "results/QC/data/merged_filtered_seurat.RData")


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
                          subset = (nUMI >= 1000) & 
                            (nGene >= 500) & 
                            (mitoRatio < 0.20))


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
                             margin = margin(t = 15, b = 10)),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_text(size = 22, margin = margin(t = 25)),
  axis.title.y = element_text(size = 22, margin = margin(r = 20)),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 24, 
                            margin = margin(b = 15)),
  strip.text = element_text(size = 20),
  plot.margin = unit(c(6, 3, 3.5, 3), "cm"),  # Increased bottom and right margins
  legend.position = "none"
)

# Create individual violin plots with enhanced theme
plot_nGene <- VlnPlot(filtered_seurat, features = "nGene", pt.size = 0.0001, 
                      group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("Genes per Cell") +
  ylab("Number of Genes") +
  xlab("Sample")

plot_nUMI <- VlnPlot(filtered_seurat, features = "nUMI", pt.size = 0.0001, 
                     group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("UMIs per Cell") +
  ylab("Number of UMIs") +
  xlab("Sample")

plot_mitoRatio <- VlnPlot(filtered_seurat, features = "mitoRatio", pt.size = 0.0001, 
                          group.by = "sample", raster = FALSE) + 
  violin_theme +
  ggtitle("Mitochondrial Ratio") +
  ylab("Mitochondrial Ratio") +
  xlab("Sample")

# Combine plots with cowplot for better alignment
library(cowplot)
violin_plot <- plot_grid(
  plot_nGene, 
  plot_nUMI, 
  plot_mitoRatio, 
  ncol = 3,
  align = "h"
)

ggsave("results/QC/figures/filtered/violin_plot.png", 
       plot = violin_plot, width = 45, height = 18, dpi = 300, limitsize = FALSE)


# Create .RData object to load at any time
#save(filtered_seurat, file = "results/QC/data/seurat_filtered.RData")
