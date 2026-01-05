# Single-cell RNA-seq - gene set enrichment analysis 
# Enhanced version with directional analysis (up/down regulated genes)
# Script ran on high performance Linux machine - access provided by Prof. Pieter De Maayer


# Clean environment
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F)
set.seed(42)


# Load libraries
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(ReactomePA)


# Get list of all significant DE results files
de_files <- list.files("results/DGEA/data/Primary_vs_Organoid/", 
                       pattern = "_signif_genes.csv$", full.names = TRUE)

cat("Found", length(de_files), "DE results files\n")
print(basename(de_files))


# Function to run enrichment analysis for a gene list
run_enrichment_analysis <- function(gene_list, analysis_name, output_dir, direction_label = "") {
  
  cat("\n--- Running enrichment for:", analysis_name, direction_label, "---\n")
  cat("Number of genes:", length(gene_list), "\n")
  
  # Skip if too few genes
  if (length(gene_list) < 5) {
    cat("Too few genes (< 5) - skipping\n")
    return(NULL)
  }
  
  # Clean gene names
  gene_list <- trimws(gene_list)
  gene_list <- gene_list[gene_list != ""]
  gene_list <- unique(gene_list)
  
  cat("After cleaning:", length(gene_list), "unique genes\n")
  
  # Test gene ID conversion
  cat("Testing gene ID mapping...\n")
  test_mapping <- bitr(gene_list, 
                       fromType = "SYMBOL",
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
  
  use_alias <- FALSE
  if (nrow(test_mapping) == 0) {
    cat("Trying ALIAS mapping...\n")
    test_mapping_alias <- bitr(gene_list, 
                               fromType = "ALIAS",
                               toType = "ENTREZID", 
                               OrgDb = org.Hs.eg.db)
    
    if (nrow(test_mapping_alias) == 0) {
      cat("ERROR: No genes could be mapped. Skipping.\n")
      return(NULL)
    } else {
      cat("Success with ALIAS mapping:", nrow(test_mapping_alias), "genes mapped\n")
      use_alias <- TRUE
    }
  } else {
    cat("Success with SYMBOL mapping:", nrow(test_mapping), "genes mapped\n")
  }
  
  from_type <- if (use_alias) "ALIAS" else "SYMBOL"
  
  # Create subdirectory if direction label exists
  if (direction_label != "") {
    output_subdir <- paste0(output_dir, "/", direction_label)
    dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  } else {
    output_subdir <- output_dir
  }
  
  results_summary <- list()
  
  # 1. GO Biological Process
  cat("Running GO Biological Process...\n")
  go_bp <- enrichGO(gene = gene_list,
                    OrgDb = org.Hs.eg.db,
                    keyType = from_type,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
  
  if (!is.null(go_bp) && nrow(go_bp) > 0) {
    write.csv(as.data.frame(go_bp), 
              file = paste0(output_subdir, "/GO_BP_results.csv"), 
              row.names = FALSE)
    
    p1 <- dotplot(go_bp, showCategory = 15, label_format = 50) + 
      ggtitle(paste("GO Biological Process -", analysis_name, direction_label)) + 
      theme(axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            plot.margin = margin(50, 20, 20, 80))
    
    ggsave(paste0(output_subdir, "/GO_BP_dotplot.png"), p1, 
           width = 14, height = 10, dpi = 300, bg = "white")
    
    cat("GO BP results saved:", nrow(go_bp), "terms\n")
    results_summary$GO_BP <- nrow(go_bp)
  } else {
    cat("No significant GO BP terms found\n")
    results_summary$GO_BP <- 0
  }
  
  # 2. GO Molecular Function
  cat("Running GO Molecular Function...\n")
  go_mf <- enrichGO(gene = gene_list,
                    OrgDb = org.Hs.eg.db,
                    keyType = from_type,
                    ont = "MF",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
  
  if (!is.null(go_mf) && nrow(go_mf) > 0) {
    write.csv(as.data.frame(go_mf), 
              file = paste0(output_subdir, "/GO_MF_results.csv"), 
              row.names = FALSE)
    
    p2 <- dotplot(go_mf, showCategory = 15, label_format = 50) + 
      ggtitle(paste("GO Molecular Function -", analysis_name, direction_label)) +
      theme(axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            plot.margin = margin(50, 20, 20, 80))
    
    ggsave(paste0(output_subdir, "/GO_MF_dotplot.png"), p2, 
           width = 14, height = 10, dpi = 300, bg = "white")
    
    cat("GO MF results saved:", nrow(go_mf), "terms\n")
    results_summary$GO_MF <- nrow(go_mf)
  } else {
    cat("No significant GO MF terms found\n")
    results_summary$GO_MF <- 0
  }
  
  # 3. GO Cellular Component
  cat("Running GO Cellular Component...\n")
  go_cc <- enrichGO(gene = gene_list,
                    OrgDb = org.Hs.eg.db,
                    keyType = from_type,
                    ont = "CC",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
  
  if (!is.null(go_cc) && nrow(go_cc) > 0) {
    write.csv(as.data.frame(go_cc), 
              file = paste0(output_subdir, "/GO_CC_results.csv"), 
              row.names = FALSE)
    
    p3 <- dotplot(go_cc, showCategory = 15, label_format = 50) + 
      ggtitle(paste("GO Cellular Component -", analysis_name, direction_label)) +
      theme(axis.text.y = element_text(size = rel(1.5)),
            axis.text.x = element_text(size = rel(1.5)),
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            plot.margin = margin(50, 20, 20, 80))
    
    ggsave(paste0(output_subdir, "/GO_CC_dotplot.png"), p3, 
           width = 14, height = 10, dpi = 300, bg = "white")
    
    cat("GO CC results saved:", nrow(go_cc), "terms\n")
    results_summary$GO_CC <- nrow(go_cc)
  } else {
    cat("No significant GO CC terms found\n")
    results_summary$GO_CC <- 0
  }
  
  # 4. KEGG Pathway Analysis
  cat("Running KEGG pathway analysis...\n")
  gene_entrez <- bitr(gene_list, 
                      fromType = from_type,
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 0) {
    cat("Converted", nrow(gene_entrez), "genes to ENTREZ IDs for KEGG\n")
    
    kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID,
                               organism = 'hsa',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.2)
    
    if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
      write.csv(as.data.frame(kegg_results), 
                file = paste0(output_subdir, "/KEGG_results.csv"), 
                row.names = FALSE)
      
      p4 <- dotplot(kegg_results, showCategory = 15, label_format = 50) + 
        ggtitle(paste("KEGG Pathways -", analysis_name, direction_label)) +
        theme(axis.text.y = element_text(size = rel(1.5)),
              axis.text.x = element_text(size = rel(1.5)),
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              plot.margin = margin(50, 20, 20, 80))
      
      ggsave(paste0(output_subdir, "/KEGG_dotplot.png"), p4, 
             width = 14, height = 10, dpi = 300, bg = "white")
      
      cat("KEGG results saved:", nrow(kegg_results), "terms\n")
      results_summary$KEGG <- nrow(kegg_results)
    } else {
      cat("No significant KEGG terms found\n")
      results_summary$KEGG <- 0
    }
  } else {
    cat("No genes could be converted to ENTREZ IDs for KEGG\n")
    results_summary$KEGG <- 0
  }
  
  # 5. Reactome Pathway Analysis
  cat("Running Reactome pathway analysis...\n")
  if (nrow(gene_entrez) > 0) {
    tryCatch({
      reactome_results <- enrichPathway(gene = gene_entrez$ENTREZID,
                                        organism = "human",
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 0.2,
                                        readable = TRUE)
      
      if (!is.null(reactome_results) && nrow(reactome_results) > 0) {
        write.csv(as.data.frame(reactome_results), 
                  file = paste0(output_subdir, "/Reactome_results.csv"), 
                  row.names = FALSE)
        
        p5 <- dotplot(reactome_results, showCategory = 15, label_format = 50) + 
          ggtitle(paste("Reactome Pathways -", analysis_name, direction_label)) +
          theme(axis.text.y = element_text(size = rel(1.5), hjust = 1, margin = margin(r = 10)),
                axis.text.x = element_text(size = rel(1.5)),
                plot.title = element_text(size = rel(1.5), hjust = 0.5),
                plot.margin = margin(50, 20, 20, 80))
        
        ggsave(paste0(output_subdir, "/Reactome_dotplot.png"), p5, 
               width = 14, height = 10, dpi = 300, bg = "white")
        
        cat("Reactome results saved:", nrow(reactome_results), "terms\n")
        results_summary$Reactome <- nrow(reactome_results)
      } else {
        cat("No significant Reactome terms found\n")
        results_summary$Reactome <- 0
      }
    }, error = function(e) {
      cat("Reactome analysis failed:", e$message, "\n")
      results_summary$Reactome <- 0
    })
  }
  
  return(results_summary)
}


# Main processing loop
all_results_summary <- data.frame()

for (i in 1:length(de_files)) {
  
  file <- de_files[i]
  file_name <- basename(file)
  analysis_name <- gsub("_signif_genes.csv$", "", file_name)
  
  cat("\n\n========================================\n")
  cat("=== Processing:", analysis_name, "===\n")
  cat("========================================\n")
  
  # Read the significant genes file
  de_results <- read.csv(file, stringsAsFactors = FALSE)
  
  cat("Total significant genes:", nrow(de_results), "\n")
  
  # Check if avg_log2FC column exists
  if (!"avg_log2FC" %in% colnames(de_results)) {
    cat("WARNING: avg_log2FC column not found. Checking for alternatives...\n")
    # Try common alternative column names
    fc_col <- NULL
    for (col_name in c("log2FoldChange", "logFC", "log2FC", "avg_logFC")) {
      if (col_name %in% colnames(de_results)) {
        fc_col <- col_name
        cat("Using", fc_col, "as fold change column\n")
        break
      }
    }
    if (is.null(fc_col)) {
      cat("ERROR: Could not find fold change column. Skipping directional analysis.\n")
      next
    }
    de_results$avg_log2FC <- de_results[[fc_col]]
  }
  
  # Extract cell type from filename for clean directory names
  if (grepl("^B Cells_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "B_Cells"
  } else if (grepl("^Endothelial Cells_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Endothelial_Cells"
  } else if (grepl("^Epithelial Cells_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Epithelial_Cells"
  } else if (grepl("^Fibroblasts_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Fibroblasts"
  } else if (grepl("^Macrophages_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Macrophages"
  } else if (grepl("^Mast Cells_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Mast_Cells"
  } else if (grepl("^Plasma_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "Plasma_Cells"
  } else if (grepl("^T Cells_", analysis_name, ignore.case = TRUE)) {
    dir_name <- "T_Cells"
  } else {
    dir_name <- gsub("[^A-Za-z0-9_]", "_", analysis_name)
  }
  
  output_dir <- paste0("results/GSEA/", dir_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Split genes by direction
  de_results_up <- de_results[de_results$avg_log2FC > 0, ]  # Upregulated in Tumour
  de_results_down <- de_results[de_results$avg_log2FC < 0, ]  # Upregulated in Organoid
  
  cat("\nUpregulated in Tumour:", nrow(de_results_up), "genes\n")
  cat("Upregulated in Organoid:", nrow(de_results_down), "genes\n")
  
  # Extract gene lists
  gene_list_all <- de_results$gene
  gene_list_up <- de_results_up$gene
  gene_list_down <- de_results_down$gene
  
  # Run enrichment for ALL genes (combined analysis)
  cat("\n*** Running COMBINED analysis (all significant genes) ***\n")
  summary_all <- run_enrichment_analysis(gene_list_all, analysis_name, output_dir, 
                                         direction_label = "Combined")
  
  # Run enrichment for UPREGULATED genes (Tumour-enriched)
  cat("\n*** Running TUMOUR-ENRICHED analysis (upregulated in tumour) ***\n")
  summary_up <- run_enrichment_analysis(gene_list_up, analysis_name, output_dir, 
                                        direction_label = "Tumour_Enriched")
  
  # Run enrichment for DOWNREGULATED genes (Organoid-enriched)
  cat("\n*** Running ORGANOID-ENRICHED analysis (upregulated in organoid) ***\n")
  summary_down <- run_enrichment_analysis(gene_list_down, analysis_name, output_dir, 
                                          direction_label = "Organoid_Enriched")
  
  # Compile summary
  if (!is.null(summary_all) || !is.null(summary_up) || !is.null(summary_down)) {
    row_summary <- data.frame(
      CellType = dir_name,
      Analysis = analysis_name,
      Total_Genes = nrow(de_results),
      Tumour_Up_Genes = nrow(de_results_up),
      Organoid_Up_Genes = nrow(de_results_down),
      Combined_GO_BP = ifelse(!is.null(summary_all), summary_all$GO_BP, 0),
      Combined_KEGG = ifelse(!is.null(summary_all), summary_all$KEGG, 0),
      Tumour_GO_BP = ifelse(!is.null(summary_up), summary_up$GO_BP, 0),
      Tumour_KEGG = ifelse(!is.null(summary_up), summary_up$KEGG, 0),
      Organoid_GO_BP = ifelse(!is.null(summary_down), summary_down$GO_BP, 0),
      Organoid_KEGG = ifelse(!is.null(summary_down), summary_down$KEGG, 0)
    )
    all_results_summary <- rbind(all_results_summary, row_summary)
  }
  
  cat("\n=== Completed:", analysis_name, "===\n")
  cat("Results saved in:", output_dir, "\n")
}

cat("\n\n========================================\n")
cat("=== ALL ANALYSES COMPLETE ===\n")
cat("========================================\n")
cat("Results saved in: results/GSEA/\n\n")

# Save comprehensive summary
if (nrow(all_results_summary) > 0) {
  write.csv(all_results_summary, "results/GSEA/comprehensive_summary.csv", row.names = FALSE)
  cat("Comprehensive summary saved: results/GSEA/comprehensive_summary.csv\n")
  
  # Print summary table
  cat("\n=== ANALYSIS SUMMARY ===\n")
  print(all_results_summary)
}

cat("\nScript completed successfully!\n")