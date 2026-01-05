# masters

# Single-Cell RNA Sequencing Analysis of Colorectal Cancer Organoids and Primary Tissues

## 📋 Project Overview

This repository contains the computational analysis pipeline and results for my Master's dissertation investigating the cellular and transcriptional differences between patient-derived colorectal cancer (CRC) organoids and matched primary tumour tissues using single-cell RNA sequencing (scRNA-seq).

### Research Objectives

1. **Establish and validate** a robust scRNA-seq analysis workflow using three independent publicly available CRC datasets
2. **Quantify compositional differences** between organoid cultures and primary tumour tissues
3. **Characterise transcriptional divergence** at the epithelial cell level between in vitro and in vivo systems
4. **Identify functional pathways** lost or altered during organoid culture adaptation

### Key Findings

- **Compositional Segregation**: Organoids are 90.56% epithelial cells vs. 31.21% in primary tissues (χ² = 31,390, p < 2.2×10⁻¹⁶)
- **Transcriptional Simplification**: 88.7% of differentially expressed genes (2,730/3,078) are upregulated in primary tissues
- **Loss of TME Signaling**: Organoids systematically suppress immune-communication, cytokine-response, and ECM-interaction pathways
- **Preserved Differentiation**: Despite simplification, organoids maintain epithelial differentiation hierarchies (stem → progenitor → differentiated)

---

## 📊 Datasets Analysed

### Validation Datasets
1. **Dataset 1 (PRJNA725335)** - Che et al. (2021)
   - 136,446 cells from 6 CRC patients
   - Primary tumours, liver metastases, PBMC
   
2. **Dataset 2 (PRJNA841584)** - Li et al. (2023) ⭐ *Main Analysis Dataset*
   - 80,380 cells from 16 treatment-naïve CRC patients
   - Organoids, primary tumours, adjacent normal colon, metastases
   
3. **Dataset 3 (PRJNA932556)** - Wu et al. (2023)
   - 66,812 cells from 6 MMR-deficient/MSI-high CRC patients
   - Response to anti-PD-1 immunotherapy
