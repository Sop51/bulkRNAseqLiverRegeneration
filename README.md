# Bulk RNA-seq Analysis: Zebrafish Liver Regeneration

This repository contains scripts for the analysis of bulk RNA-seq data generated from a zebrafish liver regeneration model. Liver samples were collected at three timepoints following ablation:

- **0 days post ablation (dpa)**
- **1.5 days post ablation**
- **4.5 days post ablation**

---

## Analysis Overview

The workflow consists of the following major steps:

1. Differential expression analysis using **edgeR**
2. Quality control and PCA visualization
3. Functional enrichment analysis (GO, KEGG, Reactome)
4. Visualization of gene expression dynamics across timepoints
5. Targeted analyses of liver markers, contamination, and key biological pathways

---

## File Structure

### Core Analysis

- **`edgeR.R`**  
  Performs differential expression analysis using edgeR and generates initial PCA plots across all samples.

- **`control_sample_plots.R`**  
  Conducts KEGG and Gene Ontology (GO) enrichment analyses and generates a PCA plot using control samples only.

- **`plotting.R`**  
  Generates downstream visualizations, including:
  - KEGG and GO enrichment plots  
  - Gene expression time-course plots  
  - Volcano plots  
  - Bar graph summarizing the number of differentially expressed genes across timepoints

- **`enrichment_analysis.R`**  
  Creates enrichment network plots for treated samples using **Reactome pathway** analysis.

---

### Targeted Gene Set Analyses

- **`markers_and_contamination_plots.R`**  
  Produces plots assessing nearby tissue contamination and expression of liver-specific marker genes.

- **`tgfBetaGenes.R`**  
  Generates differential expression plots for genes involved in the **TGF-β signaling pathway**.

- **`ecmGenes.R`**  
  Generates differential expression plots for genes related to the **extracellular matrix (ECM)**.

- **`inflammationGenes.R`**  
  Generates differential expression plots for genes associated with the **inflammatory response** (GO Biological Process).
