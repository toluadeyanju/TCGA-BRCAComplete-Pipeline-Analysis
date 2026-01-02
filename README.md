# TCGA-BRCA Complete Pipeline Analysis: Basal vs. Luminal A

This repository contains a complete bioinformatics pipeline for the differential expression and clinical survival analysis of Breast Invasive Carcinoma (BRCA) using data from The Cancer Genome Atlas (TCGA).

## Project Overview
The goal of this study is to identify molecular drivers that distinguish the **Basal-like** subtype from the **Luminal A** subtype. These subtypes represent two of the most distinct clinical entities in breast cancer, with Basal-like tumors often being more aggressive and lacking traditional hormone receptor targets.


## Key Features
- **Data Acquisition:** Automated retrieval of TCGA-BRCA RNA-Seq data using `TCGAbiolinks`.
- **Differential Expression:** Robust statistical analysis using `DESeq2` to identify significantly up- and down-regulated genes.
- **Functional Enrichment:** Pathway analysis via Gene Ontology (GO) and KEGG to understand biological shifts.
- **GSEA:** Gene Set Enrichment Analysis focused on Hallmark pathways (e.g., Estrogen Response, Cell Cycle).
- **Survival Analysis:** Kaplan-Meier plotting to link gene expression levels with patient clinical outcomes.

## Tools & Technologies
- **Language:** R
- **Bioconductor Suite:** DESeq2, TCGAbiolinks, clusterProfiler, SummarizedExperiment
- **Visualization:** EnhancedVolcano, pheatmap, survminer, ggplot2
- **AI Collaboration:** This pipeline was developed with the assistance of **Claude Code** and **Google Gemini** for code optimization and robust clinical data handling.

## Repository Structure
```text
TCGA-BRCA-Basal-vs-LumA/
├── analysis_pipeline.R     # Main R script containing the full analysis
├── README.md               # Project documentation
├── .gitignore              # Files to be excluded from GitHub (data/ and large objects)
└── results/
    ├── figures/            # Visualizations (Volcano plots, Heatmaps, Survival curves)
    └── tables/             # CSV files of DEGs and Enrichment results
