################################################################################
# TCGA-BRCA Differential Expression Analysis Pipeline
# Comparison: Basal vs LumA PAM50 Subtypes
# Author: Toluwanimi Adeyanju
# AI Partner: Claude, Google Gemini
# Date: 2025-01-01
################################################################################

# Set random seed for reproducibility
set.seed(123)

################################################################################
# 0. SETUP AND INSTALLATION
################################################################################

# Create directory structure
dir.create("data", showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Install required packages (run once)
install_if_missing <- function(packages, bioc = FALSE) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Bioconductor packages
bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "limma",
               "edgeR", "org.Hs.eg.db", "clusterProfiler", "EnhancedVolcano",
               "enrichplot", "GSEABase")
install_if_missing(bioc_pkgs, bioc = TRUE)

# CRAN packages
cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "dplyr", "stringr",
               "msigdbr", "scales")
install_if_missing(cran_pkgs, bioc = FALSE)

################################################################################
# 1. LOAD LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(EnhancedVolcano)
  library(enrichplot)
  library(dplyr)
  library(stringr)
  library(msigdbr)
  library(scales)
})

cat("✓ All libraries loaded successfully\n")

################################################################################
# 2. DOWNLOAD AND PREPARE DATA
################################################################################

data_file <- "data/TCGA_BRCA_SE.rda"

if (file.exists(data_file)) {
  cat("Loading existing data from", data_file, "\n")
  load(data_file)
} else {
  cat("Downloading TCGA-BRCA data...\n")

  query <- GDCquery(
    project       = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  GDCdownload(query)
  data <- GDCprepare(query)

  # Save for future use
  save(data, file = data_file)
  cat("✓ Data downloaded and saved\n")
}

################################################################################
# 3. QUALITY CONTROL AND DATA INSPECTION
################################################################################

cat("\n=== DATA SUMMARY ===\n")
cat("Total samples:", ncol(data), "\n")
cat("Total genes:", nrow(data), "\n")
cat("Assays available:", paste(assayNames(data), collapse = ", "), "\n")

# Check PAM50 subtype distribution
subtype_col <- "paper_BRCA_Subtype_PAM50"
cat("\nPAM50 Subtype Distribution:\n")
print(table(colData(data)[[subtype_col]], useNA = "always"))

# Check sequencing depth
seq_depth <- colSums(assay(data))
cat("\nSequencing Depth Summary:\n")
print(summary(seq_depth))

# Save QC report
sink("results/tables/QC_summary.txt")
cat("TCGA-BRCA Quality Control Report\n")
cat("Date:", as.character(Sys.Date()), "\n\n")
cat("Sample counts by subtype:\n")
print(table(colData(data)[[subtype_col]], useNA = "always"))
cat("\nSequencing depth:\n")
print(summary(seq_depth))
sink()

################################################################################
# 4. PREPARE DESeq2 OBJECT
################################################################################

# Filter samples with PAM50 annotation
data_filtered <- data[, !is.na(colData(data)[[subtype_col]])]
cat("\n✓ Filtered to", ncol(data_filtered), "samples with PAM50 annotation\n")

# Set factor levels
colData(data_filtered)[[subtype_col]] <- factor(
  colData(data_filtered)[[subtype_col]],
  levels = c("LumA", "LumB", "Her2", "Basal", "Normal")
)

# Build DESeqDataSet
dds <- DESeqDataSet(data_filtered, design = ~ paper_BRCA_Subtype_PAM50)

# Filter low-count genes (keep genes with >10 total counts)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
cat("✓ Retained", nrow(dds), "genes after filtering low counts\n")

# Run DESeq2
cat("\nRunning DESeq2 analysis...\n")
dds <- DESeq(dds)
cat("✓ DESeq2 analysis complete\n")

################################################################################
# 5. DIFFERENTIAL EXPRESSION: Basal vs LumA
################################################################################

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
res <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "LumA"))
res <- res[order(res$padj), ]

cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in Basal (L2FC > 1):",
    sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated in Basal (L2FC < -1):",
    sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")

################################################################################
# 6. GENE ANNOTATION
################################################################################

# Clean Ensembl IDs (remove version suffix)
clean_ids <- gsub("\\..*", "", rownames(res))

# Map to gene symbols
res$symbol <- mapIds(org.Hs.eg.db,
                     keys      = clean_ids,
                     column    = "SYMBOL",
                     keytype   = "ENSEMBL",
                     multiVals = "first")

# Save full results table
write.csv(as.data.frame(res),
          "results/tables/DESeq2_Basal_vs_LumA_full.csv",
          row.names = TRUE)

# Save significant genes only
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes),
          "results/tables/DESeq2_Basal_vs_LumA_significant.csv",
          row.names = TRUE)

cat("✓ Results saved to results/tables/\n")

################################################################################
# 7. VISUALIZATION: VOLCANO PLOT
################################################################################

cat("\nGenerating volcano plot...\n")

p_volcano <- EnhancedVolcano(
  res,
  lab              = res$symbol,
  x                = "log2FoldChange",
  y                = "padj",
  xlab             = bquote(~Log[2]~ "fold change"),
  ylab             = bquote(~-Log[10]~ "adjusted p-value"),
  title            = "Basal vs LumA (PAM50)",
  subtitle         = "DESeq2 results",
  pCutoff          = 0.05,
  FCcutoff         = 1.5,
  pointSize        = 2.0,
  labSize          = 4.0,
  colAlpha         = 0.8,
  legendPosition   = 'right',
  legendLabSize    = 10,
  legendIconSize   = 4.0
)

ggsave("results/figures/01_volcano_plot.png",
       plot = p_volcano,
       width = 12,
       height = 10,
       dpi = 300)

cat("✓ Volcano plot saved\n")

################################################################################
# 8. VISUALIZATION: HEATMAP OF TOP DE GENES
################################################################################

cat("Generating heatmap...\n")

# Top 50 DE genes by adjusted p-value
topGenes <- head(order(res$padj), 50)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)[topGenes, ]

# Sample annotation
annotation_col <- data.frame(
  Subtype = colData(data_filtered)$paper_BRCA_Subtype_PAM50
)
rownames(annotation_col) <- colnames(expr_matrix)

# Add gene symbols to row names
gene_ids <- gsub("\\..*", "", rownames(expr_matrix))
symbols <- mapIds(org.Hs.eg.db,
                  keys      = gene_ids,
                  column    = "SYMBOL",
                  keytype   = "ENSEMBL",
                  multiVals = "first")
rownames(expr_matrix) <- ifelse(!is.na(symbols), symbols, rownames(expr_matrix))

# Generate heatmap
png("results/figures/02_heatmap_top50.png", width = 3000, height = 2500, res = 300)
pheatmap(
  expr_matrix,
  annotation_col           = annotation_col,
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  fontsize_row             = 8,
  scale                    = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color                    = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main                     = "Top 50 DE Genes (Basal vs LumA)"
)
dev.off()

cat("✓ Heatmap saved\n")

################################################################################
# 9. FUNCTIONAL ENRICHMENT: GO AND KEGG (ORA)
################################################################################

cat("\n=== OVER-REPRESENTATION ANALYSIS ===\n")

# Get significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
gene_ids_sig <- gsub("\\..*", "", rownames(sig_genes))

# Convert to Entrez IDs
entrez <- mapIds(org.Hs.eg.db,
                 keys      = gene_ids_sig,
                 column    = "ENTREZID",
                 keytype   = "ENSEMBL",
                 multiVals = "first")
entrez <- na.omit(entrez)

cat("Genes for enrichment:", length(entrez), "\n")

## GO Biological Process
ego <- enrichGO(
  gene          = entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

if (!is.null(ego) && nrow(ego) > 0) {
  p_go <- dotplot(ego, showCategory = 15, title = "GO Biological Processes") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("results/figures/03_GO_enrichment.png", p_go,
         width = 10, height = 8, dpi = 300)
  write.csv(as.data.frame(ego), "results/tables/GO_enrichment.csv", row.names = FALSE)
  cat("✓ GO enrichment complete\n")
}

## KEGG Pathways
ekegg <- enrichKEGG(
  gene         = entrez,
  organism     = "hsa",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.5
)

if (!is.null(ekegg) && nrow(ekegg) > 0) {
  # Wrap long descriptions
  ekegg@result$Description <- str_wrap(ekegg@result$Description, width = 50)

  p_kegg <- dotplot(ekegg, showCategory = 15, title = "KEGG Pathways") +
    theme(axis.text.y = element_text(size = 9, lineheight = 0.9))
  ggsave("results/figures/04_KEGG_enrichment.png", p_kegg,
         width = 10, height = 10, dpi = 300)
  write.csv(as.data.frame(ekegg), "results/tables/KEGG_enrichment.csv", row.names = FALSE)
  cat("✓ KEGG enrichment complete\n")
}

################################################################################
# 10. GENE SET ENRICHMENT ANALYSIS (GSEA) - GO
################################################################################

cat("\n=== GENE SET ENRICHMENT ANALYSIS ===\n")

# Prepare ranked gene list
ranked_res <- data.frame(
  log2FC = res$log2FoldChange,
  symbol = res$symbol
) %>%
  filter(!is.na(log2FC) & !is.na(symbol)) %>%
  group_by(symbol) %>%
  slice_max(abs(log2FC), n = 1) %>%
  ungroup()

geneList <- ranked_res$log2FC
names(geneList) <- ranked_res$symbol
geneList <- sort(geneList, decreasing = TRUE)

cat("Ranked gene list size:", length(geneList), "\n")

# Run GSEA for GO BP
gse_go <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  eps          = 0
)

if (!is.null(gse_go) && nrow(gse_go) > 0) {
  cat("✓ GSEA found", nrow(gse_go), "enriched GO terms\n")

  # Dotplot
  p_gse <- dotplot(gse_go, showCategory = 10, split = ".sign") +
    facet_grid(.~.sign) +
    theme(axis.text.y = element_text(size = 8))
  ggsave("results/figures/05_GSEA_GO_dotplot.png", p_gse,
         width = 12, height = 8, dpi = 300)

  # Save results
  write.csv(as.data.frame(gse_go), "results/tables/GSEA_GO_results.csv", row.names = FALSE)
}

################################################################################
# 11. GSEA - HALLMARK PATHWAYS
################################################################################

cat("\nRunning Hallmark GSEA...\n")

# Get Hallmark gene sets
h_geneset <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

# Run GSEA
gse_hallmark <- GSEA(
  geneList      = geneList,
  TERM2GENE     = h_geneset,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

if (!is.null(gse_hallmark) && nrow(gse_hallmark) > 0) {
  cat("✓ Found", nrow(gse_hallmark), "enriched Hallmark pathways\n")

  # Dotplot
  p_hall <- dotplot(gse_hallmark, showCategory = 15, x = "NES") +
    ggtitle("Hallmark Pathways: Basal vs LumA") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("results/figures/06_GSEA_Hallmark_dotplot.png", p_hall,
         width = 10, height = 8, dpi = 300)

  # GSEA plot for top pathway
  p_top <- gseaplot2(gse_hallmark,
                     geneSetID = 1,
                     title = gse_hallmark$Description[1],
                     pvalue_table = TRUE)
  ggsave("results/figures/07_GSEA_top_pathway.png", p_top,
         width = 10, height = 6, dpi = 300)

  # Compare opposing pathways (if both exist)
  if (nrow(gse_hallmark) >= 2) {
    p_comp <- gseaplot2(gse_hallmark,
                        geneSetID = 1:2,
                        title = "Top 2 Enriched Pathways")
    ggsave("results/figures/08_GSEA_comparison.png", p_comp,
           width = 10, height = 6, dpi = 300)
  }

  # Save results
  gse_hallmark_readable <- setReadable(gse_hallmark, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
  write.csv(as.data.frame(gse_hallmark_readable),
            "results/tables/GSEA_Hallmark_results.csv",
            row.names = FALSE)

  cat("✓ Hallmark GSEA complete\n")
}


################################################################################
# 12. CLINICAL SURVIVAL ANALYSIS (Robust Version)
################################################################################

cat("\n=== STARTING SURVIVAL ANALYSIS ===\n")

library(survival)
library(survminer)

# 1. Get Clinical Data
clinical_data <- as.data.frame(colData(dds))

# 2. Robustly find Time and Status columns
# TCGA naming varies; we search for the most likely candidates
death_col <- intersect(c("days_to_death", "days_to_death.demographic"), colnames(clinical_data))[1]
followup_col <- intersect(c("days_to_last_followup", "days_to_last_follow_up", "days_to_last_follow_up.demographic"), colnames(clinical_data))[1]

cat("Using column for death:", death_col, "\n")
cat("Using column for follow-up:", followup_col, "\n")

# 3. Create 'time' and 'status' safely
# If columns are missing, we stop and alert
if(is.na(followup_col)) stop("Could not find follow-up columns in clinical data.")

clinical_data$time <- as.numeric(ifelse(clinical_data$vital_status == "Dead", 
                                        clinical_data[[death_col]], 
                                        clinical_data[[followup_col]]))

clinical_data$status <- ifelse(clinical_data$vital_status == "Dead", 1, 0)

# 4. Extract Expression for top gene
top_gene_id <- rownames(res)[1]
top_gene_symbol <- res$symbol[1]
gene_expr <- assay(vsd)[top_gene_id, ]

# 5. Create survival dataframe
surv_df <- data.frame(
  time = clinical_data$time,
  status = clinical_data$status,
  expression = gene_expr
)

# 6. Filter out invalid survival data (Negative days or NAs)
surv_df <- surv_df %>% filter(!is.na(time) & time > 0 & !is.na(status))

# 7. Median Split
surv_df$group <- ifelse(surv_df$expression > median(surv_df$expression), "High", "Low")

# 8. Fit and Plot
fit <- survfit(Surv(time, status) ~ group, data = surv_df)

p_surv <- ggsurvplot(
  fit, 
  data = surv_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#E41A1C", "#377EB8"),
  title = paste("Survival Impact of", top_gene_symbol),
  legend.labs = c("High Expr", "Low Expr"),
  xlab = "Time (Days)"
)

# Save
png("results/figures/09_survival_top_gene.png", width = 2400, height = 2400, res = 300)
print(p_surv)
dev.off()

cat("✓ Survival plot saved successfully\n")



################################################################################
# 13. BATCH SURVIVAL ANALYSIS (TOP 5 GENES)
################################################################################

cat("\n=== STARTING BATCH SURVIVAL ANALYSIS ===\n")

# 1. Define how many genes to test
num_genes <- 5
top_genes_indices <- 1:num_genes

# 2. Extract and clean clinical data once
clinical_data <- as.data.frame(colData(dds))
death_col <- intersect(c("days_to_death", "days_to_death.demographic"), colnames(clinical_data))[1]
followup_col <- intersect(c("days_to_last_followup", "days_to_last_follow_up"), colnames(clinical_data))[1]

clinical_data$time <- as.numeric(ifelse(clinical_data$vital_status == "Dead", 
                                        clinical_data[[death_col]], 
                                        clinical_data[[followup_col]]))
clinical_data$status <- ifelse(clinical_data$vital_status == "Dead", 1, 0)

# 3. Loop through the top genes
for (i in top_genes_indices) {
  
  gene_id <- rownames(res)[i]
  gene_symbol <- res$symbol[i]
  
  if(is.na(gene_symbol)) gene_symbol <- gene_id # Fallback to ID if symbol is missing
  
  cat(paste0("Processing [", i, "/", num_genes, "]: "), gene_symbol, "\n")
  
  # Prepare data for this specific gene
  surv_df <- data.frame(
    time = clinical_data$time,
    status = clinical_data$status,
    expression = assay(vsd)[gene_id, ]
  ) %>% filter(!is.na(time) & time > 0 & !is.na(status))
  
  # Stratify by median
  surv_df$group <- ifelse(surv_df$expression > median(surv_df$expression), "High", "Low")
  
  # Fit Survival Model
  fit <- survfit(Surv(time, status) ~ group, data = surv_df)
  
  # Generate Plot
  p <- ggsurvplot(
    fit, 
    data = surv_df,
    pval = TRUE,
    risk.table = TRUE,
    title = paste("Survival: ", gene_symbol),
    subtitle = paste("Comparison: Basal vs LumA Subset"),
    legend.labs = c("High", "Low"),
    palette = "Set1"
  )
  
  # Save plot with gene name in filename
  file_name <- paste0("results/figures/survival_", i, "_", gene_symbol, ".png")
  png(file_name, width = 2000, height = 2000, res = 300)
  print(p)
  dev.off()
  
  # Clean up memory for this iteration
  rm(surv_df, fit, p)
  gc()
}

cat("\n✓ All", num_genes, "survival plots saved to results/figures/\n")


################################################################################
# 14. COMBINED SURVIVAL PLOT (FIXED FOR COLUMN MISMATCH)
################################################################################

cat("\n=== GENERATING COMBINED SURVIVAL PLOT ===\n")

# 1. Get Clinical Data
clinical_data <- as.data.frame(colData(dds))

# 2. DYNAMICALLY FIND COLUMNS (The "Hunter" Logic)
# We search all column names for keywords to avoid the "0 rows" error
all_cols <- colnames(clinical_data)

# Find Vital Status
v_status_col <- all_cols[grep("vital_status", all_cols, ignore.case = TRUE)][1]

# Find Death Days (usually 'days_to_death')
death_days_col <- all_cols[grep("days_to_death", all_cols, ignore.case = TRUE)][1]

# Find Follow-up Days (usually 'days_to_last_followup')
followup_days_col <- all_cols[grep("days_to_last_follow", all_cols, ignore.case = TRUE)][1]

cat("Detected columns:\n", 
    "- Status:", v_status_col, "\n", 
    "- Death Days:", death_days_col, "\n", 
    "- Follow-up Days:", followup_days_col, "\n")

# 3. Handle NULLs - If a column is missing, create a dummy of NAs so the code doesn't crash
if(is.na(death_days_col)) { clinical_data$death_dummy <- NA; death_days_col <- "death_dummy" }
if(is.na(followup_days_col)) { clinical_data$follow_dummy <- NA; followup_days_col <- "follow_dummy" }

# 4. Create Time and Status
# We use 'pmax' to safely pick the available number, and as.numeric to ensure length
clinical_data$time <- as.numeric(ifelse(clinical_data[[v_status_col]] == "Dead", 
                                        clinical_data[[death_days_col]], 
                                        clinical_data[[followup_days_col]]))

clinical_data$status <- ifelse(clinical_data[[v_status_col]] == "Dead", 1, 0)

# 5. FINAL CHECK: Does 'time' have data?
if(all(is.na(clinical_data$time))) {
  stop("CRITICAL ERROR: Survival time could not be calculated. Please check colnames(colData(dds))")
}

# 6. Setup Plotting
num_genes <- 5
plot_list <- list()
top_genes_indices <- 1:num_genes

for (i in top_genes_indices) {
  gene_id <- rownames(res)[i]
  gene_symbol <- res$symbol[i]
  if(is.na(gene_symbol)) gene_symbol <- gene_id
  
  # Prepare data
  surv_df <- data.frame(
    time = clinical_data$time,
    status = clinical_data$status,
    expression = assay(vsd)[gene_id, ]
  )
  
  # Remove NAs BEFORE processing
  surv_df <- surv_df[complete.cases(surv_df), ]
  surv_df <- surv_df[surv_df$time > 0, ]
  
  if(nrow(surv_df) > 10) { # Only plot if we have enough samples
    surv_df$group <- ifelse(surv_df$expression > median(surv_df$expression), "High", "Low")
    fit <- survfit(Surv(time, status) ~ group, data = surv_df)
    
    p <- ggsurvplot(
      fit, data = surv_df, pval = TRUE, pval.size = 3,
      title = gene_symbol, legend = "none", palette = "Set1",
      ggtheme = theme_bw(base_size = 8)
    )
    plot_list[[gene_symbol]] <- p
  }
}

# 7. Stitch and Save
combined_res <- arrange_ggsurvplots(plot_list, print = FALSE, ncol = 3, nrow = 2)
ggsave("results/figures/10_combined_survival_top5.png", combined_res, width = 12, height = 8)

cat("✓ Combined plot saved successfully.\n")

################################################################################
# 15. SESSION INFO AND COMPLETION
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All results saved to results/ directory\n")
cat("- Figures: results/figures/\n")
cat("- Tables: results/tables/\n")

# Save session info
sink("results/sessionInfo.txt")
cat("TCGA-BRCA Analysis Session Information\n")
cat("Date:", as.character(Sys.Date()), "\n\n")
sessionInfo()
sink()

cat("\n✓ Session info saved to results/sessionInfo.txt\n")
cat("\n=== PIPELINE FINISHED SUCCESSFULLY ===\n")
