# pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
#                glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
#                apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
#                factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
#                genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
#                GSEAmining, ggrepel, progress, mnormt, psych, igraph, 
#                reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
#                mikropml, glmnet, scales, stats, caret, nnet, pROC)
# pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
#                glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
#                apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
#                factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
#                genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
#                clusterProfiler, gghalves, GSEABase, GSEAmining, ggrepel, progress, mnormt, psych, igraph, 
#                reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
#                mikropml, glmnet, scales, stats, caret, nnet, pROC)

pacman::p_load(
  apeglm, biomaRt, boot, car, caret, clusterProfiler, 
  coefplot, colorspace, cowplot, DESeq2, devtools, 
  dplyr, edgeR, factoextra, fmsb, genefilter, ggforce, 
  gghalves, gglasso, ggplot2, ggpubr, ggrepel, ggvenn, 
  glmnet, gplots, grid, gridExtra, GSEABase, GSEAmining, 
  GSVA, Hmisc, igraph, limma, magrittr, MatrixGenerics, 
  matrixStats, mikropml, mixOmics, mnormt, msigdbr, 
  nnet, pheatmap, plyr, pROC, progress, psych, 
  randomForest, RColorBrewer, rdist, reactome.db, 
  reshape2, ROCR, scales, stats, tidyverse, VennDiagram
)

library(patchwork)
library(tibble)
library(Seurat)
library(EnhancedVolcano)
library(stringr)
library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(gghalves)

# 1. Load the Data ----
# Organize Data
setwd("C:/Users/17343/Desktop/AntiCD3Monitoring/Data")
getwd()

# Metadata Import
meta_batch <- read.csv("C:/Users/17343/Desktop/AntiCD3Monitoring/Data/AntiCD3_metadata.csv", 
                       header = TRUE, 
                       check.names = FALSE)
meta_batch <- as.data.frame(meta_batch)

# Counts Import
counts_batch <- read.csv("C:/Users/17343/Desktop/AntiCD3Monitoring/Data/AntiCD3_gene_counts_annot.csv", 
                       header = TRUE, 
                       check.names = FALSE)
counts_batch <- as.data.frame(counts_batch)
counts_batch <- na.omit(counts_batch)

#Remove duplicate names
counts_batch <- counts_batch[!duplicated(counts_batch$external_gene_name),]
genes <- counts_batch$external_gene_name
rownames(counts_batch) <- genes
counts_batch$external_gene_name <- NULL

#Create new column with high-level analysis groups
meta_batch$Outcomes <- meta_batch$Group


# 2.Preprocessing and Cleaning ----
getwd()
setwd("C:/Users/17343/Desktop/IsletTransplantRejection/Code")

# Select columns in combined_counts that match the remaining sample names in meta_combined
counts_batch <- counts_batch[, meta_batch$Samples]  # Ensure Sample_IDs match column names in combined_counts

source("AllFunctions.R")
AntiCD3Counts <- flexiDEG.function1(counts_batch, meta_batch, # Genes in rows 
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

rows_to_remove <- grep("^Gm[0-9]", rownames(AntiCD3Counts))

# Remove those rows from case1_f1
AntiCD3Counts <- AntiCD3Counts[-rows_to_remove, ]

# Subsetting data for day 14
meta_batch <- meta_batch %>% filter(Day ==14)
AntiCD3Counts <- AntiCD3Counts[, colnames(AntiCD3Counts) %in% meta_batch$Samples]

# Saving case1_f1 dataframe as a CSV file
write.csv(AntiCD3Counts, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/Anti_CD3/AntiCD3_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)

# Saving meta_combined dataframe as a CSV file
write.csv(meta_batch, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/Anti_CD3/AntiCD3_Metadata_ITx.csv", row.names = FALSE)

# Color palettes
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
colSide <- flexiDEG.colors(meta_batch)
unique_colSide <- unique(colSide)

# 3.Create Universal DESqEQ Object ----
AntiCD3Counts <- as.matrix(AntiCD3Counts)
storage.mode(AntiCD3Counts) <- "integer"

dds_AntiCD3 <- DESeqDataSetFromMatrix(AntiCD3Counts, meta_batch,
                                              design = ~ 1)   # dummy design for now

saveRDS(dds_AntiCD3, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/RObjects/dds_AntiCD3_master.rds")

q_cut  <- 0.10
fc_cut <- 1
library(EnhancedVolcano)
make_keyvals_fdr_fc <- function(df, q = 0.10, fc = 1,
                                col_up = "#D62728", col_down = "#1F77B4", col_ns = "gray70") {
  # valid stats (for coloring); everything else becomes NS
  ok   <- !is.na(df$padj) & !is.na(df$log2FoldChange)
  
  up   <- ok & df$padj < q & df$log2FoldChange >=  fc
  down <- ok & df$padj < q & df$log2FoldChange <= -fc
  
  # default color + label
  key   <- rep(col_ns, nrow(df))
  label <- rep("Not significant", nrow(df))
  
  # overwrite where significant
  key[up]   <- col_up
  key[down] <- col_down
  
  label[up]   <- paste0("Upregulated")
  label[down] <- paste0("Downregulated")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}

library(dplyr)
pick_labels <- function(df, q = 0.10, fc = 1, topN = 30) {
  idx <- which(!is.na(df$padj) & !is.na(df$log2FoldChange) &
                 df$padj < q & abs(df$log2FoldChange) >= fc)
  if (length(idx) == 0) return(character(0))
  ord <- order(df$padj[idx], -abs(df$log2FoldChange[idx]), na.last = NA)  # tie-break by |LFC|
  labs <- df$gene[idx][ord]
  labs <- make.unique(labs)  # avoid dup labels
  labs[seq_len(min(topN, length(labs)))]
}

# 4a.Compare Anti-CD3 vs Isotype at Day 14 (Design ~ Treatment; Subset Day 14) ----
# subset samples
sel <- colData(dds_AntiCD3)$Treatment %in% c("Isotype","Anti-CD3")
dds_AntiCD3_AntiCD3VsIso <- dds_AntiCD3[, sel]

dds_AntiCD3_AntiCD3VsIso$Treatment <- factor(dds_AntiCD3_AntiCD3VsIso$Treatment,
                                             levels = c("Isotype", 
                                                        "Anti-CD3"))
levels(dds_AntiCD3_AntiCD3VsIso$Treatment)

cd <- as.data.frame(colData(dds_AntiCD3_AntiCD3VsIso))
# Basic sanity
lapply(cd[, c("Treatment")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Treatment")], function(x) any(is.na(x)))
# Model matrix rank
mm <- model.matrix(~ Treatment, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank

table(cd$Group,cd$Day)
#             14
# Isotype      5
# Resistant    8
# Sensitive    9

keep <- rowSums(counts(dds_AntiCD3_AntiCD3VsIso) >= 10) >= 5  # Smallest number of samples in a treatment group
dds_AntiCD3_AntiCD3VsIso <- dds_AntiCD3_AntiCD3VsIso[keep, ]

design(dds_AntiCD3_AntiCD3VsIso) <- ~ Treatment
dds_AntiCD3_AntiCD3VsIso <- DESeq(dds_AntiCD3_AntiCD3VsIso)
design(dds_AntiCD3_AntiCD3VsIso)
resultsNames(dds_AntiCD3_AntiCD3VsIso)
res_AntiCD3VsIso <- results(dds_AntiCD3_AntiCD3VsIso,
                            name = "Treatment_Anti.CD3_vs_Isotype")

d14_AntiCD3VsIso <- as.data.frame(res_AntiCD3VsIso); d14_AntiCD3VsIso$gene <- rownames(res_AntiCD3VsIso)
write.csv(d14_AntiCD3VsIso, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_AntiCD3VsIso.csv", row.names = TRUE)

keyvals_d14_AntiCD3VsIso  <- make_keyvals_fdr_fc(d14_AntiCD3VsIso)
selLab_d14_AntiCD3VsIso  <- pick_labels(d14_AntiCD3VsIso,  q_cut, fc_cut, 30)
xmax_d14_AntiCD3VsIso  <- max(2, ceiling(max(abs(d14_AntiCD3VsIso$log2FoldChange),  na.rm=TRUE)))
dev.new(width = 10, height = 10)
EnhancedVolcano(
  d14_AntiCD3VsIso,
  lab           = d14_AntiCD3VsIso$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Anti-CD3 vs Isotype Treatment",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_AntiCD3VsIso$padj<0.10 & abs(d14_AntiCD3VsIso$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14_AntiCD3VsIso, xmax_d14_AntiCD3VsIso),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d14_AntiCD3VsIso,
  legendPosition= "right",
  selectLab     = selLab_d14_AntiCD3VsIso
)

# 4b.Compare the Resistant vs Sensitive Group in the anti-CD3 treatment (Again for Day 14 only) ----
sel <- colData(dds_AntiCD3)$Group %in% c("Sensitive","Resistant")
dds_AntiCD3_ResistantVsSensitive <- dds_AntiCD3[, sel]

dds_AntiCD3_ResistantVsSensitive$Group <- factor(dds_AntiCD3_ResistantVsSensitive$Group,
                                                 levels = c("Sensitive", 
                                                            "Resistant"))
levels(dds_AntiCD3_ResistantVsSensitive$Group)

cd <- as.data.frame(colData(dds_AntiCD3_ResistantVsSensitive))
# Basic sanity
lapply(cd[, c("Group")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Group")], function(x) any(is.na(x)))
# Model matrix rank
mm <- model.matrix(~ Group, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank

table(cd$Group,cd$Day)
#           14
# Sensitive  9
# Resistant  8

keep <- rowSums(counts(dds_AntiCD3_ResistantVsSensitive) >= 10) >= 8  # Smallest number of samples in a treatment group
dds_AntiCD3_ResistantVsSensitive <- dds_AntiCD3_ResistantVsSensitive[keep, ]

#Design Formula
design(dds_AntiCD3_ResistantVsSensitive) <- ~ Group
dds_AntiCD3_ResistantVsSensitive <- DESeq(dds_AntiCD3_ResistantVsSensitive)
resultsNames(dds_AntiCD3_ResistantVsSensitive)
res_ResistantVsSensitive <- results(dds_AntiCD3_ResistantVsSensitive,
                            name = "Group_Resistant_vs_Sensitive")

d14_ResistantVsSensitive <- as.data.frame(res_ResistantVsSensitive); d14_ResistantVsSensitive$gene <- rownames(res_ResistantVsSensitive)
write.csv(d14_ResistantVsSensitive, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_ResistantVsSensitive.csv", row.names = TRUE)

keyvals_d14_ResistantVsSensitive  <- make_keyvals_fdr_fc(d14_ResistantVsSensitive)
selLab_d14_ResistantVsSensitive  <- pick_labels(d14_ResistantVsSensitive,  q_cut, fc_cut, 30)
xmax_d14_ResistantVsSensitive  <- max(2, ceiling(max(abs(d14_ResistantVsSensitive$log2FoldChange),  na.rm=TRUE)))
dev.new(width = 10, height = 10)
EnhancedVolcano(
  d14_ResistantVsSensitive,
  lab           = d14_ResistantVsSensitive$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Resistant vs Sensitive Group",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_ResistantVsSensitive$padj<0.10 & abs(d14_ResistantVsSensitive$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14_ResistantVsSensitive, xmax_d14_ResistantVsSensitive),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d14_ResistantVsSensitive,
  legendPosition= "right",
  selectLab     = selLab_d14_ResistantVsSensitive
)

# 5.GSEA Analysis----
# Paths to your saved results
AntiCD3VsIso_path <- "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_AntiCD3VsIso.csv"
ResistantVsSensitive_path <- "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_ResistantVsSensitive.csv"

# Import
AntiCD3VsIso <- read.csv(AntiCD3VsIso_path,  row.names = 1)
ResistantVsSensitive <- read.csv(ResistantVsSensitive_path, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_AntiCD3VsIso <- AntiCD3VsIso$stat;  names(lfc_vector_AntiCD3VsIso) <- rownames(AntiCD3VsIso)
lfc_vector_ResistantVsSensitive <- ResistantVsSensitive$stat; names(lfc_vector_ResistantVsSensitive) <- rownames(ResistantVsSensitive)

# Drop NAs
lfc_vector_AntiCD3VsIso <- lfc_vector_AntiCD3VsIso[!is.na(lfc_vector_AntiCD3VsIso)]
lfc_vector_ResistantVsSensitive <- lfc_vector_ResistantVsSensitive[!is.na(lfc_vector_ResistantVsSensitive)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_AntiCD3VsIso <- sort(lfc_vector_AntiCD3VsIso,  decreasing = TRUE)
lfc_vector_ResistantVsSensitive <- sort(lfc_vector_ResistantVsSensitive, decreasing = TRUE)

# --- Collect each set and convert into 2-column (gs_name, gene_symbol) ---

# C8
CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
mm_c8_df <- data.frame(
  gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
  gene_symbol = unlist(mm_c8_sets)
)

# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)

# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)

# Combines C8, Hallmark, and KEGG genes
mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)

# Isotype Vs AntiCD3 Treatment
gsea_results_AntiCD3VsIso <- GSEA(
  geneList      = lfc_vector_AntiCD3VsIso,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_AntiCD3VsIso_df <- as.data.frame(gsea_results_AntiCD3VsIso)

# Sensitive vs Resistant Group
gsea_results_ResistantVsSensitive <- GSEA(
  geneList      = lfc_vector_ResistantVsSensitive,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ResistantVsSensitive_df <- as.data.frame(gsea_results_ResistantVsSensitive)

# Full results Isotype Vs AntiCD3 Treatment
write.csv(gsea_results_AntiCD3VsIso_df,
          "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/GSEAResults_AntiCD3VsIso.csv",
          row.names = FALSE)

# Full results Sensitive Vs Resistant Group
write.csv(gsea_results_ResistantVsSensitive_df,
          "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/GSEAResults_ResistantVsSensitive.csv",
          row.names = FALSE)
