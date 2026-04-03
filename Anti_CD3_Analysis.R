pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

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

# 3.Create Univeral DESqEQ Object ----
AntiCD3Counts <- as.matrix(AntiCD3Counts)
storage.mode(AntiCD3Counts) <- "integer"

dds_AntiCD3 <- DESeqDataSetFromMatrix(AntiCD3Counts, meta_batch,
                                              design = ~ 1)   # dummy design for now

saveRDS(dds_AntiCD3, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/RObjects/dds_AntiCD3_master.rds")

# 4a.Compare Anti-CD3 vs Isotype at Day 14 (Design ~ Treatment; Subset Day 14) ----
# subset samples
sel <- colData(dds_AntiCD3)$Treatment %in% c("Isotype","Anti-CD3")
dds_AntiCD3_IsoVsAntiCD3 <- dds_AntiCD3[, sel]

dds_AntiCD3_IsoVsAntiCD3$Treatment <- factor(dds_AntiCD3_IsoVsAntiCD3$Treatment,
                                             levels = c("Isotype", 
                                                        "Anti-CD3"))
levels(dds_AntiCD3_IsoVsAntiCD3$Treatment)

cd <- as.data.frame(colData(dds_AntiCD3_IsoVsAntiCD3))
# Basic sanity
lapply(cd[, c("Treatment")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Treatment")], function(x) any(is.na(x)))
# Model matrix rank
mm <- model.matrix(~ Treatment, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank

# This might not be correct
table(cd$Group,cd$Day)
#             14
# Isotype      5
# Resistant    8
# Sensitive    9

keep <- rowSums(counts(dds_AntiCD3_IsoVsAntiCD3) >= 10) >= 5  # Smallest number of samples in a treatment group
dds_AntiCD3_IsoVsAntiCD3 <- dds_AntiCD3_IsoVsAntiCD3[keep, ]

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
  
  label[up]   <- paste0("Allogeneic-Upregulated")
  label[down] <- paste0("Syngeneic-Upregulated")
  
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

design(dds_AntiCD3_IsoVsAntiCD3) <- ~ Treatment
dds_AntiCD3_IsoVsAntiCD3 <- DESeq(dds_AntiCD3_IsoVsAntiCD3)
design(dds_AntiCD3_IsoVsAntiCD3)
resultsNames(dds_AntiCD3_IsoVsAntiCD3)
res_IsoVsAntiCD3 <- results(dds_AntiCD3_IsoVsAntiCD3,
                                             name = "Treatment_Anti.CD3_vs_Isotype")

d14_IsoVsAntiCD3 <- as.data.frame(res_IsoVsAntiCD3); d14_IsoVsAntiCD3$gene <- rownames(res_IsoVsAntiCD3)
write.csv(d14_IsoVsAntiCD3, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_IsoVsAntiCD3.csv", row.names = TRUE)

keyvals_d14_IsoVsAntiCD3  <- make_keyvals_fdr_fc(d14_IsoVsAntiCD3)
selLab_d14_IsoVsAntiCD3  <- pick_labels(d14_IsoVsAntiCD3,  q_cut, fc_cut, 30)
xmax_d14_IsoVsAntiCD3  <- max(2, ceiling(max(abs(d14_IsoVsAntiCD3$log2FoldChange),  na.rm=TRUE)))
dev.new(width = 10, height = 10)
EnhancedVolcano(
  d14_IsoVsAntiCD3,
  lab           = d14_IsoVsAntiCD3$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_IsoVsAntiCD3$padj<0.10 & abs(d14_IsoVsAntiCD3$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14_IsoVsAntiCD3, xmax_d14_IsoVsAntiCD3),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d14_IsoVsAntiCD3,
  legendPosition= "right",
  selectLab     = selLab_d14_IsoVsAntiCD3
)

# 4b.Compare the Resistant vs Sensitive Group in the anti-CD3 treatment (Again for Day 14 only) ----
sel <- colData(dds_AntiCD3)$Group %in% c("Sensitive","Resistant")
dds_AntiCD3_SensitiveVsResistant <- dds_AntiCD3[, sel]

dds_AntiCD3_SensitiveVsResistant$Group <- factor(dds_AntiCD3_SensitiveVsResistant$Group,
                                                 levels = c("Sensitive", 
                                                            "Resistant"))
levels(dds_AntiCD3_SensitiveVsResistant$Group)

cd <- as.data.frame(colData(dds_AntiCD3_SensitiveVsResistant))
# Basic sanity
lapply(cd[, c("Group")], function(x) table(x, useNA="ifany"))
# Check for NAs
sapply(cd[, c("Group")], function(x) any(is.na(x)))
# Model matrix rank
mm <- model.matrix(~ Group, data = cd)
qr(mm)$rank; ncol(mm)             # if rank < ncol(mm), not full rank

# This may not be correct
table(cd$Group,cd$Day)
#           14
# Sensitive  9
# Resistant  8

keep <- rowSums(counts(dds_AntiCD3_SensitiveVsResistant) >= 10) >= 8  # Smallest number of samples in a treatment group
dds_AntiCD3_SensitiveVsResistant <- dds_AntiCD3_SensitiveVsResistant[keep, ]

#Design Formula
design(dds_AntiCD3_SensitiveVsResistant) <- ~ Group
dds_AntiCD3_SensitiveVsResistant <- DESeq(dds_AntiCD3_SensitiveVsResistant)
resultsNames(dds_AntiCD3_SensitiveVsResistant)
res_SensitiveVsResistant <- results(dds_AntiCD3_SensitiveVsResistant,
                            name = "Group_Resistant_vs_Sensitive")

d14_SensitiveVsResistant <- as.data.frame(res_SensitiveVsResistant); d14_SensitiveVsResistant$gene <- rownames(res_SensitiveVsResistant)
write.csv(d14_SensitiveVsResistant, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_SensitiveVsResistant.csv", row.names = TRUE)

keyvals_d14_SensitiveVsResistant  <- make_keyvals_fdr_fc(d14_SensitiveVsResistant)
selLab_d14_SensitiveVsResistant  <- pick_labels(d14_SensitiveVsResistant,  q_cut, fc_cut, 30)
xmax_d14_SensitiveVsResistant  <- max(2, ceiling(max(abs(d14_SensitiveVsResistant$log2FoldChange),  na.rm=TRUE)))
dev.new(width = 10, height = 10)
EnhancedVolcano(
  d14_SensitiveVsResistant,
  lab           = d14_SensitiveVsResistant$gene,
  x             = "log2FoldChange",
  y             = "padj",
  pCutoff       = 0.10,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change)"),
  ylab          = expression("-log"[10]*"(FDR)"),
  title         = "Early vs Late Progessors",
  subtitle      = paste0("FDR ≤0.10 & |LFC| ≥1 (n=", sum(d14_SensitiveVsResistant$padj<0.10 & abs(d14_SensitiveVsResistant$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-xmax_d14_SensitiveVsResistant, xmax_d14_SensitiveVsResistant),
  ylim = c(0,5),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 15,
  colCustom     = keyvals_d14_SensitiveVsResistant,
  legendPosition= "right",
  selectLab     = selLab_d14_SensitiveVsResistant
)
