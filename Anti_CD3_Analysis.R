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
#setwd("C:/Users/17343/Desktop/AntiCD3Monitoring/Data")
setwd("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/")
getwd()

# Metadata Import
# meta_batch <- read.csv("C:/Users/17343/Desktop/AntiCD3Monitoring/Data/AntiCD3_metadata.csv", 
#                        header = TRUE, 
#                        check.names = FALSE)
meta_batch <- read.csv("AntiCD3_metadata.csv", 
                       header = TRUE, 
                       check.names = FALSE)
meta_batch <- as.data.frame(meta_batch)

# Counts Import
# counts_batch <- read.csv("C:/Users/17343/Desktop/AntiCD3Monitoring/Data/AntiCD3_gene_counts_annot.csv", 
#                        header = TRUE, 
#                        check.names = FALSE)
counts_batch <- read.csv("AntiCD3_gene_counts_annot.csv", 
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

# Select columns in combined_counts that match the remaining sample names in meta_combined
counts_batch <- counts_batch[, meta_batch$Samples]  # Ensure Sample_IDs match column names in combined_counts

source("AllFunctions.R")
AntiCD3Counts <- flexiDEG.function1(counts_batch, meta_batch, # Genes in rows 
                                            convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                            batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0, 0

remove_pattern <- "^Gm[0-9]|^AC[0-9]|^AL[0-9]|^AI[0-9]|^AW[0-9]|^AF[0-9]|^BB[0-9]|^BC[0-9]|^CT[0-9]|^CAAA|^BX[0-9]|^CN[0-9]|^CR[0-9]|^C[0-9]{4,}|^Olfr"
rows_to_remove <- grep(remove_pattern, rownames(AntiCD3Counts)) #Remove Gm genes
# Remove those rows from case1_f1
AntiCD3Counts <- AntiCD3Counts[-rows_to_remove, ]

# # Subsetting data for day 14
# meta_batch <- meta_batch %>% filter(Day ==14)
# AntiCD3Counts <- AntiCD3Counts[, colnames(AntiCD3Counts) %in% meta_batch$Samples]
# 
# # Saving case1_f1 dataframe as a CSV file
# write.csv(AntiCD3Counts, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/Anti_CD3/AntiCD3_Raw_Counts_ITx_Filtered.csv", row.names = TRUE)
# 
# # Saving meta_combined dataframe as a CSV file
# write.csv(meta_batch, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/Anti_CD3/AntiCD3_Metadata_ITx.csv", row.names = FALSE)

# Color palettes
# coul <- colorRampPalette(brewer.pal(11, "RdBu"))(100) # Palette for gene heatmaps
# coul_gsva <- colorRampPalette(brewer.pal(11, "PRGn"))(100) # Palette for gsva heatmaps
# colSide <- flexiDEG.colors(meta_batch)
# unique_colSide <- unique(colSide)

# 3.Create Universal DESqEQ Object ----
AntiCD3Counts <- as.matrix(AntiCD3Counts)
storage.mode(AntiCD3Counts) <- "integer"

dds_AntiCD3 <- DESeqDataSetFromMatrix(AntiCD3Counts, meta_batch,
                                              design = ~ 1)   # dummy design for now

#saveRDS(dds_AntiCD3, file = "C:/Users/17343/Desktop/AntiCD3Monitoring/Data/RObjects/dds_AntiCD3_master.rds")
saveRDS(dds_AntiCD3, file = "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Robject/dds_AntiCD3_master.rds")

q_cut  <- 0.10
fc_cut <- 1
make_keyvals_fdr_fc_aCD3_Iso <- function(df, q = 0.10, fc = 1,
                                col_up = "#1F77B4", col_down = "#D62728", col_ns = "gray70") {
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
  
  label[up]   <- paste0("Anti-CD3")
  label[down] <- paste0("Isotype")
  
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

# 4 Anti-CD3 vs Isotype at Day 14 (Design ~ Treatment; Subset Day 14) ----

# subset samples
sel <- colData(dds_AntiCD3)$Treatment %in% c("Isotype","Anti-CD3") & colData(dds_AntiCD3)$Day %in% c(14)
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

table(cd$Treatment,cd$Day)
#             14
# Isotype      5
# Resistant    8
# Sensitive    9

keep <- rowSums(counts(dds_AntiCD3_AntiCD3VsIso) >= 10) >= 5  # Smallest number of samples in a treatment group
dds_AntiCD3_AntiCD3VsIso <- dds_AntiCD3_AntiCD3VsIso[keep, ]

design(dds_AntiCD3_AntiCD3VsIso) <- ~ Treatment
dds_AntiCD3_AntiCD3VsIso <- DESeq(dds_AntiCD3_AntiCD3VsIso)
resultsNames(dds_AntiCD3_AntiCD3VsIso)
res_AntiCD3VsIso <- results(dds_AntiCD3_AntiCD3VsIso,
                            name = "Treatment_Anti.CD3_vs_Isotype")

d14_AntiCD3VsIso <- as.data.frame(res_AntiCD3VsIso); d14_AntiCD3VsIso$gene <- rownames(res_AntiCD3VsIso)
#write.csv(d14_AntiCD3VsIso, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_AntiCD3VsIso.csv", row.names = TRUE)
write.csv(d14_AntiCD3VsIso, "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/AntiCD3VsIsotype/DESEQResults_Day14_AntiCD3VsIso.csv", row.names = TRUE)

keyvals_d14_AntiCD3VsIso  <- make_keyvals_fdr_fc_aCD3_Iso(d14_AntiCD3VsIso)
selLab_d14_AntiCD3VsIso  <- c(
  # Up in Isotype
  "Ctla4", "Icos", "Foxp3", "Itk", "Cd2",
  #"Trat1", "Sit1", 
  "Rorc", "Klrg1",
  #"Klrb1a", "Klrb1c",
  #"Pou2af1", "Mzb1", 
  "Jchain",
  "Igkc", "Ighg1", 
  #"Ighm", "Iglc1", "Iglv1",
  #"Trgc1", 
  "Trdc",
  #"Blk",
  "Gzma", "Gzmd",
  #Up in Anti-CD3
  # Cytokines / chemokines
  "Il1b",  
  # Innate immune / myeloid / inflammation
  "Mpo", "Mmp9", 
  "S100a8", "Pglyrp1", "Cd177",
  "Adgre4",
  "Siglece", "Cd300e",
  # Immune regulation / checkpoint-like / metabolic immune modulators
  "Ido1"
)
xmax_d14_AntiCD3VsIso  <- max(2, ceiling(max(abs(d14_AntiCD3VsIso$log2FoldChange),  na.rm=TRUE)))
#dev.new(width = 10, height = 10)
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
  xlim          = c(-6, 6),
  ylim = c(0,13),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 5,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = Inf,
  colCustom     = keyvals_d14_AntiCD3VsIso,
  legendPosition= "right",
  selectLab     = selLab_d14_AntiCD3VsIso
)

## PCA-DE Genes----

# 1. Get significant DE genes
sig_genes_AntiCD3VsIso <- d14_AntiCD3VsIso %>%
  dplyr::filter(
    !is.na(padj),
    padj <= 0.10,
    abs(log2FoldChange) > 1
  ) %>%
  dplyr::pull(gene)

length(sig_genes_AntiCD3VsIso)
# 2. Variance-stabilized expression

vsd_AntiCD3VsIso <- vst(dds_AntiCD3_AntiCD3VsIso, blind = FALSE)
expr_AntiCD3VsIso <- assay(vsd_AntiCD3VsIso)

# Keep only DE genes
expr_AntiCD3VsIso_sig <- expr_AntiCD3VsIso[rownames(expr_AntiCD3VsIso) %in% sig_genes_AntiCD3VsIso, ]

# Remove zero-variance genes if any
expr_AntiCD3VsIso_sig <- expr_AntiCD3VsIso_sig[
  apply(expr_AntiCD3VsIso_sig, 1, var, na.rm = TRUE) > 0,
  , drop = FALSE
]

# 3. PCA on samples

mat_pca <- t(expr_AntiCD3VsIso_sig)
pca <- prcomp(mat_pca, scale. = TRUE)

# Sample metadata
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Treatment = colData(dds_AntiCD3_AntiCD3VsIso)$Treatment
)

# Make sure treatment names are exactly what you want in plot
pca_df$Treatment <- factor(pca_df$Treatment, levels = c("Isotype", "Anti-CD3"))

# 4. Define colors and shapes

group_colors <- c(
  "Isotype" = "#B23A48",
  "Anti-CD3" = "#2A6F97"
)

group_shapes <- c(
  "Isotype" = 15,
  "Anti-CD3" = 17
)

# 5. Plot PCA

ggplot(pca_df, aes(PC1, PC2, color = Treatment, shape = Treatment)) +
  geom_point(aes(fill = Treatment), size = 5, stroke = 1.2) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.2,
    aes(fill = Treatment),
    show.legend = FALSE,
    level = 0.7
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA of DE Genes: Anti-CD3 vs Isotype",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )



## GSEA Analysis ----

# Paths to your saved results
#AntiCD3VsIso_path <- "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_AntiCD3VsIso.csv"
AntiCD3VsIso_path <- "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/AntiCD3VsIsotype/DESEQResults_Day14_AntiCD3VsIso.csv"
# Import
AntiCD3VsIso <- read.csv(AntiCD3VsIso_path,  row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_AntiCD3VsIso <- AntiCD3VsIso$stat;  names(lfc_vector_AntiCD3VsIso) <- rownames(AntiCD3VsIso)

# Drop NAs
lfc_vector_AntiCD3VsIso <- lfc_vector_AntiCD3VsIso[!is.na(lfc_vector_AntiCD3VsIso)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_AntiCD3VsIso <- sort(lfc_vector_AntiCD3VsIso,  decreasing = TRUE)

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
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_AntiCD3VsIso_df <- as.data.frame(gsea_results_AntiCD3VsIso)

write.csv(gsea_results_AntiCD3VsIso_df,
          "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/AntiCD3VsIsotype/GSEAResults_AntiCD3VsIso.csv",
          row.names = FALSE)

library(enrichplot)

# Check exact pathway name
gsea_results_AntiCD3VsIso_df %>%
  dplyr::filter(ID == "TRAVAGLINI_LUNG_CD4_MEMORY_EFFECTOR_T_CELL")

# Running enrichment plot
gseaplot2(
  gsea_results_AntiCD3VsIso,
  geneSetID = "TRAVAGLINI_LUNG_CD4_MEMORY_EFFECTOR_T_CELL",
  title = "CD4 MEMORY EFFECTOR T CELL",
  color = "#B23A48",
  base_size = 16,
  subplots = 1:3,
  pvalue_table = TRUE
)


# Check exact pathway name
gsea_results_AntiCD3VsIso_df %>%
  dplyr::filter(ID == "GAUTAM_EYE_IRIS_CILIARY_BODY_CYTOTOXIC_T_CELLS")

# Running enrichment plot
gseaplot2(
  gsea_results_AntiCD3VsIso,
  geneSetID = "GAUTAM_EYE_IRIS_CILIARY_BODY_CYTOTOXIC_T_CELLS",
  title = "CYTOTOXIC_T_CELLS",
  color = "#B23A48",
  base_size = 16,
  subplots = 1:3,
  pvalue_table = TRUE
)

gseaplot2(
  gsea_results_AntiCD3VsIso,
  geneSetID = "GAUTAM_EYE_IRIS_CILIARY_BODY_CYTOTOXIC_T_CELLS",
  title = "Cytotoxic T Cells Geneset",
  base_size = 25
)
gseaplot2(
  gsea_results_AntiCD3VsIso,
  geneSetID = "GAUTAM_EYE_IRIS_CILIARY_BODY_CYTOTOXIC_T_CELLS",
  title = "Cytotoxic T Cell Signature",
  color = "#1F3A5F",   # deep navy
  base_size = 25
)

Cytotoxic_TCells <- gsea_results_AntiCD3VsIso[
  gsea_results_AntiCD3VsIso$Description == "GAUTAM_EYE_IRIS_CILIARY_BODY_CYTOTOXIC_T_CELLS",
]

# 5.Compare the Resistant vs Sensitive Group in the anti-CD3 treatment (Again for Day 14 only) ----
sel <- colData(dds_AntiCD3)$Group %in% c("Sensitive","Resistant")  & colData(dds_AntiCD3)$Day %in% c(14)
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
#write.csv(d14_ResistantVsSensitive, "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/DESEQResults_Day14_ResistantVsSensitive.csv", row.names = TRUE)
write.csv(d14_ResistantVsSensitive, "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/ResistantVsSensitive/DESEQResults_Day14_ResistantVsSensitive.csv", row.names = TRUE)

getwd()
keyvals_d14_ResistantVsSensitive  <- make_keyvals_fdr_fc(d14_ResistantVsSensitive)
selLab_d14_ResistantVsSensitive  <- pick_labels(d14_ResistantVsSensitive,  q_cut, fc_cut, 30)
xmax_d14_ResistantVsSensitive  <- max(2, ceiling(max(abs(d14_ResistantVsSensitive$log2FoldChange),  na.rm=TRUE)))
dev.off()

## Elastic Net Feature Selection ----

sig_genes_ResistantVsSensitive <- d14_ResistantVsSensitive$gene[
  !is.na(d14_ResistantVsSensitive$pvalue) &
    d14_ResistantVsSensitive$pvalue <= 0.05 &
    abs(d14_ResistantVsSensitive$log2FoldChange) >= 0.5
]

length(sig_genes_ResistantVsSensitive)

vsd <- vst(dds_AntiCD3_ResistantVsSensitive, blind = FALSE)
mat <- assay(vsd)
cd <- as.data.frame(colData(dds_AntiCD3_ResistantVsSensitive))

# expression matrix for glmnet: samples x genes
x_ResVsSen <- t(mat[sig_genes_ResistantVsSensitive, , drop = FALSE])

# binary outcome
y_ResVsSen <- ifelse(cd$Group == "Resistant", 1, 0)

# Find best alpha with LOOCV
set.seed(123)
alpha_grid <- seq(0, 1, by = 0.1)
cv_summary <- data.frame()

for (a in alpha_grid) {
  cvfit <- cv.glmnet(
    x = x_ResVsSen,
    y = y_ResVsSen,
    family = "binomial",
    alpha = a,
    foldid = 1:length(y_ResVsSen),   # LOOCV
    type.measure = "deviance",
    standardize = TRUE
  )
  
  cv_summary <- rbind(
    cv_summary,
    data.frame(
      alpha = a,
      cv_error = min(cvfit$cvm),
      lambda_min = cvfit$lambda.min,
      lambda_1se = cvfit$lambda.1se
    )
  )
}

cv_summary[order(cv_summary$cv_error), ]
cv_summary

best_alpha <- cv_summary$alpha[which.min(cv_summary$cv_error)]

# Stability selection
set.seed(123)
n_iter <- 1000
n_samples <- nrow(x_ResVsSen)


class1_idx <- which(y_ResVsSen == unique(y_ResVsSen)[1])
class2_idx <- which(y_ResVsSen == unique(y_ResVsSen)[2])

selected_list <- vector("list", n_iter)

for (i in 1:n_iter) {
  
  # Subsample ~80% of samples each time
  
  idx1 <- sample(class1_idx, size = round(0.8 * length(class1_idx)))
  idx2 <- sample(class2_idx, size = round(0.8 * length(class2_idx)))
  idx  <- c(idx1, idx2)
  
  x_sub <- x_ResVsSen[idx, ]
  y_sub <- y_ResVsSen[idx]
  
  cvfit <- cv.glmnet(
    x_sub,
    y_sub,
    family = "binomial",
    alpha = best_alpha,
    foldid = 1:length(y_sub),
    type.measure = "deviance",
    standardize = TRUE
  )
  
  coef_mat <- coef(cvfit, s = "lambda.1se")
  genes <- rownames(coef_mat)[coef_mat[,1] != 0]
  genes <- setdiff(genes, "(Intercept)")
  
  selected_list[[i]] <- genes
}

#Select Genes which appear atleast 70% of times
all_genes <- colnames(x_ResVsSen)

freq <- sapply(all_genes, function(g) {
  mean(sapply(selected_list, function(s) g %in% s))
})

freq_table <- data.frame(
  gene = names(freq),
  selection_frequency = freq
)
freq_table <- freq_table[order(freq_table$selection_frequency, decreasing = TRUE), ]
stable_genes <- subset(freq_table, selection_frequency >= 0.7) #Select genes appearing atleast 70 percent of time
stable_genes


#Plot Heatmp and PCA using Selected Genes 
stable_gene_names <- stable_genes$gene
mat_stable <- mat[stable_gene_names, , drop = FALSE]
mat_scaled <- t(scale(t(mat_stable)))
annotation_col <- data.frame(
  Group = cd$Group
)
rownames(annotation_col) <- colnames(mat_scaled)


annotation_colors <- list(
  Group = c(
    "Resistant" = "darkred",   # Scarlet
    "Sensitive" = "navyblue"     # green
  )
)
dev.new()
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 10,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101),  
  main = "Anti-CD3 Resistant Vs Sensitive Signature"
)

#PCA Analysis
mat_pca <- t(mat_stable)
pca <- prcomp(mat_pca, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = cd$Group
)

# Define colors & shapes
group_colors <- c(    "Resistant" = "darkred",   # Scarlet
                      "Sensitive" = "navyblue")     # green 
group_shapes <- c("Resistant" = 16, "Sensitive" = 24)

ggplot(pca_df, aes(PC1, PC2, color = Group,, shape = Group)) +
  geom_point(aes(fill = Group), size = 5, stroke = 1.2) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_classic(base_size = 18) +
  labs(
    title = "PCA: Anti-CD3 Resistant vs Sensitive",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1],1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2],1), "%)")
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )




## GSEA Analysis----
# Paths to your saved results
ResistantVsSensitive_path <- "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/ResistantVsSensitive/DESEQResults_Day14_ResistantVsSensitive.csv"

# Import

ResistantVsSensitive <- read.csv(ResistantVsSensitive_path, row.names = 1)

# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_ResistantVsSensitive <- ResistantVsSensitive$stat; names(lfc_vector_ResistantVsSensitive) <- rownames(ResistantVsSensitive)

# Drop NAs
lfc_vector_ResistantVsSensitive <- lfc_vector_ResistantVsSensitive[!is.na(lfc_vector_ResistantVsSensitive)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_ResistantVsSensitive <- sort(lfc_vector_ResistantVsSensitive, decreasing = TRUE)

# --- Collect each set and convert into 2-column (gs_name, gene_symbol) ---

# # C8
# CellTypeMSigDB_gene_sets <- msigdbr(species="Mus musculus", category="C8")
# mm_c8_sets <- split(CellTypeMSigDB_gene_sets$gene_symbol, CellTypeMSigDB_gene_sets$gs_name)
# mm_c8_df <- data.frame(
#   gs_name = rep(names(mm_c8_sets), sapply(mm_c8_sets, length)),
#   gene_symbol = unlist(mm_c8_sets)
# )
# 
# # Hallmark
# hallmark <- msigdbr(species = "Mus musculus", category  = "H")
# mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
# mm_h_df <- data.frame(
#   gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
#   gene_symbol = unlist(mm_h_sets)
# )
# 
# # KEGG
# kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
# mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
# mm_kegg_df <- data.frame(
#   gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
#   gene_symbol = unlist(mm_kegg_sets)
# )
# 
# # Combines C8, Hallmark, and KEGG genes
# mm_all_df <- rbind(mm_c8_df, mm_h_df, mm_kegg_df)


# Sensitive vs Resistant Group
gsea_results_ResistantVsSensitive <- GSEA(
  geneList      = lfc_vector_ResistantVsSensitive,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_ResistantVsSensitive_df <- as.data.frame(gsea_results_ResistantVsSensitive)

# Full results Isotype Vs AntiCD3 Treatment
# write.csv(gsea_results_AntiCD3VsIso_df,
#           "C:/Users/17343/Desktop/AntiCD3Monitoring/Results/GSEAResults_AntiCD3VsIso.csv",
#           row.names = FALSE)
write.csv(gsea_results_ResistantVsSensitive_df,
          "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Scaffold Bulk RNA/Results/ResistantVsSensitive/GSEAResults_ResistantVsSensitive.csv",
          row.names = FALSE)




immune_master_Rejection_vs_Sensitive <- unique(c(
  #Resistant
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_ENDOCYTOSIS",
  "KEGG_RIBOSOME",
  "FAN_OVARY_CL13_MONOCYTE_MACROPHAGE",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "TRAVAGLINI_LUNG_CD4_MEMORY_EFFECTOR_T_CELL",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "TRAVAGLINI_LUNG_CD4_NAIVE_T_CELL",
  "TRAVAGLINI_LUNG_PROLIFERATING_NK_T_CELL",
  #Sensitive
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HE_LIM_SUN_FETAL_LUNG_C2_APOE_POS_M2_MACROPHAGE_CELL",
  "DESCARTES_FETAL_LIVER_MYELOID_CELLS",
  "HALLMARK_BILE_ACID_METABOLISM",
  "KEGG_PEROXISOME",
  "KEGG_ABC_TRANSPORTERS",
  "HALLMARK_BILE_ACID_METABOLISM",
  "KEGG_LYSINE_DEGRADATION",
  "KEGG_ECM_RECEPTOR_INTERACTION"
))

# Helper to standardize clusterProfiler GSEA columns
coerce_gsea_tbl <- function(df, day_label){
  # try common column names: Description/ID/pathway/setName; p.adjust/padj
  gs  <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else if ("pathway" %in% names(df)) df$pathway else if ("setName" %in% names(df)) df$setName else rownames(df)
  pad <- if ("p.adjust"   %in% names(df)) df$p.adjust   else if ("padj" %in% names(df)) df$padj else df$pval
  tibble(
    gs_name = as.character(gs),
    NES     = as.numeric(df$NES),
    padj    = as.numeric(pad),
    Day     = day_label
  )
}

plot_df <- gsea_results_ResistantVsSensitive_df %>%
  filter(ID %in% immune_master_Rejection_vs_Sensitive) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    neglog10_padj = ifelse(is.infinite(neglog10_padj), NA, neglog10_padj)
  ) %>%
  arrange(NES) %>%
  mutate(
    ID = factor(ID, levels = ID)
  )


ggplot(plot_df, aes(x = NES, y = ID, size = neglog10_padj, color = NES)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](adjusted~italic(p))), range = c(3, 10)) +
  scale_color_gradient2(
    low = "navyblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    name = "NES"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic(base_size = 16) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Anti-CD3 Resistant vs Sensitive"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


