library(Seurat)
library(ggplot2)
library(sctransform)

# Mouse scRNA Seq----
mouse_aCD3<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/pre_analysis_mouse_pancreas_seurat.rds")
DefaultAssay(mouse_aCD3)
getwd()
mouse_aCD3[["HTO"]] <- NULL
ncol(mouse_aCD3)

length(Cells(mouse_aCD3))

dim(mouse_aCD3)

Layers(mouse_aCD3[["RNA"]])

sapply(Layers(mouse_aCD3[["RNA"]]), function(x)
  ncol(LayerData(mouse_aCD3, assay = "RNA", layer = x)))
x <- LayerData(mouse_aCD3, assay = "RNA", layer = "counts.1")


table(mouse_aCD3$HTO_classification.global, useNA = "ifany")
table(mouse_aCD3$hash.ID, useNA = "ifany")
table(mouse_aCD3@meta.data$orig.ident,mouse_aCD3$HTO_classification.global, useNA = "ifany")

# Hashed libraries: keep only HTO singlets
hashed_libs <- c(
  "GSA001", "GSA005", "GSA011", "GSA012", "GSA015", "GSA016",
  "GSA019", "GSA020", "GSA023", "GSA024"
)

# Single-sample libraries: no HTO demux available/needed
single_libs <- c("GSA002", "GSA006")

mouse_aCD3_rev <- subset(
  mouse_aCD3,
  subset =
    (
      orig.ident %in% hashed_libs &
        HTO_classification.global %in% c("Hashtag1", "Hashtag2")
    ) |
    (
      orig.ident %in% single_libs &
        nFeature_RNA > 200
    )
)
unique(mouse_aCD3_rev$orig.ident)
library(dplyr)
library(tibble)

library(dplyr)
library(tibble)

mouse_aCD3_rev$orig.ident <- as.character(mouse_aCD3_rev$orig.ident)
mouse_aCD3_rev$hash.ID <- as.character(mouse_aCD3_rev$hash.ID)

sample_map_hashed <- tribble(
  ~orig.ident, ~hash.ID,    ~MouseID,    ~Group,          ~Tissue,
  "GSA001",    "Hashtag1",  "aCD3_807",  "Sensitive",     "Pancreas",
  "GSA001",    "Hashtag2",  "aCD3_818",  "Sensitive",     "Pancreas",
  
  "GSA005",    "Hashtag1",  "aCD3_807",  "Sensitive",     "Blood",
  "GSA005",    "Hashtag2",  "aCD3_818",  "Sensitive",     "Blood",
  
  "GSA011",    "Hashtag1",  "aCD3_531",  "Sensitive",     "Pancreas",
  "GSA011",    "Hashtag2",  "aCD3_540",  "Sensitive",     "Pancreas",
  
  "GSA012",    "Hashtag1",  "aCD3_534",  "Resistant",     "Pancreas",
  "GSA012",    "Hashtag2",  "aCD3_535",  "Resistant",     "Pancreas",
  
  "GSA015",    "Hashtag1",  "aCD3_531",  "Sensitive",     "Blood",
  "GSA015",    "Hashtag2",  "aCD3_540",  "Sensitive",     "Blood",
  
  "GSA016",    "Hashtag1",  "aCD3_534",  "Resistant",     "Blood",
  "GSA016",    "Hashtag2",  "aCD3_535",  "Resistant",     "Blood",
  
  "GSA019",    "Hashtag1",  "RO_676",    "Recent-Onset",  "Pancreas",
  "GSA019",    "Hashtag2",  "RO_681",    "Recent-Onset",  "Pancreas",
  
  "GSA020",    "Hashtag1",  "RO_683",    "Recent-Onset",  "Pancreas",
  "GSA020",    "Hashtag2",  "RO_684",    "Recent-Onset",  "Pancreas",
  
  "GSA023",    "Hashtag1",  "RO_676",    "Recent-Onset",  "Blood",
  "GSA023",    "Hashtag2",  "RO_681",    "Recent-Onset",  "Blood",
  
  "GSA024",    "Hashtag1",  "RO_683",    "Recent-Onset",  "Blood",
  "GSA024",    "Hashtag2",  "RO_684",    "Recent-Onset",  "Blood"
)

sample_map_single <- tribble(
  ~orig.ident, ~MouseID,    ~Group,      ~Tissue,
  "GSA002",    "aCD3_801",  "Resistant", "Pancreas",
  "GSA006",    "aCD3_801",  "Resistant", "Blood"
)

meta_new <- mouse_aCD3_rev@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::left_join(sample_map_hashed, by = c("orig.ident", "hash.ID")) %>%
  dplyr::left_join(
    sample_map_single,
    by = "orig.ident",
    suffix = c("", ".single")
  ) %>%
  dplyr::mutate(
    MouseID = dplyr::coalesce(MouseID, MouseID.single),
    Group   = dplyr::coalesce(Group, Group.single),
    Tissue  = dplyr::coalesce(Tissue, Tissue.single)
  ) %>%
  dplyr::select(-MouseID.single, -Group.single, -Tissue.single) %>%
  tibble::column_to_rownames("cell")

mouse_aCD3_rev@meta.data <- meta_new
table(mouse_aCD3_rev$orig.ident, mouse_aCD3_rev$MouseID, useNA = "ifany")
table(mouse_aCD3_rev$Group, mouse_aCD3_rev$Tissue, useNA = "ifany")

mouse_aCD3_rev@meta.data %>%
  dplyr::filter(is.na(MouseID)) %>%
  dplyr::count(orig.ident, hash.ID)
VlnPlot(
  mouse_aCD3_rev,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "orig.ident",
  pt.size = 0
)

mouse_aCD3_rev<- subset(mouse_aCD3_rev, subset = nFeature_RNA > 200)

saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered.rds"
)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mouse_aCD3_rev[["percent.mt"]] <- PercentageFeatureSet(mouse_aCD3_rev, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(mouse_aCD3_rev, features = c("percent.mt"),pt.size = 0)
mouse_aCD3_rev <- subset(mouse_aCD3_rev, subset = percent.mt < 10)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(mouse_aCD3_rev, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_aCD3_rev, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

mouse_aCD3_rev <- SCTransform(mouse_aCD3_rev, vars.to.regress = "percent.mt", verbose = TRUE)


mouse_aCD3_rev <- RunPCA(mouse_aCD3_rev, verbose = TRUE)
mouse_aCD3_rev <- RunUMAP(mouse_aCD3_rev, dims = 1:30, verbose = TRUE)

mouse_aCD3_rev <- FindNeighbors(mouse_aCD3_rev, dims = 1:30, verbose = FALSE)
mouse_aCD3_rev <- FindClusters(mouse_aCD3_rev, verbose = FALSE)

DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T) 
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,group.by = "orig.ident") 

# Cell Type Annotation----
core_markers <- list(
  T_cells = c("Cd3d", "Cd3e", "Trac"),
  B_cells = c("Cd79a", "Cd79b", "Ms4a1", "Cd19"),
  NK_cells = c("Nkg7", "Ncr1", "Klrb1c", "Prf1"),
  DCs = c("Itgax", "Flt3", "Zbtb46", "Clec9a", "Xcr1","Siglech","Bst2"),
  Macrophages = c("Adgre1", "Csf1r", "Cd68" ,"C1qa", "C1qb", "C1qc","Mertk"),
  Monocytes = c("Ly6c2", "Lyz2"),
  Neutrophils = c("Retnlg", "Lcn2","S100a8", "S100a9"),
  Basophils = c("Mcpt8", "Fcer1a", "Ms4a2", "Cpa3"),
  RBCs = c( "Hba-a1",  "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2","Bpgm")
)
DotPlot(
  mouse_aCD3_rev,
  features = core_markers
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  )

mouse_aCD3_rev <- PrepSCTFindMarkers(mouse_aCD3_rev)

cluster18_markers <- FindMarkers(
  mouse_aCD3_rev,
  ident.1 = 18,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cluster18_markers <- cluster18_markers[order(cluster18_markers$avg_log2FC, decreasing = TRUE), ]

head(cluster18_markers, 50)
Idents(mouse_aCD3_rev)<-mouse_aCD3_rev$Idents
mouse_aCD3_rev$Idents<- Idents(mouse_aCD3_rev)
Idents(mouse_aCD3_rev)<- mouse_aCD3_rev$Idents 
mouse_aCD3_rev <- RenameIdents(mouse_aCD3_rev,
  `0` = "B cells",
  `1` = "T cells",
  `2` = "Neutrophils",
  `3` = "Neutrophils",
  `4` = "T cells",
  `5` = "B cells",
  `6` = "T cells",
  `7` = "T cells",
  `8` = "B cells",
  `9` = "T cells",
  `10` = "NK cells",
  `11` = "Monocytes/Macrophages",
  `12` = "Monocytes/Macrophages",
  `13` = "B cells",
  `14` = "T cells",
  `15` = "T cells",
  `16` = "Neutrophils",
  `17` = "Monocytes/Macrophages",
  `18` = "DCs",
  `19` = "Neutrophils",
  `20` = "B cells",
  `21` = "Monocytes/Macrophages",
  `22` = "T cells",
  `23` = "T cells",
  `24` = "B cells",
  `25` = "B cells",
  `26` = "pDCs",
  `27` = "T cells",
  `28` = "Basophils",
  `29` = "T cells",
  `30` = "RBCs", 
  `31` = "B cells",
  `32` = "T cells", 
  `33` = "Monocytes/Macrophages",
  `34` = "Neutrophils", 
  `35` = "B cells",
  `36` = "B cells",
  `37` = "Neutrophils" 
)

table(mouse_aCD3_rev$Idents,mouse_aCD3_rev$Tissue)
DimPlot(mouse_aCD3_rev, label = TRUE,label.box = T) 
saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v1.rds"
)

# replace with your Seurat object name
barcodes <- colnames(mouse_aCD3_rev)

head(barcodes)
barcodes_clean <- barcodes
# remove sample/library prefix
barcodes_clean <- sub(".*_", "", barcodes_clean)

# remove 10x suffix
barcodes_clean <- sub("-1$", "", barcodes_clean)

head(barcodes_clean)

length(barcodes_clean)
length(unique(barcodes_clean))
library(dplyr)
library(stringr)
library(purrr)

outdir <- "/Volumes/One_Touch/Pancreas_scRNAseq_TRUST4/barcodes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta_bc <- mouse_aCD3_rev@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    barcode_clean = cell %>%
      str_remove(".*_") %>%
      str_remove("-1$")
  )

barcode_tbl <- meta_bc %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    barcodes = list(unique(barcode_clean)),
    .groups = "drop"
  )

for (i in seq_len(nrow(barcode_tbl))) {
  writeLines(
    barcode_tbl$barcodes[[i]],
    file.path(outdir, paste0(barcode_tbl$orig.ident[i], "_barcodes.txt"))
  )
}
list.files(outdir)
sapply(list.files(outdir, full.names = TRUE), function(f) length(readLines(f)))

## T Cell Subtyping and Refining----
mouse_aCD3_rev<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v1.rds")
DimPlot(mouse_aCD3_rev)
colnames(mouse_aCD3_rev@meta.data)
mouse_T <- subset(
  mouse_aCD3_rev,
  idents = "T cells"
)
DefaultAssay(mouse_T) <- "RNA"
library(future)

options(
  future.globals.maxSize = 4 * 1024^3
)

mouse_T <- SCTransform(
  mouse_T,
  verbose = FALSE
)
Assays(mouse_T)
mouse_T <- RunPCA( mouse_T,  npcs = 50, verbose = FALSE)
ElbowPlot(mouse_T, ndims = 50)

mouse_T <- FindNeighbors(mouse_T, dims = 1:30)
mouse_T <- FindClusters(mouse_T,resolution = 0.5)

mouse_T <- RunUMAP( mouse_T, reduction = "pca",dims = 1:30)

DimPlot(mouse_T, label = TRUE, label.box = T) 

tcell_markers <- c(
  # General T-cell markers
  "Cd3d", "Cd3e", "Trac",
  # CD4 and CD8
  "Cd4", "Cd8a", "Cd8b1",
  # Naive / central-memory
  "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r", "Mald1",   "Bcl2",
  # Effector-memory / cytotoxic
  "Ccl5", "Nkg7", "Gzma", "Gzmk", "Gzmb", "Prf1", "Ifng", "Ctsw",
  # Regulatory T cells
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18",
  # Activated T cells
  "Cd69", "Icos", "Nr4a1", "Nr4a2", "Tnfrsf9",
  # Dysfunction / exhaustion
  "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox",
  # Proliferating
  "Mki67", "Top2a", "Stmn1",
  # Gamma-delta T cells
  "Trdc", "Trgc1", "Trgc2",
  # NK contamination
  "Ncr1", "Klrd1", "Tyrobp", "Fcerg1",
  # B Cell Contimation
  "Cd79a", "Ms4a1", "Pax5",
  #Neutrophils
  "Retnlg", "Lcn2","S100a8", "S100a9",
  #Macrophages
  "Adgre1", "Csf1r", "Cd68" ,"C1qa", "C1qb", "C1qc","Mertk",
  # ILC Like
  "Gata3","Rora",
  # APC or doublet check
  "H2-Aa", "H2-Ab1", "Cd74",
  "Lyz2", "Tyrobp", "Cd79a",
  "Ms4a1", "Epcam"
)

tcell_markers_present <- intersect(
  tcell_markers,
  rownames(mouse_T)
)

DotPlot(
  mouse_T,
  features = tcell_markers_present,
  assay = "SCT"
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  RotatedAxis()

Idents(mouse_T) <- "seurat_clusters"

cd4_cd8_genes <- intersect(
  c("Cd4", "Cd8a", "Cd8b1"),
  rownames(mouse_T)
)

cd4_cd8_genes
VlnPlot(
  mouse_T,
  features = cd4_cd8_genes,
  assay = "SCT",
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)

tbd_clusters <- c("9", "13", "15", "16", "17")
DefaultAssay(mouse_T) <- "SCT"
mouse_T <- PrepSCTFindMarkers(mouse_T)
Idents(mouse_T) <- "seurat_clusters"
tbd_markers <- lapply(
  tbd_clusters,
  function(cluster_id) {
    
    markers <- FindMarkers(
      mouse_T,
      ident.1 = cluster_id,
      assay = "SCT",
      only.pos = TRUE,
      min.pct = 0.10,
      logfc.threshold = 0.25
    )
    markers$gene <- rownames(markers)
    markers$cluster <- cluster_id
    markers
  }
)

tbd_markers <- dplyr::bind_rows(tbd_markers)
top_tbd_markers <- tbd_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  )

print(top_tbd_markers, n = Inf)
mouse_T$Idents<-Idents(mouse_T)
mouse_T <- RenameIdents(mouse_T,
                               `0` = "Cytotoxic CD8 T cells",#Cytotoxic
                               `1` = "CD4 T cells",
                               `2` = "CD4 T cells",
                               `3` = "NK cells",
                               `4` = "Memory/naive-like CD8 T cells",#Memory/naive
                               `5` = "CD4 T cells",
                               `6` = "CD4 T cells",
                               `7` = "Memory/naive-like CD8 T cells",#Memory/naive
                               `8` = "Tregs",
                               `9` = "B cells",
                               `10` = "Cytotoxic CD8 T cells",#Cytotoxic
                               `11` = "Gamma Delta T cells",
                               `12` = "Cytotoxic CD8 T cells",#Cytotoxic
                               `13` = "Memory/naive-like T cells",
                               `14` = "Cytotoxic CD8 T cells",#Cytotoxic
                               `15` = "Neutrophils",
                               `16` = "ILC2-like cells",
                               `17` = "Gamma Delta T cells"
)


DimPlot(mouse_T, label = TRUE, label.box = T) 
saveRDS(
  mouse_T,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/Tcells_mouse_aCD3_filtered_v1.rds"
)

# Change the naming in original object
tcell_annotations <- as.character(Idents(mouse_T))
names(tcell_annotations) <- colnames(mouse_T)

# Start with the original identities for all cells
mouse_aCD3_rev$celltype_refined <-
  as.character(Idents(mouse_aCD3_rev))

# Replace annotations for every cell present in mouse_T
mouse_aCD3_rev$celltype_refined[colnames(mouse_T)] <-
  tcell_annotations[colnames(mouse_T)]

# Check the updated annotations
table(
  mouse_aCD3_rev$celltype_refined,
  useNA = "ifany"
)

# Preserve original annotations
mouse_aCD3_rev$celltype_original <-
  as.character(Idents(mouse_aCD3_rev))
Idents(mouse_aCD3_rev) <- "celltype_refined"
#Save Version 2- Conatins T cell subsets and refined clusters
saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v2.rds"
)

# Verify
table(Idents(mouse_aCD3_rev))
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T)

sum(is.na(mouse_aCD3_rev$celltype_refined))
sum(table(mouse_aCD3_rev$celltype_refined)) == ncol(mouse_aCD3_rev)

## Macrophage Refining----
mouse_aCD3_rev<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v2.rds")
DimPlot(mouse_aCD3_rev)

mouse_aCD3_rev$celltype<-Idents(mouse_aCD3_rev)
mac_aCD3 <- subset(
  mouse_aCD3_rev,
  subset = celltype == "Monocytes/Macrophages"
)

options(
  future.globals.maxSize = 4 * 1024^3
)

mac_aCD3 <- SCTransform(
  mac_aCD3,
  verbose = FALSE
)
Assays(mac_aCD3)
mac_aCD3 <- RunPCA( mac_aCD3,  npcs = 50, verbose = FALSE)
ElbowPlot(mac_aCD3, ndims = 50)

mac_aCD3 <- FindNeighbors(mac_aCD3, dims = 1:30)
mac_aCD3 <- FindClusters(mac_aCD3,resolution = 0.5)

mac_aCD3 <- RunUMAP( mac_aCD3, reduction = "pca",dims = 1:30)

DimPlot(mac_aCD3, label = TRUE, label.box = T) 

cell_markers <- c(
  # General T-cell markers
  "Cd3d", "Cd3e", "Trac",
  # CD4 and CD8
  "Cd4", "Cd8a", "Cd8b1",
  # Naive / central-memory
  "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r", "Mald1",   "Bcl2",
  # Effector-memory / cytotoxic
  "Ccl5", "Nkg7", "Gzma", "Gzmk", "Gzmb", "Prf1", "Ifng", "Ctsw",
  # Regulatory T cells
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18",
  # Activated T cells
  "Cd69", "Icos", "Nr4a1", "Nr4a2", "Tnfrsf9",
  # Dysfunction / exhaustion
  "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox",
  # Proliferating
  "Mki67", "Top2a", "Stmn1",
  # Gamma-delta T cells
  "Trdc", "Trgc1", "Trgc2",
  # NK contamination
  "Ncr1", "Klrd1", "Tyrobp", "Fcerg1",
  # B Cell Contimation
  "Cd79a", "Ms4a1", "Pax5",
  #Neutrophils
  "Retnlg", "Lcn2","S100a8", "S100a9",
  #Macrophages
  "Adgre1", "Csf1r", "Cd68" ,"C1qa", "C1qb", "C1qc","Mertk",
  # ILC Like
  "Gata3","Rora",
  # APC or doublet check
  "H2-Aa", "H2-Ab1", "Cd74",
  "Lyz2"
)

cell_markers_present <- intersect(
  cell_markers,
  rownames(mac_aCD3)
)

DotPlot(
  mac_aCD3,
  features = cell_markers_present,
  assay = "SCT"
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  RotatedAxis()
mac_aCD3$Idents<-Idents(mac_aCD3)
mac_aCD3 <- RenameIdents(mac_aCD3,
                        `0` = "Monocytes/Macrophages",
                        `1` = "Monocytes/Macrophages",
                        `2` = "Monocytes/Macrophages",
                        `3` = "Monocytes/Macrophages",
                        `4` = "Monocytes/Macrophages",
                        `5` = "Monocytes/Macrophages",
                        `6` = "B cells",
                        `7` = "Monocytes/Macrophages",
                        `8` = "Monocytes/Macrophages",
                        `9` = "Monocytes/Macrophages",
                        `10` = "Monocytes/Macrophages",
                        `11` = "Neutrophils",
                        `12` = "Memory/naive-like CD8 T cells",
                        `13` = "Monocytes/Macrophages"
)

DimPlot(mac_aCD3, label = TRUE, label.box = T) 
saveRDS(
  mac_aCD3,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/Macrophages_mouse_aCD3_filtered_v1.rds"
)

# Change the naming in original object
mac_annotations <- as.character(Idents(mac_aCD3))
names(mac_annotations) <- colnames(mac_aCD3)

# Start with the original identities for all cells
mouse_aCD3_rev$celltype_refined <-
  as.character(Idents(mouse_aCD3_rev))

# Replace annotations for every cell present in mouse_T
mouse_aCD3_rev$celltype_refined[colnames(mac_aCD3)] <-
  mac_annotations[colnames(mac_aCD3)]

# Check the updated annotations
table(
  mouse_aCD3_rev$celltype_refined,
  useNA = "ifany"
)

# Preserve original annotations
mouse_aCD3_rev$celltype_original <-
  as.character(Idents(mouse_aCD3_rev))
Idents(mouse_aCD3_rev) <- "celltype_refined"
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T) 
#Save Version 3- Conatins Macrophage cell subsets and refined clusters
saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v3.rds"
)

# Verify
table(Idents(mouse_aCD3_rev))
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T)

sum(is.na(mouse_aCD3_rev$celltype_refined))
sum(table(mouse_aCD3_rev$celltype_refined)) == ncol(mouse_aCD3_rev)


## B cell Refining----
mouse_aCD3_rev<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v3.rds")

mouse_aCD3_rev$celltype<-Idents(mouse_aCD3_rev)
Bcells_aCD3 <- subset(
  mouse_aCD3_rev,
  subset = celltype == "B cells"
)

options(
  future.globals.maxSize = 4 * 1024^3
)

Bcells_aCD3 <- SCTransform(
  Bcells_aCD3,
  verbose = FALSE
)
Assays(Bcells_aCD3)
Bcells_aCD3 <- RunPCA( Bcells_aCD3,  npcs = 50, verbose = FALSE)
ElbowPlot(Bcells_aCD3, ndims = 50)

Bcells_aCD3 <- FindNeighbors(Bcells_aCD3, dims = 1:30)
Bcells_aCD3 <- FindClusters(Bcells_aCD3,resolution = 0.5)

Bcells_aCD3 <- RunUMAP( Bcells_aCD3, reduction = "pca",dims = 1:30)

DimPlot(Bcells_aCD3, label = TRUE, label.box = T) 

cell_markers <- c(
  # B Cell 
  "Cd79a", "Cd79b", "Ms4a1", "Cd19","Jchain",
  # General T-cell markers
  "Cd3d", "Cd3e", "Trac",
  # CD4 and CD8
  "Cd4", "Cd8a", "Cd8b1",
  # Naive / central-memory
  "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r", "Mald1",   "Bcl2",
  # Effector-memory / cytotoxic
  "Ccl5", "Nkg7", "Gzma", "Gzmk", "Gzmb", "Prf1", "Ifng", "Ctsw",
  # Regulatory T cells
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18",
  # Activated T cells
  "Cd69", "Icos",
  # Dysfunction / exhaustion
  "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox",
  # Proliferating
  "Mki67", "Top2a", "Stmn1",
  # Gamma-delta T cells
  "Trdc", "Trgc1", "Trgc2",
  # NK contamination
  "Ncr1", "Klrd1", "Tyrobp", "Fcerg1",
  #Neutrophils
  "Retnlg", "Lcn2","S100a8", "S100a9",
  #Macrophages
  "Adgre1", "Csf1r", "Cd68" ,"C1qa", "C1qb", "C1qc","Mertk",
  # ILC Like
  "Gata3","Rora",
  # APC or doublet check
  "H2-Aa", "H2-Ab1", "Cd74",
  "Lyz2"
)

cell_markers_present <- intersect(
  cell_markers,
  rownames(Bcells_aCD3)
)

DotPlot(
  Bcells_aCD3,
  features = cell_markers_present,
  assay = "SCT"
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  RotatedAxis()
Bcells_aCD3$Idents<-Idents(Bcells_aCD3)
Bcells_aCD3 <- RenameIdents(Bcells_aCD3,
                         `0` = "B cells",
                         `1` = "B cells",
                         `2` = "B cells",
                         `3` = "B cells",
                         `4` = "B cells",
                         `5` = "B cells",
                         `6` = "B cells",
                         `7` = "Memory/naive-like CD8 T cells",
                         `8` = "B cells",
                         `9` = "B cells",
                         `10` = "B cells",
                         `11` = "Monocytes/Macrophages",
                         `12` = "Cytotoxic CD8 T cells",
                         `13` = "Plasma cells",
                         `14` = "B cells",
                         `15` = "Plasma cells",
                         `16` = "B cells",
                         `17` = "Neutrophils",
                         `18` = "B cells",
                         `19` = "B cells"
)

DimPlot(Bcells_aCD3, label = TRUE, label.box = T) 
saveRDS(
  Bcells_aCD3,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/Bcells_mouse_aCD3_filtered_v1.rds"
)

# Change the naming in original object
Bcells_annotations <- as.character(Idents(Bcells_aCD3))
names(Bcells_annotations) <- colnames(Bcells_aCD3)

# Start with the original identities for all cells
mouse_aCD3_rev$celltype_refined <-
  as.character(Idents(mouse_aCD3_rev))

# Replace annotations for every cell present in mouse_T
mouse_aCD3_rev$celltype_refined[colnames(Bcells_aCD3)] <-
  Bcells_annotations[colnames(Bcells_aCD3)]

# Check the updated annotations
table(
  mouse_aCD3_rev$celltype_refined,
  useNA = "ifany"
)

# Preserve original annotations
mouse_aCD3_rev$celltype_original <-
  as.character(Idents(mouse_aCD3_rev))
Idents(mouse_aCD3_rev) <- "celltype_refined"
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T) 
#Save Version 4- Conatins B cell subsets and refined clusters
saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v4.rds"
)

# Verify
table(Idents(mouse_aCD3_rev))
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T)

sum(is.na(mouse_aCD3_rev$celltype_refined))
sum(table(mouse_aCD3_rev$celltype_refined)) == ncol(mouse_aCD3_rev)

##Neutrophil Refining----
mouse_aCD3_rev<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v4.rds")

mouse_aCD3_rev$celltype<-Idents(mouse_aCD3_rev)
Neutrophils_aCD3 <- subset(
  mouse_aCD3_rev,
  subset = celltype == "Neutrophils"
)

options(
  future.globals.maxSize = 4 * 1024^3
)

Neutrophils_aCD3 <- SCTransform(
  Neutrophils_aCD3,
  verbose = FALSE
)
Assays(Neutrophils_aCD3)
Neutrophils_aCD3 <- RunPCA( Neutrophils_aCD3,  npcs = 50, verbose = FALSE)
ElbowPlot(Neutrophils_aCD3, ndims = 50)

Neutrophils_aCD3 <- FindNeighbors(Neutrophils_aCD3, dims = 1:30)
Neutrophils_aCD3 <- FindClusters(Neutrophils_aCD3,resolution = 0.5)

Neutrophils_aCD3 <- RunUMAP( Neutrophils_aCD3, reduction = "pca",dims = 1:30)

DimPlot(Neutrophils_aCD3, label = TRUE, label.box = T) 

cell_markers <- c(
  #Neutrophils
  "Retnlg", "Lcn2","S100a8", "S100a9",
  # B Cell 
  "Cd79a", "Cd79b", "Ms4a1", "Cd19","Jchain",
  # General T-cell markers
  "Cd3d", "Cd3e", "Trac",
  # CD4 and CD8
  "Cd4", "Cd8a", "Cd8b1",
  # Naive / central-memory
  "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r", "Mald1",   "Bcl2",
  # Effector-memory / cytotoxic
  "Ccl5", "Nkg7", "Gzma", "Gzmk", "Gzmb", "Prf1", "Ifng", "Ctsw",
  # Regulatory T cells
  "Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18",
  # Activated T cells
  "Cd69", "Icos",
  # Dysfunction / exhaustion
  "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox",
  # Proliferating
  "Mki67", "Top2a", "Stmn1",
  # Gamma-delta T cells
  "Trdc", "Trgc1", "Trgc2",
  # NK contamination
  "Ncr1", "Klrd1", "Fcerg1",
  #Macrophages
  "Adgre1", "Csf1r", "Cd68" ,"C1qa", "C1qb", "C1qc","Mertk",
  # ILC Like
  "Gata3","Rora",
  # APC or doublet check
  "H2-Aa", "H2-Ab1", "Cd74",
  "Lyz2"
)

cell_markers_present <- intersect(
  cell_markers,
  rownames(Neutrophils_aCD3)
)

DotPlot(
  Neutrophils_aCD3,
  features = cell_markers_present,
  assay = "SCT"
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  RotatedAxis()
Neutrophils_aCD3$Idents<-Idents(Neutrophils_aCD3)
Neutrophils_aCD3 <- RenameIdents(Neutrophils_aCD3,
                            `0` = "Neutrophils",
                            `1` = "Neutrophils",
                            `2` = "Neutrophils",
                            `3` = "Neutrophils",
                            `4` = "Neutrophils",
                            `5` = "Neutrophils",
                            `6` = "Neutrophils",
                            `7` = "B cells",
                            `8` = "Neutrophils",
                            `9` = "Neutrophils",
                            `10` = "Memory/naive-like CD8 T cells",
                            `11` = "Neutrophils",
                            `12` = "Neutrophils",
                            `13` = "Monocytes/Macrophages",
                            `14` = "Neutrophils",
                            `15` = "Neutrophils",
                            `16` = "Neutrophils"
)

DimPlot(Neutrophils_aCD3, label = TRUE, label.box = T) 
saveRDS(
  Neutrophils_aCD3,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/Neutrophils_mouse_aCD3_filtered_v1.rds"
)

# Change the naming in original object
Neutrophils_annotations <- as.character(Idents(Neutrophils_aCD3))
names(Neutrophils_annotations) <- colnames(Neutrophils_aCD3)

# Start with the original identities for all cells
mouse_aCD3_rev$celltype_refined <-
  as.character(Idents(mouse_aCD3_rev))

# Replace annotations for every cell present in mouse_T
mouse_aCD3_rev$celltype_refined[colnames(Neutrophils_aCD3)] <-
  Neutrophils_annotations[colnames(Neutrophils_aCD3)]

# Check the updated annotations
table(
  mouse_aCD3_rev$celltype_refined,
  useNA = "ifany"
)

# Preserve original annotations
mouse_aCD3_rev$celltype_original <-
  as.character(Idents(mouse_aCD3_rev))
Idents(mouse_aCD3_rev) <- "celltype_refined"
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T) 
#Final Refined Version-v5
saveRDS(
  mouse_aCD3_rev,
  "/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v5.rds"
)

# Verify
table(Idents(mouse_aCD3_rev))
DimPlot(mouse_aCD3_rev, label = TRUE, label.box = T,repel = T)

sum(is.na(mouse_aCD3_rev$celltype_refined))
sum(table(mouse_aCD3_rev$celltype_refined)) == ncol(mouse_aCD3_rev)


#  UMAP and Annontation for Whole Dataset ----
p <- DimPlot(
  mouse_aCD3_rev,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  label.size = 4.8,
  pt.size = 0.18,
  raster = FALSE
) +
  theme_classic(base_size = 14) +
  theme(
    # Remove axes
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    # Remove legend
    legend.position = "none",
    
    # Clean background
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    # Remove grid
    panel.grid = element_blank()
  )

p


p1 <- DimPlot(
  mouse_aCD3_rev,
  reduction = "umap",
  group.by = 'Tissue',
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  label.size = 4.8,
  pt.size = 0.18,
  raster = FALSE
) +
  theme_classic(base_size = 14) +
  theme(
    # Remove axes
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    # Remove legend
    legend.position = "none",
    
    # Clean background
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    # Remove grid
    panel.grid = element_blank()
  )

p1


# Pancreas Macrophages ----
mac_pancreas <- subset(
  mouse_aCD3_rev,
  subset =
    celltype_refined == "Monocytes/Macrophages" &
    Tissue == "Pancreas"
)

mac_pancreas
table(mac_pancreas$Group)
table(mac_pancreas$MouseID, mac_pancreas$Group)

mac_sample_counts <- mac_pancreas@meta.data %>%
  dplyr::group_by(MouseID, Group) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Group, MouseID)

mac_sample_counts
#Keep mice with at least 10 macrophages
valid_mice <- mac_sample_counts %>%
  dplyr::filter(n_cells >= 10) %>%
  dplyr::pull(MouseID)

valid_mice
mac_pancreas_pb <- subset(
  mac_pancreas,
  subset = MouseID %in% valid_mice
)
table(
  mac_pancreas_pb$MouseID,
  mac_pancreas_pb$Group
)

#Aggregate raw RNA counts by mice
DefaultAssay(mac_pancreas_pb) <- "RNA"

mac_counts_pb <- AggregateExpression(
  mac_pancreas_pb,
  assays = "RNA",
  group.by = "MouseID",
  slot = "counts",
  return.seurat = FALSE
)$RNA
#Create mouse level metadata
mac_coldata <- mac_pancreas_pb@meta.data %>%
  dplyr::select(MouseID, Group) %>%
  dplyr::distinct() %>%
  as.data.frame()

mac_coldata
mac_coldata$PseudobulkID <- gsub(
  pattern = "_",
  replacement = "-",
  x = mac_coldata$MouseID
)

mac_coldata
rownames(mac_coldata) <- mac_coldata$PseudobulkID
mac_coldata <- mac_coldata[
  ,
  c("MouseID", "Group"),
  drop = FALSE
]

mac_coldata <- mac_coldata[
  colnames(mac_counts_pb),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(mac_coldata),
    colnames(mac_counts_pb)
  )
)

mac_coldata

#Set group factor
mac_coldata$Group <- factor(
  mac_coldata$Group,
  levels = c(
    "Recent-Onset",
    "Sensitive",
    "Resistant"
  )
)

mac_coldata

keep_genes <- rowSums(mac_counts_pb >= 10) >= 2

mac_counts_filtered <- mac_counts_pb[
  keep_genes,
  ,
  drop = FALSE
]

dds_mac <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(mac_counts_filtered),
  colData = mac_coldata,
  design = ~ Group
)

dds_mac <- DESeq2::DESeq(dds_mac)
vsd_mac <- vst(dds_mac,blind = FALSE)

## Whole Transcriptome----
pca_df <- plotPCA(
  vsd_mac,
  intgroup = "Group",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_df, "percentVar"),
  1
)

pca_df$MouseID <- mac_coldata[
  pca_df$name,
  "MouseID"
]
pca_df

pca_df$PlotGroup <- ifelse(
  pca_df$group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

pca_df$PlotGroup <- factor(
  pca_df$PlotGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

table(pca_df$group, pca_df$PlotGroup)

pca_df$MouseID <- mac_coldata[
  pca_df$name,
  "MouseID"
]

p_mac_pca_combined <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = PlotGroup,
    shape = PlotGroup
  )
) +
  stat_ellipse(
    aes(
      fill = PlotGroup,
      group = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.20,
    level = 0.70,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Recent-Onset" = 15,
      "Anti-CD3" = 16
    )
  ) +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_mac_pca_combined_ellipse <- p_mac_pca_combined +
  stat_ellipse(
    aes(
      group = PlotGroup,
      fill = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.2,
    level = 0.7,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  )

p_mac_pca_combined_ellipse

## IN Myeloid Anti-CD3 Signature----
sig_genes <- intersect(
  stable_gene_names,
  rownames(assay(vsd_mac))
)

length(sig_genes)
sig_genes
sig_mat <- assay(vsd_mac)[sig_genes, ]
pca_sig <- prcomp(
  t(sig_mat),
  center = TRUE,
  scale. = TRUE
)
percentVar_sig <- round(
  100 * pca_sig$sdev^2 / sum(pca_sig$sdev^2),
  1
)
pca_sig_df <- data.frame(
  PseudobulkID = rownames(pca_sig$x),
  PC1 = pca_sig$x[, 1],
  PC2 = pca_sig$x[, 2],
  stringsAsFactors = FALSE
)


pca_sig_df$MouseID <- mac_coldata[
  pca_sig_df$PseudobulkID,
  "MouseID"
]

pca_sig_df$Group <- mac_coldata[
  pca_sig_df$PseudobulkID,
  "Group"
]
pca_sig_df$PlotGroup <- ifelse(
  pca_sig_df$Group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

pca_sig_df$PlotGroup <- factor(
  pca_sig_df$PlotGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

p_mac_pca_combined_sig <- ggplot(
  pca_sig_df,
  aes(
    PC1,
    PC2,
    color = PlotGroup,
    shape = PlotGroup
  )
) +
  stat_ellipse(
    aes(
      fill = PlotGroup,
      group = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.20,
    level = 0.70,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Recent-Onset" = 15,
      "Anti-CD3" = 16
    )
  ) +
  labs(
    x = paste0("PC1 (", percentVar_sig[1], "%)"),
    y = paste0("PC2 (", percentVar_sig[2], "%)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_mac_pca_combined_ellipse_sig <- p_mac_pca_combined_sig +
  stat_ellipse(
    aes(
      group = PlotGroup,
      fill = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.2,
    level = 0.7,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  )

p_mac_pca_combined_ellipse_sig

## DE analyssis-----
dds_mac <- DESeq2::DESeq(dds_mac)
dds_mac_pancreas <- dds_mac
mac_counts_pb_pancreas <- mac_counts_pb
mac_coldata_pancreas <- mac_coldata

mac_coldata_pancreas$TreatmentGroup <- ifelse(
  mac_coldata_pancreas$Group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

mac_coldata_pancreas$TreatmentGroup <- factor(
  mac_coldata_pancreas$TreatmentGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

dds_mac_pancreas_combined <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds_mac_pancreas),
  colData = mac_coldata_pancreas,
  design = ~ TreatmentGroup
)

dds_mac_pancreas_combined <- DESeq2::DESeq(
  dds_mac_pancreas_combined
)

res_pancreas <- DESeq2::results(
  dds_mac_pancreas_combined,
  contrast = c(
    "TreatmentGroup",
    "Anti-CD3",
    "Recent-Onset"
  )
)

res_pancreas_df <- as.data.frame(res_pancreas) |>
  tibble::rownames_to_column("Gene") |>
  dplyr::select(
    Gene,
    pancreas_log2FC = log2FoldChange,
    pancreas_pvalue = pvalue,
    pancreas_padj = padj
  )

# Blood Monocytes ----
mac_blood <- subset(
  mouse_aCD3_rev,
  subset =
    celltype_refined == "Monocytes/Macrophages" &
    Tissue == "Blood"
)

mac_blood
table(mac_blood$Group)
table(mac_blood$MouseID, mac_blood$Group)

mac_sample_counts <- mac_blood@meta.data %>%
  dplyr::group_by(MouseID, Group) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Group, MouseID)

mac_sample_counts
#Keep mice with at least 10 monocytes
valid_mice <- mac_sample_counts %>%
  dplyr::filter(n_cells >= 10) %>%
  dplyr::pull(MouseID)

valid_mice
mac_blood_pb <- subset(
  mac_blood,
  subset = MouseID %in% valid_mice
)
table(
  mac_blood_pb$MouseID,
  mac_blood_pb$Group
)

#Aggregate raw RNA counts by mice
DefaultAssay(mac_blood_pb) <- "RNA"

mac_counts_pb <- AggregateExpression(
  mac_blood_pb,
  assays = "RNA",
  group.by = "MouseID",
  slot = "counts",
  return.seurat = FALSE
)$RNA
#Create mouse level metadata
mac_coldata <- mac_blood_pb@meta.data %>%
  dplyr::select(MouseID, Group) %>%
  dplyr::distinct() %>%
  as.data.frame()

mac_coldata
mac_coldata$PseudobulkID <- gsub(
  pattern = "_",
  replacement = "-",
  x = mac_coldata$MouseID
)

mac_coldata
rownames(mac_coldata) <- mac_coldata$PseudobulkID
mac_coldata <- mac_coldata[
  ,
  c("MouseID", "Group"),
  drop = FALSE
]

mac_coldata <- mac_coldata[
  colnames(mac_counts_pb),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(mac_coldata),
    colnames(mac_counts_pb)
  )
)

mac_coldata

#Set group factor
mac_coldata$Group <- factor(
  mac_coldata$Group,
  levels = c(
    "Recent-Onset",
    "Sensitive",
    "Resistant"
  )
)

mac_coldata

keep_genes <- rowSums(mac_counts_pb >= 10) >= 2

mac_counts_filtered <- mac_counts_pb[
  keep_genes,
  ,
  drop = FALSE
]

dds_mac <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(mac_counts_filtered),
  colData = mac_coldata,
  design = ~ Group
)

dds_mac <- DESeq2::DESeq(dds_mac)
vsd_mac <- vst(dds_mac,blind = FALSE)

## Whole Transcriptome----
pca_df <- plotPCA(
  vsd_mac,
  intgroup = "Group",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_df, "percentVar"),
  1
)

pca_df$MouseID <- mac_coldata[
  pca_df$name,
  "MouseID"
]
pca_df

pca_df$PlotGroup <- ifelse(
  pca_df$group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

pca_df$PlotGroup <- factor(
  pca_df$PlotGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

table(pca_df$group, pca_df$PlotGroup)

pca_df$MouseID <- mac_coldata[
  pca_df$name,
  "MouseID"
]

p_mac_pca_combined <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = PlotGroup,
    shape = PlotGroup
  )
) +
  stat_ellipse(
    aes(
      fill = PlotGroup,
      group = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.20,
    level = 0.70,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Recent-Onset" = 15,
      "Anti-CD3" = 16
    )
  ) +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_mac_pca_combined_ellipse <- p_mac_pca_combined +
  stat_ellipse(
    aes(
      group = PlotGroup,
      fill = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.2,
    level = 0.7,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  )

p_mac_pca_combined_ellipse

## IN Myeloid Anti-CD3 Signature----
sig_genes <- intersect(
  stable_gene_names,
  rownames(assay(vsd_mac))
)

length(sig_genes) #28
sig_genes
sig_mat <- assay(vsd_mac)[sig_genes, ]
pca_sig <- prcomp(
  t(sig_mat),
  center = TRUE,
  scale. = TRUE
)
percentVar_sig <- round(
  100 * pca_sig$sdev^2 / sum(pca_sig$sdev^2),
  1
)
pca_sig_df <- data.frame(
  PseudobulkID = rownames(pca_sig$x),
  PC1 = pca_sig$x[, 1],
  PC2 = pca_sig$x[, 2],
  stringsAsFactors = FALSE
)


pca_sig_df$MouseID <- mac_coldata[
  pca_sig_df$PseudobulkID,
  "MouseID"
]

pca_sig_df$Group <- mac_coldata[
  pca_sig_df$PseudobulkID,
  "Group"
]
pca_sig_df$PlotGroup <- ifelse(
  pca_sig_df$Group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

pca_sig_df$PlotGroup <- factor(
  pca_sig_df$PlotGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

p_mac_pca_combined_sig <- ggplot(
  pca_sig_df,
  aes(
    PC1,
    PC2,
    color = PlotGroup,
    shape = PlotGroup
  )
) +
  stat_ellipse(
    aes(
      fill = PlotGroup,
      group = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.20,
    level = 0.70,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5,
    stroke = 1
  ) +
  scale_color_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Recent-Onset" = 15,
      "Anti-CD3" = 16
    )
  ) +
  labs(
    x = paste0("PC1 (", percentVar_sig[1], "%)"),
    y = paste0("PC2 (", percentVar_sig[2], "%)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    legend.title = element_blank(),
    legend.position = "right"
  )

p_mac_pca_combined_ellipse_sig <- p_mac_pca_combined_sig +
  stat_ellipse(
    aes(
      group = PlotGroup,
      fill = PlotGroup
    ),
    geom = "polygon",
    alpha = 0.2,
    level = 0.7,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Recent-Onset" = "#D62728",
      "Anti-CD3" = "black"
    )
  )

p_mac_pca_combined_ellipse_sig

## DE analysis----
dds_mac <- DESeq2::DESeq(dds_mac)
dds_mac_blood <- dds_mac
mac_counts_pb_blood <- mac_counts_pb
mac_coldata_blood <- mac_coldata

mac_coldata_blood$TreatmentGroup <- ifelse(
  mac_coldata_blood$Group == "Recent-Onset",
  "Recent-Onset",
  "Anti-CD3"
)

mac_coldata_blood$TreatmentGroup <- factor(
  mac_coldata_blood$TreatmentGroup,
  levels = c("Recent-Onset", "Anti-CD3")
)

dds_mac_blood_combined <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds_mac_blood),
  colData = mac_coldata_blood,
  design = ~ TreatmentGroup
)

dds_mac_blood_combined <- DESeq2::DESeq(
  dds_mac_blood_combined
)

res_blood <- DESeq2::results(
  dds_mac_blood_combined,
  contrast = c(
    "TreatmentGroup",
    "Anti-CD3",
    "Recent-Onset"
  )
)

res_blood_df <- as.data.frame(res_blood) |>
  tibble::rownames_to_column("Gene") |>
  dplyr::select(
    Gene,
    blood_log2FC = log2FoldChange,
    blood_pvalue = pvalue,
    blood_padj = padj
  )

# Pan Tissue Analysis----
comparison_df <-
  d14_AntiCD3VsIso %>%
  dplyr::select(
    Gene = gene,
    scaffold_log2FC = log2FoldChange
  ) %>%
  inner_join(
    res_pancreas_df,
    by = "Gene"
  ) %>%
  inner_join(
    res_blood_df,
    by = "Gene"
  ) %>%
  filter(
    !is.na(scaffold_log2FC),
    !is.na(pancreas_log2FC),
    !is.na(blood_log2FC)
  )

nrow(comparison_df)

comparison_sig <-
  comparison_df %>%
  filter(
    Gene %in% stable_gene_names
  )

nrow(comparison_sig)

cor_pancreas <- cor.test(
  comparison_sig$scaffold_log2FC,
  comparison_sig$pancreas_log2FC,
  method = "spearman"
)

cor_blood <- cor.test(
  comparison_sig$scaffold_log2FC,
  comparison_sig$blood_log2FC,
  method = "spearman"
)

cor_pancreas_blood <- cor.test(
  comparison_sig$pancreas_log2FC,
  comparison_sig$blood_log2FC,
  method = "spearman"
)



cor_pancreas
cor_blood
cor_pancreas_blood


label_genes <- comparison_sig 

p_scaffold_pancreas <-
  
  ggplot(
    comparison_sig,
    aes(
      scaffold_log2FC,
      pancreas_log2FC
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_point(
    size = 3,
    colour = "#1f78b4"
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "black",
    se = FALSE
  ) +
  
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 4.5
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.5,
    label = paste0(
      "Spearman r = ",
      round(cor_pancreas$estimate,2),
      "\nP = ",
      signif(cor_pancreas$p.value,2)
    )
  ) +
  
  labs(
    x = expression(
      Scaffold~log[2]~FC~"(Anti-CD3 vs Isotype)"
    ),
    y = expression(
      Pancreas~log[2]~FC~"(Anti-CD3 vs Recent-Onset)"
    )
  ) +
  
  theme_classic(base_size = 16)

p_scaffold_pancreas

p_scaffold_blood <-
  
  ggplot(
    comparison_sig,
    aes(
      scaffold_log2FC,
      blood_log2FC
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_point(
    size = 3,
    colour = "#d95f02"
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "black",
    se = FALSE
  ) +
  
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 4.5
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.5,
    label = paste0(
      "Spearman r = ",
      round(cor_blood$estimate,2),
      "\nP = ",
      signif(cor_blood$p.value,2)
    )
  ) +
  
  labs(
    x =  expression(
      Scaffold~log[2]~FC~"(Anti-CD3 vs Isotype)"
    ),
    y = expression(Blood~log[2]~FC~"(Anti-CD3 vs Recent-Onset)")
  ) +
  
  theme_classic(base_size = 16)

p_scaffold_blood


p_pancreas_blood <-
  ggplot(
    comparison_sig,
    aes(
      pancreas_log2FC,
      blood_log2FC
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_point(
    size = 3,
    colour = "black"
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "black",
    se = FALSE
  ) +
  
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 4.5
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.5,
    label = paste0(
      "Spearman r = ",
      round(cor_pancreas_blood$estimate, 2),
      "\nP = ",
      signif(cor_pancreas_blood$p.value, 2)
    )
  ) +
  
  labs(
    x =  expression(
      Pancreas~log[2]~FC~"(Anti-CD3 vs Recent-Onset)"),
    y =  expression(
      Blood~log[2]~FC~"(Anti-CD3 vs Recent-Onset)")
  ) +
  
  theme_classic(base_size = 16)

p_pancreas_blood

# TCR Analysis Seq----
mouse_aCD3_rev<-readRDS("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/mouse_aCD3_filtered_v1.rds")
DimPlot(mouse_aCD3_rev)

setwd("/Users/jyotirmoyroy/Desktop/Anti-CD3 Remission/Sequencing Data/Mouse scRNA/trust4_sc_results")

library(dplyr)
library(readr)

airr.files <- list.files(
  pattern = "_barcode_airr.tsv$",
  full.names = TRUE
)

length(airr.files)
airr.list <- lapply(airr.files, read_tsv)
names(airr.list) <- gsub("_barcode_airr.tsv","",basename(airr.files))


BiocManager::install("scRepertoire")
library(scRepertoire)
packageVersion("scRepertoire")


contig.list <- loadContigs(
  input = airr.list,
  format = "AIRR"
)
names(contig.list) <- names(airr.list)
colnames(contig.list[[1]])
head(contig.list[[1]])
table(contig.list[[1]]$chain, useNA = "ifany")

combined.TCR <- combineTCR(
  input.data = contig.list,
  samples = names(contig.list),
  ID = NULL,
  filterNonproductive = TRUE
)

combined.TCR <- lapply(
  combined.TCR,
  function(x) {
    x$barcode <- ifelse(
      grepl("-1$", x$barcode),
      x$barcode,
      paste0(x$barcode, "-1")
    )
    x
  }
)

head(combined.TCR[[1]]$barcode)

head(colnames(mouse_aCD3_rev), 20)
head(mouse_aCD3_rev@meta.data)

trust4_barcodes <- unlist(
  lapply(combined.TCR, function(x) x$barcode),
  use.names = FALSE
)

seurat_barcodes <- colnames(mouse_aCD3_rev)

sum(trust4_barcodes %in% seurat_barcodes)

sum(trust4_barcodes %in% seurat_barcodes) /
  length(trust4_barcodes)

mouse_aCD3_rev <- combineExpression(
  input.data = combined.TCR,
  sc.data = mouse_aCD3_rev,
  cloneCall = "aa",
  proportion = FALSE,
  cloneSize = c(
    Singleton = 1,
    Small = 5,
    Medium = 10,
    Large = 20,
    Hyperexpanded = 1000
  )
)
grep(
  "CT|clone|clonal",
  colnames(mouse_aCD3_rev@meta.data),
  value = TRUE
)
table(
  mouse_aCD3_rev$cloneSize,
  useNA = "ifany"
)

library(dplyr)

meta <- mouse_aCD3_rev@meta.data %>%
  mutate(cell = rownames(.))

clone_counts <- meta %>%
  filter(
    !is.na(CTaa),
    CTaa != ""
  ) %>%
  count(
    MouseID,
    CTaa,
    name = "clone_size_mouse"
  ) %>%
  mutate(
    clone_key = paste(MouseID, CTaa, sep = "___")
  )

meta <- meta %>%
  mutate(
    clone_key = ifelse(
      !is.na(CTaa) & CTaa != "",
      paste(MouseID, CTaa, sep = "___"),
      NA_character_
    )
  )

idx <- match(meta$clone_key, clone_counts$clone_key)

mouse_aCD3_rev$clone_size_mouse <-
  clone_counts$clone_size_mouse[idx]

mouse_aCD3_rev$clone_size_mouse[
  is.na(mouse_aCD3_rev$clone_size_mouse)
] <- 0

mouse_aCD3_rev$clone_category_mouse <- cut(
  mouse_aCD3_rev$clone_size_mouse,
  breaks = c(-Inf, 0, 1, 5, 10, 20, Inf),
  labels = c(
    "No TCR",
    "Singleton",
    "2–5 cells",
    "6–10 cells",
    "11–20 cells",
    ">20 cells"
  )
)

DimPlot(
  mouse_aCD3_rev,
  reduction = "umap",
  group.by = "clone_category_mouse",
  shuffle = TRUE
)

mouse_aCD3_rev$log_clone_size_mouse <-
  log10(mouse_aCD3_rev$clone_size_mouse + 1)

FeaturePlot(
  mouse_aCD3_rev,
  features = "log_clone_size_mouse",
  reduction = "umap",
  order = TRUE
)
DimPlot(
  mouse_aCD3_rev
)

library(Seurat)
library(dplyr)
library(ggplot2)
Idents(mouse_aCD3_rev)
## Extract UMAP coordinates
umap_df <- as.data.frame(
  Embeddings(mouse_aCD3_rev, reduction = "umap")
)

umap_df$cell <- rownames(umap_df)

## Add metadata
plot_df <- umap_df %>%
  left_join(
    mouse_aCD3_rev@meta.data %>%
      mutate(cell = rownames(.)),
    by = "cell"
  )

clonotype_umap <- ggplot(
  plot_df,
  aes(
    x = umap_1,
    y = umap_2
  )
) +
  geom_point(
    aes(
      color = factor(seurat_clusters),
      size = clone_size_mouse
    ),
    alpha = 0.75,
    stroke = 0
  ) +
  facet_wrap(
    ~ Group
  ) +
  scale_size_continuous(
    name = "Clonotype size",
    range = c(0.3, 5),
    breaks = c(1, 5, 10, 20, 50)
  ) +
  labs(
    color = "T-cell cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold")
  )

clonotype_umap
library(dplyr)
library(ggplot2)

## Check the available cell-type labels
library(Seurat)
library(dplyr)
library(ggplot2)

## Get barcodes classified as T cells
tcell_barcodes <- WhichCells(
  mouse_aCD3_rev,
  idents = "T cells"
)

## Keep only T cells in the plotting dataframe
tcell_df <- plot_df %>%
  filter(cell %in% tcell_barcodes)

nrow(tcell_df)
table(tcell_df$Group)

clonotype_umap <- ggplot(
  tcell_df,
  aes(
    x = umap_1,
    y = umap_2
  )
) +
  geom_point(
    aes(
      color = factor(seurat_clusters),
      size = log_clone_size_mouse
    ),
    alpha = 0.75,
    stroke = 0
  ) +
  facet_wrap(
    ~ Group
  ) +
  scale_size_continuous(
    name = "Clonotype size",
    range = c(0.4, 5),
    breaks = log10(c(1, 2, 5, 10, 20, 50) + 1),
    labels = c("1", "2", "5", "10", "20", "50")
  ) +
  labs(
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold")
  )

clonotype_umap
