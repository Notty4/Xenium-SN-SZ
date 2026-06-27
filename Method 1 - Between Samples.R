# =============================================================================
# Method 2: De Novo Clustering & Within-Sample Analysis
# =============================================================================
# Project:  Astrocytic Reprogramming in Schizophrenia - Substantia Nigra
# Author:   Harrison Nott
# Degree:   Master of Brain and Mind Science, University of Sydney (2025)
# Tissue:   Post-mortem human midbrain (substantia nigra pars compacta)
# Method:   Xenium spatial transcriptomics (10x Genomics)
#
# Description:
#   This script loads raw Xenium segmentation output for two tissue samples
#   (SN1: healthy control, SN4: schizophrenia), constructs Seurat objects
#   with spatial coordinates, merges them, and runs the full clustering
#   pipeline (SCTransform, PCA, UMAP, Leiden clustering). Clusters are
#   manually annotated by marker gene expression into 9 cell types.
#   Differential gene expression (SCZ vs Control) is then run across the
#   full merged dataset.
#
# Input files (place in ./data/ directory):
#   - SN1/segmentation.csv                  : Xenium transcript assignments (Control)
#   - SN1/segmentation_cell_stats.csv       : Cell-level spatial stats (Control)
#   - SN4/baysor_output_segmentation.csv    : Baysor transcript assignments (SCZ)
#   - SN4/baysor_output_cell_stats.csv      : Cell-level spatial stats (SCZ)
#
# Note: Raw data (tissue samples) are not publicly available as they
#   originate from the NSW Brain Tissue Resource Centre and are governed
#   by institutional ethics approval. This script is provided for
#   methodological transparency.
#
# Dependencies: Seurat, data.table, dplyr, tidyr, ggplot2, patchwork,
#               tibble, EnhancedVolcano, pheatmap
# =============================================================================

# =========================================================================
# 0. SETUP: LOAD LIBRARIES
# =========================================================================
library(Seurat)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(EnhancedVolcano)
library(pheatmap)

# =========================================================================
# 1. LOAD AND PRE-PROCESS SPATIAL SAMPLES
# =========================================================================
# This function loads the raw data, creates a Seurat object, and adds the spatial info
process_dataset <- function(dataset_name, seg_path, scs_path, condition) {
  seg <- fread(seg_path)
  cat(paste0("--- Processing Sample: ", dataset_name, " ---\n"))
  
  seg <- seg %>% filter(is_noise == FALSE, assignment_confidence >= 0.9)
  
  gene_counts <- seg %>% 
    dplyr::count(gene, cell) %>%
    pivot_wider(names_from = cell, values_from = n, values_fill = list(n = 0))
  
  gene_matrix <- as.matrix(gene_counts[,-1])
  rownames(gene_matrix) <- gene_counts$gene
  
  seu_obj <- CreateSeuratObject(counts = gene_matrix, project = dataset_name)
  
  scs <- fread(scs_path)
  scs <- scs |> filter(cell %in% colnames(seu_obj))
  rownames(scs) <- scs$cell
  seu_obj <- AddMetaData(seu_obj, metadata = scs)
  
  seu_obj <- subset(seu_obj, subset = area >= 100 & nCount_RNA >= 20 & nFeature_RNA >= 20)
  
  # Add spatial coordinates to the object
  centroids <- data.frame(
    x = seu_obj@meta.data[["x"]],
    y = seu_obj@meta.data[["y"]],
    cell = colnames(seu_obj)
  )
  cents <- CreateCentroids(centroids)
  coords <- CreateFOV(coords = list("centroids" = cents), type = "centroids")
  seu_obj[[dataset_name]] <- coords # Adds the FOV to the @images slot
  
  seu_obj$condition <- condition
  
  cat(paste0("--- Finished processing sample: ", dataset_name, " ---\n\n"))
  return(seu_obj)
}

# --- Load your two spatial samples ---
# !!! IMPORTANT: Double-check these file paths are correct !!!
SN1c <- process_dataset(
  dataset_name = "SN1",
  seg_path = "data/SN1/segmentation.csv",
  scs_path = "data/SN1/segmentation_cell_stats.csv",
  condition = "Control"
)

SN4c <- process_dataset(
  dataset_name = "SN4",
  seg_path = "data/SN4/baysor_output_segmentation.csv",
  scs_path = "data/SN4/baysor_output_cell_stats.csv",
  condition = "Schizophrenia"
)

# =========================================================================
# 2. MERGE, PROCESS, AND CLUSTER
# =========================================================================

# --- Merge objects while preserving spatial data ---
# 1. Merge the assays and metadata
merged_obj <- merge(SN1c, y = SN4c, add.cell.ids = c("SN1", "SN4"), project = "ThesisProject")

# 2. Manually transfer the spatial information to the merged object
merged_obj@images$SN1 <- SN1c@images$SN1
merged_obj@images$SN4 <- SN4c@images$SN4

# --- Normalize and run the main analysis pipeline ---
merged_obj <- SCTransform(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:30)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.5)

# =========================================================================
# 3. ANNOTATE CLUSTERS
# =========================================================================

# --- Find marker genes for each cluster ---
all_markers <- FindAllMarkers(merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# --- Define cell type annotations based on markers ---
new_cluster_ids <- c("0" = "Oligodendrocytes",
                     "1" = "Astrocytes",
                     "2" = "Astrocytes (Reactive)",
                     "3" = "GABAergic Neurons",
                     "4" = "Dopaminergic Neurons",
                     "5" = "Endothelial Cells",
                     "6" = "Dopaminergic Neurons (DAT+)",
                     "7" = "Neurons",
                     "8" = "Dopaminergic Neurons (Calbindin+)")

# --- Rename the clusters in the Seurat object ---
merged_obj <- RenameIdents(merged_obj, new_cluster_ids)

# =========================================================================
# 4. ANALYZE AND VISUALIZE RESULTS
# =========================================================================

# --- Plot 1: Annotated UMAP ---
DimPlot(merged_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Annotated Cell Types in the Substantia Nigra")

# --- Plot 2: Spatial Maps of Cell Types ---
p_control <- SpatialDimPlot(merged_obj, images = "SN1", label = TRUE, label.size = 3) + 
  ggtitle("Control (SN1)")
p_schizophrenia <- SpatialDimPlot(merged_obj, images = "SN4", label = TRUE, label.size = 3) + 
  ggtitle("Schizophrenia (SN4)")
p_control + p_schizophrenia

# --- Analysis 1: Cell Proportions ---
counts_by_condition <- table(Idents(merged_obj), merged_obj$condition)
print(counts_by_condition)
chi_squared_result <- chisq.test(counts_by_condition)
print(chi_squared_result)
print("Chi-squared Residuals:")
print(chi_squared_result$residuals)

# --- Analysis 2: Differential Gene Expression (Schizophrenia vs. Control) ---
Idents(merged_obj) <- "condition"
merged_obj <- PrepSCTFindMarkers(merged_obj)
de_results <- FindMarkers(merged_obj, ident.1 = "Schizophrenia", ident.2 = "Control")

# --- Plot 3: Volcano Plot of DGE results ---
EnhancedVolcano(
  as.data.frame(de_results),
  lab = rownames(de_results),
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.5
)

# --- Save the final, fully processed object ---
saveRDS(merged_obj, file = "data/Final_Analyzed_Object.rds")

message("--- ANALYSIS COMPLETE ---")