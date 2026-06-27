# =============================================================================
# Method 1: Reference Atlas Mapping & Cell-Type-Specific DGE
# =============================================================================
# Project:  Astrocytic Reprogramming in Schizophrenia - Substantia Nigra
# Author:   Harrison Nott
# Degree:   Master of Brain and Mind Science, University of Sydney (2025)
# Tissue:   Post-mortem human midbrain (substantia nigra pars compacta)
# Method:   Xenium spatial transcriptomics (10x Genomics)
# Reference: Agarwal et al. (2020) midbrain atlas (GSE140231)
#
# Description:
#   This script performs reference atlas mapping of spatial transcriptomics
#   data onto the Agarwal midbrain reference atlas using Seurat's
#   FindTransferAnchors/MapQuery pipeline. It then runs cell-type-specific
#   differential gene expression (DGE) and GSEA across anatomical subregions
#   (SNcd, SNcv, VTA) in schizophrenia vs. control tissue.
#
# Input files (place in ./data/ directory):
#   - Method2_directSlices.rds       : Processed Seurat object from Method 2
#   - Atlas/integrated_reference_corrected.rds : Agarwal reference atlas
#   - Final_Mapped_Object.rds        : Atlas-mapped Seurat object
#   - Final_Analyzed_Object.rds      : Annotated merged Seurat object
#   - SN1/SN4 ROI subset .rds files  : Spatial subregion objects
#
# Note: Raw data (tissue samples) are not publicly available as they
#   originate from the NSW Brain Tissue Resource Centre and are governed
#   by institutional ethics approval. This script is provided for
#   methodological transparency.
#
# Dependencies: Seurat, ggplot2, dplyr, tibble, pheatmap, clusterProfiler,
#               org.Hs.eg.db, ggrepel, patchwork, readr
# =============================================================================

# =========================================================================
# 0. SETUP: LOAD LIBRARIES
# =========================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)


# =========================================================================
# 1. LOAD REQUIRED DATA
# =========================================================================
# This script requires two inputs:
# 1. Your fully analyzed object from Method 2.
# 2. The corrected reference atlas we built earlier.

# --- Load your analyzed object from Method 2 ---
# It contains your own clusters and annotations.
query_obj <- readRDS("data/Method2_directSlices.rds")

# --- Load the reference atlas ---
reference <- readRDS("data/Atlas/integrated_reference_corrected.rds")


# =========================================================================
# 2. PREPARE OBJECTS FOR MAPPING
# =========================================================================
# Set the default assay for both objects to SCT for a consistent workflow
DefaultAssay(query_obj) <- "SCT"
DefaultAssay(reference) <- "SCT"

# Defensive check: Ensure the reference has a projectable UMAP model and clusters
# This prevents errors in the MapQuery step.
if (!"umap" %in% names(reference@reductions) || is.null(reference@reductions$umap@misc$model)) {
  reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
}
if (!"seurat_clusters" %in% colnames(reference@meta.data)) {
  reference <- FindNeighbors(reference, dims = 1:30)
  reference <- FindClusters(reference, resolution = 0.5)
}

# =========================================================================
# 3. PERFORM REFERENCE MAPPING
# =========================================================================
# Find the set of genes that are common to both your data and the reference
common_genes <- intersect(rownames(reference), rownames(query_obj))

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = query_obj,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30,
  features = common_genes
)

# Map your data onto the reference atlas
# This adds the predicted cell types from the reference to your object's metadata
mapped_obj <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query_obj,
  refdata = list(predicted_celltype = "seurat_clusters"), # Transfer the reference's cluster IDs
  reference.reduction = "pca",
  reduction.model = "umap"
)


library(Seurat)
library(dplyr)

# --- 1. Load your reference atlas ---
reference <- readRDS("data/Atlas/integrated_reference_corrected.rds")

# --- 2. Define the new, detailed cluster names ---
# This is a named vector based on the marker gene analysis.
# The number on the left is the original cluster ID, the name on the right is the new annotation.
new_reference_ids <- c(
  "0" = "Astrocytes (Type 1)",
  "1" = "Neurons (Generic)",
  "2" = "Excitatory Neurons",
  "3" = "Oligodendrocytes (Type 1)",
  "4" = "Dopaminergic Neurons (Subtype 1)",
  "5" = "GABAergic Neurons (PVALB+)",
  "6" = "Astrocytes (Type 2)",
  "7" = "GABAergic Neurons (VIP+)",
  "8" = "Microglia",
  "9" = "Excitatory Neurons",
  "10" = "Dopaminergic Neurons (Subtype 2)",
  "11" = "Oligodendrocyte Precursors (OPCs)",
  "12" = "Oligodendrocytes (Type 2)",
  "13" = "GABAergic Neurons",
  "14" = "GABAergic Neurons (SST+)",
  "15" = "Dopaminergic Neurons (Subtype 3)",
  "16" = "Oligodendrocytes (Type 3)",
  "17" = "Dopaminergic Neurons (Subtype 4)",
  "18" = "Oligodendrocytes (Type 4)",
  "19" = "Dopaminergic Neurons (Subtype 5)",
  "20" = "Oligodendrocytes (Type 5)",
  "21" = "GABAergic Neurons (PVALB+)",
  "22" = "GABAergic Neurons (LAMP5+)",
  "23" = "Endothelial Cells"
)

# --- 3. Apply the new annotations ---
# First, set the active identity to the cluster numbers
Idents(reference) <- "seurat_clusters"

# Now, rename the identities using your list
reference <- RenameIdents(reference, new_reference_ids)

# --- 4. Visualize the fully annotated reference atlas ---
DimPlot(reference, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Annotated Cell Types in the Reference Atlas")

# --- 5. (Optional but Recommended) Save the annotated object ---
# This saves your work so you don't have to re-annotate every time.
saveRDS(reference, file = "data/Atlas/integrated_reference_ANNOTATED.rds")

# =========================================================================
# 4. EXTRACT KEY "METHOD 1" STATISTICS AND VISUALIZATIONS
# =========================================================================

# --- Result 1: UMAP of Predicted Cell Identities ---
# This plot shows your cells colored by the labels transferred from the reference atlas.
DimPlot(
  mapped_obj,
  reduction = "umap",
  group.by = "predicted.predicted_celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Cell Identities Predicted by Reference Mapping")

# --- Result 2: Mapping Confidence Score ---
# This histogram shows how confident the algorithm was in its predictions.
ggplot(mapped_obj@meta.data, aes(x = predicted.predicted_celltype.score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribution of Mapping Confidence Scores",
    x = "Prediction Score (0 to 1)",
    y = "Number of Cells"
  )

# --- Result 3: Focused, Cell-Type-Specific DGE (Your Key Finding) ---
# This analysis finds gene changes WITHIN a specific population of interest.
# First, you need to know which cluster ID from the REFERENCE corresponds to DA neurons.
# Let's assume it's cluster '3'. You can check your reference object to be sure.
# !!! IMPORTANT: Update the cluster number '3' if needed !!!
subset_obj <- subset(
  mapped_obj,
  subset = (predicted.predicted_celltype == 3) # Filter for cells mapped to reference cluster 3
)

Idents(subset_obj) <- "condition"
subset_obj <- PrepSCTFindMarkers(subset_obj)

de_celltype_specific <- FindMarkers(
  subset_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control"
)

# Print the top 10 genes from this focused analysis
print("--- Top Dysregulated Genes in Predicted Dopaminergic Neurons ---")
print(head(de_celltype_specific, 10))


library(Seurat)
library(dplyr)
library(tibble)

# --- 1. Load your reference atlas ---
reference <- readRDS("data/Atlas/integrated_reference_corrected.rds")

# --- 2. Set the cluster identities ---
Idents(reference) <- "seurat_clusters"

# --- 3. !! THE DEFINITIVE FIX !! Prepare the integrated object for DE ---
# This command is required before running FindMarkers on an integrated SCT object.
message("--- Preparing the integrated reference for marker analysis... ---")
reference <- PrepSCTFindMarkers(reference)

# --- 4. Manually find markers for each cluster using a loop ---
message("--- Starting manual marker finding loop... ---")

# Get the list of all unique cluster IDs
cluster_ids <- levels(Idents(reference))
all_markers_list <- list()

# Loop through each cluster ID
for (cluster_id in cluster_ids) {
  
  message(paste("Finding markers for Cluster", cluster_id))
  
  # This will now work correctly
  cluster_markers <- FindMarkers(
    reference, 
    ident.1 = cluster_id, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )
  
  cluster_markers <- cluster_markers %>%
    rownames_to_column("gene") %>%
    mutate(cluster = cluster_id)
  
  all_markers_list[[cluster_id]] <- cluster_markers
}

# --- 5. Combine all results into a single data frame ---
reference_markers <- bind_rows(all_markers_list)
message("--- Manual marker finding complete. ---")

# --- 6. Get the Top 10 Markers ---
top10_ref_markers <- reference_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# Print the final list
print(top10_ref_markers, n = Inf)





library(Seurat)
library(dplyr)
library(tibble)

# Load your mapped object
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")

# --- Let's quickly confirm what's in the prediction column ---
# This will print the list of unique cluster IDs (they will be numbers)
print("Unique predicted cluster IDs found in your data:")
print(unique(mapped_obj$predicted.predicted_celltype))

# --- Perform DGE within a specific predicted cell type ---
# We use the CLUSTER NUMBER (14) instead of the text label.
# !!! You can change '14' to any other reference cluster number you want to analyze !!!
subset_obj <- subset(
  mapped_obj,
  subset = (predicted.predicted_celltype == 14) # Use the number, not the name
)

Idents(subset_obj) <- "condition"
subset_obj <- PrepSCTFindMarkers(subset_obj)

de_celltype_specific <- FindMarkers(
  subset_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control"
)

# --- View the results ---
# This table will show you the gene changes specific to only that cell type.
print("--- Top Dysregulated Genes in Predicted Cluster 14 (GABAergic Neurons, SST+) ---")
print(head(de_celltype_specific, 10))

library(dplyr)
library(tibble)

# First, convert the rownames to a 'gene' column
de_celltype_final <- de_celltype_specific %>%
  rownames_to_column("gene")

# Now, filter out the placeholder names
de_celltype_final <- de_celltype_final %>%
  filter(!grepl("DeprecatedCodeword", gene)) %>%
  arrange(desc(abs(avg_log2FC))) # Sort by the magnitude of the change

# Print the final, clean, and sorted table
print("--- Final Cleaned DGE Table for SST+ Neurons ---")
print(head(de_celltype_final, 10))
# =========================================================================
# 5. SAVE THE FINAL MAPPED OBJECT
# =========================================================================
saveRDS(mapped_obj, file = "data/Final_Mapped_Object.rds")

message("--- METHOD 1 ANALYSIS IS COMPLETE ---")



#POST


library(Seurat)
library(dplyr)
library(pheatmap)

# ===================================================================
# 1. LOAD YOUR FINAL MAPPED OBJECT
# ===================================================================
# This object contains both your own annotations and the predicted ones.
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")

# ===================================================================
# 2. ENSURE YOUR OWN ANNOTATIONS ARE PRESENT
# ===================================================================
# Let's quickly re-apply your own annotations to be safe.
new_cluster_ids <- c("0" = "Oligodendrocytes",
                     "1" = "Astrocytes",
                     "2" = "Astrocytes (Reactive)",
                     "3" = "GABAergic Neurons",
                     "4" = "Dopaminergic Neurons",
                     "5" = "Endothelial Cells",
                     "6" = "Dopaminergic Neurons (DAT+)",
                     "7" = "Neurons",
                     "8" = "Dopaminergic Neurons (Calbindin+)")

# First, set the active identity to the original cluster numbers
Idents(mapped_obj) <- "seurat_clusters"
# Now, rename them with your annotations
mapped_obj <- RenameIdents(mapped_obj, new_cluster_ids)

# ===================================================================
# 3. CREATE A CONTINGENCY TABLE OF CLUSTER ASSIGNMENTS
# ===================================================================
# This table compares your annotations (Method 2) against the reference's predictions (Method 1)
# Rows: Your annotations (from Method 2)
# Columns: The reference's predicted cluster numbers (from Method 1)
contingency_table <- table(
  `Your Annotation` = Idents(mapped_obj),
  `Reference Prediction` = mapped_obj$predicted.predicted_celltype
)

# Print the raw counts to see the numbers
print("--- Raw Cell Counts: Your Clusters vs. Reference Predictions ---")
print(contingency_table)

# ===================================================================
# 4. VISUALIZE THE OVERLAP WITH A HEATMAP
# ===================================================================
# For a fair comparison, we'll convert the raw counts to percentages.
# This answers the question: "Of the cells I called 'Astrocytes', what percentage
# did the reference map to each of its own clusters?"
percentage_table <- prop.table(contingency_table, margin = 1) * 100

# Generate the heatmap
pheatmap(
  percentage_table,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = FALSE, # Keep your annotations in their original order
  cluster_cols = FALSE, # Keep the reference clusters in numerical order
  main = "Overlap Between Your Clusters and Reference Predictions (%)",
  fontsize = 8
)







library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# =========================================================================
# FIGURE FOR SLIDE 2: Log2 Fold Change in CELL TYPE PROPORTIONS
# =========================================================================
message("--- Generating Figure for Slide 2: Cell Proportions ---")

# --- 1. Load your main analyzed object from Method 2 ---
# This object contains your 9 annotated cell types.
merged_obj <- readRDS("data/Final_Analyzed_Object.rds")

# --- 2. Create the cell count table ---
cell_counts_table <- table(Idents(merged_obj), merged_obj$condition)
cell_counts_df <- as.data.frame(cell_counts_table)
colnames(cell_counts_df) <- c("CellType", "Condition", "Count")

# --- 3. Calculate Proportions and Log2 Fold Change ---
proportion_data <- cell_counts_df %>%
  group_by(Condition) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  select(CellType, Condition, Proportion) %>%
  pivot_wider(names_from = Condition, values_from = Proportion) %>%
  mutate(
    # Add a small number (pseudocount) to avoid division by zero
    Control = Control + 1e-9,
    Schizophrenia = Schizophrenia + 1e-9,
    log2FoldChange = log2(Schizophrenia / Control)
  )

# --- 4. Create the Plot ---
cell_proportion_plot <- ggplot(proportion_data, aes(x = reorder(CellType, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flips the axes to make it a horizontal bar chart
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = "none") +
  labs(
    title = "Log2 Fold Change in Cell Type Proportions (Schizophrenia vs. Control)",
    x = "Cell Type",
    y = "Log2 Fold Change"
  ) +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_minimal()

# --- 5. Display the Plot ---
print(cell_proportion_plot)


# =========================================================================
# FIGURE FOR SLIDE 3: Log2 Fold Change of KEY GENES
# =========================================================================
message("--- Generating Figure for Slide 3: Key Genes ---")

# --- 1. Load your final MAPPED object from Method 1 ---
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")

# --- 2. Perform the focused DGE within a specific cell type ---
# This is where your most biologically relevant gene list came from.
# We will use the SST+ GABAergic Neuron cluster (Reference Cluster 14) as the example.
# !!! IMPORTANT: Update the cluster number if you want to analyze a different population !!!
subset_obj <- subset(
  mapped_obj,
  subset = (predicted.predicted_celltype == 14) # Filter for cells mapped to reference cluster 14
)

Idents(subset_obj) <- "condition"
subset_obj <- PrepSCTFindMarkers(subset_obj)

de_celltype_specific <- FindMarkers(
  subset_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control"
)

# --- 3. Create the final, curated list of top genes ---
# This cleans the table and selects the top genes by fold change.
final_genes_table <- de_celltype_specific %>%
  rownames_to_column("gene") %>%
  filter(!grepl("DeprecatedCodeword", gene) & !is.na(gene)) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(10) # Select the top 10 most changed genes

# --- 4. Create the Plot ---
key_genes_plot <- ggplot(final_genes_table, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), guide = "none") +
  labs(
    title = "Top Dysregulated Genes in SST+ Neurons (Schizophrenia vs. Control)",
    x = "Gene",
    y = "Log2 Fold Change"
  ) +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_minimal()

# --- 5. Display the Plot ---
print(key_genes_plot)






# =========================================================================
# FAILSAFE DGE ON CLUSTERS 6 & 10 (USING CLASSIC WORKFLOW)
# =========================================================================
library(Seurat)
library(dplyr)
library(tibble)

# --- 1. Load your final MAPPED object from Method 1 ---
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")


# =========================================================================
# ANALYSIS 1: CELLS MAPPED TO REFERENCE CLUSTER 6 (ASTROCYTE-LIKE)
# =========================================================================

# --- Subset for cells mapped to Cluster 6 ---
subset_c6 <- subset(
  mapped_obj,
  subset = (predicted.predicted_celltype == 6)
)

# --- !! THE FAILSAFE FIX !! ---
# Switch to the raw RNA assay and run the classic workflow
DefaultAssay(subset_c6) <- "RNA"
subset_c6 <- NormalizeData(subset_c6)
subset_c6 <- FindVariableFeatures(subset_c6)
subset_c6 <- ScaleData(subset_c6)

# --- Now, FindMarkers will work correctly on the RNA assay ---
Idents(subset_c6) <- "condition"
de_cluster_6 <- FindMarkers(
  subset_c6,
  ident.1 = "Schizophrenia",
  ident.2 = "Control"
)

# --- View the results for Cluster 6 ---
print("--- Top Dysregulated Genes in Cells Mapped to Reference Cluster 6 (Astrocytes) ---")
print(head(de_cluster_6, 10))


# =========================================================================
# ANALYSIS 2: CELLS MAPPED TO REFERENCE CLUSTER 10 (NEURON-LIKE)
# =========================================================================

# --- Subset for cells mapped to Cluster 10 ---
subset_c10 <- subset(
  mapped_obj,
  subset = (predicted.predicted_celltype == 10)
)

# --- !! THE FAILSAFE FIX !! ---
# Switch to the raw RNA assay and run the classic workflow
DefaultAssay(subset_c10) <- "RNA"
subset_c10 <- NormalizeData(subset_c10)
subset_c10 <- FindVariableFeatures(subset_c10)
subset_c10 <- ScaleData(subset_c10)

# --- Now, FindMarkers will work correctly on the RNA assay ---
Idents(subset_c10) <- "condition"
de_cluster_10 <- FindMarkers(
  subset_c10,
  ident.1 = "Schizophrenia",
  ident.2 = "Control"
)

# --- View the results for Cluster 10 ---
print("--- Top Dysregulated Genes in Cells Mapped to Reference Cluster 10 (Neurons) ---")
print(head(de_cluster_10, 10))



# This assumes 'merged_obj' is your final annotated object from Method 2

# --- Final failsafe script to highlight a specific cell type on the spatial map ---
library(ggplot2)
library(dplyr)

# 1. Extract all the necessary data into a clean data frame
plot_data <- data.frame(
  x = merged_obj@meta.data$x,
  y = merged_obj@meta.data$y,
  cell_type = Idents(merged_obj),
  condition = merged_obj@meta.data$condition
)

# 2. Create a new column to specify which cells to highlight
# We will set the color to "grey" for all cells, then change it for our target.
plot_data$highlight_color <- "grey80" # Default color for non-target cells
plot_data$highlight_color[plot_data$cell_type == "Astrocytes (Reactive)"] <- "yellow"

# 3. Create a size column to make the highlighted cells slightly larger
plot_data$highlight_size <- 0.5 # Default size
plot_data$highlight_size[plot_data$cell_type == "Astrocytes (Reactive)"] <- 1.0

# 4. Create the final, styled plot
# We sort the data frame so the highlighted cells are plotted on top
ggplot(plot_data %>% arrange(highlight_color == "grey80"), aes(x = x, y = y, color = highlight_color, size = highlight_size)) +
  geom_point() +
  # Use facet_wrap to create two separate panels for Control and Schizophrenia
  facet_wrap(~condition, ncol = 2) +
  # Manually set the colors
  scale_color_identity() +
  # Manually set the sizes
  scale_size_identity() +
  coord_fixed() + # Ensures the aspect ratio of the tissue is correct
  labs(title = "Spatial Location of Reactive Astrocytes") +
  # This detailed theme section replicates the professional, dark aesthetic
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black", color = "black"),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey40", size = 0.25),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(color = "white", hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none" # Remove the legend as the title is self-explanatory
  )






# Extract coordinates and TH expression data
plot_data <- data.frame(
  x = merged_obj@meta.data$x,
  y = merged_obj@meta.data$y,
  TH = GetAssayData(merged_obj, assay = "SCT")["TH",]
)

# Create the plot with only the vertical flip
ggplot(plot_data, aes(x = x, y = y, color = TH)) +
  geom_point(size = 1.2) +
  scale_color_gradient(low = "grey", high = "yellow") +
  coord_fixed() +
  theme_void() +
  ggtitle("Correctly Oriented TH Expression") +
  
  # --- ONLY THE VERTICAL FLIP REMAINS ---
  scale_y_reverse()






# --- Coordinates for the ENTIRE SNc ---
roi_coords_sn1_snc <- data.frame(
  x = c(6090, 5684, 5441, 4441, 3468, 2833, 3765, 4535, 5198, 6523, 7667, 8392, 9069, 9368, 8516, 7603, 7101, 6780),
  y = c(4748, 4343, 3897, 3532, 2640, 1518, 1599, 1856, 2329, 2517, 2452, 2546, 3232, 3760, 4433, 4433, 4678, 4356)
)

# --- Coordinates for the DIVIDING LINE (Updated) ---
dividing_line_coords <- data.frame(
  x = c(3605, 3971, 4767, 5436, 5943, 6760, 7189, 7661, 8034, 8630), # <-- Added 8569
  y = c(2629, 2622, 2932, 3291, 3509, 3692, 3671, 3411, 3045, 2800)  # <-- Added 2833
)

# --- Coordinates for the VTA ---
roi_coords_sn1_vta <- data.frame(
  x = c(10562, 9992, 9605, 9189, 8619, 8745, 9288, 9978, 10245, 10548, 10844, 10851),
  y = c(3023, 3185, 3129, 2967, 2566, 2340, 2002, 1756, 1869, 2108, 2446, 2714)
)

# =========================================================================
# 4. CREATE SNcd AND SNcv POLYGONS
# =========================================================================
# --- SNcd Polygon ---
sncd_polygon <- rbind(
  roi_coords_sn1_snc[6:13, ],
  dividing_line_coords[nrow(dividing_line_coords):1, ]
)
# --- SNcv Polygon ---
sncv_polygon <- rbind(
  roi_coords_sn1_snc[1:5, ],
  dividing_line_coords,
  roi_coords_sn1_snc[14:18, ]
)

# =========================================================================
# 5. EXTRACT SN1 CELL COORDINATES
# =========================================================================
sn1_cell_coords <- subset(merged_obj@meta.data, subset = orig.ident == "SN1")[, c("x", "y")]

# =========================================================================
# 6. ISOLATE CELLS FOR EACH REGION (SNcd, SNcv, VTA)
# =========================================================================

# --- SNcd ---
points_in_sncd <- point.in.polygon(point.x = sn1_cell_coords$x, point.y = sn1_cell_coords$y, pol.x = sncd_polygon$x, pol.y = sncd_polygon$y)
cells_in_sn1_sncd <- rownames(sn1_cell_coords)[points_in_sncd > 0]
sn1_sncd_obj <- subset(merged_obj, cells = cells_in_sn1_sncd)
message(paste("Isolated", length(cells_in_sn1_sncd), "cells in SNcd for SN1."))
saveRDS(sn1_sncd_obj, file = "data/SN1_SNcd_Only_Object.rds")

# --- SNcv ---
points_in_sncv <- point.in.polygon(point.x = sn1_cell_coords$x, point.y = sn1_cell_coords$y, pol.x = sncv_polygon$x, pol.y = sncv_polygon$y)
cells_in_sn1_sncv <- rownames(sn1_cell_coords)[points_in_sncv > 0]
sn1_sncv_obj <- subset(merged_obj, cells = cells_in_sn1_sncv)
message(paste("Isolated", length(cells_in_sn1_sncv), "cells in SNcv for SN1."))
saveRDS(sn1_sncv_obj, file = "data/SN1_SNcv_Only_Object.rds")

# --- VTA ---
points_in_vta <- point.in.polygon(point.x = sn1_cell_coords$x, point.y = sn1_cell_coords$y, pol.x = roi_coords_sn1_vta$x, pol.y = roi_coords_sn1_vta$y)
cells_in_sn1_vta <- rownames(sn1_cell_coords)[points_in_vta > 0]
sn1_vta_obj <- subset(merged_obj, cells = cells_in_sn1_vta)
message(paste("Isolated", length(cells_in_sn1_vta), "cells in VTA for SN1."))
saveRDS(sn1_vta_obj, file = "data/SN1_VTA_Only_Object.rds")

# =========================================================================
# 7. (Optional) VISUALIZE ALL SELECTIONS
# =========================================================================
plot_data <- sn1_cell_coords
plot_data$subregion <- "Outside"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn1_sncd] <- "SNcd"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn1_sncv] <- "SNcv"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn1_vta] <- "VTA"

plot_data$subregion <- factor(plot_data$subregion, levels = c("Outside", "VTA", "SNcv", "SNcd"))

ggplot(plot_data, aes(x = x, y = y, color = subregion, size = subregion == "Outside")) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("SNcd" = "magenta", "SNcv" = "green", "VTA" = "cyan", "Outside" = "grey40")) +
  scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 1.5)) +
  coord_fixed() +
  ggtitle("SN1 SNcd, SNcv, and VTA Cell Selections") +
  scale_y_reverse() +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(color = "white", hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )


##SN4##

roi_coords_sn4_snc <- data.frame(
  x = c(3690, 3902, 4051, 7268, 9387, 9371, 9285, 9026, 8571, 6256, 4805),
  y = c(2850, 1720, 943, 1147, 896, 1971, 2685, 2991, 3235, 3282, 3133)
)

# --- Coordinates for the DIVIDING LINE (SN4) ---
dividing_line_coords_sn4 <- data.frame(
  x = c(3902, 6713, 8384, 9371),
  y = c(1720, 2030, 2077, 1971)
)

# --- Coordinates for the VTA (SN4) ---
roi_coords_sn4_vta <- data.frame(
  x = c(9662, 9301, 9167, 9202, 9406, 9851, 10048, 10478, 10425, 10070),
  y = c(2387, 2168, 1924, 1915, 1598, 1394, 1447, 1628, 1998, 2459)
)

# =========================================================================
# 4. CREATE SNcd AND SNcv POLYGONS (SN4)
# =========================================================================
# --- SNcd Polygon (SN4) ---
# Assumes points 2-6 form the top curve of the SNc outline
sncd_polygon_sn4 <- rbind(
  roi_coords_sn4_snc[2:6, ],
  dividing_line_coords_sn4[nrow(dividing_line_coords_sn4):1, ] # Reversed dividing line
)
# --- SNcv Polygon (SN4) ---
# Assumes points 1 and 7-11 form the bottom curve of the SNc outline
sncv_polygon_sn4 <- rbind(
  roi_coords_sn4_snc[1, ], # Point 1
  dividing_line_coords_sn4, # Original dividing line
  roi_coords_sn4_snc[7:11, ] # Points 7-11
)

# =========================================================================
# 5. EXTRACT SN4 CELL COORDINATES
# =========================================================================
sn4_cell_coords <- subset(merged_obj@meta.data, subset = orig.ident == "SN4")[, c("x", "y")]

# =========================================================================
# 6. ISOLATE CELLS FOR EACH REGION (SN4)
# =========================================================================

# --- SNcd (SN4) ---
points_in_sncd_sn4 <- point.in.polygon(point.x = sn4_cell_coords$x, point.y = sn4_cell_coords$y, pol.x = sncd_polygon_sn4$x, pol.y = sncd_polygon_sn4$y)
cells_in_sn4_sncd <- rownames(sn4_cell_coords)[points_in_sncd_sn4 > 0]
sn4_sncd_obj <- subset(merged_obj, cells = cells_in_sn4_sncd)
message(paste("Isolated", length(cells_in_sn4_sncd), "cells in SNcd for SN4."))
saveRDS(sn4_sncd_obj, file = "data/SN4_SNcd_Only_Object.rds")

# --- SNcv (SN4) ---
points_in_sncv_sn4 <- point.in.polygon(point.x = sn4_cell_coords$x, point.y = sn4_cell_coords$y, pol.x = sncv_polygon_sn4$x, pol.y = sncv_polygon_sn4$y)
cells_in_sn4_sncv <- rownames(sn4_cell_coords)[points_in_sncv_sn4 > 0]
sn4_sncv_obj <- subset(merged_obj, cells = cells_in_sn4_sncv)
message(paste("Isolated", length(cells_in_sn4_sncv), "cells in SNcv for SN4."))
saveRDS(sn4_sncv_obj, file = "data/SN4_SNcv_Only_Object.rds")

# --- VTA (SN4) ---
points_in_vta_sn4 <- point.in.polygon(point.x = sn4_cell_coords$x, point.y = sn4_cell_coords$y, pol.x = roi_coords_sn4_vta$x, pol.y = roi_coords_sn4_vta$y)
cells_in_sn4_vta <- rownames(sn4_cell_coords)[points_in_vta_sn4 > 0]
sn4_vta_obj <- subset(merged_obj, cells = cells_in_sn4_vta)
message(paste("Isolated", length(cells_in_sn4_vta), "cells in VTA for SN4."))
saveRDS(sn4_vta_obj, file = "data/SN4_VTA_Only_Object.rds")

# =========================================================================
# 7. (Optional) VISUALIZE ALL SELECTIONS (SN4)
# =========================================================================
plot_data <- sn4_cell_coords
plot_data$subregion <- "Outside"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn4_sncd] <- "SNcd"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn4_sncv] <- "SNcv"
plot_data$subregion[rownames(plot_data) %in% cells_in_sn4_vta] <- "VTA"

plot_data$subregion <- factor(plot_data$subregion, levels = c("Outside", "VTA", "SNcv", "SNcd"))

ggplot(plot_data, aes(x = x, y = y, color = subregion, size = subregion == "Outside")) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("SNcd" = "magenta", "SNcv" = "green", "VTA" = "cyan", "Outside" = "grey40")) +
  scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 1.5)) +
  coord_fixed() +
  ggtitle("SN4 SNcd, SNcv, and VTA Cell Selections") +
  scale_y_reverse() + # Keep if needed for orientation
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black", color = "black"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(color = "white", hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )





###ROI ANALYSIS###

# --- 1. Load Required Libraries ---
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr) # For data manipulation

# =========================================================================
# 2. LOAD ALL SUBSET OBJECTS
# =========================================================================
# We will load the six separate .rds files you created.
# !!! Update these paths if you saved them in a different location !!!

sn1_vta_obj <- readRDS("data/SN1_VTA_Only_Object.rds")
sn1_sncd_obj <- readRDS("data/SN1_SNcd_Only_Object.rds")
sn1_sncv_obj <- readRDS("data/SN1_SNcv_Only_Object.rds")

sn4_vta_obj <- readRDS("data/SN4_VTA_Only_Object.rds")
sn4_sncd_obj <- readRDS("data/SN4_SNcd_Only_Object.rds")
sn4_sncv_obj <- readRDS("data/SN4_SNcv_Only_Object.rds")

# =========================================================================
# 3. ADD REGION-SPECIFIC LABELS
# =========================================================================
# We add new metadata columns to label each object BEFORE merging.

sn1_vta_obj$subregion <- "VTA"
sn1_sncd_obj$subregion <- "SNcd"
sn1_sncv_obj$subregion <- "SNcv"

sn4_vta_obj$subregion <- "VTA"
sn4_sncd_obj$subregion <- "SNcd"
sn4_sncv_obj$subregion <- "SNcv"

# =========================================================================
# 4. MERGE INTO ONE FINAL ANALYSIS OBJECT
# =========================================================================
# This merges all six objects into a single, combined Seurat object.
analysis_obj <- merge(
  x = sn1_vta_obj,
  y = list(sn1_sncd_obj, sn1_sncv_obj, sn4_vta_obj, sn4_sncd_obj, sn4_sncv_obj)
)

# --- Create a combined metadata column for easy grouping ---
# This gives us labels like "Control_SNcd" and "Schizophrenia_VTA"
analysis_obj$region_condition <- paste(analysis_obj$condition, analysis_obj$subregion, sep = "_")

# Set the factor levels for correct plot ordering
analysis_obj$region_condition <- factor(analysis_obj$region_condition, levels = c(
  "Control_SNcd", "Schizophrenia_SNcd",
  "Control_SNcv", "Schizophrenia_SNcv",
  "Control_VTA", "Schizophrenia_VTA"
))

# =========================================================================
# 5. START ANALYSIS: COMPARE "DIFFERENT CLUSTERS" (CELL PROPORTIONS)
# =========================================================================

# --- !! ADD THIS FIX !! ---
# We must re-set the active Idents to the cell type clusters from Method 2.
# The original cluster IDs (0-8) are stored in the 'seurat_clusters' column.
Idents(analysis_obj) <- "seurat_clusters"

# Re-apply the cell type names from your Method 2 script
new_cluster_ids <- c("0" = "Oligodendrocytes",
                     "1" = "Astrocytes",
                     "2" = "Astrocytes (Reactive)",
                     "3" = "GABAergic Neurons",
                     "4" = "Dopaminergic Neurons",
                     "5" = "Endothelial Cells",
                     "6" = "Dopaminergic Neurons (DAT+)",
                     "7" = "Neurons",
                     "8" = "Dopaminergic Neurons (Calbindin+)")

# Rename the Idents
analysis_obj <- RenameIdents(analysis_obj, new_cluster_ids)
# --- END FIX ---


# --- a. Calculate proportions ---
# This will now correctly use the cell type names as the rows of the table
prop_table <- table(Idents(analysis_obj), analysis_obj$region_condition)

# Convert to proportions (percentages)
prop_table_percent <- prop.table(prop_table, margin = 2) * 100

# Convert to a data frame for plotting
plot_data <- as.data.frame(prop_table_percent)
colnames(plot_data) <- c("CellType", "Region_Condition", "Percentage")

# --- b. Create the plot ---
# This plot will now be correct
ggplot(plot_data, aes(x = Region_Condition, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Cell Type Proportions by Region and Condition",
    x = "Region",
    y = "Percentage of Cells"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_viridis_d(option = "D")



# =========================================================================
# 6. ANALYZE THE "BLURRY IDENTITY" CLUSTER (Dopaminergic Neurons)
# =========================================================================

# --- a. Subset the main object for only "Dopaminergic Neurons" ---
# This pulls this cell type from all 6 of your anatomical regions
da_neurons_obj <- subset(analysis_obj, idents = "Dopaminergic Neurons")

message(paste("Created a new object with", ncol(da_neurons_obj), "Dopaminergic Neurons."))

# --- b. Prepare the subset object for DGE ---
# Set the identity to 'condition' to compare Schizophrenia vs. Control
Idents(da_neurons_obj) <- "condition"
da_neurons_obj <- PrepSCTFindMarkers(da_neurons_obj)

# --- c. Run the Differential Gene Expression test ---
# This finds what's different between Schizophrenia and Control
# *only within this specific cell population*
da_neurons_dge <- FindMarkers(
  da_neurons_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control",
  logfc.threshold = 0.25, # You can adjust this threshold
  min.pct = 0.1 # Find genes expressed in at least 10% of cells in one group
)

# --- d. View the top results ---
print("--- Top Dysregulated Genes in 'Dopaminergic Neurons' (Schizophrenia vs. Control) ---")
print(head(da_neurons_dge, 20))







# =========================================================================
# ANALYSIS 1: DIFFERENTIAL GENE EXPRESSION (DGE) BETWEEN REGIONS
# =========================================================================
# This assumes your 'analysis_obj' is already loaded and contains
# all 6 subregions.

# --- a. Set the active identity ---
# We tell Seurat to group cells by our 'region_condition' labels
# (e.g., "Control_SNcd", "Schizophrenia_SNcd")
Idents(analysis_obj) <- "region_condition"

# --- b. Prepare the object for DGE (if using SCT) ---
# This step is good practice when working with SCTransform data
analysis_obj <- PrepSCTFindMarkers(analysis_obj)

# --- c. Run DGE for the SNcd ---
# We compare the schizophrenia SNcd directly against the control SNcd
message("--- Running DGE for SNcd (Schizophrenia vs. Control) ---")
dge_sncd <- FindMarkers(
  analysis_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25, # Find genes with at least a 0.25 log-fold change
  min.pct = 0.1           # Find genes expressed in at least 10% of cells in one group
)

# --- d. View the top results for SNcd ---
# This prints the 20 most significant genes that are different in the SNcd
print("--- Top Dysregulated Genes in SNcd ---")
print(head(dge_sncd, 20))


# =========================================================================
# ANALYSIS 1 (Continued): RUN DGE FOR SNcd and VTA
# =========================================================================

# --- Run DGE for the SNcd ---
message("--- Running DGE for SNcd (Schizophrenia vs. Control) ---")
dge_sncd <- FindMarkers(
  analysis_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd ---")
print(head(dge_sncd, 20))


# --- Run DGE for the VTA ---
message("--- Running DGE for VTA (Schizophrenia vs. Control) ---")
dge_vta <- FindMarkers(
  analysis_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA ---")
print(head(dge_vta, 20))



# =========================================================================
# ANALYSIS 1 (Continued): COMPARE CELLS OUTSIDE THE ROIs
# =========================================================================
# This assumes 'merged_obj' (all cells) and 'analysis_obj' (the 6 ROIs)
# are loaded in your R session.

# --- a. Get the names (barcodes) of all cells in your ROIs ---
cells_in_rois <- colnames(analysis_obj)

# --- b. Get the names of all cells in the entire dataset ---
all_cells <- colnames(merged_obj)

# --- c. Find the cells that are OUTSIDE the ROIs ---
cells_outside_rois <- setdiff(all_cells, cells_in_rois)
message(paste("Found", length(cells_outside_rois), "cells outside the defined ROIs."))

# --- d. Create a new Seurat object for 'Outside' cells ---
outside_obj <- subset(merged_obj, cells = cells_outside_rois)

# --- e.Set the identity to 'condition' for DGE ---
Idents(outside_obj) <- "condition"
outside_obj <- PrepSCTFindMarkers(outside_obj)

# --- f. Run DGE for the 'Outside' cells ---
message("--- Running DGE for 'Outside' cells (Schizophrenia vs. Control) ---")
dge_outside <- FindMarkers(
  outside_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# --- g. View the top results for 'Outside' cells ---
print("--- Top Dysregulated Genes in 'Outside' Cells ---")
print(head(dge_outside, 20))




# =========================================================================
# ANALYSIS 2: FOCUSED DGE *WITHIN* A SINGLE CELL TYPE
# =========================================================================

# --- a. Subset for all Astrocytes (Corrected) ---
astrocyte_obj <- subset(analysis_obj, subset = seurat_clusters %in% c(1, 2))

message(paste("Created a new object with", ncol(astrocyte_obj), "total astrocytes."))


# --- !! THE FAILSAFE FIX IS HERE !! ---
# --- b. Prepare the astrocyte object for DGE using the RNA assay ---

# Switch the default assay from "SCT" to "RNA"
DefaultAssay(astrocyte_obj) <- "RNA"

# Run the classic normalization workflow on the raw RNA data
# This bypasses all the SCT model errors.
astrocyte_obj <- NormalizeData(astrocyte_obj, verbose = FALSE)
astrocyte_obj <- FindVariableFeatures(astrocyte_obj, verbose = FALSE)
astrocyte_obj <- ScaleData(astrocyte_obj, verbose = FALSE)

# Now, set the identity to 'region_condition'
Idents(astrocyte_obj) <- "region_condition"

# =========================================================================
# Run DGE for Astrocytes in each region
# =========================================================================
# This will now run on the RNA assay without errors.

# --- c. DGE for Astrocytes in SNcd ---
message("--- Running DGE for Astrocytes in SNcd (Schiz vs. Control) ---")
dge_astro_sncd <- FindMarkers(
  astrocyte_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Astrocytes ---")
print(head(dge_astro_sncd, 20))

# --- d. DGE for Astrocytes in SNcv ---
message("--- Running DGE for Astrocytes in SNcv (Schiz vs. Control) ---")
dge_astro_sncv <- FindMarkers(
  astrocyte_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Astrocytes ---")
print(head(dge_astro_sncv, 20))

# --- e. DGE for Astrocytes in VTA ---
message("--- Running DGE for Astrocytes in VTA (Schiz vs. Control) ---")
dge_astro_vta <- FindMarkers(
  astrocyte_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Astrocytes ---")
print(head(dge_astro_vta, 20))




# =========================================================================
# ANALYSIS 3: FOCUSED DGE USING ATLAS MAPPING (SST+ NEURONS)
# =========================================================================
# This assumes the cell lists (e.g., 'cells_in_sn1_sncd') from your
# previous script are still in your R environment.

# --- a. Load the ONE Mapped Object ---
# This object has ALL cells and ALL metadata (including atlas predictions)
message("Loading atlas mapping data...")
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")

# --- b. Combine all your ROI cell lists ---
# We make one big list of all the cell barcodes you want to keep
all_roi_cells <- c(
  cells_in_sn1_sncd,
  cells_in_sn1_sncv,
  cells_in_sn1_vta,
  cells_in_sn4_sncd,
  cells_in_sn4_sncv,
  cells_in_sn4_vta
)

# --- c. Create the final analysis object ---
# We subset the main mapped object ONE time.
# This keeps all the cell names and metadata perfectly intact.
analysis_obj <- subset(mapped_obj, cells = all_roi_cells)

# --- d. Add the 'region_condition' labels ---
# We need to re-create the labels for this new object
analysis_obj$region_condition <- "Unknown" # Initialize column
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_sncd] <- "Control_SNcd"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_sncv] <- "Control_SNcv"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_vta] <- "Control_VTA"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_sncd] <- "Schizophrenia_SNcd"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_sncv] <- "Schizophrenia_SNcv"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_vta] <- "Schizophrenia_VTA"

# Set the factor levels for correct plot ordering
analysis_obj$region_condition <- factor(analysis_obj$region_condition, levels = c(
  "Control_SNcd", "Schizophrenia_SNcd",
  "Control_SNcv", "Schizophrenia_SNcv",
  "Control_VTA", "Schizophrenia_VTA"
))

message("Final analysis object created successfully.")
# --- END FIX ---


# --- e. Subset for Atlas-Defined SST+ Neurons ---
# Now this command will work because 'predicted.predicted_celltype' exists
sst_neuron_obj <- subset(analysis_obj, subset = predicted.predicted_celltype == 14)

message(paste("Created a new object with", ncol(sst_neuron_obj), "total SST+ Neurons."))

# --- f. Prepare the SST+ object for DGE using the RNA assay ---
DefaultAssay(sst_neuron_obj) <- "RNA"
sst_neuron_obj <- NormalizeData(sst_neuron_obj, verbose = FALSE)
sst_neuron_obj <- FindVariableFeatures(sst_neuron_obj, verbose = FALSE)
sst_neuron_obj <- ScaleData(sst_neuron_obj, verbose = FALSE)

# Set the identity to 'region_condition'
Idents(sst_neuron_obj) <- "region_condition"

# =========================================================================
# Run DGE for SST+ Neurons in each region
# =========================================================================

# --- g. DGE for SST+ Neurons in SNcd ---
# This is the only comparison that will work
message("--- Running DGE for SST+ Neurons in SNcd (Schiz vs. Control) ---")
dge_sst_sncd <- FindMarkers(
  sst_neuron_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd SST+ Neurons ---")
print(head(dge_sst_sncd, 20))

# --- h. DGE for SST+ Neurons in SNcv ---
# This will fail because 'Schizophrenia_SNcv' is empty
# We will comment it out
# message("--- Running DGE for SST+ Neurons in SNcv (Schiz vs. Control) ---")
# dge_sst_sncv <- FindMarkers(...)

# --- i. DGE for SST+ Neurons in VTA ---
# This will fail because 'Control_VTA' is empty
# We will comment it out
# message("--- Running DGE for SST+ Neurons in VTA (Schiz vs. Control) ---")
# dge_sst_vta <- FindMarkers(...)





# =========================================================================
# ANALYSIS 4: FOCUSED DGE *WITHIN* "GABAergic Neurons" (Cluster 3)
# =========================================================================

# --- a. Subset for "GABAergic Neurons" (Cluster 3) ---
# We use the 'subset' argument to find cluster 3
gaba_obj <- subset(analysis_obj, subset = seurat_clusters == 3)

message(paste("Created a new object with", ncol(gaba_obj), "total GABAergic Neurons."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the GABA object for DGE using the RNA assay ---

# Switch the default assay from "SCT" to "RNA"
DefaultAssay(gaba_obj) <- "RNA"

# Run the classic normalization workflow
gaba_obj <- NormalizeData(gaba_obj, verbose = FALSE)
gaba_obj <- FindVariableFeatures(gaba_obj, verbose = FALSE)
gaba_obj <- ScaleData(gaba_obj, verbose = FALSE)

# Now, set the identity to 'region_condition'
Idents(gaba_obj) <- "region_condition"

# =========================================================================
# Run DGE for GABAergic Neurons in each region
# =========================================================================

# --- c. DGE for GABAergic Neurons in SNcd ---
message("--- Running DGE for GABAergic Neurons in SNcd (Schiz vs. Control) ---")
dge_gaba_sncd <- FindMarkers(
  gaba_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd GABAergic Neurons ---")
print(head(dge_gaba_sncd, 20))

# --- d. DGE for GABAergic Neurons in SNcv ---
message("--- Running DGE for GABAergic Neurons in SNcv (Schiz vs. Control) ---")
dge_gaba_sncv <- FindMarkers(
  gaba_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv GABAergic Neurons ---")
print(head(dge_gaba_sncv, 20))

# --- e. DGE for GABAergic Neurons in VTA ---
message("--- Running DGE for GABAergic Neurons in VTA (Schiz vs. Control) ---")
dge_gaba_vta <- FindMarkers(
  gaba_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA GABAergic Neurons ---")
print(head(dge_gaba_vta, 20))




# =========================================================================
# ANALYSIS 5: FOCUSED DGE *WITHIN* "Dopaminergic Neurons" (Cluster 4)
# =========================================================================

# --- a. Subset for "Dopaminergic Neurons" (Cluster 4) ---
# We use the 'subset' argument to find cluster 4
da_obj <- subset(analysis_obj, subset = seurat_clusters == 4)

message(paste("Created a new object with", ncol(da_obj), "total 'Dopaminergic Neurons'."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the DA object for DGE using the RNA assay ---

# Switch the default assay from "SCT" to "RNA"
DefaultAssay(da_obj) <- "RNA"

# Run the classic normalization workflow
da_obj <- NormalizeData(da_obj, verbose = FALSE)
da_obj <- FindVariableFeatures(da_obj, verbose = FALSE)
da_obj <- ScaleData(da_obj, verbose = FALSE)

# Now, set the identity to 'region_condition'
Idents(da_obj) <- "region_condition"

# =========================================================================
# Run DGE for Dopaminergic Neurons in each region
# =========================================================================

# --- c. DGE for DA Neurons in SNcd ---
message("--- Running DGE for DA Neurons in SNcd (Schiz vs. Control) ---")
dge_da_sncd <- FindMarkers(
  da_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd DA Neurons ---")
print(head(dge_da_sncd, 20))

# --- d. DGE for DA Neurons in SNcv ---
message("--- Running DGE for DA Neurons in SNcv (Schiz vs. Control) ---")
dge_da_sncv <- FindMarkers(
  da_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv DA Neurons ---")
print(head(dge_da_sncv, 20))

# --- e. DGE for DA Neurons in VTA ---
message("--- Running DGE for DA Neurons in VTA (Schiz vs. Control) ---")
dge_da_vta <- FindMarkers(
  da_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA DA Neurons ---")
print(head(dge_da_vta, 20))





# =========================================================================
# ANALYSIS 6: FOCUSED DGE *WITHIN* "Dopaminergic Neurons (DAT+)" (Cluster 6)
# =========================================================================

# --- a. Subset for "Dopaminergic Neurons (DAT+)" (Cluster 6) ---
dat_pos_obj <- subset(analysis_obj, subset = seurat_clusters == 6)

message(paste("Created a new object with", ncol(dat_pos_obj), "total 'Dopaminergic Neurons (DAT+)'."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the DAT+ object for DGE using the RNA assay ---

# Switch the default assay from "SCT" to "RNA"
DefaultAssay(dat_pos_obj) <- "RNA"

# Run the classic normalization workflow
dat_pos_obj <- NormalizeData(dat_pos_obj, verbose = FALSE)
dat_pos_obj <- FindVariableFeatures(dat_pos_obj, verbose = FALSE)
dat_pos_obj <- ScaleData(dat_pos_obj, verbose = FALSE)

# Now, set the identity to 'region_condition'
Idents(dat_pos_obj) <- "region_condition"

# =========================================================================
# Run DGE for DAT+ Neurons in each region
# =========================================================================

# --- c. DGE for DAT+ Neurons in SNcd ---
message("--- Running DGE for DAT+ Neurons in SNcd (Schiz vs. Control) ---")
dge_dat_sncd <- FindMarkers(
  dat_pos_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd DAT+ Neurons ---")
print(head(dge_dat_sncd, 20))

# --- d. DGE for DAT+ Neurons in SNcv ---
message("--- Running DGE for DAT+ Neurons in SNcv (Schiz vs. Control) ---")
dge_dat_sncv <- FindMarkers(
  dat_pos_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv DAT+ Neurons ---")
print(head(dge_dat_sncv, 20))

# --- e. DGE for DAT+ Neurons in VTA ---
message("--- Running DGE for DAT+ Neurons in VTA (Schiz vs. Control) ---")
dge_dat_vta <- FindMarkers(
  dat_pos_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA DAT+ Neurons ---")
print(head(dge_dat_vta, 20))




# =========================================================================
# ANALYSIS 7: FOCUSED DGE *WITHIN* "Dopaminergic Neurons (Calbindin+)" (Cluster 8)
# =========================================================================

# --- a. Subset for "Dopaminergic Neurons (Calbindin+)" (Cluster 8) ---
calbindin_pos_obj <- subset(analysis_obj, subset = seurat_clusters == 8)

message(paste("Created a new object with", ncol(calbindin_pos_obj), "total 'Dopaminergic Neurons (Calbindin+)'."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Calbindin+ object for DGE using the RNA assay ---

# Switch the default assay from "SCT" to "RNA"
DefaultAssay(calbindin_pos_obj) <- "RNA"

# Run the classic normalization workflow
calbindin_pos_obj <- NormalizeData(calbindin_pos_obj, verbose = FALSE)
calbindin_pos_obj <- FindVariableFeatures(calbindin_pos_obj, verbose = FALSE)
calbindin_pos_obj <- ScaleData(calbindin_pos_obj, verbose = FALSE)

# Now, set the identity to 'region_condition'
Idents(calbindin_pos_obj) <- "region_condition"

# =========================================================================
# Run DGE for Calbindin+ Neurons in each region
# =========================================================================

# --- c. DGE for Calbindin+ Neurons in SNcd ---
#message("--- Running DGE for Calbindin+ Neurons in SNcd (Schiz vs. Control) ---")
#dge_calbindin_sncd <- FindMarkers(
  #calbindin_pos_obj,
  #ident.1 = "Schizophrenia_SNcd",
  #ident.2 = "Control_SNcd",
  #logfc.threshold = 0.25,
 # min.pct = 0.1
#)
#print("--- Top Dysregulated Genes in SNcd Calbindin+ Neurons ---")
#print(head(dge_calbindin_sncd, 20))

# --- d. DGE for Calbindin+ Neurons in SNcv ---
#message("--- Running DGE for Calbindin+ Neurons in SNcv (Schiz vs. Control) ---")
#dge_calbindin_sncv <- FindMarkers(
 ## calbindin_pos_obj,
  #ident.1 = "Schizophrenia_SNcv",
  #ident.2 = "Control_SNcv",
  #logfc.threshold = 0.25,
  #min.pct = 0.1
#)
#print("--- Top Dysregulated Genes in SNcv Calbindin+ Neurons ---")
#print(head(dge_calbindin_sncv, 20))

# --- e. DGE for Calbindin+ Neurons in VTA ---
#message("--- Running DGE for Calbindin+ Neurons in VTA (Schiz vs. Control) ---")
#dge_calbindin_vta <- FindMarkers(
 # calbindin_pos_obj,
  #ident.1 = "Schizophrenia_VTA",
  #ident.2 = "Control_VTA",
 # logfc.threshold = 0.25,
 # min.pct = 0.1
#)
#print("--- Top Dysregulated Genes in VTA Calbindin+ Neurons ---")
#print(head(dge_calbindin_vta, 20))



# =========================================================================
# ANALYSIS 8: FOCUSED DGE USING ATLAS-DEFINED ASTROCYTES
# =========================================================================
# This assumes 'analysis_obj' is still loaded, which contains
# the 'predicted.predicted_celltype' column.

# --- a. Subset for Atlas-Defined Astrocytes (Clusters 0 and 6) ---
atlas_astro_obj <- subset(analysis_obj, subset = predicted.predicted_celltype %in% c(0, 6))

message(paste("Created a new object with", ncol(atlas_astro_obj), "total Atlas-Defined Astrocytes."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Atlas-Astrocyte object for DGE using the RNA assay ---

DefaultAssay(atlas_astro_obj) <- "RNA"
atlas_astro_obj <- NormalizeData(atlas_astro_obj, verbose = FALSE)
atlas_astro_obj <- FindVariableFeatures(atlas_astro_obj, verbose = FALSE)
atlas_astro_obj <- ScaleData(atlas_astro_obj, verbose = FALSE)

Idents(atlas_astro_obj) <- "region_condition"

# =========================================================================
# Run DGE for Atlas-Defined Astrocytes in each region
# =========================================================================

# --- c. DGE for Atlas Astrocytes in SNcd ---
message("--- Running DGE for Atlas Astrocytes in SNcd (Schiz vs. Control) ---")
dge_atlas_astro_sncd <- FindMarkers(
  atlas_astro_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Atlas Astrocytes ---")
print(head(dge_atlas_astro_sncd, 20))

# --- d. DGE for Atlas Astrocytes in SNcv ---
message("--- Running DGE for Atlas Astrocytes in SNcv (Schiz vs. Control) ---")
dge_atlas_astro_sncv <- FindMarkers(
  atlas_astro_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Atlas Astrocytes ---")
print(head(dge_atlas_astro_sncv, 20))

# --- e. DGE for Atlas Astrocytes in VTA ---
message("--- Running DGE for Atlas Astrocytes in VTA (Schiz vs. Control) ---")
dge_atlas_astro_vta <- FindMarkers(
  atlas_astro_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Atlas Astrocytes ---")
print(head(dge_atlas_astro_vta, 20))





# =========================================================================
# ANALYSIS 9: FOCUSED DGE USING ATLAS-DEFINED PVALB+ NEURONS
# =========================================================================
# This assumes 'analysis_obj' is still loaded, which contains
# the 'predicted.predicted_celltype' column.

# --- a. Subset for Atlas-Defined PVALB+ Neurons (Clusters 5 and 21) ---
atlas_pvalb_obj <- subset(analysis_obj, subset = predicted.predicted_celltype %in% c(5, 21))

message(paste("Created a new object with", ncol(atlas_pvalb_obj), "total Atlas-Defined PVALB+ Neurons."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Atlas-PVALB object for DGE using the RNA assay ---

DefaultAssay(atlas_pvalb_obj) <- "RNA"
atlas_pvalb_obj <- NormalizeData(atlas_pvalb_obj, verbose = FALSE)
atlas_pvalb_obj <- FindVariableFeatures(atlas_pvalb_obj, verbose = FALSE)
atlas_pvalb_obj <- ScaleData(atlas_pvalb_obj, verbose = FALSE)

Idents(atlas_pvalb_obj) <- "region_condition"

# =========================================================================
# Run DGE for Atlas-Defined PVALB+ Neurons in each region
# =========================================================================

# --- c. DGE for Atlas PVALB+ Neurons in SNcd ---
message("--- Running DGE for Atlas PVALB+ Neurons in SNcd (Schiz vs. Control) ---")
dge_atlas_pvalb_sncd <- FindMarkers(
  atlas_pvalb_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Atlas PVALB+ Neurons ---")
print(head(dge_atlas_pvalb_sncd, 20))

# --- d. DGE for Atlas PVALB+ Neurons in SNcv ---
message("--- Running DGE for Atlas PVALB+ Neurons in SNcv (Schiz vs. Control) ---")
dge_atlas_pvalb_sncv <- FindMarkers(
  atlas_pvalb_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Atlas PVALB+ Neurons ---")
print(head(dge_atlas_pvalb_sncv, 20))

# --- e. DGE for Atlas PVALB+ Neurons in VTA ---
message("--- Running DGE for Atlas PVALB+ Neurons in VTA (Schiz vs. Control) ---")
dge_atlas_pvalb_vta <- FindMarkers(
  atlas_pvalb_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Atlas PVALB+ Neurons ---")
print(head(dge_atlas_pvalb_vta, 20))



# =========================================================================
# ANALYSIS 10: FOCUSED DGE *WITHIN* "Oligodendrocytes" (Cluster 0)
# =========================================================================
# This assumes 'analysis_obj' is still loaded.

# --- a. Subset for "Oligodendrocytes" (Cluster 0) ---
oligo_obj <- subset(analysis_obj, subset = seurat_clusters == 0)

message(paste("Created a new object with", ncol(oligo_obj), "total Oligodendrocytes."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Oligo object for DGE using the RNA assay ---

DefaultAssay(oligo_obj) <- "RNA"
oligo_obj <- NormalizeData(oligo_obj, verbose = FALSE)
oligo_obj <- FindVariableFeatures(oligo_obj, verbose = FALSE)
oligo_obj <- ScaleData(oligo_obj, verbose = FALSE)

Idents(oligo_obj) <- "region_condition"

# =========================================================================
# Run DGE for Oligodendrocytes in each region
# =========================================================================

# --- c. DGE for Oligodendrocytes in SNcd ---
message("--- Running DGE for Oligodendrocytes in SNcd (Schiz vs. Control) ---")
dge_oligo_sncd <- FindMarkers(
  oligo_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Oligodendrocytes ---")
print(head(dge_oligo_sncd, 20))

# --- d. DGE for Oligodendrocytes in SNcv ---
message("--- Running DGE for Oligodendrocytes in SNcv (Schiz vs. Control) ---")
dge_oligo_sncv <- FindMarkers(
  oligo_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Oligodendrocytes ---")
print(head(dge_oligo_sncv, 20))

# --- e. DGE for Oligodendrocytes in VTA ---
message("--- Running DGE for Oligodendrocytes in VTA (Schiz vs. Control) ---")
dge_oligo_vta <- FindMarkers(
  oligo_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Oligodendrocytes ---")
print(head(dge_oligo_vta, 20))




# =========================================================================
# ANALYSIS 11: FOCUSED DGE *WITHIN* "Endothelial Cells" (Cluster 5)
# =========================================================================
# This assumes 'analysis_obj' is still loaded.

# --- a. Subset for "Endothelial Cells" (Cluster 5) ---
endo_obj <- subset(analysis_obj, subset = seurat_clusters == 5)

message(paste("Created a new object with", ncol(endo_obj), "total Endothelial Cells."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Endothelial object for DGE using the RNA assay ---

DefaultAssay(endo_obj) <- "RNA"
endo_obj <- NormalizeData(endo_obj, verbose = FALSE)
endo_obj <- FindVariableFeatures(endo_obj, verbose = FALSE)
endo_obj <- ScaleData(endo_obj, verbose = FALSE)

Idents(endo_obj) <- "region_condition"

# =========================================================================
# Run DGE for Endothelial Cells in each region
# =========================================================================

# --- c. DGE for Endothelial Cells in SNcd ---
message("--- Running DGE for Endothelial Cells in SNcd (Schiz vs. Control) ---")
dge_endo_sncd <- FindMarkers(
  endo_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Endothelial Cells ---")
print(head(dge_endo_sncd, 20))

# --- d. DGE for Endothelial Cells in SNcv ---
message("--- Running DGE for Endothelial Cells in SNcv (Schiz vs. Control) ---")
dge_endo_sncv <- FindMarkers(
  endo_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Endothelial Cells ---")
print(head(dge_endo_sncv, 20))

# --- e. DGE for Endothelial Cells in VTA ---
message("--- Running DGE for Endothelial Cells in VTA (Schiz vs. Control) ---")
dge_endo_vta <- FindMarkers(
  endo_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Endothelial Cells ---")
print(head(dge_endo_vta, 20))






# =========================================================================
# ANALYSIS 12: FOCUSED DGE USING *ATLAS-DEFINED* DA NEURONS
# =========================================================================
# This assumes 'analysis_obj' is still loaded.

# --- a. Subset for Atlas-Defined DA Neurons (Clusters 4, 10, 11, 15, 17, 19) ---
# This list is based on the Agarwal et al. paper (GSE140231) and your TH FeaturePlot
atlas_da_obj <- subset(analysis_obj, subset = predicted.predicted_celltype %in% c(4, 10, 11, 15, 17, 19))

message(paste("Created a new object with", ncol(atlas_da_obj), "total Atlas-Defined DA Neurons."))

# --- !! THE FAILSAFE FIX !! ---
# --- b. Prepare the Atlas-DA object for DGE using the RNA assay ---

DefaultAssay(atlas_da_obj) <- "RNA"
atlas_da_obj <- NormalizeData(atlas_da_obj, verbose = FALSE)
atlas_da_obj <- FindVariableFeatures(atlas_da_obj, verbose = FALSE)
atlas_da_obj <- ScaleData(atlas_da_obj, verbose = FALSE)

Idents(atlas_da_obj) <- "region_condition"

# =========================================================================
# Run DGE for Atlas-Defined DA Neurons in each region
# =========================================================================

# --- c. DGE for Atlas DA Neurons in SNcd ---
message("--- Running DGE for Atlas DA Neurons in SNcd (Schiz vs. Control) ---")
dge_atlas_da_sncd <- FindMarkers(
  atlas_da_obj,
  ident.1 = "Schizophrenia_SNcd",
  ident.2 = "Control_SNcd",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcd Atlas DA Neurons ---")
print(head(dge_atlas_da_sncd, 20))

# --- d. DGE for Atlas DA Neurons in SNcv ---
message("--- Running DGE for Atlas DA Neurons in SNcv (Schiz vs. Control) ---")
dge_atlas_da_sncv <- FindMarkers(
  atlas_da_obj,
  ident.1 = "Schizophrenia_SNcv",
  ident.2 = "Control_SNcv",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in SNcv Atlas DA Neurons ---")
print(head(dge_atlas_da_sncv, 20))

# --- e. DGE for Atlas DA Neurons in VTA ---
message("--- Running DGE for Atlas DA Neurons in VTA (Schiz vs. Control) ---")
dge_atlas_da_vta <- FindMarkers(
  atlas_da_obj,
  ident.1 = "Schizophrenia_VTA",
  ident.2 = "Control_VTA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
print("--- Top Dysregulated Genes in VTA Atlas DA Neurons ---")
print(head(dge_atlas_da_vta, 20))




message("Loading libraries (clusterProfiler, org.Hs.eg.db)...")
library(clusterProfiler)
library(org.Hs.eg.db) # Annotation database for Homo sapiens (human)


# 2. HELPER FUNCTION (Base R)
# ---------------------------------------------------------------------
# This function prepares your DGE data for GSEA using Base R
# to avoid package conflicts.
prepare_gsea_list_base_r <- function(dge_table) {
  # Convert rownames to a 'gene' column
  dge_df <- as.data.frame(dge_table)
  dge_df$gene <- rownames(dge_df)
  
  # Filter out 'Codeword' genes
  dge_df <- dge_df[!grepl("Codeword", dge_df$gene), ]
  
  # Select only the gene and avg_log2FC columns
  dge_df <- dge_df[, c("gene", "avg_log2FC")]
  
  # Create the ranked list
  gene_list <- dge_df$avg_log2FC
  names(gene_list) <- dge_df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(
    names(gene_list),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Match the ranked list to the newly mapped Entrez IDs
  ranked_gene_list <- gene_list[entrez_ids$SYMBOL]
  names(ranked_gene_list) <- entrez_ids$ENTREZID
  
  # Remove any NA values that failed to map
  ranked_gene_list <- ranked_gene_list[!is.na(names(ranked_gene_list))]
  
  return(ranked_gene_list)
}

# This function extracts top 5 results using Base R
extract_top_5_base_r <- function(gsea_result) {
  gsea_df <- as.data.frame(gsea_result)
  
  # Sort all results by the absolute value of NES
  order_indices <- order(abs(gsea_df$NES), decreasing = TRUE)
  sorted_results <- gsea_df[order_indices, ]
  
  # Extract the top 5 rows
  top_5_overall <- head(sorted_results, 5)
  return(top_5_overall)
}


# 3. ANALYSIS 1: SNcd DA NEURONS
# ---------------------------------------------------------------------
message("--- Running GSEA for SNcd DA Neurons ---")
ranked_list_da_sncd <- prepare_gsea_list_base_r(dge_atlas_da_sncd)
gsea_da_sncd <- gseKEGG(
  geneList = ranked_list_da_sncd,
  organism = 'hsa', minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.5, verbose = FALSE, seed = 42
)
top_5_da_sncd <- extract_top_5_base_r(gsea_da_sncd)
print("--- TOP 5 OVERALL GSEA PATHWAYS (SNcd Atlas DA Neurons) ---")
print(top_5_da_sncd)


# 4. ANALYSIS 2: SNcv DA NEURONS
# ---------------------------------------------------------------------
message("--- Running GSEA for SNcv DA Neurons ---")
ranked_list_da_sncv <- prepare_gsea_list_base_r(dge_atlas_da_sncv)
gsea_da_sncv <- gseKEGG(
  geneList = ranked_list_da_sncv,
  organism = 'hsa', minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.5, verbose = FALSE, seed = 42
)
top_5_da_sncv <- extract_top_5_base_r(gsea_da_sncv)
print("--- TOP 5 OVERALL GSEA PATHWAYS (SNcv Atlas DA Neurons) ---")
print(top_5_da_sncv)


# 5. ANALYSIS 3: SNcd ASTROCYTES
# ---------------------------------------------------------------------
message("--- Running GSEA for SNcd Astrocytes ---")
ranked_list_astro_sncd <- prepare_gsea_list_base_r(dge_atlas_astro_sncd)
gsea_astro_sncd <- gseKEGG(
  geneList = ranked_list_astro_sncd,
  organism = 'hsa', minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.5, verbose = FALSE, seed = 42
)
top_5_astro_sncd <- extract_top_5_base_r(gsea_astro_sncd)
print("--- TOP 5 OVERALL GSEA PATHWAYS (SNcd Atlas Astrocytes) ---")
print(top_5_astro_sncd)


# 6. ANALYSIS 4: SNcv ASTROCYTES
# ---------------------------------------------------------------------
message("--- Running GSEA for SNcv Astrocytes ---")
ranked_list_astro_sncv <- prepare_gsea_list_base_r(dge_atlas_astro_sncv)
gsea_astro_sncv <- gseKEGG(
  geneList = ranked_list_astro_sncv,
  organism = 'hsa', minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.5, verbose = FALSE, seed = 42
)
top_5_astro_sncv <- extract_top_5_base_r(gsea_astro_sncv)
print("--- TOP 5 OVERALL GSEA PATHWAYS (SNcv Atlas Astrocytes) ---")
print(top_5_astro_sncv)





















library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(merged_obj) <- "RNA"
# If not already normalized/scaled, uncomment:
# merged_obj <- NormalizeData(merged_obj, verbose = FALSE) |>
#   FindVariableFeatures(verbose = FALSE) |>
#   ScaleData(verbose = FALSE)

# <<< pick ONE: >>>
cluster_col <- "seurat_clusters"   # or "cluster"

Idents(merged_obj) <- cluster_col
clusters <- sort(unique(Idents(merged_obj)))

min.cells <- 10
dge_list <- list()

for (cl in clusters) {
  cl_obj <- subset(merged_obj, idents = cl)
  Idents(cl_obj) <- "condition"
  # guard for very small groups
  n_per_cond <- table(Idents(cl_obj))
  if (all(n_per_cond >= min.cells)) {
    res <- FindMarkers(
      cl_obj, ident.1 = "Schizophrenia", ident.2 = "Control",
      logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox"
    )
    if (nrow(res) > 0) {
      res$gene <- rownames(res)
      res$cluster <- as.character(cl)
      res$padj <- p.adjust(res$p_val, method = "BH")
      dge_list[[as.character(cl)]] <- res
    }
  }
}

dge_all <- if (length(dge_list)) bind_rows(dge_list) else tibble()
# summary: how many sig genes per cluster
sig_summary <- dge_all %>%
  group_by(cluster) %>%
  summarise(n_sig = sum(padj < 0.05), .groups = "drop") %>%
  arrange(desc(n_sig))
sig_summary




plot_top10 <- function(df, cl) {
  df %>%
    filter(cluster == cl, padj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    slice_head(n = 10) %>%
    mutate(direction = ifelse(avg_log2FC >= 0, "Up (SZ)", "Down (SZ)"),
           gene = factor(gene, levels = rev(gene))) %>%
    ggplot(aes(x = gene, y = avg_log2FC, fill = direction)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Top dysregulated genes • cluster", cl),
         y = "log2 FC (SZ vs Ctrl)", x = NULL) +
    scale_fill_manual(values = c("Down (SZ)" = "#377eb8", "Up (SZ)" = "#e41a1c")) +
    theme_bw(base_size = 12)
}
# Example: plot_top10(dge_all, cl = sig_summary$cluster[1])







library(ggplot2)
library(dplyr)
library(tidyr)

# 1) Stacked % composition per condition
comp_df <- as.data.frame.matrix(table(merged_obj[[cluster_col]][,1], merged_obj$condition))
comp_df$cluster <- rownames(comp_df)
comp_long <- comp_df |>
  pivot_longer(cols = -cluster, names_to = "condition", values_to = "n") |>
  group_by(condition) |>
  mutate(percent = 100 * n / sum(n)) |>
  ungroup()

ggplot(comp_long, aes(x = condition, y = percent, fill = cluster)) +
  geom_col() +
  labs(title = "Cell-type composition by condition",
       y = "Percent of cells", x = NULL) +
  theme_bw(base_size = 12)

# 2) Chi-square + standardized residuals heatmap & sorted bars
tab <- table(merged_obj[[cluster_col]][,1], merged_obj$condition)
cs  <- suppressWarnings(chisq.test(tab))

stdres <- as.data.frame(as.table(cs$stdres)) |>
  rename(cluster = Var1, condition = Var2, stdres = Freq)

ggplot(stdres, aes(condition, cluster, fill = stdres)) +
  geom_tile() + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(title = sprintf("Standardized residuals (chisq p = %.2e, df = %d)", cs$p.value, cs$parameter),
       x = NULL, y = NULL, fill = "Std. resid") +
  theme_bw(base_size = 12)

# Sorted bar of SZ column (over/under-represented in SZ)
std_sz <- stdres %>% filter(condition == "Schizophrenia") %>% arrange(desc(stdres))
ggplot(std_sz, aes(x = reorder(cluster, stdres), y = stdres, fill = stdres > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac")) +
  labs(title = "Over-/under-representation in SZ (Std. residuals, SZ column)",
       x = NULL, y = "Standardized residual") +
  theme_bw(base_size = 12)














# ============================================================
# Rebuild ROI DGE (Atlas clusters), plus Outside-ROI DGE
# and generate volcano plots (SNcd, SNcv, VTA + Outside)
# ============================================================

# ---- 0) Libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
})

set.seed(42)

# ---- 1) Paths / output dirs ----

# ---- 2) Load base objects ----
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")
merged_obj <- readRDS("data/Final_Analyzed_Object.rds")

# ---- 3) Load ROI cell lists from your saved ROI objects ----
get_cells <- function(path) {
  if (file.exists(path)) {
    colnames(readRDS(path))
  } else {
    warning(sprintf("Missing ROI object: %s", path))
    NULL
  }
}

cells_in_sn1_sncd <- get_cells("data/SN1_SNcd_Only_Object.rds")
cells_in_sn1_sncv <- get_cells("data/SN1_SNcv_Only_Object.rds")
cells_in_sn1_vta  <- get_cells("data/SN1_VTA_Only_Object.rds")
cells_in_sn4_sncd <- get_cells("data/SN4_SNcd_Only_Object.rds")
cells_in_sn4_sncv <- get_cells("data/SN4_SNcv_Only_Object.rds")
cells_in_sn4_vta  <- get_cells("data/SN4_VTA_Only_Object.rds")

all_roi_cells <- c(
  cells_in_sn1_sncd, cells_in_sn1_sncv, cells_in_sn1_vta,
  cells_in_sn4_sncd, cells_in_sn4_sncv, cells_in_sn4_vta
)
all_roi_cells <- unique(all_roi_cells)

# ---- 4) Build analysis_obj (only ROI cells) + region_condition ----
analysis_obj <- subset(mapped_obj, cells = all_roi_cells)
analysis_obj$region_condition <- "Unknown"

# tag membership (exclusive)
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_sncd] <- "Control_SNcd"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_sncv] <- "Control_SNcv"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn1_vta]  <- "Control_VTA"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_sncd] <- "Schizophrenia_SNcd"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_sncv] <- "Schizophrenia_SNcv"
analysis_obj$region_condition[colnames(analysis_obj) %in% cells_in_sn4_vta]  <- "Schizophrenia_VTA"

analysis_obj$region_condition <- factor(
  analysis_obj$region_condition,
  levels = c("Control_SNcd","Schizophrenia_SNcd",
             "Control_SNcv","Schizophrenia_SNcv",
             "Control_VTA","Schizophrenia_VTA")
)

# quick sanity
message("ROI cells: ", ncol(analysis_obj))
stopifnot("region_condition" %in% colnames(analysis_obj@meta.data))
stopifnot("predicted.predicted_celltype" %in% colnames(analysis_obj@meta.data))

# ---- 5) Helper to prep a Seurat object for DGE ----
prep_for_dge <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj
}

# ---- 6) Run DGE for Atlas clusters by region (SNcd/SNcv/VTA) ----
atlas_dge_results <- list()

atlas_clusters <- sort(unique(analysis_obj$predicted.predicted_celltype))
message("Atlas clusters present in ROI: {", paste(atlas_clusters, collapse = ","), "}")

for (cl in atlas_clusters) {
  # subset ROI cells that belong to this atlas cluster
  sub <- subset(analysis_obj, subset = predicted.predicted_celltype == cl)
  if (ncol(sub) < 10) {
    message(sprintf("Skip cluster %s (too few cells: %d)", cl, ncol(sub)))
    next
  }
  sub <- prep_for_dge(sub)
  Idents(sub) <- "region_condition"
  present <- names(table(Idents(sub)))
  
  do_comp <- function(id1, id2) {
    if (all(c(id1,id2) %in% present)) {
      res <- tryCatch({
        FindMarkers(sub, ident.1 = id1, ident.2 = id2,
                    logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox")
      }, error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        res$gene <- rownames(res)
        res$padj <- p.adjust(res$p_val, method = "BH")
      }
      return(res)
    }
    NULL
  }
  
  dge_sncd <- do_comp("Schizophrenia_SNcd", "Control_SNcd")
  dge_sncv <- do_comp("Schizophrenia_SNcv", "Control_SNcv")
  dge_vta  <- do_comp("Schizophrenia_VTA",  "Control_VTA")
  
  atlas_dge_results[[as.character(cl)]] <- list(SNcd = dge_sncd, SNcv = dge_sncv, VTA = dge_vta)
}

# Save full DGE results list
saveRDS(atlas_dge_results, "results/dge/atlas_dge_results_final.rds")

# ---- 7) Make a tidy data.frame of counts of significant genes ----
sig_counts <- lapply(names(atlas_dge_results), function(cl){
  e <- atlas_dge_results[[cl]]
  tibble(
    cluster = cl,
    SNcd = if (!is.null(e$SNcd)) sum(e$SNcd$padj < 0.05, na.rm = TRUE) else 0L,
    SNcv = if (!is.null(e$SNcv)) sum(e$SNcv$padj < 0.05, na.rm = TRUE) else 0L,
    VTA  = if (!is.null(e$VTA )) sum(e$VTA$padj  < 0.05, na.rm = TRUE) else 0L
  )
}) %>% bind_rows()

write_csv(sig_counts, "results/dge/atlas_sig_counts_per_region.csv")
print(sig_counts)

# ---- 8) Pick TOP cluster per region (max # of sig genes) ----
pick_top <- function(region_col) {
  sig_counts %>%
    mutate(n_sig = .data[[region_col]]) %>%
    arrange(desc(n_sig)) %>%
    slice(1) %>%
    pull(cluster) %>%
    as.character()
}

top_SNcd <- pick_top("SNcd")
top_SNcv <- pick_top("SNcv")
top_VTA  <- pick_top("VTA")

message(sprintf("Top clusters -> SNcd: %s | SNcv: %s | VTA: %s",
                top_SNcd, top_SNcv, top_VTA))

# ---- 9) Volcano plot helper ----
is_deprecated <- function(x) grepl("^DeprecatedCodeword|^Codeword", x, ignore.case = TRUE)

make_volcano <- function(res, title_text) {
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  
  df <- res %>%
    mutate(
      padj = ifelse(is.na(padj), 1, padj),
      neglog10p = -log10(padj),
      is_sig = padj < 0.05,
      dir = case_when(
        is_sig & avg_log2FC > 0 ~ "Up",
        is_sig & avg_log2FC < 0 ~ "Down",
        TRUE ~ "NotSig"
      )
    ) %>%
    mutate(gene = rownames(.)) %>%
    filter(!is_deprecated(gene))   # 👈 remove codewords
  
  # label top genes only
  top_labs <- df %>%
    filter(is_sig) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(15)
  
  ggplot(df, aes(x = avg_log2FC, y = neglog10p)) +
    geom_point(aes(color = dir), alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
    geom_text_repel(data = top_labs, aes(label = gene), size = 3, max.overlaps = 15) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey70")) +
    labs(
      title = title_text,
      x = "log2 Fold Change",
      y = "-log10(adj p-value)"
    ) +
    theme_minimal(base_size = 13)
}


# ---- 10) Build one volcano per region (using the top cluster we picked) ----
get_region_dge <- function(cluster_id, region_name){
  if (is.null(cluster_id) || is.na(cluster_id)) return(NULL)
  entry <- atlas_dge_results[[as.character(cluster_id)]]
  if (is.null(entry)) return(NULL)
  entry[[region_name]]
}

vol_SNcd <- make_volcano(get_region_dge(top_SNcd, "SNcd"),
                         sprintf("SNcd — Atlas Cluster %s (top by #sig genes)", top_SNcd))
vol_SNcv <- make_volcano(get_region_dge(top_SNcv, "SNcv"),
                         sprintf("SNcv — Atlas Cluster %s (top by #sig genes)", top_SNcv))
vol_VTA  <- make_volcano(get_region_dge(top_VTA,  "VTA"),
                         sprintf("VTA — Atlas Cluster %s (top by #sig genes)",  top_VTA))

# ---- 11) Outside-ROI DGE + volcano ----
# define outside cells using your ROI membership:
cells_in_rois   <- colnames(analysis_obj)
all_cells_align <- intersect(colnames(merged_obj), colnames(mapped_obj))
cells_outside   <- setdiff(all_cells_align, cells_in_rois)

outside_plot <- NULL
if (length(cells_outside) >= 20) {
  outside_obj <- subset(merged_obj, cells = cells_outside)
  Idents(outside_obj) <- "condition"
  outside_obj <- prep_for_dge(outside_obj)
  
  dge_outside <- tryCatch({
    FindMarkers(outside_obj, ident.1 = "Schizophrenia", ident.2 = "Control",
                logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox")
  }, error = function(e) NULL)
  
  if (!is.null(dge_outside) && nrow(dge_outside) > 0) {
    dge_outside$gene <- rownames(dge_outside)
    dge_outside$padj <- p.adjust(dge_outside$p_val, method = "BH")
    
    # Save outside table
    write_csv(dge_outside, "results/dge/DGE_Outside_vs_ROI.csv")
    
    outside_plot <- make_volcano(dge_outside, "Outside ROI — All Cells (SZ vs Ctrl)")
  } else {
    message("Outside-ROI DGE had no rows or failed; skipping outside volcano.")
  }
} else {
  message("No/too few Outside cells found; skipping Outside-ROI DGE.")
}

# ---- 12) Save plots (individual + 2x2 grid including Outside if present) ----
if (!is.null(vol_SNcd)) ggsave("results/figs/"Volcano_SNcd_TOP.png"), vol_SNcd, width = 6, height = 5, dpi = 300)
if (!is.null(vol_SNcv)) ggsave("results/figs/"Volcano_SNcv_TOP.png"), vol_SNcv, width = 6, height = 5, dpi = 300)
if (!is.null(vol_VTA )) ggsave("results/figs/"Volcano_VTA_TOP.png"),  vol_VTA,  width = 6, height = 5, dpi = 300)
if (!is.null(outside_plot)) ggsave("results/figs/"Volcano_Outside.png"), outside_plot, width = 6, height = 5, dpi = 300)

plots <- list(vol_SNcd, vol_SNcv, vol_VTA, outside_plot)
plots <- Filter(Negate(is.null), plots)

if (length(plots) > 0) {
  # 2x2 layout when we have 4; sensible wrap otherwise
  combo <- patchwork::wrap_plots(plots, ncol = 2)
  ggsave("results/figs/"Volcano_AllRegions_2x2.png"), combo, width = 12, height = 10, dpi = 300)
  message("Saved combined volcano: Volcano_AllRegions_2x2.png")
} else {
  message("No volcanoes to combine.")
}

message("Done ✅  Outputs are in:")
message(" - DGE list (RDS): ", "results/dge/atlas_dge_results_final.rds")
message(" - Sig counts CSV: ", "results/dge/atlas_sig_counts_per_region.csv")
message(" - Volcano PNGs   : ", fig_dir)









library(Seurat)
library(dplyr)

# Load mapped and merged objects
mapped_obj <- readRDS("data/Final_Mapped_Object.rds")
merged_obj <- readRDS("data/Final_Analyzed_Object.rds")

# Get all ROI cells again (from saved objects)
get_cells <- function(path) {
  obj <- readRDS(path)
  colnames(obj)
}

roi_files <- c(
  "data/SN1_SNcd_Only_Object.rds",
  "data/SN1_SNcv_Only_Object.rds",
  "data/SN1_VTA_Only_Object.rds",
  "data/SN4_SNcd_Only_Object.rds",
  "data/SN4_SNcv_Only_Object.rds",
  "data/SN4_VTA_Only_Object.rds"
)

roi_cells <- unlist(lapply(roi_files, get_cells))

# Identify OUTSIDE cells
all_cells <- intersect(colnames(merged_obj), colnames(mapped_obj))
outside_cells <- setdiff(all_cells, roi_cells)

outside_obj <- subset(merged_obj, cells = outside_cells)
Idents(outside_obj) <- "condition"
DefaultAssay(outside_obj) <- "RNA"

outside_obj <- NormalizeData(outside_obj) %>% 
  FindVariableFeatures() %>%
  ScaleData()

outside_dge <- FindMarkers(
  outside_obj,
  ident.1 = "Schizophrenia",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

outside_dge$gene <- rownames(outside_dge)
outside_dge$padj <- p.adjust(outside_dge$p_val, "BH")

saveRDS(outside_dge, "outside_dge_final.rds")
message("✅ Outside DGE saved as outside_dge_final.rds")









# ============================
# Volcanoes: Top cluster per region + Outside
# ============================

# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
})

# --- 1) Inputs/outputs ---
atlas_rds   <- "atlas_dge_results_final.rds"    # <- change if needed
outside_rds <- "outside_dge_final.rds"          # <- optional; skip if missing
outdir      <- "Volcano_TopClusters"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- 2) Load DGE objects (atlas list + outside df) ---
if (exists("atlas_dge_results") && is.list(atlas_dge_results)) {
  message("Using in-memory 'atlas_dge_results'.")
} else {
  stopifnot(file.exists(atlas_rds))
  atlas_dge_results <- readRDS(atlas_rds)
  message("Loaded: ", atlas_rds)
}

outside_dge <- NULL
if (exists("outside_dge") && is.data.frame(outside_dge)) {
  message("Using in-memory 'outside_dge'.")
} else if (file.exists(outside_rds)) {
  outside_dge <- readRDS(outside_rds)
  message("Loaded: ", outside_rds)
} else {
  message("No outside_dge RDS found. Outside panel will be skipped if no object present.")
}

# --- 3) Helpers ---
is_deprecated <- function(x) grepl("^Deprecated|Codeword", x, ignore.case = TRUE)

# normalize a Seurat FindMarkers-like table to consistent columns
normalize_dge <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  df <- as.data.frame(df)
  
  # gene column
  if (!"gene" %in% names(df)) df <- tibble::rownames_to_column(df, "gene")
  
  # harmonize log2FC column
  if (!"avg_log2FC" %in% names(df)) {
    if ("avg_logFC" %in% names(df)) df$avg_log2FC <- df$avg_logFC
    else if ("log2FoldChange" %in% names(df)) df$avg_log2FC <- df$log2FoldChange
    else stop("No log2FC column found.")
  }
  
  # harmonize adjusted p
  if (!("padj" %in% names(df))) {
    if ("p_val_adj" %in% names(df)) df$padj <- df$p_val_adj
    else if ("p_val" %in% names(df)) df$padj <- p.adjust(df$p_val, method = "BH")
    else stop("No p-value column found.")
  }
  df$padj[is.na(df$padj)] <- 1
  
  # filter out codewords
  df <- df %>% filter(!is_deprecated(gene))
  
  df
}

# scoring function to pick the “best” cluster for a region
score_dge <- function(df, fdr_cut = 0.05) {
  if (is.null(df) || !nrow(df)) return(list(n_sig = 0L, score = 0))
  df2 <- normalize_dge(df)
  if (is.null(df2)) return(list(n_sig = 0L, score = 0))
  sig <- df2 %>% filter(padj < fdr_cut)
  n_sig <- nrow(sig)
  # score: number of sig genes + a tie-breaker based on sum of |log2FC|
  score <- n_sig + 0.001 * sum(abs(sig$avg_log2FC), na.rm = TRUE)
  list(n_sig = n_sig, score = score)
}

# volcano plot maker
make_volcano <- function(res, title_text, fdr_cut = 0.05, label_top = 20) {
  if (is.null(res) || !nrow(res)) return(NULL)
  df <- normalize_dge(res)
  if (is.null(df) || !nrow(df)) return(NULL)
  
  df <- df %>%
    mutate(
      neglog10p = -log10(pmax(padj, 1e-300)),
      is_sig    = padj < fdr_cut,
      dir = case_when(
        is_sig & avg_log2FC > 0 ~ "Up",
        is_sig & avg_log2FC < 0 ~ "Down",
        TRUE ~ "NotSig"
      )
    )
  
  lab <- df %>%
    filter(is_sig) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(label_top)
  
  ggplot(df, aes(x = avg_log2FC, y = neglog10p, color = dir)) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    geom_point(size = 1.6, alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                             max.overlaps = Inf, box.padding = 0.3) +
    scale_color_manual(values = c(Down = "#2b6cb0", Up = "#e53e3e", NotSig = "grey70")) +
    labs(title = title_text, x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_bw(base_size = 12)
}

# --- 4) Pick top cluster per region ---
regions <- c("SNcd", "SNcv", "VTA")
top_by_region <- tibble(region = character(), cluster = character(),
                        n_sig = integer(), score = double())

get_region_df <- function(cluster_entry, region_name) {
  # atlas_dge_results[[cluster]] is expected to be a list with $SNcd, $SNcv, $VTA entries
  if (!is.list(cluster_entry)) return(NULL)
  cand <- cluster_entry[[region_name]]
  if (is.character(cand)) return(NULL)  # "Comparison not run..." strings
  cand
}

for (reg in regions) {
  best_cluster <- NA_character_
  best_score   <- -Inf
  best_nsig    <- 0L
  
  for (cl_name in names(atlas_dge_results)) {
    df <- get_region_df(atlas_dge_results[[cl_name]], reg)
    sc <- score_dge(df)
    if (sc$score > best_score) {
      best_score <- sc$score
      best_cluster <- cl_name
      best_nsig <- sc$n_sig
    }
  }
  
  top_by_region <- bind_rows(
    top_by_region,
    tibble(region = reg, cluster = best_cluster, n_sig = best_nsig, score = best_score)
  )
}

# write a small summary CSV
readr::write_csv(top_by_region, file.path(outdir, "TopCluster_per_Region.csv"))

message("Top picks:")
print(top_by_region)

# --- 5) Build plots ---
plots <- list()

# regions (SNcd, SNcv, VTA)
for (i in seq_len(nrow(top_by_region))) {
  reg <- top_by_region$region[i]
  cl  <- top_by_region$cluster[i]
  entry <- atlas_dge_results[[cl]][[reg]]
  p <- make_volcano(entry, paste0(reg, " — Top Atlas Cluster: ", cl))
  if (!is.null(p)) {
    plots[[reg]] <- p
    ggsave(filename = file.path(outdir, paste0("Volcano_", reg, "_TopCluster_", cl, ".png")),
           plot = p, width = 6, height = 5, dpi = 300)
  }
}

# Outside ROI (optional)
if (!is.null(outside_dge) && is.data.frame(outside_dge) && nrow(outside_dge) > 0) {
  p_out <- make_volcano(outside_dge, "Outside ROI — All Cells (Schizophrenia vs Control)")
  if (!is.null(p_out)) {
    plots[["Outside"]] <- p_out
    ggsave(filename = file.path(outdir, "Volcano_Outside.png"),
           plot = p_out, width = 6, height = 5, dpi = 300)
  }
}

# --- 6) Combined 2×2 panel (SNcd | SNcv) / (VTA | Outside) ---
# if Outside missing, we still combine the 3 we have
panel <- NULL
have_out <- "Outside" %in% names(plots)

if (all(c("SNcd","SNcv","VTA") %in% names(plots))) {
  if (have_out) {
    panel <- (plots$SNcd | plots$SNcv) / (plots$VTA | plots$Outside)
    ggsave(file.path(outdir, "Volcano_TopClusters_plus_Outside.png"),
           panel, width = 12, height = 9, dpi = 300)
  } else {
    panel <- (plots$SNcd | plots$SNcv) / (plots$VTA | patchwork::plot_spacer())
    ggsave(file.path(outdir, "Volcano_TopClusters.png"),
           panel, width = 12, height = 9, dpi = 300)
  }
} else {
  message("⚠️ Could not build the 2×2 panel (missing one or more region plots).")
}

message("Done. PNGs and CSV are in: ", normalizePath(outdir))

















# --- SCRIPT TO EXTRACT DGE FOR ALL 24 CLUSTERS ---
#
# This script loads your final DGE results file and
# saves the DGE data for EVERY atlas-defined cluster
# to a single, comprehensive CSV file. This is the
# definitive check.
# ---

# 1. Load libraries
library(dplyr)
library(tibble)

# 2. Load the correct DGE results
atlas_results <- readRDS("results/dge/atlas_dge_results_final.rds")

# 3. Define ALL 24 cluster IDs
all_cluster_ids <- as.character(0:23)

# 4. Define the full cluster name map
# (From your master script)
cluster_name_map <- c(
  "0"="Astrocytes (Type 1)", "1"="Neurons (Generic)", "2"="Excitatory Neurons",
  "3"="Oligodendrocytes (Type 1)", "4"="Dopaminergic Neurons (Subtype 1)",
  "5"="GABAergic Neurons (PVALB+)", "6"="Astrocytes (Type 2)", "7"="GABAergic Neurons (VIP+)",
  "8"="Microglia", "9"="Excitatory Neurons", "10"="Dopaminergic Neurons (Subtype 2)",
  "11"="Oligodendrocyte Precursors (OPCs)", "12"="Oligodendrocytes (Type 2)",
  "13"="GABAergic Neurons", "14"="GABAergic Neurons (SST+)", "15"="Dopaminergic Neurons (Subtype 3)",
  "16"="Oligodendrocytes (Type 3)", "17"="Dopaminergic Neurons (Subtype 4)",
  "18"="Oligodendrocytes (Type 4)", "19"="Dopaminergic Neurons (Subtype 5)",
  "20"="Oligodendrocytes (Type 5)", "21"="GABAergic Neurons (PVALB+)",
  "22"="GABAergic Neurons (LAMP5+)", "23"="Endothelial Cells"
)

# 5. Create a list to hold all the data frames
all_dge_list <- list()

# 6. Loop through every cluster (0-23)
for (cluster_id in all_cluster_ids) {
  
  # Get the friendly name
  cluster_name <- ifelse(cluster_id %in% names(cluster_name_map), 
                         cluster_name_map[cluster_id], 
                         paste("Cluster", cluster_id))
  
  message(paste("Extracting data for:", cluster_name, "(Cluster", cluster_id, ")"))
  
  # Get the DGE data for this cluster
  cluster_dge_list <- atlas_results[[cluster_id]]
  
  if (is.null(cluster_dge_list)) {
    message(paste("  -> No data found for cluster", cluster_id, "in .rds file. Skipping."))
    next
  }
  
  # Extract each region (SNcd, SNcv, VTA)
  # --- SNcd ---
  df_sncd <- cluster_dge_list$SNcd
  if (!is.null(df_sncd) && nrow(df_sncd) > 0) {
    df_sncd <- df_sncd %>%
      mutate(cluster_name = cluster_name,
             region = "SNcd",
             cluster_id = cluster_id)
    all_dge_list[[paste0(cluster_id, "_SNcd")]] <- df_sncd
  }
  
  # --- SNcv ---
  df_sncv <- cluster_dge_list$SNcv
  if (!is.null(df_sncv) && nrow(df_sncv) > 0) {
    df_sncv <- df_sncv %>%
      mutate(cluster_name = cluster_name,
             region = "SNcv",
             cluster_id = cluster_id)
    all_dge_list[[paste0(cluster_id, "_SNcv")]] <- df_sncv
  }
  
  # --- VTA ---
  df_vta <- cluster_dge_list$VTA
  if (!is.null(df_vta) && nrow(df_vta) > 0) {
    df_vta <- df_vta %>%
      mutate(cluster_name = cluster_name,
             region = "VTA",
             cluster_id = cluster_id)
    all_dge_list[[paste0(cluster_id, "_VTA")]] <- df_vta
  }
}

# 7. Combine them all into one data frame
final_all_clusters_dge <- bind_rows(all_dge_list)

# 8. Save the final, comprehensive CSV
write.csv(final_all_clusters_dge, "ALL_CLUSTERS_DGE_CORRECTED.csv", row.names = FALSE)

message("\nSuccess! All 24 cluster results have been saved to:")
message("ALL_CLUSTERS_DGE_CORRECTED.csv")