# Spatial Transcriptomics Pipeline: Astrocytic Reprogramming in Schizophrenia

**Master of Brain and Mind Science Thesis**  
Harrison Nott | University of Sydney | 2025  
Thesis awarded High Distinction

---

## Overview

This repository contains the R analysis pipeline for my postgraduate thesis investigating astrocytic reprogramming in the substantia nigra pars compacta (SNpc) in schizophrenia using Xenium spatial transcriptomics on post-mortem human midbrain tissue.

The central finding was a reactive astrocyte signature — upregulation of **SERPINA3**, **HILPDA**, and **IGFBP4** — localised to the SNpc in schizophrenia tissue, alongside a spatial reduction in SST+ interneuron density across anatomical subregions (SNcd, SNcv, VTA). GSEA identified significant enrichment of neuroinflammatory pathways including TNF signalling and oxidative stress.

---

## Repository Structure

```
├── 01_Within_Sample_Clustering.R    # De novo clustering, spatial processing, ROI analysis
├── 02_Atlas_Integration_DGE.R       # Reference atlas mapping and cell-type-specific DGE
├── data/                            # Input data directory (not included — see note below)
│   ├── SN1/                         # Control tissue Xenium output
│   └── SN4/                         # Schizophrenia tissue Xenium output
└── results/                         # Output directory for figures and DGE tables
    ├── figs/
    └── dge/
```

---

## Analysis Pipeline

### Method 1 — De Novo Clustering (`01_Within_Sample_Clustering.R`)

Loads raw Xenium segmentation output, builds Seurat objects with spatial coordinates, and runs the full clustering and within-sample analysis pipeline:

- Quality filtering (area ≥ 100, nCount_RNA ≥ 20, nFeature_RNA ≥ 20)
- Normalisation with SCTransform
- Dimensionality reduction (PCA → UMAP)
- Leiden clustering and manual cell-type annotation into 9 populations
- Cell type proportion analysis (chi-squared)
- Global DGE between schizophrenia and control (FindMarkers, SCT)
- Spatial ROI definition for SNcd, SNcv, and VTA subregions using polygon coordinates
- Cell-type-specific DGE within each anatomical subregion
- Volcano plot visualisation (EnhancedVolcano)

**Cell types identified:** Oligodendrocytes, Astrocytes, Reactive Astrocytes, GABAergic Neurons, Dopaminergic Neurons (3 subtypes), Endothelial Cells, Generic Neurons

### Method 2 — Atlas Integration (`02_Atlas_Integration_DGE.R`)

Maps the clustered data onto the Agarwal et al. (2020) human midbrain reference atlas (GSE140231) using Seurat's anchor-based transfer learning:

- FindTransferAnchors / MapQuery for cell-type label transfer
- Validation of de novo annotations against reference predictions (contingency heatmap)
- Cell-type-specific DGE using atlas-defined populations across SNcd, SNcv, and VTA
- GSEA using clusterProfiler (KEGG pathways)
- Volcano plot generation and combined figure output

---

## Dependencies

```r
install.packages(c("Seurat", "data.table", "dplyr", "tidyr", "ggplot2",
                   "patchwork", "tibble", "pheatmap", "ggrepel", "readr"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
```

---

## Data Availability

Raw tissue data is not publicly available. Samples were obtained from the **NSW Brain Tissue Resource Centre (NSWBTRC)** under institutional ethics approval from the University of Sydney Human Research Ethics Committee. This code is provided for methodological transparency.

The reference atlas used for mapping is publicly available: [Agarwal et al. (2020), GSE140231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140231).

---

## Related Project

An independent computational replication study using a publicly available post-mortem midbrain bulk RNA-seq dataset (GSE311978, n=97) to validate the astrocytic signature identified here is available at: [link to second repo]

---

## Contact

Harrison Nott  
harrisonnott.au@gmail.com  
[ResearchGate](https://www.researchgate.net)
