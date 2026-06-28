# =============================================================================
# TH Expression Visualisation for ROI Delineation
# =============================================================================
# Project:  Astrocytic Reprogramming in Schizophrenia - Substantia Nigra
# Author:   Harrison Nott
# Degree:   Master of Brain and Mind Science, University of Sydney (2025)
#
# Description:
#   This script generates an interactive spatial plot of tyrosine hydroxylase
#   (TH) transcript expression across segmented cells in the Xenium dataset.
#   The resulting visualisation was used to manually delineate anatomical
#   subregions of the substantia nigra pars compacta (SNcd, SNcv) and VTA
#   by recording X/Y coordinates from the interactive plot window.
#   These coordinates were then used in 01_Within_Sample_Clustering.R to
#   define polygon-based regions of interest (ROIs) for spatial subsetting.
#
# Input:
#   - Baysor segmentation output CSV (baysor_output_segmentation.csv)
#     Place in ./data/SN4/ directory
#
# Output:
#   - Interactive matplotlib window showing TH expression across all cells
#   - X/Y coordinates recorded manually from this plot were used as ROI
#     polygon vertices in the R analysis pipeline
#
# Usage:
#   Run this script for each sample (SN1/SN4) by updating the file_path
#   variable. Use the interactive plot to record coordinates at anatomical
#   boundaries of interest.
#
# Dependencies:
#   pip install pandas matplotlib
# =============================================================================

import pandas as pd
import matplotlib
matplotlib.use('TkAgg')  # Interactive backend for coordinate reading
import matplotlib.pyplot as plt

# =========================================================================
# 1. LOAD DATA
# =========================================================================
# !!! Update this path to your local data directory !!!
file_path = "data/SN4/baysor_output_segmentation.csv"
print(f"Loading data from: {file_path}")

df = pd.read_csv(file_path)
print(f"Loaded {len(df)} transcripts across {df['cell'].nunique()} cells.")

# =========================================================================
# 2. PROCESS DATA
# =========================================================================
# Calculate cell centroids from transcript coordinates
cell_coords = df.groupby('cell')[['x', 'y']].mean().reset_index()

# Count TH transcripts per cell (marker of dopaminergic neurons)
th_transcripts = df[df['gene'] == 'TH']
th_counts = th_transcripts['cell'].value_counts().reset_index()
th_counts.columns = ['cell', 'TH_count']

# Merge coordinates with TH counts
plot_data = pd.merge(cell_coords, th_counts, on='cell', how='left').fillna(0)
print(f"Processing complete. {int((plot_data['TH_count'] > 0).sum())} cells express TH.")

# =========================================================================
# 3. GENERATE INTERACTIVE PLOT
# =========================================================================
# Hover over points in the interactive window to read X/Y coordinates.
# These coordinates define the polygon vertices used for ROI delineation
# in the R analysis pipeline.

fig, ax = plt.subplots(figsize=(10, 7))

scatter = ax.scatter(
    plot_data['x'],
    plot_data['y'],
    c=plot_data['TH_count'],
    cmap='hot',               # Black -> Red -> Yellow: highlights DA neurons
    s=5,
    alpha=0.8
)

ax.set_aspect('equal', adjustable='box')
ax.set_facecolor('black')
ax.invert_yaxis()             # Matches standard tissue image orientation

plt.title("SN4 (Schizophrenia): TH Expression — Use for ROI Delineation",
          color='white', pad=10)
fig.patch.set_facecolor('black')
ax.tick_params(colors='white')
ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('TH Transcript Count', color='white')
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

print("Displaying interactive plot.")
print("Hover over anatomical boundaries to read X/Y coordinates for ROI delineation.")
plt.tight_layout()
plt.show()

print("Done.")

