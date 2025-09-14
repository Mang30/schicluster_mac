import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# Load the h5ad file
adata = sc.read_h5ad('31_svd_embedding/output/decomp/embedding_data.h5ad')

# Compute neighbors and UMAP embedding if not already computed
if 'X_umap' not in adata.obsm:
    print("Computing neighbors and UMAP embedding...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata)

# Get unique stages and cell types
stages = sorted(adata.obs['stage'].unique())
cell_types = sorted(adata.obs['celltype'].unique())
print(f"Stages: {stages}")
print(f"Cell types: {len(cell_types)} unique types")

# Create a color map for cell types
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
colors = cm.tab20(np.linspace(0, 1, len(cell_types)))
cell_type_colors = dict(zip(cell_types, colors))

# Create UMAP plots for each stage, colored by cell type
fig, axes = plt.subplots(1, len(stages), figsize=(5*len(stages), 5))
if len(stages) == 1:
    axes = [axes]
elif len(stages) == 0:
    print("No stages found in data")
    exit()

# Plot for each stage
for i, stage in enumerate(stages):
    # Create a mask for the current stage
    mask = adata.obs['stage'] == stage

    # Plot each cell type with its specific color
    for cell_type in cell_types:
        cell_mask = (adata.obs['stage'] == stage) & (adata.obs['celltype'] == cell_type)
        if np.sum(cell_mask) > 0:  # Only plot if there are cells of this type
            axes[i].scatter(
                adata.obsm['X_umap'][cell_mask, 0],
                adata.obsm['X_umap'][cell_mask, 1],
                c=[cell_type_colors[cell_type]],
                s=8, alpha=0.7, label=cell_type
            )

    axes[i].set_title(f'Stage {stage}')
    axes[i].set_xlabel('UMAP1')
    axes[i].set_ylabel('UMAP2')

# Add legend to the figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.99, 0.95), fontsize=8)

plt.tight_layout()
output_path = '/home/duxuyan/Projects/schicluster_mac/31_svd_embedding/output/decomp/stage_umap_subplots_by_celltype.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Stage-based UMAP subplots colored by cell type saved to: {output_path}")
plt.show()