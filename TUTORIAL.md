# ALARMIST Tutorial

ALARMIST (**A**nalysis of **L**igand-receptor interactions **A**cross **R**egions using **M**atrix factorization for **I**dentifying **S**patial **T**ranscriptomic motifs) is a Python package for spatial ligand-receptor interaction (LRI) analysis with BPTF matrix factorization.

## Installation

```bash
cd /path/to/alarmist
pip install -e .

# Required: Install BPTF
pip install git+https://github.com/aschein/bptf.git
```

## Package Overview

```python
import alarmist as al
```

All main functions are accessible from the top-level import:

| Category | Functions |
|----------|-----------|
| LRI Analysis | `PatchLRIAnalyzer`, `NeighborhoodLRIAnalyzer` |
| Factorization | `run_bptf`, `extract_factors`, `save_bptf_results`, `project_cell_loadings` |
| Data I/O | `load_patch_lri_results`, `load_cell_lri_results`, `load_bptf_results` |
| GLM Analysis | `run_poisson_glm_analysis`, `glm_volcano`, `glm_forest` |
| Marker Genes | `compute_exclusion_mask`, `load_exclusion_mask` |
| Single Cell | `gmm_binarize_all_motifs`, `weighted_celltypes_by_motif` |

### Saving Convention

All analysis functions follow a consistent pattern:
- **`output_dir=None`** (default): Results kept in memory only, nothing saved
- **`output_dir="path/to/dir"`**: Results saved to disk AND returned in memory

This allows flexible workflows—quick exploration without saving, or persistent checkpoints for long analyses.

---

## Single-Sample Workflow

A complete analysis pipeline from raw data to visualization.

```python
import alarmist as al
import anndata
import numpy as np
import pandas as pd

# ============================================================
# 1. Load Data
# ============================================================
adata = anndata.read_h5ad("data/sample.h5ad")
print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

# ============================================================
# 2. Patch-Based LRI Analysis
# ============================================================
patch_analyzer = al.PatchLRIAnalyzer(
    patch_size=50.0,              # µm
    resource_name='cellchatdb',
    cell_type_column='cell_type'
)

# Run analysis (save to disk for later use)
patch_results = patch_analyzer.run_patchify(
    adata,
    output_dir="results/patch_lri"  # omit to skip saving
)

# Results available immediately in memory
patch_lri_matrix = patch_results['patch_lri_matrix']
column_names = patch_results['column_names']
print(f"Patch LRI matrix: {patch_lri_matrix.shape}")

# ============================================================
# 3. BPTF Factorization
# ============================================================
np.random.seed(0)

model = al.run_bptf(
    patch_lri_matrix,
    n_components=15,
    max_iter=10000,
    random_state=0
)

patch_loadings, lri_factors = al.extract_factors(model)
print(f"Patch loadings: {patch_loadings.shape}")
print(f"LRI factors: {lri_factors.shape}")

# Save BPTF results
al.save_bptf_results(
    model=model,
    patch_loadings=patch_loadings,
    lri_factors=lri_factors,
    column_names=column_names,
    patch_metadata_df=patch_results['patch_tma_df'],
    patch_lri_matrix=patch_lri_matrix,
    output_dir="results/bptf"  # omit to skip saving
)

# ============================================================
# 4. Cell-Level LRI Analysis
# ============================================================
cell_analyzer = al.NeighborhoodLRIAnalyzer(
    neighborhood_size=50.0,       # µm
    resource_name='cellchatdb'
)

# Align columns to patch matrix for projection
cell_results = cell_analyzer.run_neighborhood(
    adata,
    output_dir="results/cell_lri",
    required_columns=column_names  # ensures column alignment
)

cell_lri_matrix = cell_results['cell_lri_matrix']
print(f"Cell LRI matrix: {cell_lri_matrix.shape}")

# ============================================================
# 5. Project Cell Loadings
# ============================================================
cell_loadings = al.project_cell_loadings(
    model=model,
    cell_lri_matrix=cell_lri_matrix,
    model_lri_columns=column_names,
    cell_lri_columns=cell_results['column_names'],
    max_iter=200,
    chunk_size=50000,
    verbose=True
)

print(f"Cell loadings: {cell_loadings.shape}")
np.save("results/cell_loadings.npy", cell_loadings)

# ============================================================
# 6. GMM Binarization (ON/OFF Calling)
# ============================================================
gmm_summary = al.gmm_binarize_all_motifs(
    cell_loadings=cell_loadings,
    adata=adata,
    random_state=42
)

print(gmm_summary)
# Results now in adata.obs: motif_0_state, motif_0_posprob, etc.

# ============================================================
# 7. GLM Differential Expression
# ============================================================
glm_results = al.run_poisson_glm_analysis(
    cell_loadings=cell_loadings,
    adata=adata,
    lri_column_names=pd.Series(column_names),
    output_dir="results/glm",  # omit to skip saving
    random_state=42
)

print(f"GLM completed for {len(glm_results)} motif-celltype combinations")

# ============================================================
# 8. Compute Marker Gene Exclusion Mask
# ============================================================
# This can be slow—save for reuse
cell_types, genes, exclusion_mask = al.compute_exclusion_mask(
    adata,
    marker_lfc=1.0,
    marker_pvalue=1e-5,
    output_dir="results/markers"  # omit to skip saving
)

# ============================================================
# 9. Visualization
# ============================================================
import matplotlib.pyplot as plt

# Volcano plots for specific motifs
fig = al.glm_volcano(
    adata=adata,
    de_results=glm_results,
    cell_types=cell_types,
    all_genes=genes,
    exclusion_mask=exclusion_mask,
    motif_id=[0, 1, 2]
)
plt.show()

# Generate all volcano plots (multi-page PDF)
al.glm_volcano(
    adata=adata,
    de_results=glm_results,
    cell_types=cell_types,
    all_genes=genes,
    exclusion_mask=exclusion_mask,
    output_dir="results/volcano"  # required when motif_id=None
)

# Forest plots
al.glm_forest(
    adata=adata,
    de_results=glm_results,
    cell_types=cell_types,
    all_genes=genes,
    exclusion_mask=exclusion_mask,
    output_dir="results/forest"
)

print("Single-sample analysis complete!")
```

---

## Multi-Sample Workflow

When analyzing multiple samples, ALARMIST automatically:
- Uses the **intersection of genes** across all samples
- Creates a **unified column structure** (cell type pairs × LRI pairs × signaling modes)
- Assigns **globally unique** patch/cell indices
- Fits **GMM on combined data** for consistent thresholds

```python
import alarmist as al
import anndata
import numpy as np
import pandas as pd

# ============================================================
# 1. Load Multiple Samples
# ============================================================
adata_dict = {
    'sample_A': anndata.read_h5ad("data/sample_A.h5ad"),
    'sample_B': anndata.read_h5ad("data/sample_B.h5ad"),
    'sample_C': anndata.read_h5ad("data/sample_C.h5ad"),
}

for name, ad in adata_dict.items():
    print(f"{name}: {ad.n_obs} cells, {ad.n_vars} genes")

# ============================================================
# 2. Patch-Based LRI Analysis (Multi-Sample)
# ============================================================
patch_analyzer = al.PatchLRIAnalyzer(
    patch_size=50.0,
    resource_name='cellchatdb',
    cell_type_column='cell_type'
)

# Pass dict of AnnData objects
patch_results = patch_analyzer.run_patchify(
    adata_dict,  # dict, not single AnnData
    output_dir="results/patch_lri"
)

patch_lri_matrix = patch_results['patch_lri_matrix']
column_names = patch_results['column_names']
print(f"Combined patch matrix: {patch_lri_matrix.shape}")

# ============================================================
# 3. BPTF Factorization
# ============================================================
np.random.seed(0)

model = al.run_bptf(
    patch_lri_matrix,
    n_components=15,
    max_iter=10000,
    random_state=0
)

patch_loadings, lri_factors = al.extract_factors(model)

al.save_bptf_results(
    model=model,
    patch_loadings=patch_loadings,
    lri_factors=lri_factors,
    column_names=column_names,
    patch_metadata_df=patch_results['patch_metadata_df'],
    patch_lri_matrix=patch_lri_matrix,
    output_dir="results/bptf"
)

# ============================================================
# 4. Cell-Level LRI Analysis (Multi-Sample)
# ============================================================
cell_analyzer = al.NeighborhoodLRIAnalyzer(
    neighborhood_size=50.0,
    resource_name='cellchatdb'
)

cell_results = cell_analyzer.run_neighborhood(
    adata_dict,  # dict of AnnData
    output_dir="results/cell_lri",
    required_columns=column_names
)

cell_lri_matrix = cell_results['cell_lri_matrix']
print(f"Combined cell matrix: {cell_lri_matrix.shape}")

# ============================================================
# 5. Project Cell Loadings
# ============================================================
cell_loadings = al.project_cell_loadings(
    model=model,
    cell_lri_matrix=cell_lri_matrix,
    model_lri_columns=column_names,
    cell_lri_columns=cell_results['column_names'],
    max_iter=200,
    chunk_size=50000
)

np.save("results/cell_loadings.npy", cell_loadings)

# Split by sample if needed
cell_meta = cell_results['cell_metadata_df']
loadings_dict = {}
for sample_id in adata_dict.keys():
    mask = cell_meta['sample_id'] == sample_id
    loadings_dict[sample_id] = cell_loadings[mask]
    print(f"{sample_id}: {loadings_dict[sample_id].shape}")

# ============================================================
# 6. GMM Binarization (Multi-Sample)
# ============================================================
# GMM fitted on combined data for consistent thresholds
gmm_summary = al.gmm_binarize_all_motifs(
    cell_loadings=loadings_dict,  # dict of arrays
    adata=adata_dict,             # dict of AnnData
    random_state=42
)

print(gmm_summary)
# Results added to each adata in adata_dict

# ============================================================
# 7. Downstream Analysis Per Sample
# ============================================================
# After GMM, analyze each sample separately or combined

# Example: GLM on combined data
combined_adata = anndata.concat(
    adata_dict.values(),
    keys=adata_dict.keys(),
    label='sample'
)

glm_results = al.run_poisson_glm_analysis(
    cell_loadings=cell_loadings,
    adata=combined_adata,
    lri_column_names=pd.Series(column_names),
    output_dir="results/glm"
)

# Compute exclusion mask on combined data
cell_types, genes, exclusion_mask = al.compute_exclusion_mask(
    combined_adata,
    output_dir="results/markers"
)

# Visualization
al.glm_volcano(
    adata=combined_adata,
    de_results=glm_results,
    cell_types=cell_types,
    all_genes=genes,
    exclusion_mask=exclusion_mask,
    output_dir="results/volcano"
)

print("Multi-sample analysis complete!")
```

---

## Loading Saved Results

Resume analysis from saved checkpoints:

```python
import alarmist as al
import numpy as np

# Load patch LRI results
patch_results = al.load_patch_lri_results("results/patch_lri")
patch_lri_matrix = patch_results['patch_lri_matrix']
column_names = patch_results['column_names']

# Load BPTF results
bptf_results = al.load_bptf_results("results/bptf")
patch_loadings = bptf_results['patch_loadings']
lri_factors = bptf_results['lri_factors']
lri_motifs = bptf_results['lri_motifs']  # DataFrame with factor_norm

# Load cell LRI results
cell_results = al.load_cell_lri_results("results/cell_lri")
cell_lri_matrix = cell_results['cell_lri_matrix']

# Load cell loadings
cell_loadings = np.load("results/cell_loadings.npy")

# Load exclusion mask
cell_types, genes, exclusion_mask = al.load_exclusion_mask(
    "results/markers/exclusion_matrix.csv"
)
```

---

## Key Concepts

### Signaling Modes

ALARMIST distinguishes three signaling modes:
- **Autocrine**: Same cell expresses both ligand and receptor
- **Paracrine**: Different cells (secreted signaling), excludes autocrine
- **Juxtacrine**: Cell-cell contact signaling, excludes autocrine

### Ligand/Receptor Complexes

Both ligand and receptor complexes are supported:
- `"IL12A_IL12B"` or `"IL12A,IL12B"` → both genes required (AND logic)
- A cell only counts as expressing the complex if ALL genes are expressed

### Column Alignment

When projecting cell loadings from a patch model:
- Use `required_columns` in `run_neighborhood()` to pre-align columns
- Or let `project_cell_loadings()` handle alignment automatically

### Memory Efficiency

For large datasets:
- `project_cell_loadings()` processes cells in chunks (default 50,000)
- `compute_exclusion_mask()` subsamples cells for marker detection

---

## Visualization Functions

### BPTF Diagnostics

```python
from alarmist.plotting import plot_bptf_diagnostics

plot_bptf_diagnostics(
    patch_loadings=patch_loadings,
    lri_factors=lri_factors,
    output_dir="results/bptf_plots"  # or omit to display
)
```

### Single Cell Plots

```python
from alarmist.plotting import (
    plot_motif_celltype_composition,
    plot_motif_state_counts,
    plot_motif_spatial
)

# Cell type composition per motif
tidy_df = al.weighted_celltypes_by_motif(
    cell_loadings, adata.obs[['cell_type']],
    normalize=True, top_n_per_motif=8
)
plot_motif_celltype_composition(tidy_df, save_path="composition.pdf")

# ON/OFF counts
counts_df = al.compute_motif_state_counts(adata)
plot_motif_state_counts(counts_df, save_path="state_counts.pdf")

# Spatial distribution
plot_motif_spatial(adata, motif_id=0, save_path="motif_0_spatial.pdf")
```

### Motif Visualization

```python
from alarmist.plotting import (
    plot_top_lri_interactions_dot,
    plot_celltype_communication_by_motif,
    add_lri_components,
    annotate_pathways
)

# Prepare data
lri_motifs = bptf_results['lri_motifs']
lri_motifs = add_lri_components(lri_motifs)

# Top LRI interactions
unique_ct = sorted(adata.obs['cell_type'].unique())
plot_top_lri_interactions_dot(
    lri_motifs, unique_ct,
    use_normalized=True,
    top_n=35,
    save_path="top_lri.pdf"
)

# Communication heatmaps
plot_celltype_communication_by_motif(
    lri_factors, column_names,
    save_path="comm_heatmap.pdf"
)
```

---

## Tips

1. **Check BPTF availability**:
   ```python
   print(f"BPTF available: {al.BPTF_AVAILABLE}")
   ```

2. **Reproducibility**: Set random seeds before BPTF and GMM:
   ```python
   np.random.seed(0)
   model = al.run_bptf(...)
   ```

3. **Save intermediate results**: For large datasets, save after each major step to avoid recomputation.

4. **Exclusion mask reuse**: Computing marker genes is slow. Save once, load for subsequent visualizations:
   ```python
   # First run
   ct, genes, mask = al.compute_exclusion_mask(adata, output_dir="markers")

   # Later
   ct, genes, mask = al.load_exclusion_mask("markers/exclusion_matrix.csv")
   ```

5. **Visualization modes**:
   - `motif_id=[0, 1, 2]`: View specific motifs (returns figure)
   - `motif_id=None` + `output_dir`: Generate all motifs (multi-page PDF)
