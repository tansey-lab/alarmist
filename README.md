# ALARMIST

**A**ssessment of **L**igand **A**nd **R**eceptor Interaction **M**otifs **I**n **S**patial **T**ranscriptomics

ALARMIST identifies recurring ligand-receptor interaction (LRI) patterns in spatial transcriptomics data using Bayesian Poisson Tensor Factorization (BPTF).

## Overview

ALARMIST discovers **microenvironment motifs** - coordinated patterns of cell-cell communication that recur across tissue regions. The pipeline:

1. **Patchifies** tissue into spatial units and quantifies LRI activity per patch
2. **Factorizes** the patch-LRI matrix using BPTF to discover latent motifs
3. **Projects** motifs to single-cell resolution
4. **Analyzes** downstream impact via Poisson GLM differential expression

## Installation

```bash
# Clone the repository
git clone https://github.com/tansey-lab/alarmist.git
cd alarmist

# Install dependencies
pip install -r requirements.txt

# Install ALARMIST
pip install -e .
```

## Quick Start

```python
import alarmist as al
import scanpy as sc

# Load your spatial transcriptomics data
adata = sc.read_h5ad('your_data.h5ad')
# Required: adata.obs['cell_type'], adata.obsm['spatial']

# === Step 1: Patchify & Count LRI ===
analyzer = al.PatchLRIAnalyzer(
    patch_size=50.0,  # micrometers
    resource_name='cellchatdb'
)
results = analyzer.run_patchify(adata, output_dir='results/')

# === Step 2: BPTF Factorization ===
model = al.run_bptf(results['patch_lri_matrix'], n_components=15)
bptf_results = al.process_bptf_results(model, results, output_dir='results/bptf')

# === Step 3: Single-Cell Projection ===
cell_analyzer = al.NeighborhoodLRIAnalyzer(neighborhood_size=50.0)
cell_results = cell_analyzer.run_neighborhood(
    adata,
    required_columns=results['column_names'],
    output_dir='results/single_cell'
)
cell_loadings = al.project_cell_loadings(
    model=model,
    cell_lri_matrix=cell_results['cell_lri_matrix']
)

# === Step 4: Motif Calling (GMM Binarization) ===
gmm_summary = al.gmm_binarize_all_motifs(cell_loadings, adata)
# Now adata.obs contains motif_{k}_state ('positive'/'negative')

# === Step 5: Downstream Analysis ===
glm_results = al.run_poisson_glm_analysis(
    cell_loadings=cell_loadings,
    adata=adata,
    lri_column_names=results['column_name'],
    prefilter_spearman=True  # faster
)
```

## Input Requirements

AnnData (`.h5ad`) with:
- `adata.obsm['spatial']`: Spatial coordinates (microns)
- `adata.obs['cell_type']`: Cell type annotations
- Raw counts in `adata.X` or `adata.layers['counts']`

## Multi-Sample Support

ALARMIST supports multiple samples via:

```python
# Option 1: Merged AnnData with sample column
results = analyzer.run_patchify(
    adata,
    multi_sample=True,
    sample_column='sample_id'
)

# Option 2: Dictionary of AnnData objects
results = analyzer.run_patchify({
    'sample_A': adata_a,
    'sample_B': adata_b
})
```

## Tutorials

See [`tutorials/GBM.ipynb`](tutorials/GBM.ipynb) for a complete walkthrough with visualizations.

## Citation

If you use ALARMIST in your research, please cite:

```
[Citation pending]
```

## License

MIT License
