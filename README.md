# ALARMIST

**A**ssessment of **L**igand **A**nd **R**eceptor Interaction **M**otifs **I**n **S**patial **T**ranscriptomics

ALARMIST identifies recurring ligand-receptor interaction (LRI) patterns in spatial transcriptomics data using Bayesian Poisson Tensor Factorization (BPTF).

## Overview

ALARMIST discovers **microenvironment motifs** - coordinated patterns of cell-cell communication that recur across tissue regions. The pipeline:

1. **Patchify** - Divide tissue into spatial patches and quantify LRI activity
2. **BPTF** - Factorize the patch-LRI matrix to discover latent motifs
3. **Project** - Map motifs to single-cell resolution
4. **GLM** - Analyze downstream impact via Poisson GLM differential expression
5. **Visualize** - Generate publication-ready figures

## Installation

### Using uv (recommended)

```bash
# Install uv if you don't have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and install
git clone https://github.com/tansey-lab/alarmist.git
cd alarmist
uv sync
```

### Using pip

```bash
git clone https://github.com/tansey-lab/alarmist.git
cd alarmist
pip install -e .
```

### Development Setup

```bash
# Full dev environment with pre-commit hooks
make dev
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

## CLI Usage

ALARMIST provides CLI commands for each pipeline step:

```bash
# Step 1: Patchify tissue and count LRI
alarmist-patchify --adata data.h5ad --output-dir results/patchify

# Step 2: Run BPTF factorization
alarmist-bptf --input-dir results/patchify --output-dir results/bptf --n-components 15

# Step 3: Project to single cells
alarmist-project --adata data.h5ad --bptf-dir results/bptf --output-dir results/project

# Step 4: Run GLM analysis
alarmist-glm --input-dir results/project --adata data.h5ad --output-dir results/glm

# Step 5: Generate visualizations
alarmist-visualize --glm-dir results/glm --bptf-dir results/bptf --output-dir results/plots
```

## Nextflow Pipeline

For production use, run the full nf-core pipeline:

```bash
cd nextflow

# Run with Docker
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile docker

# Run with Singularity (HPC)
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile singularity
```

### Samplesheet Format

Create a CSV file with your samples:

```csv
sample_id,adata_path
sample_A,/path/to/sample_A.h5ad
sample_B,/path/to/sample_B.h5ad
```

### Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--patch_size` | 50 | Patch size in micrometers |
| `--n_components` | 15 | Number of BPTF motifs |
| `--cell_type_column` | cell_type | Column name for cell types |
| `--resource` | cellphonedb | LRI database (cellphonedb, cellchatdb) |

See `nextflow/nextflow.config` for all parameters.

## Tutorials

See [`tutorials/GBM.ipynb`](tutorials/GBM.ipynb) for a complete walkthrough with visualizations.

## Citation

If you use ALARMIST in your research, please cite:

```
[Citation pending]
```

## License

MIT License
