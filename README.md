# ALARMIST: Assessment of Ligand And Receptor Interaction Motifs in Spatial Transcriptomics

A focused pipeline for spatial transcriptomics analysis using patch-based ligand-receptor interaction (LRI) mapping and Bayesian Poisson Tensor Factorization (BPTF).

## Overview

This pipeline implements:

1. **Patch-based LRI analysis** - Spatial tissue segmentation and LRI quantification
2. **BPTF matrix factorization** - Identification of latent microenvironment programs
3. **BPTF visualization** - Comprehensive plots of motifs, networks, and spatial distributions
4. **Poisson GLM analysis** - Differential expression analysis using BPTF factors
5. **GLM results visualization** - Volcano plots and forest plots with marker gene filtering

## Workflow

### Step 1: Patch-LRI Analysis

```bash
python scripts/01_run_patch_lri_analysis.py --data-file <input.h5ad> --output-dir results/patch_lri
```

### Step 2: BPTF Matrix Factorization

```bash
python scripts/02_bptf_matrix_factorization.py --input-dir results/patch_lri --output-dir results/bptf --n-components 15
```

### Step 3: BPTF Visualization

```bash
python scripts/03_bptf_visualization.py --bptf-dir results/bptf --patch-dir results/patch_lri --data-file <input.h5ad> --output-dir results/plots/bptf_plots
```

### Step 4: Poisson GLM Analysis

```bash
python scripts/04_poisson_glm.py --data-file <input.h5ad> --results-dir results/bptf --patch-lri-dir results/patch_lri --output-dir results/glm
```

### Step 5: GLM Results Visualization

```bash
python scripts/05_glm_results.py --data-file <input.h5ad> --results-dir results/glm --output-dir results/plots/glm_plots
```

### Complete Pipeline (might not be runnable)

```bash
python run_pipeline.py --data-file <input.h5ad> --output-dir results --n-components 15
```

### Pipeline with Conda Environment Management (TODO: unify into one env)

The pipeline automatically manages conda environments for different steps:

```bash
# Default environments (bptf for BPTF step, tf2.10 for others)
python run_pipeline.py --data-file <input.h5ad> --output-dir results --n-components 15

# Custom environment names
python run_pipeline.py --data-file <input.h5ad> --output-dir results \
  --bptf-conda-env my_bptf_env --main-conda-env my_main_env

# Run without conda environment switching
python run_pipeline.py --data-file <input.h5ad> --output-dir results --no-conda
```

## Requirements

### Conda Environments (TODO: unify into one env)

The pipeline uses two conda environments:

**Main Environment (`tf2.10` by default):**

- Python 3.8+
- Standard scientific Python stack (numpy, pandas, scipy, matplotlib)
- Scanpy for single-cell analysis
- Liana for accessing LRI databases
- scikit-learn for statistical modeling
- Used for steps 1, 3, 4, and 5

**BPTF Environment (`bptf` by default):**

- Python 3.8+
- BPTF package: `https://github.com/aschein/bptf`
- Numpy, scipy (compatible versions with BPTF)
- Used for step 2 only

### Input Data Format

Expected input is AnnData (.h5ad) format with:

- `adata.obsm['spatial']`: Spatial coordinates (μm)
- `adata.obs['cell_type']`: Cell type annotations
- `adata.obs['tma_id']`: Sample/TMA identifiers

## Key Parameters

- **Patch size**: 50μm × 50μm (adjustable)
- **LRI database**: CellChatDB (via liana)
- **BPTF components**: 15 (adjust iteratively based on results)

## Output Structure

```
results/
├── patch_lri/          # Sparse patch-LRI matrices
├── bptf/              # BPTF factorization results
├── glm_results/               # GLM differential expression results
└── plots/             # All visualization outputs
    ├── bptf_plots/    # BPTF analysis plots
    └── glm_plots/     # GLM results plots
```
