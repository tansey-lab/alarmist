# Samplesheets

This directory contains example samplesheets for the ALARMIST Nextflow pipeline.

## Format

CSV file with columns:
- `sample_id`: Unique identifier for each sample
- `adata_path`: Path to the AnnData (.h5ad) file

## Example

```csv
sample_id,adata_path
patient_001,/data/spatial/patient_001.h5ad
patient_002,/data/spatial/patient_002.h5ad
```

## Requirements

Each AnnData file must contain:
- `adata.obsm['spatial']`: Spatial coordinates (in micrometers)
- `adata.obs['cell_type']`: Cell type annotations
- Gene expression in `adata.X` or `adata.layers['counts']`
