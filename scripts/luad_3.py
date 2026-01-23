#!/usr/bin/env python
"""
GLM Differential Expression Analysis for ALARMIST

Analyzes ligand-receptor interaction motifs using Poisson GLM
to identify genes associated with each motif in different cell types.
"""

import numpy as np
import pandas as pd
import anndata as ad
import alarmist as al

# =============================================================================
# Configuration
# =============================================================================

DATA_DIR = "data/linghua"
RESULTS_DIR = "results/AIS_LUAD"

FILES = {
    "P17_AIS": f"{DATA_DIR}/P17_AIS_Xenium.h5ad",
    "P17_LUAD": f"{DATA_DIR}/P17_LUAD_Xenium.h5ad",
    "P21_AIS": f"{DATA_DIR}/P21_AIS_Xenium.h5ad",
    "P21_LUAD": f"{DATA_DIR}/P21_LUAD_Xenium.h5ad",
}

CELL_TYPE_COL = "annotation_coarse"  # Column to rename to 'cell_type'
RANDOM_STATE = 42


# =============================================================================
# Load and preprocess data
# =============================================================================

def load_adata_dict(files: dict, cell_type_col: str = "annotation_coarse") -> dict:
    """
    Load multiple h5ad files and preprocess cell type annotations.
    
    Parameters
    ----------
    files : dict
        Mapping of sample names to file paths
    cell_type_col : str
        Column name to rename to 'cell_type'
    
    Returns
    -------
    dict
        Dictionary of preprocessed AnnData objects
    """
    adata_dict = {}
    
    for name, path in files.items():
        adata = ad.read_h5ad(path)
        
        # Rename cell type column
        adata.obs.rename(columns={cell_type_col: "cell_type"}, inplace=True)
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
        
        # Remove cells with missing cell type
        mask = adata.obs["cell_type"].notna()
        n_removed = (~mask).sum()
        print(f"{name}: {adata.n_obs} -> {mask.sum()} cells ({n_removed} NaN removed)")
        
        adata_dict[name] = adata[mask].copy()
    
    return adata_dict

import pandas as pd
from pathlib import Path

def load_glm_results(results_dir: str) -> dict:
    """
    Load GLM DE results from saved CSV files.
    
    Parameters
    ----------
    results_dir : str
        Directory containing *_de_results.csv files
    
    Returns
    -------
    dict
        Dictionary mapping motif-celltype keys to DataFrames
    """
    results_path = Path(results_dir)
    glm_results = {}
    
    for csv_file in results_path.glob("motif_*_celltype_*_de_results.csv"):
        # Extract key from filename: motif_0_celltype_Tcell_de_results.csv -> motif_0_celltype_Tcell
        key = csv_file.stem.replace("_de_results", "")
        glm_results[key] = pd.read_csv(csv_file)
    
    print(f"Loaded {len(glm_results)} motif-celltype combinations")
    return glm_results


def main():
    # Load individual samples
    print("=" * 60)
    print("Loading data...")
    print("=" * 60)
    adata_dict = load_adata_dict(FILES, cell_type_col=CELL_TYPE_COL)
    
    # Combine samples
    print("\nCombining samples...")
    adata_combined = ad.concat(adata_dict, label="sample", join="outer")
    print(f"Combined shape: {adata_combined.shape}")
    print(adata_combined.obs["sample"].value_counts())
    
    # # Load BPTF results
    # print("\n" + "=" * 60)
    # print("Loading BPTF results...")
    # print("=" * 60)
    # lri_results = al.load_patch_lri_results(RESULTS_DIR)
    # cell_loadings = np.load(f"{RESULTS_DIR}/cell_loadings.npy")
    # print(f"Cell loadings shape: {cell_loadings.shape}")
    
    # # Verify alignment
    # assert cell_loadings.shape[0] == adata_combined.shape[0], (
    #     f"Mismatch: cell_loadings has {cell_loadings.shape[0]} cells, "
    #     f"adata has {adata_combined.shape[0]} cells"
    # )
    
    # # Run GLM analysis
    # print("\n" + "=" * 60)
    # print("Running Poisson GLM differential expression analysis...")
    # print("=" * 60)
    
    # glm_results = al.run_poisson_glm_analysis(
    #     cell_loadings=cell_loadings,
    #     adata=adata_combined,
    #     lri_column_names=pd.Series(lri_results["column_names"]),
    #     output_dir=f"{RESULTS_DIR}/impact",
    #     count_layer="X",
    #     random_state=RANDOM_STATE,
    # )

    # print("\nGLM analysis complete.")
    # print(f"Results saved to: {RESULTS_DIR}/impact")

    glm_results = load_glm_results(f"{RESULTS_DIR}/impact")

    # ========== Step 1: Compute Marker Gene Exclusion Mask ==========
    # print("1. Computing marker gene exclusion mask...")

    # # First time: compute and optionally save
    # cell_types, genes, exclusion_mask = al.compute_exclusion_mask(
    #     adata_combined,
    #     marker_lfc=1.0,
    #     marker_pvalue=1e-5,
    #     marker_subsample=50000,
    #     output_dir=f"{RESULTS_DIR}/markers"  # Save for reuse
    # )
    # print(f"Computed exclusion mask for {len(cell_types)} cell types, {len(genes)} genes")

    cell_types, genes, exclusion_mask = al.load_exclusion_mask(f"{RESULTS_DIR}/markers/exclusion_matrix.csv")

    # ========== Step 2: Volcano Plots ==========
    print("\n2. Creating volcano plots...")
    al.glm_volcano(
        adata=adata_combined,
        de_results=glm_results,
        cell_types=cell_types,
        all_genes=genes,
        exclusion_mask=exclusion_mask,
        fdr_threshold=0.05,
        lfc_threshold=0.2,
        output_dir = f"{RESULTS_DIR}/impact"
    )

    # ========== Step 3: Forest Plots ==========
    print("\n3. Creating forest plots...")
    al.glm_forest(
        adata=adata_combined,
        de_results=glm_results,
        cell_types=cell_types,
        all_genes=genes,
        exclusion_mask=exclusion_mask,
        output_dir=f"{RESULTS_DIR}/impact"  # Required when motif_id=None
    )
    print("\nAll analyses complete.")

if __name__ == "__main__":
    main()