"""
Core analysis modules for alarmist package

Includes LRI analysis, factorization, and GLM functionality.
"""

# LRI Analysis
# Spatial contours
from .contours import (
    motif_loading_contours,
    motif_loading_contours_from_adata,
)

# Factorization
from .factorization import (
    BPTF_AVAILABLE,
    extract_factors,
    get_top_motifs,
    process_bptf_results,
    project_cell_loadings,
    run_bptf,
    save_bptf_results,
)

# GLM Analysis
from .glm import (
    analyze_glm_results,
    check_memory_usage,
    compute_exclusion_mask,
    differential_expression,
    extract_lri_genes,
    glm_forest,
    glm_volcano,
    load_exclusion_mask,
    load_glm_results,
    prepare_cell_data_from_adata,
    prepare_cell_data_memory_efficient,
    run_poisson_glm_analysis,
    run_univariate_de_sklearn_by_celltype,
    save_de_results,
    spearman_corr_chunked,
    spearman_prefilter_genes,
)
from .lri import (
    NeighborhoodLRIAnalyzer,
    PatchLRIAnalyzer,
    SingleCellLRIAnalyzer,
    load_database_genes,
    load_database_resource,
)

# Single Cell Analysis
from .single_cell import (
    compute_motif_state_counts,
    compute_positive_motifs_per_cell,
    gmm_binarize_all_motifs,
    weighted_celltypes_by_motif,
)

__all__ = [
    # LRI
    "PatchLRIAnalyzer",
    "NeighborhoodLRIAnalyzer",
    "SingleCellLRIAnalyzer",
    "load_database_genes",
    "load_database_resource",
    # Spatial contours
    "motif_loading_contours",
    "motif_loading_contours_from_adata",
    # Factorization
    "run_bptf",
    "extract_factors",
    "get_top_motifs",
    "save_bptf_results",
    "process_bptf_results",
    "project_cell_loadings",
    "BPTF_AVAILABLE",
    # GLM
    "run_poisson_glm_analysis",
    "analyze_glm_results",
    "glm_volcano",
    "glm_forest",
    "differential_expression",
    "compute_exclusion_mask",
    "load_exclusion_mask",
    "extract_lri_genes",
    "run_univariate_de_sklearn_by_celltype",
    "prepare_cell_data_from_adata",
    "prepare_cell_data_memory_efficient",
    "check_memory_usage",
    "save_de_results",
    "load_glm_results",
    "spearman_corr_chunked",
    "spearman_prefilter_genes",
    # Single Cell
    "weighted_celltypes_by_motif",
    "gmm_binarize_all_motifs",
    "compute_motif_state_counts",
    "compute_positive_motifs_per_cell",
]
