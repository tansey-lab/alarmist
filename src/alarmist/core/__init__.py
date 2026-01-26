"""
Core analysis modules for alarmist package

Includes LRI analysis, factorization, and GLM functionality.
"""

# LRI Analysis
from .lri import (
    PatchLRIAnalyzer,
    NeighborhoodLRIAnalyzer,
    SingleCellLRIAnalyzer
)

# Factorization
from .factorization import (
    run_bptf,
    extract_factors,
    get_top_motifs,
    save_bptf_results,
    process_bptf_results,
    project_cell_loadings,
    BPTF_AVAILABLE
)

# GLM Analysis
from .glm import (
    run_poisson_glm_analysis,
    analyze_glm_results,
    glm_volcano,
    glm_forest,
    differential_expression,
    compute_exclusion_mask,
    load_exclusion_mask,
    extract_lri_genes,
    run_univariate_de_sklearn_by_celltype,
    prepare_cell_data_from_adata,
    prepare_cell_data_memory_efficient,
    check_memory_usage,
    save_de_results,
    load_glm_results,
    spearman_corr_chunked,
    spearman_prefilter_genes
)

# Single Cell Analysis
from .single_cell import (
    weighted_celltypes_by_motif,
    gmm_binarize_all_motifs,
    compute_motif_state_counts,
    compute_positive_motifs_per_cell
)

__all__ = [
    # LRI
    'PatchLRIAnalyzer',
    'NeighborhoodLRIAnalyzer',
    'SingleCellLRIAnalyzer',
    # Factorization
    'run_bptf',
    'extract_factors',
    'get_top_motifs',
    'save_bptf_results',
    'process_bptf_results',
    'project_cell_loadings',
    'BPTF_AVAILABLE',
    # GLM
    'run_poisson_glm_analysis',
    'analyze_glm_results',
    'glm_volcano',
    'glm_forest',
    'differential_expression',
    'compute_exclusion_mask',
    'load_exclusion_mask',
    'extract_lri_genes',
    'run_univariate_de_sklearn_by_celltype',
    'prepare_cell_data_from_adata',
    'prepare_cell_data_memory_efficient',
    'check_memory_usage',
    'save_de_results',
    'load_glm_results',
    'spearman_corr_chunked',
    'spearman_prefilter_genes',
    # Single Cell
    'weighted_celltypes_by_motif',
    'gmm_binarize_all_motifs',
    'compute_motif_state_counts',
    'compute_positive_motifs_per_cell',
]
