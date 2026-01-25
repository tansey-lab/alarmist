"""
Alarmist: Spatial Ligand-Receptor Interaction Analysis with BPTF

Main package for spatial transcriptomics LRI analysis and matrix factorization.
"""

__version__ = "0.1.0"

# Import main classes and functions for Python API
from alarmist.core import (
    # LRI Analysis
    PatchLRIAnalyzer,
    SingleCellLRIAnalyzer,
    NeighborhoodLRIAnalyzer,
    # BPTF Factorization
    run_bptf,
    extract_factors,
    get_top_motifs,
    save_bptf_results,
    process_bptf_results,
    project_cell_loadings,
    BPTF_AVAILABLE,
    # GLM Analysis
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
    # Single Cell Analysis
    weighted_celltypes_by_motif,
    gmm_binarize_all_motifs,
    compute_motif_state_counts,
    compute_positive_motifs_per_cell,
)

from alarmist.data.loaders import load_patch_lri_results, load_cell_lri_results, load_bptf_results

# Import commonly used plotting functions
from alarmist.plotting import (
    # Color management
    set_celltype_colors,
    get_celltype_colors,
    clear_celltype_colors,
    # Plots
    plot_bptf_diagnostics,
    plot_motif_activities,
    plot_factor_distributions,
    plot_lri_factor_scatter,
    volcano_plot,
    forest_plot,
    plot_cells_per_patch,
)

__all__ = [
    # LRI analysis
    'PatchLRIAnalyzer',
    'SingleCellLRIAnalyzer',
    'NeighborhoodLRIAnalyzer',

    # Factorization
    'run_bptf',
    'extract_factors',
    'get_top_motifs',
    'save_bptf_results',
    'process_bptf_results',
    'project_cell_loadings',
    'BPTF_AVAILABLE',

    # Data I/O
    'load_patch_lri_results',
    'load_cell_lri_results',
    'load_bptf_results',

    # Color management
    'set_celltype_colors',
    'get_celltype_colors',
    'clear_celltype_colors',

    # Plotting
    'plot_bptf_diagnostics',
    'plot_motif_activities',
    'plot_factor_distributions',
    'plot_lri_factor_scatter',
    'plot_cells_per_patch',

    # GLM analysis
    'run_poisson_glm_analysis',
    'analyze_glm_results',
    'glm_volcano',
    'glm_forest',
    'volcano_plot',
    'forest_plot',
    'differential_expression',
    'compute_exclusion_mask',
    'load_exclusion_mask',
    'extract_lri_genes',
    'run_univariate_de_sklearn_by_celltype',
    'prepare_cell_data_from_adata',
    'prepare_cell_data_memory_efficient',
    'check_memory_usage',
    'save_de_results',

    # Single Cell Analysis
    'weighted_celltypes_by_motif',
    'gmm_binarize_all_motifs',
    'compute_motif_state_counts',
    'compute_positive_motifs_per_cell',
]
