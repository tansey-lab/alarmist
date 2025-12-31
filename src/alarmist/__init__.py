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
    project_cell_loadings,
    # GLM Analysis
    run_poisson_glm_analysis,
    analyze_glm_results,
    glm_volcano,
    glm_forest,
    # Single Cell Analysis
    weighted_celltypes_by_motif,
    gmm_binarize_all_motifs,
    compute_motif_state_counts,
    compute_positive_motifs_per_cell,
)

from alarmist.data.loaders import load_patch_lri_results, load_cell_lri_results, load_bptf_results

# Import commonly used plotting functions
from alarmist.plotting import (
    plot_bptf_diagnostics,
    plot_motif_activities,
    plot_factor_distributions,
    volcano_plot,
    forest_plot,
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
    'project_cell_loadings',

    # Data I/O
    'load_patch_lri_results',
    'load_cell_lri_results',
    'load_bptf_results',

    # Plotting
    'plot_bptf_diagnostics',
    'plot_motif_activities',
    'plot_factor_distributions',

    # GLM analysis
    'run_poisson_glm_analysis',
    'analyze_glm_results',
    'glm_volcano',
    'glm_forest',
    'volcano_plot',
    'forest_plot',

    # Single Cell Analysis
    'weighted_celltypes_by_motif',
    'gmm_binarize_all_motifs',
    'compute_motif_state_counts',
    'compute_positive_motifs_per_cell',
]
