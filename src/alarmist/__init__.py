"""
Alarmist: Spatial Ligand-Receptor Interaction Analysis with BPTF

Main package for spatial transcriptomics LRI analysis and matrix factorization.
"""

__version__ = "0.1.0"

# Import main classes and functions for Python API
# Column name constants
from alarmist.constants import (
    COLUMN_NAME_CELL_TYPE,
    COLUMN_NAME_GENE,
    COLUMN_NAME_LOGFC,
    COLUMN_NAME_MOTIF,
    COLUMN_NAME_PVAL,
    COLUMN_NAME_QVAL,
    COLUMN_NAME_SE,
    COLUMN_NAME_SIGNIFICANT,
)
from alarmist.core import (
    BPTF_AVAILABLE,
    NeighborhoodLRIAnalyzer,
    # LRI Analysis
    PatchLRIAnalyzer,
    SingleCellLRIAnalyzer,
    analyze_glm_results,
    check_memory_usage,
    compute_exclusion_mask,
    compute_motif_state_counts,
    compute_positive_motifs_per_cell,
    differential_expression,
    extract_factors,
    extract_lri_genes,
    get_top_motifs,
    glm_forest,
    glm_volcano,
    gmm_binarize_all_motifs,
    load_exclusion_mask,
    load_glm_results,
    motif_loading_contours,
    motif_loading_contours_from_adata,
    prepare_cell_data_from_adata,
    prepare_cell_data_memory_efficient,
    process_bptf_results,
    project_cell_loadings,
    # BPTF Factorization
    run_bptf,
    # GLM Analysis
    run_poisson_glm_analysis,
    run_univariate_de_sklearn_by_celltype,
    save_bptf_results,
    save_de_results,
    spearman_corr_chunked,
    spearman_prefilter_genes,
    # Single Cell Analysis
    weighted_celltypes_by_motif,
)
from alarmist.data.loaders import (
    load_bptf_results,
    load_cell_lri_results,
    load_patch_lri_results,
)

# Import commonly used plotting functions
from alarmist.plotting import (
    # Single cell plots
    analyze_motif_celltype_composition,
    analyze_motif_celltype_counts,
    analyze_motif_state_counts,
    clear_celltype_colors,
    forest_plot,
    get_celltype_colors,
    # Plots
    plot_bptf_diagnostics,
    plot_cells_per_patch,
    plot_factor_distributions,
    plot_lr_pair_overlap,
    plot_lri_database_overlap,
    plot_lri_factor_scatter,
    plot_lri_networks,
    plot_lri_networks_html,
    plot_motif_activities,
    plot_top_lri_interactions_dot,
    # Color management
    set_celltype_colors,
    volcano_plot,
)

__all__ = [
    # LRI analysis
    "PatchLRIAnalyzer",
    "SingleCellLRIAnalyzer",
    "NeighborhoodLRIAnalyzer",
    # Factorization
    "run_bptf",
    "extract_factors",
    "get_top_motifs",
    "save_bptf_results",
    "process_bptf_results",
    "project_cell_loadings",
    "BPTF_AVAILABLE",
    # Data I/O
    "load_patch_lri_results",
    "load_cell_lri_results",
    "load_bptf_results",
    # Color management
    "set_celltype_colors",
    "get_celltype_colors",
    "clear_celltype_colors",
    # Plotting
    "plot_bptf_diagnostics",
    "plot_motif_activities",
    "plot_factor_distributions",
    "plot_lri_database_overlap",
    "plot_lr_pair_overlap",
    "plot_lri_factor_scatter",
    "plot_lri_networks",
    "plot_lri_networks_html",
    "plot_top_lri_interactions_dot",
    "plot_cells_per_patch",
    # GLM analysis
    "run_poisson_glm_analysis",
    "analyze_glm_results",
    "glm_volcano",
    "glm_forest",
    "volcano_plot",
    "forest_plot",
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
    # Spatial contours
    "motif_loading_contours",
    "motif_loading_contours_from_adata",
    # Single Cell Analysis
    "weighted_celltypes_by_motif",
    "gmm_binarize_all_motifs",
    "compute_motif_state_counts",
    "compute_positive_motifs_per_cell",
    "analyze_motif_celltype_composition",
    "analyze_motif_celltype_counts",
    "analyze_motif_state_counts",
    # Column name constants
    "COLUMN_NAME_GENE",
    "COLUMN_NAME_LOGFC",
    "COLUMN_NAME_SE",
    "COLUMN_NAME_PVAL",
    "COLUMN_NAME_QVAL",
    "COLUMN_NAME_SIGNIFICANT",
    "COLUMN_NAME_MOTIF",
    "COLUMN_NAME_CELL_TYPE",
]
