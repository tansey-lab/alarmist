"""
Visualization and plotting module
"""

from alarmist.plotting.bptf_plots import (
    plot_bptf_diagnostics,
    plot_factor_distributions,
    plot_factor_sparsity,
    plot_lri_factor_scatter,
    plot_motif_activities,
)
from alarmist.plotting.colors import (
    clear_celltype_colors,
    get_celltype_colors,
    set_celltype_colors,
)
from alarmist.plotting.glm_plots import (
    forest_plot,
    generate_forest_plots,
    generate_volcano_plots,
    volcano_plot,
)
from alarmist.plotting.lri_overlap import plot_lri_database_overlap
from alarmist.plotting.motif_plots import (
    # Data preprocessing
    add_lri_components,
    annotate_pathways,
    build_master_edge_gate,
    get_cell_type_colors,
    # Utility functions
    parse_lri_full,
    plot_celltype_communication_by_motif,
    # Visualization functions
    plot_lri_clustermap,
    plot_lri_networks,
    plot_single_motif_cellpair_lollipop,
    plot_single_motif_lri_lollipop,
    plot_top_lri_interactions_by_pathway,
    plot_top_lri_interactions_dot,
)
from alarmist.plotting.single_cell_plots import (
    analyze_motif_celltype_composition,
    analyze_motif_celltype_counts,
    analyze_motif_state_counts,
    plot_motif_celltype_composition,
    plot_motif_spatial,
    plot_motif_state_counts,
    plot_positive_motifs_distribution,
)
from alarmist.plotting.spatial_plots import plot_cells_per_patch

__all__ = [
    # Color management
    "set_celltype_colors",
    "get_celltype_colors",
    "clear_celltype_colors",
    # BPTF plots
    "plot_bptf_diagnostics",
    "plot_motif_activities",
    "plot_factor_distributions",
    "plot_factor_sparsity",
    "plot_lri_factor_scatter",
    "volcano_plot",
    "forest_plot",
    "generate_volcano_plots",
    "generate_forest_plots",
    "plot_motif_celltype_composition",
    "plot_motif_state_counts",
    "plot_positive_motifs_distribution",
    "plot_motif_spatial",
    "analyze_motif_celltype_composition",
    "analyze_motif_celltype_counts",
    "analyze_motif_state_counts",
    # Spatial plots
    "plot_cells_per_patch",
    # Motif plots
    "parse_lri_full",
    "get_cell_type_colors",
    "add_lri_components",
    "annotate_pathways",
    "plot_lri_clustermap",
    "plot_celltype_communication_by_motif",
    "plot_top_lri_interactions_dot",
    "plot_single_motif_lri_lollipop",
    "plot_single_motif_cellpair_lollipop",
    "plot_top_lri_interactions_by_pathway",
    "build_master_edge_gate",
    "plot_lri_networks",
    "plot_lri_database_overlap",
]
