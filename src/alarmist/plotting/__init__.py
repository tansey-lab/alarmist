"""
Visualization and plotting module
"""
from alarmist.plotting.colors import (
    set_celltype_colors,
    get_celltype_colors,
    clear_celltype_colors,
    _get_colors_for_plotting,
)

from alarmist.plotting.bptf_plots import (
    plot_bptf_diagnostics,
    plot_motif_activities,
    plot_factor_distributions,
    plot_factor_sparsity,
    plot_lri_factor_scatter
)

from alarmist.plotting.glm_plots import (
    volcano_plot,
    forest_plot,
    generate_volcano_plots,
    generate_forest_plots
)

from alarmist.plotting.single_cell_plots import (
    plot_motif_celltype_composition,
    plot_motif_state_counts,
    plot_positive_motifs_distribution,
    plot_motif_spatial
)

from alarmist.plotting.spatial_plots import (
    plot_cells_per_patch
)

from alarmist.plotting.motif_plots import (
    # Utility functions
    parse_lri_full,
    get_cell_type_colors,
    # Data preprocessing
    add_lri_components,
    annotate_pathways,
    # Visualization functions
    plot_lri_clustermap,
    plot_celltype_communication_by_motif,
    plot_top_lri_interactions_dot,
    plot_single_motif_lri_lollipop,
    plot_single_motif_cellpair_lollipop,
    plot_top_lri_interactions_by_pathway,
    build_master_edge_gate,
    plot_lri_networks
)

__all__ = [
    # Color management
    'set_celltype_colors',
    'get_celltype_colors',
    'clear_celltype_colors',
    # BPTF plots
    'plot_bptf_diagnostics',
    'plot_motif_activities',
    'plot_factor_distributions',
    'plot_factor_sparsity',
    'plot_lri_factor_scatter',
    'volcano_plot',
    'forest_plot',
    'generate_volcano_plots',
    'generate_forest_plots',
    'plot_motif_celltype_composition',
    'plot_motif_state_counts',
    'plot_positive_motifs_distribution',
    'plot_motif_spatial',
    # Spatial plots
    'plot_cells_per_patch',
    # Motif plots
    'parse_lri_full',
    'get_cell_type_colors',
    'add_lri_components',
    'annotate_pathways',
    'plot_lri_clustermap',
    'plot_celltype_communication_by_motif',
    'plot_top_lri_interactions_dot',
    'plot_single_motif_lri_lollipop',
    'plot_single_motif_cellpair_lollipop',
    'plot_top_lri_interactions_by_pathway',
    'build_master_edge_gate',
    'plot_lri_networks',
]
