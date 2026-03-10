"""
CLI command: alarmist-visualize

Generate visualizations from ALARMIST pipeline results.
"""

import argparse
import sys

from alarmist.cli import common, log_config
from alarmist.constants import (
    COLUMN_NAME_GENE,
    COLUMN_NAME_LOGFC,
    COLUMN_NAME_NEG_LOG10_Q,
    COLUMN_NAME_QVAL,
)


def get_parser():
    """Create argument parser for visualize command"""
    parser = argparse.ArgumentParser(
        prog="alarmist-visualize",
        description="Generate visualizations from ALARMIST pipeline results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - generate all plots
  alarmist-visualize --glm-dir results/glm --bptf-dir results/bptf --output-dir results/plots

  # Generate only specific plot types
  alarmist-visualize --glm-dir results/glm --bptf-dir results/bptf --output-dir results/plots \\
      --plot-types volcano forest heatmap
        """,
    )

    # Input arguments
    parser.add_argument(
        "--glm-dir",
        "-g",
        type=str,
        required=True,
        help="Directory containing GLM results (from alarmist-glm)",
    )
    parser.add_argument(
        "--bptf-dir",
        "-b",
        type=str,
        required=True,
        help="Directory containing BPTF results (from alarmist-bptf)",
    )
    parser.add_argument(
        "--project-dir",
        "-p",
        type=str,
        default=None,
        help="Directory containing project results with projected_adata.h5ad (from alarmist-project)",
    )
    parser.add_argument(
        "--patchify-dir",
        type=str,
        default=None,
        help="Directory containing patchify results (from alarmist-patchify)",
    )

    common.add_output_arguments(parser)

    # Plot options
    parser.add_argument(
        "--plot-types",
        type=str,
        nargs="+",
        default=["volcano", "forest", "heatmap", "motif_summary", "spatial"],
        choices=["volcano", "forest", "heatmap", "motif_summary", "spatial", "all"],
        help="Types of plots to generate (default: volcano forest heatmap motif_summary spatial)",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="png",
        choices=["png", "pdf", "svg"],
        help="Output format for plots (default: png)",
    )
    parser.add_argument(
        "--dpi", type=int, default=150, help="DPI for raster outputs (default: 150)"
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold for highlighting (default: 0.05)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top features to show in summary plots (default: 20)",
    )

    log_config.add_logging_args(parser)

    return parser


def main():
    """Main entry point for visualize command"""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging
    logger = log_config.configure_logging(args)

    # Import heavy dependencies only after argument parsing
    logger.info("Loading dependencies...")
    from pathlib import Path

    import matplotlib
    import pandas as pd
    import yaml

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt

    import alarmist as al

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve 'all' plot type
    plot_types = args.plot_types
    if "all" in plot_types:
        plot_types = ["volcano", "forest", "heatmap", "motif_summary", "spatial"]

    # Load GLM results
    logger.info(f"Loading GLM results from {args.glm_dir}")
    try:
        glm_results = al.load_glm_results(args.glm_dir)
        logger.info(f"Loaded GLM results with {len(glm_results)} entries")
    except Exception as e:
        logger.warning(f"Could not load GLM results: {e}")
        glm_results = None

    # Load BPTF results
    logger.info(f"Loading BPTF results from {args.bptf_dir}")
    try:
        bptf_results = al.load_bptf_results(args.bptf_dir)
        n_motifs = bptf_results.get("n_components", "unknown")
        logger.info(f"Loaded BPTF results with {n_motifs} motifs")
    except Exception as e:
        logger.warning(f"Could not load BPTF results: {e}")
        bptf_results = None

    # Load projected adata if project-dir provided
    adata = None
    if args.project_dir:
        import os

        import anndata

        adata_path = os.path.join(args.project_dir, "projected_adata.h5ad")
        if os.path.exists(adata_path):
            logger.info(f"Loading projected adata from {adata_path}")
            adata = anndata.read_h5ad(adata_path)
            logger.info(f"Loaded adata with {adata.n_obs} cells, {adata.n_vars} genes")
        else:
            logger.warning(f"projected_adata.h5ad not found in {args.project_dir}")

    # Generate plots
    plots_generated = []

    # Helper function to write MultiQC metadata YAML
    def write_mqc_yaml(plot_path, section_name, description, plot_type="image"):
        """Write companion MultiQC YAML file for a plot"""
        yaml_path = plot_path.with_suffix(".yaml")
        mqc_data = {
            "id": plot_path.stem,
            "section_name": section_name,
            "description": description,
            "plot_type": plot_type,
        }
        with open(yaml_path, "w") as f:
            yaml.dump(mqc_data, f, default_flow_style=False)

    # Volcano plots
    if "volcano" in plot_types and glm_results is not None:
        logger.info("Generating volcano plots...")
        import numpy as np

        for motif_key, df in glm_results.items():
            if isinstance(df, pd.DataFrame) and COLUMN_NAME_LOGFC in df.columns:
                # Compute -log10(qval) for y-axis
                df = df.copy()
                df[COLUMN_NAME_NEG_LOG10_Q] = -np.log10(
                    df[COLUMN_NAME_QVAL].clip(1e-300)
                )

                fig, ax = plt.subplots(figsize=(10, 10))
                al.volcano_plot(
                    df,
                    x_col=COLUMN_NAME_LOGFC,
                    y_col=COLUMN_NAME_NEG_LOG10_Q,
                    label_col=COLUMN_NAME_GENE,
                    fdr=args.alpha,
                    ax=ax,
                )
                ax.set_title(f"Volcano Plot - {motif_key}")
                outpath = output_dir / f"volcano_{motif_key}_mqc.{args.format}"
                fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
                plt.close(fig)
                write_mqc_yaml(
                    outpath,
                    f"Volcano Plot: {motif_key}",
                    f"Differential expression volcano plot for {motif_key}",
                )
                plots_generated.append(str(outpath))
        logger.info(
            f"Generated {len([p for p in plots_generated if 'volcano' in p])} volcano plots"
        )

    # Forest plots
    if "forest" in plot_types and glm_results is not None:
        logger.info("Generating forest plots...")
        for motif_key, df in glm_results.items():
            if isinstance(df, pd.DataFrame) and COLUMN_NAME_LOGFC in df.columns:
                fig, ax = plt.subplots(figsize=(8, 10))
                al.forest_plot(df, ax=ax, n_top=args.top_n)
                ax.set_title(f"Forest Plot - {motif_key}")
                outpath = output_dir / f"forest_{motif_key}_mqc.{args.format}"
                fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
                plt.close(fig)
                write_mqc_yaml(
                    outpath,
                    f"Forest Plot: {motif_key}",
                    f"Top differentially expressed genes forest plot for {motif_key}",
                )
                plots_generated.append(str(outpath))
        logger.info(
            f"Generated {len([p for p in plots_generated if 'forest' in p])} forest plots"
        )

    # Heatmap of motif activities
    if "heatmap" in plot_types and bptf_results is not None:
        logger.info("Generating motif heatmap...")
        patch_loadings = bptf_results["patch_loadings"]
        fig = al.plot_motif_activities(patch_loadings)
        if fig is not None:
            outpath = output_dir / f"motif_heatmap_mqc.{args.format}"
            fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
            write_mqc_yaml(
                outpath,
                "Motif Activity Heatmap",
                "Heatmap showing motif activities across spatial patches",
            )
            plots_generated.append(str(outpath))
            logger.info("Generated motif heatmap")

    # Motif summary plots
    if "motif_summary" in plot_types and bptf_results is not None:
        logger.info("Generating motif summary plots...")
        patch_loadings = bptf_results["patch_loadings"]
        lri_factors = bptf_results["lri_factors"]

        fig = al.plot_factor_distributions(patch_loadings, lri_factors)
        if fig is not None:
            outpath = output_dir / f"motif_distributions_mqc.{args.format}"
            fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
            write_mqc_yaml(
                outpath,
                "Motif Distributions",
                "Distribution of motif loadings across patches and LRI factors",
            )
            plots_generated.append(str(outpath))

        fig = al.plot_bptf_diagnostics(patch_loadings, lri_factors)
        if fig is not None:
            outpath = output_dir / f"bptf_diagnostics_mqc.{args.format}"
            fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
            write_mqc_yaml(
                outpath,
                "BPTF Diagnostics",
                "Diagnostic plots for Bayesian Poisson Tensor Factorization",
            )
            plots_generated.append(str(outpath))

        logger.info("Generated motif summary plots")

    # Spatial plots
    if "spatial" in plot_types and adata is not None and bptf_results is not None:
        logger.info("Generating spatial plots...")
        import numpy as np

        # Check if spatial coordinates exist
        if "spatial" not in adata.obsm:
            logger.warning(
                "No spatial coordinates found in adata.obsm['spatial'], skipping spatial plots"
            )
        else:
            coords = adata.obsm["spatial"][:, :2]
            n_motifs = bptf_results["patch_loadings"].shape[1]

            # Get cell type colors if available
            cell_type_col = "cell_type"
            has_celltypes = cell_type_col in adata.obs.columns

            # Plot 1: Cells colored by cell type (reference)
            if has_celltypes:
                logger.info("  Generating cell type spatial plot...")
                fig, ax = plt.subplots(figsize=(10, 10))
                cell_types = adata.obs[cell_type_col].astype("category")
                n_types = len(cell_types.cat.categories)

                # Use a qualitative colormap
                cmap = plt.cm.get_cmap("tab20", n_types)
                colors = [cmap(i) for i in range(n_types)]
                color_map = dict(zip(cell_types.cat.categories, colors))
                cell_colors = [color_map[ct] for ct in cell_types]

                scatter = ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=cell_colors,
                    s=1,
                    alpha=0.8,
                    rasterized=True,
                )
                ax.set_aspect("equal")
                ax.set_xlabel("X")
                ax.set_ylabel("Y")
                ax.set_title("Cells by Cell Type")

                # Add legend (top N cell types by count)
                ct_counts = cell_types.value_counts()
                top_cts = ct_counts.head(15).index
                handles = [
                    plt.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=color_map[ct],
                        markersize=8,
                        label=ct,
                    )
                    for ct in top_cts
                ]
                ax.legend(
                    handles=handles,
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                    fontsize=8,
                    title="Cell Type",
                )

                outpath = output_dir / f"spatial_celltypes_mqc.{args.format}"
                fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
                plt.close(fig)
                write_mqc_yaml(
                    outpath,
                    "Spatial Cell Types",
                    "Spatial distribution of cells colored by cell type",
                )
                plots_generated.append(str(outpath))

            # Plot 2: Cells colored by motif loading (one per motif)
            # Check if cell loadings are in adata.obsm
            loading_key = None
            for key in ["cell_motif_loadings", "X_motif", "motif_loadings"]:
                if key in adata.obsm:
                    loading_key = key
                    break

            if loading_key is None:
                # Try to load from parquet file
                import os

                loadings_path = os.path.join(
                    args.project_dir, "cell_motif_loadings.parquet"
                )
                if os.path.exists(loadings_path):
                    logger.info(f"  Loading cell loadings from {loadings_path}")
                    cell_loadings = pd.read_parquet(loadings_path)
                    # Assume columns are motif_0, motif_1, etc.
                    loading_cols = [
                        c for c in cell_loadings.columns if c.startswith("motif_")
                    ]
                    if loading_cols:
                        cell_loadings_arr = cell_loadings[sorted(loading_cols)].values
                    else:
                        cell_loadings_arr = None
                else:
                    cell_loadings_arr = None
            else:
                cell_loadings_arr = adata.obsm[loading_key]

            if cell_loadings_arr is not None:
                n_motifs = cell_loadings_arr.shape[1]
                logger.info(f"  Generating {n_motifs} motif loading spatial plots...")

                for k in range(n_motifs):
                    loadings = cell_loadings_arr[:, k]

                    fig, ax = plt.subplots(figsize=(10, 10))
                    scatter = ax.scatter(
                        coords[:, 0],
                        coords[:, 1],
                        c=loadings,
                        cmap="viridis",
                        s=1,
                        alpha=0.8,
                        rasterized=True,
                    )
                    ax.set_aspect("equal")
                    ax.set_xlabel("X")
                    ax.set_ylabel("Y")
                    ax.set_title(f"Motif {k} Loading")
                    plt.colorbar(scatter, ax=ax, label="Loading")

                    outpath = (
                        output_dir / f"spatial_motif_{k}_loading_mqc.{args.format}"
                    )
                    fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
                    plt.close(fig)
                    write_mqc_yaml(
                        outpath,
                        f"Spatial Motif {k} Loading",
                        f"Spatial distribution of cells colored by Motif {k} loading",
                    )
                    plots_generated.append(str(outpath))

                logger.info(f"  Generated {n_motifs} motif loading plots")
            else:
                logger.warning(
                    "  No cell motif loadings found, skipping motif loading plots"
                )

            logger.info("Generated spatial plots")

    # Report results
    logger.info(f"Generated {len(plots_generated)} plots")
    logger.info(f"Results saved to: {args.output_dir}")

    if plots_generated:
        logger.info("Generated files:")
        for p in plots_generated:
            logger.info(f"  - {p}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
