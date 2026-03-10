"""
CLI command: alarmist-visualize

Generate visualizations from ALARMIST pipeline results.
"""

import argparse
import sys

from alarmist.cli import common, log_config


def get_parser():
    """Create argument parser for visualize command"""
    parser = argparse.ArgumentParser(
        prog='alarmist-visualize',
        description='Generate visualizations from ALARMIST pipeline results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - generate all plots
  alarmist-visualize --glm-dir results/glm --bptf-dir results/bptf --output-dir results/plots

  # Generate only specific plot types
  alarmist-visualize --glm-dir results/glm --bptf-dir results/bptf --output-dir results/plots \\
      --plot-types volcano forest heatmap
        """
    )

    # Input arguments
    parser.add_argument(
        '--glm-dir', '-g',
        type=str,
        required=True,
        help='Directory containing GLM results (from alarmist-glm)'
    )
    parser.add_argument(
        '--bptf-dir', '-b',
        type=str,
        required=True,
        help='Directory containing BPTF results (from alarmist-bptf)'
    )
    parser.add_argument(
        '--patch-lri-dir',
        type=str,
        default=None,
        help='Directory containing patch-LRI results (optional)'
    )

    common.add_output_arguments(parser)

    # Plot options
    parser.add_argument(
        '--plot-types',
        type=str,
        nargs='+',
        default=['volcano', 'forest', 'heatmap', 'motif_summary'],
        choices=['volcano', 'forest', 'heatmap', 'motif_summary', 'spatial', 'all'],
        help='Types of plots to generate (default: volcano forest heatmap motif_summary)'
    )
    parser.add_argument(
        '--format',
        type=str,
        default='png',
        choices=['png', 'pdf', 'svg'],
        help='Output format for plots (default: png)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='DPI for raster outputs (default: 150)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.05,
        help='Significance threshold for highlighting (default: 0.05)'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=20,
        help='Number of top features to show in summary plots (default: 20)'
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
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import alarmist as al

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve 'all' plot type
    plot_types = args.plot_types
    if 'all' in plot_types:
        plot_types = ['volcano', 'forest', 'heatmap', 'motif_summary']

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
        n_motifs = bptf_results.get('n_components', 'unknown')
        logger.info(f"Loaded BPTF results with {n_motifs} motifs")
    except Exception as e:
        logger.warning(f"Could not load BPTF results: {e}")
        bptf_results = None

    # Generate plots
    plots_generated = []

    # Volcano plots
    if 'volcano' in plot_types and glm_results is not None:
        logger.info("Generating volcano plots...")
        try:
            for motif_key, df in glm_results.items():
                if isinstance(df, pd.DataFrame) and 'log2fc' in df.columns.str.lower():
                    fig = al.volcano_plot(
                        df,
                        alpha=args.alpha,
                        title=f'Volcano Plot - {motif_key}'
                    )
                    if fig is not None:
                        outpath = output_dir / f'volcano_{motif_key}.{args.format}'
                        fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
                        plt.close(fig)
                        plots_generated.append(str(outpath))
            logger.info(f"Generated {len([p for p in plots_generated if 'volcano' in p])} volcano plots")
        except Exception as e:
            logger.warning(f"Could not generate volcano plots: {e}")

    # Forest plots
    if 'forest' in plot_types and glm_results is not None:
        logger.info("Generating forest plots...")
        try:
            for motif_key, df in glm_results.items():
                if isinstance(df, pd.DataFrame) and 'coef' in df.columns.str.lower():
                    fig = al.forest_plot(
                        df,
                        top_n=args.top_n,
                        title=f'Forest Plot - {motif_key}'
                    )
                    if fig is not None:
                        outpath = output_dir / f'forest_{motif_key}.{args.format}'
                        fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
                        plt.close(fig)
                        plots_generated.append(str(outpath))
            logger.info(f"Generated {len([p for p in plots_generated if 'forest' in p])} forest plots")
        except Exception as e:
            logger.warning(f"Could not generate forest plots: {e}")

    # Heatmap of motif activities
    if 'heatmap' in plot_types and bptf_results is not None:
        logger.info("Generating motif heatmap...")
        try:
            fig = al.plot_motif_activities(bptf_results)
            if fig is not None:
                outpath = output_dir / f'motif_heatmap.{args.format}'
                fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
                plt.close(fig)
                plots_generated.append(str(outpath))
                logger.info("Generated motif heatmap")
        except Exception as e:
            logger.warning(f"Could not generate heatmap: {e}")

    # Motif summary plots
    if 'motif_summary' in plot_types and bptf_results is not None:
        logger.info("Generating motif summary plots...")
        try:
            fig = al.plot_factor_distributions(bptf_results)
            if fig is not None:
                outpath = output_dir / f'motif_distributions.{args.format}'
                fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
                plt.close(fig)
                plots_generated.append(str(outpath))

            fig = al.plot_bptf_diagnostics(bptf_results)
            if fig is not None:
                outpath = output_dir / f'bptf_diagnostics.{args.format}'
                fig.savefig(outpath, dpi=args.dpi, bbox_inches='tight')
                plt.close(fig)
                plots_generated.append(str(outpath))

            logger.info("Generated motif summary plots")
        except Exception as e:
            logger.warning(f"Could not generate motif summary plots: {e}")

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
