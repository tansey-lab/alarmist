"""
Cell type color management for consistent plotting across alarmist.

Usage:
    import alarmist as al

    # Set colors from list (auto-generate using palette)
    al.set_celltype_colors(['Tumor', 'T cell', 'Macrophage'], palette='tab20')

    # Set colors from adata
    al.set_celltype_colors(adata, column='cell_type', palette='Set2')

    # Set custom colors
    al.set_celltype_colors({'Tumor': '#E41A1C', 'T cell': '#377EB8'})

    # Get current colors
    colors = al.get_celltype_colors()

    # Clear colors
    al.clear_celltype_colors()
"""

from typing import Any, Union

import anndata
import matplotlib.pyplot as plt

# Module-level registry for cell type colors
_CELLTYPE_COLORS: dict[str, Any] = {}


def set_celltype_colors(
    source: Union[list[str], dict[str, Any], "anndata.AnnData"],
    column: str | None = None,
    palette: str = "tab20",
) -> dict[str, Any]:
    """
    Set global cell type colors for consistent plotting.

    Parameters
    ----------
    source : list, dict, or AnnData
        - list of str: Cell type names, colors auto-generated from palette
        - dict: Mapping of cell type -> color (hex or named color)
        - AnnData: Extract unique values from specified column
    column : str, optional
        Column name in adata.obs to extract cell types from.
        Required if source is AnnData.
    palette : str, default "tab20"
        Matplotlib colormap name for auto-generating colors.
        Common options: 'tab20', 'tab10', 'Set1', 'Set2', 'Set3', 'Paired'

    Returns
    -------
    dict
        The color mapping that was set

    Examples
    --------
    >>> al.set_celltype_colors(['Tumor', 'T cell', 'Macrophage'])
    >>> al.set_celltype_colors(['Tumor', 'T cell'], palette='Set1')
    >>> al.set_celltype_colors({'Tumor': 'red', 'T cell': 'blue'})
    >>> al.set_celltype_colors(adata, column='cell_type')
    """
    global _CELLTYPE_COLORS

    if isinstance(source, dict):
        # Direct color mapping
        _CELLTYPE_COLORS = source.copy()
        print(f"Set {len(_CELLTYPE_COLORS)} cell type colors (custom)")
    elif isinstance(source, list):
        # Generate colors from palette
        _CELLTYPE_COLORS = _generate_colors(source, palette)
        print(f"Set {len(_CELLTYPE_COLORS)} cell type colors using palette '{palette}'")
    else:
        # Assume AnnData
        if column is None:
            raise ValueError("column parameter required when source is AnnData")
        celltypes = sorted(source.obs[column].unique().tolist())
        _CELLTYPE_COLORS = _generate_colors(celltypes, palette)
        print(
            f"Set {len(_CELLTYPE_COLORS)} cell type colors from adata.obs['{column}'] using palette '{palette}'"
        )
    return _CELLTYPE_COLORS.copy()


def get_celltype_colors() -> dict[str, Any]:
    """
    Get current global cell type color mapping.

    Returns
    -------
    dict
        Mapping of cell type -> color. Empty dict if not set.
    """
    return _CELLTYPE_COLORS.copy()


def clear_celltype_colors() -> None:
    """Clear the global cell type color registry."""
    global _CELLTYPE_COLORS
    _CELLTYPE_COLORS = {}
    print("Cleared cell type colors")


def _generate_colors(celltypes: list[str], palette: str = "tab20") -> dict[str, tuple]:
    """Generate color mapping from a matplotlib palette."""
    n = len(celltypes)
    cmap = plt.get_cmap(palette, n)
    return {ct: cmap(i) for i, ct in enumerate(celltypes)}


def _get_colors_for_plotting(
    ct_colors: dict[str, Any] | None = None,
    df_celltypes: list[str] | None = None,
    palette: str = "tab20",
) -> dict[str, Any]:
    """
    Internal function to get colors for plotting functions.

    Priority:
    1. ct_colors parameter (if provided)
    2. Global registry (if set)
    3. Auto-generate from df_celltypes (if provided)

    Parameters
    ----------
    ct_colors : dict, optional
        Explicit color mapping passed to plotting function
    df_celltypes : list, optional
        Cell types extracted from DataFrame (fallback)
    palette : str, default "tab20"
        Palette to use if auto-generating

    Returns
    -------
    dict
        Color mapping to use
    """
    # Priority 1: explicit parameter
    if ct_colors is not None:
        return ct_colors

    # Priority 2: global registry
    if _CELLTYPE_COLORS:
        return _CELLTYPE_COLORS.copy()

    # Priority 3: auto-generate from DataFrame
    if df_celltypes is not None:
        return _generate_colors(sorted(set(df_celltypes)), palette)

    # No colors available
    return {}
