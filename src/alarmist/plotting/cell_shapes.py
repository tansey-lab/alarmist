"""
Cell boundary polygons from Xenium Ranger output, for drawing cells as their
segmented shapes instead of scatter points.

Requires the optional ``xenium`` extra (``pip install 'alarmist[xenium]'``),
which provides ``spatialdata-io`` and ``shapely``.
"""

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection

logger = logging.getLogger(__name__)

_XENIUM_EXTRA_HINT = "install the xenium extra: pip install 'alarmist[xenium]'"


def load_xenium_cell_shapes(xenium_dir: str):
    """
    Read cell boundary polygons from a Xenium Ranger output directory.

    Only the cell boundaries are parsed; images, transcripts, labels and the
    cell table are skipped to keep the read fast.

    Parameters
    ----------
    xenium_dir : str
        Xenium Ranger output directory (the one containing ``experiment.xenium``).

    Returns
    -------
    geopandas.GeoDataFrame
        One row per cell, indexed by Xenium ``cell_id``. Geometry is in the
        element's intrinsic coordinates (µm), the same frame as the
        ``x_centroid``/``y_centroid`` that usually populate ``obsm['spatial']``.
    """
    try:
        import spatialdata_io
    except ImportError as e:
        raise ImportError(
            f"Reading Xenium cell shapes requires spatialdata-io; {_XENIUM_EXTRA_HINT}"
        ) from e

    logger.info(f"Reading Xenium cell boundaries from {xenium_dir}")
    sdata = spatialdata_io.xenium(
        xenium_dir,
        cells_boundaries=True,
        nucleus_boundaries=False,
        cells_as_circles=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
        aligned_images=False,
        cells_table=False,
    )
    if "cell_boundaries" not in sdata.shapes:
        raise ValueError(
            f"No 'cell_boundaries' shapes found in {xenium_dir}; "
            f"available shapes: {list(sdata.shapes.keys())}"
        )
    shapes = sdata.shapes["cell_boundaries"]
    shapes.index = shapes.index.astype(str)
    logger.info(f"Loaded {len(shapes):,} cell boundary polygons")
    return shapes


def align_cell_shapes(
    shapes, cell_ids: Sequence[str] | pd.Index
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reorder cell polygons to match a list of cell IDs.

    Parameters
    ----------
    shapes : geopandas.GeoDataFrame
        Polygons indexed by cell ID (e.g. from :func:`load_xenium_cell_shapes`).
    cell_ids : sequence of str
        Target order, typically ``adata.obs_names``.

    Returns
    -------
    geoms : np.ndarray of object
        Geometry per cell ID, ``None`` where the ID has no polygon.
    matched : np.ndarray of bool
        Whether each cell ID was found in ``shapes``.
    """
    ids = pd.Index(np.asarray(cell_ids).astype(str))
    if not shapes.index.is_unique:
        raise ValueError("Cell shape index contains duplicate cell IDs")
    geoms = shapes.geometry.reindex(ids).to_numpy(dtype=object)
    matched = ids.isin(shapes.index)
    geoms[~matched] = None
    return geoms, np.asarray(matched)


def check_shape_alignment(
    geoms: np.ndarray, coords: np.ndarray, matched: np.ndarray
) -> float:
    """
    Median distance between polygon centroids and the cells' recorded coordinates.

    A value far above a cell radius means the polygons and ``obsm['spatial']``
    are in different frames (e.g. pixels vs µm) or the IDs are mis-paired.

    Returns ``nan`` when no cell is matched.
    """
    import shapely

    if not matched.any():
        return float("nan")
    centroids = shapely.get_coordinates(shapely.centroid(geoms[matched]))
    dist = np.linalg.norm(centroids - np.asarray(coords)[matched, :2], axis=1)
    return float(np.median(dist))


def _polygon_vertices(geoms: np.ndarray) -> tuple[list[np.ndarray], np.ndarray]:
    """Exterior-ring vertex arrays for every polygon part, plus each part's owning geom index."""
    import shapely

    parts, owner = shapely.get_parts(geoms, return_index=True)
    rings = shapely.get_exterior_ring(parts)
    xy, ring_idx = shapely.get_coordinates(rings, return_index=True)
    splits = np.flatnonzero(np.diff(ring_idx)) + 1
    return np.split(xy, splits), owner


def draw_cell_shapes(
    ax: Axes,
    geoms: np.ndarray,
    *,
    color=None,
    facecolors: Sequence | None = None,
    values: np.ndarray | None = None,
    cmap=None,
    norm=None,
    alpha: float | None = None,
    rasterized: bool = True,
) -> PolyCollection:
    """
    Draw cells as filled polygons in a single PolyCollection.

    Give exactly one of ``color`` (a single colour for every cell),
    ``facecolors`` (one colour per geom) or ``values`` (one scalar per geom,
    mapped through ``cmap``/``norm``). ``None`` geoms are skipped. MultiPolygons contribute one patch
    per part, all sharing the cell's colour.

    Returns the collection, which can be passed to ``plt.colorbar`` when
    ``values`` was used.
    """
    if sum(x is not None for x in (color, facecolors, values)) != 1:
        raise ValueError("Pass exactly one of color, facecolors or values")
    geoms = np.asarray(geoms, dtype=object)
    present = np.array([g is not None and not g.is_empty for g in geoms], dtype=bool)
    idx = np.flatnonzero(present)

    if idx.size:
        verts, owner = _polygon_vertices(geoms[idx])
        owner = idx[owner]
    else:
        verts, owner = [], np.array([], dtype=int)

    coll = PolyCollection(
        verts,
        edgecolors="none",
        linewidths=0,
        alpha=alpha,
        rasterized=rasterized,
    )
    if values is not None:
        coll.set_array(np.asarray(values)[owner])
        coll.set_cmap(cmap)
        if norm is not None:
            coll.set_norm(norm)
        else:
            vals = np.asarray(values, dtype=float)[present]
            if vals.size:
                coll.set_clim(np.nanmin(vals), np.nanmax(vals))
    elif facecolors is not None:
        if len(facecolors) != len(geoms):
            raise ValueError(
                f"facecolors has {len(facecolors)} entries for {len(geoms)} geoms"
            )
        coll.set_facecolor([facecolors[i] for i in owner])
    elif color is not None:
        coll.set_facecolor(color)

    ax.add_collection(coll)
    ax.autoscale_view()
    return coll
