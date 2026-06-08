"""Spatial contour polygons for motif loadings, exported as vector geometry.

For each motif its per-cell loading is thresholded at a percentile (default the
95th, computed globally across all cells) and the iso-contour of the smoothed,
interpolated loading field is traced as polygons. Interpolation is performed
*per spatial group* (e.g. one TMA core) so the field is never interpolated
across the empty gaps between disjoint tissue regions. The result is a
GeoDataFrame that can be written to GeoJSON (or any other OGR vector format).

Coordinates are passed through unchanged (typically spatial microns), so the
output geometry is planar; set ``crs`` if the coordinates carry one.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _require_geo():
    """Import the optional geo stack, raising a friendly error if missing."""
    try:
        import contourpy  # noqa: F401
        import geopandas  # noqa: F401
        import shapely  # noqa: F401
    except ImportError as e:  # pragma: no cover - exercised only without deps
        raise ImportError(
            "motif_loading_contours requires the optional geo dependencies "
            "'geopandas', 'shapely' and 'contourpy'. Install them with "
            "`pip install geopandas shapely contourpy`."
        ) from e


def _rings_to_polygons(lines, Polygon):
    """Turn closed contour polylines into shapely polygons.

    A contour that sits directly inside another is treated as a hole (one level
    of nesting, covering the common 'donut' hotspot). Deeper nesting is rare for
    upper-percentile loading fields and is not resolved.
    """
    polys = []
    for ln in lines:
        if len(ln) >= 4:
            p = Polygon(ln)
            if not p.is_valid:
                p = p.buffer(0)
            if (not p.is_empty) and p.geom_type == "Polygon":
                polys.append(p)
    if not polys:
        return []
    order = sorted(range(len(polys)), key=lambda i: polys[i].area, reverse=True)
    used: set[int] = set()
    out = []
    for i in order:
        if i in used:
            continue
        ext = list(polys[i].exterior.coords)
        holes = []
        for j in order:
            if j == i or j in used:
                continue
            if polys[i].contains(polys[j].representative_point()):
                holes.append(list(polys[j].exterior.coords))
                used.add(j)
        out.append(Polygon(ext, holes))
    return out


def motif_loading_contours(
    coords,
    loadings,
    group_labels=None,
    *,
    percentile: float = 95.0,
    grid_res: float = 20.0,
    smooth_sigma: float = 1.5,
    min_area: float = 2000.0,
    pad: float = 60.0,
    min_cells: int = 50,
    motif_names: dict | None = None,
    crs=None,
):
    """Trace per-motif percentile contours of the spatial loading field.

    Parameters
    ----------
    coords : array-like, shape (n_cells, 2)
        Spatial coordinates (x, y) per cell.
    loadings : array-like, shape (n_cells, n_motifs)
        Per-cell motif loading matrix.
    group_labels : array-like, shape (n_cells,), optional
        Spatial group per cell (e.g. TMA core id). Interpolation/contouring is
        done independently within each group so the field never bridges the gaps
        between disjoint tissue regions. If ``None`` all cells form one group
        (the field is then interpolated over the whole convex hull).
    percentile : float, default 95.0
        Threshold percentile, computed globally per motif across all cells.
    grid_res : float, default 20.0
        Interpolation grid spacing, in coordinate units.
    smooth_sigma : float, default 1.5
        Gaussian smoothing of the gridded field, in grid cells.
    min_area : float, default 2000.0
        Drop contour polygons smaller than this (coordinate-unit^2).
    pad : float, default 60.0
        Padding added around each group's bounding box so edge contours close.
    min_cells : int, default 50
        Skip a group with fewer cells than this.
    motif_names : dict, optional
        Mapping motif index -> display name, added as a ``motif_name`` field.
    crs : optional
        CRS for the output GeoDataFrame (default ``None`` -> planar).

    Returns
    -------
    geopandas.GeoDataFrame
        One polygon feature per (motif, group, contour) with attributes
        ``motif``, [``motif_name``], ``group``, ``threshold``, ``pctile`` and
        ``area_um2``.
    """
    _require_geo()
    import geopandas as gpd
    from contourpy import LineType, contour_generator
    from scipy.interpolate import LinearNDInterpolator
    from scipy.ndimage import gaussian_filter
    from scipy.spatial import Delaunay
    from shapely.geometry import Polygon

    coords = np.asarray(coords, dtype=float)
    loadings = np.asarray(loadings, dtype=float)
    if loadings.ndim != 2:
        raise ValueError("loadings must be 2D (n_cells, n_motifs)")
    n_cells, n_motifs = loadings.shape
    if coords.shape[0] != n_cells:
        raise ValueError("coords and loadings must have the same number of rows")

    if group_labels is None:
        group_labels = np.zeros(n_cells, dtype=object)
    else:
        group_labels = np.asarray(group_labels)

    thresholds = {
        k: float(np.nanpercentile(loadings[:, k], percentile))
        for k in range(n_motifs)
    }

    def _is_missing(g):
        return g is None or (isinstance(g, float) and np.isnan(g))

    groups = [g for g in pd.unique(group_labels) if not _is_missing(g)]
    has_names = motif_names is not None
    records = []

    for grp in groups:
        mask = group_labels == grp
        if int(mask.sum()) < min_cells:
            continue
        pts = coords[mask]
        try:
            tri = Delaunay(pts)
        except Exception as e:  # collinear / too few unique points
            logger.warning("contours: group %r could not triangulate (%s)", grp, e)
            continue

        x0, y0 = pts[:, 0].min() - pad, pts[:, 1].min() - pad
        x1, y1 = pts[:, 0].max() + pad, pts[:, 1].max() + pad
        nx = max(8, int(np.ceil((x1 - x0) / grid_res)))
        ny = max(8, int(np.ceil((y1 - y0) / grid_res)))
        GX, GY = np.meshgrid(np.linspace(x0, x1, nx), np.linspace(y0, y1, ny))

        for k in range(n_motifs):
            thr = thresholds[k]
            Z = LinearNDInterpolator(tri, loadings[mask, k])(GX, GY)
            valid = ~np.isnan(Z)
            if int(valid.sum()) < 10:
                continue
            # NaN-aware gaussian smoothing (normalized convolution).
            zf = gaussian_filter(np.where(valid, Z, 0.0), smooth_sigma, mode="nearest")
            wf = gaussian_filter(valid.astype(float), smooth_sigma, mode="nearest")
            with np.errstate(invalid="ignore", divide="ignore"):
                Zs = np.where(wf > 1e-6, zf / wf, np.nan)
            Zs[~valid] = np.nan
            finite = Zs[np.isfinite(Zs)]
            if finite.size == 0 or float(finite.max()) < thr:
                continue
            # Fill outside-tissue with a sub-threshold value so the iso-contour
            # closes along the tissue boundary instead of running off-grid.
            Zc = np.where(np.isnan(Zs), float(finite.min()) - 1.0, Zs)

            cg = contour_generator(GX, GY, Zc, line_type=LineType.Separate)
            for g in _rings_to_polygons(cg.lines(thr), Polygon):
                if g.is_empty or g.area < min_area:
                    continue
                rec = {
                    "motif": int(k),
                    "group": str(grp),
                    "threshold": float(thr),
                    "pctile": float(percentile),
                    "area_um2": float(g.area),
                    "geometry": g,
                }
                if has_names:
                    rec["motif_name"] = str(motif_names.get(k, f"Motif {k}"))[:80]
                records.append(rec)

    cols = ["motif"]
    if has_names:
        cols.append("motif_name")
    cols += ["group", "threshold", "pctile", "area_um2", "geometry"]

    if not records:
        logger.warning("motif_loading_contours: no contour polygons produced")
        return gpd.GeoDataFrame(columns=cols, geometry="geometry", crs=crs)
    return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)[cols]


def motif_loading_contours_from_adata(
    adata,
    *,
    spatial_key: str = "spatial",
    loading_key: str | None = None,
    group_column: str | None = None,
    **kwargs,
):
    """Convenience wrapper: build motif loading contours from an AnnData.

    Pulls coordinates from ``adata.obsm[spatial_key]`` and the loading matrix
    from ``adata.obsm[loading_key]`` (auto-detected among ``X_motif``,
    ``cell_motif_loadings``, ``motif_loadings`` when ``loading_key`` is None).
    ``group_column`` names an ``adata.obs`` column to group by (e.g. ``tma_id``).
    Remaining keyword arguments are forwarded to :func:`motif_loading_contours`.
    """
    if spatial_key not in adata.obsm:
        raise KeyError(f"adata.obsm[{spatial_key!r}] not found")
    coords = np.asarray(adata.obsm[spatial_key])[:, :2]

    if loading_key is None:
        for key in ("X_motif", "cell_motif_loadings", "motif_loadings"):
            if key in adata.obsm:
                loading_key = key
                break
        if loading_key is None:
            raise KeyError(
                "no motif loading matrix found in adata.obsm "
                "(looked for X_motif, cell_motif_loadings, motif_loadings)"
            )
    loadings = np.asarray(adata.obsm[loading_key])

    group_labels = None
    if group_column is not None:
        if group_column not in adata.obs.columns:
            raise KeyError(f"group column {group_column!r} not in adata.obs")
        group_labels = adata.obs[group_column].to_numpy()

    return motif_loading_contours(coords, loadings, group_labels, **kwargs)
