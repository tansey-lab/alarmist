"""
Tests for drawing cells as Xenium boundary polygons (alarmist-visualize
--xenium-ranger-dir) without needing a real Xenium Ranger output.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pytest

gpd = pytest.importorskip("geopandas")
shapely = pytest.importorskip("shapely")

from matplotlib.collections import PolyCollection  # noqa: E402

from alarmist.cli import visualize  # noqa: E402
from alarmist.plotting import cell_shapes  # noqa: E402
from alarmist.plotting.cell_shapes import (  # noqa: E402
    align_cell_shapes,
    check_shape_alignment,
    draw_cell_shapes,
)


def _square(cx, cy, half=2.0):
    return shapely.box(cx - half, cy - half, cx + half, cy + half)


def _shapes_for(adata, shift=0.0):
    """One square per cell centred on obsm['spatial'], indexed like Xenium cell_id."""
    coords = adata.obsm["spatial"]
    geoms = [_square(x + shift, y + shift) for x, y in coords]
    return gpd.GeoDataFrame(
        {"geometry": geoms}, index=adata.obs_names.astype(str).rename("cell_id")
    )


def test_align_cell_shapes_reorders_and_flags_missing():
    shapes = gpd.GeoDataFrame(
        {"geometry": [_square(0, 0), _square(10, 10)]}, index=["a", "b"]
    )
    geoms, matched = align_cell_shapes(shapes, ["b", "x", "a"])
    assert matched.tolist() == [True, False, True]
    assert geoms[1] is None
    assert geoms[0].equals(shapes.geometry["b"])
    assert geoms[2].equals(shapes.geometry["a"])


def test_check_shape_alignment(adata_single):
    coords = adata_single.obsm["spatial"]
    matched = np.ones(adata_single.n_obs, dtype=bool)
    centred, _ = align_cell_shapes(_shapes_for(adata_single), adata_single.obs_names)
    shifted, _ = align_cell_shapes(
        _shapes_for(adata_single, shift=30.0), adata_single.obs_names
    )
    assert check_shape_alignment(centred, coords, matched) == pytest.approx(0, abs=1e-9)
    assert check_shape_alignment(shifted, coords, matched) == pytest.approx(
        30 * np.sqrt(2)
    )


def test_draw_cell_shapes_values_multipolygon_and_missing():
    multi = shapely.MultiPolygon([_square(0, 0, 1), _square(5, 5, 1)])
    geoms = np.array([_square(10, 10), None, multi], dtype=object)
    fig, ax = plt.subplots()
    coll = draw_cell_shapes(ax, geoms, values=np.array([1.0, 99.0, 3.0]))
    # 1 polygon + 2 MultiPolygon parts; the None cell is skipped
    assert isinstance(coll, PolyCollection)
    assert len(coll.get_paths()) == 3
    assert coll.get_array().tolist() == [1.0, 3.0, 3.0]
    # colour limits ignore the undrawn cell
    assert coll.get_clim() == (1.0, 3.0)
    plt.close(fig)


def test_draw_cell_shapes_per_cell_facecolors():
    geoms = np.array([_square(0, 0), _square(5, 5)], dtype=object)
    fig, ax = plt.subplots()
    coll = draw_cell_shapes(ax, geoms, facecolors=["red", "blue"])
    fc = coll.get_facecolor()
    assert tuple(fc[0][:3]) == (1.0, 0.0, 0.0)
    assert tuple(fc[1][:3]) == (0.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        draw_cell_shapes(ax, geoms, color="red", values=np.ones(2))
    plt.close(fig)


def _run_visualize(monkeypatch, tmp_path, adata, shapes, extra_args=()):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    adata.write_h5ad(project_dir / "projected_adata.h5ad")
    out_dir = tmp_path / "plots"

    monkeypatch.setattr(cell_shapes, "load_xenium_cell_shapes", lambda _dir: shapes)
    n_motifs = adata.obsm["cell_motif_loadings"].shape[1]
    monkeypatch.setattr(
        "alarmist.load_bptf_results",
        lambda _dir: {
            "patch_loadings": np.zeros((4, n_motifs)),
            "n_components": n_motifs,
        },
    )
    drawn = []
    real_draw = cell_shapes.draw_cell_shapes
    monkeypatch.setattr(
        cell_shapes,
        "draw_cell_shapes",
        lambda *a, **kw: drawn.append(real_draw(*a, **kw)) or drawn[-1],
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "alarmist-visualize",
            "--glm-dir",
            str(tmp_path / "no_glm"),
            "--bptf-dir",
            str(tmp_path / "no_bptf"),
            "--project-dir",
            str(project_dir),
            "--output-dir",
            str(out_dir),
            "--plot-types",
            "spatial",
            "motif_states",
            "--xenium-ranger-dir",
            str(tmp_path / "xenium"),
            *extra_args,
        ],
    )
    visualize.main()
    return out_dir, drawn


def test_visualize_draws_polygons(monkeypatch, tmp_path, adata_single):
    adata = adata_single
    adata.obsm["cell_motif_loadings"] = adata.obs[
        [f"motif_{k}_loading" for k in range(3)]
    ].to_numpy()
    # Xenium-style IDs held in an obs column; obs_names don't match
    adata.obs["xenium_id"] = [f"cell{i:05d}-1" for i in range(adata.n_obs)]
    shapes = _shapes_for(adata)
    shapes.index = adata.obs["xenium_id"].to_numpy()

    out_dir, drawn = _run_visualize(
        monkeypatch,
        tmp_path,
        adata,
        shapes,
        extra_args=["--xenium-cell-id-column", "xenium_id"],
    )

    assert (out_dir / "spatial_celltypes_S1_mqc.png").exists()
    for k in range(3):
        assert (out_dir / f"spatial_motif_{k}_loading_S1_mqc.png").exists()
        assert (out_dir / f"spatial_motif_{k}_state_S1_mqc.png").exists()
    # 1 cell-type map + 3 loading maps + 3 motifs × (OFF, ON)
    assert len(drawn) == 1 + 3 + 3 * 2
    assert all(isinstance(c, PolyCollection) for c in drawn)
    assert sum(len(c.get_paths()) for c in drawn[:1]) == adata.n_obs


def test_visualize_rejects_unmatched_ids(monkeypatch, tmp_path, adata_single):
    adata_single.obsm["cell_motif_loadings"] = np.ones((adata_single.n_obs, 3))
    shapes = _shapes_for(adata_single)
    shapes.index = [f"other{i}" for i in range(adata_single.n_obs)]
    with pytest.raises(ValueError, match="match a Xenium cell_id"):
        _run_visualize(monkeypatch, tmp_path, adata_single, shapes)


def test_visualize_rejects_multiple_samples(monkeypatch, tmp_path, adata_merged):
    adata_merged.obsm["cell_motif_loadings"] = np.ones((adata_merged.n_obs, 3))
    with pytest.raises(ValueError, match="single Xenium Ranger output"):
        _run_visualize(monkeypatch, tmp_path, adata_merged, _shapes_for(adata_merged))
