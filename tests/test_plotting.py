"""
Tests for plotting functions moved out of the analysis notebooks, and for the
optional kwargs folded into existing plotting functions.
"""

import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from alarmist.core.similarity import match_motifs, row_similarity
from alarmist.core.single_cell import compute_motif_celltype_enrichment
from alarmist.plotting import (
    gsea_barplot,
    plot_cells_per_patch,
    plot_lri_networks,
    plot_motif_celltype_enrichment,
    plot_motif_score_spatial,
    plot_motif_similarity,
    plot_motif_spatial,
    plot_single_motif_lri_lollipop,
    plot_spatial_categorical,
    plot_total_counts_per_cell,
    plot_two_motif_spatial,
)


def _lollipop_rows(ax):
    """(ligand, receptor) label pairs drawn on a lollipop axes, top to bottom."""
    texts = [t.get_text() for t in ax.texts]
    return list(zip(texts[0::2], texts[1::2], strict=True))


# --------------------------------------------------------------------------
# Lollipop: include / merge_celltypes / hide_celltypes / dedup / xlabel
# --------------------------------------------------------------------------


def test_lollipop_defaults_unchanged(lri_motifs_df):
    fig = plot_single_motif_lri_lollipop(lri_motifs_df, motif_idx=0, top_n=10)
    ax = fig.axes[0]
    assert ax.get_xlabel() == "factor"
    # Default dedup keeps the top ligand per (signaling, receptor, sender, receiver)
    dfp = lri_motifs_df[lri_motifs_df.motif_idx == 0].sort_values(
        "factor", ascending=False
    )
    expected = dfp.drop_duplicates(
        ["signaling_type", "receptor", "celltype1", "celltype2"]
    ).nlargest(10, "factor")
    assert _lollipop_rows(ax) == list(
        zip(expected.ligand, expected.receptor, strict=True)
    )


def test_lollipop_include_single_and_pair(lri_motifs_df):
    fig = plot_single_motif_lri_lollipop(
        lri_motifs_df, motif_idx=0, include="Tumor", top_n=200, dedup=False
    )
    n_involving = (
        (lri_motifs_df.motif_idx == 0)
        & ((lri_motifs_df.celltype1 == "Tumor") | (lri_motifs_df.celltype2 == "Tumor"))
    ).sum()
    assert len(_lollipop_rows(fig.axes[0])) == n_involving
    assert "involving Tumor" in fig.axes[0].get_title()

    fig = plot_single_motif_lri_lollipop(
        lri_motifs_df,
        motif_idx=0,
        include=["Tumor", "Macrophage"],
        top_n=200,
        dedup=False,
    )
    # Two directions (Tumor→Mac, Mac→Tumor) × 6 LR pairs; same-type pairs dropped
    assert len(_lollipop_rows(fig.axes[0])) == 12


def test_lollipop_include_conflicts_with_sender(lri_motifs_df):
    with pytest.raises(ValueError, match="include"):
        plot_single_motif_lri_lollipop(
            lri_motifs_df, motif_idx=0, include="Tumor", sender_type="Tumor"
        )


def test_lollipop_merge_and_hide(lri_motifs_df):
    merge = {"Macrophage": "Immune", "T cell": "Immune"}
    fig = plot_single_motif_lri_lollipop(
        lri_motifs_df,
        motif_idx=0,
        merge_celltypes=merge,
        hide_celltypes=["Fibroblast"],
        sender_type="Immune",
        receiver_type="Tumor",
        top_n=50,
        dedup=False,
        xlabel="Loading",
    )
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Loading"
    # Macrophage→Tumor and T cell→Tumor collapse onto one Immune→Tumor row per LR
    assert len(_lollipop_rows(ax)) == 6
    sub = lri_motifs_df[
        (lri_motifs_df.motif_idx == 0)
        & lri_motifs_df.celltype1.isin(["Macrophage", "T cell"])
        & (lri_motifs_df.celltype2 == "Tumor")
    ]
    top_sum = sub.groupby(["ligand", "receptor"]).factor.sum().max()
    endpoints = [
        c.get_offsets()[0][0] for c in ax.collections if c.get_offsets()[0][0] > 0
    ]
    assert max(endpoints) == pytest.approx(top_sum)
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Immune" in legend_labels and "Fibroblast" not in legend_labels


# --------------------------------------------------------------------------
# plot_lri_networks: return_celltypes
# --------------------------------------------------------------------------


@pytest.mark.skipif(shutil.which("dot") is None, reason="Graphviz not installed")
def test_lri_networks_return_celltypes(lri_motifs_df):
    fig = plot_lri_networks(lri_motifs_df, top_n=50)
    assert isinstance(fig, plt.Figure)
    fig, cts = plot_lri_networks(lri_motifs_df, top_n=50, return_celltypes=True)
    assert isinstance(fig, plt.Figure)
    assert set(cts) == {0, 1, 2}
    assert all(c <= set(lri_motifs_df.celltype1) for c in cts.values())


# --------------------------------------------------------------------------
# plot_motif_spatial: new kwargs
# --------------------------------------------------------------------------


def test_motif_spatial_ax_title_and_legend(adata_single):
    fig, ax = plt.subplots()
    out = plot_motif_spatial(
        adata_single,
        motif_idx=0,
        ax=ax,
        title="S1",
        show_celltype_legend=False,
        ct_colors={"Tumor": "red"},
        invert_y=True,
    )
    assert out is fig
    assert ax.get_title() == "S1"
    assert ax.yaxis_inverted()
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert len(labels) == 2 and labels[0].startswith("OFF")


def test_motif_spatial_ax_rejects_grid(adata_dict):
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="single-panel"):
        plot_motif_spatial(adata_dict, motif_idx=0, ax=ax, color_by_celltype=False)


def test_motif_spatial_save_path_and_per_panel(adata_merged, tmp_path):
    out = tmp_path / "m0.png"
    plot_motif_spatial(
        adata_merged,
        motif_idx=0,
        sample_column="sample_id",
        color_by_celltype=False,
        save_path=str(out),
    )
    assert out.exists()

    paths = plot_motif_spatial(
        adata_merged,
        motif_idx=1,
        sample_column="sample_id",
        color_by_celltype=False,
        per_panel_files=True,
        output_dir=str(tmp_path / "panels"),
        file_format="svg",
    )
    assert sorted(p.split("/")[-1] for p in paths) == [
        "motif_1_S1_spatial.svg",
        "motif_1_S2_spatial.svg",
    ]


# --------------------------------------------------------------------------
# plot_cells_per_patch: DataFrame input
# --------------------------------------------------------------------------


def test_cells_per_patch_dataframe():
    df = pd.DataFrame({"n_cells": [3, 5, 5, 8]})
    fig = plot_cells_per_patch(df, show=False)
    stats = fig.axes[0].texts[0].get_text()
    assert "Total patches: 4" in stats and "Total cells: 21" in stats


# --------------------------------------------------------------------------
# New spatial plots
# --------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", ["adata_single", "adata_dict", "adata_merged"])
def test_motif_score_spatial_inputs(fixture, request, tmp_path):
    adata = request.getfixturevalue(fixture)
    kwargs = {"sample_column": "sample_id"} if fixture == "adata_merged" else {}
    out = tmp_path / "score.png"
    fig = plot_motif_score_spatial(
        adata, motif_idx=0, log=True, save_path=str(out), **kwargs
    )
    assert isinstance(fig, plt.Figure) and out.exists()
    n_panels = 1 if fixture == "adata_single" else 2
    # Scatter panels plus one shared colourbar axes
    panels = [a for a in fig.axes if a.get_xlabel() == "X"]
    assert len(panels) == n_panels
    assert len(fig.axes) == n_panels + 1


def test_motif_score_spatial_missing_column(adata_single):
    with pytest.raises(ValueError, match="not found"):
        plot_motif_score_spatial(adata_single, motif_idx=0, value_key="nope")


def test_two_motif_spatial_counts(adata_single):
    fig = plot_two_motif_spatial(adata_single, motif_a=0, motif_b=1, cell_type="Tumor")
    obs = adata_single.obs
    tumor = obs.cell_type == "Tumor"
    a = tumor & (obs.motif_0_state == "positive")
    b = tumor & (obs.motif_1_state == "positive")
    title = fig.axes[0].get_title()
    assert f"A&B={int((a & b).sum()):,}" in title
    assert f"(n={int(tumor.sum()):,})" in title


def test_spatial_categorical_split_and_grid(adata_merged, tmp_path):
    fig = plot_spatial_categorical(
        adata_merged,
        split_by="sample_id",
        grid_size=50,
        categories=["Tumor", "Macrophage"],
        colors={"Tumor": "red", "Macrophage": "blue"},
        invert_y=True,
        save_path=str(tmp_path / "cat.png"),
    )
    panels = [a for a in fig.axes if a.get_title()]
    assert [a.get_title() for a in panels] == ["S1", "S2"]
    assert all(a.yaxis_inverted() for a in panels)
    # Same limits in every panel (shared extent) and grid lines drawn
    assert panels[0].get_xlim() == panels[1].get_xlim()
    assert len(panels[0].collections) == 3  # scatter + vlines + hlines
    labels = [t.get_text() for t in fig.legends[0].get_texts()]
    assert labels == ["Tumor", "Macrophage"]


def test_spatial_categorical_obs_coords(adata_single):
    adata_single.obs["x"] = adata_single.obsm["spatial"][:, 0]
    adata_single.obs["y"] = adata_single.obsm["spatial"][:, 1]
    fig = plot_spatial_categorical(adata_single, x_col="x", y_col="y")
    assert fig.axes[0].get_xlabel() == "x"


# --------------------------------------------------------------------------
# Enrichment
# --------------------------------------------------------------------------


def test_enrichment_math():
    import anndata

    obs = pd.DataFrame(
        {
            "cell_type": ["A", "A", "A", "B"],
            "motif_0_state": ["positive", "positive", "negative", "positive"],
        },
        index=[f"c{i}" for i in range(4)],
    )
    enr = compute_motif_celltype_enrichment(anndata.AnnData(obs=obs))
    e = enr.set_index("cell_type")["enrichment"]
    # p(A)=3/4, p_0(A)=2/3 -> 8/9 ; p(B)=1/4, p_0(B)=1/3 -> 4/3
    assert e["A"] == pytest.approx(8 / 9)
    assert e["B"] == pytest.approx(4 / 3)


def test_enrichment_dict_pools_samples(adata_dict, adata_merged):
    from_dict = compute_motif_celltype_enrichment(adata_dict)
    from_merged = compute_motif_celltype_enrichment(adata_merged)
    pd.testing.assert_frame_equal(
        from_dict.reset_index(drop=True), from_merged.reset_index(drop=True)
    )
    fig = plot_motif_celltype_enrichment(adata_dict, log2=True)
    assert fig.axes[0].get_ylabel() == "Motif"


# --------------------------------------------------------------------------
# Similarity
# --------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["cosine", "pearson", "spearman"])
def test_similarity_recovers_permutation(method):
    rng = np.random.default_rng(0)
    A = rng.gamma(1.0, 1.0, (5, 200))
    perm = [3, 0, 4, 1, 2]
    B = A[perm] + rng.normal(0, 0.01, A.shape)
    S = row_similarity(A, B, method=method)
    assert S.shape == (5, 5)
    # Row j of B is row perm[j] of A
    assert sorted(match_motifs(S)) == sorted((p, j) for j, p in enumerate(perm))


def test_similarity_plot_boxes_matches():
    S = np.eye(4) * 0.9 + 0.05
    fig = plot_motif_similarity(S, matches="hungarian", min_similarity=0.5)
    assert len(fig.axes[0].patches) == 4
    fig = plot_motif_similarity(S, matches=[(0, 1)], annotate_min=0.5)
    ax = fig.axes[0]
    assert len(ax.patches) == 1 and len(ax.texts) == 4


# --------------------------------------------------------------------------
# GSEA and QC
# --------------------------------------------------------------------------


def test_gsea_barplot_dataframe_and_result_object():
    df = pd.DataFrame(
        {
            "Term": [f"HALLMARK_SET_{i}" for i in range(6)],
            "NES": [2.1, 1.5, 0.4, -0.8, -1.9, -2.5],
            "FDR q-val": [0.01, 0.1, 0.9, 0.5, 0.02, 0.001],
        }
    )
    fig = gsea_barplot(df, n_top=2, pos_label="Up in B", neg_label="Up in A")
    ax = fig.axes[0]
    assert len(ax.patches) == 4
    assert [t.get_text() for t in ax.get_legend().get_texts()][:2] == [
        "Up in B",
        "Up in A",
    ]

    class Result:
        res2d = df

    assert len(gsea_barplot(Result(), n_top=1).axes[0].patches) == 2


def test_total_counts_per_cell(adata_single):
    fig = plot_total_counts_per_cell(adata_single, layer="counts", threshold=20)
    ax = fig.axes[0]
    expected = np.asarray(adata_single.layers["counts"]).sum(axis=1)
    heights = [p.get_height() for p in ax.patches]
    assert sum(heights) == len(expected)
    assert "layers['counts']" in ax.get_title()
