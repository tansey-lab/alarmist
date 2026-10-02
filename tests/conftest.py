"""
Shared pytest fixtures: small synthetic ALARMIST outputs for plotting tests.
"""

import matplotlib

matplotlib.use("Agg")

import anndata  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

CELL_TYPES = ["Tumor", "Macrophage", "T cell", "Fibroblast"]
N_MOTIFS = 3


def _make_adata(n_cells: int, seed: int, sample: str) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                rng.choice(CELL_TYPES, n_cells), categories=CELL_TYPES
            ),
            "sample_id": sample,
        },
        index=[f"{sample}_c{i}" for i in range(n_cells)],
    )
    for k in range(N_MOTIFS):
        loading = rng.gamma(1.0, 1.0, n_cells)
        obs[f"motif_{k}_loading"] = loading
        obs[f"motif_{k}_state"] = pd.Categorical(
            np.where(loading > 1.5, "positive", "negative"),
            categories=["negative", "positive"],
        )
    X = rng.poisson(2.0, (n_cells, 10)).astype(np.float32)
    adata = anndata.AnnData(X=X, obs=obs)
    adata.layers["counts"] = X.copy()
    adata.obsm["spatial"] = rng.uniform(0, 500, (n_cells, 2))
    return adata


@pytest.fixture
def adata_single():
    return _make_adata(300, seed=0, sample="S1")


@pytest.fixture
def adata_dict():
    return {
        "S1": _make_adata(200, seed=1, sample="S1"),
        "S2": _make_adata(150, seed=2, sample="S2"),
    }


@pytest.fixture
def adata_merged(adata_dict):
    return anndata.concat(list(adata_dict.values()))


@pytest.fixture
def lri_motifs_df():
    """Long LRI × motif table as produced by process_bptf_results."""
    rows = []
    pairs = [
        ("GRN", "SORT1", "Secreted Signaling"),
        ("ANXA1", "FPR1", "Secreted Signaling"),
        ("CD99", "CD99", "Cell-Cell Contact"),
        ("SPP1", "CD44", "Secreted Signaling"),
        ("SPP1", "ITGB1", "Secreted Signaling"),
        ("COL1A1", "CD44", "ECM-Receptor"),
    ]
    rng = np.random.default_rng(3)
    for k in range(N_MOTIFS):
        for s in CELL_TYPES:
            for r in CELL_TYPES:
                for lig, rec, sig in pairs:
                    signaling = "autocrine" if s == r else sig
                    rows.append(
                        {
                            "motif_idx": k,
                            "celltype1": s,
                            "celltype2": r,
                            "ligand": lig,
                            "receptor": rec,
                            "signaling_type": signaling,
                            "factor": float(rng.gamma(1.0, 1.0)),
                        }
                    )
    return pd.DataFrame(rows)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")
