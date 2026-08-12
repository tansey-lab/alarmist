#!/usr/bin/env python
"""Export the GBM/LGG Xenium TMA as CellChat inputs — one directory per CONDITION.

Unlike the NICHES/CytoSignal exports (one directory per TMA core), CellChat wants all
sections of ONE condition inside ONE object, distinguished by a `samples` column, and a
SEPARATE object per condition:

    "a column named `samples` should be provided for spatial transcriptomics analysis,
     which is useful for analyzing cell-cell communication by aggregating multiple
     samples/replicates. Of note, for comparison analysis across different conditions,
     users still need to create a CellChat object separately for each condition."
        -- CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.Rmd:82
           and FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.Rmd:27

So this writes `high/` (7 cores) and `low/` (6 cores), each carrying every core's cells with
`samples` = tma_id. CellChat then computes cell-group distances PER SAMPLE and averages them
(modeling.R:1194-1231), which is what keeps physically separate punches from being treated as
spatial neighbours.

NORMALIZED data is exported, not raw counts: CellChat requires "Normalized data (e.g.,
library-size normalization and then log-transformed with a pseudocount of 1)" as
`data.input` (V1:49) and the vignettes feed it Seurat's `slot="data"`. adata.X in this file
is already log-normalized (CLAUDE.md), so X is exported as-is and the normalization is
ASSERTED here rather than assumed — see check_normalization().

Per condition:
    data.mtx            genes x cells, log-normalized (from adata.X)
    genes.tsv           one gene symbol per line   (rows of data.mtx)
    barcodes.tsv        one cell id per line       (cols of data.mtx)
    meta.csv            cell_id,labels,samples,x,y,grade   (labels = cell_type)
    spatial_factors.csv sample,ratio,tol -- ONE ROW PER SAMPLE, in levels() order

`spatial_factors.csv` row order matters: computeRegionDistance indexes ratio[k]/tol[k]
positionally over levels(meta$samples) (modeling.R:1165,1212).

Xenium spatial factors are the FAQ's own values for single-cell-resolution platforms:
ratio = 1 (coordinates already in microns), tol = spot.size/2 = 10/2 = 5, where spot.size is
"the typical human cell size" (FAQ:78-83).

Env: bptf.  Usage:
    python prepare_gbm_input.py --h5ad data/xenium_mm_final_cell_id.h5ad \
        --out-dir results/comparators/cellchat/GBM/input
"""
import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def check_normalization(X, n=500):
    """Assert adata.X is log1p(library-size-normalized), and report the scale factor.

    CellChat's own normalizeData() is `log1p(x / colSums(x) * scale.factor)`, so feeding it
    raw counts, or data that has been logged twice, silently changes every probability.
    """
    sub = X[:n]
    dense = np.asarray(sub.todense()) if sp.issparse(sub) else np.asarray(sub)
    is_integer = np.allclose(dense % 1, 0)
    if is_integer:
        raise SystemExit(
            "adata.X looks like RAW INTEGER COUNTS. CellChat needs normalized data "
            "(V1:49). Point --count-layer at the log-normalized layer, or normalize first."
        )
    # invert log1p and look at the per-cell library size it implies
    lib = np.expm1(dense).sum(axis=1)
    return {
        "x_max": float(dense.max()),
        "implied_libsize_median": float(np.median(lib)),
        "implied_libsize_iqr": [float(np.percentile(lib, 25)), float(np.percentile(lib, 75))],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-column", default="tma_id")
    p.add_argument("--cell-type-column", default="cell_type")
    p.add_argument("--condition-column", default="grade")
    p.add_argument("--ratio", type=float, default=1.0, help="pixels->microns; Xenium = 1")
    p.add_argument("--spot-size", type=float, default=10.0, help="typical human cell size (um)")
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    print(f"reading {a.h5ad}")
    adata = ad.read_h5ad(a.h5ad)
    print(adata)

    X = adata.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    norm = check_normalization(X)
    print(f"  X normalization check: {norm}")

    xy = adata.obsm["spatial"]
    genes = adata.var_names.to_numpy().astype(str)
    samples = adata.obs[a.sample_column].astype(str)
    labels = adata.obs[a.cell_type_column].astype(str)
    grade = adata.obs[a.condition_column].astype(str)

    tol = a.spot_size / 2.0
    conditions = sorted(grade.unique())
    cond_rows, core_rows = [], []

    for cond in conditions:
        mask = (grade == cond).to_numpy()
        d = os.path.join(a.out_dir, cond)
        os.makedirs(d, exist_ok=True)

        sub = X[mask]  # cells x genes
        sio.mmwrite(os.path.join(d, "data.mtx"), sub.T.tocoo())  # genes x cells
        with open(os.path.join(d, "genes.tsv"), "w") as fh:
            fh.write("\n".join(genes) + "\n")
        cell_ids = adata.obs_names[mask].to_numpy().astype(str)
        with open(os.path.join(d, "barcodes.tsv"), "w") as fh:
            fh.write("\n".join(cell_ids) + "\n")

        meta = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "labels": labels.to_numpy()[mask],
                "samples": samples.to_numpy()[mask],
                "x": xy[mask, 0],
                "y": xy[mask, 1],
                "grade": cond,
            }
        )
        meta.to_csv(os.path.join(d, "meta.csv"), index=False)

        # one row per sample, in the SAME order R's levels(factor(samples)) will produce.
        # R sorts character factor levels with the C locale; sort the strings the same way.
        sample_levels = sorted(meta["samples"].unique())
        pd.DataFrame(
            {"sample": sample_levels, "ratio": a.ratio, "tol": tol}
        ).to_csv(os.path.join(d, "spatial_factors.csv"), index=False)

        for s in sample_levels:
            m = meta["samples"] == s
            core_rows.append(
                {
                    "tma_id": s,
                    "grade": cond,
                    "n_cells": int(m.sum()),
                    "n_cell_types": int(meta.loc[m, "labels"].nunique()),
                    "x_span_um": round(float(np.ptp(meta.loc[m, "x"])), 1),
                    "y_span_um": round(float(np.ptp(meta.loc[m, "y"])), 1),
                }
            )
        cond_rows.append(
            {
                "condition": cond,
                "n_cells": int(mask.sum()),
                "n_samples": len(sample_levels),
                "samples": sample_levels,
                "n_cell_types": int(meta["labels"].nunique()),
                "cell_types": sorted(meta["labels"].unique()),
            }
        )
        print(
            f"  {cond}: {mask.sum()} cells, {len(sample_levels)} cores "
            f"({', '.join(sample_levels)}), {meta['labels'].nunique()} cell types"
        )

    # CellChat's functional-similarity comparison requires identical cell-type composition
    ct_sets = [set(c["cell_types"]) for c in cond_rows]
    same_composition = all(s == ct_sets[0] for s in ct_sets)
    print(f"  identical cell-type composition across conditions: {same_composition}")

    pd.DataFrame(core_rows).to_csv(os.path.join(a.out_dir, "cores.csv"), index=False)
    with open(os.path.join(a.out_dir, "input_manifest.json"), "w") as fh:
        json.dump(
            {
                "h5ad": a.h5ad,
                "n_cells_total": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "data_source": "adata.X (log-normalized; raw counts are in layers['counts'])",
                "normalization_check": norm,
                "sample_column": a.sample_column,
                "cell_type_column": a.cell_type_column,
                "condition_column": a.condition_column,
                "coordinate_source": "obsm['spatial'] (microns)",
                "spatial_factors": {"ratio": a.ratio, "tol": tol, "spot_size_um": a.spot_size},
                "identical_celltype_composition": same_composition,
                "conditions": cond_rows,
                "cores": core_rows,
            },
            fh,
            indent=2,
        )
    print(f"\nwrote {len(cond_rows)} conditions + cores.csv + input_manifest.json to {a.out_dir}")


if __name__ == "__main__":
    main()
