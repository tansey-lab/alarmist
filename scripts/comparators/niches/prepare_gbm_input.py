#!/usr/bin/env python
"""Export the GBM/LGG Xenium TMA as per-TMA-core NICHES inputs.

NICHES' ComputeEdgelist builds a DENSE n_cells x n_cells euclidean distance matrix
(R/ComputeEdgelist.R:36) -- at the full 100,197 cells that is ~80 GB, so the slide MUST be
split. obs['tma_id'] gives 13 physically separate tissue punches with disjoint spatial
extents, so a per-core split is also the only correct neighbour graph: cells in different
punches are not spatial neighbours. This matches how SpatialDM (13 cores) and COMMOT
(11/13 cores) were already run in this repo.

Writes, per core, the same input contract the cytosignal runner uses:
    counts.mtx      genes x cells, RAW integer counts (from layers['counts'])
    genes.tsv       one gene symbol per line   (rows of counts.mtx)
    barcodes.tsv    one cell id per line       (cols of counts.mtx)
    meta.csv        cell_id,x,y,celltype,grade,tma_id

RAW counts are exported on purpose: the NICHES vignettes all call Seurat::NormalizeData()
themselves, so normalisation must happen in R, not here. adata.X in this file is already
log-normalised (CLAUDE.md), which is why layers['counts'] is the layer we read.

Also emits the CellChatDB v2 human resource in NICHES' custom_LR_database format: a 2-column
data.frame whose first column is the ligand subunits and second the receptor subunits, each
'_'-separated (R/LoadCustom.R:19). Our CSV export already stores subunits '_'-joined, so this
is a straight column selection -- no complex re-encoding.

Env: bptf.  Usage:
    python prepare_gbm_input.py --h5ad data/xenium_mm_final_cell_id.h5ad \
        --cellchatdb data/LRdatabase/CellChatDBv2.0.human.csv \
        --out-dir results/comparators/niches/GBM/input
"""
import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cellchatdb", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-column", default="tma_id")
    p.add_argument("--cell-type-column", default="cell_type")
    p.add_argument("--condition-column", default="grade")
    p.add_argument("--count-layer", default="counts")
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    print(f"reading {a.h5ad}")
    adata = ad.read_h5ad(a.h5ad)
    print(adata)

    counts = adata.layers[a.count_layer]
    counts = sp.csr_matrix(counts) if not sp.issparse(counts) else counts.tocsr()
    dense_sample = np.asarray(counts[:200].todense())
    assert np.allclose(dense_sample % 1, 0), "layers['counts'] is not integer -- wrong layer?"

    xy = adata.obsm["spatial"]
    genes = adata.var_names.to_numpy().astype(str)

    # --- LR resource: CellChatDB v2 -> NICHES custom_LR_database (2 cols, '_'-separated) ---
    db = pd.read_csv(a.cellchatdb)
    lr = db[["ligand", "receptor"]].drop_duplicates().reset_index(drop=True)
    lr_path = os.path.join(a.out_dir, "custom_LR_database.csv")
    lr.to_csv(lr_path, index=False)
    panel = set(genes)

    def in_panel(row):
        return all(g in panel for g in row.ligand.split("_")) and all(
            g in panel for g in row.receptor.split("_")
        )

    n_panel = int(lr.apply(in_panel, axis=1).sum())
    print(
        f"wrote {lr_path}: {len(lr)} unique LR pairs from {len(db)} DB rows; "
        f"{n_panel} have every subunit on the {len(genes)}-gene panel "
        f"(NICHES::FilterGroundTruth keeps exactly these)"
    )

    # --- per-core export ---
    samples = adata.obs[a.sample_column].astype(str)
    index_rows = []
    for core in sorted(samples.unique(), key=lambda s: int(s)):
        mask = (samples == core).to_numpy()
        n = int(mask.sum())
        d = os.path.join(a.out_dir, f"core{core}")
        os.makedirs(d, exist_ok=True)

        sub = counts[mask]  # cells x genes
        # NICHES/Seurat want genes x cells
        sio.mmwrite(os.path.join(d, "counts.mtx"), sub.T.tocoo(), field="integer")
        with open(os.path.join(d, "genes.tsv"), "w") as fh:
            fh.write("\n".join(genes) + "\n")
        cell_ids = adata.obs_names[mask].to_numpy().astype(str)
        with open(os.path.join(d, "barcodes.tsv"), "w") as fh:
            fh.write("\n".join(cell_ids) + "\n")

        meta = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "x": xy[mask, 0],
                "y": xy[mask, 1],
                "celltype": adata.obs[a.cell_type_column].astype(str).to_numpy()[mask],
                "grade": adata.obs[a.condition_column].astype(str).to_numpy()[mask],
                "tma_id": core,
            }
        )
        meta.to_csv(os.path.join(d, "meta.csv"), index=False)

        grade = meta.grade.unique()
        assert len(grade) == 1, f"core {core} spans multiple grades: {grade}"
        index_rows.append(
            {
                "core": f"core{core}",
                "tma_id": core,
                "n_cells": n,
                "grade": grade[0],
                "n_cell_types": meta.celltype.nunique(),
                "x_span_um": round(float(np.ptp(xy[mask, 0])), 1),
                "y_span_um": round(float(np.ptp(xy[mask, 1])), 1),
            }
        )
        print(f"  core{core}: {n} cells, grade={grade[0]}, {meta.celltype.nunique()} cell types")

    idx = pd.DataFrame(index_rows).sort_values("n_cells").reset_index(drop=True)
    idx.to_csv(os.path.join(a.out_dir, "cores.csv"), index=False)

    with open(os.path.join(a.out_dir, "input_manifest.json"), "w") as fh:
        json.dump(
            {
                "h5ad": a.h5ad,
                "n_cells_total": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "count_layer": a.count_layer,
                "sample_column": a.sample_column,
                "cell_type_column": a.cell_type_column,
                "condition_column": a.condition_column,
                "coordinate_source": "obsm['spatial'] (microns)",
                "cellchatdb": a.cellchatdb,
                "n_lr_pairs": int(len(lr)),
                "n_lr_pairs_on_panel": n_panel,
                "cores": index_rows,
            },
            fh,
            indent=2,
        )
    print(f"\nwrote {len(index_rows)} cores + cores.csv + input_manifest.json to {a.out_dir}")


if __name__ == "__main__":
    main()
