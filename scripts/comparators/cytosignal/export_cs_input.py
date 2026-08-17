#!/usr/bin/env python
"""Export ONE prepped LUAD/AIS section to the CytoSignal input contract.

Replaces ``results/comparators/cytosignal/bundle_bignode/export_p21_full.py``, which lives
in the OUTPUT tree (against SKILL.md's code/output separation) and, more importantly, reads
``adata.X``. That was correct when ``X`` held raw counts; after
``_common/prepare_luad_input.py`` runs, ``X`` is log1p(CP10K) and the raw counts live in
``layers['counts']``. Reading ``X`` here would silently truncate every count to 0 or 1 when
cast to int32.

  *** THIS SCRIPT READS layers['counts'], NEVER X. ***

Output contract, consumed by ``run_cytosignal.R`` (its lines 25-31) and by
``run_nebula_stage.R``:

    counts.mtx       genes x cells, MatrixMarket, field=integer
    genes.tsv        one gene symbol per line   (rows of counts.mtx)
    barcodes.tsv     one cell id per line       (cols of counts.mtx)
    meta.csv         cell_id,x,y,celltype,sample,stage,patient,alarmist_row
    provenance.json

``run_cytosignal.R`` reads only ``cell_id,x,y,celltype`` by name; the extra columns are
inert there and are what ``run_nebula_stage.R`` uses to build its per-section metadata.

Coordinates are already in microns (Xenium), so downstream ``inferEpsParams`` uses
``scale.factor = 1``.

Env: bptf.  Usage:
    python export_cs_input.py <prepped_section.h5ad> <out_dir>
"""

from __future__ import annotations

import collections
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

LAYER = "counts"


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)

    a = ad.read_h5ad(src)
    if LAYER not in a.layers:
        sys.exit(f"ERROR: layers['{LAYER}'] missing from {src}; have {list(a.layers)}. "
                 f"Run scripts/comparators/_common/prepare_luad_input.py first.")
    for col in ("cell_type", "sample", "stage", "patient"):
        if col not in a.obs:
            sys.exit(f"ERROR: obs['{col}'] missing from {src}; prep was not run or is stale.")

    X = sp.csr_matrix(a.layers[LAYER]).astype(np.float64)
    X.eliminate_zeros()

    # CytoSignal wants integer counts. Assert rather than assume -- a silent float->int cast
    # of a log-normalized matrix is exactly the failure this script exists to prevent.
    sample = X.data[: min(X.data.size, 200000)]
    if not np.all(np.mod(sample, 1) == 0):
        sys.exit(f"ERROR: layers['{LAYER}'] is not integral -- refusing to cast to int32. "
                 f"Check that prepare_luad_input.py wrote raw counts there.")

    # Drop zero-count cells (prep already dropped unannotated ones; on this data the two
    # filters coincide, but run_cytosignal.R's own QC assumes no all-zero column).
    nz = np.asarray(X.sum(axis=1)).ravel() > 0
    if not nz.all():
        print(f"dropping {int((~nz).sum())} zero-count cells")
        a = a[nz].copy()
        X = sp.csr_matrix(a.layers[LAYER]).astype(np.float64)
        X.eliminate_zeros()

    genes = np.asarray(a.var_names, dtype=str)
    cells = np.asarray(a.obs["cell_id"].astype(str))
    if len(set(cells)) != len(cells):
        sys.exit("ERROR: cell ids are not unique; prepare_luad_input.py should have prefixed them.")
    coords = np.asarray(a.obsm["spatial"], dtype=float)
    cts = a.obs["cell_type"].astype(str).to_numpy()

    print(f"{os.path.basename(src)}: {a.n_obs} cells x {a.n_vars} genes, nnz {X.nnz}")
    print(collections.Counter(cts).most_common())

    sio.mmwrite(os.path.join(out, "counts.mtx"), X.T.tocsc().astype(np.int32), field="integer")
    pd.Series(genes).to_csv(os.path.join(out, "genes.tsv"), index=False, header=False)
    pd.Series(cells).to_csv(os.path.join(out, "barcodes.tsv"), index=False, header=False)
    pd.DataFrame(
        {
            "cell_id": cells,
            "x": coords[:, 0],
            "y": coords[:, 1],
            "celltype": cts,
            "sample": a.obs["sample"].astype(str).to_numpy(),
            "stage": a.obs["stage"].astype(str).to_numpy(),
            "patient": a.obs["patient"].astype(str).to_numpy(),
            "alarmist_row": a.obs["alarmist_row"].to_numpy(),
        }
    ).to_csv(os.path.join(out, "meta.csv"), index=False)

    json.dump(
        {
            "source": os.path.abspath(src),
            "count_source": f"layers['{LAYER}'] (raw integer counts) -- NOT X",
            "sample": str(a.obs["sample"].iloc[0]),
            "stage": str(a.obs["stage"].iloc[0]),
            "patient": str(a.obs["patient"].iloc[0]),
            "n_cells": int(a.n_obs),
            "n_genes": int(a.n_vars),
            "nnz": int(X.nnz),
            "coordinate_units": "microns (Xenium) -> scale.factor = 1",
        },
        open(os.path.join(out, "provenance.json"), "w"),
        indent=2,
    )
    print("wrote", out)


if __name__ == "__main__":
    main()
