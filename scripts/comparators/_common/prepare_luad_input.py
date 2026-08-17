#!/usr/bin/env python
"""Build the shared comparator input for the LUAD/AIS four-section dataset.

WHY THIS EXISTS
---------------
The four source h5ads (``P{17,21}_{AIS,LUAD}_Xenium.h5ad``) are NOT directly runnable by
the comparator runners:

  1. There is no ``cell_type`` column (only ``annotation_coarse``), and no sample / stage /
     patient column at all -- the sample key exists solely as the filename.
  2. ``X`` is raw counts with NO ``layers`` and NO ``.raw``. SpatialDM needs log-scale ``.X``
     AND raw counts reachable as ``.raw``; ``run_spatialdm.py`` builds ``.raw`` from
     ``layers[<name>]`` and passing ``--count-layer X`` leaves ``.raw = None``, which dies
     inside ``spatialdm_local`` (``spatialdm/utils.py:181``) only AFTER the expensive global
     pass has already run.
  3. ``obs_names`` is a plain 0..n-1 integer index in every file, and 348 ``cell_id`` strings
     occur in two different sections -- a concatenation collides both ways.
  4. 17-23% of cells have a NaN ``annotation_coarse``. Five of the seven runners do
     ``.astype(str)`` on the cell-type column, which turns those into a literal ``'nan'``
     cell type that then enters every chord plot, every cell-type pair and every p-value.

So this script emits ONE canonical set of inputs that every method reads.

WHAT IT WRITES
--------------
  <out-dir>/P17_AIS.prepped.h5ad          per-section, for the runners with no split flag
  <out-dir>/P17_LUAD.prepped.h5ad         (stLearn, LIANA inflow, CytoSignal export)
  <out-dir>/P21_AIS.prepped.h5ad
  <out-dir>/P21_LUAD.prepped.h5ad
  <out-dir>/AIS_LUAD_4sections.h5ad       concatenated, for CellChat's prepare step
  <out-dir>/prep_manifest.json            provenance + the stLearn grid geometry

In every file: ``X`` = log1p(CP10K) float32, ``layers['counts']`` = raw counts float32,
``obsm['spatial']`` = the ORIGINAL per-section micron coordinates, unchanged.

THE CONCATENATED FILE IS FOR METADATA-LEVEL CONSUMERS ONLY. Its four coordinate frames
overlap numerically. Anything that builds a spatial graph over the whole object -- stLearn's
``grid()``, LIANA's ``spatial_neighbors`` -- MUST read the per-section files instead, or it
will fabricate neighbours between two patients. CellChat is the exception: it computes
distances strictly within ``meta$samples`` (``CellChat/R/modeling.R:1194-1228``), so the
overlap is inert there.

obs columns added
-----------------
  sample        P17_AIS | P17_LUAD | P21_AIS | P21_LUAD
  stage         AIS | LUAD                (the condition contrast)
  patient       P17 | P21                 (the only blocking factor available)
  cell_type     copy of annotation_coarse (19 levels, identical vocabulary in all four)
  orig_row      0-based row in the SOURCE h5ad, so any result can be mapped back
  alarmist_row  0-based row in ALARMIST's own cell ordering -- see below

``alarmist_row`` is what makes the comparator outputs joinable to the ALARMIST fit.
``results/AIS_LUAD/single_cell/cell_loadings.npy`` is (1676162, 25) and its row order is
"annotated cells in file order, sections in the order P17_AIS, P17_LUAD, P21_AIS, P21_LUAD".
Dropping the unannotated cells here reproduces exactly that ordering.

  *** Join cells by ``alarmist_row``. Do NOT join by LRI name. ***
  ``results/AIS_LUAD/`` was fitted with the PRE-2026-07-28 CellChatDB export; these
  comparator runs use the current one, and 1,120 of 3,218 LR keys changed between them
  (complex subunits were reordered, e.g. RAMP2_CALCR -> CALCR_RAMP2). See CLAUDE.md.

Env: bptf.  Read-only with respect to the source h5ads.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Section order is load-bearing: it defines alarmist_row. Do not reorder.
SECTIONS = ["P17_AIS", "P17_LUAD", "P21_AIS", "P21_LUAD"]

# stLearn's grid rule (scripts/comparators/stlearn/DEVIATIONS.md row 11): preserve the
# tutorial's SPOT AREA of 2,637 um^2 rather than its 125x125 bin count, which is specific
# to that tutorial's section extent. Equivalent square edge:
STLEARN_SPOT_AREA_UM2 = 2637.0
STLEARN_SPOT_EDGE_UM = float(np.sqrt(STLEARN_SPOT_AREA_UM2))  # 51.352 um


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def normalize_cp10k_log1p(X):
    """CP10K + log1p, done sparsely. Mirrors sc.pp.normalize_total(1e4) + sc.pp.log1p."""
    X = X.tocsr().astype(np.float32)
    libsize = np.asarray(X.sum(axis=1)).ravel()
    libsize[libsize == 0] = 1.0  # no annotated cell has zero counts, but be explicit
    inv = (1e4 / libsize).astype(np.float32)
    Xn = sp.diags(inv, dtype=np.float32) @ X
    Xn.data = np.log1p(Xn.data, dtype=np.float32)
    return Xn.tocsr()


def substring_collisions(levels):
    """stLearn maps cell types to deconvolution columns by SUBSTRING, first hit
    (run_stlearn.py:160-169 guards this and sys.exit()s when the mapping is wrong).
    In the LUAD vocabulary 'T' is a substring of 'Tumor_epi'. Report every such pair so a
    late abort inside stLearn is never a surprise."""
    out = []
    for a in levels:
        for b in levels:
            if a != b and a in b:
                out.append({"needle": a, "haystack": b})
    return out


def prep_one(src, section, celltype_col, drop_unannotated):
    log(f"--- {section}: reading {src}")
    adata = ad.read_h5ad(src)
    n_src = adata.n_obs

    if celltype_col not in adata.obs:
        sys.exit(f"ERROR: '{celltype_col}' not in obs of {src}; have {list(adata.obs.columns)}")

    adata.obs["orig_row"] = np.arange(n_src, dtype=np.int32)

    ct = adata.obs[celltype_col]
    keep = ct.notna().to_numpy()
    n_drop = int((~keep).sum())
    if drop_unannotated:
        adata = adata[keep].copy()
        log(f"    dropped {n_drop} unannotated cells ({100 * n_drop / n_src:.1f}%) -> {adata.n_obs}")
    elif n_drop:
        log(f"    WARNING: keeping {n_drop} unannotated cells; five runners will render them "
            f"as a literal 'nan' cell type")

    patient, stage = section.split("_", 1)
    adata.obs["sample"] = pd.Categorical([section] * adata.n_obs, categories=SECTIONS)
    adata.obs["stage"] = pd.Categorical([stage] * adata.n_obs, categories=["AIS", "LUAD"])
    adata.obs["patient"] = pd.Categorical([patient] * adata.n_obs, categories=["P17", "P21"])
    # Keep the source category order so the 19 levels are identical across sections --
    # CellChat's functional-similarity comparison requires identical composition.
    adata.obs["cell_type"] = adata.obs[celltype_col].astype("category")

    # obs_names must be unique across sections. Overwriting cell_id is load-bearing:
    # the CytoSignal/NICHES exporters PREFER obs['cell_id'] over obs_names, so leaving the
    # bare Xenium id there would reintroduce the 348 cross-section duplicates in
    # barcodes.tsv and trip the stopifnot() in run_cytosignal.R / run_cellchat.R.
    base = adata.obs["cell_id"].astype(str) if "cell_id" in adata.obs else pd.Series(
        adata.obs_names.astype(str), index=adata.obs.index)
    adata.obs_names = pd.Index(f"{section}_" + base.to_numpy().astype(str))
    adata.obs["cell_id"] = adata.obs_names.astype(str)
    assert adata.obs_names.is_unique, f"{section}: obs_names not unique after prefixing"

    # X: raw counts in, log1p(CP10K) out; the raw counts survive in layers['counts'].
    X = adata.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    X.eliminate_zeros()
    frac_int = float(np.mean(np.mod(X.data[: min(X.data.size, 100000)], 1) == 0))
    if frac_int < 0.999:
        sys.exit(f"ERROR: {section}: X does not look like raw counts "
                 f"(only {frac_int:.3%} of a 100k-value sample are integers)")
    adata.layers["counts"] = X.astype(np.float32)
    adata.X = normalize_cp10k_log1p(X)

    # SpatialData leftovers that mean nothing here and confuse downstream readers.
    if "region" in adata.obs:
        del adata.obs["region"]
    adata.uns.pop("spatialdata_attrs", None)
    adata.raw = None

    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    x_ext = float(xy[:, 0].max() - xy[:, 0].min())
    y_ext = float(xy[:, 1].max() - xy[:, 1].min())
    # Computed on the ANNOTATED subset on purpose: run_stlearn.py drops unannotated cells
    # (:98-99) BEFORE it reads xs/ys (:122), so the grid it builds is over this extent.
    n_col = int(round(x_ext / STLEARN_SPOT_EDGE_UM))
    n_row = int(round(y_ext / STLEARN_SPOT_EDGE_UM))

    info = {
        "section": section,
        "source": os.path.abspath(src),
        "source_md5": md5(src),
        "n_cells_source": int(n_src),
        "n_cells_unannotated_dropped": n_drop if drop_unannotated else 0,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "nnz": int(adata.layers["counts"].nnz),
        "stage": stage,
        "patient": patient,
        "cell_types": [str(c) for c in adata.obs["cell_type"].cat.categories],
        "extent_um": {"x": round(x_ext, 1), "y": round(y_ext, 1)},
        "stlearn_grid": {
            "spot_area_um2": STLEARN_SPOT_AREA_UM2,
            "spot_edge_um": round(STLEARN_SPOT_EDGE_UM, 3),
            "n_col": n_col,
            "n_row": n_row,
            "realised_spot_x_um": round(x_ext / n_col, 2),
            "realised_spot_y_um": round(y_ext / n_row, 2),
        },
    }
    log(f"    {adata.n_obs} cells x {adata.n_vars} genes | extent {x_ext:.0f} x {y_ext:.0f} um "
        f"| stLearn grid --n-col {n_col} --n-row {n_row}")
    return adata, info


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", required=True, help="directory holding P{17,21}_{AIS,LUAD}_Xenium.h5ad")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cell-type-column", default="annotation_coarse")
    p.add_argument("--keep-unannotated", action="store_true",
                   help="do NOT drop cells with a NaN cell type (default is to drop; see docstring)")
    p.add_argument("--no-concat", action="store_true",
                   help="skip the concatenated file (CellChat needs it; nothing else does)")
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    drop = not a.keep_unannotated

    infos, out_paths = [], []
    running_offset = 0
    for s in SECTIONS:
        src = os.path.join(a.in_dir, f"{s}_Xenium.h5ad")
        if not os.path.exists(src):
            sys.exit(f"ERROR: {src} not found. Set --in-dir to the directory holding the four "
                     f"P{{17,21}}_{{AIS,LUAD}}_Xenium.h5ad files.")
        adata, info = prep_one(src, s, a.cell_type_column, drop)

        # alarmist_row: position in the concatenated ANNOTATED stream, sections in SECTIONS
        # order. Reproduces the row order of results/AIS_LUAD/single_cell/cell_loadings.npy.
        adata.obs["alarmist_row"] = np.arange(
            running_offset, running_offset + adata.n_obs, dtype=np.int32)
        info["alarmist_row_range"] = [running_offset, running_offset + adata.n_obs - 1]
        running_offset += adata.n_obs

        out = os.path.join(a.out_dir, f"{s}.prepped.h5ad")
        log(f"    writing {out}")
        adata.write_h5ad(out, compression="gzip")
        info["prepped"] = os.path.abspath(out)
        infos.append(info)
        out_paths.append(out)
        del adata

    ct_sets = [set(i["cell_types"]) for i in infos]
    identical_ct = all(c == ct_sets[0] for c in ct_sets)
    collisions = substring_collisions(sorted(ct_sets[0]))
    log(f"identical cell-type vocabulary across the four sections: {identical_ct}")
    if collisions:
        log("cell_type SUBSTRING COLLISIONS (stLearn matches deconvolution columns by "
            "substring, first hit):")
        for c in collisions:
            log(f"    '{c['needle']}' is a substring of '{c['haystack']}'")
        log("    -> if run_stlearn.py aborts on the cell-type mapping guard, that is why.")

    concat_path = None
    if not a.no_concat:
        concat_path = os.path.join(a.out_dir, "AIS_LUAD_4sections.h5ad")
        log(f"concatenating -> {concat_path}")
        try:
            from anndata.experimental import concat_on_disk
            concat_on_disk(out_paths, concat_path, label=None, index_unique=None)
            log("    used anndata.experimental.concat_on_disk (low memory)")
        except Exception as e:  # noqa: BLE001 -- fall back rather than fail the whole prep
            log(f"    concat_on_disk unavailable/failed ({e}); falling back to in-memory concat")
            merged = ad.concat([ad.read_h5ad(p_) for p_ in out_paths], axis=0,
                               join="outer", index_unique=None)
            merged.write_h5ad(concat_path, compression="gzip")
            del merged
        log(f"    total {running_offset} cells")

    manifest = {
        "generated_by": os.path.abspath(__file__),
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "in_dir": os.path.abspath(a.in_dir),
        "out_dir": os.path.abspath(a.out_dir),
        "section_order": SECTIONS,
        "cell_type_column": a.cell_type_column,
        "unannotated_dropped": drop,
        "n_cells_total": running_offset,
        "X": "log1p(counts / libsize * 1e4), float32",
        "layers_counts": "raw integer counts, float32",
        "spatial": "obsm['spatial'], original per-section microns, UNCHANGED",
        "obs_added": ["sample", "stage", "patient", "cell_type", "orig_row", "alarmist_row"],
        "identical_celltype_vocabulary": identical_ct,
        "cell_type_substring_collisions": collisions,
        "concatenated": os.path.abspath(concat_path) if concat_path else None,
        "concatenated_warning": (
            "The four coordinate frames OVERLAP numerically. Use this file only for "
            "metadata-level consumers (CellChat's prepare step, which splits by "
            "meta$samples before computing any distance). Never build a single spatial "
            "graph over it."
        ),
        "alarmist_join": (
            "Join cells to results/AIS_LUAD/single_cell/cell_loadings.npy by obs['alarmist_row']. "
            "Do NOT join by LRI name: that ALARMIST fit used the pre-2026-07-28 CellChatDB "
            "export and 1,120 of 3,218 LR keys differ from the DB these comparator runs use."
        ),
        "sections": infos,
    }
    mpath = os.path.join(a.out_dir, "prep_manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)

    log(f"\nwrote {len(infos)} per-section files"
        f"{' + 1 concatenated' if concat_path else ''} + prep_manifest.json to {a.out_dir}")
    log("stLearn grid geometry (pass these to run_stlearn.py):")
    for i in infos:
        g = i["stlearn_grid"]
        log(f"    {i['section']:10} --n-col {g['n_col']:4d} --n-row {g['n_row']:4d}"
            f"   ({g['realised_spot_x_um']} x {g['realised_spot_y_um']} um spots)")


if __name__ == "__main__":
    main()
