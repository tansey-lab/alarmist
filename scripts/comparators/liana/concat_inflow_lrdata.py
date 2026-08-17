#!/usr/bin/env python
"""Concatenate the four per-section LIANA+ ``inflow`` outputs into ONE object for a joint
MOFA-Flex fit, with the coordinate frames laid out side by side.

WHY
---
``run_inflow.py`` is per-section (its ``li.ut.spatial_neighbors`` graph is built over the
whole object, so four superimposed coordinate frames would fabricate neighbours between two
patients). That leaves four independent factorisations whose factor *k* have no relation to
each other. MOFA-Flex is the one LIANA branch that can put all four sections into a single
comparable factor space -- ``mofaflex`` shares factors across the whole fit -- so it is the
only place a cross-section comparison is available at all for LIANA (open deviation CD-2:
LIANA has no native multi-sample mode for the spatial per-pair branches).

THE COORDINATE TRAP THIS SCRIPT EXISTS TO CLOSE
-----------------------------------------------
``run_mofaflex.py:433-434`` builds ``mfl.priors.GaussianProcess(covariates_mkey=spatial_key)``
-- the factor prior is a Gaussian process **over the spatial coordinates**. The four sections
all span roughly 0-11,000 x 0-13,000 um, so a naive concatenation would place four different
tissues on top of each other and the GP would model them as one continuous field.

``liana.utils.expand_coordinates`` is LIANA's own answer to exactly this. Its 1.8.1 CHANGELOG
entry names the use case verbatim -- "enabling multi-sample spatial analyses (e.g. a joint
`spatial_neighbors` graph)" -- and the transform is a pure per-sample translation onto a grid,
so within-sample distances are preserved exactly and the originals are kept in
``obsm['spatial_original']``.

  *** RECORD THIS AS A DEVIATION. *** No LIANA tutorial calls ``expand_coordinates``
  (verified: zero hits across docs/notebooks/). It is an author-provided utility used for its
  documented purpose, not a tutorial step. The alternative -- four separate MOFA-Flex fits --
  yields four incomparable factor spaces and answers nothing.

FEATURE JOIN
------------
Inflow feature sets differ per section because ``nz_prop`` (and the optional SVI filter) are
evaluated per object. The default here is an INNER join: only features present in all four
sections enter the joint fit. An outer join would zero-fill absent features, and MOFA-Flex
cannot distinguish a structural zero from a measured one -- it would read "this interaction
was filtered out of P17_AIS" as "this interaction is off in P17_AIS". The dropped features
are written out, never silently discarded.

Env: comp-liana.  Usage:
    python concat_inflow_lrdata.py --in-dirs <d1> <d2> <d3> <d4> --out <joint_lrdata.h5ad>
where each <d> is a run_inflow.py output dir (the file read is <d>/data/inflow_lrdata.h5ad).
"""

from __future__ import annotations

import argparse
import json
import os
import time

import anndata as ad
import numpy as np
import pandas as pd

import liana as li


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dirs", nargs="+", required=True,
                   help="run_inflow.py output dirs, one per section, in section order")
    p.add_argument("--out", required=True, help="joint lrdata .h5ad for run_mofaflex.py --lrdata")
    p.add_argument("--join", choices=["inner", "outer"], default="inner")
    p.add_argument("--sample-key", default="sample")
    p.add_argument("--spatial-key", default="spatial")
    p.add_argument("--margin", type=float, default=0.1,
                   help="expand_coordinates grid margin, fraction of the largest extent")
    p.add_argument("--n-cols", type=int, default=2,
                   help="expand_coordinates grid columns; 2 gives a 2x2 layout for 4 sections")
    a = p.parse_args()

    parts, infos = [], []
    for d in a.in_dirs:
        f = os.path.join(d, "data", "inflow_lrdata.h5ad")
        if not os.path.exists(f):
            f2 = d if d.endswith(".h5ad") else None
            if f2 and os.path.exists(f2):
                f = f2
            else:
                raise SystemExit(f"ERROR: {f} not found (and {d} is not an .h5ad either)")
        log(f"reading {f}")
        x = ad.read_h5ad(f)
        if a.sample_key not in x.obs:
            raise SystemExit(
                f"ERROR: obs['{a.sample_key}'] missing from {f}. run_inflow.py must have been "
                f"given a prepped h5ad from _common/prepare_luad_input.py.")
        samp = str(pd.unique(x.obs[a.sample_key].astype(str))[0])
        log(f"    {samp}: {x.shape[0]} cells x {x.shape[1]} features")
        infos.append({"dir": os.path.abspath(d), "sample": samp,
                      "n_cells": int(x.shape[0]), "n_features": int(x.shape[1])})
        parts.append(x)

    feat_sets = [set(map(str, x.var_names)) for x in parts]
    inter = set.intersection(*feat_sets)
    union = set.union(*feat_sets)
    log(f"features: union {len(union)}, intersection {len(inter)}, join={a.join}")

    joint = ad.concat(parts, axis=0, join=a.join, index_unique=None, merge="first",
                      uns_merge="first")
    log(f"concatenated: {joint.shape[0]} cells x {joint.shape[1]} features")
    assert joint.obs_names.is_unique, "obs_names collided; prepped ids should be prefixed"

    # Make sample a proper categorical in section order so expand_coordinates lays the grid
    # out deterministically.
    order = [i["sample"] for i in infos]
    joint.obs[a.sample_key] = pd.Categorical(joint.obs[a.sample_key].astype(str), categories=order)

    log(f"li.ut.expand_coordinates(sample_key='{a.sample_key}', n_cols={a.n_cols}, "
        f"margin={a.margin})")
    before = np.asarray(joint.obsm[a.spatial_key], dtype=float)
    joint = li.ut.expand_coordinates(joint, sample_key=a.sample_key,
                                     spatial_key=a.spatial_key,
                                     n_cols=a.n_cols, margin=a.margin)
    after = np.asarray(joint.obsm[a.spatial_key], dtype=float)
    assert f"{a.spatial_key}_original" in joint.obsm, "expand_coordinates did not keep originals"

    # Prove the translation is per-sample and rigid: within a sample the coordinate spread is
    # unchanged, and across samples the bounding boxes are now disjoint.
    boxes = {}
    for s in order:
        m = (joint.obs[a.sample_key].astype(str) == s).to_numpy()
        b0, b1 = before[m], after[m]
        spread_ok = np.allclose(b0.max(0) - b0.min(0), b1.max(0) - b1.min(0), atol=1e-6)
        boxes[s] = {"xmin": float(b1[:, 0].min()), "xmax": float(b1[:, 0].max()),
                    "ymin": float(b1[:, 1].min()), "ymax": float(b1[:, 1].max()),
                    "extent_preserved": bool(spread_ok)}
        log(f"    {s}: x[{boxes[s]['xmin']:.0f},{boxes[s]['xmax']:.0f}] "
            f"y[{boxes[s]['ymin']:.0f},{boxes[s]['ymax']:.0f}] "
            f"extent preserved: {spread_ok}")
    overlaps = []
    for i, s in enumerate(order):
        for t in order[i + 1:]:
            A, B = boxes[s], boxes[t]
            if not (A["xmax"] < B["xmin"] or B["xmax"] < A["xmin"]
                    or A["ymax"] < B["ymin"] or B["ymax"] < A["ymin"]):
                overlaps.append([s, t])
    log(f"    overlapping bounding-box pairs after expansion: {len(overlaps)} {overlaps}")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    joint.write_h5ad(a.out, compression="gzip")
    log(f"wrote {a.out}")

    man = {
        "generated_by": os.path.abspath(__file__),
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "inputs": infos,
        "join": a.join,
        "n_features_union": len(union),
        "n_features_intersection": len(inter),
        "n_features_used": int(joint.shape[1]),
        "features_dropped": sorted(union - inter) if a.join == "inner" else [],
        "n_cells": int(joint.shape[0]),
        "expand_coordinates": {"function": "liana.utils.expand_coordinates",
                               "sample_key": a.sample_key, "n_cols": a.n_cols,
                               "margin": a.margin,
                               "originals_kept_in": f"obsm['{a.spatial_key}_original']",
                               "post_expansion_boxes": boxes,
                               "overlapping_box_pairs": overlaps},
        "deviation": ("expand_coordinates is a liana 1.8.1 utility used for its documented "
                      "purpose, but NO liana tutorial calls it. Log this in "
                      "liana/DEVIATIONS.md before reporting any joint MOFA-Flex result."),
    }
    with open(os.path.splitext(a.out)[0] + "_concat_manifest.json", "w") as fh:
        json.dump(man, fh, indent=2)
    log("wrote concat manifest")


if __name__ == "__main__":
    main()
