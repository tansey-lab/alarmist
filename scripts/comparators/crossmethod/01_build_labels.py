"""Build the per-section label + ALARMIST-feature table for the 1-D classifier benchmark.

Reads ONLY `obs` and `obsm['spatial']` from the LUAD h5ad -- nothing is
recomputed, nothing is written back.  Produces one .npz per section holding the
join key, the two binarised pathologist labels, the two ALARMIST motif
loadings, and the reference features.

JOIN KEY.  `obs_names` in this file are literally '0', '1', '2', ... and carry
no identity (verified 2026-09-10), so the key is `sample_id + '_' + cell_id`,
which reproduces the comparator cell ids exactly (182,378/182,378 matched for
P17_AIS, zero on either side unmatched).  Never use obs_names here.

REFERENCE FEATURES.  Beyond the comparators, three features are scored in the
same 1-D framework so the AUCs can be read against something:
  * `local_density_80um` -- number of cells within 80 um (ALARMIST's patch
    size).  This is the "tissue architecture alone" control: any spatially
    smoothed feature inherits some of this, and a reviewer will ask.
  * `total_counts` -- per-cell transcript depth.
  * cell-type indicators are built at scoring time from `cell_type`.
These are NOT part of the meeting takeaway; they are labelled `reference` in
the output and can be dropped without touching anything else.

Usage:
    python 01_build_labels.py [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import ADATA, ALARMIST_FEATURES, ALL_MOTIF_FEATURES, LABEL_DEFS, OUT, \
    SECTIONS, TISSUE_LAYERS

DENSITY_RADIUS_UM = 80.0   # = ALARMIST's LUAD patch size (CLAUDE.md Fig 4)


def read_obs_column(f: h5py.File, name: str) -> np.ndarray:
    """Read one obs column, decoding anndata's categorical encoding."""
    o = f['obs'][name]
    if isinstance(o, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in o['categories'][:]])
        return cats[o['codes'][:]]
    return o[:]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=OUT)
    args = ap.parse_args()
    out_dir = args.out / 'labels'
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f'reading obs from {ADATA}')
    with h5py.File(ADATA, 'r') as f:
        sample_id = read_obs_column(f, 'sample_id').astype(str)
        cell_id = read_obs_column(f, 'cell_id').astype(str)
        tissue = read_obs_column(f, 'tissue_layer').astype(str)
        cell_type = read_obs_column(f, 'cell_type').astype(str)
        total_counts = read_obs_column(f, 'total_counts').astype(np.float64)
        # ALARMIST_FEATURES (motifs 24 and 10) are the claim; ALL_MOTIF_FEATURES
        # is the all-25 control.  The two dicts overlap by value, not by name.
        motifs = {name: read_obs_column(f, col).astype(np.float64)
                  for name, col in {**ALARMIST_FEATURES, **ALL_MOTIF_FEATURES}.items()}
        xy = f['obsm']['spatial'][:].astype(np.float64)
    key = np.char.add(np.char.add(sample_id, '_'), cell_id)
    print(f'  {len(key):,} cells, {time.time() - t0:.0f}s')

    observed = sorted(set(tissue))
    unexpected = [t for t in observed if t not in TISSUE_LAYERS]
    if unexpected:
        raise SystemExit(f'unexpected tissue_layer levels on disk: {unexpected}')

    from scipy.spatial import cKDTree

    manifest = {'adata': str(ADATA), 'n_cells_total': int(len(key)),
                'density_radius_um': DENSITY_RADIUS_UM, 'sections': {}}

    for sec in SECTIONS:
        m = sample_id == sec
        n = int(m.sum())
        tl = tissue[m]

        labels = {}
        for lname, d in LABEL_DEFS.items():
            lab = np.full(n, -1, dtype=np.int8)          # -1 = excluded
            lab[np.isin(tl, d['pos'])] = 1
            lab[np.isin(tl, d['neg'])] = 0
            labels[lname] = lab

        # Local density within 80 um -- self included, so the minimum is 1.
        t1 = time.time()
        pts = xy[m]
        tree = cKDTree(pts)
        density = np.asarray(
            tree.query_ball_point(pts, r=DENSITY_RADIUS_UM, return_length=True),
            dtype=np.float64)

        np.savez_compressed(
            out_dir / f'{sec}.npz',
            key=key[m], tissue_layer=tl, cell_type=cell_type[m],
            x=pts[:, 0], y=pts[:, 1],
            total_counts=total_counts[m], local_density_80um=density,
            **{f'label_{k}': v for k, v in labels.items()},
            **{name: motifs[name][m] for name in motifs},
        )

        sec_info = {'n_cells': n,
                    'tissue_layer_counts': {t: int((tl == t).sum()) for t in TISSUE_LAYERS},
                    'density_median': float(np.median(density)),
                    'density_seconds': round(time.time() - t1, 1)}
        for lname, lab in labels.items():
            sec_info[f'label_{lname}'] = {'n_pos': int((lab == 1).sum()),
                                          'n_neg': int((lab == 0).sum()),
                                          'n_excluded': int((lab == -1).sum())}
            if sec_info[f'label_{lname}']['n_pos'] == 0 or sec_info[f'label_{lname}']['n_neg'] == 0:
                raise SystemExit(f'{sec}/{lname}: degenerate label')
        manifest['sections'][sec] = sec_info
        print(f'  {sec:10s} n={n:>9,}  clean {sec_info["label_clean"]["n_pos"]:>7,}+/'
              f'{sec_info["label_clean"]["n_neg"]:>7,}-  full '
              f'{sec_info["label_full"]["n_pos"]:>7,}+/{sec_info["label_full"]["n_neg"]:>7,}-'
              f'  density_med={sec_info["density_median"]:.0f}  ({sec_info["density_seconds"]}s)')

    manifest['wall_min'] = round((time.time() - t0) / 60, 2)
    (args.out / 'labels' / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'wrote {out_dir}  ({manifest["wall_min"]} min)')


if __name__ == '__main__':
    main()
