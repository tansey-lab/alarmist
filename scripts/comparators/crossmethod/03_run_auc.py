"""Score every feature of every method with a 1-D classifier against the pathologist label.

For each section x label x method, every individual feature (one LR pair, one
factor, one motif) is scored on its own by AUC -- which, for a one-dimensional
logistic regression, IS that model's AUC (proved numerically, see auc.py).
Output: one tidy row per (section, label, method, feature).

WHAT EACH METHOD CONTRIBUTES, AND ON WHAT UNIT
----------------------------------------------
  alarmist    cell   motif 24 + motif 10 per-cell loadings              [the claim]
  spatialdm   cell   local Moran z per LR pair                          [comparator]
  liana       cell   inflow score per `celltype^ligand^receptor`        [comparator]
  cytosignal  cell   LR score per CCI (diffusion + contact)             [comparator]
  stlearn     spot   LR co-expression score per pair, on ITS OWN grid   [comparator]
  cellchat    ct*    communication probability of the cell's own type   [comparator]
  mofaflex    cell   MOFA-Flex factor scores                    [aggregate baseline]
  reference   cell   local density / total counts / cell-type indicator [control]

*CellChat has NO per-cell score.  Its unit is (LR pair, sender type, receiver
type) within a GROUP (AIS or LUAD), so within a section it carries no spatial
variation at all.  The most generous per-cell feature derivable from it is the
communication probability attached to a cell's OWN cell type, built here both
as incoming (cell is receiver) and outgoing (cell is sender) so CellChat gets
twice the features and direction-agnostic credit on each.  Such a feature takes
at most 19 distinct values, so its AUC is bounded by what cell-type identity
alone achieves -- which the `best_cell_type_indicator` reference measures.
This is a fact about CellChat's unit of analysis, not a defect in the run.

stLearn scores a ~51 um GRID, not cells.  Its spots are labelled by majority
vote of the pathologist labels of the cells falling in them, and ALARMIST's
motif loadings are averaged onto the SAME spots so the two are read on one
unit.  That aggregation is applied to both sides or neither.

SCALE IS DELIBERATELY NOT HARMONISED (comparator-benchmark SKILL.md:45-46):
each method is compared only against ITS OWN distribution, and feature counts
are never compared across methods.

Usage:
    python 03_run_auc.py [--sections ...] [--methods ...] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import ALARMIST_FEATURES, ALL_MOTIF_FEATURES, LABEL_DEFS, \
    METHOD_INPUTS, OUT, RESULTS, SECTIONS
from auc import auc_abs, auc_dense, auc_sparse_nonneg, compare_implementations, \
    validate_logistic_identity

FEAT_CHUNK = 256
ROW_CHUNK = 100_000


# ------------------------------------------------------------------ helpers
def auc_featmajor(Z, y, chunk=FEAT_CHUNK):
    """AUC per ROW of a (n_features, n_units) array, without transposing it whole."""
    out = np.empty(Z.shape[0], dtype=np.float64)
    for a in range(0, Z.shape[0], chunk):
        b = min(a + chunk, Z.shape[0])
        out[a:b] = auc_dense(np.ascontiguousarray(Z[a:b].T), y)
    return out


def align(cells, want, method, section, *, subset_ok=False):
    """Row index into `cells` for each entry of `want`; -1 when absent."""
    order = {c: i for i, c in enumerate(cells)}
    idx = np.fromiter((order.get(k, -1) for k in want), dtype=np.int64, count=len(want))
    n_missing = int((idx < 0).sum())
    if n_missing and not subset_ok:
        raise SystemExit(f'{method}/{section}: {n_missing:,} labelled cells absent from its output')
    if n_missing:
        print(f'      {method}: {n_missing:,} of {len(want):,} labelled cells absent '
              f'(the method\'s own QC subset) -- dropped')
    return idx


def result(names, auc, *, unit, kind, n_units, n_pos, n_neg, note=''):
    return {'names': np.asarray(names, dtype=str), 'auc': np.asarray(auc, dtype=np.float64),
            'unit': unit, 'kind': kind, 'n_units': int(n_units),
            'n_pos': int(n_pos), 'n_neg': int(n_neg), 'note': note}


def alarmist_on(L, lab, y, sub=None, *, note=''):
    """ALARMIST motif AUC on EXACTLY the cells a comparator retained.

    Methods drop cells for their own reasons -- CytoSignal's removeLowQuality
    takes out 13% of them, SpatialDM leaves a handful with non-finite local z --
    so scoring the motif on all labelled cells and the comparator on a subset
    would compare two numbers with different denominators.  Every comparator
    therefore gets its own paired motif row computed on its own units.  AUC does
    not depend on row order, so a boolean mask over the labelled cells is
    enough.
    """
    names = np.array(list(ALARMIST_FEATURES), dtype=str)
    V = np.column_stack([L[n][lab] for n in names])
    if sub is not None:
        V, y = V[sub], y[sub]
    return result(names, auc_dense(V, y), unit='cell', kind='alarmist_paired',
                  n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum()), note=note)


# ------------------------------------------------------------------ scorers
# Every scorer receives the LABELLED cells only: `key` (join keys, in table
# order), `y` (bool, positive class) and `L` (the whole label table + mask).
def score_alarmist(section, key, y, L, lab):
    names = np.array(list(ALARMIST_FEATURES), dtype=str)
    V = np.column_stack([L[n][lab] for n in names])
    return result(names, auc_dense(V, y), unit='cell', kind='alarmist',
                  n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum())), V


def score_all_motifs(section, key, y, L, lab):
    """All 25 motifs of the K=25 fit -- the control for motif SELECTION.

    Motifs 24 and 10 were named from their cell-type composition, and Wesley
    named them for this test, so a reader is entitled to ask whether they are
    simply the best two of 25.  Scoring all of them answers that directly.
    """
    names = np.array(list(ALL_MOTIF_FEATURES), dtype=str)
    V = np.column_stack([L[n][lab] for n in names])
    return result(names, auc_dense(V, y), unit='cell', kind='alarmist_all_motifs',
                  n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum()),
                  note='control: motif selection, not part of the headline comparison')


def score_reference(section, key, y, L, lab):
    feats = {'local_density_80um': L['local_density_80um'][lab],
             'total_counts': L['total_counts'][lab]}
    ct = L['cell_type'][lab].astype(str)
    for t in np.unique(ct):
        feats[f'celltype_is_{t}'] = (ct == t).astype(np.float64)
    names = list(feats)
    V = np.column_stack([feats[n] for n in names])
    return result(names, auc_dense(V, y), unit='cell', kind='reference',
                  n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum()))


def score_spatialdm(section, key, y, L, lab):
    """SpatialDM local Moran z, (pairs x cells).

    `local_z.npz` was written from a bare ndarray, so its `pairs`/`cells` fields
    degenerated to positional indices '0','1',... (run_spatialdm.py:258-260 --
    `getattr(v, "index", np.arange(...))`).  The real names must come from the
    sibling files, and which pairs are present is NOT all of them: SpatialDM
    computes local statistics only for the pairs its global stage SELECTED
    (579 of 1,692 for P17_AIS).  `local_n_spots.csv` is indexed by exactly that
    selected set, in `global_res` order -- verified identical to
    `global_res[selected]` and to `selected_spots.csv.gz` -- so it is the name
    source, and the lengths are asserted rather than assumed.
    """
    d = METHOD_INPUTS['spatialdm']['path'](section).parent
    z = np.load(METHOD_INPUTS['spatialdm']['path'](section), allow_pickle=True)
    Z_all = z['values']

    pairs = pd.read_csv(d / 'local_n_spots.csv', index_col=0).index.to_numpy().astype(str)
    gr = pd.read_csv(d / 'global_res.csv', index_col=0)
    sel = gr.index[gr['selected'].astype(bool)].to_numpy().astype(str)
    if not np.array_equal(pairs, sel):
        raise SystemExit(f'spatialdm/{section}: local_n_spots index != global_res[selected]')
    cells = pd.read_csv(d / 'cell_meta.csv', usecols=['cell'])['cell'].to_numpy().astype(str)
    if Z_all.shape != (len(pairs), len(cells)):
        raise SystemExit(f'spatialdm/{section}: local_z {Z_all.shape} != '
                         f'({len(pairs)}, {len(cells)}) from the name files')

    idx = align(cells, key, 'spatialdm', section)
    Z = Z_all[:, idx]

    # Non-finite local z is a property of a HANDFUL OF CELLS, not of pairs: every
    # affected pair is NaN on exactly the same cells (10 of 182,378 in P17_AIS,
    # 3 in P21_AIS -- isolated cells with no usable spatial neighbourhood), and no
    # pair is all-NaN.  Dropping those cells keeps all 579 pairs; dropping pairs
    # instead would have thrown away 440 of them over 10 cells.
    good = np.isfinite(Z).all(axis=0)
    n_bad = int((~good).sum())
    if n_bad:
        print(f'      spatialdm: dropping {n_bad} cells with non-finite local z '
              f'(shared across all affected pairs); all {len(pairs)} pairs kept')
        Z, yy = Z[:, good], y[good]
    else:
        yy = y
    if not np.isfinite(Z).all():
        raise SystemExit(f'spatialdm/{section}: non-finite values remain after the cell drop')
    return (result(pairs, auc_featmajor(Z, yy), unit='cell', kind='comparator',
                   n_units=int(good.sum()), n_pos=int(yy.sum()), n_neg=int((~yy).sum()),
                   note=f'globally selected pairs only ({len(pairs)} of {len(gr)})'
                        + (f'; {n_bad} cells dropped for non-finite local z' if n_bad else '')),
            alarmist_on(L, lab, y, good if n_bad else None,
                        note='paired to SpatialDM cells'))


def score_liana(section, key, y, L, lab):
    z = np.load(METHOD_INPUTS['liana']['path'](section), allow_pickle=True)
    idx = align(z['cells'].astype(str), key, 'liana', section)
    V = z['values']                              # dense (n_cells, n_features)
    blocks = [sparse.csr_matrix(V[idx[a:a + ROW_CHUNK]])
              for a in range(0, idx.size, ROW_CHUNK)]
    del V, z
    M = sparse.vstack(blocks).tocsc()
    del blocks
    if M.data.size and M.data.min() < 0:
        raise SystemExit('liana inflow has negative values; the sparse path assumes nonneg')
    feats = np.load(METHOD_INPUTS['liana']['path'](section),
                    allow_pickle=True)['features'].astype(str)
    # CLAUDE.md: when a package primitive and a hand path both exist, prefer the
    # package one and RECORD the numerical comparison.  LIANA is the only matrix
    # here that is sparse and nonnegative, so it is the only place both apply.
    score_liana.impl_check = compare_implementations(M, y, n_features=200)
    return (result(feats, auc_sparse_nonneg(M, y), unit='cell', kind='comparator',
                   n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum())),
            alarmist_on(L, lab, y, note='paired to LIANA cells'))


def score_cytosignal(section, key, y, L, lab):
    d = METHOD_INPUTS['cytosignal']['path'](section)
    cells = np.array((d / 'cells.tsv').read_text().split('\n')[:-1], dtype=str)
    idx = align(cells, key, 'cytosignal', section, subset_ok=True)
    present = idx >= 0
    rows, yy = idx[present], y[present]
    names_all, auc_all = [], []
    for slot in ('diffusion_Raw_smooth', 'contact_Raw_smooth'):
        shp = json.loads((d / f'shape_{slot}.json').read_text())
        intr = np.array((d / f'intr_{slot}.tsv').read_text().split('\n')[:-1], dtype=str)
        M = np.memmap(d / f'score_{slot}.f32', dtype=np.float32, mode='r',
                      shape=(shp['n_cells'], shp['n_intr']), order='F')
        a = np.empty(shp['n_intr'], dtype=np.float64)
        for i in range(0, shp['n_intr'], FEAT_CHUNK):
            j = min(i + FEAT_CHUNK, shp['n_intr'])
            a[i:j] = auc_dense(np.asarray(M[:, i:j])[rows].astype(np.float64), yy)
        names_all.append(np.char.add(f'{slot.split("_")[0]}:', intr))
        auc_all.append(a)
        del M
    return (result(np.concatenate(names_all), np.concatenate(auc_all), unit='cell',
                   kind='comparator', n_units=int(present.sum()), n_pos=int(yy.sum()),
                   n_neg=int((~yy).sum()),
                   note='CytoSignal removeLowQuality subset' if not present.all() else ''),
            alarmist_on(L, lab, y, present,
                        note='paired to CytoSignal post-QC cells'))


def score_mofaflex(section, key, y, L, lab):
    fs = pd.read_csv(RESULTS / 'liana/LUAD/mofaflex_inflow_joint/data/factor_scores.csv.gz')
    fs = fs.rename(columns={fs.columns[0]: 'cell'}).set_index('cell')
    idx = align(fs.index.to_numpy().astype(str), key, 'mofaflex', section)
    V = fs.to_numpy(dtype=np.float64)[idx]
    return (result(np.array(fs.columns, dtype=str), auc_dense(V, y), unit='cell',
                   kind='aggregate', n_units=len(y), n_pos=int(y.sum()), n_neg=int((~y).sum())),
            alarmist_on(L, lab, y, note='paired to MOFA-Flex cells'))


def score_cellchat(section, key, y, L, lab):
    group = 'AIS' if section.endswith('AIS') else 'LUAD'
    net = pd.read_csv(RESULTS / f'cellchat/LUAD/cellchatdb2/quant/{group}_net_full.csv.gz')
    ct = L['cell_type'][lab].astype(str)
    names, aucs, valid_any = [], [], None
    for role, col in (('incoming', 'target'), ('outgoing', 'source')):
        agg = net.groupby(['interaction_name', col])['prob'].sum().unstack(fill_value=0.0)
        codes = pd.Categorical(ct, categories=list(agg.columns))
        valid = np.asarray(~pd.isna(codes))
        valid_any = valid if valid_any is None else (valid_any | valid)
        V = agg.to_numpy(dtype=np.float64).T[codes.codes[valid]]
        names.append(np.char.add(f'{role}:', np.array(agg.index, dtype=str)))
        aucs.append(auc_dense(V, y[valid]))
    sub = None if valid_any.all() else valid_any
    return (result(np.concatenate(names), np.concatenate(aucs), unit='cell_type_constant',
                   kind='comparator', n_units=int(valid_any.sum()), n_pos=int(y[valid_any].sum()),
                   n_neg=int((~y[valid_any]).sum()),
                   note='piecewise-constant on cell type; bounded by cell-type identity'),
            alarmist_on(L, lab, y, sub, note='paired to CellChat-scorable cells'))


def score_stlearn(section, key, y, L, lab):
    """stLearn's own ~51 um grid: spots labelled by majority vote of their cells."""
    from scipy.spatial import cKDTree
    d = METHOD_INPUTS['stlearn']['path'](section)
    meta = pd.read_csv(d / 'spot_meta.csv')
    scores = pd.read_csv(d / 'spot_lr_scores.csv.gz', index_col=0)
    centres = meta[['imagecol', 'imagerow']].to_numpy(dtype=np.float64)
    spot_of = {s: i for i, s in enumerate(meta['spot'])}
    row_of_score = np.array([spot_of[s] for s in scores.index], dtype=np.int64)

    xy = np.column_stack([L['x'][lab], L['y'][lab]])
    dist, nearest = cKDTree(centres).query(xy, k=1)
    # A cell further than one spot diagonal from every centre is off-grid; drop it.
    diag = float(np.hypot(*(meta[['imagecol', 'imagerow']].diff().abs().max())))
    on = dist <= max(diag, 60.0)
    pos = np.bincount(nearest[on & y], minlength=len(meta)).astype(np.float64)
    neg = np.bincount(nearest[on & ~y], minlength=len(meta)).astype(np.float64)
    tot = pos + neg
    usable = (tot > 0)[row_of_score]
    yp = (pos > neg)[row_of_score][usable]
    ties = int(((pos == neg) & (tot > 0)).sum())
    if ties:
        print(f'      stlearn: {ties} spots tied on the majority vote -- counted negative')
    # How much evidence is behind each spot's label?  A spot carried by a single
    # labelled cell is a noisy label.  That noise hits the comparator row and the
    # paired ALARMIST row equally (they share these spots), so it costs power
    # rather than fairness -- but the reader should be able to see it.
    lc = tot[row_of_score][usable]
    print(f'      stlearn: labelled cells per usable spot -- median {np.median(lc):.0f}, '
          f'{int((lc == 1).sum()):,} ({100 * (lc == 1).mean():.1f}%) rest on one cell, '
          f'{int((lc >= 5).sum()):,} ({100 * (lc >= 5).mean():.1f}%) on >=5')

    # ALARMIST onto the SAME spots, so both sides get the same aggregation.
    extra = {}
    for name in ALARMIST_FEATURES:
        v = L[name][lab]
        s = np.bincount(nearest[on], weights=v[on], minlength=len(meta))
        n = np.bincount(nearest[on], minlength=len(meta)).astype(np.float64)
        extra[name] = np.divide(s, n, out=np.zeros_like(s), where=n > 0)[row_of_score][usable]

    V = scores.to_numpy(dtype=np.float64)[usable]
    res = result(np.array(scores.columns, dtype=str), auc_dense(V, yp), unit='spot',
                 kind='comparator', n_units=int(usable.sum()), n_pos=int(yp.sum()),
                 n_neg=int((~yp).sum()),
                 note=f'spot label = majority vote of its cells; median '
                      f'{np.median(lc):.0f} labelled cells/spot')
    alar = result(np.array(list(extra), dtype=str),
                  auc_dense(np.column_stack([extra[k] for k in extra]), yp),
                  unit='spot', kind='alarmist_paired', n_units=int(usable.sum()),
                  n_pos=int(yp.sum()), n_neg=int((~yp).sum()),
                  note='motif loadings averaged onto stLearn\'s grid')
    return res, alar


SCORERS = {'alarmist': score_alarmist, 'all_motifs': score_all_motifs,
           'reference': score_reference,
           'spatialdm': score_spatialdm, 'liana': score_liana,
           'cytosignal': score_cytosignal, 'mofaflex': score_mofaflex,
           'cellchat': score_cellchat, 'stlearn': score_stlearn}


def to_rows(section, label, method, r):
    a = r['auc']
    return pd.DataFrame({
        'section': section, 'label': label, 'method': method, 'unit': r['unit'],
        'kind': r['kind'], 'feature': r['names'], 'n_units': r['n_units'],
        'n_pos': r['n_pos'], 'n_neg': r['n_neg'], 'auc': a, 'auc_abs': auc_abs(a),
        'direction': np.where(a >= 0.5, 'up_in_tumor', 'down_in_tumor'), 'note': r['note'],
    })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sections', nargs='*', default=SECTIONS)
    ap.add_argument('--methods', nargs='*', default=list(SCORERS))
    ap.add_argument('--out', type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    all_rows, checks = [], {}
    for section in args.sections:
        d = np.load(args.out / 'labels' / f'{section}.npz', allow_pickle=False)
        L = {k: d[k] for k in d.files}
        key_all = L['key'].astype(str)
        for label in LABEL_DEFS:
            raw = L[f'label_{label}']
            lab = raw >= 0
            y = raw[lab] == 1
            key = key_all[lab]
            print(f'{section} / {label}: {lab.sum():,} labelled cells  '
                  f'{y.sum():,}+ / {(~y).sum():,}-')
            for method in args.methods:
                t0 = time.time()
                try:
                    out = SCORERS[method](section, key, y, L, lab)
                except FileNotFoundError as e:
                    print(f'      {method:11s} SKIPPED -- {e}')
                    continue
                extras = []
                if method == 'alarmist':
                    out, V = out
                    if 'logistic_identity' not in checks:
                        checks['logistic_identity'] = validate_logistic_identity(V, y, n_features=2)
                        print(f'      [check] 1-D logistic AUC vs rank AUC: max|diff| = '
                              f'{checks["logistic_identity"]["max_abs_diff"]:.2e}')
                elif method in ('all_motifs', 'reference'):
                    pass                       # standalone; no paired row to emit
                else:
                    out, alar = out
                    extras = [(f'alarmist@{method}', alar)]
                for nm, r in [(method, out)] + extras:
                    all_rows.append(to_rows(section, label, nm, r))
                print(f'      {method:11s} {len(out["names"]):>6,} features  n={out["n_units"]:>8,}'
                      f'  max|AUC-.5|={np.abs(out["auc"] - .5).max():.3f}  ({time.time()-t0:.0f}s)')

    if hasattr(score_liana, 'impl_check'):
        checks['sparse_vs_dense_auc'] = score_liana.impl_check
        print(f'[check] package sparse U vs dense rankdata AUC on 200 LIANA features: '
              f'max|diff| = {score_liana.impl_check["max_abs_diff"]:.2e}')
    df = pd.concat(all_rows, ignore_index=True)
    df.to_csv(args.out / 'feature_auc.csv.gz', index=False)
    (args.out / 'checks.json').write_text(json.dumps(checks, indent=2))
    print(f'\nwrote {args.out / "feature_auc.csv.gz"}  ({len(df):,} rows)')


if __name__ == '__main__':
    main()
