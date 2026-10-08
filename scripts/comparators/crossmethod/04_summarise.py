"""Turn per-feature AUCs into the answer: where does ALARMIST's motif fall?

For each section x label x method we report the whole distribution of that
method's single-feature AUCs and the position of ALARMIST motif 24 and motif 10
in it.

HOW THE POSITION IS REPORTED, AND WHAT IT IS NOT
------------------------------------------------
`n_ge` counts the method's own features whose direction-agnostic AUC matches or
beats the motif's, out of `n_features`.  `empirical_p = (1 + n_ge) / (1 + n)` is
the rank of the motif in the OBSERVED distribution, in the style of a
permutation p-value.  It is NOT a test against a sampling null: the features are
not exchangeable draws, and the cells they are scored on are spatially
autocorrelated, so a nominal p from any parametric test would be meaningless
here.  The defensible statement is the rank itself -- "k of N single LR pairs
predict the pathologist's label at least as well as the motif does" -- and that
is what the write-up should say.

DIRECTION-AGNOSTIC ON PURPOSE.  Everything is ranked on `auc_abs =
max(auc, 1 - auc)`, so an LR pair that is strongly DOWN in tumour gets full
credit.  That makes the comparison harder for ALARMIST, not easier.

Usage:  python 04_summarise.py [--out DIR]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import OUT, RESULTS, SECTIONS

# Every comparator is compared against ITS OWN paired ALARMIST row, written by
# 03_run_auc.py as `alarmist@<method>` and computed on exactly the units that
# method retained -- CytoSignal's post-QC cells, SpatialDM's finite-z cells,
# stLearn's grid spots.  Comparing against the all-cells motif row instead would
# put a different denominator on each side of the comparison.
COMPARE_KINDS = ('comparator', 'aggregate')


def cytosignal_lut():
    """CCI-xxxxx -> 'LIGAND - RECEPTOR', from the runs' own signif_summary files.

    CytoSignal names interactions by opaque id in the score matrix; the readable
    name lives only in `quant/signif_summary_<slot>.csv`.  Without this the top
    features are unreadable and the collagen observation that motivated the
    analysis cannot be checked.
    """
    lut = {}
    for sec in SECTIONS:
        for slot in ('diffusion', 'contact', 'Raw'):
            f = (RESULTS / f'cytosignal/LUAD/cellchatdb2/{sec}/quant/'
                 f'signif_summary_{slot}_Raw_smooth.csv')
            if f.exists():
                t = pd.read_csv(f)
                lut.update(dict(zip(t['interaction_id'], t['name'])))
    return lut


def readable(feature, lut):
    if isinstance(feature, str) and ':' in feature and feature.split(':')[1] in lut:
        pre, cid = feature.split(':', 1)
        return f'{pre}:{lut[cid]}'
    return feature


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=OUT)
    args = ap.parse_args()

    df = pd.read_csv(args.out / 'feature_auc.csv.gz')
    LUT = cytosignal_lut()
    df['feature_label'] = [readable(f, LUT) for f in df['feature']]
    rows, tops = [], []

    for (section, label), g in df.groupby(['section', 'label'], sort=False):
        alarm = {k: v.set_index('feature')['auc_abs'].to_dict()
                 for k, v in g[g.kind == 'alarmist_paired'].groupby('method')}
        paired_n = {k: int(v['n_units'].iloc[0])
                    for k, v in g[g.kind == 'alarmist_paired'].groupby('method')}
        for method, m in g[g.kind.isin(COMPARE_KINDS)].groupby('method', sort=False):
            src = f'alarmist@{method}'
            if src not in alarm:
                print(f'  [warn] no paired ALARMIST row for {method} in {section}/{label}')
                continue
            if paired_n[src] != int(m['n_units'].iloc[0]):
                raise SystemExit(
                    f'{section}/{label}/{method}: paired ALARMIST was scored on '
                    f'{paired_n[src]:,} units but the method on {int(m["n_units"].iloc[0]):,}')
            a = m['auc_abs'].to_numpy()
            n = a.size
            for motif, mv in alarm[src].items():
                n_ge = int((a >= mv).sum())
                rows.append({
                    'section': section, 'label': label, 'method': method,
                    'unit': m['unit'].iloc[0], 'kind': m['kind'].iloc[0],
                    'n_units': int(m['n_units'].iloc[0]), 'n_pos': int(m['n_pos'].iloc[0]),
                    'n_neg': int(m['n_neg'].iloc[0]), 'n_features': n,
                    'motif': motif, 'motif_auc_abs': mv,
                    'method_median_auc_abs': float(np.median(a)),
                    'method_p90_auc_abs': float(np.quantile(a, 0.90)),
                    'method_max_auc_abs': float(a.max()),
                    'method_best_feature': m.loc[m['auc_abs'].idxmax(), 'feature_label'],
                    'n_ge': n_ge, 'frac_ge': n_ge / n,
                    'empirical_p': (1 + n_ge) / (1 + n),
                })
            t = m.nlargest(5, 'auc_abs')[['feature', 'feature_label', 'auc', 'auc_abs']].copy()
            t.insert(0, 'method', method); t.insert(0, 'label', label)
            t.insert(0, 'section', section); t['rank'] = np.arange(1, len(t) + 1)
            tops.append(t)

    S = pd.DataFrame(rows)
    S.to_csv(args.out / 'summary_motif_vs_methods.csv', index=False)
    pd.concat(tops, ignore_index=True).to_csv(args.out / 'top_features_per_method.csv', index=False)

    ref = df[df.kind == 'reference'].copy()
    ref.to_csv(args.out / 'reference_features.csv', index=False)

    # ---------------------------------------------------------------- printout
    for label in S['label'].unique():
        print(f'\n{"="*104}\nLABEL = {label}\n{"="*104}')
        sub = S[S.label == label]
        for motif in sorted(sub['motif'].unique()):
            print(f'\n--- {motif} ---')
            print(f'{"section":10s} {"method":11s} {"unit":10s} {"n_units":>9s} '
                  f'{"n_feat":>7s} {"motif":>7s} {"med":>6s} {"max":>6s} {"n>=":>6s} '
                  f'{"emp_p":>9s}  best feature')
            for _, r in sub[sub.motif == motif].iterrows():
                print(f'{r.section:10s} {r.method:11s} {r.unit:10s} {r.n_units:>9,} '
                      f'{r.n_features:>7,} {r.motif_auc_abs:>7.3f} '
                      f'{r.method_median_auc_abs:>6.3f} {r.method_max_auc_abs:>6.3f} '
                      f'{r.n_ge:>6,} {r.empirical_p:>9.2e}  {str(r.method_best_feature)[:34]}')

    print(f'\n--- reference features (context, not part of the comparison) ---')
    for (section, label), g in ref.groupby(['section', 'label'], sort=False):
        best_ct = g[g.feature.str.startswith('celltype_is_')].nlargest(1, 'auc_abs')
        dens = g[g.feature == 'local_density_80um']
        line = f'{section:10s} {label:6s} '
        if len(dens):
            line += f'local_density_80um={dens.auc_abs.iloc[0]:.3f}  '
        if len(best_ct):
            line += f'best_cell_type={best_ct.feature.iloc[0][12:]} ({best_ct.auc_abs.iloc[0]:.3f})'
        print(line)

    print(f'\nwrote {args.out / "summary_motif_vs_methods.csv"}')


if __name__ == '__main__':
    main()
