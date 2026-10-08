"""Figures for the 1-D classifier benchmark.

Figure 1 (per label): four panels, one per section.  Within a panel, one row per
method showing the FULL distribution of that method's single-feature AUCs, with
ALARMIST's two motifs marked in the row.  The motif marker in each row is the
one measured on THAT row's unit -- per-cell for the cell-level methods, averaged
onto stLearn's grid for the stLearn row -- so the comparison inside a row is
always like-for-like.  stLearn is drawn below a divider because its unit is a
~51 um spot, not a cell; `comparator-benchmark` SKILL.md:45-46 forbids
harmonising that away, so it is separated rather than merged.

Figure 2: how many of each method's features match or beat the motif, as a
fraction, per section.

Colour encodes section (Okabe-Ito, the fixed CLAUDE.md mapping) and is
redundant with panel position; motif identity is encoded by MARKER SHAPE, so
nothing here is colour-alone.

Usage:  python 05_plot.py [--out DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/comparators
from _common.plotting import apply_publication_style, save_all_formats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import OUT, SECTIONS

SECCOL = {'P17_AIS': '#0072b2', 'P17_LUAD': '#e69f00',
          'P21_AIS': '#009e73', 'P21_LUAD': '#cc79a7'}
MOTIF_MARK = {'ALARMIST_motif_24': ('o', 'motif 24 (healthy vasculature)'),
              'ALARMIST_motif_10': ('D', 'motif 10 (tumor vasculature)')}
CELL_METHODS = ['spatialdm', 'cytosignal', 'liana', 'cellchat', 'mofaflex']
SPOT_METHODS = ['stlearn']
PRETTY = {'spatialdm': 'SpatialDM', 'cytosignal': 'CytoSignal', 'liana': 'LIANA inflow',
          'cellchat': 'CellChat*', 'mofaflex': 'MOFA-Flex (factors)', 'stlearn': 'stLearn'}


def panel(ax, g, section, alarm, rng):
    methods = [m for m in CELL_METHODS + SPOT_METHODS if m in set(g.method)]
    ypos = {m: i for i, m in enumerate(methods)}

    for m in methods:
        a = g.loc[g.method == m, 'auc_abs'].to_numpy()
        y = ypos[m] + rng.uniform(-0.16, 0.16, a.size)
        ax.scatter(a, y, s=1.4, c='0.55', alpha=0.30, linewidths=0, rasterized=True,
                   zorder=2)
        q1, med, q3 = np.quantile(a, [0.25, 0.5, 0.75])
        ax.plot([q1, q3], [ypos[m]] * 2, c='0.15', lw=1.6, solid_capstyle='butt', zorder=3)
        ax.plot([med], [ypos[m]], marker='|', c='k', ms=9, mew=1.6, zorder=4)
        ax.text(0.998, ypos[m] - 0.26, f'n={a.size:,}', fontsize=4.6, color='0.4',
                va='center', ha='right')

    for m in methods:
        # Always the PAIRED motif row -- the one 03_run_auc.py scored on exactly
        # this method's own units (its post-QC cells, its grid spots).  Using the
        # all-cells motif row here would draw a marker measured on a different
        # denominator from the distribution behind it.
        src = f'alarmist@{m}'
        for motif, (mk, _) in MOTIF_MARK.items():
            v = alarm.get((src, motif))
            if v is None:
                continue
            ax.plot([v], [ypos[m]], marker=mk, ms=5.2, mfc=SECCOL[section],
                    mec='k', mew=0.7, ls='none', zorder=6)

    if SPOT_METHODS and set(SPOT_METHODS) & set(methods):
        ax.axhline(len(CELL_METHODS) - 0.5, c='0.75', lw=0.7, ls=(0, (4, 3)), zorder=1)

    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels([PRETTY[m] for m in methods], fontsize=6)
    ax.set_ylim(-0.6, len(methods) - 0.4)
    ax.invert_yaxis()
    ax.set_xlim(0.5, 1.0)
    ax.set_title(section, fontsize=7, color=SECCOL[section], pad=3, loc='left')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=6)


def figure_distributions(df, ref, label, out):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(0)
    d = df[df.label == label]
    secs = [s for s in SECTIONS if s in set(d.section)]
    fig, axes = plt.subplots(len(secs), 1, figsize=(6.9, 2.0 * len(secs)), sharex=True)
    axes = np.atleast_1d(axes)

    for ax, section in zip(axes, secs):
        g = d[(d.section == section) & (d.kind.isin(['comparator', 'aggregate']))]
        alarm = {(r.method, r.feature): r.auc_abs
                 for r in d[(d.section == section)
                            & (d.kind == 'alarmist_paired')].itertuples()}
        panel(ax, g, section, alarm, rng)

        r = ref[(ref.section == section) & (ref.label == label)]
        ct = r[r.feature.str.startswith('celltype_is_')]
        for feat, style, lab in (
                ('local_density_80um', (0, (1, 1.6)), 'local density (80 um)'),
                (None, (0, (5, 2)), 'best cell-type indicator')):
            v = (r.loc[r.feature == feat, 'auc_abs'].max() if feat
                 else (ct['auc_abs'].max() if len(ct) else np.nan))
            if np.isfinite(v):
                ax.axvline(v, c='0.45', lw=0.7, ls=style, zorder=0)

    axes[-1].set_xlabel('AUC of a one-dimensional logistic regression on that single feature\n'
                        '(direction-agnostic: max(AUC, 1 - AUC))', fontsize=6.5)

    handles = [Line2D([], [], marker=mk, ls='none', mfc='0.4', mec='k', mew=0.7, ms=5.2,
                      label=lab) for mk, lab in MOTIF_MARK.values()]
    handles += [Line2D([], [], c='0.45', lw=0.7, ls=(0, (1, 1.6)), label='local density (80 um)'),
                Line2D([], [], c='0.45', lw=0.7, ls=(0, (5, 2)), label='best cell-type indicator'),
                Line2D([], [], c='0.55', marker='o', ls='none', ms=2, label='one LR pair / factor'),
                Line2D([], [], c='0.15', lw=1.6, label='IQR, | = median')]
    fig.legend(handles=handles, fontsize=5.6, loc='upper center', bbox_to_anchor=(0.5, 0.962),
               frameon=False, ncol=3, handlelength=1.8, columnspacing=1.6)
    fig.suptitle(f'Can a single LR pair predict the pathologist label as well as an '
                 f'ALARMIST motif?   (label = {label})', fontsize=7.6, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    save_all_formats(fig, out / f'motif_vs_single_lr_{label}', dpi=450, close=True, verbose=True)


def figure_fraction(S, label, out):
    import matplotlib.pyplot as plt

    d = S[S.label == label]
    methods = [m for m in CELL_METHODS + SPOT_METHODS if m in set(d.method)]
    motifs = sorted(d.motif.unique())
    fig, axes = plt.subplots(1, len(motifs), figsize=(3.4 * len(motifs), 2.9), sharey=True)
    axes = np.atleast_1d(axes)

    for ax, motif in zip(axes, motifs):
        for i, m in enumerate(methods):
            for section in SECTIONS:
                r = d[(d.method == m) & (d.motif == motif) & (d.section == section)]
                if not len(r):
                    continue
                r = r.iloc[0]
                f = max(r.frac_ge, 1.0 / (r.n_features * 4))     # so 0 is still visible
                ax.plot([f], [i], marker='o', ms=4.4, mfc=SECCOL[section], mec='k',
                        mew=0.5, ls='none', alpha=0.95)
        ax.set_xscale('log')
        ax.axvline(0.05, c='0.7', lw=0.7, ls=(0, (4, 3)))
        ax.set_yticks(range(len(methods)))
        ax.set_yticklabels([PRETTY[m] for m in methods], fontsize=6)
        ax.invert_yaxis()
        ax.set_title(motif.replace('ALARMIST_motif_', 'motif '), fontsize=7)
        ax.set_xlabel('fraction of that method\'s features\nthat match or beat the motif',
                      fontsize=6.2)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=6)
    handles = [plt.Line2D([], [], marker='o', ls='none', mfc=SECCOL[s], mec='k', mew=0.5,
                          ms=4.4, label=s) for s in SECTIONS]
    axes[-1].legend(handles=handles, fontsize=5.6, frameon=False, loc='lower right')
    fig.tight_layout()
    save_all_formats(fig, out / f'fraction_beating_motif_{label}', dpi=450, close=True, verbose=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=OUT)
    args = ap.parse_args()
    fig_dir = args.out / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)

    apply_publication_style(**{'font.size': 6.5, 'axes.linewidth': 0.6,
                               'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
    df = pd.read_csv(args.out / 'feature_auc.csv.gz')
    S = pd.read_csv(args.out / 'summary_motif_vs_methods.csv')
    ref = pd.read_csv(args.out / 'reference_features.csv')

    for label in df.label.unique():
        figure_distributions(df, ref, label, fig_dir)
        figure_fraction(S, label, fig_dir)
    print(f'figures in {fig_dir}')


if __name__ == '__main__':
    main()
