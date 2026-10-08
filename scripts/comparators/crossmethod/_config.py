"""Shared paths / constants for the cross-method label-prediction benchmark.

The question (Wesley, 2026-09-10 meeting): using the pathologist's tissue
annotation as the label, fit the SIMPLEST possible classifier -- a
one-dimensional logistic regression -- on ONE feature at a time.  ALARMIST
contributes motif 24 and motif 10 loadings; every comparator contributes each
of its individual LR pairs.  That yields, per method, the whole DISTRIBUTION of
"how well can a single LR pair predict this label", and the question is where
ALARMIST's motif falls in it.

Nothing here re-runs any comparator.  Every score matrix is read as it was
written by its own run; `comparator-benchmark` SKILL.md:45-46 forbids
harmonising spatial scale across methods, so each method is compared only
against ITS OWN distribution and LR counts are never put side by side.
"""
from __future__ import annotations

from pathlib import Path

# ------------------------------------------------------------------ host paths
ADATA = Path('/data1/tanseyw/projects/fanj2/alarmist_luad/adata_with_region.h5ad')
RESULTS = Path('/data1/tanseyw/projects/fanj2/results')
OUT = RESULTS / '_crossmethod' / 'label_auc'

# ------------------------------------------------------------------ the design
# Canonical order (ALARMIST's own row ordering, per _common/luad_config.sh).
SECTIONS = ['P17_AIS', 'P17_LUAD', 'P21_AIS', 'P21_LUAD']

# obs['tissue_layer'] is a PATHOLOGIST annotation, not model-derived -- which is
# the whole point of using it as the label here.  Levels verified on disk.
TISSUE_LAYERS = ['Normal_distal', 'Normal_peri', 'Interface_normal',
                 'Interface_tumor', 'Tumor_region']

# Two binarisations.  `clean` is the main-figure contrast; `full` is the
# sensitivity analysis that keeps every cell.  User decision 2026-09-10.
LABEL_DEFS = {
    'clean': {'pos': ['Tumor_region'],
              'neg': ['Normal_peri', 'Normal_distal'],
              'drop': ['Interface_tumor', 'Interface_normal']},
    'full':  {'pos': ['Tumor_region', 'Interface_tumor'],
              'neg': ['Normal_peri', 'Normal_distal', 'Interface_normal'],
              'drop': []},
}

# ALARMIST features.  Motif 24 = healthy vasculature, motif 10 = tumor
# vasculature (CLAUDE.md, Fig 4).  Wesley asked for these two only.
ALARMIST_FEATURES = {'ALARMIST_motif_24': 'motif_24_loading',
                     'ALARMIST_motif_10': 'motif_10_loading'}

# All 25 motifs of the K=25 LUAD fit, scored as a CONTROL only.  Wesley asked for
# motifs 24 and 10 and those two carry the claim; this set exists so the obvious
# objection -- "you picked the two motifs that happen to work" -- can be answered
# with the whole distribution rather than an assurance.  Rows are tagged
# `kind == 'alarmist_all_motifs'` and are not part of the headline comparison.
N_MOTIFS = 25
ALL_MOTIF_FEATURES = {f'motif_{k}': f'motif_{k}_loading' for k in range(N_MOTIFS)}

# ------------------------------------------------------------ comparator inputs
# Each entry: how to reach that method's per-unit x per-LR score matrix.
# `unit` is the row unit the method natively scores -- NOT harmonised.
METHOD_INPUTS = {
    'spatialdm': {
        'unit': 'cell',
        'path': lambda s: RESULTS / f'spatialdm/LUAD/cellchatdb2/{s}/data/local_z.npz',
        'scale': 'SpatialDM default cutoff (l selected by the method)',
    },
    'liana': {
        'unit': 'cell',
        'path': lambda s: RESULTS / f'liana/LUAD/{s}/cellchatdb2_inflow/data/inflow_scores.npz',
        'scale': 'LIANA bandwidth 13.1454 um (deviation CD-1)',
    },
    'cytosignal': {
        'unit': 'cell',
        'path': lambda s: OUT / 'cytosignal_export' / s,
        'scale': "CytoSignal's own contact / diffusion kernels (defaults)",
    },
    'stlearn': {
        'unit': 'spot',
        'path': lambda s: RESULTS / f'stlearn/LUAD/cellchatdb2/{s}/data',
        'scale': 'stLearn grid ~51.4 x 51.2 um, distance=250 um',
    },
}

SEED = 0
