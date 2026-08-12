#!/usr/bin/env python
"""Cosine similarity between LIANA inflow+MOFA-Flex factors and ALARMIST GBM motifs.

This follows the matching procedure in .claude/skills/alarmist (cosine on the LRI factors,
restricted to SHARED cell-type-LRI combinations). Three things must be handled explicitly or
the heatmap is an artefact:

1. THE FEATURE SPACES ARE NOT THE SAME OBJECT.
   ALARMIST  : (sender, receiver, ligand, receptor, contact mode)   25,271 features
   MOFA-Flex : (sender, ligand, receptor)                              779 features
   The only common space is (sender, ligand, receptor), so ALARMIST must be COLLAPSED over
   receiver and contact mode. That discards the coordinates ALARMIST has and LIANA does not,
   so this comparison is deliberately biased in LIANA's favour -- it is an upper bound on the
   agreement, not a neutral measurement. The fraction of ALARMIST mass lost to the collapse
   is reported.

2. THE TWO RUNS USED DIFFERENT CellChatDB EXPORTS.
   ALARMIST used the pre-2026-07-28 export, LIANA the re-export; 1,120 of 3,218 LR keys differ
   by heteromeric SUBUNIT ORDER alone (TGFBR2_TGFBR1 vs TGFBR1_TGFBR2). Keys are therefore
   canonicalised by sorting subunits before joining, and the raw-vs-canonical overlap is
   reported so the reader can see how much the canonicalisation recovered.

3. MOFA-FLEX IS SIGNED, BPTF IS NOT.
   MOFA-Flex weights are signed and each factor's overall sign is arbitrary; BPTF factors are
   non-negative. Cosine between a signed and a non-negative vector is not meaningful as a
   signed quantity, so |weight| is used on the MOFA-Flex side and |cosine| is reported.

Both ALARMIST value columns are compared, because they behave very differently:
   `factor` = raw V   -- dominated by high-prevalence adhesion pairs
   `score`  = V*      -- V / (mean_LR + 1), the prevalence-normalised column the skill says to
                         rank on. LIANA has no equivalent normalisation.

Run:
    /Users/jiayifan/anaconda3/envs/bptf/bin/python \\
        scripts/comparators/liana/cosine_factors_vs_motifs.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")

p = argparse.ArgumentParser()
p.add_argument("--mofaflex-loadings",
               default=str(ROOT / "results/comparators/liana/GBM/mofaflex_inflow_reachnorm/data/mofaflex_loadings.csv"))
p.add_argument("--alarmist-motifs", default=str(ROOT / "results/GBM/bptf/lri_motifs.csv"))
p.add_argument("--out-dir", default=str(ROOT / "results/comparators/liana/GBM/vs_alarmist"))
p.add_argument("--label", default="mofaflex_reachnorm")
p.add_argument("--n-perms", type=int, default=200,
               help="permutations for the null on max-cosine per motif")
p.add_argument("--sign-mode", choices=["poles", "abs", "signed"], default="poles",
               help="how to reconcile MOFA-Flex signed weights with non-negative BPTF "
                    "factors; see the comment at the loading step. Default 'poles'.")
a = p.parse_args()
OUT = Path(a.out_dir); (OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(exist_ok=True)

canon = lambda s: "_".join(sorted(str(s).split("_")))

# ------------------------------------------------------------------ MOFA-Flex side
M = pd.read_csv(a.mofaflex_loadings)
fac = [c for c in M.columns if c.startswith("Factor")]
M["key"] = (M["source"].astype(str) + "|" + M["ligand_complex"].map(canon)
            + "|" + M["receptor_complex"].map(canon))

# --- how to reconcile a SIGNED factorisation with a NON-NEGATIVE one -------------------
# A BPTF motif is non-negative: it means "these interactions are high TOGETHER".
# A MOFA-Flex factor is an AXIS: features at the + pole and the - pole are ANTI-correlated,
# and the global sign of a factor is arbitrary.
#
#   abs    -- take |weight|. Treats the factor as one unsigned feature set. This MERGES the
#             two anti-correlated poles, so it manufactures similarity: a factor contrasting
#             A against B is scored as if it contained A and B together. On this run 57.2% of
#             weights are negative and the minor pole holds a median 38.5% of the mass, so
#             this is a large distortion, not a rounding. Kept only for reproducibility.
#   poles  -- split each factor into  Factor_k(+) = max(w,0)  and  Factor_k(-) = max(-w,0).
#             Each pole is then a genuine non-negative set, directly comparable to a motif,
#             and the contrast between poles is preserved. THIS IS THE DEFAULT.
#   signed -- keep signed weights and report |cosine|; handles the arbitrary global sign but
#             a signed-vs-non-negative cosine is still not a like-for-like quantity.
if a.sign_mode == "abs":
    Mg = M.groupby("key")[fac].apply(lambda d: d.abs().sum())
    col_labels = fac
elif a.sign_mode == "signed":
    Mg = M.groupby("key")[fac].sum()
    col_labels = fac
else:  # poles
    parts = {}
    for f in fac:
        parts[f"{f} (+)"] = M[f].clip(lower=0)
        parts[f"{f} (-)"] = (-M[f]).clip(lower=0)
    Mp = pd.concat([M[["key"]], pd.DataFrame(parts, index=M.index)], axis=1)
    Mg = Mp.groupby("key").sum()
    col_labels = list(parts)
    # a pole with almost no mass is not a program; drop it rather than score noise
    keep_cols = [c for c in col_labels if Mg[c].sum() > 0]
    dropped = [c for c in col_labels if c not in keep_cols]
    if dropped:
        print(f"  dropped {len(dropped)} empty pole(s): {dropped}")
    Mg, col_labels = Mg[keep_cols], keep_cols
print(f"MOFA-Flex : {M.shape[0]} features -> {Mg.shape[0]} canonical (sender,L,R) keys, "
      f"{len(col_labels)} columns  [sign_mode={a.sign_mode}]")
fac = col_labels

# ------------------------------------------------------------------ ALARMIST side
A = pd.read_csv(a.alarmist_motifs)
A["key"] = (A["celltype1"].astype(str) + "|" + A["ligand"].map(canon)
            + "|" + A["receptor"].map(canon))
n_motifs = A.motif_idx.nunique()
raw_overlap = len(set(M["source"].astype(str) + "|" + M.ligand_complex.astype(str) + "|"
                      + M.receptor_complex.astype(str))
                  & set(A.celltype1.astype(str) + "|" + A.ligand.astype(str) + "|"
                        + A.receptor.astype(str)))

results = {}
for valcol, vlabel in [("factor", "raw V"), ("score", "V* = V/(mean_LR+1)")]:
    Ag = A.pivot_table(index="key", columns="motif_idx", values=valcol, aggfunc="sum")
    shared = Mg.index.intersection(Ag.index)
    # how much ALARMIST mass survives the collapse onto (sender, L, R) and the join?
    kept = A[A.key.isin(shared)][valcol].sum() / A[valcol].sum()
    Mx = Mg.loc[shared].values                      # keys x factors
    Ax = Ag.loc[shared].fillna(0).values            # keys x motifs
    Mn = Mx / np.maximum(np.linalg.norm(Mx, axis=0, keepdims=True), 1e-12)
    An = Ax / np.maximum(np.linalg.norm(Ax, axis=0, keepdims=True), 1e-12)
    C = pd.DataFrame(Mn.T @ An, index=fac, columns=[f"motif{m}" for m in Ag.columns])
    C.to_csv(OUT / "data" / f"cosine_{a.label}_{a.sign_mode}_vs_alarmist_{valcol}.csv")
    # --- permutation null -------------------------------------------------------------
    # `poles` doubles the number of candidate vectors (40 vs 20), so a higher max-cosine per
    # motif could be pure multiple comparisons. Permuting each MOFA column over the shared
    # keys destroys the correspondence while preserving that column's sparsity and magnitude
    # distribution, which is the right null for exactly that concern.
    rng = np.random.default_rng(0)
    null = np.empty((a.n_perms, C.shape[1]))
    for i in range(a.n_perms):
        Mp_ = Mx.copy()
        for j in range(Mp_.shape[1]):
            Mp_[:, j] = rng.permutation(Mp_[:, j])
        Mpn = Mp_ / np.maximum(np.linalg.norm(Mp_, axis=0, keepdims=True), 1e-12)
        null[i] = (Mpn.T @ An).max(0)
    obs_max = C.max(0).values
    emp_p = (null >= obs_max[None, :]).mean(0)
    print(f"  null ({a.n_perms} perms): median max-cosine {np.median(null):.3f}, "
          f"95th pct {np.percentile(null, 95):.3f}  |  observed median {np.median(obs_max):.3f}")
    print(f"  motifs beating the null at p<0.05 : {int((emp_p < 0.05).sum())}/{C.shape[1]}")
    pd.DataFrame({"motif": C.columns, "best_pole": C.idxmax(0).values,
                  "cosine": obs_max, "null_median": np.median(null, 0),
                  "null_p95": np.percentile(null, 95, axis=0), "emp_p": emp_p}).to_csv(
        OUT / "data" / f"bestmatch_{a.label}_{a.sign_mode}_{valcol}.csv", index=False)

    results[valcol] = dict(
        null_median_max_cosine=float(np.median(null)),
        null_p95=float(np.percentile(null, 95)),
        n_motifs_beating_null_p05=int((emp_p < 0.05).sum()),
        value_column=valcol, description=vlabel,
        n_shared_keys=int(len(shared)), n_mofaflex_keys=int(Mg.shape[0]),
        n_alarmist_keys=int(Ag.shape[0]),
        alarmist_mass_retained=float(kept),
        max_cosine=float(C.values.max()), mean_cosine=float(C.values.mean()),
        n_motifs_with_match_above_0p5=int((C.max(0) > 0.5).sum()),
        best_per_motif={c: [C[c].idxmax(), float(C[c].max())] for c in C.columns},
    )
    print(f"\n=== ALARMIST value column: {valcol}  ({vlabel}) ===")
    print(f"  shared (sender,L,R) keys : {len(shared)}  "
          f"[MOFA-Flex {Mg.shape[0]}, ALARMIST {Ag.shape[0]}]")
    print(f"  ALARMIST mass retained after collapse+join : {100*kept:.1f}%")
    print(f"  max cosine {C.values.max():.3f} | mean {C.values.mean():.3f} | "
          f"motifs with a match > 0.5: {int((C.max(0)>0.5).sum())}/{C.shape[1]}")

    VMAX = max(0.6, float(C.values.max()))
    if a.sign_mode == "signed":
        _lim = float(np.abs(C.values).max())
        CMAP, VMIN, VMAX = "RdBu_r", -_lim, _lim      # diverging, white at 0
    else:
        CMAP, VMIN = "Reds", 0.0

    # (a) sorted heatmap -- rows/cols by their own best match, so the strongest pairings sit
    #     top-left and the reader can see whether ANY factor claims a motif
    order_f = C.max(1).sort_values(ascending=False).index
    order_m = C.max(0).sort_values(ascending=False).index
    fig, ax = plt.subplots(figsize=(1 + 0.42 * C.shape[1], 1 + 0.40 * C.shape[0]))
    sns.heatmap(C.loc[order_f, order_m], cmap=CMAP, vmin=VMIN, vmax=VMAX,
                ax=ax, cbar_kws={"label": "cosine similarity"},
                linewidths=0.3, linecolor="white")
    ax.set_title(f"MOFA-Flex factors vs ALARMIST motifs  (sorted by best match)\n"
                 f"{len(shared)} shared (sender, ligand, receptor) keys | "
                 f"ALARMIST scored on {vlabel}", fontsize=9)
    ax.set_xlabel("ALARMIST motif"); ax.set_ylabel("MOFA-Flex factor pole" if a.sign_mode=="poles" else "MOFA-Flex factor")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / "figures" / f"cosine_{a.label}_{a.sign_mode}_vs_alarmist_{valcol}.{ext}", dpi=200)
    plt.close(fig)

    # (b) clustermap -- hierarchical clustering on BOTH axes. If the two methods recovered the
    #     same programs there should be a permutation of rows/cols with a bright diagonal;
    #     clustering finds that permutation if it exists, so a clustermap with no diagonal
    #     block structure is stronger evidence of non-correspondence than a sorted heatmap.
    g = sns.clustermap(C, cmap=CMAP, vmin=VMIN, vmax=VMAX,
                       figsize=(2.2 + 0.44 * C.shape[1], 2.2 + 0.42 * C.shape[0]),
                       linewidths=0.3, linecolor="white",
                       cbar_kws={"label": "cosine similarity"},
                       dendrogram_ratio=(0.14, 0.14),
                       method="average", metric="euclidean")
    g.ax_heatmap.set_xlabel("ALARMIST motif"); g.ax_heatmap.set_ylabel("MOFA-Flex factor pole" if a.sign_mode=="poles" else "MOFA-Flex factor")
    g.fig.suptitle(f"MOFA-Flex factors vs ALARMIST motifs  (clustered both axes)\n"
                   f"{len(shared)} shared (sender, ligand, receptor) keys | "
                   f"ALARMIST scored on {vlabel}", fontsize=9, y=1.0)
    for ext in ("png", "pdf", "svg"):
        g.fig.savefig(OUT / "figures" / f"clustermap_{a.label}_{a.sign_mode}_vs_alarmist_{valcol}.{ext}",
                      dpi=200, bbox_inches="tight")
    plt.close(g.fig)
    print(f"  wrote figures/{{cosine,clustermap}}_{a.label}_{a.sign_mode}_vs_alarmist_{valcol}.png/pdf/svg")

results["raw_key_overlap_before_canonicalisation"] = int(raw_overlap)
json.dump(results, open(OUT / f"cosine_{a.label}_{a.sign_mode}_vs_alarmist.json", "w"), indent=2, default=str)
print(f"\nraw (sender,L,R) key overlap BEFORE canonicalising subunit order: {raw_overlap}")
print(f"wrote {OUT}/cosine_{a.label}_vs_alarmist.json")
