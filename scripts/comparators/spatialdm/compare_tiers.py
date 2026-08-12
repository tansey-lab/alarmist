#!/usr/bin/env python
"""Compare SpatialDM's two LR-database tiers on the same cores.

`default`     = bundled CellChatDB **v1** (1,939 interactions, 3 signaling categories)
`cellchatdb2` = CellChatDB **v2** injected by run_spatialdm.py (3,233 rows / 3,218 unique pairs),
                including 994 `Non-protein Signaling` rows that v1 does not have and that we had
                to remap onto the long-range RBF kernel for SpatialDM to run at all.

The point of the two tiers (comparator-benchmark skill) is to remove the LR resource as a
confounder. This script answers: how much of the cellchatdb2 result survives when the database
is the authors' own, and how much of it was the Non-protein remap?

Pairs are matched across tiers by the **resolved subunit sets**, not by interaction name --
v1 names complexes its own way (`TGFB1_TGFBR1_TGFBR2`) while our v2 names are built as
`ligand_receptor` from the complex table, so string joins would silently under-match.

Usage: python compare_tiers.py --root results/comparators/spatialdm/GBM \
           --tier-a default --tier-b cellchatdb2 --out-dir results/comparators/spatialdm/GBM/tier_comparison
"""
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common.plotting import apply_publication_style

apply_publication_style()

p = argparse.ArgumentParser()
p.add_argument("--root", required=True, help="dir holding the tier subdirs")
p.add_argument("--tier-a", default="default")
p.add_argument("--tier-b", default="cellchatdb2")
p.add_argument("--out-dir", required=True)
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
p.add_argument("--v2-csv", default="data/LRdatabase/CellChatDBv2.0.human.csv")
a = p.parse_args()
OUT = a.out_dir
os.makedirs(OUT, exist_ok=True)


def save_all(stem, fig=None):
    fig = fig or plt.gcf()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def load(tier, core):
    """global_res keyed by the resolved (ligand subunits | receptor subunits) signature."""
    d = os.path.join(a.root, tier, core, "data")
    gr = pd.read_csv(os.path.join(d, "global_res.csv"), index_col=0)
    lig = pd.read_csv(os.path.join(d, "lr_ligand_subunits.csv"), index_col=0)
    rec = pd.read_csv(os.path.join(d, "lr_receptor_subunits.csv"), index_col=0)
    subs = lambda f: f.reindex(gr.index).apply(
        lambda r: "_".join(sorted(r.dropna().astype(str))), axis=1)
    sig = subs(lig) + "|" + subs(rec)
    gr = gr.assign(sig=sig, name=gr.index.astype(str))
    # a signature can appear twice if the DB lists the same resolved pair under two names;
    # keep the most significant so the join stays 1:1
    return gr.sort_values("fdr").drop_duplicates("sig").set_index("sig")


def nrows(tier, core):
    """rows in global_res BEFORE collapsing duplicate resolved subunit pairs"""
    return len(pd.read_csv(os.path.join(a.root, tier, core, "data", "global_res.csv"), index_col=0))


cores = sorted([c for c in os.listdir(os.path.join(a.root, a.tier_a))
                if c.isdigit() and os.path.exists(
                    os.path.join(a.root, a.tier_a, c, "data", "global_res.csv"))],
               key=int)
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]

# Which tier-B pairs came from the remapped Non-protein category. Match on the v2
# INTERACTION NAME, which is exactly the key run_spatialdm.py built the DB with -- matching on
# resolved subunits would silently miss every pair whose complex has an off-panel subunit,
# because the stored subunits are filtered to the 5,119-gene panel and the raw DB's are not.
np_names = set()
if os.path.exists(a.v2_csv):
    v2 = pd.read_csv(a.v2_csv)
    v2["interaction_name"] = v2.ligand.astype(str) + "_" + v2.receptor.astype(str)
    np_names = set(v2.loc[v2.signaling_type == "Non-protein Signaling", "interaction_name"])

rows, per_core_shared = [], {}
for c in cores:
    A, B, rawA, rawB = (*[load(t, c) for t in (a.tier_a, a.tier_b)],
                        *[nrows(t, c) for t in (a.tier_a, a.tier_b)])
    shared = A.index.intersection(B.index)
    selA, selB = set(A.index[A.selected.astype(bool)]), set(B.index[B.selected.astype(bool)])
    sA, sB = selA & set(shared), selB & set(shared)
    both, onlyA, onlyB = len(sA & sB), len(sA - sB), len(sB - sA)
    r = np.nan
    if len(shared) > 2:
        r = float(np.corrcoef(A.loc[shared, "z"], B.loc[shared, "z"])[0, 1])
    npB = sum(B.loc[s, "name"] in np_names for s in selB)
    rows.append(dict(
        core=int(c),
        # rows in the tier's global_res, and how many DISTINCT resolved subunit pairs they
        # collapse to on this gene panel -- v2 lists several DB rows that become the same
        # measurable pair once off-panel subunits are dropped.
        n_rows_A=rawA, n_rows_B=rawB, n_valid_A=len(A), n_valid_B=len(B),
        pct_redundant_A=round(100 * (1 - len(A) / rawA), 1),
        pct_redundant_B=round(100 * (1 - len(B) / rawB), 1),
        n_shared=len(shared), pct_B_shared=round(100 * len(shared) / len(B), 1),
        n_sel_A=len(selA), n_sel_B=len(selB),
        n_sel_B_non_protein=npB,
        pct_sel_B_non_protein=round(100 * npB / max(len(selB), 1), 1),
        sel_shared_both=both, sel_only_A=onlyA, sel_only_B=onlyB,
        jaccard_on_shared=round(both / max(both + onlyA + onlyB, 1), 3),
        z_corr_on_shared=round(r, 3) if r == r else None))
    per_core_shared[c] = (A, B, shared)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "per_core_tier_comparison.csv"), index=False)
pd.set_option("display.width", 250)
print(df.to_string(index=False))

# ---- requested LRs side by side -------------------------------------------------
req_rows = []
for c in cores:
    A, B, _ = per_core_shared[c]
    an, bn = A.reset_index().set_index("name"), B.reset_index().set_index("name")
    for lr in requested:
        req_rows.append(dict(
            core=int(c), lr=lr,
            v1_testable=lr in an.index, v2_testable=lr in bn.index,
            v1_fdr=float(an.loc[lr, "fdr"]) if lr in an.index else None,
            v2_fdr=float(bn.loc[lr, "fdr"]) if lr in bn.index else None,
            v1_selected=bool(an.loc[lr, "selected"]) if lr in an.index else False,
            v2_selected=bool(bn.loc[lr, "selected"]) if lr in bn.index else False))
rq = pd.DataFrame(req_rows)
rq.to_csv(os.path.join(OUT, "requested_lr_by_tier.csv"), index=False)
print("\nRequested LRs — selected in how many of the 13 cores:")
print(rq.groupby("lr")[["v1_testable", "v2_testable", "v1_selected", "v2_selected"]].sum().to_string())

# ---- figures ---------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
ax[0].scatter(df.n_sel_A, df.n_sel_B, s=45, c="#4C72B0")
lim = [0, max(df.n_sel_A.max(), df.n_sel_B.max()) * 1.08]
ax[0].plot(lim, lim, ls="--", c="grey", lw=1); ax[0].set_xlim(lim); ax[0].set_ylim(lim)
for _, r in df.iterrows():
    ax[0].annotate(str(r.core), (r.n_sel_A, r.n_sel_B), fontsize=7,
                   xytext=(3, 3), textcoords="offset points")
ax[0].set_xlabel(f"significant pairs — {a.tier_a} (v1)")
ax[0].set_ylabel(f"significant pairs — {a.tier_b} (v2)")
ax[0].set_title("per-core significant counts")

ax[1].bar(df.core.astype(str), df.jaccard_on_shared, color="#55A868")
ax[1].set_ylim(0, 1); ax[1].set_xlabel("core")
ax[1].set_ylabel("Jaccard of selected sets\n(pairs present in both DBs)")
ax[1].set_title("agreement on the shared pairs")

ax[2].bar(df.core.astype(str), df.pct_sel_B_non_protein, color="#C44E52")
ax[2].set_xlabel("core"); ax[2].set_ylabel("% of v2 significant pairs\nfrom Non-protein Signaling")
ax[2].set_ylim(0, 100); ax[2].set_title("what the remap contributes")
fig.tight_layout()
save_all("tier_comparison_overview", fig)

# z-score concordance, pooled over cores
xs = np.concatenate([per_core_shared[c][0].loc[per_core_shared[c][2], "z"].values for c in cores])
ys = np.concatenate([per_core_shared[c][1].loc[per_core_shared[c][2], "z"].values for c in cores])
fig, ax = plt.subplots(figsize=(4.6, 4.4))
ax.scatter(xs, ys, s=4, alpha=.25, c="#4C72B0", edgecolors="none")
lo, hi = float(min(xs.min(), ys.min())), float(max(xs.max(), ys.max()))
ax.plot([lo, hi], [lo, hi], ls="--", c="grey", lw=1)
ax.set_xlabel(f"Moran z — {a.tier_a} (v1)"); ax.set_ylabel(f"Moran z — {a.tier_b} (v2)")
ax.set_title(f"shared pairs, all cores (n={len(xs):,})\nPearson r = {np.corrcoef(xs, ys)[0,1]:.3f}")
fig.tight_layout()
save_all("z_concordance_shared_pairs", fig)

summary = dict(
    tier_a=a.tier_a, tier_b=a.tier_b, cores=[int(c) for c in cores],
    total_rows_A=int(df.n_rows_A.sum()), total_rows_B=int(df.n_rows_B.sum()),
    total_valid_A=int(df.n_valid_A.sum()), total_valid_B=int(df.n_valid_B.sum()),
    pct_redundant_A=round(100 * (1 - df.n_valid_A.sum() / df.n_rows_A.sum()), 1),
    pct_redundant_B=round(100 * (1 - df.n_valid_B.sum() / df.n_rows_B.sum()), 1),
    total_shared=int(df.n_shared.sum()),
    pooled_sel_A=int(df.n_sel_A.sum()), pooled_sel_B=int(df.n_sel_B.sum()),
    pooled_sel_B_non_protein=int(df.n_sel_B_non_protein.sum()),
    pct_pooled_sel_B_non_protein=round(100 * df.n_sel_B_non_protein.sum() / df.n_sel_B.sum(), 1),
    median_jaccard_on_shared=float(df.jaccard_on_shared.median()),
    pooled_z_corr=round(float(np.corrcoef(xs, ys)[0, 1]), 3),
    median_pct_of_v2_pairs_present_in_v1=float(df.pct_B_shared.median()))
json.dump(summary, open(os.path.join(OUT, "tier_comparison_summary.json"), "w"), indent=2)
print("\n" + json.dumps(summary, indent=2))
