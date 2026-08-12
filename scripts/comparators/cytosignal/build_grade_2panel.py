#!/usr/bin/env python
# 2-panel grade comparison (styled like al.glm_volcano):
#  (1) CytoSignal runNEBULA (native NB mixed model) per-LRI volcano
#  (2) ALARMIST plot_motif_fraction_by_grade (motif-ON fraction per core, RAW Mann-Whitney p)
# Panel 1 labels: top-10 up + top-10 down (by significance) + the SIGNIFICANT motif-1 LRIs (orange).
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common.plotting import apply_publication_style

apply_publication_style()

NG = "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM/cellchatdb2/nebula_grade"
PD = "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM/cellchatdb2/run_full/plots"
RED, BLUE, ORANGE = "#E64B35", "#3C5488", "#E69F00"
clean = lambda s: str(s).replace("-ligand", "").replace("-receptor", "")

neb = pd.read_csv(NG + "/nebula_grade_results.csv").dropna(subset=["p", "logFC"]).copy()
m1  = pd.read_csv(NG + "/motif1_top_lris_grade.csv")
m1_top100 = set(m1.loc[m1["alarmist_rank"] <= 100, "cs_name"].dropna())

FDR, XT, NTOP = 0.05, 0.2, 10
fig, ax = plt.subplots(1, 2, figsize=(24, 12))

# ---------- Panel 1: CytoSignal runNEBULA volcano ----------
a = ax[0]
x = neb["logFC"].values; y = -np.log10(np.clip(neb["padj"].values, 1e-300, None))
sig = y >= -np.log10(FDR); eff = np.abs(x) >= XT
mid = ~eff; blk = eff & ~sig; up = sig & (x >= XT); dn = sig & (x <= -XT)
a.scatter(x[mid], y[mid], c="lightgray", alpha=0.25, s=16, linewidths=0)
a.scatter(x[blk], y[blk], c="black", alpha=0.55, s=16, linewidths=0)
a.scatter(x[up], y[up], c=RED, alpha=0.85, s=26, linewidths=0)
a.scatter(x[dn], y[dn], c=BLUE, alpha=0.85, s=26, linewidths=0)
neb["ism1"] = neb["interaction"].isin(m1_top100)
mm = neb["ism1"].values
a.scatter(x[mm], y[mm], facecolors="none", edgecolors=ORANGE, s=70, linewidths=1.3, label="motif-1 top-100 LRI")
a.axvline(XT, ls="--", c="k", alpha=0.4); a.axvline(-XT, ls="--", c="k", alpha=0.4); a.axhline(-np.log10(FDR), ls="--", c="k", alpha=0.4)

sig_all = sig & eff                          # for red/blue colouring (q<0.05 & |logFC|>XT)
n_sig_fdr = int(sig.sum())                   # runNEBULA-official significant count: q<0.05 (per-mode BH)
n_m1_tested = int(neb.ism1.sum())            # motif-1 top-100 LRIs actually tested (nebula's cpc filter drops some)
m1_sig = neb[neb.ism1 & sig]; n_m1_sig = len(m1_sig)   # significant motif-1 LRIs (q<0.05)
# labels: top-10 up + top-10 down by significance (padj), + significant motif-1 LRIs
up_lab = neb[neb.logFC > 0].nsmallest(NTOP, "padj")
dn_lab = neb[neb.logFC < 0].nsmallest(NTOP, "padj")
lab = pd.concat([up_lab, dn_lab, m1_sig]).drop_duplicates("interaction")
texts = []
for _, r in lab.iterrows():
    yi = -np.log10(max(r["padj"], 1e-300)); col = ORANGE if r["ism1"] else "black"
    texts.append(a.text(r["logFC"], yi, clean(r["interaction"]), fontsize=8, ha="center", va="center", color=col,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=col if r["ism1"] else "gray", alpha=0.9)))
adjust_text(texts, force_text=(0.5, 0.6), arrowprops=dict(arrowstyle="-", color="gray", lw=0.5), ax=a)
a.set_xlabel("logFC (high vs low grade)", fontsize=15); a.set_ylabel(r"-$\log_{10}(q)$", fontsize=15)
a.set_title(f"(1) CytoSignal runNEBULA — native NB mixed model (cpc=0.001, per-mode BH)\n"
            f"{len(neb)} LRIs tested · {n_sig_fdr} significant (q<{FDR}) · "
            f"motif-1: {n_m1_sig} sig of {n_m1_tested} tested (top-100)", fontsize=12)
a.legend(loc="upper left", fontsize=10)

# ---------- Panel 2: ALARMIST motif fraction (raw p) ----------
al = pd.read_csv(NG + "/alarmist_motif_fraction_grade.csv"); al["d"] = al.frac_high - al.frac_low
b = ax[1]; yv = -np.log10(al.p.clip(1e-300))
sigm = al.p < 0.05
b.scatter(al.d[~sigm], yv[~sigm], c="lightgray", s=90, linewidths=0, label="motif p≥0.05")
b.scatter(al.d[sigm], yv[sigm], c=RED, s=110, linewidths=0, label="motif p<0.05")
mo1 = al[al.motif == 1].iloc[0]
b.scatter([mo1.d], [-np.log10(mo1.p)], marker="*", s=420, c=ORANGE, edgecolors="k", linewidths=0.6, zorder=6, label="motif 1 (mGAM)")
b.axhline(-np.log10(0.05), ls="--", c="k", alpha=0.4)
texts2 = [b.text(r.d, -np.log10(r.p), f"motif {int(r.motif)}", fontsize=9, ha="center", va="center",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="gray", alpha=0.85)) for _, r in al.iterrows()]
adjust_text(texts2, force_text=(0.4, 0.5), arrowprops=dict(arrowstyle="-", color="gray", lw=0.5), ax=b)
b.set_xlabel("Δ ON-cell fraction (high − low grade)", fontsize=15); b.set_ylabel(r"-$\log_{10}(p)$ (raw Mann-Whitney)", fontsize=15)
b.set_title(f"(2) ALARMIST plot_motif_fraction_by_grade\n{len(al)} motifs · {int(sigm.sum())} with p<0.05 · motif 1 p={mo1.p:.3f}", fontsize=13)
b.legend(loc="lower right", fontsize=10)

fig.suptitle("High- vs low-grade GBM (same 13 TMA cores). CytoSignal's native runNEBULA finds grade-differential LRIs — mostly generic "
             "(CDH2, junctions, SLC-glutamate→GRM3), not mGAM;\nof motif-1's top-100 LRIs (%d tested by nebula) %d are significant. "
             "ALARMIST flags the mGAM motif itself (motif 1, p=0.014)." % (n_m1_tested, n_m1_sig),
             fontsize=12.5)
out = f"{PD}/grade_comparison_2panel.png"
plt.tight_layout(rect=[0, 0, 1, 0.95]); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
print("wrote", out)
print(f"runNEBULA: {len(neb)} LRIs, {n_sig_fdr} sig (q<{FDR}); motif-1 top-100: {n_m1_sig} sig of {n_m1_tested} tested:",
      list(m1_sig["interaction"]))
print("up-labeled:", [clean(s) for s in up_lab.interaction.head(10)])
print("down-labeled:", [clean(s) for s in dn_lab.interaction.head(10)])
