#!/usr/bin/env python
# Reconstruct ALARMIST motif 1 (mGAM) from CytoSignal: for motif 1's top-N LRIs, sum each LRI's
# CytoSignal per-cell score weighted by its ALARMIST factor (two versions: factor, factor_lrnorm),
# and plot side-by-side with the real motif-1 cell loading. Tests whether CytoSignal's scores
# CONTAIN the motif signal — recoverable only if you already have ALARMIST's motif weights.
import numpy as np, pandas as pd, scipy.io as sio, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

H5 = "/Users/jiayifan/tansey_lab/alarmist/data/xenium_mm_final_cell_id.h5ad"
LOADINGS = "/Users/jiayifan/tansey_lab/alarmist/results/GBM/single_cell/cell_loadings.npy"
LM = "/Users/jiayifan/tansey_lab/alarmist/results/GBM/bptf/lri_motifs.csv"
ROOT = "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM"
QD, PD = ROOT+"/cellchatdb2/run_full/quant", ROOT+"/cellchatdb2/run_full/plots"
TOPN, MGAM_COL = 100, 1

m1 = pd.read_csv(LM); m1 = m1[m1["motif_idx"] == 1].copy()
adata = ad.read_h5ad(H5, backed="r")
cid = adata.obs["cell_id"].astype(str).values; xy = adata.obsm["spatial"][:]
mgam = np.load(LOADINGS)[:, MGAM_COL]
pos = {c: i for i, c in enumerate(cid)}

def load_mode(slot):
    cells = [l.strip() for l in open(f"{QD}/score_{slot}.cells.tsv")]
    intr  = [l.strip() for l in open(f"{QD}/score_{slot}.intr.tsv")]
    S = sio.mmread(f"{QD}/score_{slot}.mtx").tocsc()
    summ = pd.read_csv(f"{QD}/signif_summary_{slot}.csv")
    id2col = {cid_: i for i, cid_ in enumerate(intr)}
    name2col = {nm: id2col[iid] for nm, iid in zip(summ.name, summ.interaction_id) if iid in id2col}
    return cells, S, name2col

dcells, Sd, nd = load_mode("diffusion_Raw_smooth")
ccells, Sc, nc = load_mode("contact_Raw_smooth")
cells = dcells                                   # canonical order
creorder = np.array([ccells.index(c) for c in cells]) if ccells != dcells else np.arange(len(cells))
srows = np.array([pos[c] for c in cells])
xy_s = xy[srows]; mgam_s = mgam[srows]

def cyto(name, mode):
    diff = (name, nd, Sd, None); cont = (name, nc, Sc, creorder)
    order = [cont, diff] if mode == "juxtacrine" else [diff, cont]   # try native mode first, then fallback
    for nm, tbl, mat, reorder in order:
        j = tbl.get(nm)
        if j is not None:
            v = np.asarray(mat[:, j].todense()).ravel()
            return v if reorder is None else v[reorder]
    return None

def reconstruct(wcol, standardize=False):
    top = m1.sort_values(wcol, ascending=False).head(TOPN)
    recon = np.zeros(len(cells)); matched = 0; mw = 0.0; tw = 0.0; miss = []
    for _, r in top.iterrows():
        tw += r[wcol]; name = f"{r.ligand} - {r.receptor}"
        v = cyto(name, r.signaling_type)
        if v is None: miss.append(name); continue
        if standardize:
            sd = v.std(); v = (v - v.mean()) / sd if sd > 0 else v * 0.0
        recon += r[wcol] * v; matched += 1; mw += r[wcol]
    return recon, matched, mw / tw, miss

def pearson(a, b): return float(np.corrcoef(a, b)[0, 1])

recons = {}
for wcol in ["factor", "factor_lrnorm"]:
    for std in [False, True]:
        recon, matched, wfrac, miss = reconstruct(wcol, std)
        r = pearson(recon, mgam_s); recons[(wcol, std)] = (recon, r)
        print(f"[{wcol} {'z-scored' if std else 'raw':8s}] matched {matched}/{TOPN} "
              f"({100*wfrac:.0f}% weight); Pearson r = {r:.3f}")

def spatial(ax, val, title):
    lo, hi = np.percentile(val, 2), np.percentile(val, 98)
    if (val >= 0).all(): lo = 0.0
    o = np.argsort(val)
    s = ax.scatter(xy_s[o, 0], xy_s[o, 1], c=val[o], s=1.0, cmap="magma", vmin=lo, vmax=hi, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off"); ax.set_title(title, fontsize=9.5)
    plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
for row, std in enumerate([False, True]):
    tag = "z-scored per-LRI" if std else "raw score"
    spatial(axes[row, 0], mgam_s, "ALARMIST motif 1 (mGAM)\ncell loading (ground truth)")
    spatial(axes[row, 1], recons[("factor", std)][0],
            f"CytoSignal recon · weight=factor · {tag}\ntop{TOPN} LRIs; r={recons[('factor', std)][1]:.2f}")
    spatial(axes[row, 2], recons[("factor_lrnorm", std)][0],
            f"CytoSignal recon · weight=factor_lrnorm · {tag}\ntop{TOPN} LRIs; r={recons[('factor_lrnorm', std)][1]:.2f}")
fig.suptitle(f"Reconstruct ALARMIST motif 1 from CytoSignal: Σ (ALARMIST factor × CytoSignal LRI score) over motif 1's "
             f"top-{TOPN} LRIs.\nThe LRIs' scores exist in CytoSignal, but the motif is recoverable only WITH ALARMIST's "
             f"weights — CytoSignal alone never groups them.", fontsize=12)
out = f"{PD}/reconstruct_motif1_from_cytosignal.png"
plt.tight_layout(); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
print("wrote", out)
