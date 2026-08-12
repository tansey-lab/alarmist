#!/usr/bin/env python
# Spatial CytoSignal score maps for ALARMIST motif 1's top-25 LR pairs (by factor_lrnorm),
# de-duplicated to unique ligand-receptor (CytoSignal is cell-type-agnostic, so one map per L-R).
import numpy as np, pandas as pd, scipy.io as sio, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

H5 = "/Users/jiayifan/tansey_lab/alarmist/data/xenium_mm_final_cell_id.h5ad"
LM = "/Users/jiayifan/tansey_lab/alarmist/results/GBM/bptf/lri_motifs.csv"
ROOT = "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM"
QD, PD = ROOT+"/cellchatdb2/run_full/quant", ROOT+"/cellchatdb2/run_full/plots"
WCOL, TOPN = "factor_lrnorm", 25

m1 = pd.read_csv(LM); m1 = m1[m1["motif_idx"] == 1]
idx = m1.groupby(["ligand", "receptor"])[WCOL].idxmax()          # best row per unique L-R
uniq = m1.loc[idx].sort_values(WCOL, ascending=False).head(TOPN).reset_index(drop=True)

adata = ad.read_h5ad(H5, backed="r")
cid = adata.obs["cell_id"].astype(str).values; xy = adata.obsm["spatial"][:]
pos = {c: i for i, c in enumerate(cid)}

def load_mode(slot):
    cells = [l.strip() for l in open(f"{QD}/score_{slot}.cells.tsv")]
    intr  = [l.strip() for l in open(f"{QD}/score_{slot}.intr.tsv")]
    S = sio.mmread(f"{QD}/score_{slot}.mtx").tocsc()
    summ = pd.read_csv(f"{QD}/signif_summary_{slot}.csv").sort_values("n_hq", ascending=False).reset_index(drop=True)
    id2col = {c: i for i, c in enumerate(intr)}
    name2col = {nm: id2col[iid] for nm, iid in zip(summ.name, summ.interaction_id) if iid in id2col}
    name2rank = {nm: i + 1 for i, nm in enumerate(summ.name)}
    name2n = dict(zip(summ.name, summ.n_hq))
    return cells, S, name2col, name2rank, name2n, len(summ)

dcells, Sd, nd, rd, hd, Nd = load_mode("diffusion_Raw_smooth")
ccells, Sc, nc, rc, hc, Nc = load_mode("contact_Raw_smooth")
cells = dcells
creorder = np.array([ccells.index(c) for c in cells]) if ccells != dcells else np.arange(len(cells))
srows = np.array([pos[c] for c in cells]); xy_s = xy[srows]

def cyto(name, mode):   # returns (score, cytosignal-rank-label) — native mode first, then fallback
    seq = [("cont", nc, Sc, creorder, rc, hc, Nc), ("diff", nd, Sd, None, rd, hd, Nd)] if mode == "juxtacrine" \
          else [("diff", nd, Sd, None, rd, hd, Nd), ("cont", nc, Sc, creorder, rc, hc, Nc)]
    for tag, tbl, mat, reorder, rk, hh, N in seq:
        j = tbl.get(name)
        if j is not None:
            v = np.asarray(mat[:, j].todense()).ravel()
            v = v if reorder is None else v[reorder]
            return v, f"CytoSignal #{rk[name]}/{N} ({tag}), {int(hh[name])} cells"
    return None, None

ncol = 5; nrow = int(np.ceil(TOPN/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.3*ncol, 5.0*nrow)); axes = axes.ravel()
for i, (ax, (_, r)) in enumerate(zip(axes, uniq.iterrows())):
    name = f"{r.ligand} - {r.receptor}"; v, cslabel = cyto(name, r.signaling_type)
    if v is None:
        ax.axis("off"); ax.set_title(f"ALARMIST #{i+1}  {name}\n(not scored by CytoSignal)", fontsize=8); continue
    vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1.0
    o = np.argsort(v)
    s = ax.scatter(xy_s[o, 0], xy_s[o, 1], c=v[o], s=0.8, cmap="magma", vmin=0, vmax=vmax, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    ax.set_title(f"{name}\nALARMIST #{i+1} (motif1): {r.celltype1}→{r.celltype2} {r.signaling_type[:4]} lrnorm={r[WCOL]:.0f}\n"
                 f"{cslabel}", fontsize=7.2)
    plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02)
for ax in axes[len(uniq):]: ax.axis("off")
fig.suptitle(f"CytoSignal spatial score for ALARMIST motif 1's top-{TOPN} LR pairs (by {WCOL}, unique L-R).\n"
             f"ALARMIST groups + directions these into one mGAM motif; CytoSignal gives them as {TOPN} separate undirected maps.",
             fontsize=13, y=0.998)
out = f"{PD}/motif1_top{TOPN}_lris_cytosignal.png"
plt.tight_layout(rect=[0, 0, 1, 0.975]); plt.savefig(out, dpi=145, bbox_inches="tight"); plt.close()
print("wrote", out, "| panels:", len(uniq))
print(uniq[["ligand", "receptor", "celltype1", "celltype2", "signaling_type", WCOL]].to_string())
