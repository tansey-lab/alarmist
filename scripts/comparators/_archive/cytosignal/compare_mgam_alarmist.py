#!/usr/bin/env python
# Head-to-head: ALARMIST motif 1 (mGAM) program  vs  CytoSignal's per-LRI scores, same GBM cells.
import os, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

H5   = "/Users/jiayifan/tansey_lab/alarmist/data/xenium_mm_final_cell_id.h5ad"
ROOT = "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM"
LOADINGS = "/Users/jiayifan/tansey_lab/alarmist/results/GBM/single_cell/cell_loadings.npy"
QD, PD = ROOT+"/cellchatdb2/run_full/quant", ROOT+"/cellchatdb2/run_full/plots"; os.makedirs(PD, exist_ok=True)
MGAM_COL = 1  # 0-based column in cell_loadings = mGAM motif (notebook motif_idx=1; verified: mGAM's top-loading motif)
GENES = ["GRN", "SORT1", "ANXA1", "FPR1"]
CTS = ["mGAM", "non-mGAM", "MES-like", "AC-like", "NPC-like", "OPC-like", "Glial-Neuronal", "Vascular", "Lymphoid"]

A = ad.read_h5ad(H5); obs = A.obs
cid_all = obs["cell_id"].astype(str).values
xy_all = A.obsm["spatial"]; ct_all = obs["cell_type"].astype(str).values
cell_loadings = np.load(LOADINGS)                 # (n_cells x 20) projected motif loadings, same order as obs
assert cell_loadings.shape[0] == len(cid_all), (cell_loadings.shape, len(cid_all))
def gene_vec(g):
    col = A.X[:, A.var_names.get_loc(g)]
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
expr = {g: gene_vec(g) for g in GENES}

# CytoSignal diffusion scores aligned to score cells
cells = [l.strip() for l in open(f"{QD}/score_diffusion_Raw_smooth.cells.tsv")]
intr  = [l.strip() for l in open(f"{QD}/score_diffusion_Raw_smooth.intr.tsv")]
S = sio.mmread(f"{QD}/score_diffusion_Raw_smooth.mtx").tocsc()
summ = pd.read_csv(f"{QD}/signif_summary_diffusion_Raw_smooth.csv").sort_values("n_hq", ascending=False).reset_index(drop=True)
summ["rank"] = np.arange(1, len(summ)+1); NTOT = len(summ)
n2id = dict(zip(summ.name, summ.interaction_id)); n2rank = dict(zip(summ.name, summ["rank"]))
row_of = {c: i for i, c in enumerate(cid_all)}
srows = np.array([row_of[c] for c in cells])
xy_s = xy_all[srows]; ct_s = ct_all[srows]
def score(nm):
    return np.asarray(S[:, intr.index(n2id[nm])].todense()).ravel()
sc_grn = score("GRN - SORT1"); sc_anx = score("ANXA1 - FPR1")

def spatial(ax, x, y, val, title, cmap="magma", cbar=True):
    vmax = np.percentile(val[val > 0], 99) if (val > 0).any() else 1.0
    o = np.argsort(val)
    s = ax.scatter(x[o], y[o], c=val[o], s=1.0, cmap=cmap, vmin=0, vmax=vmax, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off"); ax.set_title(title, fontsize=10)
    if cbar: plt.colorbar(s, ax=ax, fraction=0.046, pad=0.02)

fig = plt.figure(figsize=(18, 11.5)); gs = GridSpec(2, 3, figure=fig, hspace=0.30, wspace=0.26)

# (a) ALARMIST motif-1 cell loading (real projected loadings; NOT the obs['motif'] argmax label)
spatial(fig.add_subplot(gs[0, 0]), xy_all[:, 0], xy_all[:, 1], cell_loadings[:, MGAM_COL],
        "ALARMIST motif 1 (mGAM) — cell loading\n(one program; al.project_cell_loadings)")

# (b,c) CytoSignal per-LRI score maps
spatial(fig.add_subplot(gs[0, 1]), xy_s[:, 0], xy_s[:, 1], sc_grn,
        f"CytoSignal: GRN→SORT1 score\nrank {n2rank['GRN - SORT1']}/{NTOT}  (one isolated pair)")
spatial(fig.add_subplot(gs[0, 2]), xy_s[:, 0], xy_s[:, 1], sc_anx,
        f"CytoSignal: ANXA1→FPR1 score\nrank {n2rank['ANXA1 - FPR1']}/{NTOT}  (one isolated pair)")

# (d) directionality heatmap: mean expr genes x cell types
ax = fig.add_subplot(gs[1, 0])
M = np.array([[expr[g][ct_all == c].mean() for c in CTS] for g in GENES])
im = ax.imshow(M, aspect="auto", cmap="viridis")
ax.set_xticks(range(len(CTS))); ax.set_xticklabels(CTS, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontsize=9)
for i in range(len(GENES)):
    j = int(np.argmax(M[i])); ax.text(j, i, "★", ha="center", va="center", color="w", fontsize=9)
ax.set_title("Expression by cell type (★=max): GRN & FPR1 are mGAM-specific —\n"
             "mGAM sends GRN, receives ANXA1 (the anchors ALARMIST gives direction)", fontsize=8.5)

# (e) ranking: mGAM LRIs buried
ax = fig.add_subplot(gs[1, 1]); top = summ.head(15)
ax.barh(range(len(top)), top.n_hq, color="#95a5a6"); ax.invert_yaxis()
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.name, fontsize=6.5)
ax.set_xlabel("significant cells")
ax.set_title(f"CytoSignal ranking: top 15 of {NTOT} (all generic WNT/BMP)\n"
             f"mGAM LRIs buried at #{n2rank['GRN - SORT1']} and #{n2rank['ANXA1 - FPR1']}", fontsize=8.5)

# (f) concordance: CytoSignal score by receiver cell type (the two directions land on different cells)
ax = fig.add_subplot(gs[1, 2])
grp = ["mGAM", "MES-like"]
g_grn = [sc_grn[ct_s == c].mean() for c in grp]; g_anx = [sc_anx[ct_s == c].mean() for c in grp]
w = 0.35; x = np.arange(2)
ax.bar(x - w/2, g_grn, w, label="GRN→SORT1", color="#e67e22")
ax.bar(x + w/2, g_anx, w, label="ANXA1→FPR1", color="#8e44ad")
ax.set_xticks(x); ax.set_xticklabels(grp); ax.set_ylabel("mean CytoSignal score")
ax.legend(fontsize=8)
ax.set_title("Both mGAM-LRI scores peak on mGAM cells (the shared hub) —\n"
             "yet delivered as 2 separate buried maps, never grouped into one motif", fontsize=8.5)

fig.suptitle("mGAM⇄MES-like signaling: ALARMIST recovers it as one directional motif (motif 1); "
             "CytoSignal scatters it into 2 undirected maps ranked #66 & #255 of 895", fontsize=13, y=0.99)
out = f"{PD}/comparison_mGAM_vs_cytosignal.png"
plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
print("wrote", out)
print("score by receiver:  GRN-SORT1  mGAM=%.3f MES=%.3f | ANXA1-FPR1 mGAM=%.3f MES=%.3f"
      % (g_grn[0], g_grn[1], g_anx[0], g_anx[1]))
