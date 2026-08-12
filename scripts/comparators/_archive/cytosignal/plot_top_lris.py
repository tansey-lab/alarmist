#!/usr/bin/env python
# Spatial score-map grids of CytoSignal's TOP-ranked LRIs (per mode) from a run's quant/ dir.
# Usage: python plot_top_lris.py [run_dir] [input_dir]
#   run_dir   : has quant/ (score_<slot>.mtx/.cells.tsv/.intr.tsv, signif_summary_<slot>.csv)
#   input_dir : has meta.csv (cell_id,x,y,celltype)
import os, sys, numpy as np, pandas as pd, scipy.io as sio
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RUN = sys.argv[1] if len(sys.argv) > 1 else "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM/cellchatdb2/run_full"
INP = sys.argv[2] if len(sys.argv) > 2 else "/Users/jiayifan/tansey_lab/alarmist/results/comparators/cytosignal/GBM/input_full"
QD, PD = RUN+"/quant", RUN+"/plots"; os.makedirs(PD, exist_ok=True)
meta = pd.read_csv(INP+"/meta.csv", dtype={"cell_id": str}).set_index("cell_id")

def load_slot(slot):
    cells = [l.strip() for l in open(f"{QD}/score_{slot}.cells.tsv")]
    intr  = [l.strip() for l in open(f"{QD}/score_{slot}.intr.tsv")]
    S = sio.mmread(f"{QD}/score_{slot}.mtx").tocsc()
    summ = pd.read_csv(f"{QD}/signif_summary_{slot}.csv").sort_values("n_hq", ascending=False).reset_index(drop=True)
    id2col = {cid: i for i, cid in enumerate(intr)}
    s = meta.loc[cells]
    return dict(S=S, summ=summ, id2col=id2col, x=s.x.values, y=s.y.values)

def col(slot, cid):
    return np.asarray(slot["S"][:, slot["id2col"][cid]].todense()).ravel()

def grid(slot, rows_df, title, fname, ncol=4):
    n = len(rows_df); nrow = int(np.ceil(n/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3*ncol, 4.3*nrow)); axes = np.atleast_1d(axes).ravel()
    for ax, (_, r) in zip(axes, rows_df.iterrows()):
        cid = r["interaction_id"]
        if cid not in slot["id2col"]: ax.axis("off"); continue
        v = col(slot, cid); vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1.0
        o = np.argsort(v)
        sc = ax.scatter(slot["x"][o], slot["y"][o], c=v[o], s=1.0, cmap="magma", vmin=0, vmax=vmax, linewidths=0)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
        ax.set_title(f"#{r['rank']}  {r['name']}\n{int(r['n_hq'])} cells", fontsize=8)
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    for ax in axes[n:]: ax.axis("off")
    fig.suptitle(title, fontsize=13)
    plt.tight_layout(); plt.savefig(f"{PD}/{fname}", dpi=140, bbox_inches="tight"); plt.close()
    print("wrote", fname, f"({n} panels)")

for slot_name, nm_top, ncol in [("diffusion_Raw_smooth", 16, 4), ("contact_Raw_smooth", 12, 4)]:
    slot = load_slot(slot_name)
    summ = slot["summ"].copy(); summ["rank"] = np.arange(1, len(summ)+1)
    mode = slot_name.split("_")[0]
    grid(slot, summ.head(nm_top), f"CytoSignal GBM — top {nm_top} {mode} LRIs by significant-cell count",
         f"top_{mode}_grid.png", ncol=ncol)
print("DONE ->", PD)
