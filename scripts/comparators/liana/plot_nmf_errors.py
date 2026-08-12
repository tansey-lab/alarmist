#!/usr/bin/env python
"""A3: NMF reconstruction error vs k, from the saved nmf_errors.csv of both branches.

Plots only what li.multi.nmf already computed and stored -- no stability selection, no
resampling, no rank test of any kind that LIANA does not itself perform. The dashed line is
the rank kneedle chose; both branches are shown on a common k axis so they are comparable.
"""
import json, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = "results/comparators/liana/GBM"
runs = [("bivariate", f"{G}/nmf_bivariate"), ("inflow", f"{G}/nmf_inflow")]
fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
combined = {}
for i, (tag, d) in enumerate(runs):
    e = pd.read_csv(f"{d}/data/nmf_errors.csv")
    r = json.load(open(f"{d}/run_manifest.json"))["nmf_rank"]
    combined[tag] = (e, r)
    ax = axes[i]
    ax.plot(e.k, e.error, "-o", ms=3, lw=1.2, color="#4C72B0")
    ax.axvline(r, ls="--", c="red", lw=1.2)
    ax.annotate(f"rank = {r}", (r, e.error.max()), color="red", fontsize=9,
                ha="left", va="top", xytext=(4, -2), textcoords="offset points")
    ax.set_xlabel("Component number (k)"); ax.set_ylabel("Reconstruction error")
    ax.set_title(f"{tag}  (k_range {int(e.k.min())}..{int(e.k.max())})", fontsize=10)
    ax.grid(alpha=.3)
ax = axes[2]
for tag, (e, r) in combined.items():
    n = (e.error - e.error.min()) / (e.error.max() - e.error.min())
    ax.plot(e.k, n, "-o", ms=3, lw=1.2, label=f"{tag} (rank {r})")
    ax.axvline(r, ls="--", lw=1, alpha=.6,
               color=ax.get_lines()[-1].get_color())
ax.set_xlabel("Component number (k)"); ax.set_ylabel("error, min-max scaled")
ax.set_title("both branches, scaled to a common axis", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.tight_layout()
out = f"{G}/nmf_error_vs_k.png"
fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close("all")
print("wrote", out)
for tag, (e, r) in combined.items():
    print(f"  {tag}: k {int(e.k.min())}..{int(e.k.max())}, rank {r}, "
          f"error {e.error.iloc[0]:.6f} -> {e.error.iloc[-1]:.6f}")
