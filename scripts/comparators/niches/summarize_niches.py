#!/usr/bin/env python
"""Roll the per-core NICHES manifests + quant tables into the summary CSVs used in METHODS.md.

Emits, into <root>/cellchatdb2/:
    summary_per_core.csv        one row per (core, imputation): cells, edges, density, ALRA rank,
                                mechanisms detected, wall time, peak R heap
    summary_per_tier.csv        one row per imputation tier: totals and grade means
    summary_confound.csv        the density ~ core-size ~ grade confound, quantified
    summary_requested_lr.csv    GRN-SORT1 / ANXA1-FPR1 per core and tier

Env: bptf.  Usage:
    python summarize_niches.py --root results/comparators/niches/GBM
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd

REQUESTED = ["GRN—SORT1", "ANXA1—FPR1"]
TIERS = {"noimpute": "none", "alra": "alra"}


def alra_rank(root, core):
    """ALRA prints 'Rank k = N'; RunALRA picks it per object, so it varies core to core."""
    log = os.path.join(root, "cellchatdb2", "logs", f"{core}_alra.log")
    if not os.path.exists(log):
        return np.nan
    m = re.search(r"^Rank k = (\d+)", open(log).read(), re.M)
    return int(m.group(1)) if m else np.nan


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default="results/comparators/niches/GBM")
    a = p.parse_args()
    base = os.path.join(a.root, "cellchatdb2")

    rows, req_rows = [], []
    for tier in TIERS:
        for mf in sorted(glob.glob(os.path.join(base, tier, "core*", "run_manifest.json"))):
            m = json.load(open(mf))
            cdir = os.path.dirname(mf)
            core = f"core{m['tma_id']}"
            summ = pd.read_csv(os.path.join(cdir, "quant", "NeighborhoodToCell_mechanism_summary.csv"))
            n2c, c2cs = m["organizations"]["NeighborhoodToCell"], m["organizations"]["CellToCellSpatial"]
            rows.append(
                dict(
                    core=core, tier=tier, grade=m["grade"], n_cells=m["n_cells"],
                    n_edges=c2cs["n_columns"], edges_per_cell=round(c2cs["n_columns"] / m["n_cells"], 2),
                    n_mechanisms=n2c["n_mechanisms"],
                    n_detected=int((summ.n_nonzero > 0).sum()),
                    n2c_density_pct=round(n2c["density"] * 100, 3),
                    c2cs_density_pct=round(c2cs["density"] * 100, 3),
                    alra_rank_k=alra_rank(a.root, core) if tier == "alra" else np.nan,
                    wall_min=m["wall_minutes"], peak_r_heap_gb=m["r_heap_peak_gb"],
                )
            )
            for lr in REQUESTED:
                hit = summ[summ.mechanism == lr]
                req_rows.append(
                    dict(
                        core=core, tier=tier, grade=m["grade"], lr=lr,
                        present=bool(len(hit)),
                        n_nonzero=int(hit.n_nonzero.iloc[0]) if len(hit) else 0,
                        frac_nonzero=round(float(hit.frac_nonzero.iloc[0]), 4) if len(hit) else 0.0,
                        mean_score=round(float(hit.mean_score.iloc[0]), 4) if len(hit) else 0.0,
                        rank_by_frac=int((summ.frac_nonzero > hit.frac_nonzero.iloc[0]).sum()) + 1
                        if len(hit) else np.nan,
                    )
                )

    per_core = pd.DataFrame(rows).sort_values(["tier", "n_cells"])

    # True peak RSS from /usr/bin/time -l. The `peak_r_heap_gb` in the manifest is R's own
    # gc() "max used" accounting (Ncells*56B + Vcells*8B), which counts allocations that were
    # freed and reused, so it overstates the process high-water mark by ~2x. Prefer the RSS.
    tpath = os.path.join(base, "run_timings.csv")
    if os.path.exists(tpath):
        t = pd.read_csv(tpath)
        t["tier"] = t["impute"]
        t = t[["core", "tier", "wall_sec", "peak_rss_gb"]].drop_duplicates(["core", "tier"], keep="last")
        per_core = per_core.merge(t, on=["core", "tier"], how="left")
    per_core.to_csv(os.path.join(base, "summary_per_core.csv"), index=False)

    per_tier = (
        per_core.groupby("tier")
        .agg(n_cores=("core", "size"), n_cells=("n_cells", "sum"), n_edges=("n_edges", "sum"),
             mean_detected=("n_detected", "mean"), min_detected=("n_detected", "min"),
             max_detected=("n_detected", "max"), mean_density_pct=("n2c_density_pct", "mean"),
             total_wall_min=("wall_min", "sum"),
             max_peak_rss_gb=("peak_rss_gb", "max"), max_r_heap_gb=("peak_r_heap_gb", "max"))
        .round(2).reset_index()
    )
    per_tier.to_csv(os.path.join(base, "summary_per_tier.csv"), index=False)

    conf = []
    for tier, s in per_core.groupby("tier"):
        r = np.corrcoef(np.log10(s.n_cells), s.n2c_density_pct)[0, 1]
        g = s.groupby("grade").agg(n_cores=("core", "size"), mean_cells=("n_cells", "mean"),
                                   mean_density_pct=("n2c_density_pct", "mean"),
                                   mean_detected=("n_detected", "mean"))
        conf.append(dict(
            tier=tier, pearson_r_logcells_vs_density=round(float(r), 3),
            high_n_cores=int(g.loc["high", "n_cores"]), low_n_cores=int(g.loc["low", "n_cores"]),
            high_mean_cells=round(float(g.loc["high", "mean_cells"]), 0),
            low_mean_cells=round(float(g.loc["low", "mean_cells"]), 0),
            high_mean_density_pct=round(float(g.loc["high", "mean_density_pct"]), 2),
            low_mean_density_pct=round(float(g.loc["low", "mean_density_pct"]), 2),
            cell_count_ratio_high_over_low=round(
                float(g.loc["high", "mean_cells"] / g.loc["low", "mean_cells"]), 2),
        ))
    pd.DataFrame(conf).to_csv(os.path.join(base, "summary_confound.csv"), index=False)

    pd.DataFrame(req_rows).sort_values(["lr", "tier", "core"]).to_csv(
        os.path.join(base, "summary_requested_lr.csv"), index=False)

    print(per_core.to_string(index=False))
    print("\n--- per tier ---")
    print(per_tier.to_string(index=False))
    print("\n--- confound ---")
    print(pd.DataFrame(conf).to_string(index=False))
    print("\n--- requested LRs (mean over cores) ---")
    rq = pd.DataFrame(req_rows)
    print(rq.groupby(["lr", "tier"]).agg(
        cores_present=("present", "sum"), mean_frac_nonzero=("frac_nonzero", "mean"),
        mean_rank=("rank_by_frac", "mean")).round(4).to_string())
    print(f"\nwrote 4 summary CSVs to {base}")


if __name__ == "__main__":
    main()
