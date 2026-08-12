#!/usr/bin/env python
"""Assemble the LIANA+ GBM figure report: COPY existing PNGs, draw nothing.

Scope, as requested:
  * `cellchatdb2` tier only. (That is in fact every LIANA+ GBM run -- the `default`
    tier was never produced -- but the restriction is stated so the report cannot be
    mistaken for a two-tier comparison.)
  * ONLY figures the authors' own tutorials demonstrate, plus the two standing
    required interactions (GRN^SORT1, ANXA1^FPR1).

Covers FOUR of the authors' notebooks, all run on the cellchatdb2 resource:
  bivariate.ipynb | inflow_score.ipynb | inflow_mofaflex.ipynb | LRIC_tutorial.ipynb | misty.ipynb

What is therefore EXCLUDED, and why (all of it is real output, none of it is deleted --
it stays under results/comparators/liana/GBM/):

  default/, default_inflow/, nmf_*_default/   the `default` tier (LIANA's own `consensus`
                                    resource). Real, complete, and 115 figures -- but the
                                    request for this report was cellchatdb2 ONLY. Its result is
                                    summarised in the README instead: per-interaction statistics
                                    are bit-identical across the two resources.
  cellchatdb2_morans/               a metric VARIANT (local_name='morans'), not the tutorial's
                                    default cosine. Same 131 pairs, same calls. Summarised in
                                    the README; figures left in place.
  factor_annotation/                our PROGENy / CellChatDB-pathway annotation of the NMF
                                    factors. No tutorial analogue.
  nmf_inflow/punch_level/           our 7-vs-6 punch-level grade test.
  mofaflex_inflow/**/factor_by_punch_grade.png   same -- our punch-level test, not tutorial.
                                    (The tutorial's own group summary, cell 63, IS included as
                                    focus_factors_by_region.png.)
  lric_percore/figures/{grade_comparison_*, support_*_per_punch, target_availability,
                        direction_asymmetry_paired}   our cross-punch aggregation and grade
                                    analysis, not outputs of the tutorial's own calls.
  _benchmarks/                      timing probes, not results.

  nmf_*/plots_full/factors/*        our own construction. bivariate.ipynb's NMF section
                                    (cells 50-57) draws exactly two figures -- the elbow
                                    and the per-factor spatial maps. Cell 57 is a
                                    `.head(10)` TABLE, not a plot. The identity panels,
                                    the three factor x {sender, receiver, LR} heatmaps,
                                    celltype_communication_by_factor and top_lri_dot_by_factor
                                    are ALARMIST-style summaries we added; two of them are
                                    marked "OUR plots -- LIANA demonstrates neither" in
                                    plot_liana_full.py.
  nmf_*/plots_full/{global,interactions}/  tutorial FIGURE TYPES, but drawn a second time on
                                    a different LR selection (top-N by total inflow rather
                                    than the tutorial's own criterion). The canonical
                                    reproduction of the same tutorial cells already lives in
                                    cellchatdb2*/plots/. Kept out to avoid two versions of
                                    the same tutorial figure in one report.
  nmf_inflow/punch_level/*          entirely our own (the 7-vs-6 punch-level grade test).
  GBM/nmf_error_vs_k.png            our cross-branch elbow comparison; no tutorial analogue.

Usage:
  python build_report_figures.py --results-dir results/comparators/liana/GBM \
                                 --out-dir reports/liana_plus_GBM_cellchatdb2
"""
import argparse, hashlib, json, os, shutil, subprocess, time

p = argparse.ArgumentParser()
p.add_argument("--results-dir", default="results/comparators/liana/GBM")
p.add_argument("--out-dir", default="reports/liana_plus_GBM_cellchatdb2")
p.add_argument("--required-lrs", default="GRN-SORT1,ANXA1-FPR1",
               help="filename tags (^ rendered as -) of the standing required interactions")
p.add_argument("--dry-run", action="store_true")
a = p.parse_args()

R, OUT = a.results_dir, a.out_dir
REQ = [s.strip() for s in a.required_lrs.split(",") if s.strip()]
log = lambda *m: print(*m, flush=True)

# ---------------------------------------------------------------- the whitelist
# (source subdir, filename predicate, destination subdir, tutorial provenance)
# Predicates take the basename. `None` destination means "route by required-LR match".
BIV = "cellchatdb2/plots"
INF = "cellchatdb2_inflow/plots"
MFX = "mofaflex_inflow/plots"
MFXS = "mofaflex_inflow/sensitivity_nzf0.001/plots"
LRIC = "lric_percore/figures"
MSTY = "misty/linear_fullslide/plots"

# MISTy names its per-target figures after the RECEPTOR alone (importances_SORT1.png), because
# its unit is the target, not the LR pair. So the required-LR test has to accept both the pair
# tag and the bare receptor token, or the two motif-1 figures would land in the top-hits dir.
REQ_TOKENS = REQ + [t.split("-")[-1] for t in REQ]


def is_req(name):
    return any(t in name for t in REQ_TOKENS)

PLAN = [
    # ---- bivariate branch: bivariate.ipynb -------------------------------------
    (BIV, lambda n: n.startswith("bandwidth_query"),      "bivariate/00_calibration", "cell 19  li.ut.query_bandwidth"),
    (BIV, lambda n: n.startswith("connectivity_idx"),     "bivariate/00_calibration", "cell 23  li.pl.connectivity"),
    (BIV, lambda n: n.startswith("top_"),                 "bivariate/01_overview",    "cells 32-36  sc.pl.spatial, top pairs by Moran's R"),
    (BIV, lambda n: n.startswith("local_") and n.endswith("_pvals.png"), "bivariate/02_top6_lri", "cell 42  permutation p-value map"),
    (BIV, lambda n: n.startswith("local_") and n.endswith("_cats.png"),  "bivariate/02_top6_lri", "cell 44  local category map (high-high / low-low)"),
    (BIV, lambda n: n.startswith("local_"),               "bivariate/02_top6_lri",    "cells 32-40  local bivariate score map"),
    (BIV, lambda n: n.startswith("genes_"),               "bivariate/02_top6_lri",    "cell 37  constituent ligand / receptor expression"),
    ("cellchatdb2/plots/requested", lambda n: True,       "bivariate/03_requested",   "same calls, forced onto the two required LRs"),
    ("nmf_bivariate/plots", lambda n: n == "elbow.png",       "bivariate/04_nmf", "cell 50  li.multi.nmf rank selection"),
    ("nmf_bivariate/plots", lambda n: n == "factor_maps.png", "bivariate/04_nmf", "cell 55  sc.pl.spatial over NMF_W"),

    # ---- inflow branch: inflow_score.ipynb -------------------------------------
    (INF + "/global", lambda n: n.startswith("bandwidth_query"),  "inflow/00_calibration", "cells 19/20  li.ut.query_bandwidth"),
    (INF + "/global", lambda n: n.startswith("connectivity_idx"), "inflow/00_calibration", "cell 24  li.pl.connectivity"),
    (INF + "/global", lambda n: n.startswith("spatial_"),         "inflow/01_global",      "cells 11/43  sc.pl.embedding(cell_type / region)"),
    (INF + "/global", lambda n: n == "dotplot_global.png",        "inflow/01_global",      "cell 59  li.pl.dotplot(uns_key='global_interactions')"),
    (INF + "/global", lambda n: n == "pair_proximity.png",        "inflow/01_global",      "cells 75-78  li.ut.spatial_pair_proximity"),
    (INF + "/global", lambda n: n == "rank_aggregate_dotplot.png","inflow/01_global",      "cells 83-88  li.mt.rank_aggregate (spatially constrained)"),
    (INF + "/global", lambda n: n.startswith("rank_aggregate_"),  None,                    "cells 83-88  rank_aggregate, per interaction"),
    (INF + "/interactions", lambda n: True,                       None,                    "cells 45-71  per-interaction suite"),
    ("nmf_inflow/plots", lambda n: n == "elbow.png",              "inflow/04_nmf",         "cell 50 of bivariate.ipynb, applied to inflow output"),
    ("nmf_inflow/plots", lambda n: n == "factor_maps.png",        "inflow/04_nmf",         "cell 55 of bivariate.ipynb, applied to inflow output"),

    # ---- MOFA-Flex branch: inflow_mofaflex.ipynb -------------------------------
    # The authors' PRESCRIBED unsupervised route for single-cell spatial. Both fits are kept:
    # the tutorial-faithful one (nzf>0.01), in which the authors' own QC deletes both arms of
    # ALARMIST motif 1, and the sensitivity fit (nzf>0.001) in which they survive.
    (MFX,  lambda n: n == "qc_inflow_distributions.png",   "mofaflex/00_qc",        "cell 18  inflow score + nonzero_fraction QC"),
    (MFX,  lambda n: n == "data_overview.png",             "mofaflex/00_qc",        "cell 25  MOFAFLEX.fit(plot_data_overview=True)"),
    (MFX,  lambda n: n == "factor_correlation.png",        "mofaflex/01_factors",   "cell 31  mfl.pl.factor_correlation"),
    (MFX,  lambda n: n.startswith("variance_explained"),   "mofaflex/01_factors",   "cells 33/35  mfl.pl.variance_explained (overall and by view = by sender)"),
    (MFX,  lambda n: n == "top_weights.png",               "mofaflex/01_factors",   "cell 43  mfl.pl.top_weights"),
    (MFX,  lambda n: n.startswith("umap_"),                "mofaflex/02_embedding", "cells 48/50  UMAP over factor space, leiden + focus factors"),
    (MFX,  lambda n: n.startswith("spatial_"),             "mofaflex/02_embedding", "cells 52/61  spatial maps, leiden + focus factors"),
    (MFX,  lambda n: n == "focus_factors_by_region.png",   "mofaflex/02_embedding", "cell 63  group summary (grade, in place of major_brain_region)"),
    (MFX,  lambda n: n == "dotplot_focus_factors.png",     "mofaflex/03_liana_bridge", "cell 57  li.pl.dotplot over MOFA-Flex weights"),
    (MFX,  lambda n: n.startswith("circle_plot_"),         "mofaflex/03_liana_bridge", "cell 59  li.pl.circle_plot per focus factor"),
    (MFXS, lambda n: n != "factor_by_punch_grade.png",     "mofaflex/04_sensitivity_nzf0.001", "identical calls at nzf>0.001, where both motif-1 arms survive QC"),

    # ---- LRIC / cross-PCF: LRIC_tutorial.ipynb ---------------------------------
    # Run PER PUNCH: both functions normalise by density against the bounding box, and the
    # TMA's global bbox is 42.2% occupied, so a whole-slide run inflates every g(r) ~2.4x.
    (LRIC, lambda n: n == "gr_curves_per_punch.png",       "lric", "cells 12-30  li.mt.lric g(r), cell-type-informed, per punch"),
    (LRIC, lambda n: n == "gr_curves_agnostic.png",        "lric", "cells 12-30  li.mt.lric g(r), cell-type-agnostic"),
    (LRIC, lambda n: n == "cross_pcf_curves_per_punch.png","lric", "cells 12-20  li.mt.cross_pcf g(r) baseline"),
    (LRIC, lambda n: n == "ratio_vs_distance.png",         "lric", "LRIC / cross-PCF per distance bin -- does expression add over co-location?"),
    (LRIC, lambda n: n == "whole_slide_vs_perpunch.png",   "lric", "control: the density artefact that forces the per-punch design"),

    # ---- LR-MISTy: misty.ipynb -------------------------------------------------
    (MSTY, lambda n: n.startswith("importances_"),         None,   "li.pl.interactions, forced onto the two required receptors"),
    (MSTY, lambda n: n.startswith("target_metrics_"),      "misty/01_metrics", "li.pl.target_metrics -- gain_R2 and multi_R2 by target"),
    (MSTY, lambda n: n == "gain_R2_hist.png",              "misty/01_metrics", "distribution of gain_R2 over all 382 receptor targets"),
    (MSTY, lambda n: n == "contributions_top.png",         "misty/01_metrics", "li.pl.contributions -- per-view contribution"),
    (MSTY, lambda n: n == "interactions_extra_top200.png", "misty/01_metrics", "li.pl.interactions -- ligand predictor importances, extra view"),
]

# destination when the rule routes by required-LR membership
ROUTE = {INF + "/global": ("inflow/01_global", "inflow/03_requested"),
         INF + "/interactions": ("inflow/02_top6_lri", "inflow/03_requested"),
         MSTY: ("misty/01_metrics", "misty/03_requested")}

# ------------------------------------------------------------------------ copy
copied, seen, manifest = 0, set(), []
for src_sub, pred, dst_sub, prov in PLAN:
    src_dir = os.path.join(R, src_sub)
    if not os.path.isdir(src_dir):
        log(f"  MISSING source dir, skipped: {src_dir}"); continue
    for name in sorted(os.listdir(src_dir)):
        if not name.endswith(".png"): continue
        key = (src_sub, name)
        if key in seen or not pred(name): continue          # first matching rule wins
        seen.add(key)
        dst_rel = dst_sub if dst_sub else ROUTE[src_sub][1 if is_req(name) else 0]
        dst_dir = os.path.join(OUT, dst_rel)
        src, dst = os.path.join(src_dir, name), os.path.join(dst_dir, name)
        if not a.dry_run:
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copy2(src, dst)                          # copy2 preserves mtime
        h = (hashlib.md5(open(src, "rb").read()).hexdigest()[:8] if not a.dry_run else "")
        manifest.append({"figure": f"{dst_rel}/{name}",
                         "source": os.path.relpath(src),
                         "tutorial": prov, "md5_8": h,
                         "required_lr": bool(is_req(name))})
        copied += 1

n_req = sum(1 for m in manifest if m["required_lr"])
log(f"copied {copied} PNG ({n_req} of them required-LR figures) -> {OUT}")

if not a.dry_run:
    try: sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except Exception: sha = None
    json.dump({"script": "build_report_figures.py", "results_dir": R, "out_dir": OUT,
               "tier": "cellchatdb2", "scope": "authors' tutorial figures + required LRs",
               "required_lrs": REQ, "n_png": copied, "n_required_lr_png": n_req,
               "git_sha": sha, "built": time.strftime("%Y-%m-%d %H:%M:%S"),
               "figures": manifest},
              open(os.path.join(OUT, "figure_manifest.json"), "w"), indent=2)
    log(f"wrote {OUT}/figure_manifest.json")
