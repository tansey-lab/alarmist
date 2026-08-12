#!/usr/bin/env python
"""SpatialDM completeness pass — replots from the PERSISTED spatialdm.h5ad (no recomputation).

Fixes three distinct bugs from the first pass, all of which produced BLANK 7,544-byte PNGs
because the code did `fn(); plt.gcf().savefig(...)` for functions that do not leave a drawn
figure current:

  1. `pl.plot_pairs` ends every iteration with `plt.show(); plt.close()`, so `gcf()` afterwards
     is a NEW empty figure. Fix: call `plot_selected_pair` directly (it does not close), or use
     the authors' own `pdf=` argument.
  2. `pl.ligand_ct` / `pl.receptor_ct` are NOT plot functions at all -- they `return` a
     cells x cell-type DataFrame. Fix: save the frame and plot it ourselves.
  3. `pl.chord_celltype` / `pl.chord_LR` render with holoviews+bokeh, not matplotlib, and take
     their own `save=` path. Fix: pass `save=`.

Every write goes through _save(), which REFUSES to write a figure with no drawn content -- so a
blank plot is reported as a failure instead of silently landing on disk.

Usage: python plot_spatialdm_full.py --run-dir DIR [--cores 1,13] [--out-sub plots_full]
"""
import argparse, json, os, sys, time, warnings
# The env's bin/ is NOT on PATH when the interpreter is invoked by absolute path -- which is how
# this repo must run things, because conda activate is broken on this box (CLAUDE.md). bokeh's
# PNG/SVG export shells out to `geckodriver`/`firefox`, which are installed in exactly that bin/,
# so without this every chord diagram dies with "Neither firefox and geckodriver nor a variant of
# chromium browser and chromedriver are available on system PATH".
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spatialdm as sdm
import spatialdm.plottings as pl
warnings.filterwarnings("ignore")
# bokeh/holoviews (used by chord_celltype / chord_LR) writes a default output_file named
# after the CALLING SCRIPT into the script's own directory -- i.e. a 600 KB .html dropped into
# scripts/, violating the code/output separation rule. Pin it to a temp path.
try:
    import tempfile as _tf
    from bokeh.io import output_file as _bokeh_output_file
    _bokeh_output_file(os.path.join(_tf.gettempdir(), "bokeh_scratch.html"))
except Exception:
    pass

p = argparse.ArgumentParser()
p.add_argument("--run-dir", required=True)
p.add_argument("--cores", default="", help="comma list; default = all with a saved h5ad")
p.add_argument("--out-sub", default="plots_full")
p.add_argument("--n-top", type=int, default=3)
p.add_argument("--requested-lrs", default="GRN_SORT1,ANXA1_FPR1")
a = p.parse_args()
t0 = time.time()
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
NBLANK = [0]

def _has_content(fig):
    for ax in fig.get_axes():
        if (ax.collections or ax.images or ax.lines or ax.patches or
                ax.texts or len(ax.get_children()) > 10):
            return True
    return False

def _save(path, fig=None):
    """Write only if the figure actually has drawn content."""
    fig = fig or plt.gcf()
    if not fig.get_axes() or not _has_content(fig):
        plt.close("all"); NBLANK[0] += 1
        raise RuntimeError("figure is BLANK — wrong save mechanism for this function")
    fig.savefig(path, dpi=170, bbox_inches="tight"); plt.close("all")

def guard(fn, nm, d):
    try:
        fn(); log(f"    ok: {nm}")
    except Exception as e:
        plt.close("all"); log(f"    FAILED {nm}: {type(e).__name__}: {e}")

cores = ([c.strip() for c in a.cores.split(",") if c.strip()] or
         sorted([d for d in os.listdir(a.run_dir)
                 if os.path.exists(os.path.join(a.run_dir, d, "data", "spatialdm.h5ad"))],
                key=lambda x: int(x) if x.isdigit() else 0))
requested = [s.strip() for s in a.requested_lrs.split(",") if s.strip()]
log(f"cores with a persisted object: {cores}")

made = 0
for core in cores:
    ddir = os.path.join(a.run_dir, core, "data")
    out = os.path.join(a.run_dir, core, a.out_sub); os.makedirs(out, exist_ok=True)
    log(f"--- core {core} ---")
    sub = sdm.read_spatialdm_h5ad(os.path.join(ddir, "spatialdm.h5ad"))
    gr = sub.uns["global_res"]
    ranked = gr.sort_values("fdr").index.tolist()
    top = [x for x in ranked if gr.loc[x, "selected"]][:a.n_top]
    focus = top + [r for r in requested if r in gr.index]
    xy = np.asarray(sub.obsm["spatial"], float)
    ASPECT = (np.ptp(xy[:, 1]) / np.ptp(xy[:, 0])) if np.ptp(xy[:, 0]) else 1.0

    if "celltypes" not in sub.obsm:
        sub.obsm["celltypes"] = pd.get_dummies(sub.obs["cell_type"]).astype(float)

    # --- v1-column shim -------------------------------------------------------------
    # compute_pathway does  geneInter.groupby('pathway_name').interaction_name  and chord_LR
    # does  geneInter.interaction_name_2.str.split('-')  -- both are CellChatDB *v1* column
    # conventions that our injected v2 frame does not carry (ours has the name in the INDEX).
    # Add them rather than leave two author plots unavailable.
    gi = sub.uns["geneInter"]
    if "interaction_name" not in gi.columns:
        gi["interaction_name"] = gi.index.astype(str)
    if "interaction_name_2" not in gi.columns:
        L = sub.uns["ligand"].apply(lambda r: "_".join(r.dropna().astype(str)), axis=1)
        R = sub.uns["receptor"].apply(lambda r: "_".join(r.dropna().astype(str)), axis=1)
        gi["interaction_name_2"] = (L.reindex(gi.index).astype(str) + "-" +
                                    R.reindex(gi.index).astype(str))
    sub.uns["geneInter"] = gi

    # --- 1. per-pair spatial detail: call plot_selected_pair directly (it does NOT close) ---
    meth = sub.uns["local_stat"]["local_method"]
    key = "local_z_p" if meth == "z-score" else "local_perm_p"
    spots = 1 - sub.uns[key]; sel_ind = sub.uns[key].index
    for lr in focus:
        def _pp(lr=lr):
            pl.plot_selected_pair(sub, lr, spots, sel_ind, figsize=(22, 4.2),
                                  cmap="Greens", cmap_l="Reds", cmap_r="Reds", s=0.6)
            _save(os.path.join(out, f"selected_pair_{lr}.png"))
        guard(_pp, f"selected_pair_{lr}", out); made += 1
    # authors' own multi-pair saver -> real PDF.
    # `plot_selected_pair` indexes `spots`/`n_spots` by the GLOBALLY SELECTED pairs only, so a
    # non-selected requested LR raises KeyError -- and because PdfPages.__exit__ still flushes
    # the pages written so far, the file lands on disk silently truncated at that point (and the
    # os.path.exists check below passes). Filter to what SpatialDM can actually draw.
    pdf_focus = [lr for lr in focus if lr in sel_ind]
    dropped = [lr for lr in focus if lr not in sel_ind]

    def _pdf():
        pl.plot_pairs(sub, pdf_focus, pdf=os.path.join(out, "plot_pairs"), figsize=(22, 4.2), s=0.6)
        plt.close("all")
        if not os.path.exists(os.path.join(out, "plot_pairs.pdf")):
            raise RuntimeError("pdf not written")
    guard(_pdf, f"plot_pairs.pdf ({len(pdf_focus)} pages)", out)
    if dropped:
        log(f"    note: {dropped} not globally selected in this core -> no local statistic "
            f"exists, so SpatialDM cannot plot them (excluded from the pdf, by design)")

    # --- 2. ligand_ct / receptor_ct RETURN DATAFRAMES -> save + plot ourselves ---
    for lr in focus:
        def _ct(lr=lr):
            L = pl.ligand_ct(sub, lr); R = pl.receptor_ct(sub, lr)
            Ls, Rs = np.asarray(L).sum(0), np.asarray(R).sum(0)
            cts = list(sub.obsm["celltypes"].columns)
            df = pd.DataFrame({"cell_type": cts, "ligand_contrib": Ls, "receptor_contrib": Rs})
            df.to_csv(os.path.join(out, f"celltype_contrib_{lr}.csv"), index=False)
            o = np.argsort(-(Ls + Rs)); idx = np.arange(len(cts))
            fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(cts) + 1.8))
            ax.barh(idx - .2, Ls[o], height=.4, label="ligand", color="#4C72B0")
            ax.barh(idx + .2, Rs[o], height=.4, label="receptor", color="#C44E52")
            ax.set_yticks(idx); ax.set_yticklabels([cts[i] for i in o], fontsize=8)
            ax.invert_yaxis(); ax.legend(fontsize=8)
            ax.set_xlabel("summed local contribution"); ax.set_title(f"{lr}: cell-type attribution",
                                                                    fontsize=10)
            _save(os.path.join(out, f"celltype_contrib_{lr}.png"), fig)
        guard(_ct, f"celltype_contrib_{lr}", out); made += 1

    # --- 3. chord diagrams render via holoviews/bokeh -> use their own save= ---
    for lr in focus[:2]:
        guard(lambda lr=lr: pl.chord_celltype(sub, pairs=[lr], min_quantile=0.5, ncol=1,
                                              save=os.path.join(out, f"chord_celltype_{lr}.png")),
              f"chord_celltype_{lr}", out)
    # chord_LR does `adata.obs.loc[:, sender]` -- it wants ONE obs COLUMN PER CELL TYPE, because
    # the intestine tutorial's obs holds deconvolution proportions (`A1.obsm['celltypes'] =
    # A1.obs[A1.obs.columns[:-1]]`). Passing category names raises KeyError. Single-cell labels
    # give the same thing as a one-hot, which is what obsm['celltypes'] already is -- copy it out.
    # Note senders/receivers are ZIPPED, not crossed: this plots the self-self diagonal, one
    # panel per cell type.
    cts = list(sub.obsm["celltypes"].columns)
    for ct in cts:
        sub.obs[ct] = np.asarray(sub.obsm["celltypes"][ct], dtype=float)
    guard(lambda: pl.chord_LR(sub, senders=cts, receivers=cts, min_quantile=0.5,
                              save=os.path.join(out, "chord_LR_all.png")), "chord_LR_all", out)

    # --- 4. global volcano + pathway enrichment (matplotlib; these do work with gcf) ---
    def _gp():
        pl.global_plot(sub, pairs=top[:4], figsize=(6, 5), cmap="RdGy_r", vmin=-1.5, vmax=2)
        _save(os.path.join(out, "global_plot.png"))
    guard(_gp, "global_plot", out); made += 1
    def _dp():
        sel = [x for x in ranked if gr.loc[x, "selected"]]
        res = sdm.compute_pathway(sub, dic={f"core{core}": pd.Index(sel)})
        res.to_csv(os.path.join(out, "pathway_enrichment.csv"))
        pl.dot_path(res, figsize=(6, 10)); _save(os.path.join(out, "dot_path.png"))
    guard(_dp, "dot_path", out); made += 1

n_png = sum(len([f for f in fs if f.endswith(".png")]) for _, _, fs in os.walk(a.run_dir))
blanks = 0
for root, _, fs in os.walk(a.run_dir):
    for f in fs:
        if f.endswith(".png") and os.path.getsize(os.path.join(root, f)) == 7544: blanks += 1
json.dump({"script": "plot_spatialdm_full.py", "cores": cores, "n_png_total": n_png,
           "blank_7544_remaining": blanks, "blanks_refused": NBLANK[0],
           "wall_min": round((time.time() - t0) / 60, 1)},
          open(os.path.join(a.run_dir, "plot_full_manifest.json"), "w"), indent=2)
log(f"DONE in {(time.time()-t0)/60:.1f} min | refused {NBLANK[0]} blank figures | "
    f"{blanks} legacy 7544-byte blanks still on disk")
