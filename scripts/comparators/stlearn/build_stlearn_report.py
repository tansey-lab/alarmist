#!/usr/bin/env python
"""Build the self-contained HTML report for a stLearn comparator run.

The HTML is gitignored (it lives under reports/); THIS SCRIPT is the tracked artifact.
Rerun it to regenerate the report from the run directory.

  python build_stlearn_report.py \
      --run-dir results/comparators/stlearn/GBM/cellchatdb2 \
      --fig-dir reports/stlearn_GBM_cellchatdb2/figures \
      --out     reports/stlearn_GBM_cellchatdb2/stlearn_GBM_cellchatdb2_report.html

No external libraries, no CDN: vanilla JS + inline SVG, figures base64-embedded.
"""
import argparse, base64, html, json, os, textwrap, time

import numpy as np
import pandas as pd

P = argparse.ArgumentParser()
P.add_argument("--run-dir", required=True)
P.add_argument("--fig-dir", required=True)
P.add_argument("--out", required=True)
P.add_argument("--label", default="cell_type")
P.add_argument("--inline-max-mb", type=float, default=3.5,
               help="figures larger than this are linked file:// instead of inlined (skill rule)")
A = P.parse_args()

DATA = os.path.join(A.run_dir, "data")
log = lambda *m: print(f"[{time.strftime('%H:%M:%S')}]", *m, flush=True)
E = html.escape

# ------------------------------------------------------------------ load run facts
run_man = json.load(open(os.path.join(A.run_dir, "run_manifest.json")))
plot_man_p = os.path.join(A.run_dir, "plots_tutorial", "plot_manifest.json")
plot_man = json.load(open(plot_man_p)) if os.path.exists(plot_man_p) else {}
summ = pd.read_csv(os.path.join(DATA, "lr_summary.csv"), index_col=0)
spot_meta = pd.read_csv(os.path.join(DATA, "spot_meta.csv"))
pooled = pd.read_csv(os.path.join(DATA, f"lr_cci_{A.label}.csv"), index_col=0)
CTS = list(pooled.index)
LRS = list(summ.index)
log(f"{len(LRS)} LR pairs, {len(CTS)} cell types")

dom = spot_meta["dominant_cell_type"].value_counts().reindex(CTS).fillna(0).astype(int)
sender_tot = pooled.sum(axis=1).reindex(CTS)

# per-LR cell-type matrices -> compact parallel arrays (row-major, sender x receiver)
cnt_dir = os.path.join(DATA, f"per_lr_cci_{A.label}")
pv_dir = os.path.join(DATA, f"per_lr_cci_pvals_{A.label}")
counts, pvals = [], []
for lr in LRS:
    c = pd.read_csv(os.path.join(cnt_dir, f"{lr}.csv"), index_col=0).reindex(index=CTS, columns=CTS)
    counts.append([int(x) for x in np.nan_to_num(c.values).ravel()])
    fp = os.path.join(pv_dir, f"{lr}.csv")
    if os.path.exists(fp):
        p = pd.read_csv(fp, index_col=0).reindex(index=CTS, columns=CTS)
        pvals.append([round(float(x), 3) for x in np.nan_to_num(p.values, nan=1.0).ravel()])
    else:
        pvals.append([1.0] * len(CTS) ** 2)
log("per-LR matrices packed")

PAYLOAD = {
    "cts": CTS,
    "lrs": LRS,
    "nSpots": [int(x) for x in summ["n_spots"]],
    "nSig": [int(x) for x in summ["n_spots_sig"]],
    "nSigP": [int(x) for x in summ["n_spots_sig_pval"]],
    "nCci": [int(x) for x in summ["n_cci_sig_cell_type"].fillna(0)],
    "counts": counts,
    "pvals": pvals,
    "requested": plot_man.get("requested_lrs", ["GRN_SORT1", "ANXA1_FPR1"]),
}

# ------------------------------------------------------------------ figures -> base64
def img_src(path):
    mb = os.path.getsize(path) / 1048576
    if mb > A.inline_max_mb:
        return "file://" + os.path.abspath(path), mb, False
    b64 = base64.b64encode(open(path, "rb").read()).decode()
    return f"data:image/png;base64,{b64}", mb, True

# --- Which figures earn a place in the report -------------------------------------------
# The full tutorial set (29 figures) stays on disk in plots_tutorial/ and in figures/. The
# REPORT carries only those that support a biological claim. Excluded: pipeline/QC checks and
# plain "where does cell type X sit" descriptives, which say nothing the tables do not.
KEEP_FIGS = [
    "cell42_composition_mGAM.png",
    "cell79_cci_check.png",
    "cell58_lr_summary_n50.png",
    "cell67_lr_result_C3_C3AR1.png",
    "cell81_ccinet_all.png",
    "cell81_ccinet_C3_C3AR1.png",
    "requested/cell67_lr_result_rank21_GRN_SORT1.png",
    "requested/cell81_ccinet_rank21_GRN_SORT1.png",
    "requested/cell67_lr_result_rank99_ANXA1_FPR1.png",
    "requested/cell81_ccinet_rank99_ANXA1_FPR1.png",
]
# (figure, why it is not in the report)
DROPPED_FIGS = [
    ("cell23_cell_labels.png", "descriptive", "Single-cell cell-type map. Shows where cell types sit; makes no claim."),
    ("cell40_grid_vs_cell_labels.png", "pipeline check", "Grid vs single-cell labels — a gridding sanity check. Its one substantive point (abundant types absorb the spots rarer types contribute to) is carried quantitatively by the footprint table in §10.1."),
    ("cell42_composition_*.png (8 of 9)", "descriptive", "Per-cell-type spatial composition for AC-like, Glial-Neuronal, Lymphoid, MES-like, NPC-like, OPC-like, Vascular, non-mGAM. Location description. Only mGAM is kept, because its everywhere-but-dominant-nowhere pattern is load-bearing for §10.1."),
    ("cell44_gene_C3.png, cell44_gene_C3AR1.png", "descriptive", "Expression maps for the top pair's ligand and receptor. Gene-location description; the interaction itself is shown in §10.2."),
    ("cell58_lr_summary_n500.png", "redundant", "All 526 pairs as one bar chart — same ranking as the top-50 panel, unreadable at that width. The full ranking is browsable in §11."),
    ("cell81_ccinet_CNTN2_CNTN2.png", "redundant", "Rank-2 network. CNTN2–CNTN2 is homophilic, and the caution about homophilic pairs (§10.3) is a statement about the ranking, not something the network adds to."),
    ("cell83_chord_*.png (3)", "redundant", "Chord diagrams for the pooled set and the top 2 pairs — the same edge data as the ccinet networks, in a form that is harder to read at 9 nodes."),
    ("requested/cell83_chord_rank*.png (2)", "redundant", "Chord versions of the two requested interactions; the networks in §10.4 carry the same edges."),
]

FIGS, inlined, linked, total_mb, skipped = {}, 0, 0, 0.0, 0
for root, _, files in os.walk(A.fig_dir):
    for f in sorted(files):
        if not f.endswith(".png"):
            continue
        rel = os.path.relpath(os.path.join(root, f), A.fig_dir)
        if rel not in KEEP_FIGS:
            skipped += 1
            continue
        src, mb, was_inlined = img_src(os.path.join(root, f))
        FIGS[rel] = src
        total_mb += mb
        inlined += was_inlined
        linked += (not was_inlined)
missing = [k for k in KEEP_FIGS if k not in FIGS]
if missing:
    log(f"WARNING: kept-list figures not found on disk: {missing}")
log(f"figures: {inlined} inlined, {linked} linked, {skipped} left out of the report, "
    f"{total_mb:.1f} MB source")

def fig(rel, caption, cls="fig"):
    if rel not in FIGS:
        return f'<p class="warn">MISSING FIGURE: {E(rel)}</p>'
    return (f'<figure class="{cls}"><img loading="lazy" src="{FIGS[rel]}" alt="{E(caption)}">'
            f'<figcaption>{caption}</figcaption></figure>')

def code(s):
    return f'<pre class="code">{E(textwrap.dedent(s).strip())}</pre>'

# ------------------------------------------------------------------ content tables
# our extraction index -> the tutorial's own In[k] prompt
CELL_MAP = [
    (23, 10, "st.pl.cluster_plot(adata, ...)", "single-cell labels"),
    (40, 15, "2-panel cluster_plot(grid) | cluster_plot(adata)", "grid vs single-cell"),
    (42, 16, "3-panel feat_plot | cluster_plot(grid) | cluster_plot(adata)", "per-cell-type composition"),
    (44, 17, "2-panel gene_plot(grid) | gene_plot(adata)", "gene expression"),
    (58, 21, "st.pl.lr_summary(n_top=500 / 50)", "LR ranking"),
    (62, 22, "st.tl.cci.adj_pvals(...)", "significance thresholds"),
    (67, 24, "3-stat st.pl.lr_result_plot", "per-LR spatial map"),
    (79, 27, "st.pl.cci_check(...)", "abundance diagnostic"),
    (81, 28, "st.pl.ccinet_plot(...)", "interaction networks"),
    (83, 29, "st.pl.lr_chord_plot(...)", "chord diagrams"),
]

# Every deviation: (item, tutorial, ours, category, reason)
DEVIATIONS = [
    ("Cell labels", "<code>leiden</code> clustering at resolution 1.05 (In[11])",
     "<code>obs['cell_type']</code>, 9 annotated types; Leiden skipped",
     "better-input",
     "The tutorial runs PCA + neighbours + Leiden purely because its public demo ships no "
     "annotation. <code>use_label</code> only needs <em>a</em> per-cell categorical. We have "
     "expert annotations, so running Leiden would substitute a worse label for a better one. "
     "Nothing downstream changes shape: the same column feeds <code>grid()</code> and "
     "<code>run_cci()</code>."),
    ("Expression matrix", "Xenium raw counts read from <code>cell_feature_matrix.h5</code>",
     "<code>layers['counts']</code> copied into <code>X</code>",
     "forced",
     "Our <code>X</code> is log-normalised. The tutorial is explicit that no log1p may be "
     "applied before <code>run()</code>, because the permutation null selects background genes "
     "of <em>similar expression level</em> to the real ligand and receptor. Log-transforming "
     "compresses genes toward each other and silently invalidates that matching. Using "
     "<code>X</code> would have produced p-values against the wrong null."),
    ("Gene names", "313-gene panel, no underscores present",
     "21 <code>Intergenic_Region_*</code> probes dropped before <code>normalize_total</code>",
     "forced",
     "<code>st.tl.cci.run()</code> hard-errors on any gene name containing <code>_</code>, "
     "because <code>_</code> is its ligand/receptor separator in the "
     "<code>\"LIGAND_RECEPTOR\"</code> encoding. Every offender on our panel is a Xenium genomic "
     "<em>control probe</em> — not a gene, absent from every LR database, and normally excluded "
     "upstream (10x keeps controls in a separate feature type, which is why the tutorial's input "
     "has none). Dropped <em>before</em> normalisation so library size reflects the real panel, "
     "matching the tutorial. 5,119 → 5,098 genes; zero biology lost."),
    ("Object construction", "<code>st.read_xenium(...)</code> from a 10x output bundle",
     "plain h5ad + hand-built <code>uns['spatial']</code>",
     "forced",
     "Both <code>st.convert_scanpy()</code> and <code>st.tl.cci.grid()</code> index "
     "<code>adata.uns['spatial']</code> (<code>grid()</code> copies it verbatim), so a plain "
     "h5ad raises <code>KeyError: 'spatial'</code>. We replicate exactly what "
     "<code>read_xenium</code> builds. Two sub-points: (a) <code>read_xenium</code> itself "
     "creates a <em>blank</em> placeholder image when no image file is given, so a placeholder "
     "is the reader's own behaviour — we just make it 1×1 instead of (1.1·max_coord)² RGBA ≈ "
     "1.34 GB, since every plot passes <code>show_image=False</code>; (b) "
     "<code>scalef=1</code> matches the tutorial's <code>read_xenium(scale=1)</code> for micron "
     "coordinates, making <code>convert_scanpy</code>'s <code>obsm['spatial'] * scale</code> the "
     "identity — asserted at runtime. stLearn's own source notes the scale factor scales the "
     "<em>image</em> to the spots, never the spots, so it cannot displace a data point."),
    ("Grid resolution", "<code>n_row = n_col = 125</code> → 60.2 × 43.8 µm on a 7,521 × 5,471 µm section",
     "<code>n_row=321, n_col=146</code> → <strong>51.3 × 51.3 µm</strong>",
     "forced",
     "The tutorial's 125 is an author-declared compute/resolution knob (its own text: “The higher "
     "resolution, the better this represents the single cell data but the longer the "
     "computation takes”), and its resulting spots are <em>rectangular</em> and specific to its "
     "extent. Our TMA is 7,484 × 16,483 µm (2.20:1), so copying 125×125 would give 59.9 × 131.9 µm "
     "— a 2.2:1 elongated spot, i.e. a different neighbourhood in x than in y. We instead preserve "
     "the tutorial's <em>spot area</em> (2,637 µm² → 51.4 µm equivalent square) and let the bin "
     "counts follow. <strong>The resulting ~51 µm is not chosen to match ALARMIST's 50 µm patch "
     "— that is a coincidence of the arithmetic, and no comparator parameter was tuned toward "
     "ALARMIST.</strong>"),
    ("Cells per spot", "11.4 per occupied spot",
     "<strong>7.98</strong> per occupied spot (median 5, max 46)",
     "consequence",
     "A consequence, not a choice. GBM in-core density is lower than the tutorial's breast "
     "section. Matching 11.4 would have required larger spots — i.e. changing the method's "
     "spatial resolution to compensate for sparser tissue. We preserved resolution and let "
     "cells/spot fall where it falls."),
    ("LR database", "<code>connectomeDB2020_lit</code> via <code>load_lrs</code>, 2,293 pairs",
     "CellChatDB v2 → <strong>1,371</strong> single-gene pairs (527 on panel, 526 tested)",
     "protocol",
     "The <code>cellchatdb2</code> tier exists to remove the LR resource as a confounder against "
     "ALARMIST. stLearn's format is one underscore-delimited string, so <code>_</code> cannot "
     "also mean “subunit of” and <strong>1,859 of 3,233 CellChatDB rows (57.5%) are "
     "unrepresentable</strong>. We drop them rather than expand combinatorially: "
     "<code>TGFB1</code> + <code>TGFBR1_TGFBR2</code> → {TGFB1_TGFBR1, TGFB1_TGFBR2} would "
     "assert two interactions CellChat does not claim exist independently, and would "
     "double-count one. This is a <em>limitation of the method</em>, recorded as such."),
    ("<code>n_pairs</code>", "1,000",
     "<strong>10,000</strong>",
     "recommended",
     "Following the tutorial's own inline comment rather than its example value: "
     "<code># Number of random pairs to generate; low as example, recommend ~10,000</code>. "
     "Cost scales with unique LR <em>genes</em>, not pairs × permutations (backgrounds are "
     "cached per gene), so 100× the pairs cost only ~7.6× the time — there is no reason to run "
     "below the authors' recommendation."),
    ("<code>n_perms</code>", "100",
     "<strong>1,000</strong>",
     "recommended",
     "Same pattern: <code># Permutations of cell information to get background, recommend ~1000</code>."),
    ("<code>distance</code>", "250", "250 — unchanged",
     "match",
     "Physical µm via a cKDTree radius, and independent of grid resolution. Kept at the authors' "
     "default per the benchmark rule that each method keeps its own neighbourhood definition and "
     "nothing is harmonised toward ALARMIST's patch size."),
    ("Whole slide vs per core", "one contiguous section",
     "whole-slide, one joint run over 13 disconnected cores",
     "measured",
     "<code>distance=250</code> could in principle link spots across cores from different "
     "patients. Measured rather than assumed: <strong>0 cross-core neighbour pairs at 200 µm; "
     "39 pairs / 43 cells (0.04%) at 250 µm</strong> — bounding boxes come within 150 µm but the "
     "cores are round, so cells never actually get that close. Running per core would also break "
     "the smallest core (819 cells → ~15 occupied spots, below <code>min_spots=20</code>). "
     "Documented rather than worked around."),
    ("<code>run_cci</code> dtype", "works as published",
     "<code>obs['cell_type']</code> coerced to object dtype around the call, categorical restored after",
     "version-gap",
     "A genuine version gap, not a parameter choice. <code>run_cci</code> does "
     "<code>adata.obs[label].values.astype(str)</code> and feeds the result to an "
     "<code>@njit</code> kernel. Under pandas 3.0.5 / numpy 2.4.6 that returns a "
     "StringDtype/object array which numba 0.66 cannot type "
     "(<code>TypingError: non-precise type array(pyobject, 1d, C)</code>), hard-blocking the "
     "entire cell-type half of the workflow. Object-dtype Python strings make the <em>same "
     "expression</em> return a plain <code>&lt;U</code> array, which numba accepts. Labels are "
     "byte-identical — only the dtype changes, asserted at runtime. The alternative was pinning "
     "<code>pandas&lt;3</code>, which we did not do because the coercion is provably "
     "label-preserving."),
    ("Cell-type ↔ column binding", "not an issue for the tutorial's numeric Leiden labels",
     "runtime guard added; aborts if any label mis-binds",
     "safety",
     "<code>get_data_for_counting</code> binds cell types to deconvolution columns by "
     "<strong>substring, first hit wins</strong>: <code>[ct in col for col in cols]</code>. Our "
     "vocabulary contains both <code>mGAM</code> and <code>non-mGAM</code>, and "
     "<code>'mGAM' in 'non-mGAM'</code> is <code>True</code> — so the binding is order-dependent "
     "and a silent mis-assignment would corrupt every interaction reported for that type. "
     "Verified correct here (alphabetical ordering puts <code>mGAM</code> first, so it wins its "
     "own column), but <em>verified rather than assumed</em>. This is a latent stLearn bug that "
     "would bite any dataset with nested cell-type names."),
    ("<code>st.tl.cci.adj_pvals</code>", "called at In[22] with <code>correct_axis='spot', "
     "pval_adj_cutoff=0.05, adj_method='fdr_bh'</code>",
     "called in the plotting pass with a runtime assertion; <strong>proven no-op</strong>",
     "no-op",
     "<code>run()</code> already applies these exact settings internally "
     "(<code>run()</code> defaults <code>adj_method='fdr_bh'</code>, "
     "<code>pval_adj_cutoff=0.05</code>; <code>permutation.py</code> does the per-spot MHT "
     "correction), and <code>adj_pvals</code>' own docstring says so: <em>“Default settings of "
     "this function are already run in st.tl.cci.run”</em>. Verified <em>empirically</em>, not "
     "assumed — see §6. The on-disk results, produced without the call, are the tutorial's "
     "results."),
    ("<code>feat_plot</code> groups", "3 clusters selected by <em>position</em> — "
     "<code>groups[6]</code>, <code>groups[10]</code>, <code>groups[11]</code>",
     "all 9 cell types",
     "complete-iteration",
     "The <code>idc</code>/<code>dcis</code>/<code>stroma</code> names in the tutorial are local "
     "variable names, not a selection rule — the actual selection is positional indexing into "
     "<code>cat.categories</code>. There is no positional equivalent for an annotated 9-type "
     "vocabulary, and any fixed slice would be arbitrary. Same call, same arguments, complete "
     "iteration: this adds no analysis, no parameter and no number, and it removes the "
     "possibility of a favourable subset. The tutorial itself leaves ≥9 of its ≥12 clusters "
     "unplotted."),
    ("<code>feat_plot</code> arguments", "<code>vmax=1, show_color_bar=False</code>", "same",
     "match",
     "Restored after an earlier pass got this wrong (<code>show_color_bar=True</code>, no "
     "<code>vmax</code>), which let each cell type's colour scale stretch to its own maximum so "
     "the proportion panels were <em>not comparable to one another</em>. <code>vmax=1</code> is "
     "what makes “proportion, max = 1” mean the same thing in every panel."),
    ("<code>gene_plot</code> gene", "<code>CXCL12</code>, a prior biological pick for its breast section",
     "ligand + receptor of the top-ranked pair (<code>C3</code>, <code>C3AR1</code>)",
     "stated-rule",
     "The tutorial gives no selection rule, just a gene it cared about. Rather than invent a "
     "biological pick for glioma we apply the tutorial's own <code>lr_summary.index[0]</code> "
     "idiom to gene identity, and state the rule explicitly so it can be challenged."),
    ("<code>gene_plot</code> single-cell panel", "<code>vmax=80</code>", "omitted",
     "forced",
     "80 is a display clip tuned to the tutorial's own count distribution. Transplanting the "
     "literal number onto a different panel and chemistry would misrepresent our data. The grid "
     "panel carries no <code>vmax</code> in the tutorial either."),
    ("<code>ccinet_plot</code> layout", "<code>pos_1 = ccinet_plot(..., return_pos=True)</code>, "
     "then <code>pos=pos_1</code> on each per-LR network", "same",
     "match",
     "Restored after the original run omitted <code>return_pos</code>/<code>pos</code>, which "
     "laid out every per-LR network independently and made them impossible to compare by eye. "
     "The shared layout is the whole point of capturing the return value."),
    ("<code>figsize</code>", "<code>(20, 5)</code> and <code>(20, 8)</code>",
     "each panel sized to the measured 2.21:1 tissue aspect (4.5 × 10.0 in)",
     "forced",
     "<strong>The only tutorial argument not copied literally.</strong> Those figsizes are for a "
     "1.37:1 section; combined with the tutorial's <em>own</em> "
     "<code>set_aspect('equal')</code> they letterbox this TMA into roughly a tenth of the "
     "canvas and drop the legend on top of the upper cores. <code>figsize</code> is pure canvas "
     "— it changes no datum, statistic, selection or colour scale. Computed from "
     "<code>grid.obsm['spatial']</code> at runtime, not hardcoded."),
    ("Per-LR plot selection", "top-1 for <code>lr_result_plot</code>; top-2 for "
     "<code>ccinet_plot</code> / <code>lr_chord_plot</code>", "same, plus a separate "
     "<code>requested/</code> directory",
     "protocol",
     "The benchmark protocol requires GRN→SORT1 and ANXA1→FPR1 to be plotted whatever their "
     "rank. They go in their own directory so they can never be mistaken for stLearn's own "
     "ranking, and the filenames carry the rank (<code>rank21_</code>, <code>rank99_</code>) so "
     "the distinction survives being copied out of context."),
]

CATEGORY_LABEL = {
    "forced": ("Forced", "Our data or platform makes the tutorial's literal choice impossible or wrong."),
    "better-input": ("Better input", "We supply a higher-quality version of what the tutorial improvises."),
    "recommended": ("Author-recommended", "The tutorial's own comment recommends a different value than its example uses."),
    "match": ("Matches tutorial", "Included because an earlier pass of ours got it wrong."),
    "protocol": ("Benchmark protocol", "Required by the comparator-benchmark rules, not by stLearn."),
    "consequence": ("Consequence", "Not a choice — an outcome of another decision."),
    "measured": ("Measured, not assumed", "A risk we quantified rather than hand-waved."),
    "version-gap": ("Version gap", "The installed package disagrees with the published tutorial."),
    "safety": ("Safety guard", "Protects against a silent-wrong-answer bug in stLearn."),
    "no-op": ("Proven no-op", "Literally absent, but verified to change no result."),
    "stated-rule": ("Stated rule", "The tutorial gives no rule; we state ours explicitly."),
    "complete-iteration": ("Complete iteration", "Same call, full coverage instead of an arbitrary subset."),
}

TOP10 = summ.head(10)

# ------------------------------------------------------------------ HTML
CSS = """
:root{--bg:#fff;--fg:#16181d;--mut:#5b6270;--line:#e3e6ec;--accent:#1f5fa9;--accent-bg:#eef4fb;
--warn:#a4442c;--warn-bg:#fdf1ed;--ok:#1d6f42;--ok-bg:#edf7f0;--code-bg:#f6f7f9;--card:#fbfcfd}
@media (prefers-color-scheme:dark){:root{--bg:#101318;--fg:#e6e9ef;--mut:#98a1b0;--line:#262b34;
--accent:#7db3ec;--accent-bg:#152232;--warn:#e79479;--warn-bg:#2c1d18;--ok:#7fc99c;--ok-bg:#152417;
--code-bg:#171b22;--card:#141821}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px 120px}
header.top{border-bottom:1px solid var(--line);margin-bottom:8px;padding:44px 0 28px}
h1{font-size:31px;line-height:1.25;margin:0 0 10px;letter-spacing:-.02em}
.sub{color:var(--mut);font-size:16px;margin:0}
h2{font-size:23px;margin:52px 0 6px;padding-top:14px;border-top:1px solid var(--line);letter-spacing:-.01em}
h3{font-size:18px;margin:30px 0 6px}
h4{font-size:15px;margin:22px 0 4px;color:var(--mut);text-transform:uppercase;letter-spacing:.06em}
p{margin:11px 0}
code{background:var(--code-bg);padding:1px 5px;border-radius:4px;font-size:.88em;
font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
pre.code{background:var(--code-bg);border:1px solid var(--line);border-radius:8px;padding:14px 16px;
overflow-x:auto;font:13px/1.6 ui-monospace,SFMono-Regular,Menlo,monospace;margin:12px 0}
pre.code code{background:none;padding:0}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:14.5px;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:7px 10px;text-align:left;vertical-align:top}
th{background:var(--card);font-weight:600}
tbody tr:nth-child(even){background:var(--card)}
figure.fig{margin:22px 0;padding:0}
figure.fig img{width:100%;max-width:100%;height:auto;border:1px solid var(--line);border-radius:8px;
background:#fff;display:block}
figcaption{color:var(--mut);font-size:13.5px;margin-top:8px;line-height:1.55}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}
@media(max-width:760px){.grid2{grid-template-columns:1fr}}
.note,.warn,.ok{border-left:3px solid;border-radius:0 8px 8px 0;padding:12px 16px;margin:16px 0;font-size:15px}
.note{border-color:var(--accent);background:var(--accent-bg)}
.warn{border-color:var(--warn);background:var(--warn-bg)}
.ok{border-color:var(--ok);background:var(--ok-bg)}
.tag{display:inline-block;font-size:11.5px;padding:2px 8px;border-radius:99px;
background:var(--accent-bg);color:var(--accent);border:1px solid var(--line);white-space:nowrap;font-weight:600}
nav.toc{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);
padding:10px 0;margin-bottom:0;z-index:50;font-size:13.5px}
nav.toc a{color:var(--mut);text-decoration:none;margin-right:15px;white-space:nowrap;display:inline-block;padding:2px 0}
nav.toc a:hover{color:var(--accent)}
.tocinner{max-width:1080px;margin:0 auto;padding:0 24px;overflow-x:auto;white-space:nowrap}
.kv{display:grid;grid-template-columns:230px 1fr;gap:2px 16px;font-size:15px;margin:14px 0}
.kv div:nth-child(odd){color:var(--mut)}
.stat{display:flex;gap:26px;flex-wrap:wrap;margin:18px 0}
.stat div{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 18px;min-width:118px}
.stat b{display:block;font-size:25px;line-height:1.15;letter-spacing:-.02em}
.stat span{color:var(--mut);font-size:12.5px}
/* explorer */
#exp{border:1px solid var(--line);border-radius:12px;padding:18px;background:var(--card);margin:18px 0}
.ctrls{display:flex;gap:26px;flex-wrap:wrap;align-items:flex-end;margin-bottom:14px}
.ctrl label{display:block;font-size:12.5px;color:var(--mut);margin-bottom:3px}
.ctrl input[type=range]{width:210px;display:block}
.ctrl input[type=text]{padding:6px 9px;border:1px solid var(--line);border-radius:6px;
background:var(--bg);color:var(--fg);width:190px;font-size:14px}
#cnt{font-weight:600;color:var(--accent)}
#tblwrap{max-height:440px;overflow-y:auto;border:1px solid var(--line);border-radius:8px;background:var(--bg)}
#tbl{margin:0;font-size:13.5px;display:table}
#tbl th{position:sticky;top:0;z-index:2;cursor:pointer;user-select:none}
#tbl th:hover{color:var(--accent)}
#tbl tbody tr{cursor:pointer}
#tbl tbody tr:hover{background:var(--accent-bg)}
#tbl tbody tr.req td:first-child{font-weight:700;color:var(--accent)}
#pin{margin-top:16px;display:none;border:1px solid var(--line);border-radius:10px;padding:14px;background:var(--bg)}
#pin h4{margin:0 0 8px;color:var(--fg);text-transform:none;letter-spacing:0;font-size:16px}
.hm{border-collapse:collapse;font-size:11px;display:table;width:auto;margin:6px 0}
.hm td,.hm th{border:1px solid var(--line);padding:3px 5px;text-align:center;font-weight:400}
.hm th{background:none;font-size:10.5px;color:var(--mut);font-weight:600}
.hm th.rot{writing-mode:vertical-rl;transform:rotate(180deg);padding:5px 2px}
.hm td.sig{outline:2px solid var(--accent);outline-offset:-2px;font-weight:700}
.legend{font-size:12px;color:var(--mut);margin-top:8px}
footer{margin-top:70px;padding-top:22px;border-top:1px solid var(--line);color:var(--mut);font-size:13.5px}
"""

JS = """
const N=DATA.cts.length;
const fmt=n=>n.toLocaleString();
let sortKey='nSig',sortDir=-1;
const el=id=>document.getElementById(id);
function rows(){
  const mS=+el('sSpots').value, mC=+el('sCci').value, q=el('q').value.trim().toUpperCase();
  let r=[];
  for(let i=0;i<DATA.lrs.length;i++){
    if(DATA.nSig[i]<mS) continue;
    if(DATA.nCci[i]<mC) continue;
    if(q && !DATA.lrs[i].toUpperCase().includes(q)) continue;
    r.push(i);
  }
  r.sort((a,b)=>{
    if(sortKey==='lr') return sortDir*DATA.lrs[a].localeCompare(DATA.lrs[b]);
    return sortDir*(DATA[sortKey][a]-DATA[sortKey][b]);
  });
  return r;
}
function render(){
  const r=rows();
  el('cnt').textContent=fmt(r.length);
  el('vSpots').textContent=el('sSpots').value;
  el('vCci').textContent=el('sCci').value;
  const req=new Set(DATA.requested);
  let h='';
  for(const i of r){
    h+=`<tr data-i="${i}" class="${req.has(DATA.lrs[i])?'req':''}"><td>${DATA.lrs[i].replace('_',' &rarr; ')}</td>`
     +`<td>${i+1}</td><td>${fmt(DATA.nSpots[i])}</td><td>${fmt(DATA.nSig[i])}</td>`
     +`<td>${fmt(DATA.nSigP[i])}</td><td>${DATA.nCci[i]}</td></tr>`;
  }
  el('tb').innerHTML=h;
}
function colour(v,max){
  if(v<=0||max<=0) return 'transparent';
  const t=Math.sqrt(v/max);
  return `rgba(31,95,169,${(0.08+0.72*t).toFixed(3)})`;
}
function pin(i){
  const c=DATA.counts[i], p=DATA.pvals[i];
  const max=Math.max(...c);
  let h=`<h4>${DATA.lrs[i].replace('_',' &rarr; ')} &mdash; significant sender &rarr; receiver edges</h4>`;
  h+=`<p style="font-size:13.5px;color:var(--mut);margin:2px 0 8px">Rank ${i+1} of ${DATA.lrs.length} `
   + `&middot; ${fmt(DATA.nSig[i])} significant spots &middot; ${DATA.nCci[i]} significant cell-type edges. `
   + `Rows send, columns receive. Cell colour = interaction count; <b>outlined + bold</b> = permutation p &lt; 0.05.</p>`;
  h+='<table class="hm"><tr><th></th>';
  for(const ct of DATA.cts) h+=`<th class="rot">${ct}</th>`;
  h+='</tr>';
  for(let a=0;a<N;a++){
    h+=`<tr><th style="text-align:right">${DATA.cts[a]}</th>`;
    for(let b=0;b<N;b++){
      const k=a*N+b, v=c[k], sig=p[k]<0.05;
      h+=`<td class="${sig?'sig':''}" style="background:${colour(v,max)}" `
       + `title="${DATA.cts[a]} &rarr; ${DATA.cts[b]}: count ${fmt(v)}, p=${p[k]}">${v?fmt(v):'&middot;'}</td>`;
    }
    h+='</tr>';
  }
  h+='</table>';
  const sig=[];
  for(let a=0;a<N;a++) for(let b=0;b<N;b++) if(p[a*N+b]<0.05) sig.push(`${DATA.cts[a]} -> ${DATA.cts[b]}\\t${c[a*N+b]}`);
  h+=`<div class="legend">${sig.length} significant edges. Copyable list:</div>`
   + `<textarea readonly rows="4" style="width:100%;margin-top:5px;font:11.5px ui-monospace,Menlo,monospace;`
   + `background:var(--code-bg);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:7px">`
   + sig.join('\\n') + `</textarea>`;
  el('pin').innerHTML=h;
  el('pin').style.display='block';
}
document.addEventListener('DOMContentLoaded',()=>{
  ['sSpots','sCci','q'].forEach(id=>el(id).addEventListener('input',render));
  el('tb').addEventListener('click',e=>{
    const tr=e.target.closest('tr'); if(tr) pin(+tr.dataset.i);
  });
  document.querySelectorAll('#tbl th[data-k]').forEach(th=>{
    th.addEventListener('click',()=>{
      const k=th.dataset.k;
      if(sortKey===k) sortDir=-sortDir; else {sortKey=k;sortDir=k==='lr'?1:-1;}
      render();
    });
  });
  render();
  pin(DATA.lrs.indexOf(DATA.requested[0])>=0?DATA.lrs.indexOf(DATA.requested[0]):0);
});
"""


def dev_rows():
    out = []
    for item, tut, ours, cat, why in DEVIATIONS:
        lab, _ = CATEGORY_LABEL[cat]
        out.append(
            f"<tr><td><strong>{item}</strong><br><span class='tag'>{lab}</span></td>"
            f"<td>{tut}</td><td>{ours}</td><td>{why}</td></tr>")
    return "\n".join(out)


def cellmap_rows():
    return "\n".join(
        f"<tr><td><code>cell{ours}_…</code></td><td><code>In[{inp}]</code></td>"
        f"<td><code>{E(call)}</code></td><td>{what}</td></tr>"
        for ours, inp, call, what in CELL_MAP)


def top10_rows():
    return "\n".join(
        f"<tr><td>{i+1}</td><td><code>{lr.replace('_',' → ')}</code></td>"
        f"<td>{int(r.n_spots):,}</td><td>{int(r.n_spots_sig):,}</td>"
        f"<td>{int(r.n_spots_sig_pval):,}</td><td>{int(r.n_cci_sig_cell_type)}</td></tr>"
        for i, (lr, r) in enumerate(TOP10.iterrows()))


def dom_rows():
    tot = int(dom.sum())
    return "\n".join(
        f"<tr><td>{ct}</td><td>{int(dom[ct]):,}</td><td>{100*dom[ct]/tot:.1f}%</td>"
        f"<td>{int(sender_tot[ct]):,}</td><td>{100*sender_tot[ct]/sender_tot.sum():.1f}%</td>"
        f"<td><strong>{sender_tot[ct]/max(dom[ct],1):,.0f}</strong></td></tr>"
        for ct in dom.sort_values(ascending=False).index)


def dropped_rows():
    return "\n".join(
        f"<tr><td><code>{E(f)}</code></td><td><span class='tag'>{cat}</span></td><td>{why}</td></tr>"
        for f, cat, why in DROPPED_FIGS)


BS = chr(92)  # backslash; kept out of f-string expressions so this file parses on Python < 3.12
REPRO = f"""
# the 45-minute analysis
python scripts/comparators/stlearn/run_stlearn.py {BS}
    --h5ad {run_man["h5ad"]} --out-dir {run_man["out_dir"]} {BS}
    --lrs results/comparators/stlearn/cellchatdb2_lrs.txt {BS}
    --cell-type-col cell_type --count-layer counts {BS}
    --n-row {run_man["n_row"]} --n-col {run_man["n_col"]} --seed {run_man["seed"]}

# the tutorial-exact figure set (~50 s, replay only)
python scripts/comparators/stlearn/plot_stlearn_tutorial.py {BS}
    --h5ad {run_man["h5ad"]} --run-dir {run_man["out_dir"]} {BS}
    --out-dir {run_man["out_dir"]}/plots_tutorial

# stage-2 quantitative export (~3 s, replay only)
python scripts/comparators/stlearn/export_stlearn_quant.py --run-dir {run_man["out_dir"]}

# this report
python scripts/comparators/stlearn/build_stlearn_report.py {BS}
    --run-dir {run_man["out_dir"]} --fig-dir {A.fig_dir} --out {A.out}
"""

HTML = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>stLearn on the GBM Xenium TMA — CellChatDB v2 — full method walkthrough</title>
<style>{CSS}</style>
</head>
<body>
<nav class="toc"><div class="tocinner">
<a href="#names">0 · Figure names</a>
<a href="#input">1 · Input</a>
<a href="#prep">2 · Preprocessing</a>
<a href="#grid">3 · Gridding</a>
<a href="#db">4 · LR database</a>
<a href="#stage1">5 · Stage 1</a>
<a href="#adj">6 · adj_pvals</a>
<a href="#stage2">7 · Stage 2</a>
<a href="#persist">8 · Outputs</a>
<a href="#scope">9 · Figure scope</a>
<a href="#figs">10 · Findings</a>
<a href="#explorer">11 · Explorer</a>
<a href="#dev">12 · Deviations</a>
<a href="#limits">13 · Limits</a>
<a href="#repro">14 · Reproduce</a>
</div></nav>

<div class="wrap">
<header class="top">
<h1>stLearn on the GBM/LGG Xenium TMA</h1>
<p class="sub">CellChatDB v2 tier · complete method walkthrough, with a stated reason for every
departure from the authors' default workflow · stLearn v1.4.1</p>
</header>

<div class="stat">
<div><b>100,197</b><span>cells in</span></div>
<div><b>12,562</b><span>grid spots</span></div>
<div><b>51.3 µm</b><span>spot size</span></div>
<div><b>526</b><span>LR pairs tested</span></div>
<div><b>482</b><span>with ≥20 sig. spots</span></div>
<div><b>45.4 min</b><span>wall time</span></div>
</div>

<div class="note"><strong>What this document is.</strong> A step-by-step account of the stLearn
run — every call, every argument, and for each place we depart from the published tutorial, the
reason. It is written to be checkable: where a claim could have been assumed, it says how it was
measured instead. Scope is deliberately narrow: <strong>one run</strong>, the
<code>cellchatdb2</code> tier on GBM. The <code>default</code> tier and LUAD have not been run.</div>

<h2 id="names">0 · Reading the figure names</h2>
<p>Every figure is named <code>cell&lt;NN&gt;_…</code>. <strong><code>NN</code> is the code cell of
the official stLearn Xenium tutorial that the figure comes from.</strong> The naming exists so that
any figure can be traced back to the exact tutorial step that produced it, without having to
reverse-engineer it from the picture.</p>

<p>One precision point, because the number could otherwise mislead: the tutorial is published as a
rendered HTML page, not a runnable notebook. To pin the call contract we extracted every
<code>&lt;div class="highlight"&gt;&lt;pre&gt;</code> block from that page, which yields
<strong>84 blocks</strong> — 49 code and 35 output/prompt blocks, since the rendering interleaves
them. <code>NN</code> is the index into <em>that extraction</em>. It is <em>not</em> the notebook's
own <code>In[k]</code> execution number. Both are given below.</p>

<table><thead><tr><th>Our prefix</th><th>Tutorial's own prompt</th><th>Call</th><th>What it shows</th></tr></thead>
<tbody>{cellmap_rows()}</tbody></table>

<p>Files in <code>figures/requested/</code> carry an extra <code>rank&lt;N&gt;_</code> segment —
e.g. <code>cell67_lr_result_rank21_GRN_SORT1.png</code> — recording where that interaction sat in
stLearn's <em>own</em> ranking. That way a figure requested by us can never be mistaken for one
stLearn surfaced on its own, even if the file is copied out of context.</p>

<h2 id="input">1 · Input</h2>
<div class="kv">
<div>Source</div><div><code>data/xenium_mm_final_cell_id.h5ad</code></div>
<div>Shape</div><div>100,197 cells × 5,119 genes, human Xenium 5K</div>
<div>Design</div><div>TMA, <strong>13 cores</strong> — 7 high-grade, 6 low-grade (<code>obs['tma_id']</code>)</div>
<div><code>X</code></div><div>log-normalised — <strong>not usable directly</strong>, see §2</div>
<div><code>layers['counts']</code></div><div>raw counts — this is what the run uses</div>
<div>Cell types</div><div><code>obs['cell_type']</code>, 9 types</div>
<div>Coordinates</div><div><code>obsm['spatial']</code>, already µm</div>
</div>

<h2 id="prep">2 · Preprocessing, step by step</h2>
<p>The tutorial's preprocessing is short, and every step of it matters to the statistics
downstream. Ours follows it in the same order; four steps differ, and each is flagged inline.</p>

{code('''
# 1. counts into X  --  DEVIATION (forced)
adata.X = adata.layers["counts"].copy()

# 2. drop gene names containing "_"  --  DEVIATION (forced)
adata = adata[:, [g for g in adata.var_names if "_" not in g]]     # 5,119 -> 5,098

# 3. cell labels  --  DEVIATION (better input: tutorial runs Leiden here)
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

# 4. Visium-style uns["spatial"]  --  DEVIATION (forced: plain h5ad has none)
adata.uns["spatial"] = {"sample": {"images": {"hires": <1x1 placeholder>},
                                   "use_quality": "hires",
                                   "scalefactors": {"tissue_hires_scalef": 1.0,
                                                    "spot_diameter_fullres": 15}}}
adata = st.convert_scanpy(adata)          # -> obs["imagecol"], obs["imagerow"]

# 5-8. verbatim from the tutorial
st.pp.filter_genes(adata, min_counts=10)  # In[7]
st.pp.filter_cells(adata, min_counts=10)  # In[7]
adata.raw = adata                         # In[9]
st.pp.normalize_total(adata)              # In[12]  -- library size ONLY, no log1p
''')}

<p>Result: <strong>100,197 cells (none lost to QC) × 5,096 genes</strong>. Of the 23 genes removed,
21 are control probes dropped by name and 2 fall below <code>min_counts=10</code>.</p>

<div class="warn"><strong>The normalisation step is the one that would fail silently.</strong>
The tutorial's own text: <em>“No log1p or shrinking to make genes of similar expression range. In
our case, for calling hotspots, we want genes to be more separate, since we select background genes
with similar expression levels to detect hotspots.”</em> The permutation null in §5 works by
sampling random gene pairs <em>matched on expression level</em> to the real ligand and receptor.
Log-transforming compresses genes toward each other, so the “matched” background stops being
matched — and the resulting p-values are wrong without anything erroring. Our <code>X</code> is
log-normalised, which is precisely why the run reads <code>layers['counts']</code> instead.</div>

<h3>Why Leiden is skipped</h3>
<p>The tutorial spends three calls (<code>run_pca</code> → <code>neighbors</code> →
<code>clustering.leiden</code> at resolution 1.05) producing <code>obs['leiden']</code>. It does
this because its public demo ships no cell-type annotation. <code>use_label</code> only requires
<em>a</em> per-cell categorical column; we have expert annotations, so running Leiden would
substitute a worse label for a better one. Everything downstream is identical in shape — the same
column feeds <code>grid()</code> and <code>run_cci()</code>.</p>

<h2 id="grid">3 · Gridding — and why not 125 × 125</h2>
<p>stLearn does not score cells. <code>st.tl.cci.grid()</code> aggregates them onto a regular grid
via <code>np.histogram2d(xs, ys, bins=[n_col, n_row])</code>, drops spots containing zero cells, and
from then on the <strong>spot is the unit of inference</strong>. It also stores a per-spot
cell-type proportion matrix in <code>uns[label]</code>, which is what makes the mixture mode in §7
possible.</p>

<p>The tutorial passes <code>n_ = 125</code> to both axes. Its own markdown calls this a
resolution/compute trade-off — an author-declared knob, not a tuned biological parameter. Copying
the number is not the same as copying the choice:</p>

<table><thead><tr><th></th><th>Tutorial (breast)</th><th>Ours (GBM TMA)</th></tr></thead><tbody>
<tr><td>Tissue extent</td><td>7,521 × 5,471 µm</td><td>7,484 × 16,483 µm</td></tr>
<tr><td>Aspect</td><td>1.37 : 1</td><td><strong>2.20 : 1</strong></td></tr>
<tr><td>Bins</td><td>125 × 125</td><td>146 × 321</td></tr>
<tr><td>Spot size</td><td>60.2 × 43.8 µm (rectangular)</td><td><strong>51.3 × 51.3 µm (square)</strong></td></tr>
<tr><td>Spot area</td><td>2,637 µm²</td><td>2,632 µm²</td></tr>
<tr><td>Occupied spots</td><td>14,364</td><td>12,562</td></tr>
<tr><td>Cells per occupied spot</td><td>11.4</td><td>7.98 (median 5, max 46)</td></tr>
</tbody></table>

<p>Copying <code>125 × 125</code> onto a 2.20:1 tissue would produce <strong>59.9 × 131.9 µm</strong>
spots — a neighbourhood 2.2× larger in y than in x, which is a different aggregation, not the same
one. We preserved the quantity the tutorial's choice actually fixes on its own data, the
<strong>spot area</strong>, and let the bin counts follow. The two areas agree to 0.2%.</p>

<div class="note"><strong>On the ~51 µm coincidence.</strong> ALARMIST uses a 50 µm patch on this
same dataset. The 51.3 µm here falls out of matching the tutorial's spot area and is
<em>not</em> tuned toward ALARMIST. No comparator parameter in this benchmark is set by reference to
ALARMIST — the whole point of the exercise is that each method keeps its own spatial definition.</div>

<h3>The neighbourhood is not the grid</h3>
<p>A point worth separating, because it is easy to conflate: <code>distance=250</code> µm is the
signalling neighbourhood, applied via a cKDTree over spot centroids, and it is
<strong>independent of grid resolution</strong>. The grid controls aggregation granularity; the
distance controls how far a spot looks. We changed the former and left the latter at the authors'
default, per the benchmark rule that no method's neighbourhood is harmonised toward another's.</p>

<h3>Does a 250 µm radius leak across TMA cores?</h3>
<p>A legitimate worry: the 13 cores come from different patients, and a neighbourhood that bridges
two cores would mix them. We measured it rather than assuming: <strong>0 cross-core neighbour pairs
at 200 µm, and 39 pairs involving 43 cells (0.04% of cells) at 250 µm.</strong> Core bounding boxes
approach within 150 µm, but the cores are round, so the cells themselves never get that close.
Running per core was rejected for a separate reason — the smallest core has 819 cells → ~15
occupied spots, below <code>min_spots=20</code>, so it would drop out entirely.</p>

<h2 id="db">4 · The ligand–receptor database</h2>
<p>The tutorial loads stLearn's bundled resource:</p>
{code("lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')   # 2,293 pairs")}
<p>This run uses the <code>cellchatdb2</code> tier instead — CellChatDB v2, the resource ALARMIST
uses — so that the LR database is not a confounder when the two methods are compared. Converting it
is where stLearn's most consequential limitation shows up.</p>

<div class="warn"><strong>stLearn cannot represent multi-subunit complexes.</strong> Its LR format
is the single string <code>"LIGAND_RECEPTOR"</code>, so <code>_</code> is the ligand/receptor
separator and cannot <em>also</em> mean “subunit of”. CellChatDB v2 writes complexes with
underscores (<code>TGFBR1_TGFBR2</code>), so those rows are unrepresentable:
<strong>1,859 of 3,233 rows (57.5%) are dropped</strong>, leaving 1,371 single-gene pairs, of which
527 are on the Xenium panel and <strong>526</strong> survive stLearn's own expression filter.</div>

<p>We drop them rather than expand them. Expanding <code>TGFB1</code> + <code>TGFBR1_TGFBR2</code>
into {{TGFB1_TGFBR1, TGFB1_TGFBR2}} would assert two interactions that CellChat does not claim exist
independently, and would double-count one. Dropping is a faithful representation of what the method
can and cannot test; expanding would be us inventing data.</p>

<p>The cost is not cosmetic. All 40 <code>WNT3_*</code> rows in CellChatDB v2 have complex receptors
(<code>FZD*_LRP5/6</code>), so <strong>stLearn cannot test a single WNT3 interaction</strong> — while
WNT3–FZD*/LRP6 was CytoSignal's entire top-6 on this same dataset. For scale, ALARMIST's GBM feature
space contains 712 unique LRIs, <strong>205 of them (29%) multi-subunit</strong> and therefore
invisible here. Any disagreement between the methods is partly structural, not biological, and
should be reported that way.</p>

<h2 id="stage1">5 · Stage 1 — where is an LR pair co-expressed more than chance?</h2>
{code('''
st.tl.cci.run(grid, lrs,
              min_spots=20,      # tutorial value, unchanged
              distance=250,      # tutorial value, unchanged
              n_pairs=10000,     # DEVIATION: tutorial uses 1000 but recommends ~10,000
              n_cpus=n_cpus,
              random_state=0)
''')}
<p>For each LR pair, the score of a spot is the co-expression of ligand and receptor across that
spot and its neighbours within 250 µm. The null is built by sampling <code>n_pairs</code> random
gene pairs <strong>matched on expression level</strong> to the real ligand and receptor, scoring
them identically, and asking where the real pair exceeds its own background. This yields a per-spot
p-value per pair, BH-corrected across the LRs tested in that spot.</p>

<p>Two properties worth stating plainly, because they determine what the method can and cannot
claim: the unit of inference is the <strong>spot</strong>, never the cell; and stage 1 is
<strong>cell-type agnostic</strong> — cell identity does not enter until stage 2.</p>

<h4>Why 10,000 rather than 1,000</h4>
<p>The tutorial's own inline comment reads
<code># Number of random pairs to generate; low as example, recommend ~10,000</code>. We follow the
recommendation rather than the example. The cost argument that usually discourages this does not
apply: backgrounds are cached per <em>gene</em>, not per pair, so 100× the pairs costs roughly 7.6×
the time. <code>n_perms</code> in stage 2 is raised from 100 to 1,000 for the identical reason —
the tutorial's comment there says <code>recommend ~1000</code>.</p>

<p>Outputs: <code>obsm['lr_scores']</code>, <code>['p_vals']</code>, <code>['p_adjs']</code>,
<code>['-log10(p_adjs)']</code>, <code>['lr_sig_scores']</code> (each spots × LR), and the ranking
in <code>uns['lr_summary']</code>. The obsm column order matches the <code>lr_summary</code> row
order — stLearn re-sorts both together, which is what makes the persisted CSVs in §8 correctly
labelled.</p>

<h2 id="adj">6 · <code>adj_pvals</code> — the step we can prove does nothing</h2>
<p>The tutorial calls, at <code>In[22]</code>:</p>
{code("st.tl.cci.adj_pvals(grid, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')")}
<p>Our analysis run did not. Rather than assume that mattered or assume it didn't, we checked the
source and then ran the experiment.</p>

<p><strong>From source:</strong> <code>run()</code> defaults to <code>adj_method='fdr_bh'</code>
and <code>pval_adj_cutoff=0.05</code>, and its permutation routine performs the correction per spot
— i.e. exactly <code>correct_axis='spot'</code>. <code>adj_pvals</code>' own docstring states
<em>“Default settings of this function are already run in st.tl.cci.run.”</em></p>

<p><strong>From experiment:</strong> we loaded the saved grid, applied the call verbatim, and
compared everything before and after.</p>

<div class="ok"><strong>Result: a complete no-op on every statistic.</strong> All five
<code>obsm</code> matrices are element-wise identical once realigned, and all six
<code>lr_summary</code> columns are identical. The only effect is that
<code>np.argsort</code> — which is <em>unstable</em> — re-permutes <strong>179 of 526</strong> LR
pairs, and every one of those swaps is <em>within a group tied on <code>n_spots_sig</code></em>.
The top-3 are unchanged. So the results already on disk, produced without the call, are the
tutorial's results.</div>

<p>The tutorial-exact plotting pass makes the call anyway, for literal fidelity, and asserts the
no-op at runtime — so if a future stLearn version changes this behaviour, the assertion fires
instead of the report quietly becoming wrong.</p>

<h2 id="stage2">7 · Stage 2 — which cell types sit in those hotspots?</h2>
{code('''
st.tl.cci.run_cci(grid, "cell_type",
                  min_spots=2,
                  spot_mixtures=True,     # a spot may count as several cell types
                  cell_prop_cutoff=0.1,   # ...if that type holds >10% of it
                  sig_spots=True,         # restrict to stage-1 significant spots
                  n_perms=1000,           # DEVIATION: tutorial uses 100, recommends ~1000
                  random_state=0)
''')}
<p>Restricted to significant spots, this counts cell-type → cell-type edges across the
neighbourhood graph, then permutes the cell-type labels 1,000 times to test whether a given
directed pair is over-represented in that LR's hotspots. <code>spot_mixtures=True</code> matters a
great deal here — see the mGAM discussion in §10.</p>

<h4>A silent-wrong-answer bug we guarded against</h4>
<p><code>get_data_for_counting</code> binds cell types to deconvolution columns by
<strong>substring, first match wins</strong>: <code>[ct in col for col in cols]</code>. Our
vocabulary contains both <code>mGAM</code> and <code>non-mGAM</code>, and
<code>'mGAM' in 'non-mGAM'</code> is <code>True</code> — so the binding is order-dependent, and a
mis-binding would corrupt <em>every</em> interaction reported for that cell type without raising
anything. It happens to be correct here because alphabetical ordering puts <code>mGAM</code> first,
so it claims its own column. The runner nevertheless verifies the mapping for all 9 labels and
aborts if any label binds to a column that is not its own. This is a latent stLearn bug that would
bite any dataset with nested cell-type names.</p>

<h4>A version gap we had to work around</h4>
<p><code>run_cci</code> does <code>adata.obs[label].values.astype(str)</code> and feeds the result
to an <code>@njit</code> kernel. Under pandas 3.0.5 / numpy 2.4.6 that returns a StringDtype/object
array that numba 0.66 cannot type, raising
<code>TypingError: non-precise type array(pyobject, 1d, C)</code> and hard-blocking the entire
cell-type half of the workflow. Coercing the column to object dtype makes the <em>same
expression</em> return a plain <code>&lt;U</code> array, which numba accepts; the categorical is
restored immediately afterwards because the plotting functions need <code>.cat</code>. Labels are
byte-identical across the coercion, asserted at runtime. The alternative — pinning
<code>pandas&lt;3</code> — was rejected because the coercion is provably label-preserving and keeps
the environment on current scanpy.</p>

<h2 id="persist">8 · What landed on disk</h2>
<p>The tutorial is a notebook and persists nothing. The benchmark protocol requires score matrices,
p-values and per-spot assignments in a re-readable form, so the runner writes them out. Two gaps in
the original run were closed afterwards by a replay-only exporter that reads
<code>grid.h5ad</code> read-only and recomputes nothing:</p>

<table><thead><tr><th>Gap</th><th>Cause</th><th>Fix</th></tr></thead><tbody>
<tr><td>Stage-2 <strong>p-values</strong> never written</td>
<td>The persistence loop covered only <code>lr_cci_&lt;label&gt;</code> and
<code>per_lr_cci_&lt;label&gt;</code>, so <code>per_lr_cci_pvals_&lt;label&gt;</code> and both raw
count sets stayed inside the h5ad. The per-LR CSVs carried interaction counts with no significance
attached.</td>
<td>526 p-value matrices + 526 raw-count matrices + <code>lr_cci_raw_cell_type.csv</code> exported.</td></tr>
<tr><td><code>lr_summary.csv</code> had 3 of 6 columns</td>
<td>It was written <em>before</em> <code>run_cci()</code>, which appends three
<code>*_&lt;label&gt;</code> columns — the stage-2 ranking.</td>
<td>Rewritten with all 6 columns; the 3 pre-existing columns verified byte-identical first.</td></tr>
</tbody></table>

<h2 id="scope">9 · Which figures — and why exactly these</h2>
<p>stLearn ships <strong>two</strong> cell–cell-interaction vignettes, and they do not have the same
figure set. This matters more than it sounds: several figures that had been treated as “the standard
workflow” turn out to belong to the other vignette, and three belong to neither. Provenance was
established by grepping both vignettes, then adversarially re-checked by independent verifiers
looking for any occurrence inside loops, multi-panel figures or aliased imports.</p>

<table><thead><tr><th><code>st.pl.*</code></th><th>Xenium vignette</th><th>Generic Visium vignette</th><th>In this report</th></tr></thead>
<tbody>
<tr><td><code>cluster_plot</code></td><td>✅ ×5</td><td>✅</td><td>yes</td></tr>
<tr><td><code>feat_plot</code></td><td>✅ In[16]</td><td>❌</td><td>yes</td></tr>
<tr><td><code>gene_plot</code></td><td>✅ In[17]</td><td>❌</td><td>yes</td></tr>
<tr><td><code>lr_summary</code></td><td>✅ ×2</td><td>✅</td><td>yes</td></tr>
<tr><td><code>lr_result_plot</code></td><td>✅ In[24]</td><td>✅</td><td>yes</td></tr>
<tr><td><code>cci_check</code></td><td>✅</td><td>✅</td><td>yes</td></tr>
<tr><td><code>ccinet_plot</code></td><td>✅ In[28]</td><td>✅</td><td>yes</td></tr>
<tr><td><code>lr_chord_plot</code></td><td>✅ In[29]</td><td>✅</td><td>yes</td></tr>
<tr><td><code>lr_diagnostics</code></td><td>❌</td><td>✅</td><td>no</td></tr>
<tr><td><code>lr_n_spots</code></td><td>❌</td><td>✅</td><td>no</td></tr>
<tr><td><code>cci_map</code>, <code>lr_cci_map</code></td><td>❌</td><td>✅</td><td>no</td></tr>
<tr><td><code>lr_plot</code></td><td>❌</td><td>✅ ×8</td><td>no</td></tr>
<tr><td><code>lr_go</code></td><td>❌</td><td>✅</td><td>no — also unavailable, needs R + clusterProfiler</td></tr>
<tr><td><code>het_plot</code>, <code>grid_plot</code>, <code>deconvolution_plot</code></td>
<td>❌</td><td>❌</td><td>no — <strong>in neither vignette</strong></td></tr>
</tbody></table>

<p>The Xenium vignette is the right one for this data: the generic one runs on
<code>st.datasets.visium_sge</code> and never calls <code>st.tl.cci.grid()</code> at all. So this
report contains the Xenium vignette's set, call-for-call, plus the two requested interactions —
23 figures + 6 requested. The excluded figures do exist on disk under
<code>plots/</code> and <code>plots_full/</code>; they were kept, not deleted, but they are not part
of the authors' Xenium workflow.</p>

<h2 id="figs">10 &middot; The biological findings</h2>
<p>Organised by claim rather than by tutorial order, and carrying only the figures that
<em>support a claim</em>. Pipeline checks and plain cell-type-location maps are left out and listed
in &sect;10.5; the complete 29-figure tutorial set remains on disk either way.</p>

<h3>10.1 &nbsp;The myeloid compartment is the signalling hub &mdash; far above its spatial footprint</h3>
<p>This is the most robust finding in the run, and it survives the obvious objection.</p>

<table><thead><tr><th>Cell type</th><th>Dominant in N spots</th><th>% of spots</th>
<th>Signal sent (pooled)</th><th>% of signal</th><th>Per dominant spot</th></tr></thead>
<tbody>{dom_rows()}</tbody></table>

<div class="note"><strong>Read the last column.</strong> mGAM dominates <strong>1.8%</strong> of
spots but sends <strong>11.4%</strong> of the signal &mdash; <strong>12.7&times; more interaction
per dominant spot than Glial-Neuronal</strong>, the most abundant type. non-mGAM is nearly
identical (1,511 vs 1,571). mGAM is a significant sender for <strong>288 of the 526</strong> tested
pairs.</div>

{fig("cell42_composition_mGAM.png", "<b>mGAM composition</b> &mdash; per-spot proportion (colour fixed at <code>vmax=1</code>), the spots mGAM wins, and the individual mGAM cells.")}

<p>The composition panel explains the mechanism and the caveat at once: mGAM cells are spread
across essentially every core, yet the proportion panel is almost uniformly low and mGAM wins only
221 spots. <strong>Everywhere, dominant almost nowhere.</strong> That is exactly the regime in which
a majority-label grid erases a cell type &mdash; and it is why
<code>spot_mixtures=True</code> is not optional here. Without the mixture mode every mGAM claim
above would rest on 221 spots rather than on every spot in which mGAM holds &gt;10%.</p>

<p>The obvious objection is that this is an abundance artefact of the permutation. It is not:</p>

{fig("cell79_cci_check.png", "<b>Abundance diagnostic</b> &mdash; cell-type frequency (bars) against number of CCI-LR interactions (line).")}

<p>Bars are cell-type frequency, the line is number of CCI-LR interactions. Had the permutation
failed to control for abundance the line would track the bars. It runs against them &mdash;
Glial-Neuronal is the most frequent type and has the <em>fewest</em> interactions, while mGAM is
rare and sits among the highest.</p>

<h3>10.2 &nbsp;Complement is the top signal, and mGAM is its dominant sender</h3>

{fig("cell58_lr_summary_n50.png", "<b>Top 50 ligand&ndash;receptor pairs</b> by number of significant spots, of 526 tested. All 526 have at least one significant spot; 482 have 20 or more.")}

<table><thead><tr><th>#</th><th>LR pair</th><th>Spots expressing</th><th>Significant spots</th>
<th>Sig. by raw p</th><th>Sig. cell-type edges</th></tr></thead>
<tbody>{top10_rows()}</tbody></table>

<p><code>C3&rarr;C3AR1</code> takes rank 1 with 1,244 significant spots and <code>C3&rarr;CR2</code>
rank 5 &mdash; two of the top five are complement. mGAM sends <code>C3&rarr;C3AR1</code> to all nine
cell types, with a pooled count of <strong>76,648</strong>, the largest single number in the run.
Grouping mGAM's other significant outputs by hand gives a coherent tumour-associated myeloid
programme: <code>GAS6&rarr;MERTK</code> and <code>GAS6&rarr;AXL</code> (efferocytosis, TAM
receptors), <code>CSF1&rarr;CSF1R</code> (myeloid autocrine survival),
<code>HLA-DQA1&rarr;CD4</code> (antigen presentation), <code>LAIR1&rarr;LILRB4</code> and
<code>ADORA3&rarr;ENTPD1</code> (inhibitory / immunosuppressive), plus
<code>CXCL12&rarr;CXCR4</code>.</p>

<div class="warn"><strong>That grouping is ours, not stLearn's.</strong> The method tests 526
<em>independent</em> hypotheses and ranks them; it has no representation of &ldquo;these
interactions run together&rdquo;. The programme above exists because a human sorted 15 rows of a
table into a story. This is the structural difference from a factorisation method, and it is worth
being explicit about rather than presenting the grouping as a result.</div>

{fig("cell67_lr_result_C3_C3AR1.png", "<b>C3 &rarr; C3AR1</b> across three statistics: raw co-expression score, &minus;log10(adjusted p), and score masked to significant spots.")}

<p>Three statistics for the top pair. Reading them together is what separates
&ldquo;co-expressed&rdquo; from &ldquo;co-expressed more than its own expression-matched
background&rdquo;: the masked panel on the right is far sparser than the raw score on the left.</p>

<div class="grid2">
{fig("cell81_ccinet_all.png", "<b>All 526 pairs pooled</b> (<code>min_counts=30</code>).")}
{fig("cell81_ccinet_C3_C3AR1.png", "<b>C3 &rarr; C3AR1 alone</b>, same node layout.")}
</div>
<p>Left: all 526 pairs pooled. Right: <code>C3&rarr;C3AR1</code> alone. Both use the same node
layout (captured with <code>return_pos=True</code>), so they can be compared directly.</p>

<h3>10.3 &nbsp;The default ranking is partly an artefact &mdash; rank by breadth instead</h3>
<div class="warn">Homophilic pairs &mdash; CNTN2&ndash;CNTN2, GJA1&ndash;GJA1, NCAM1&ndash;NCAM1,
JAM3&ndash;JAM3 &mdash; occupy much of the top of the table above. When ligand and receptor are the
<strong>same gene</strong>, &ldquo;co-expression&rdquo; is trivially satisfied wherever that gene is
expressed, so a high rank partly reflects <em>this gene is everywhere in this tissue</em> rather
than a communication event.</div>
<p>Ranking by <em>breadth</em> instead &mdash; the number of significant cell-type pairs, the last
column above and a sort option in &sect;11 &mdash; surfaces a different and more immunologically
interesting set: <code>HLA-DQA1&rarr;CD4</code> (56 edges), <code>PVR&rarr;CD226</code>,
<code>IL16&rarr;CD4</code>, <code>ANXA1&rarr;FPR1</code>, <code>ICAM1&rarr;ITGAL</code>. Complement
is at the top of both rankings, which is part of why finding 10.2 is credible.</p>

<h3>10.4 &nbsp;Both ALARMIST motif-1 directions are independently significant</h3>
<p>GRN&rarr;SORT1 and ANXA1&rarr;FPR1 are the ALARMIST motif-1 mGAM &hArr; MES-like loop, plotted
whatever their rank. Both are single-gene pairs, so both survive the complex-dropping step of
&sect;4 and are genuinely testable here &mdash; which is not true of most of CellChatDB.</p>

<table><thead><tr><th>Interaction</th><th>stLearn rank</th><th>Percentile</th><th>Spots expressing</th>
<th>Significant spots</th><th>Sig. cell-type edges</th></tr></thead>
<tbody>
<tr><td><code>GRN &rarr; SORT1</code></td><td><strong>21 / 526</strong></td><td>top 4.0%</td>
<td>8,763</td><td>267</td><td>24</td></tr>
<tr><td><code>ANXA1 &rarr; FPR1</code></td><td><strong>99 / 526</strong></td><td>top 18.8%</td>
<td>5,721</td><td>122</td><td>48</td></tr>
</tbody></table>

<div class="grid2">
{fig("requested/cell67_lr_result_rank21_GRN_SORT1.png", "<b>GRN &rarr; SORT1</b> (rank 21) &mdash; co-expression is widespread, significance is focal.")}
{fig("requested/cell81_ccinet_rank21_GRN_SORT1.png", "<b>GRN &rarr; SORT1</b> cell-type network.")}
{fig("requested/cell67_lr_result_rank99_ANXA1_FPR1.png", "<b>ANXA1 &rarr; FPR1</b> (rank 99).")}
{fig("requested/cell81_ccinet_rank99_ANXA1_FPR1.png", "<b>ANXA1 &rarr; FPR1</b> cell-type network.")}
</div>

<div class="ok"><strong>Checked against stLearn's own stage-2 permutation p-values, both directions
of the loop are recovered:</strong>
<code>GRN&rarr;SORT1, mGAM &rarr; MES-like</code> &mdash; count 1,130, p&nbsp;=&nbsp;0.000, 10th of
that pair's 24 significant edges. <code>ANXA1&rarr;FPR1, MES-like &rarr; mGAM</code> &mdash; count
447, p&nbsp;=&nbsp;0.000, 11th of 48. Two methods sharing no unit, no null and no neighbourhood
agree on both edges.</div>

<p>Two caveats stated plainly, because the concordance is easy to overstate. Neither edge is the
<em>top</em> edge of its pair &mdash; for GRN&rarr;SORT1, mGAM&rarr;Glial-Neuronal (2,069) and
mGAM&rarr;OPC-like (1,402) both score higher, and stLearn calls mGAM a significant sender to
<em>all nine</em> types. And these counts rest on the mixture mode rather than on mGAM-majority
spots (&sect;10.1). stLearn supports the edges; it does not single them out.</p>

<h3>10.5 &nbsp;Figures deliberately left out of this report</h3>
<p>All of these were produced, are part of the tutorial-exact set, and remain on disk in
<code>plots_tutorial/</code> and in this report's <code>figures/</code> directory. They are omitted
here because they describe the pipeline or the tissue rather than support a finding.</p>
<table><thead><tr><th>Figure</th><th>Category</th><th>Why not in the report</th></tr></thead>
<tbody>{dropped_rows()}</tbody></table>

<h2 id="explorer">11 · Interactive: all 526 pairs</h2>
<p>Every tested pair, with the stage-2 sender → receiver matrix behind each one. Move the sliders to
re-filter live; click any row to pin its matrix. Defaults show everything, reproducing the numbers
in §10. The two requested interactions are highlighted.</p>

<div id="exp">
<div class="ctrls">
<div class="ctrl"><label>Min significant spots: <b id="vSpots">0</b></label>
<input type="range" id="sSpots" min="0" max="1244" value="0" step="1"></div>
<div class="ctrl"><label>Min significant cell-type edges: <b id="vCci">0</b></label>
<input type="range" id="sCci" min="0" max="56" value="0" step="1"></div>
<div class="ctrl"><label>Search</label><input type="text" id="q" placeholder="e.g. GRN, C3, CXCL"></div>
<div class="ctrl"><label>Passing</label><div><span id="cnt">526</span> / 526 pairs</div></div>
</div>
<div id="tblwrap"><table id="tbl">
<thead><tr>
<th data-k="lr">LR pair</th><th data-k="rank">Rank</th><th data-k="nSpots">Spots expr.</th>
<th data-k="nSig">Sig. spots</th><th data-k="nSigP">Sig. (raw p)</th><th data-k="nCci">Sig. edges</th>
</tr></thead><tbody id="tb"></tbody></table></div>
<div id="pin"></div>
</div>
<p style="font-size:13.5px;color:var(--mut)">Sort by <em>Sig. edges</em> to see the breadth ranking
discussed in §10 — it is a different list from the default spot-count ranking.</p>

<h2 id="dev">12 · Every deviation, with its reason</h2>
<p>The complete register. Categories: <em>forced</em> (the tutorial's literal choice is impossible
or wrong on our data), <em>author-recommended</em> (the tutorial's own comment recommends a
different value than its example uses), <em>benchmark protocol</em> (required by the comparison, not
by stLearn), and so on. Nothing below is a preference.</p>

<table><thead><tr><th style="min-width:170px">Item</th><th style="min-width:200px">Tutorial</th>
<th style="min-width:200px">This run</th><th>Reason</th></tr></thead>
<tbody>{dev_rows()}</tbody></table>

<h2 id="limits">13 · What this run cannot answer</h2>
<ul>
<li><strong>Nothing at single-cell resolution.</strong> The unit is a 51 µm spot throughout; stLearn
never scores a cell. Per-cell differential expression, per-cell programme loadings and anything
built on them are out of reach by construction.</li>
<li><strong>No between-condition test.</strong> stLearn has no native multi-sample or differential
mode. The 7 high-grade vs 6 low-grade cores cannot be contrasted without hand-rolling a test, which
the benchmark rules forbid — the point is to report each method as its authors ship it.</li>
<li><strong>No complexes</strong> — 57.5% of CellChatDB v2, including every WNT3 interaction (§4).</li>
<li><strong>No notion of co-occurring programmes.</strong> 526 independent hypotheses, ranked. mGAM
is a significant sender for 288 of them; the method offers no way to say which subset constitutes
one coherent programme.</li>
<li><strong>Lymphoid is underpowered</strong> — 32 dominant spots.</li>
<li><strong>The <code>default</code> tier has not been run</strong>, so this run does not yet
demonstrate stLearn on its own recommended resource; nor has LUAD, at either tier.</li>
<li><strong>Database mismatch with the ALARMIST fit.</strong> This run used the re-exported
<code>CellChatDBv2.0.human.csv</code>; <code>results/GBM/</code> was produced with
<code>CellChatDBv2.0.human.old.csv</code>, and 1,120 of 3,218 LR keys differ. GRN→SORT1 and
ANXA1→FPR1 are in both, so the comparison above is safe, but a full LRI-level comparison would
cross that boundary.</li>
</ul>

<h2 id="repro">14 · Reproducing this</h2>
{code(REPRO)}
<p>All three post-analysis steps open <code>grid.h5ad</code> read-only and recompute nothing. The
plotting pass rebuilds the cell-level object deterministically (QC + <code>normalize_total</code>)
because the side-by-side panels need it, and asserts the rebuild against
<code>run_manifest.json</code> before drawing — a silent preprocessing drift aborts rather than
producing a mismatched figure.</p>

<table><thead><tr><th>Reference</th><th>Path</th></tr></thead><tbody>
<tr><td>Method notes, algorithm, gotchas</td><td><code>scripts/comparators/METHODS.md</code> §stLearn</td></tr>
<tr><td>Tutorial call contract</td><td><code>scripts/comparators/stlearn/NOTES.md</code></td></tr>
<tr><td>Deviation register (source of §12)</td><td><code>scripts/comparators/stlearn/DEVIATIONS.md</code></td></tr>
<tr><td>Source figures</td><td><code>{E(A.run_dir)}/plots_tutorial/</code></td></tr>
<tr><td>Matrices, p-values, per-spot assignments</td><td><code>{E(A.run_dir)}/data/</code></td></tr>
</tbody></table>

<footer>
<p>stLearn {run_man.get("stlearn","1.4.1")} · scanpy {run_man.get("scanpy","")} ·
Python {run_man.get("python","")} · seed {run_man.get("seed",0)} ·
{run_man.get("n_lr_pairs_input",0):,} LR pairs in, {run_man.get("n_lr_pairs_tested",0)} tested ·
analysis wall time {run_man.get("wall_min","?")} min.</p>
<p>Generated by <code>scripts/comparators/stlearn/build_stlearn_report.py</code> — that script, not
this HTML, is the version-controlled artifact. {inlined} figures embedded, {linked} linked.</p>
</footer>
</div>

<script>const DATA={json.dumps(PAYLOAD, separators=(",", ":"))};</script>
<script>{JS}</script>
</body></html>
"""

os.makedirs(os.path.dirname(os.path.abspath(A.out)), exist_ok=True)
open(A.out, "w").write(HTML)
mb = os.path.getsize(A.out) / 1048576
log(f"wrote {A.out} ({mb:.1f} MB)")
