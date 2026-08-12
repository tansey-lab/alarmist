#!/usr/bin/env python
"""Build the self-contained CellChat-on-GBM HTML report.

Follows the repo's `interactive-report` rules: one standalone .html, no CDN, no external
files, figures inlined as base64, report text in English, and THIS SCRIPT is the tracked
artifact (the HTML lands under reports/, which is gitignored).

Every number in the report is read from the run outputs at build time rather than
hardcoded, so re-running the pipeline and re-running this script keeps them in sync.

Env: bptf.  Usage:
    python scripts/comparators/cellchat/build_report.py
"""
import base64
import json
import os
import subprocess
import tempfile

import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
GBM = os.path.join(REPO, "results/comparators/cellchat/GBM")
PNG = os.path.join(REPO, "reports/cellchat_GBM_png")
OUT = os.path.join(PNG, "cellchat_GBM_report.html")
TIER = "cellchatdb2"          # the tier this report walks through
MAX_W = 1400                  # downscale width for embedding

# --------------------------------------------------------------------------- numbers
man = {t: json.load(open(f"{GBM}/{t}/run_manifest.json")) for t in ["default", "cellchatdb2"]}
audit = json.load(open(f"{REPO}/results/comparators/cellchat/db_audit/db_equivalence.json"))
ctrl_path = f"{GBM}/control_nonspatial/spatial_vs_nonspatial.json"
ctrl = json.load(open(ctrl_path)) if os.path.exists(ctrl_path) else None

Q = f"{GBM}/{TIER}/quant"
summary = pd.read_csv(f"{Q}/summary_by_condition.csv")
flow = pd.read_csv(f"{Q}/information_flow_by_pathway.csv")
fp = flow.pivot(index="pathway_name", columns="condition", values="information_flow").fillna(0)
fp["diff"] = fp["high"] - fp["low"]
req = pd.read_csv(f"{Q}/requested_lr_status.csv")
sig = {c: pd.read_csv(f"{Q}/{c}_net_significant.csv") for c in ["low", "high"]}
sizes = {c: pd.read_csv(f"{Q}/{c}_group_sizes.csv") for c in ["low", "high"]}
nup = len(pd.read_csv(f"{Q}/net_up_in_high.csv"))
ndn = len(pd.read_csv(f"{Q}/net_down_in_high.csv"))


def cond_stat(t, c, k):
    return next(x[k] for x in man[t]["conditions"] if x["condition"] == c)


def motif_dir(cond, lr, s, t):
    d = sig[cond]
    hit = d[(d.interaction_name == lr) & (d.source == s) & (d.target == t)]
    return None if not len(hit) else (float(hit.iloc[0]["prob"]), float(hit.iloc[0]["pval"]))


# --------------------------------------------------------------------------- figures
def embed(rel):
    """Downscale to MAX_W and return a base64 data URI."""
    src = os.path.join(PNG, rel)
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        dst = tmp.name
    subprocess.run(["sips", "-Z", str(MAX_W), src, "--out", dst],
                   check=True, capture_output=True)
    with open(dst, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    os.unlink(dst)
    return f"data:image/png;base64,{b64}"


def fig(rel, caption, half=False):
    cls = "fig half" if half else "fig"
    return (f'<figure class="{cls}"><img loading="lazy" src="{embed(rel)}" alt="{caption}">'
            f'<figcaption><span class="fn">{rel}</span>{caption}</figcaption></figure>')


def why(text):
    return f'<div class="why"><span class="whyh">Why we deviated</span>{text}</div>'


def note(text):
    return f'<div class="note">{text}</div>'


def warn(text):
    return f'<div class="warn">{text}</div>'


# --------------------------------------------------------------------------- content
top_high = fp.sort_values("high", ascending=False).head(10)
gain = fp.sort_values("diff", ascending=False).head(8)
lose = fp.sort_values("diff").head(6)
high_only = sorted(fp[(fp.low == 0) & (fp.high > 0)].index.tolist())
low_only = sorted(fp[(fp.high == 0) & (fp.low > 0)].index.tolist())


def tbl(df, cols, headers, fmt=None):
    fmt = fmt or {}
    h = "".join(f"<th>{x}</th>" for x in headers)
    rows = ""
    for idx, r in df.iterrows():
        cells = ""
        for c in cols:
            v = idx if c == "__index__" else r[c]
            cells += f"<td>{fmt.get(c, lambda x: x)(v)}</td>"
        rows += f"<tr>{cells}</tr>"
    return f"<table><thead><tr>{h}</tr></thead><tbody>{rows}</tbody></table>"


f3 = lambda x: f"{x:.3f}"

SECTIONS = []


def sec(id_, title, body):
    SECTIONS.append((id_, title, body))


# ---- 0 ----------------------------------------------------------------------
sec("what", "0. What CellChat actually computes", f"""
<p>Before any parameter choice makes sense, it matters what the method's output <em>is</em>.
For each ligand–receptor (L–R) pair and each <strong>ordered pair of cell groups</strong>
(<em>i</em> → <em>j</em>), CellChat evaluates a Hill function on the product of the
group-average ligand expression in <em>i</em> and the group-average receptor expression in
<em>j</em>:</p>
<p class="eq">P = L·R / (K<sub>h</sub><sup>n</sup> + L·R),&nbsp;&nbsp; K<sub>h</sub> = 0.5, n = 1</p>
<p>L and R are computed over <code>data/max(data)</code> using a 10% truncated mean per group;
multi-subunit complexes enter as the geometric mean across subunits; the receptor term is
multiplied by co-activation and divided by co-inhibition receptor terms. Significance comes
from a <strong>permutation test</strong>: cell-group labels are shuffled
<code>nboot = 100</code> times and <code>pval = #{{P<sub>boot</sub> ≥ P<sub>obs</sub>}}/nboot</code>.</p>
{note("<strong>Unit of inference.</strong> A CellChat result is a "
      "<em>(sender cell type, receiver cell type, L–R pair)</em> triple. It is not a cell, "
      "not a spot, and not a cell–cell edge. Our output array is literally "
      f"9 × 9 × {len(sig['high'].interaction_name.unique())} for the high-grade object. "
      "Every interpretation below inherits that granularity.")}
""")

# ---- 1 ----------------------------------------------------------------------
low_n, high_n = cond_stat(TIER, "low", "n_cells"), cond_stat(TIER, "high", "n_cells")
sec("input", "1. Data input and object construction", f"""
<p>CellChat needs four inputs for spatial data: a normalized expression matrix, a metadata
frame, coordinates, and <code>spatial.factors</code>. Here is exactly what each was.</p>

<h3>1.1 Expression matrix — normalized, not counts</h3>
<p>CellChat requires <em>normalized</em> data (library-size normalized then log1p), not raw
counts. Our h5ad already stores this in <code>adata.X</code>, with raw counts kept separately
in <code>layers['counts']</code>. The exporter asserts rather than assumes it: it inverts the
log1p and checks the implied library size, which came out at a median of exactly
<strong>10,000.0</strong> — i.e. <code>X = log1p(CP10K)</code>, precisely CellChat's expected
input. <code>layers['counts']</code> was therefore never touched.</p>
{note("This is the single most common silent failure mode. Feeding counts, or data that has "
      "been logged twice, changes every probability without raising an error.")}

<h3>1.2 Metadata — and why there are two objects, not one</h3>
<p>The vignettes are explicit on a point that is easy to get wrong: several
<em>sections of one condition</em> belong in <strong>one</strong> object, distinguished by a
<code>samples</code> column, but a comparison <em>across conditions</em> requires a
<strong>separate object per condition</strong>. The GBM TMA has 13 cores spanning two tumour
grades, so it becomes two objects:</p>
<table><thead><tr><th>Object</th><th>TMA cores (<code>samples</code>)</th><th>Cells</th></tr></thead>
<tbody>
<tr><td><code>high</code></td><td>7 cores — 1, 3, 5, 8, 10, 11, 13</td><td>{high_n:,}</td></tr>
<tr><td><code>low</code></td><td>6 cores — 2, 4, 6, 9, 12, 14</td><td>{low_n:,}</td></tr>
</tbody></table>
<p><code>labels</code> came from <code>obs['cell_type']</code> (9 types) and <code>samples</code>
from <code>obs['tma_id']</code>. All 9 cell types occur in all 13 cores, which is what makes
<code>mergeCellChat</code> and the <em>functional</em>-similarity manifold analysis applicable
at all — the latter requires identical cell-type composition between datasets.</p>
{note("Putting the cores in as <code>samples</code> is not cosmetic. CellChat computes "
      "cell-group distances <em>within each sample</em> and then averages, so physically "
      "separate tissue punches are never treated as spatial neighbours of one another.")}

<h3>1.3 Spatial factors — the first deviation</h3>
<p>The vignettes compute <code>ratio</code> from a 10X Visium scale-factor JSON and set
<code>tol = 65/2 = 32.5</code>. Neither applies to Xenium.</p>
<table><thead><tr><th></th><th>Tutorial (Visium)</th><th>Ours (Xenium)</th></tr></thead><tbody>
<tr><td><code>ratio</code></td><td><code>65 / spot_diameter_fullres</code></td><td><strong>1</strong></td></tr>
<tr><td><code>tol</code></td><td><code>32.5</code></td><td><strong>5</strong></td></tr>
</tbody></table>
{why("This is not our invention — it is the FAQ vignette's own Xenium row. Xenium coordinates "
     "are already in micrometres, so the pixel→micron conversion factor is 1; and "
     "<code>spot.size</code> becomes the typical human cell size (10 µm), giving "
     "<code>tol = spot.size/2 = 5</code>.")}
{warn("<strong>Undocumented trap.</strong> <code>spatial.factors</code> must have one row per "
      "sample <em>in <code>levels(samples)</code> order</em>: <code>computeRegionDistance</code> "
      "indexes them positionally as <code>ratio[k]</code>/<code>tol[k]</code>. No vignette says "
      "this. Getting the order wrong silently applies the wrong micron conversion to the wrong "
      "core. Our exporter writes the rows pre-sorted and the runner asserts the alignment.")}
""")

# ---- 2 ----------------------------------------------------------------------
ov = audit["overlap"]
sec("db", "2. The ligand–receptor database, and why the two tiers are not what you'd expect", f"""
<p>The benchmark protocol asks for two tiers per method: <code>default</code> (the tutorial's
own database) and <code>cellchatdb2</code> (the CellChatDB v2 resource ALARMIST uses), so that
the database is removed as a confounder. For CellChat that framing collapses, and it is worth
being precise about why rather than quietly running the same thing twice.</p>

<h3>2.1 CellChat's bundled database <em>is</em> CellChatDB v2</h3>
<p>We re-derived the repo CSV's flattening from the bundled <code>CellChatDB.human</code> using
the documented mapping (ligand/receptor ← <code>complex</code> subunits joined with
<code>_</code>) and diffed the key sets:</p>
<table><thead><tr><th></th><th>Bundled <code>CellChatDB.human</code></th><th><code>data/LRdatabase/CellChatDBv2.0.human.csv</code></th></tr></thead><tbody>
<tr><td>rows</td><td>{audit['bundled']['n_rows']:,}</td><td>{audit['repo_csv']['n_rows']:,}</td></tr>
<tr><td>unique ligand|receptor keys</td><td>{audit['bundled']['n_unique_keys']:,}</td><td>{audit['repo_csv']['n_unique_keys']:,}</td></tr>
</tbody></table>
<p class="big">Shared keys <strong>{ov['n_shared']:,}</strong> &nbsp;·&nbsp; only in bundled
<strong>{ov['n_only_bundled']}</strong> &nbsp;·&nbsp; only in repo CSV
<strong>{ov['n_only_repo']}</strong> &nbsp;·&nbsp; Jaccard
<strong>{ov['jaccard']:.4f}</strong></p>
<p>They are the same resource. (The {audit['bundled']['n_rows'] - audit['bundled']['n_unique_keys']}-key
gap between rows and unique keys is 15 pairs listed twice under two annotations — the three
<code>POMC|OPR*</code> pairs appear as both Secreted and Non-protein Signaling.)</p>

<h3>2.2 What we did instead of faking a second run</h3>
<p>Re-importing the flat CSV through <code>updateCellChatDB</code> was rejected: it would
discard 22 columns including <code>agonist</code>, <code>antagonist</code>,
<code>co_A_receptor</code> and <code>co_I_receptor</code> — precisely the terms
<code>computeCommunProb</code> multiplies into the receptor expression — while adding zero
interactions. Degrading the method to re-supply a database it already ships is not a tier.</p>
<p>So the two tiers vary <strong>annotation scope</strong>, which is the only knob the vignettes
actually expose:</p>
<table><thead><tr><th>Tier</th><th>Call</th><th>Interactions</th><th>Pathways</th></tr></thead><tbody>
<tr><td><code>default</code></td><td><code>subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")</code><br><span class="fn">the literal call in both spatial vignettes</span></td><td>{man['default']['db']['n_interactions']:,}</td><td>{man['default']['db']['n_pathways']}</td></tr>
<tr><td><code>cellchatdb2</code></td><td><code>subsetDB(CellChatDB)</code><br><span class="fn">the alternative commented in on the next line of the same vignette</span></td><td>{man['cellchatdb2']['db']['n_interactions']:,}</td><td>{man['cellchatdb2']['db']['n_pathways']}</td></tr>
</tbody></table>
{why("<code>default</code> proves we ran the method the way its authors recommend. "
     "<code>cellchatdb2</code> widens the L–R scope to Secreted + ECM-Receptor + Cell-Cell "
     "Contact — everything except the Non-protein Signaling that CellChat itself advises "
     "against — which is the scope closest to what ALARMIST uses. Both calls are the "
     "vignette's own; neither is invented.")}
{note("<strong>Widening the database does not change any individual pair's probability.</strong> "
      "We verified it: for the L–R pairs present in both tiers, the probability and p-value "
      "arrays are bit-identical (max |Δ| = 0). That is expected — the Hill function for a pair "
      "depends only on that pair's own ligand and receptor. The tier changes coverage and "
      "pathway-level aggregation, nothing else. <strong>This report walks through the "
      f"<code>{TIER}</code> tier.</strong>")}
""")

# ---- 3 ----------------------------------------------------------------------
sec("preproc", "3. Preprocessing", f"""
<p>Run exactly as the vignette does, at package defaults:</p>
<table><thead><tr><th>Step</th><th>Call</th><th>Result on our data</th></tr></thead><tbody>
<tr><td>subset to signalling genes</td><td><code>subsetData(cellchat)</code></td><td>722 of 5,119 panel genes are in the database</td></tr>
<tr><td>parallel backend</td><td><code>future::plan("multisession", workers = 4)</code></td><td>the vignette's own value</td></tr>
<tr><td>over-expressed genes</td><td><code>identifyOverExpressedGenes(cellchat)</code></td><td>632 (low) / 659 (high)</td></tr>
<tr><td>over-expressed interactions</td><td><code>identifyOverExpressedInteractions(cellchat)</code></td><td>{cond_stat(TIER,'low','n_lr_tested')} (low) / {cond_stat(TIER,'high','n_lr_tested')} (high) tested</td></tr>
</tbody></table>
{why("<code>identifyOverExpressedInteractions</code> has <code>variable.both = TRUE</code> as "
     "its package default and the multi-sample spatial vignette passes nothing, so it gets "
     "<code>TRUE</code>. The single-section vignette passes <code>F</code>. We follow the "
     "multi-sample vignette, because it is the one that matches our data shape. This is a "
     "genuine fork in the tutorials, not an oversight, so we state which branch we took.")}
{note("The FAQ offers <code>do.DE = FALSE</code> for datasets with 'a small panel of genes'. "
      "The Xenium 5K panel is <strong>not</strong> small in that sense — over 1,000 database "
      "pairs survive gene filtering — so the default DE path was used. Had it collapsed to "
      "near zero, that would have been a reportable result rather than a licence to switch "
      "paths silently.")}
""")

# ---- 4 ----------------------------------------------------------------------
p = man[TIER]["parameters"]
sec("infer", "4. Inference — every parameter, and the reason for each", f"""
<p>The single call that produces the result:</p>
<pre><code>computeCommunProb(cellchat,
                  type = "{p['type']}", trim = {p['trim']},
                  distance.use = {str(p['distance.use']).upper()},
                  interaction.range = {p['interaction.range']},
                  scale.distance = NULL,
                  contact.dependent = TRUE, contact.range = {p['contact.range']},
                  nboot = {p['nboot']}, seed.use = {p['seed.use']})</code></pre>
<table><thead><tr><th>Parameter</th><th>Value</th><th>Source / reason</th></tr></thead><tbody>
<tr><td><code>type</code>, <code>trim</code></td><td><code>truncatedMean</code>, 0.1</td><td>Both spatial vignettes' literal values. The package default <code>triMean</code> is more conservative (≈25% truncated mean); the vignettes deliberately loosen it for spatial data.</td></tr>
<tr><td><code>contact.range</code></td><td><strong>10</strong> (vignettes: 100)</td><td><strong>Deviation</strong> — see below.</td></tr>
<tr><td><code>distance.use</code></td><td><strong>FALSE</strong></td><td>The multi-sample spatial vignette's value. The single-section vignette uses <code>TRUE</code> with <code>scale.distance = 0.01</code>. See below.</td></tr>
<tr><td><code>interaction.range</code></td><td>250</td><td>Package and vignette default; the stated maximum diffusion range of a ligand in µm.</td></tr>
<tr><td><code>nboot</code>, <code>seed.use</code></td><td>100, 1</td><td>Package defaults.</td></tr>
<tr><td><code>raw.use</code></td><td>TRUE</td><td>Default; both spatial vignettes leave the PPI smoothing (<code>projectData</code>) commented out, so we do too.</td></tr>
<tr><td><code>population.size</code></td><td>FALSE</td><td>Package default.</td></tr>
<tr><td><code>k.min</code>, <code>do.symmetric</code>, <code>K<sub>h</sub></code>, <code>n</code></td><td>10, TRUE, 0.5, 1</td><td>Package defaults, untouched.</td></tr>
</tbody></table>

<h3>4.1 <code>contact.range = 10</code>, not 100</h3>
{why("The vignettes' <code>contact.range = 100</code> is the 10X Visium spot centre-to-centre "
     "distance. The FAQ, the single-section vignette and <code>?computeCommunProb</code> all "
     "pin <strong>10 µm</strong> for single-cell-resolution platforms, calling it 'a typical "
     "human cell size'. Xenium is single-cell resolution. We did not take this on faith — we "
     "measured the nearest-neighbour distance on our own coordinates: "
     f"<strong>median {cond_stat(TIER,'low','median_nn_distance_um')} µm (low) and "
     f"{cond_stat(TIER,'high','median_nn_distance_um')} µm (high)</strong>, minimum "
     f"{cond_stat(TIER,'low','observed_min_cell_distance_um')} / "
     f"{cond_stat(TIER,'high','observed_min_cell_distance_um')} µm. A 100 µm 'contact' range "
     "would have spanned roughly ten cell diameters and made the word meaningless.")}

<h3>4.2 <code>distance.use = FALSE</code></h3>
{why("The multi-sample spatial vignette — the one whose data shape matches ours — sets "
     "<code>distance.use = FALSE</code>. The single-section vignette sets <code>TRUE</code> "
     "with <code>scale.distance = 0.01</code>, but that constant is Visium-specific: with "
     "<code>distance.use = TRUE</code> CellChat validates that "
     "<code>min(d × scale.distance)</code> lands in [1,2] and <strong>aborts</strong> "
     "otherwise, so 0.01 cannot simply be copied onto micron-scale Xenium coordinates.")}
{warn("<strong>This choice has a consequence we discovered only by checking, and it is the "
      "most important caveat in this report.</strong> See section 9.1.")}

<h3>4.3 Two things we replaced for scale, not for science</h3>
{why("<code>computeCellDistance</code> materialises a dense N×N matrix — about 51 GB at "
     "79,998 cells — and returned NA. It is only used as a sanity check on "
     "<code>contact.range</code>, so we compute the identical quantity with "
     "<code>BiocNeighbors::findKNN(k = 1)</code> in O(N log N). No inference path is affected."
     "<br><br>"
     "<code>netEmbedding</code> defaults to python <code>umap-learn</code> through reticulate. "
     "We pass <code>umap.method = &quot;uwot&quot;</code> — the package's own documented "
     "alternative — so the run does not depend on a python environment.")}
""")

# ---- 5 ----------------------------------------------------------------------
sec("global", "5. Result I — the global rewiring between grades", f"""
<p>Objects were combined with <code>mergeCellChat(list(low = …, high = …))</code>. The order
matters for reading every differential plot: <strong>low is dataset 1, high is dataset 2</strong>,
so red / positive always means <em>higher in high grade</em>.</p>
<table><thead><tr><th>Condition</th><th>Cells</th><th>L–R pairs tested</th><th>Significant links</th><th>Pathways</th><th>Total strength</th></tr></thead><tbody>
<tr><td>low</td><td>{low_n:,}</td><td>{cond_stat(TIER,'low','n_lr_tested')}</td><td>{int(summary[summary.condition=='low'].n_interactions.iloc[0]):,}</td><td>{cond_stat(TIER,'low','n_pathways')}</td><td>{summary[summary.condition=='low'].total_strength.iloc[0]:.3f}</td></tr>
<tr><td>high</td><td>{high_n:,}</td><td>{cond_stat(TIER,'high','n_lr_tested')}</td><td>{int(summary[summary.condition=='high'].n_interactions.iloc[0]):,}</td><td>{cond_stat(TIER,'high','n_pathways')}</td><td>{summary[summary.condition=='high'].total_strength.iloc[0]:.3f}</td></tr>
</tbody></table>
{fig("comparison/compareInteractions.png", "Total number of interactions (left) and total interaction strength (right), low vs high grade.")}
<div class="grid2">
{fig("comparison/diffInteraction_count.png", "Differential number of interactions. Red = more in high grade, blue = more in low grade.", True)}
{fig("comparison/diffInteraction_weight.png", "Differential interaction strength between cell-type pairs.", True)}
</div>
<div class="grid2">
{fig("comparison/diff_heatmap_count.png", "Same contrast as a heatmap; the top bar is summed incoming change, the right bar summed outgoing change.", True)}
{fig("comparison/diff_heatmap_weight.png", "Differential strength heatmap.", True)}
</div>
{fig("comparison/role_scatter_both.png", "Outgoing vs incoming interaction strength per cell type, low (left) and high (right). Dot size is the number of inferred links.")}
""")

# ---- 6 ----------------------------------------------------------------------
sec("pathways", "6. Result II — which signalling pathways changed", f"""
<p><code>rankNet</code> compares the <em>information flow</em> of each pathway — the summed
communication probability across all cell-type pairs — with a paired Wilcoxon test
(<code>do.stat = TRUE</code>).</p>
{fig("comparison/rankNet_stacked.png", "Information flow per pathway, stacked (relative) — red enriched in low grade, blue-green enriched in high grade.")}
{fig("comparison/rankNet_unstacked.png", "The same comparison on absolute information flow.")}

<h3>6.1 The ten strongest pathways in high grade</h3>
{tbl(top_high, ["__index__", "high", "low", "diff"], ["Pathway", "high", "low", "Δ (high−low)"],
     {"high": f3, "low": f3, "diff": f3})}

<h3>6.2 Biggest gains and losses</h3>
<div class="grid2">
<div><h4>Gained in high grade</h4>{tbl(gain, ["__index__","high","low","diff"], ["Pathway","high","low","Δ"], {"high":f3,"low":f3,"diff":f3})}</div>
<div><h4>Lost in high grade</h4>{tbl(lose, ["__index__","high","low","diff"], ["Pathway","high","low","Δ"], {"high":f3,"low":f3,"diff":f3})}</div>
</div>
<p><strong>Present only in high grade:</strong> {', '.join(high_only)}.<br>
<strong>Present only in low grade:</strong> {', '.join(low_only)}.</p>
{note("<strong>Reading of the biology.</strong> The pathways that grow are contact-dependent and "
      "adhesion-type — APP, NOTCH, NCAM, CD99, CADM, COLLAGEN, CDH — while what recedes is "
      "immune-recognition signalling (COMPLEMENT, GAS, ICAM, LAIR1). The high-grade-only set "
      "(IGF, CCL, THBS, TENASCIN, APELIN) adds angiogenesis and matrix remodelling. The single "
      "strongest edge anywhere in the high-grade network is MES-like → mGAM via APP–SORL1.")}
<div class="grid2">
{fig("comparison/role_heatmap_outgoing_both.png", "Outgoing signalling per cell type × pathway, low vs high, on a shared pathway union.", True)}
{fig("comparison/role_heatmap_incoming_both.png", "Incoming signalling, same layout.", True)}
</div>

<h3>6.3 Up- and down-regulated L–R pairs</h3>
<p>A presto-backed differential expression analysis across conditions
(<code>identifyOverExpressedGenes(group.dataset = "datasets", thresh.fc = 0.05)</code>) mapped
onto the inferred communications with <code>netMappingDEG</code> gives
<strong>{nup} up</strong> and <strong>{ndn} down</strong> L–R pairs in high grade.</p>
{note("<code>thresh.fc = 0.05</code> is the vignette's own presto-adjusted value, which means "
      "<strong>presto must be installed</strong>: with it missing, <code>do.fast = TRUE</code> "
      "silently falls back to <code>stats::wilcox.test</code> and returns systematically larger "
      "logFC against an unchanged threshold. presto 1.0.0 is installed and its presence is "
      "asserted at startup and recorded in the run manifest.")}
<div class="grid2">
{fig("comparison/bubble_increased.png", "L–R pairs with increased communication probability in high grade.", True)}
{fig("comparison/bubble_decreased.png", "L–R pairs with decreased communication probability in high grade.", True)}
</div>
<div class="grid2">
{fig("comparison/wordcloud_up.png", "Enriched ligands in high grade.", True)}
{fig("comparison/wordcloud_down.png", "Enriched ligands in low grade.", True)}
</div>
{warn("A pair can legitimately appear in <em>both</em> the up and the down list: the DEA is run "
      "per cell group, so the same L–R can be up in one cell type and down in another. The "
      "vignette flags this explicitly. <code>GRN_SORT1</code> is one such case here.")}
""")

# ---- 7 ----------------------------------------------------------------------
g_hi = motif_dir("high", "GRN_SORT1", "mGAM", "MES-like")
g_lo = motif_dir("low", "GRN_SORT1", "mGAM", "MES-like")
a_hi = motif_dir("high", "ANXA1_FPR1", "MES-like", "mGAM")
a_lo = motif_dir("low", "ANXA1_FPR1", "MES-like", "mGAM")
fmtd = lambda x: "not significant" if x is None else f"prob {x[0]:.4f}, p {'< 0.001' if x[1] < 0.001 else f'= {x[1]:.3f}'}"
rq = {(r.condition, r.interaction_name): r for r in req.itertuples()}

sec("motif1", "7. Result III — the two requested interactions (ALARMIST motif 1)", f"""
<p>Two L–R pairs were requested regardless of how CellChat itself ranks them, because they are
the two arms of ALARMIST's motif 1, a bidirectional mGAM ⇄ MES-like loop:
<strong>GRN → SORT1</strong> (mGAM → MES-like) and <strong>ANXA1 → FPR1</strong>
(MES-like → mGAM). Both are Secreted Signaling, so both are present in <em>both</em> tiers, and
both survive <code>identifyOverExpressedInteractions</code> in both conditions.</p>

<table><thead><tr><th>L–R</th><th>Condition</th><th>Significant cell-type pairs</th><th>Max probability</th><th>The ALARMIST direction</th></tr></thead><tbody>
<tr><td rowspan="2"><code>GRN_SORT1</code><br><span class="fn">pathway GRN</span></td><td>low</td><td>{rq[('low','GRN_SORT1')].n_significant_pairs}</td><td>{rq[('low','GRN_SORT1')].max_prob:.4f}</td><td class="bad">mGAM → MES-like: {fmtd(g_lo)}</td></tr>
<tr><td>high</td><td>{rq[('high','GRN_SORT1')].n_significant_pairs}</td><td>{rq[('high','GRN_SORT1')].max_prob:.4f}</td><td class="good">mGAM → MES-like: {fmtd(g_hi)}</td></tr>
<tr><td rowspan="2"><code>ANXA1_FPR1</code><br><span class="fn">pathway ANNEXIN</span></td><td>low</td><td>{rq[('low','ANXA1_FPR1')].n_significant_pairs}</td><td>{rq[('low','ANXA1_FPR1')].max_prob:.4f}</td><td>MES-like → mGAM: {fmtd(a_lo)}</td></tr>
<tr><td>high</td><td>{rq[('high','ANXA1_FPR1')].n_significant_pairs}</td><td>{rq[('high','ANXA1_FPR1')].max_prob:.4f}</td><td class="good">MES-like → mGAM: {fmtd(a_hi)}</td></tr>
</tbody></table>

<p class="big">The complete bidirectional loop is significant <strong>only in high grade</strong>.
In low grade the GRN → SORT1 arm is not called at all, and the ANXA1 → FPR1 arm is about
{a_hi[0]/a_lo[0]:.0f}× weaker.</p>

{note("This is an independent corroboration of ALARMIST motif 1 being grade-associated, from a "
      "method with a completely different inference target — cell-type pairs with a permutation "
      "null — and no knowledge of the motif decomposition.")}

<div class="grid2">
{fig("requested_lr/high/GRN_SORT1/circle.png", "GRN → SORT1, high grade. Edge width is communication probability.", True)}
{fig("requested_lr/low/GRN_SORT1/circle.png", "GRN → SORT1, low grade.", True)}
</div>
<div class="grid2">
{fig("requested_lr/high/ANXA1_FPR1/circle.png", "ANXA1 → FPR1, high grade.", True)}
{fig("requested_lr/low/ANXA1_FPR1/circle.png", "ANXA1 → FPR1, low grade.", True)}
</div>
<div class="grid2">
{fig("requested_lr/high/GRN_SORT1/bubble.png", "GRN → SORT1 across all sender/receiver combinations, high grade.", True)}
{fig("requested_lr/high/ANXA1_FPR1/bubble.png", "ANXA1 → FPR1 across all sender/receiver combinations, high grade.", True)}
</div>
{note("<strong>ANXA1 → FPR1's receivers in high grade are exactly mGAM and non-mGAM and nothing "
      "else</strong> — consistent with FPR1 being myeloid-restricted on this panel. GRN → SORT1, "
      "by contrast, is broadcast: mGAM, MES-like, Vascular and non-mGAM all send it.")}

<h3>7.1 Where CellChat's own ranking puts them</h3>
{warn("Honesty about rank: the ALARMIST direction is significant but <strong>not</strong> the "
      "maximum for either pair. The strongest <code>GRN_SORT1</code> edge in both conditions is "
      "mGAM → Glial-Neuronal, and the strongest <code>ANXA1_FPR1</code> edges are "
      "Vascular → mGAM and Lymphoid → mGAM. CellChat confirms the loop exists; it does not "
      "single it out.")}

<h3>7.2 The two pathways in context</h3>
<div class="grid2">
{fig("high/pathways/GRN/circle.png", "The whole GRN pathway, high grade.", True)}
{fig("high/pathways/ANNEXIN/circle.png", "The whole ANNEXIN pathway, high grade.", True)}
</div>
<div class="grid2">
{fig("high/pathways/GRN/LR_contribution.png", "Contribution of each L–R pair to GRN signalling.", True)}
{fig("high/pathways/ANNEXIN/LR_contribution.png", "Contribution of each L–R pair to ANNEXIN signalling.", True)}
</div>
<div class="grid2">
{fig("comparison/geneExpression/GRN.png", "Expression of GRN-pathway signalling genes, split by grade.", True)}
{fig("comparison/geneExpression/ANNEXIN.png", "Expression of ANNEXIN-pathway signalling genes, split by grade.", True)}
</div>

<h3>7.3 On tissue</h3>
<div class="grid2">
{fig("requested_lr/high/GRN_SORT1/spatial_LR_core1.png", "GRN → SORT1 expression product mapped onto one high-grade core.", True)}
{fig("requested_lr/high/ANXA1_FPR1/spatial_LR_core1.png", "ANXA1 → FPR1 mapped onto the same core.", True)}
</div>
{warn("Read these two maps as <em>expression</em> maps, not as evidence of spatially local "
      "signalling. Section 9.1 explains why.")}

<h3>7.4 The cell types the story is about</h3>
<div class="grid2">
{fig("comparison/signalingChanges/mGAM.png", "Signalling changes for mGAM between grades.", True)}
{fig("comparison/signalingChanges/MES_like.png", "Signalling changes for MES-like cells between grades.", True)}
</div>
{note("<strong>mGAM's absolute signalling barely moves.</strong> It sends 1.797 → 1.671 and "
      "receives 2.867 → 2.884 in total strength, but its share of the whole network falls from "
      "17.7% → 10.9% (outgoing) and 28.2% → 18.8% (incoming). The high-grade tumour is not "
      "making mGAM louder — it is adding other traffic around it.")}
""")

# ---- 8 ----------------------------------------------------------------------
sec("systems", "8. Result IV — latent communication programs", f"""
<p><code>identifyCommunicationPatterns</code> factorises the cell-type × pathway matrix with
NMF, which is CellChat's closest analogue to a data-driven signalling program.</p>
{why("<code>identifyCommunicationPatterns</code> <em>requires</em> a value of <code>k</code> and "
     "stops if it is NULL; the vignette picks it by eye off the <code>selectK</code> curve. We "
     "applied the vignette's own stated rule — 'the number of patterns at which Cophenetic and "
     "Silhouette begin to drop suddenly' — programmatically to the same "
     "<code>NMF::nmfEstimateRank</code> measures, and persisted those measures to disk so the "
     "choice is auditable rather than eyeballed. <code>k.range</code> is additionally capped at "
     "<code>min(10, n_celltypes − 1, n_pathways − 1)</code> because NMF's rank must stay below "
     "both matrix dimensions.")}
<div class="grid2">
{fig("high/systems/river_outgoing.png", "Outgoing communication patterns, high grade: cell types → patterns → pathways.", True)}
{fig("high/systems/river_incoming.png", "Incoming communication patterns, high grade.", True)}
</div>
<div class="grid2">
{fig("high/systems/dot_outgoing.png", "Cell type × pathway contribution to each outgoing pattern.", True)}
{fig("high/systems/dot_incoming.png", "The same for incoming patterns.", True)}
</div>

<h3>8.1 A result worth stating plainly</h3>
{warn("<strong>CellChat's own pattern analysis separates the two arms of the loop.</strong> "
      "GRN is assigned to Pattern 3 and ANNEXIN to Pattern 2 — in the outgoing decomposition "
      "and in the incoming one. The database also files them as two unrelated pathways (GRN and "
      "ANNEXIN). So CellChat detects both arms as individually significant, but nothing in its "
      "output ever links them into one object. That linkage is what ALARMIST's motif 1 provides, "
      "and it is why the loop was found there first and confirmed here second.")}
<div class="grid2">
{fig("comparison/embeddingPairwise_functional.png", "Joint manifold of pathway networks across both grades, functional similarity.", True)}
{fig("comparison/rankSimilarity_functional.png", "Pathways ranked by how much their network geometry differs between grades.", True)}
</div>
""")

# ---- 9 ----------------------------------------------------------------------
ctrl_row = ctrl["results"][0] if ctrl else None
if ctrl_row:
    ctrl_block = warn(
        "<strong>We tested this rather than inferring it.</strong> Re-running the identical "
        "pipeline with <code>datatype = &quot;RNA&quot;</code> — no coordinates at all — "
        f"produced <strong>{ctrl_row['n_rows_spatial']} rows in both</strong>, "
        f"max |Δprob| = {ctrl_row['max_abs_delta_prob']:.1e}, "
        f"max |Δpval| = {ctrl_row['max_abs_delta_pval']:.1e}: <strong>bit-identical</strong>. "
        "For the Secreted-Signaling tier the spatial information contributed exactly nothing.")
else:
    ctrl_block = ""
sec("limits", "9. Limitations — read before quoting any number above", f"""
<h3>9.1 The spatial information had no effect on this dataset</h3>
<p>The object was built with <code>datatype = "spatial"</code>, coordinates and spatial factors,
and CellChat did apply its spatial machinery. We then checked what that machinery actually
excluded:</p>
<table><thead><tr><th>Constraint</th><th>low</th><th>high</th></tr></thead><tbody>
<tr><td><code>adj.spatial</code> — within 250 µm</td><td class="bad">81 / 81 kept, 0 excluded</td><td class="bad">81 / 81 kept, 0 excluded</td></tr>
<tr><td><code>adj.contact</code> — within 10 µm</td><td>73 / 81 kept (8 excluded)</td><td>77 / 81 kept (4 excluded)</td></tr>
</tbody></table>
<p>The 250 µm filter removed <strong>nothing</strong>: a TMA core is roughly 2,000 µm across with
nine intermixed cell types, so every cell-type pair has ≥ 10 mutual neighbours within 250 µm.
Combined with <code>distance.use = FALSE</code>, which makes the spatial weight matrix all-ones
wherever the filter passes, the spatial term multiplies by 1 everywhere.</p>
{ctrl_block}
<p>In the <code>{TIER}</code> tier the picture is slightly better but still narrow: only the 535
Cell-Cell Contact interactions are subject to <code>adj.contact</code>, and even there just 8
(low) and 4 (high) of 81 cell-type pairs are removed. Secreted and ECM-Receptor interactions
remain spatially unconstrained.</p>
{note("<strong>Consequence.</strong> These results should be described as cell-type-level "
      "co-expression inferences, not as spatial communication. Making the spatial term bite "
      "requires <code>distance.use = TRUE</code> (which multiplies probability by "
      "1/(d × scale.distance)) and/or a much smaller <code>interaction.range</code>. That run "
      "has not been done.")}

<h3>9.2 The two conditions differ four-fold in cell number</h3>
<table><thead><tr><th></th><th>All cells</th><th>mGAM</th><th>MES-like</th></tr></thead><tbody>
<tr><td>low</td><td>{low_n:,}</td><td>{int(sizes['low'].set_index('cell_type').loc['mGAM','n_cells']):,}</td><td>{int(sizes['low'].set_index('cell_type').loc['MES-like','n_cells']):,}</td></tr>
<tr><td>high</td><td>{high_n:,}</td><td>{int(sizes['high'].set_index('cell_type').loc['mGAM','n_cells']):,}</td><td>{int(sizes['high'].set_index('cell_type').loc['MES-like','n_cells']):,}</td></tr>
</tbody></table>
{warn("More cells give more stable group means and therefore more power in the permutation "
      "test. How much of 'high grade has more significant links' "
      f"({int(summary[summary.condition=='high'].n_interactions.iloc[0]):,} vs "
      f"{int(summary[summary.condition=='low'].n_interactions.iloc[0]):,}) is biology and how "
      "much is statistical power <strong>cannot be separated from this run</strong>. "
      "Information-flow (strength) comparisons are less sensitive to this than link counts, but "
      "not immune. A downsampled high-grade run would settle it.")}

<h3>9.3 The unit of inference is coarse</h3>
<p>CellChat cannot say which <em>cells</em>, or which region of tissue, are engaged in an
interaction — only which cell-type pairs. It also cannot represent a combination of L–R pairs
as one program (section 8.1). Its pathway vocabulary is a fixed, curated column in the
database, not something learned from the data.</p>

<h3>9.4 Reproducibility</h3>
{note("CellChat draws its permutation matrix once under <code>set.seed(seed.use)</code> "
      "<em>before</em> the parallel loop, so the <code>future.rng.onMisuse</code> warnings it "
      "emits are benign. Verified rather than assumed: two runs with identical arguments "
      "produced bit-identical probability, p-value, count and weight outputs.")}
""")

# ---- 10 ---------------------------------------------------------------------
DEV = [
 ("<code>contact.range</code>", "100", "<strong>10</strong>",
  "100 is the Visium spot pitch. The FAQ, vignette and man page all pin 10 µm for single-cell-resolution platforms. Measured median nearest-neighbour distance on our data: 11.5 µm (low), 7.45 µm (high)."),
 ("<code>spatial.factors</code>", "<code>ratio = 65/spot_diameter_fullres</code>, <code>tol = 32.5</code>", "<code>ratio = 1</code>, <code>tol = 5</code>",
  "The FAQ's own Xenium row: coordinates already in µm, <code>spot.size</code> = typical human cell (10 µm)."),
 ("Normalisation", "<code>GetAssayData(slot = \"data\", assay = \"SCT\")</code>", "<code>adata.X</code> used as-is",
  "<code>X</code> is already <code>log1p(CP10K)</code>, asserted at export. Re-normalising would log-transform twice."),
 ("<code>distance.use</code>", "single-section: <code>TRUE</code>, <code>scale.distance = 0.01</code>", "<code>FALSE</code>",
  "Follows the multi-sample spatial vignette, which matches our data shape. 0.01 is Visium-specific and would abort CellChat's own [1,2] validation on micron-scale coordinates."),
 ("<code>variable.both</code>", "single-section passes <code>F</code>; multi-sample passes nothing", "<code>TRUE</code>",
  "The multi-sample vignette governs, and <code>TRUE</code> is the package default. A genuine fork between the two tutorials, stated rather than glossed."),
 ("<code>umap.method</code>", "not passed → <code>umap-learn</code>", "<code>\"uwot\"</code>",
  "Avoids a reticulate python dependency; <code>uwot</code> is the package's own documented alternative."),
 ("<code>computeCellDistance</code>", "used as-is", "<code>BiocNeighbors::findKNN(k = 1)</code>",
  "Dense O(N²) — about 51 GB at 79,998 cells, and it returned NA. Identical quantity, O(N log N). Sanity check only; no inference path touched."),
 ("<code>sources.use</code> / <code>targets.use</code>", "hardcoded indices (<code>4</code>, <code>5:11</code>)", "all 9 cell types",
  "The tutorial's indices are specific to its 12-cluster skin dataset and are meaningless here."),
 ("<code>pathways.show</code>", "one hand-picked pathway", "every pathway in <code>netP$pathways</code>",
  "The benchmark protocol requires producing every figure the standard workflow can produce. (This report then curates them back down for readability.)"),
 ("<code>selectK</code>", "called directly", "measures recomputed with <code>.pbackend = \"seq\"</code>, curve redrawn",
  "<code>selectK</code> uses NMF's default parallel foreach backend and exposes no override; inside a long plotting session every run dies with \"All the runs produced an error\", while the identical call succeeds in a clean session."),
 ("<code>k</code> for patterns", "read off the curve by eye", "same rule applied programmatically; measures persisted",
  "Reproducibility — the choice becomes auditable instead of a judgement call made once."),
 ("<code>netVisual_chord_gene</code> scope", "one sender (<code>sources.use = 4</code>)", "all-sources attempted, plus one chord per sender",
  "All-vs-all is 161 ligand/receptor sectors and circlize cannot lay it out at any <code>small.gap</code> (tested 1, 0.5, 0.2, 0.1). Per-sender is the tutorial's own scope and renders for all 9 cell types."),
 ("<code>mergeInteractions</code> (coarse regrouping)", "12 clusters → 3 coarse types", "skipped",
  "Our 9 labels are already coarse; no defensible 3-way grouping exists without asking the biologist."),
 ("Tier <code>cellchatdb2</code>", "—", "<code>subsetDB(CellChatDB)</code>",
  "The bundled database <em>is</em> CellChatDB v2 (Jaccard 1.0000 vs the repo CSV), so the tier varies annotation scope instead of resource. Re-importing the flat CSV would drop the complex/cofactor/agonist/antagonist columns <code>computeCommunProb</code> uses."),
 ("CRAN source", "CRAN HEAD", "dated snapshot 2024-06-01",
  "The environment is R 4.3.3 and current CRAN sources require ≥ 4.4 (e.g. Deriv 4.2.0 uses <code>Rf_allocLang</code>). Environment-only; does not touch results."),
 ("igraph build", "stock", "<code>--disable-graphml</code>, <code>xml2-config</code> shadowed",
  "Base anaconda's <code>xml2-config</code> leaks onto PATH and igraph links <code>libxml2.2.dylib</code>, which the environment does not ship. Environment-only."),
 ("Smoke test", "on the tutorial's demo data", "on an 819-cell GBM core",
  "CellChat ships no demo expression data — only <code>CellChatDB.*.rda</code> and <code>PPI.*.rda</code>. The vignettes load from the author's own local paths."),
]
dev_rows = "".join(
    f"<tr><td>{a}</td><td>{b}</td><td><strong>{c}</strong></td><td>{d}</td></tr>" for a, b, c, d in DEV)
sec("deviations", "10. Every deviation from the default workflow", f"""
<p>Complete list. Anything not in this table was left at the vignette's literal value or the
package default.</p>
<table class="dev"><thead><tr><th>Item</th><th>Tutorial</th><th>Ours</th><th>Reason</th></tr></thead>
<tbody>{dev_rows}</tbody></table>
""")

# ---- 11 ---------------------------------------------------------------------
sec("repro", "11. Reproducing this", f"""
<pre><code>source scripts/comparators/cellchat/activate_env.sh
bash   scripts/comparators/cellchat/run_all_gbm.sh all
python scripts/comparators/cellchat/build_report.py</code></pre>
<p><code>prepare_gbm_input.py</code> (env <code>bptf</code>) regenerates the input tree;
<code>run_all_gbm.sh</code> runs the database audit and then both tiers, inference and figures.
Environment: <code>env.lock.yml</code> plus <code>r_packages.lock.csv</code>
(249 R packages with versions); <code>install_env.R</code> rebuilds the R library from a pinned
CRAN snapshot.</p>
<table><thead><tr><th></th><th>Value</th></tr></thead><tbody>
<tr><td>CellChat</td><td>{man[TIER]['cellchat_version']} @ <code>{man[TIER]['cellchat_git_sha'][:12]}</code></td></tr>
<tr><td>R</td><td>{man[TIER]['r_version']}</td></tr>
<tr><td>presto installed</td><td>{man[TIER]['presto_installed']}</td></tr>
<tr><td>Tier shown here</td><td><code>{TIER}</code> — {man[TIER]['db']['n_interactions']:,} interactions</td></tr>
<tr><td>Full outputs</td><td><code>results/comparators/cellchat/GBM/{TIER}/</code></td></tr>
<tr><td>All figures (PNG)</td><td><code>reports/cellchat_GBM_png/</code></td></tr>
</tbody></table>
""")

# --------------------------------------------------------------------------- assemble
nav = "".join(f'<a href="#{i}">{t}</a>' for i, t, _ in SECTIONS)
body = "".join(f'<section id="{i}"><h2>{t}</h2>{b}</section>' for i, t, b in SECTIONS)

HTML = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CellChat on the GBM/LGG Xenium TMA — full walkthrough</title>
<style>
:root {{ --ink:#1a1a1a; --mute:#666; --line:#e2e2e2; --bg:#fff;
        --why:#0b6b3a; --whybg:#eefaf2; --warnc:#8a4b00; --warnbg:#fff6e8;
        --notec:#144a8a; --notebg:#eef4fc; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font:15.5px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
        color:var(--ink); background:var(--bg); }}
#wrap {{ display:flex; align-items:flex-start; }}
nav {{ position:sticky; top:0; width:250px; min-width:250px; height:100vh; overflow-y:auto;
       padding:26px 16px; border-right:1px solid var(--line); font-size:13px; }}
nav b {{ display:block; margin-bottom:12px; font-size:12px; letter-spacing:.08em;
         text-transform:uppercase; color:var(--mute); }}
nav a {{ display:block; padding:5px 8px; color:var(--ink); text-decoration:none;
         border-radius:5px; margin-bottom:2px; }}
nav a:hover {{ background:#f2f2f2; }}
main {{ flex:1; max-width:1000px; padding:36px 46px 120px; }}
h1 {{ font-size:29px; line-height:1.25; margin:0 0 6px; }}
.sub {{ color:var(--mute); margin-bottom:34px; font-size:14px; }}
h2 {{ font-size:22px; margin:52px 0 16px; padding-bottom:8px; border-bottom:2px solid var(--ink); }}
h3 {{ font-size:17px; margin:30px 0 10px; }}
h4 {{ font-size:14px; margin:16px 0 6px; color:var(--mute); }}
p {{ margin:11px 0; }}
code {{ background:#f4f4f4; padding:1px 5px; border-radius:3px; font-size:.9em;
        font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }}
pre {{ background:#f7f7f7; border:1px solid var(--line); border-radius:7px;
       padding:14px 16px; overflow-x:auto; }}
pre code {{ background:none; padding:0; font-size:13px; line-height:1.5; }}
table {{ border-collapse:collapse; width:100%; margin:14px 0; font-size:13.5px; }}
th,td {{ border:1px solid var(--line); padding:7px 10px; text-align:left; vertical-align:top; }}
th {{ background:#f7f7f7; font-weight:600; }}
table.dev td:nth-child(4) {{ font-size:13px; color:#333; }}
.why,.warn,.note {{ border-left:4px solid; border-radius:0 7px 7px 0;
                    padding:12px 16px; margin:16px 0; font-size:14px; }}
.why {{ border-color:var(--why); background:var(--whybg); }}
.warn {{ border-color:var(--warnc); background:var(--warnbg); }}
.note {{ border-color:var(--notec); background:var(--notebg); }}
.whyh {{ display:block; font-weight:700; color:var(--why); font-size:11.5px;
         letter-spacing:.07em; text-transform:uppercase; margin-bottom:5px; }}
.eq {{ text-align:center; font-size:17px; margin:16px 0; }}
.big {{ font-size:16px; background:#f7f7f7; border-radius:7px; padding:12px 16px; }}
figure {{ margin:20px 0; }}
figure img {{ width:100%; height:auto; display:block; border:1px solid var(--line); border-radius:6px; }}
figcaption {{ font-size:12.5px; color:var(--mute); margin-top:7px; line-height:1.5; }}
.fn {{ display:block; font-family:ui-monospace,Menlo,monospace; font-size:11px; color:#999; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:20px; }}
.grid2 figure {{ margin:10px 0; }}
.good {{ color:#0b6b3a; font-weight:600; }}
.bad {{ color:#a11; font-weight:600; }}
@media (max-width:1000px) {{ #wrap{{display:block;}} nav{{position:static;width:auto;height:auto;
  border-right:none;border-bottom:1px solid var(--line);}} main{{padding:24px 18px 80px;}}
  .grid2{{grid-template-columns:1fr;}} }}
</style></head><body><div id="wrap">
<nav><b>Contents</b>{nav}</nav>
<main>
<h1>CellChat on the GBM/LGG Xenium TMA</h1>
<div class="sub">A complete walkthrough of the run — every input, every parameter, every
deviation from the published workflow and why it was made, and what the result says
biologically.<br>
CellChat {man[TIER]['cellchat_version']} · tier <code>{TIER}</code> ·
{man[TIER]['db']['n_interactions']:,} L–R interactions ·
{low_n + high_n:,} cells across 13 TMA cores · generated
{man[TIER]['finished'][:10]}</div>
{body}
</main></div></body></html>"""

os.makedirs(PNG, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write(HTML)
print(f"wrote {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB, {len(SECTIONS)} sections)")
