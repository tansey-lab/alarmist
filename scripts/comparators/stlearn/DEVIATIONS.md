# stLearn — deviations from the tutorial

The tutorial is a **single contiguous breast section**; ours is a **TMA of 13 disconnected
cores**. Every deviation below follows from that, from the panel, or from an author-declared
"example value, recommend higher".

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| cell labels | `leiden` clustering (steps 4) | `obs['cell_type']` (9 types), Leiden skipped | Leiden is the tutorial's stand-in because its demo has no annotation. We have real annotations, and `use_label` only needs *a* per-cell label. Running Leiden anyway would substitute a worse label for a better one. |
| expression matrix | Xenium raw counts from `cell_feature_matrix.h5` | `layers['counts']` of our h5ad | Our `X` is log-normalised. The tutorial explicitly forbids log1p before `run()` because the permutation null matches background genes by expression level. Using `X` would silently break the null. |
| grid | `n_row=125, n_col=125` | **`n_row=321, n_col=146`** (51.3 × 51.3 µm) | The tutorial's 125 gives **60.2 × 43.8 µm** spots on *its* 7,521 × 5,471 µm section — rectangular, and specific to that extent. Our TMA is 7,484 × 16,483 µm (2.20:1), so copying 125×125 would give **59.9 × 131.9 µm**, a 2.2:1 elongated spot. We instead preserve the tutorial's **spot area** (2,637 µm² → 51.4 µm equivalent square). Measured outcome: **12,562 occupied spots** vs the tutorial's **14,364** — same order, our TMA being mostly empty slide between cores. **The resulting ~51 µm is *not* chosen to match ALARMIST's 50 µm patch size — that is a coincidence of the arithmetic, and no comparator parameter was tuned toward ALARMIST.** |
| cells per spot | 11.4 per occupied spot | **7.98 per occupied spot** (measured; `run_manifest.json`) | Consequence, not a choice: GBM in-core density is lower than the tutorial's breast section. Matching 11.4 would have required larger spots, i.e. changing the method's spatial resolution to compensate for our tissue being sparser. We preserved resolution and let cells/spot fall where it falls. |
| LR database | `connectomeDB2020_lit`, 2,293 pairs | **CellChatDB v2, complex rows dropped**: 1,371 pairs, 527 on the GBM 5K panel | The `cellchatdb2` tier requires CellChatDB v2. stLearn's LR format is a flat `L_R` gene-pair string with no complex support, so the **1,859** rows (57.5%) with a multi-subunit ligand or receptor cannot be represented. Dropping them is a **limitation of the method**, not a modelling choice — expanding complexes combinatorially would invent LR pairs CellChat does not assert. Consequence to report: stLearn tests less than half the panel-available interactions that complex-aware methods do (527 vs CytoSignal's 1,088). |
| `n_pairs` | 1,000 | **10,000** | The tutorial's own comment: `# Number of random pairs to generate; low as example, recommend ~10,000`. Following the recommendation, not the example. |
| `n_perms` (`run_cci`) | 100 | **1,000** | Same: `# Permutations of cell information to get background, recommend ~1000`. |
| `distance` | 250 | 250 (unchanged) | Physical µm and grid-independent. Kept at the authors' default per the "do not harmonise neighbourhoods" rule. |
| whole-slide vs per-core | single section | **whole-slide, one joint run** | `distance=250` could in principle link spots across TMA cores (different patients). Measured: **0 cross-core neighbour pairs at 200 µm; 39 pairs / 43 cells (0.04%) at 250 µm** — the cores' bounding boxes come within 150 µm but the cores are round, so cells never actually get that close. Per-core running would also break the smallest cores (core 2 has 819 cells → ~15 occupied spots, below `min_spots=20`). Documented rather than worked around. |
| per-LR plots | top LR only | top LRs **plus** GRN→SORT1 and ANXA1→FPR1 | Standing request; the two go in a separate `requested/` directory so they are never confused with stLearn's own ranking. |
| `st.tl.cci.adj_pvals` | called (CELL 62) with `correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh'` | **not called during the run; called in `plot_stlearn_tutorial.py` and asserted to be a no-op** | `st.tl.cci.run()` already applies these exact settings internally (`run()` defaults `adj_method='fdr_bh'`, `pval_adj_cutoff=0.05`; `permutation.py:172-178` does the per-spot MHT correction), and `adj_pvals`' own docstring says so: *"Default settings of this function are already run in st.tl.cci.run"* (`analysis.py:382`). **Verified empirically, not assumed**: after calling it on the saved grid, all five `obsm` matrices are element-wise identical once realigned, and all six `lr_summary` columns are identical. The only effect is that `np.argsort` (unstable) re-permutes **179 / 526** LRs *within groups tied on `n_spots_sig`*; the top-3 are unchanged. So the on-disk results, produced without the call, are the tutorial's results. |
| `gene_plot` gene | `CXCL12`, a prior biological pick for the tutorial's breast section | **ligand + receptor of the top-ranked LR pair** (`C3`, `C3AR1`) | The tutorial gives no rule, just a gene it cared about. Rather than invent a biological pick we apply the tutorial's own `lr_summary.index[0]` idiom to gene identity, and state the rule. |
| `gene_plot` single-cell panel | `vmax=80` | **omitted** | 80 is a display clip tuned to the tutorial's own count distribution; transplanting the literal number onto a different panel and chemistry would misrepresent our data. The grid panel has no `vmax` in the tutorial either. |
| `feat_plot` group selection | 3 clusters chosen by **position** — `groups[6]`, `groups[10]`, `groups[11]` of `grid.obs['leiden'].cat.categories` (the `idc`/`dcis`/`stroma` names are local variables, not a selection rule) | **all 9 cell types** | There is no positional equivalent for an annotated 9-type vocabulary, and any fixed slice would be arbitrary. Same call, same arguments, complete iteration — this adds no analysis, no parameter and no number, and it removes the possibility of a favourable subset. The tutorial itself leaves ≥9 of its ≥12 clusters unplotted. |
| `feat_plot` arguments | `vmax=1, show_color_bar=False` | same | **This was wrong in the earlier `plots_full/` pass** (`show_color_bar=True`, no `vmax`), which let each cell type's colour scale stretch to its own maximum so the proportion panels were not comparable to each other. `vmax=1` is what makes "proportion, max = 1" mean the same thing in every panel. |
| `ccinet_plot` node layout | `pos_1 = ccinet_plot(..., return_pos=True)`, then `pos=pos_1` on every per-LR network | same | `run_stlearn.py` omitted `return_pos`/`pos`, so each per-LR network was laid out independently and the networks could not be compared by eye. Restored in `plot_stlearn_tutorial.py`. |
| `figsize` | `(20, 5)` for the 2- and 3-panel figures, `(20, 8)` for `gene_plot` / `lr_result_plot` | **each panel sized to the tissue aspect** (4.5 × 10.0 in, from the measured 2.21:1) | The **only** tutorial argument not copied literally, and the one place our TMA geometry forces a display change. Those figsizes are for the tutorial's 1.37:1 breast section; combined with the tutorial's *own* `set_aspect('equal')` they letterbox our 2.20:1 TMA into roughly a tenth of the canvas and drop the legend on top of the upper cores. `figsize` is pure canvas — it changes no datum, no statistic, no selection, no colour scale — so every other argument stays verbatim and only the box is matched to the data. Computed from `grid.obsm['spatial']` at runtime, not hardcoded. |
| plot scope | the Xenium vignette's own set | **`plots_tutorial/` = exactly that set; `plots/` and `plots_full/` are broader and are kept as-is** | See "Plot provenance" below. |
| `run_cci` dtype | works as published | **`obs['cell_type']` coerced to object dtype around the call**, categorical restored afterwards | **Genuine version gap, not a parameter choice.** `run_cci` does `adata.obs[label].values.astype(str)` and feeds the result to an `@njit` kernel. Under **pandas 3.0.5 / numpy 2.4.6** that returns a `StringDtype`/object array; numba 0.66 cannot type it and raises `TypingError: non-precise type array(pyobject, 1d, C)` at `het.py:227`, hard-blocking the entire cell-type CCI half of the workflow. Object-dtype Python strings make the *same expression* return a plain `<U` array, which numba accepts. **Labels are byte-identical — only the dtype changes**, asserted at runtime. `run_cci` never uses `.cat`; the plot functions do, so the categorical is restored immediately after. The alternative fix is pinning `pandas<3`, which we did not do because the coercion is provably label-preserving and keeps the env on current scanpy. |
| cell-type ↔ column mapping | not an issue for the tutorial's `leiden` labels | **runtime guard added** | `get_data_for_counting` binds cell types to deconvolution columns by **substring**, taking the first hit: `[ct in col for col in cols]`. Our labels include both `mGAM` and `non-mGAM`, and `'mGAM' in 'non-mGAM'` is `True`, so the binding is order-dependent — a silent mis-assignment would corrupt every interaction reported for that type. **Verified correct here** (alphabetical ordering puts `mGAM` first, so it wins its own column), but verified rather than assumed, and the runner now aborts if any label mis-binds. This is a latent stLearn bug that would bite any dataset with nested cell-type names. |
| control probes | tutorial's 313 genes are real panel only | **21 `Intergenic_Region_*` genes dropped** before `normalize_total` | `st.tl.cci.run()` hard-errors on any gene name containing `_`, because `_` is its L/R separator in the `"LIGAND_RECEPTOR"` encoding ("Recommend to rename adata.var_names or remove these genes"). On our panel every offender is an `Intergenic_Region_*` Xenium **genomic control probe** — not a real gene, absent from every LR database, and normally excluded upstream (10x keeps controls in separate feature types, which is why the tutorial's input has none). Dropped *before* `normalize_total` so library size reflects the real panel, as in the tutorial. 5,119 → 5,098 genes; **zero biology lost**. Note the earlier CytoSignal GBM run kept them, but CytoSignal scores only DB genes so they were inert there. |
| object construction | `st.read_xenium(...)` from a 10x bundle | plain h5ad + **hand-built `uns['spatial']`** | `st.convert_scanpy()` and `st.tl.cci.grid()` both index `adata.uns['spatial']` (grid copies it verbatim, `analysis.py:172`) — a plain h5ad has no such key, so `grid()` raises `KeyError: 'spatial'`. We replicate exactly what `read_xenium` builds: `{images:{hires:arr}, use_quality:'hires', scalefactors:{tissue_hires_scalef:1, spot_diameter_fullres:15}}`. Two sub-points: (a) `read_xenium` *itself* creates a **blank** placeholder image when no image file is given, so a placeholder is the reader's own behaviour, not our invention — we just make it 1×1 instead of `(1.1·max_coord)²` = 18,299² RGBA ≈ **1.34 GB**, since every plot passes `show_image=False` and the array is read only on the `_add_image()` path; (b) `scalef=1` matches the tutorial's `read_xenium(scale=1)` for micron coordinates, which makes `convert_scanpy`'s `image_coor = obsm['spatial'] * scale` the identity — asserted at runtime. stLearn's own source notes the scalefactor "scales the IMAGE to match the spot spatial coordinates, not the spots", so it cannot displace a data point. |

## Plot provenance — which vignette each figure comes from

Verified 2026-08-04 by grepping the two official vignettes (local mirror
`/Users/jiayifan/tansey_lab/stLearn/stlearn.readthedocs.io/en/latest/tutorials/`), then
adversarially re-checked. The **Xenium** vignette is the correct one for our data; the generic
one targets Visium (`st.datasets.visium_sge`) and does not use `st.tl.cci.grid()` at all.

| `st.pl.*` call | Xenium vignette | generic Visium vignette |
|---|---|---|
| `cluster_plot` | ✅ ×5 | ✅ |
| `feat_plot` | ✅ (CELL 42, `vmax=1`) | ❌ |
| `gene_plot` | ✅ (CELL 44) | ❌ |
| `lr_summary` | ✅ ×2 | ✅ |
| `lr_result_plot` | ✅ (CELL 67, top-1 only) | ✅ |
| `cci_check` | ✅ | ✅ |
| `ccinet_plot` | ✅ (CELL 81, shared `pos`) | ✅ |
| `lr_chord_plot` | ✅ (CELL 83) | ✅ |
| `lr_diagnostics` | ❌ | ✅ |
| `lr_n_spots` | ❌ | ✅ |
| `cci_map` | ❌ | ✅ |
| `lr_cci_map` | ❌ | ✅ |
| `lr_plot` | ❌ | ✅ ×8 |
| `lr_go` | ❌ | ✅ |
| `het_plot` | ❌ | ❌ |
| `grid_plot` | ❌ | ❌ |
| `deconvolution_plot` | ❌ | ❌ |

### The three plot directories

| Directory | Produced by | Scope |
|---|---|---|
| `plots_tutorial/` | `plot_stlearn_tutorial.py` | **The Xenium vignette's set, call-for-call — 23 PNGs + 6 in `requested/`.** This is the figure set to cite as "the authors' default workflow". |
| `plots/` | `run_stlearn.py` | The original pass: the Xenium set *plus* `lr_diagnostics`, `lr_n_spots`, `cci_map`, `lr_cci_map` (generic-vignette calls), and top-6 rather than top-1/top-2 where the tutorial subsets. 29 PNGs. |
| `plots_full/` | `plot_stlearn_full.py` | A completeness pass that added `lr_plot` (generic vignette) and `het_plot`, `grid_plot`, `deconvolution_plot` (**in neither vignette**). 26 PNGs. Retained for reference; `deconvolution_plot` never produced output — see below. |

Nothing was deleted; `plots_tutorial/` was added alongside. Note that `plot_stlearn_full.py`
also **rewrites `data/grid.h5ad`** (adding `obsm['het']`, `obs['prop_*']`, `uns['deconvolution']`),
which is why that file's mtime is later than the run. `plot_stlearn_tutorial.py` and
`export_stlearn_quant.py` both open it read-only.

### Known defects in `plot_stlearn_full.py` (left in place, not fixed)

Since that script's output is out of scope for the tutorial-fidelity deliverable, its bugs are
recorded rather than repaired:

1. **`deconvolution_plot` never ran.** The script writes `grid.uns["deconvolution"]`, but
   `stlearn/pl/deconvolution_plot.py:86` reads `adata.obsm["deconvolution"]` — all 4 calls died
   with `KeyError: 'deconvolution'`, caught by its `guard()` and visible only in the log.
   (`stlearn/adds/add_deconvolution.py:35` shows the intended field.) The call is not in either
   vignette, so under the fidelity rule the correct resolution is removal, not repair.
2. **`feat_plot` covered 4 of 9 cell types** (`[:4]` on alphabetically-ordered categories),
   dropping NPC-like, OPC-like, Vascular, mGAM and non-mGAM — mGAM being the one of interest.
   Superseded by `plots_tutorial/cell42_composition_*.png`, which covers all 9.
3. **`feat_plot` colour scales were not comparable** (`show_color_bar=True`, no `vmax`); the
   tutorial's `vmax=1` is restored in `plots_tutorial/`.
4. `gene_plot` covered only `top[:2]` + requested, so GJA1 / APP / SORL1 were never drawn.

### `lr_go` — genuinely unavailable, not skipped by choice

`st.tl.cci.run_lr_go(adata, r_path, ...)` requires a path to an R installation with
`clusterProfiler`, `org.Hs.eg.db` and `org.Mm.eg.db`. Checked 2026-08-04: none of the three
comp envs carrying R (`comp-cellchat`, `comp-cytosignal`, `comp-niches`) has `clusterProfiler`,
and `comp-stlearn` is Python-only. It is also a generic-vignette call, so it is out of scope for
`plots_tutorial/` regardless.
