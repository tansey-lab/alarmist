# stLearn — tutorial call contract

Source: `stlearn.readthedocs.io/en/latest/tutorials/cell_cell_interaction_xenium.html`
(**the Xenium-specific CCI tutorial** — the right one for our data; the generic
`cell_cell_interaction.html` targets Visium). Local mirror at
`/Users/jiayifan/tansey_lab/stLearn/`. Source read from `_modules/stlearn/tl/cci/analysis.html`.

Tutorial dataset: 10x `Xenium_FFPE_Human_Breast_Cancer_Rep1`, 167,780 cells → **164,000**
after `min_counts=10`, 313 genes, extent **7,521 × 5,471 µm**.

## One row per tutorial call

| # | Call (tutorial values) | Produces | Ours |
|---|---|---|---|
| 1 | `st.pp.filter_genes(adata, min_counts=10)` | gene QC | same |
| 2 | `st.pp.filter_cells(adata, min_counts=10)` | cell QC | same |
| 3 | `adata.raw = adata` | store raw for PSTS | same |
| 4 | `st.em.run_pca(n_comps=50, random_state=0)`, `st.pp.neighbors(n_neighbors=25, use_rep='X_pca')`, `st.tl.clustering.leiden(resolution=1.05)` | `obs['leiden']` cell labels | **skipped** — we have real annotations (`obs['cell_type']`, 9 types). Leiden is only the tutorial's stand-in for cell types. |
| 5 | `st.pp.normalize_total(adata)` | library-size normalisation **only** | same |
| 6 | `st.tl.cci.grid(adata, n_row=125, n_col=125, use_label='leiden', n_cpus=n)` | gridded AnnData; `use_label` → per-spot deconvolution in `uns` | `n_row=321, n_col=146`, `use_label='cell_type'` (see DEVIATIONS) |
| 7 | `st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')` | 2,293 LR pairs, `L_R` strings | **CellChatDB v2**, complexes dropped (see DEVIATIONS) |
| 8 | `st.tl.cci.run(grid, lrs, min_spots=20, distance=250, n_pairs=1000, n_cpus, random_state=0)` | `obsm['lr_scores','p_vals','p_adjs','-log10(p_adjs)','lr_sig_scores']`, `uns['lr_summary']` | same, but `n_pairs=10000` (tutorial says "recommend ~10,000"; 1000 is "low as example") |
| 9 | `st.tl.cci.adj_pvals(grid, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')` (CELL 62) | re-adjusts p-values and re-ranks `lr_summary` | **verified no-op** — `run()` already applies these exact defaults (`analysis.py:382`). Called in `plot_stlearn_tutorial.py` with a runtime assertion; not called during the run itself. Its only effect is an unstable re-sort that permutes LRs tied on `n_spots_sig`. |
| 10 | `st.tl.cci.run_cci(grid, 'leiden', min_spots=2, spot_mixtures=True, cell_prop_cutoff=0.1, sig_spots=True, n_perms=100, random_state=0, n_cpus)` | `uns['lr_cci_<label>']`, `uns['per_lr_cci_<label>']`, `uns['per_lr_cci_pvals_<label>']`, `uns['per_lr_cci_raw_<label>']`, `uns['lr_cci_raw_<label>']`, **+3 columns on `uns['lr_summary']`** | `'cell_type'`, `n_perms=1000` (tutorial says "recommend ~1000"; 100 is the example value) |

**Crucially — step 5 is the whole normalisation.** The tutorial is explicit: *"No log1p or
shrinking to make genes of similar expression range. In our case, for calling hotspots, we want
genes to be more separate, since we select background genes with similar expression levels to
detect hotspots."* Log-transforming would break the background-gene matching that the
permutation null depends on. Our GBM `X` is log-normalised, so we **must** load from
`layers['counts']`.

## Plots — **the Xenium vignette's set only**

These are the *only* `st.pl.*` calls the Xenium vignette makes, in order, with its own argument
values and figure compositions. `plot_stlearn_tutorial.py` reproduces this table one-to-one.

| CELL | Call (tutorial values) | Shows |
|---|---|---|
| 23 | `st.pl.cluster_plot(adata, use_label, size=1, show_image=False, bbox_to_anchor=(1.2,1))` | single-cell labels |
| 40 | **2-panel**: `cluster_plot(grid, size=10, ...)` \| `cluster_plot(adata, ...)` | grid dominant spots vs single-cell labels — the gridding sanity check |
| 42 | **3-panel × group**: `feat_plot(grid, feature='group', vmax=1, show_color_bar=False)` \| `cluster_plot(grid, list_clusters=[g])` \| `cluster_plot(adata, list_clusters=[g])` | per-spot proportion, grid subset, cell subset. Loop is over **3 clusters chosen by position** (`groups[6/10/11]`) — the `idc`/`dcis`/`stroma` names are local variables, not a rule |
| 44 | **2-panel**: `gene_plot(grid, 'CXCL12')` \| `gene_plot(adata, 'CXCL12', vmax=80)` | gene expression, gridded vs original |
| 58 | `lr_summary(grid, n_top=500)`; `lr_summary(grid, n_top=50, figsize=(10,3))` | LR ranking by number of significant spots |
| 67 | **3-stat panel** `lr_result_plot(grid, use_result=stat, use_lr=best_lr)` for `stat` ∈ `lr_scores`, `-log10(p_adjs)`, `lr_sig_scores`, **`best_lr` = top-1 only** | per-LR spatial map ← *this is where the requested LRIs go* |
| 79 | `cci_check(grid, label, figsize=(16,5))` | diagnostic: interaction vs cell-type-frequency dependence (should be ~none) |
| 81 | `pos_1 = ccinet_plot(grid, label, return_pos=True, min_counts=30)`, then **top-2** with `min_counts=2, figsize=(10,7.5), pos=pos_1` | cell-type interaction network; `pos=pos_1` gives every per-LR network the **same node layout** |
| 83 | `lr_chord_plot(grid, label)`, then **top-2** | chord diagram of cell-type interactions |

### NOT in the Xenium vignette

- **Generic (Visium) vignette only**, so out of scope here: `lr_diagnostics`, `lr_n_spots`,
  `cci_map`, `lr_cci_map`, `lr_plot`, `lr_go`. That vignette uses `st.datasets.visium_sge` and
  never calls `st.tl.cci.grid()`.
- **In neither vignette**: `het_plot`, `grid_plot`, `deconvolution_plot`.

Both groups exist in the earlier `plots/` and `plots_full/` trees; see `DEVIATIONS.md`
§"Plot provenance" for the full matrix and why they were kept rather than deleted.

## Facts confirmed from source, not assumed

- `grid()`: `np.histogram2d(xs, ys, bins=[n_col, n_row])` → **`n_col` bins x, `n_row` bins y**.
  Coordinates are cast to `int`, and **spots with 0 cells are dropped**, so "cells per spot"
  must be computed over *occupied* spots.
- `grid()` docstring: *"intended use is for single cell spatial data, not Visium"* — gridding is
  the authors' intended Xenium path, not a workaround we invented.
- **The neighbourhood is `distance=250` (µm, physical, via cKDTree), independent of the grid.**
  Grid resolution controls aggregation granularity and how many spots fall inside that 250 µm
  radius — it is *not* the neighbourhood parameter.
- The tutorial's `n_ = 125` is fed to both axes and its own markdown calls it a
  resolution/compute trade-off (*"The higher resolution, the better this represents the single
  cell data but the longer the computation takes"*) — an author-declared compute knob, **not** a
  tuned biological parameter.
