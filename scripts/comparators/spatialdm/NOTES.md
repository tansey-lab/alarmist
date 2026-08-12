# SpatialDM — tutorial call contract

Source: `/Users/jiayifan/tansey_lab/SpatialDM/tutorial/melanoma.ipynb` (main),
`differential_test_intestine.ipynb` (multi-sample), `misc_notes.ipynb`.
Package source read directly: `spatialdm/main.py`, `spatialdm/stats.py`.

**Status: PREPARED, NOT RUN.** Written while the stLearn GBM run was in flight. The run waits
on stLearn sign-off, and on the one open parameter below.

## One row per tutorial call

| # | Call (tutorial values) | Produces |
|---|---|---|
| 1 | `sdm.weight_matrix(adata, l=1.2, cutoff=0.2, single_cell=False)` | `obsp['weight']` (secreted, RBF), `obsp['nearest_neighbors']` (adjacent) |
| 2 | `sdm.extract_lr(adata, 'human', min_cell=3)` | `uns['ligand']`, `uns['receptor']`, `uns['geneInter']` — valid LR pairs |
| 3 | `sdm.spatialdm_global(adata, method='z-score', nproc=1)` | global bivariate Moran's R per LR pair |
| 4 | `sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.1)` | `uns['global_res']` with `.selected` |
| 5 | `sdm.spatialdm_local(adata, method='z-score', nproc=1)` | `uns['local_stat']['local_I']`, `local_z_p` / `local_perm_p` |
| 6 | `sdm.sig_spots(adata, method='z-score', fdr=False, threshold=0.1)` | `uns['selected_spots']`, `uns['local_stat']['n_spots']` |

Downstream (optional in the tutorial): SparseAEH `MixedGaussian` clustering of local spot
patterns → `sdm.compute_pathway(adata, dic=...)` → pathway enrichment.

## Plots the standard workflow produces

| Call | Shows |
|---|---|
| `pl.global_plot(adata, pairs=[...])` | global selection volcano; named pairs highlighted |
| `pl.plot_pairs(adata, [...], marker='s')` | **per-LR spatial maps** ← where GRN_SORT1 / ANXA1_FPR1 go |
| `plot_clusters(gaussian)` | spatial patterns from SparseAEH clustering |
| `pl.dot_path(pathway_res)` | pathway enrichment dot plot per pattern |
| `pl.chord_celltype(adata, pairs=[...])` | chord diagram of aggregated cell types |

## Facts confirmed from source, not assumed

- **Input expects log-transformed `X` with raw in `.raw`** (opposite of stLearn, which forbids
  log1p). The melanoma tutorial states: log-transformed in `adata.X`, raw in `adata.raw`.
  Our GBM h5ad already matches this natively — no `layers['counts']` swap needed.
- **Default DB is CellChatDB v1**, **1,939 human interactions + 157 complexes**, three signaling
  categories only (Secreted 1,199 / ECM-Receptor 421 / Cell-Cell Contact 319) — *not* v2, so
  unlike CellChat the `default` and `cellchatdb2` tiers genuinely differ here. It ships inside
  the wheel at `spatialdm/datasets/LR_data/human-{interaction,complex}_input_CellChatDB.csv{.gz,}`
  (the repo checkout also has it at `0_CellChatDB/{1_human,2_mouse}/`).
- **`extract_lr(datahost=…)` is named backwards.** `datahost='package'` reads those bundled
  files; `datahost='builtin'` — the function's own **default**, hence what the tutorial gets —
  **downloads from figshare on every call**. Use `--datahost package`: same DB, offline, pinned,
  and not re-fetched once per split.
- **`datahost='package'` imports `pkg_resources`, which setuptools ≥ 81 no longer ships.** The env
  had setuptools 83 → `ModuleNotFoundError: No module named 'pkg_resources'`. Pinned to
  `setuptools<81` (80.10.2); `env.lock.yml` now records it.
- **Chord plots shell out to `geckodriver`/`firefox` via bokeh, and the env's `bin/` is not on
  `PATH`** when the interpreter is invoked by absolute path (which is how this repo must run —
  conda activate is broken). Both binaries are installed in `comp-spatialdm/bin/`. Both runners
  now prepend `dirname(sys.executable)` to `PATH` at import.
- **`chord_LR` wants one `obs` COLUMN per cell type**, not category names — `adata.obs.loc[:,
  sender]`. The intestine tutorial's obs holds deconvolution proportions. Copy the one-hot
  `obsm['celltypes']` into `obs` first. `senders`/`receivers` are **zipped, not crossed**.
- **Complexes are supported natively** via `complex_input` (`subunit_1..subunit_4`), so the
  `cellchatdb2` tier loses nothing — unlike stLearn, which can only take simple `L_R` pairs and
  drops 57.5% of CellChatDB v2. Building the v2 tier means emitting our CSV in these two files'
  format.
- Kernel: `w(d) = exp(-d² / (2·l²))`, then weights below `cutoff` are zeroed
  (`stats.py:_Euclidean_to_RBF`). So the **effective radius is `R = l·sqrt(-2·ln cutoff)`**.
- `single_cell=True` zeroes the diagonal ("no within-spot communication"). Xenium is
  single-cell, so this must be `True` — the tutorial's `False` is for 200 µm ST spots.

### TRAP: `eff_dist` is a *squared* distance despite its name

`stats.py:71` sets `l = np.sqrt(-eff_dist / (2*np.log(cutoff)))`. Solving `w(R) = cutoff` for
the actual kernel gives `l = R / sqrt(-2·ln cutoff)`, i.e. the code omits the square on
`eff_dist`. Verified against the tutorial's own numbers: `l=1.2, cutoff=0.2` →
`R = 1.2·sqrt(-2·ln 0.2) = 2.153` units, and `2.153² = 4.635 = eff_dist`. So passing a physical
radius to `eff_dist` silently yields a kernel `sqrt(R)`-wide. **Set `l` directly from
`l = R / sqrt(-2·ln cutoff)`; do not pass `eff_dist` a distance.**

## SETTLED — kernel scale, measured not assumed

Loaded `dataset.melanoma()` (293 spots × 16,148 genes) and measured it:

| quantity | measured |
|---|---|
| coordinate extent | x 3–28, y 7–27 — **integers** |
| nearest-neighbour spacing | **exactly 1.0000** → coordinates are array indices, 1 unit = 1 spot pitch |
| platform pitch | 200 µm center-to-center (Thrane 2018 ST; the tutorial text says so too) |
| tutorial effective radius | `l·sqrt(-2·ln cutoff)` = `1.2·1.794` = **2.153 units = 431 µm** |
| **→ `l` for micron coordinates** | **240.0** at `cutoff = 0.2` |

`X.max() = 9.98` also confirms the tutorial's input is log-transformed, as documented.

### SUPERSEDED — `l = 240` was the wrong conclusion

The above measurement is correct but the inference from it was not. There is a **second**
calibration point: the intestine tutorial states *"spot-spot distance of 100 micrometers, `l`
will be set to 75 here… The parameters here should be determined to match the context of CCC."*

| dataset | unit pitch | `l` | effective radius | radius / pitch |
|---|---|---|---|---|
| melanoma | 200 µm | 240 µm (1.2 units) | 431 µm | 2.15 |
| intestine | 100 µm | 75 µm | 135 µm | 1.35 |

The two disagree **3× in absolute radius**, so neither absolute value is an author claim about
signalling range — both follow from data resolution. Taking melanoma's 431 µm as "the authors'
physical range" was unjustified.

### DECIDED — `l = 75`, `cutoff = 0.2`, `n_neighbors` per split, `n_nearest_neighbors = 6`, `single_cell = True`

Effective radius **134.6 µm**, i.e. the intestine tutorial's absolute radius preserved. It is
the only place the authors state their reasoning, and it sits in a plausible secreted-signalling
range. Deliberately **not** harmonised with CytoSignal's 200 µm or stLearn's 250 µm.

**`n_neighbors` and `n_nearest_neighbors` are independent** — `rbfweight` derives
`n_neighbors = n_neighbor_layers * 31` **only when `n_neighbors is None`**. `n_nearest_neighbors`
is passed as `n_neighbor_layers` and used *directly* as the size of a separate small graph
(`nnbrs0` → `obsp['nearest_neighbors']`, adjacent signalling). Passing `n_neighbors` explicitly
leaves it untouched.

**The kNN cap must not bind.** In both tutorials it never does (≈14.6 units inside melanoma's
radius, ≈5.7 inside intestine's, against a 186 default) — `cutoff` is the truncator and the cap
is a loose ceiling. Letting it bind, or reverse-engineering `l` from it, would invert how the
method operates. On GBM the required size is per-core (see DEVIATIONS.md).

**Once the cap does not bind, `W` is numerically identical for any larger `n_neighbors`** — the
extra neighbours are all zeroed by `cutoff`. Only stored `nnz` differs, because `rbfweight` sets
sub-cutoff entries to 0 without calling `eliminate_zeros()`, so stored `nnz = n_cells ×
n_neighbors` regardless. The runner calls `eliminate_zeros()` immediately after `weight_matrix`;
nothing downstream touches `.nonzero()`/`.nnz` on `W` (all uses are `.multiply`/`.sum`/`@`).

## Differential (multi-sample) call contract — `differential_test_intestine.ipynb`

Runner: `run_diff_spatialdm.py`. Consumes the persisted `<split>/data/spatialdm.h5ad` objects;
**recomputes nothing** — every z, p and `selected` flag comes off disk.

| # | Call (tutorial values) | Produces |
|---|---|---|
| 1 | `concat_obj(samples, names, 'human', 'z-score', fdr=False)` | `cdata.uns['p_df'/'zscore_df'/'tf_df']` (union-pairs × samples), `uns['ligand'/'receptor'/'geneInter']` |
| 2 | `sns.clustermap(1 - concat.uns['p_df'])` | overview of every pair × sample |
| 3 | `differential_test(cdata, subset, conditions)` | `uns['p_val']`, `uns['diff']`, `uns['diff_fdr']` — OLS likelihood-ratio of z on one covariate |
| 4 | `group_differential_pairs(cdata, c1, c2, 0.7, 0.3, fdr_co=0.1)` | `uns['<c>_specific']`, `uns['<c>_only']`, `uns['q1'/'q2'/'fdr_co']` |
| 5 | `pl.differential_dendrogram(cdata)` | clustermap of `1-p` over pairs with **raw** `p_val < 0.1` |
| 6 | `pl.differential_volcano(cdata, legend=[...])` (+ a second call with `pairs=[...]`) | `diff` vs `-log10(diff_fdr)` |
| 7 | `pl.dot_path(...)` on each `_specific` set | pathway enrichment of the condition-specific pairs |

### Facts confirmed from source, not assumed

- **`concat_db` hard-codes CellChatDB v1** (figshare) and then `geneInter.loc[ligand.index]`. Our
  injected v2 names are not in the v1 index → `KeyError`. Monkeypatched; see DEVIATIONS.md.
- **`dot_path`'s installed signature is `dot_path(pathway_res, cut_off=1, ...)`.** The tutorial's
  `pl.dot_path(concat, 'adult_specific', cut_off=2, ...)` would bind the string to `cut_off` —
  a version gap. Call `compute_pathway(dic=...)` first and pass the frame.
- **`compute_pathway` needs `geneInter.interaction_name` as a COLUMN.** `run_spatialdm.py` set it,
  but anndata wrote it out *as the index* (`uns/geneInter` attrs: `_index: interaction_name`) and
  dropped the column, so it is absent after `read_spatialdm_h5ad`. Both replot scripts re-add it.
- **`concat_obj` zero-fills, it does not mark missing.** A pair a split never tested (dropped by
  `min_cell`) gets `z = 0`, `p = 1`, `selected = False` — i.e. it is modelled as a *null* pair.
  On GBM that is 12.2% of the 1,662 × 13 grid, and it lets a pure **testability** difference
  produce a tiny likelihood-ratio p. `group_differential_pairs` guards the `_specific` sets with
  `tf_df.sum(1) ∈ [1, n_sub-1]`, but the raw `diff_fdr` ranking and the volcano cloud do not.
  The runner therefore emits `tested_mask.csv`, `tested_in_all_splits` /
  `confounded_by_testability` columns, and a `differential_results_interpretable.csv`.
- `_range = np.arange(1, n_sub)` also **excludes pairs selected in every sample** from the
  `_specific` sets and from the volcano's coloured points.
- `differential_test` accepts a **continuous** covariate (`x` goes straight into the design
  matrix); only `group_differential_pairs` requires 0/1. Used for the density sensitivity.

## Other deviations already identified

- `single_cell=True` (tutorial: `False`) — Xenium is single-cell; this zeroes the diagonal so a
  cell cannot signal to itself. The tutorial's `False` is for 200 µm multi-cell ST spots.
- `sig_spots(fdr=False)` — the installed **default is `fdr=True`**, but the tutorial passes
  `False`. We follow the tutorial's workflow, not the function default; switching to FDR here
  would silently change which spots are called significant.
