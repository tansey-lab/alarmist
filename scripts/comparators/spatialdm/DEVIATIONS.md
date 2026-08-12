# SpatialDM — deviations from the tutorials

Tutorials are **multi-cell spot** data (melanoma ST 200 µm pitch; intestine 100 µm pitch), ours
is **single-cell Xenium** on a **13-core TMA**. Every deviation below follows from that, from
CellChatDB v2 being newer than the bundled v1, or from a version gap.

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `single_cell` | `False` | **`True`** | Xenium is single-cell; this zeroes the diagonal so a cell cannot signal to itself. `False` is for multi-cell spots. |
| `.X` / `.raw` | log in `.X`, raw in `.raw` | same, wired explicitly | Both are genuinely used: `extract_lr` + global Moran read `.X`, but **local** spot selection reads `.raw` (`utils.spot_selection_matrix` → `adata.raw.to_adata()`, max-normalise, `sc.pp.scale`). Leaving `.raw` as the log matrix does not error — it silently changes which spots are called significant. |
| `l` | 1.2 (melanoma, = 240 µm) / 75 (intestine) | **`l = 75`**, `cutoff = 0.2` → **134.6 µm** | The two tutorials disagree **3×** on absolute radius (431 vs 135 µm) because `l` is set relative to spot pitch, so neither absolute value is an author claim about signalling range. We take the intestine value: it is the only place the authors state their reasoning (*"spot-spot distance of 100 micrometers, `l` will be set to 75… parameters should be determined to match the context of CCC"*) and 135 µm is a plausible secreted-signalling range. Deliberately **not** harmonised with CytoSignal's 200 µm or stLearn's 250 µm. |
| `eff_dist` | not used | **never used** | `stats.py:71` computes `l = sqrt(-eff_dist/(2·ln cutoff))` while the kernel is `exp(-d²/2l²)`, so `eff_dist` is a **squared** distance despite its name — `eff_dist=135, cutoff=0.2` yields `l=6.48` and a weight of 4e-95 at 135 µm instead of 0.2. We pass `l` directly. |
| `n_neighbors` | left `None` → 186 | **per-core, = (max cells within 134.6 µm in that core) + 1** (94–709) | `n_neighbors` and `n_nearest_neighbors` are **independent** — `rbfweight` derives `n_neighbors = n_neighbor_layers*31` *only when `n_neighbors is None`*. In both tutorials the cap never binds (≈14.6 and ≈5.7 units inside the radius vs 186), so **`cutoff` is the truncator and the cap is a loose ceiling**; letting it bind would invert how the method operates. GBM core density spans **19×** (423–7,993 cells/mm²), so no single value works: at 400 the cap binds for **34.7%** of cells (83% in core 1). Set per core so binding is **0.00% everywhere**. Note that once the cap does not bind, `W` is *numerically identical* for any larger value — only stored `nnz` differs. |
| `eliminate_zeros()` | n/a | **called after `weight_matrix`** | `rbfweight` zeroes sub-cutoff entries via `rbf_d[rows, cols] = 0` and never compacts, so stored `nnz = n_cells × n_neighbors` regardless of `cutoff`. Per-core stored would be 47.3M (core 1 alone 18.8M); compacting brings core 1 to 13.3M. Safe: nothing downstream touches `.nonzero()`/`.nnz` on `W` — all uses are `.multiply`/`.sum`/`@`. |
| `n_nearest_neighbors` | 6 | **6 (unchanged)** | Separate small graph for adjacent signalling. At Xenium resolution the 6 nearest cells approximate physically adjacent cells, which is what that graph is for. |
| LR database | bundled **CellChatDB v1** | **CellChatDB v2**, injected | `extract_lr` has no custom-DB hook (`datahost` only picks builtin vs figshare, both v1), so we build `uns['ligand'/'receptor'/'geneInter'/'num_pairs']` with identical logic. **Lossless** — SpatialDM handles multi-subunit complexes natively, so no v2 row is dropped for being a complex (contrast stLearn, which must drop 57.5%). |
| `Non-protein Signaling` | does not exist in v1 | **remapped to `Secreted Signaling`** | **Hard v1/v2 incompatibility.** CellChatDB v2 adds a fourth `signaling_type` (**994/3,233 rows, 31%**) that SpatialDM's `globle_st_compute` does not enumerate: it builds `st` as `hstack(repeat(nm0, #{ECM-Receptor, Cell-Cell Contact}), repeat(nm, #{Secreted Signaling}))`, so those rows get no variance term, `st` comes out short, and `spatialdm_global` dies with an `IndexError` against the length-N `idx_use` mask. It is also internally inconsistent beforehand: `n_short_lri = (annotation != 'Secreted Signaling')` counts them as short-range while `st` allocates them nothing. We remap to `Secreted Signaling` (long-range RBF), matching this benchmark's CytoSignal mapping (Secreted/ECM/Non-protein → diffusion) and the fact that non-protein ligands are diffusible neurotransmitters/metabolites — which matters in glioma. Dropping them instead would cost 31% of the DB. |
| whole-tissue vs per-core | single section | **one run per TMA core (13)** | See the note below — the original justification (cross-core weights) was **measured and is wrong** for these parameters; the real reasons are that a pooled fit returns one tissue-wide Moran's R per pair instead of 13, that `concat_obj` requires separately fitted objects, and that `extract_lr`'s `min_cell` filter would stop adapting per core (valid pairs 1,133–1,661). Also avoids a single 71M-nnz graph. |
| `sig_spots` | `fdr=False` | `fdr=False` | The installed **default is `fdr=True`**; the tutorial passes `False`. We follow the tutorial's workflow, not the function default. |
| per-LR plots | top pairs only | top pairs **plus** GRN→SORT1 and ANXA1→FPR1 | Standing request; in a separate `plots/requested/` directory. |

## `default` tier (CellChatDB v1, 2026-08-07)

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `extract_lr` DB source | `datahost='builtin'` (the default) → **downloads CellChatDB v1 from figshare, once per call** | **`datahost='package'`** → the identical CellChatDB v1 shipped inside the wheel (`spatialdm/datasets/LR_data/`) | The argument is named backwards from what it does. `'package'` gives the same 1,939-interaction v1 DB but offline, pinned to the installed version, reproducible, and not re-fetched 13 times (once per core). No content difference is intended or expected. |
| `setuptools` | n/a | **pinned `<81` (80.10.2)** in `comp-spatialdm` | `datahost='package'` imports `pkg_resources`; setuptools ≥ 81 no longer ships it, so the env's setuptools 83 made the authors' own bundled-DB branch raise `ModuleNotFoundError`. Frozen in `env.lock.yml`. |
| `PATH` | n/a | runners prepend `dirname(sys.executable)` | bokeh exports chord PNGs by driving headless `geckodriver`/`firefox`. Both are installed in `comp-spatialdm/bin/`, but this repo invokes the interpreter by absolute path (conda activate is broken), which leaves that `bin/` off `PATH` — every chord plot failed with "Neither firefox and geckodriver … available on system PATH". |
| `chord_LR` inputs | `obs` columns = deconvolution proportions per cell type | **one-hot `obsm['celltypes']` copied into `obs`** | `chord_LR` does `adata.obs.loc[:, sender]`; passing cell-type *names* raises `KeyError`. Single-cell labels are a one-hot version of the tutorial's proportions. Same accommodation already used for `chord_celltype`. |
| `Non-protein Signaling` remap | n/a | **not needed** | v1 has only the three categories SpatialDM enumerates, so this tier is free of the remap entirely — which is precisely what makes it the control for it. |

## Differential run (`run_diff_spatialdm.py`, 2026-08-06)

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `diff_utils.concat_db` | downloads **CellChatDB v1** from figshare, then `geneInter.loc[ligand.index]` | **monkeypatched** to build the union frame from the samples' own injected **v2** `geneInter` | Hard incompatibility, same class as `extract_lr`'s missing custom-DB hook: our v2 interaction names are not in the v1 index, so the tutorial call dies with a `KeyError`. The patch keeps the authors' structure verbatim (`pd.concat` → `~index.duplicated()` → `.loc[ligand.index]`) and only changes where `geneInter` comes from. It also removes the network dependency. `concat_obj`, `differential_test` and `group_differential_pairs` run **untouched**. |
| `geneInter.interaction_name` | present (v1 CSV column) | **re-added after load** | `run_spatialdm.py` set it as a column, but anndata wrote it out *as the index* (`uns/geneInter` attrs `_index: interaction_name`) and dropped the column. `compute_pathway` does `groupby('pathway_name').interaction_name`, so without the shim `dot_path` raises `AttributeError`. Same shim already present in `plot_spatialdm_full.py`. |
| `pl.dot_path(concat, '<c>_specific', cut_off=2, ...)` | called with an AnnData + a uns key | **`compute_pathway(dic={...})` first, then `dot_path(pathway_res, cut_off=2, ...)`** | **Version gap**, not a choice: installed 0.3.1 is `dot_path(pathway_res, cut_off=1, groups=None, ...)`, so the tutorial's second positional argument would bind the string to `cut_off`. The 0.3.1 route yields one panel per `pattern`, which is what the tutorial's two calls produce between them. |
| object passed to `concat_obj` | full per-sample objects | **light copies**: `X` + `obs[grade, cell_type]` + `obsm['spatial']` + the `uns` frames | `anndata.concat` drops `obsp` and `.raw` anyway, and nothing in the differential path reads `X`/`obsm`. We additionally drop the source h5ad's `pca`/`umap`/`neighbors`/`rank_genes_groups` and the per-core `obsm['celltypes']` dummy frame, whose **columns differ between cores** and would break the concat. Peak RSS 3.6 GB against 6.4 GB of persisted objects. |
| condition covariate | fixed `conditions` vector | grade **plus** a density-split sensitivity | See "Confound" above — reported alongside, never as the primary result. Uses the authors' own `differential_test`/`group_differential_pairs` with a different `conditions` vector; nothing about the method changes. |
| figure formats | inline / `pdf=` | **png + pdf + svg through one saver**, Arial, `pdf.fonttype=42`, `svg.fonttype='none'` | CLAUDE.md plotting rule. Blank figures are refused rather than written. |

**`dot_path` was not produced for the grade contrast** — both `_specific` sets came out empty, so
there was nothing to enrich. Recorded in `run_manifest.json` under `skipped_no_input`, not
`failed`.

### Correction (2026-08-07): why per-core, measured

The justification originally given here — "Moran's I over disconnected cores would create
cross-core weights and break the null" — **does not hold at our parameters, and I checked rather
than assumed**:

- Minimum **cell-to-cell** distance between any two cores: **222.9 µm** (cores 11 ↔ 14),
  independently reproducing the figure measured for the LIANA+ section.
- Our effective radius is **134.6 µm**, and the RBF weight at 222.9 µm is **0.0121**, far below
  the `cutoff = 0.2` that `rbfweight` zeroes at. A pooled run's `W` would be **exactly
  block-diagonal by construction**. It stays safe up to `l = 124` (radius 222.5 µm).
- The second confound one would expect from pooling — global standardisation mixing cores — is
  also small here: the between-core share of expression variance is a median of **0.9%** across
  all 5,119 genes (90th pct 4.3%; only 15 genes exceed 20%).

The reasons that **do** hold:

1. **`spatialdm_global` produces one Moran's R per LR pair per object.** Pooling gives a single
   tissue-wide statistic, not 13 — the per-core result table simply would not exist.
2. **`diff_utils.concat_obj(samples, names, …)` takes a list of separately fitted objects.** The
   grade differential, and the density-confound diagnosis that came out of it, are not
   constructible from a pooled fit at all.
3. **`extract_lr`'s `min_cell` filter would become global**, so a pair expressed in one core only
   would be tested tissue-wide. Per-core it adapts (1,133–1,661 valid pairs), which also changes
   the BH denominator.
4. The authors' own multi-sample tutorial fits each sample separately in a loop and only then
   concatenates — per-sample fitting **is** their design, not our departure from it.

Reason 3 also means the two schemes are not just a re-batching of the same computation: they test
different hypothesis sets. Points 1–2 are the decisive ones.

## Confound to report with any cross-grade result

**Core neighbourhood density correlates with grade: r = 0.659, p = 0.014** (point-biserial,
n=13 cores; Mann-Whitney p = 0.005). Median cells within 134.6 µm: **243 in high-grade cores
(range 78–532) vs 60 in low-grade (35–118)** — a **4× difference**.

SpatialDM's analytical z-score null derives from the structure of `W`, so statistical power
varies with neighbourhood size. High cellularity is itself a WHO glioma grading criterion, so
this is expected biology, not an artifact of our pipeline — but it means **any cross-grade
comparison of significant-LR counts is partly a power artifact**, for this method and for any
other whose power scales with neighbourhood size (including the existing CytoSignal grade
result, whose 200 µm ball has the same property).

## Unresolved

**Core 7 is absent** (`tma_id` runs 1–6, 8–14; 13 cores from 7 patients). `tma_id` is `int64`,
not categorical, so a dropped level leaves no trace, and `obs` carries no QC/filter column —
whether this was a deliberate QC exclusion or a silent upstream filter **cannot be determined
from the h5ad** and needs the provenance of whoever built it.

### Patient structure — corrected, and the design is largely PAIRED

**The h5ad has no patient column.** The `tma_id → patient` map exists only as a hardcoded dict in
`scripts/research/gbm_supplemental.py:29-38`; that script is the sole provenance for everything
below, and it is not verifiable against the data file.

Taking it at face value: **13 cores from 7 patients, not 8** — an off-by-one that had propagated
through METHODS.md (3 places) and this file; corrected 2026-08-07.

| patient | high | low | cores |
|---|---|---|---|
| 9736 | 1 | 1 | 3, 4 |
| 23184 | 1 | 1 | 5, 6 |
| 19882 | 1 | 1 | 8, 9 |
| 14007 | 2 | 2 | 11, 12, 13, 14 |
| 10097 | 1 | 0 | 1 |
| 7341 | 1 | 0 | 10 |
| 67927 | 0 | 1 | 2 |

**4 of the 7 patients contribute both a high-grade and a low-grade core, covering 10 of the 13
cores.** So the grade contrast is mostly a *within-patient* comparison, which is a considerably
stronger design than the between-subject one the analysis currently assumes — and nobody had
noticed.

`differential_test` cannot use it: it fits `y ~ 1 + conditions` by OLS with **no patient term**
(`diff_utils.py`, `x1 = np.vstack((np.ones(n_sub), x)).T`), so the within-patient pairing is
thrown away and patient 14007's 4 cores are treated as 4 independent observations. That cuts both
ways — it is pseudoreplication *and* lost power. A paired/mixed model on `zscore_df` (which
`run_diff_spatialdm.py` already persists) would be the natural follow-up, but it is **not**
SpatialDM's own workflow, so it was not run.
