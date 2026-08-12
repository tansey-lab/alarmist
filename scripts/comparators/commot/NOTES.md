# COMMOT — tutorial call contract

Source: `/Users/jiayifan/tansey_lab/COMMOT/docs/notebooks/{Basic_usage,visium-mouse_brain}.ipynb`.
Package source read directly (`commot/tools/_spatial_communication.py`,
`commot/preprocessing/_ligand_receptor_database.py`). Installed: **commot 0.0.3**, env `comp-commot`.

## ⚠️ The local clone is NOT the installed package (verified 2026-08-10)

`/Users/jiayifan/tansey_lab/COMMOT/` is **ahead of** the installed PyPI release 0.0.3. Any line
number cited from the clone must be re-checked against
`/Users/jiayifan/anaconda3/envs/comp-commot/lib/python3.10/site-packages/commot/` before it is
relied on — the skill requires verifying against the *installed* package.

| File | clone vs installed | Consequence |
|---|---|---|
| `tools/_spatial_communication.py` | **identical** | every citation below (incl. `:390`, the `obsp['spatial_distance']` escape hatch) is valid for both |
| `_optimal_transport/_cot.py` | **differs** | see the two commits below |

Clone `git log`: `d117445` *"Update cost construction for new behaviors of np.where on sparse
matrices"* adds `_cost_matrix_within_cutoff()`, which accepts a **sparse** cost matrix;
`ff8f1b7` *"Fix np.Inf to np.inf in _usot.py"* — i.e. **upstream has already fixed the numpy≥2
blocker that `DEVIATIONS.md` records as a hard version gap.** Neither is in 0.0.3, so neither
was in effect for any run in `results/`.

Consequence worth knowing: `spatial_communication` in *both* versions honours a precomputed
`adata.obsp['spatial_distance']` (installed `:389-393`), so with the clone's sparse support a
whole-slide run might be feasible via a radius-neighbour graph. **Untested**, and it has two
known traps — the sparse path treats implicit zeros as absent edges (so self-distance = 0 is
dropped, silently removing autocrine signalling, which the dense path includes), and the memory
bottleneck merely moves from the distance matrix to the OT cost matrix. Do not treat this as a
supported route without a smoke test.

## Signatures verified against the INSTALLED package (2026-08-10)

- `cluster_communication(adata, database_name, pathway_name=None, lr_pair=None, clustering,
  n_permutations=100, random_seed=1, copy=False)` — **`random_seed` exists and defaults to 1**;
  runs before 2026-08-10 never passed it, so their permutation p-values were not pinned to the
  run seed. Also accepts `lr_pair=`, which is the only way to get a significance test for a
  *specific* requested LR pair.
- `communication_direction(adata, database_name, pathway_name=None, lr_pair=None, k=5,
  pos_idx=None)` — writes `obsm['commot_{sender,receiver}_vf-<db>-<tag>']` for exactly one key
  (lr_pair, else pathway, else `total-total`).
- `plot_cell_communication(..., plot_method='cell', background='summary', summary='sender',
  filename=...)` — **`background` defaults to `'summary'` and needs no H&E image.**
- `plot_cluster_communication_network(adata, uns_names=[...], clustering, filename=...)` — with
  the default `nx_node_pos='cluster'` it reads `uns['cluster_pos-<clustering>']`, so
  `ct.tl.cluster_position(adata, clustering=...)` must run first.
- `plot_cluster_communication_dotplot(adata, database_name, pathway_name, clustering,
  filename=...)`.

## One row per tutorial call

| # | Call (tutorial values) | Produces |
|---|---|---|
| 1 | `adata.raw = adata`; `sc.pp.normalize_total(adata)`; `sc.pp.log1p(adata)` | non-negative abundance-like values |
| 2 | `ct.pp.ligand_receptor_database(species, signaling_type='Secreted Signaling', database='CellChat')` | 3-column `df_ligrec` = ligand, receptor, pathway |
| 3 | `ct.pp.filter_lr_database(df, adata, min_cell_pct=0.05)` | pairs where both sides are expressed in ≥5% of spots |
| 4 | `ct.tl.spatial_communication(adata, database_name, df_ligrec, dis_thr=500, heteromeric=True, pathway_sum=True)` | the OT communication inference |
| 5 | `ct.tl.communication_direction(adata, database_name, pathway_name/lr_pair, k=5)` | interpolated signalling vector field |
| 6 | `ct.tl.cluster_communication(adata, database_name, pathway_name, clustering, n_permutations=100)` | cell-type × cell-type matrix + permutation p-values |
| 7 | *(downstream)* `ct.tl.communication_deg_detection` → `communication_deg_clustering` | signalling-dependent DE genes |

## Core algorithm

COMMOT solves a **collective optimal transport** problem: ligand "mass" at sender cells is
transported to receptor "mass" at receiver cells, minimising a spatial cost, subject to
`dis_thr` forbidding any coupling beyond that distance. It is *competitive* — multiple ligands
and receptors compete for the same mass, which is the property the authors emphasise over
pairwise scoring. Output is a **cell × cell transport plan per LR pair**, summarised to per-cell
sent/received amounts. There is **no per-pair significance test**; significance appears only at
the cell-type level via `cluster_communication`'s permutations.

## Outputs

| Key | Shape | Meaning |
|---|---|---|
| `obsp['commot-<db>-<lig>-<rec>']` | cells × cells | **transport plan per LR pair** (sparse) |
| `obsp['commot-<db>-<pathway>']` | cells × cells | per-pathway sum |
| `obsp['commot-<db>-total-total']` | cells × cells | total |
| `obsm['commot-<db>-sum-sender']` | cells × (pairs+pathways) | amount sent, `s-<tag>` columns |
| `obsm['commot-<db>-sum-receiver']` | cells × (pairs+pathways) | amount received, `r-<tag>` columns |
| `uns['commot-<db>-info']` | — | `df_ligrec` + `distance_threshold` |
| `uns['commot_cluster-<clust>-<db>-<pw>']` | types × types | `communication_matrix` + `communication_pvalue` |

## The downstream chain needs R — see `SETUP_tradeSeq.md`

`ct.tl.communication_deg_detection` is the only function in this pure-Python package that requires
R. Its `import rpy2` / `import anndata2ri` sit **inside the function body** (`:100-104`), so
`import commot` and every other call work without R and the dependency stays invisible until you
call it. It drives `tradeSeq::fitGAM` + `associationTest` and `clusterExperiment::clusterExpressionPatterns`.

- **Setup, five known breakage points, data-contract traps, runtime plan** → `SETUP_tradeSeq.md`
- **Runner for the full authors' chain** (deg detection → clustering → dependent-genes heatmap →
  `communication_impact` → impact heatmap → the tutorial's 3-panel example figure) →
  `run_commot_deg.py`. **Never executed** — written against the installed source; treat the first
  run as a test.
- Runner for `communication_impact` alone, with a substitute gene list → `run_commot_impact.py`.

Two contract facts worth repeating here because they cost time:
`communication_deg_detection` reads **`adata.layers['counts']`** (the docstring says `'count'`,
singular — the code at `:139` wins), and `communication_impact` reads **`adata.raw`**. Our saved
`adata_commot.h5ad` has neither, so both are rebuilt from the source h5ad's counts. Neither needs
an OT re-run.

## Plots the standard workflow produces

`ct.pl.plot_cell_communication` (vector field over sender/receiver, needs
`communication_direction` first), `plot_cluster_communication_network`,
`plot_cluster_communication_dotplot`, `plot_communication_impact`,
`plot_communication_dependent_genes`.

## Facts confirmed from source, not assumed

- **`dis_thr` is a plain Euclidean cutoff in whatever units `obsm['spatial']` carries**
  (`cost_type='euc'` → `self.M = dmat`, `self.cutoff = dis_thr`). It is *not* in microns unless
  the coordinates are.
- **Heteromeric complexes are supported natively** — `heteromeric=True`,
  `heteromeric_delimiter='_'`, `heteromeric_rule='min'` (subunits aggregated by minimum).
- `signaling_type=None` returns **all** categories; the tutorial deliberately restricts to
  `'Secreted Signaling'`.
- `df_ligrec` is literally 3 columns (ligand, receptor, pathway) with `_`-joined subunits — the
  **same encoding our CellChatDB v2 CSV already uses**, so the v2 handover is direct and lossless.
- **A dense N×N distance matrix is materialised** (`scipy.spatial.distance_matrix`,
  `_spatial_communication.py:390`). This, not the OT solver, sets the feasibility ceiling.

### Calibration: what `dis_thr = 500` means physically

The tutorial runs on `V1_Mouse_Brain_Sagittal_Posterior`, whose `obsm['spatial']` is in
**full-resolution image pixels**, not microns — so the notebook's prose ("500 µm") does not match
its own units. Measured on the actual dataset (3,355 spots):

| quantity | measured |
|---|---|
| nearest-neighbour spacing | **137.0 units** |
| Visium hex pitch | 100 µm (hard geometric constant) |
| → scale | **0.7299 µm per coordinate unit** |
| → `dis_thr = 500` units | **365 µm** |

A cross-check via `spot_diameter_fullres` (89.51 units for a 55 µm spot) gives 0.6144 µm/unit,
a 19% disagreement; the hex pitch is the more reliable calibration and is what we use.

**→ `dis_thr = 365` for micron coordinates.** Deliberately not harmonised with CytoSignal's
200 µm, stLearn's 250 µm, or SpatialDM's 135 µm.
