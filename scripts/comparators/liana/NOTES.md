# LIANA+ — tutorial call contract

Source: `/Users/jiayifan/tansey_lab/liana-py/docs/notebooks/bivariate.ipynb` (the vignette the
benchmark's methods paragraph describes: *"identified spatial neighbors using the
spatial_neighbors function and computed bivariate scores using the bivariate function"*).
Package source read directly (`src/liana/utils/spatial_neighbors.py`, `query_bandwidth.py`,
`method/sp/_bivariate/_spatial_bivariate.py`). Installed: **liana 1.8.1**, env `comp-liana`.

## One row per tutorial call

| # | Call (tutorial values) | Produces |
|---|---|---|
| 1 | `adata.layers['counts'] = adata.X.copy()`; `sc.pp.normalize_total(target_sum=1e4)`; `sc.pp.log1p` | log-normalised input |
| 2 | `li.ut.query_bandwidth(coordinates, start=0, end=500, interval_n=20)` | **bandwidth calibration** — median neighbours (self excluded) vs a hard query **radius**; see the warning below, the x-axis is *not* the Gaussian σ |
| 3 | `li.ut.spatial_neighbors(adata, bandwidth=200, cutoff=0.1, kernel='gaussian', set_diag=True)` | `obsp['spatial_connectivities']` |
| 4 | `li.mt.bivariate(adata, resource_name='consensus', local_name='cosine', global_name='morans', n_perms=100, mask_negatives=False, add_categories=True, nz_prop=0.2, use_raw=False)` | `lrdata`: cells × LR pairs |

## Core algorithm

LIANA+'s spatial mode computes **local bivariate similarity** between a ligand and a receptor,
weighted by spatial proximity — a spatially-weighted cosine (or one of six metrics; see
`li.mt.bivariate.show_functions()`). Two levels come out at once:

- **Local**: a score *per cell per LR pair* (`lrdata.X`), with permutation p-values
  (`lrdata.layers['pvals']`) from shuffling spot/cell labels, and optional categorical labels
  (`high-high`, `low-low`, …) in `layers['cats']`.
- **Global**: per-pair summaries in `lrdata.var` — `mean`, `std`, and **bivariate Moran's R**
  (Lee's statistic), an extension of univariate Moran's I to two variables.

So LIANA+ is the only method here that natively returns **both** a per-cell score and a
per-pair global statistic with p-values, from one call.

## Spatial model

Gaussian kernel `exp(-d² / (2·bandwidth²))`, weights below `cutoff` zeroed, built on a
`max_neighbours`-capped KNN graph. Support radius = `bandwidth × sqrt(-2·ln cutoff)` = **bandwidth × 2.146** at
`cutoff = 0.1` → **28.2 µm** at our bandwidth of 13.1454.

### Calibration — SUPERSEDED, then replaced

**SUPERSEDED (do not use).** The first pass took the *bivariate* tutorial's criterion —
*"a bandwidth of 150-200 roughly includes 6 neighbours i.e. the first ring of neighbours in the
hexagonal grid of 10x Visium"* — and re-derived it with `li.ut.query_bandwidth`, giving
**bandwidth 18.75 µm / support 40.2 µm**. That value is **void**: "first ring" is a topological
property of a hexagonal lattice, and irregularly packed Xenium cells have no lattice, so the
criterion has no referent on this data. The *inflow* tutorial — the branch that actually applies
here — states no numeric rule at all, only that bandwidth "should reflect the typical range of
molecular signaling in the tissue", traded against resolution.

**In force: equal-area correspondence to an s x s patch (s = 50 µm).**

| step | value |
|---|---|
| `k = sqrt(-2·ln cutoff)`, cutoff 0.1 | 2.14597 |
| support radius `s/sqrt(pi)` | **28.2095 µm** |
| `bandwidth = support / k = s/3.804` | **13.1454 µm** |

Realised on GBM (all 100,197 cells):

| quantity | value |
|---|---|
| median neighbours / cell | **14** |
| max neighbours / cell | **52** |
| `max_neighbours` cap | 100 (LIANA default) |
| cap binding? | **no — 0.0000%** (52 < 100) |

The cap was left at the default precisely because it is already non-binding; raising it would be
a gratuitous deviation. See `../METHODS.md` for the provenance note on where the 50 µm comes from.

Both the median and the max were **re-verified independently 2026-08-04** — a plain sklearn
101-NN query on the `x`/`y` columns of `cellchatdb2/data/cell_meta.csv`, thresholded at the 28.2096 µm
support radius, reproduces median **14** and max **52** (including self, since `set_diag=True`),
with the 101-neighbour cap binding for **0** cells.

> ⚠️ **`connectivity_nnz` in the manifest is not a neighbour count.** `cellchatdb2/run_manifest.json`
> records `connectivity_nnz: 10,119,897`, which is exactly 100,197 × **101** — i.e. `max_neighbours`
> + self, for every cell. LIANA builds the capped KNN and then zeroes weights below `cutoff`
> *without eliminating them from the sparse structure*, so `nnz` counts **stored entries including
> explicit zeros**. The real degree comes from `(conn > 0).sum(axis=1)`, which is what
> `run_liana.py:126` uses and what gives median 14. Reading `nnz / n_cells` as "101 neighbours per
> cell" would contradict the max of 52 — and be wrong.

### ⚠️ `query_bandwidth`'s x-axis is a RADIUS, not the bandwidth

`li.ut.query_bandwidth` does **not** sweep the Gaussian σ. Reading the installed source
(`liana/utils/query_bandwidth.py`), each point is a `BallTree.query_radius(r=max_distance,
count_only=True)` — a **hard cutoff radius** — and the y value is **`ceil(median(count)) − 1`**
(`query_bandwidth.py:71-72`), i.e. the *median* neighbour count with self removed. The local variable
is named `avg_nn`, but it is a **median**, not a mean — the name is wrong in the source.

That one line reconciles three numbers that otherwise look inconsistent (re-verified 2026-08-06 by a
direct `BallTree.query_radius` on `cellchatdb2/data/cell_meta.csv` at r = 28.209581801983987):

| statistic | value | why it differs |
|---|---|---|
| **mean** neighbours, self excluded | **14.65** | a mean, not a median — pulled up by dense cells (max 52) |
| **median** neighbours, self **included** | **14** | the realised whole-run median, because `set_diag=True` |
| what the curve plots at that radius | **13** | `ceil(14) − 1` — the same median with self removed |

So the two axes are on different scales, and comparing a σ against them is a category error:

| point | value | curve value (`bandwidth_choice/data/query_bandwidth_tutorial_5_35.csv`) |
|---|---|---|
| Gaussian σ (our bandwidth) | 13.1454 µm | 3 — **meaningless as a neighbour count** |
| **support radius** `σ × 2.14597` | **28.2095 µm** | **13** (interp. 13.17; exactly 13 at r = 28.077) |
| realised median incl. self, whole run | — | **14** = curve value + 1 ✅ |

> The older row for this table read "13.46, interpolated from
> `cellchatdb2_inflow/data/bandwidth_query.csv`". That file's 2.95 µm grid actually interpolates to
> **13.61**, and 13.46 reproduces from neither file — **corrected 2026-08-06** to the 0.77 µm sweep
> above (see the two-sweeps table below). The substantive point is unchanged, but state it as
> *curve 13, realised median 14, and the +1 is the self that `query_bandwidth` subtracts* — not as
> an approximate match.

**Converting between the two axes has a direction, and it is easy to get backwards.** A value read
off this curve is a **radius**; `spatial_neighbors(bandwidth=…)` wants a **σ**. So

```
sigma = radius_read_off_the_curve / 2.145966      # NOT the other way round
```

Passing a curve reading straight through as `bandwidth=` inflates the support radius by 2.146× and
the neighbourhood **area by 4.6×**. This is an inconsistency in the tutorial itself — it plots one
quantity and then asks you to supply another — not something we introduced.

> ⚠️ **Retracted 2026-08-07 — the figure was never wrong.** Both this note and the 2026-08-06 pass
> asserted that the dashed vline in `cellchatdb2_inflow/plots/global/bandwidth_query.png` is drawn
> at σ = 13.1454 on a radius axis, and listed "regenerate it" as an open item. **It is not.**
> `run_inflow_downstream.py:163-167` computes `R = a.bandwidth * np.sqrt(-2 * np.log(a.cutoff))`
> and draws `geom_vline(xintercept=float(R))` — i.e. at the **support radius 28.2096 µm**, the
> correct scale — with the label `support radius = 28.2 um (gaussian sigma = 13.1454)`. Reading the
> PNG on disk confirms it: the guide sits at ≈28 on the x-axis and crosses the curve at ≈14
> neighbours, exactly where it should. The claim was never checked against either the code or the
> image; it is recorded here because a plausible wrong claim survived three documents and a review.

The one real blemish in that figure is cosmetic: the rotated annotation is anchored at
`y = dfb.neighbours.max()` with `va="bottom"`, so it is pushed above the panel and clipped — only
the first two characters (`su`) are visible. The guide line, the axis and the kernel are all
correct, and **nothing needs regenerating for correctness**. The comparable point remains the
**support radius 28.2095 µm**, where the curve gives 13.

Note also that our calls use a rescaled sweep rather than the tutorial's `0 / 500 / 20` recorded in
the contract table above — the tutorial's range is Visium-scaled and would put our entire kernel
inside the first bin. Two sweeps exist, and they are **not interchangeable** (**corrected
2026-08-06**; this note previously gave a single `start=5, end=60, interval_n=40`, which matches
neither file on disk):

| sweep | call | grid | file |
|---|---|---|---|
| inflow run | `start=5, end=120, interval_n=40` (`run_inflow_downstream.py:157`) | 2.95 µm steps | `cellchatdb2_inflow/data/bandwidth_query.csv` |
| bandwidth study | `start=5, end=35, interval_n=40` (`choose_bandwidth.py:127,144`) | 0.77 µm steps | `bandwidth_choice/data/query_bandwidth_tutorial_5_35.csv` |

Read fine-grained values off the 5–35 sweep; the 5–120 one is for the shape of the curve.

### Whole-slide is safe here — verified

| | |
|---|---|
| LIANA support radius | **28.2 µm** |
| minimum inter-core cell–cell distance | **222.9 µm** |
| cross-core pairs at 28.2 µm | **0** (still 0 at 200 µm) |
| `max_neighbours=100` cap | binds for **0.0000%** (median 14 neighbours, max 52) |

A **7.9×** margin, so the 13 TMA cores are analysed **together, in one run** — no splitting.

## Outputs

| Key | Shape | Meaning |
|---|---|---|
| `lrdata.X` | cells × LR pairs | local bivariate score |
| `lrdata.layers['pvals']` | cells × LR pairs | permutation p-value per cell |
| `lrdata.layers['cats']` | cells × LR pairs | local category (high-high / low-low / …); written only when `add_categories=True` (`_spatial_bivariate.py:260-261`) |
| `lrdata.var` | pairs × **10** | `ligand`, `receptor`, `ligand_means`, `ligand_props`, `receptor_means`, `receptor_props`, **`morans`**, `morans_pvals`, `mean`, `std` |
| `adata.obsp['spatial_connectivities']` | cells × cells | the kernel graph |

The `lrdata.var` row is **corrected 2026-08-04** — it previously listed only
`mean, std, morans, morans_pvals`. The persisted table
`cellchatdb2/data/global_scores.csv` is (131, 10) with the ten columns above; the four
`*_means`/`*_props` columns are the per-gene expression summaries `nz_prop` filters on
(`_spatial_bivariate.py:183-195`), and they are what makes the file self-contained enough to
re-derive that filter without the source h5ad.

Pair naming is `<ligand>^<receptor>`, with `_` separating heteromeric subunits
(e.g. `VTN^ITGAV_ITGB5`) — so `^` and `_` mean different things. **The subunit ORDER inside a
complex is database-export-specific** (`TGFBR1_TGFBR2` vs `TGFBR2_TGFBR1`) and differs between
these runs and `results/GBM/` — see the LR-database provenance section of `DEVIATIONS.md` before
joining LR strings across the two.

## What the MOFA-Flex figure calls actually draw

`inflow_mofaflex.ipynb`'s factor-inspection cells (31–43) and its circle-plot cell (59) call three
plotting functions whose output does not match the obvious reading of the axes. Read from the
installed **`mofaflex 0.1.0.post2.dev179+g9792b435f`** (env `comp-liana`) on 2026-08-06, not assumed.

### `mfl.pl.top_weights` — `_plotting.py:1118-1171`

- The x aesthetic is **`weightabs`**, labelled `"| Weight |"` — an **absolute** loading, so a large
  bar can be a large *negative* weight.
- Sign is carried **only by the glyph**: `p9.aes(shape="weightsgn")` with
  `scale_shape_manual(values=("$\\oplus$", "$\\ominus$"), breaks=(True, False), guide=None)`.
  `weightsgn` is `weight >= 0`, and **`guide=None` means no legend is drawn for the shape** — the
  ⊕/⊖ distinction is undocumented on the figure itself.
- Selection is `groupby("factor").apply(… argsort()[-n_features:])`, then
  `sort_values(["factor", "weightabs"], ascending=True)` — top *n* by **`|weight|`**, plotted
  ascending, so the **largest is at the TOP** of each facet.
- `facet_wrap("factor", scales="free")` — **each panel has its own x-scale**. Bar lengths are **not**
  comparable across factors.
- The weights are raw model weights. There is **no** prevalence or sparsity normalisation anywhere
  in this call.

### `mfl.pl.variance_explained` — `_plotting.py:477-526`

- It is a **heatmap** (`p9.geom_tile()`), not a bar chart:
  `p9.aes(x=x, y="factor", fill="R2")` with `scale_fill_distiller(palette="OrRd")`.
- With the default `group_by="group"`, **x = view** (for an inflow MuData, one view per **sender**
  cell type) and **y = factor**; `facet_wrap(group_by)` gives **one facet per group**, so a
  single-group model is one panel.
- `get_r2(..., ordered=True)` sorts factors by total R² descending and the order is frozen into a
  `Categorical`; plotnine's discrete y-scale then puts the **largest at the BOTTOM**.
- Column sums over factors are exactly the per-view totals we persist as `r2_per_view.csv` (verified
  exact on the reach-norm fit, max residual 1e-15).
- `group_by="view"` swaps the roles (x = group, faceted by view); with one group that degenerates to
  *n*-view one-column facets.
- **Corrected 2026-08-06:** an earlier note described this figure as bars labelled by sender. It is
  a factor × sender heatmap.

### `li.pl.circle_plot` as the tutorial wires it (cell 59) — `run_mofaflex.py:624-676`

`li.pl.circle_plot` reads everything from the `liana_res` DataFrame it is handed; **it never sees
the model.** In this call:

| what you see | where it comes from |
|---|---|
| `source` node | the feature name `<sender>^<lig>^<rec>` — **from the model** |
| *which* edges are drawn | top *n* by absolute loading — **from the model** |
| `target` node | `inflow_means`, i.e. `lrdata.to_df().groupby(obs[cell_type]).mean()` — **the receiving cell's own annotation, not the model** |
| edge weight (`score_key="inflow"`) | the same **raw** mean inflow — **not factor-weighted, and sign-blind** |

Consequence, verified on the reach-norm fit: **Factor 1 and Factor 6 share 108 identical
`(source, LR, target)` rows whose edge weights differ by exactly `0.000e+00`, while their loadings
differ by up to 1.90.** Two factors that select the same feature therefore render the *identical*
sub-network. Read the figure as *"the top-n interactions this factor selects, and where those
generally go"* — **not** as *"this factor's sender→receiver structure"*.

## Facts confirmed from source, not assumed

Re-verified against installed **liana 1.8.1** in `comp-liana` on 2026-08-04; all of the following
still hold, with two additions.

- Kernel is `exp(-d²/(2·bandwidth²))` (`_gaussian`, `spatial_neighbors.py:16`); alternatives are
  `exponential`, `linear`, `misty_rbf` (note `misty_rbf` omits the factor 2 —
  `spatial_neighbors.py:19-25`). ✅ confirmed verbatim.
- `max_neighbours` defaults to **100** (`spatial_neighbors.py:51`) — a KNN ceiling on the graph,
  exactly the class of cap that bit SpatialDM. Verified non-binding here. ✅
- `li.mt.bivariate` accepts `resource: pd.DataFrame` with `['ligand','receptor']` columns, so
  CellChatDB v2 hands over directly and **losslessly** — LIANA handles heteromeric complexes
  natively (`_`-joined subunits), unlike stLearn. ✅
- `nz_prop` defaults to **0.05** in the function signature (`_spatial_bivariate.py:58`); the
  tutorial overrides it to 0.2. ✅ (and `_inflow.py:41` defaults to **0.001** — see
  `DEVIATIONS.md` on the four values in play).
- **`set_diag` defaults to `False`** (`spatial_neighbors.py:53`). Both the tutorial and our runs
  pass `set_diag=True` explicitly, so a cell is its own neighbour — which is why the realised
  median of 14 above includes self.
- **`query_bandwidth` sweeps a hard radius, not the bandwidth.** See the warning in the
  Calibration section; the parameter name is doubly misleading because the returned column is also
  misspelled.

### Package bug

`li.ut.query_bandwidth` returns a DataFrame whose column is spelled **`bandwith`** (missing the
second `d`) — confirmed at `query_bandwidth.py:61`, and the saved
`cellchatdb2_inflow/data/bandwidth_query.csv` carries that spelling in its header. Code indexing
`'bandwidth'` raises `KeyError`.
