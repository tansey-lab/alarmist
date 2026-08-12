# NICHES — tutorial-call contract

Source: **local clone** `/Users/jiayifan/tansey_lab/NICHES`, v1.2.4,
git `d698e37b8c38ebd103c34acd0b35b03b48d3c5a3` (2026-01-29, merge of PR #69 from `dev`).
No WebFetch — every call below was read out of the vendored `vignettes/*.Rmd` and verified
against the installed `R/*.R` source.

The dataset is **spatial, single-cell, multi-sample**, so the governing vignette is
**`07 Spatiotemporal NICHES.Rmd`** (per-sample spatial NICHES → tag Condition → merge →
differential), with `01 NICHES Spatial.Rmd` for the single-section spatial idioms and
`04 Differential CellToCell.Rmd` for the differential-analysis idioms. Vignettes 02/03/05/06/08/09
are non-spatial or interoperability demos and are not used.

## Which NICHES "organization" we compute

`RunNICHES` can emit six atlases. Vignette 07 — the multi-sample **spatial** one — runs

```r
CellToCellSpatial = T, NeighborhoodToCell = T, CellToCell = F
```

and that is what we run. `CellToNeighborhood` (the outgoing mirror of `NeighborhoodToCell`)
and the three non-spatial organizations (`CellToCell`, `CellToSystem`, `SystemToCell`) are
left at their package defaults (`FALSE`); `CellToCell` in particular crosses every cell pair
regardless of distance and discards the spatial structure that motivates this dataset.

| Organization | Unit of a column | Default | Ours |
| --- | --- | --- | --- |
| `CellToCell` | one sending–receiving cell pair, non-spatial | `TRUE` | **`FALSE`** (vignette 07) |
| `CellToSystem` | one sending cell → whole system | `FALSE` | `FALSE` |
| `SystemToCell` | whole system → one receiving cell | `FALSE` | `FALSE` |
| `CellToCellSpatial` | one **spatially adjacent** cell pair | `FALSE` | **`TRUE`** (vignette 07) |
| `CellToNeighborhood` | one cell → its spatial neighbours | `FALSE` | `FALSE` |
| `NeighborhoodToCell` | a cell's **niche** (neighbours → cell) | `FALSE` | **`TRUE`** (vignettes 01 + 07) |

## Call contract — one row per tutorial call

`run_niches.R` implements exactly this table, in this order. "Tutorial value" is the literal
argument in the vignette; where we differ the cell says so and the row is repeated in
`DEVIATIONS.md`.

### Stage A — per sample (per TMA core)

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| A1 | `CreateSeuratObject(counts)` | `LoadData("stxBrain")` / `load(*.Robj)` | our `counts.mtx` (genes × cells, raw integer) | 01:40, 07:36 |
| A2 | metadata `x` / `y` | `brain@meta.data$x <- ...coordinates$row` | `x`/`y` ← `obsm['spatial']` (µm) | 01:61-62 |
| A3 | `NormalizeData(object)` | `NormalizeData(brain)` — all defaults | identical (`LogNormalize`, `scale.factor = 1e4`) | 01:65, 04:46 |
| A4 | `SeuratWrappers::RunALRA(object)` | `brain <- SeuratWrappers::RunALRA(brain)` | run for the `alra` sub-run; skipped for the `noimpute` sub-run | 01:72 |
| A5 | `object@meta.data$cell_types <- Idents(object)` | same | `cell_type` column from the h5ad | 07:58 |
| A6 | **`RunNICHES(...)`** | see below | see below | 07:60-71 |

`RunNICHES` arguments — every one either the vignette's literal value or the package default:

| Argument | Package default | Vignette 07 | Vignette 01 | Ours |
| --- | --- | --- | --- | --- |
| `assay` | `"RNA"` | `'alra'` | `"alra"` | `"alra"` (alra sub-run) / `"RNA"` (noimpute sub-run) |
| `LR.database` | `"fantom5"` | `'fantom5'` | `"fantom5"` | **`"custom"`** — the `cellchatdb2` tier |
| `custom_LR_database` | `NULL` | — | — | CellChatDB v2.0 human, 2 cols, `_`-separated |
| `species` | (required) | `'mouse'` | `"mouse"` | `"human"` (ignored when `LR.database = "custom"`, `LoadCustom.R:12`) |
| `cell_types` | `NULL` | `"cell_types"` | `"seurat_clusters"` | `"cell_type"` |
| `position.x` / `position.y` | `NULL` | `"x"` / `"y"` | `'x'` / `'y'` | `"x"` / `"y"` |
| `k` | `4` | *(not passed → 4)* | `4` | **`4`** (default; skill invariant: keep the method's own neighbourhood definition) |
| `rad.set` | `NULL` | *(not passed)* | *(not passed)* | `NULL` (ignored when `k` is set, `ComputeEdgelist.R:44`) |
| `blend` | `"mean"` | *(not passed)* | *(not passed)* | `"mean"` |
| `min.cells.per.ident` | `NULL` | *(not passed)* | `0` | `NULL` (default) |
| `min.cells.per.gene` | `NULL` | *(not passed)* | `NULL` | `NULL` (default) |
| `meta.data.to.map` | all meta cols | *(not passed)* | `c('orig.ident','seurat_clusters')` | `c('cell_type','grade','tma_id','x','y')` |
| `output_format` | `"seurat"` | *(not passed)* | *(not passed)* | `"seurat"` |
| `CellToCellSpatial` | `FALSE` | `T` | `F` | `TRUE` |
| `NeighborhoodToCell` | `FALSE` | `T` | `T` | `TRUE` |
| `CellToCell` | `TRUE` | `F` | `F` | `FALSE` |

### Stage B — merge across samples, embed

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| B1 | extract organization | `niches.list[[i]][["NeighborhoodToCell"]]` | same, per core | 07:79 |
| B2 | tag condition | `menv[[i]]$Condition <- names(niches.list)[i]` | `Condition <- grade` (`high`/`low`), plus `Core <- tma_id` | 07:80 |
| B3 | `merge()` | `merge(menv[[1]], menv[2])` | `merge(x[[1]], x[-1])` — 13 cores | 07:82 |
| B4 | low-information filter | `subset(scc.merge, nFeature_CellToCell > 5)` | same idiom on `nFeature_<organization> > 5`; **CellToCellSpatial only**, as in vignette 04 | 04:92 |
| B5 | `ScaleData` | `ScaleData(menv)` | identical | 07:86 |
| B6 | `FindVariableFeatures(selection.method = "disp")` | identical | identical | 07:88 |
| B7 | `RunPCA(npcs = 100)` | `RunPCA(menv, npcs = 100)` | identical | 07:90 |
| B8 | `ElbowPlot(ndims = 100)` | identical | identical | 07:91 |
| B9 | `PCHeatmap(dims = 40:48, cells = 100, balanced = T)` | identical | identical | 07:92 |
| B10 | `RunUMAP(dims = 1:50)` | `RunUMAP(menv, dims = 1:50)` | identical | 07:93 |
| B11 | `DimPlot` × {ident, Condition} | `DimPlot(menv)`, `DimPlot(menv, group.by='Condition')` | + `group.by = 'Core'`, `'SendingType'`, `'ReceivingType'` (04:110-113) | 07:94-96 |

### Stage C — differential signalling, high vs low grade

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| C1 | isolate a receiving population | `ROI <- "Neural tube and notochord"; subset(menv, idents = ROI)` | loop over **all 9** `ReceivingType`s, one sub-analysis each | 07:103-104 |
| C2 | `ScaleData` → `FindVariableFeatures("disp")` → `RunPCA(npcs=50)` → `ElbowPlot` → `RunUMAP(dims=1:40)` | identical | identical | 07:106-113 |
| C3 | `Idents(subs) <- subs[['Condition']]` | identical | identical | 07:127 |
| C4 | `FindAllMarkers(min.pct = 0.25, only.pos = T, test.use = "roc")` | identical | identical | 07:128 |
| C5 | drop infinite-differential markers | `marker$ratio <- pct.1/pct.2; marker[ratio < Inf,]` | identical | 07:130-131 |
| C6 | `top_n(20, myAUC)` | identical | identical | 07:133 |
| C7 | `DoHeatmap(group.by="ident", features=GOI)` | identical | identical | 07:135 |
| C8 | `FeaturePlot(subs, '<LR>')` | `FeaturePlot(subs,'Dll1—Notch4')` | method's own top LRs **and** the requested pair, in separate dirs | 07:139 |
| C9 | map back as an assay | `CreateAssayObject(data = niches.data)` + `ScaleData` | identical, per core | 07:145-157 |
| C10 | `VlnPlot(..., split.by = 'orig.ident')` | `VlnPlot(merge,'Dll1—Notch4', assay="NeighborhoodToCell", split.by='orig.ident')` | `split.by = 'grade'` | 07:161-162 |
| C11 | ligand/receptor expression check | `VlnPlot(merge, c('Dll1','Notch4'), assay='alra', split.by='orig.ident')` | same, for each requested LR's genes | 07:165 |
| C12 | spatial LR map | `SpatialFeaturePlot(brain, features=c(...), slot='scale.data')` | no Seurat `@images` for Xenium → equivalent `ggplot` scatter on x/y | 01:145-147 |

## Per-LR plots — two separate sets (skill invariant)

* `plots/top_lr/` — the LRs **NICHES itself** ranks highest (`FindAllMarkers` `myAUC`, ROC test).
* `plots/requested_lr/` — always **`GRN—SORT1`** and **`ANXA1—FPR1`** (the ALARMIST motif-1
  mGAM ⇄ MES-like loop), whatever their rank.

Both LRs are present in CellChatDB v2.0 human (`GRN—SORT1`, pathway GRN; `ANXA1—FPR1`, pathway
ANNEXIN, both "Secreted Signaling") and all four genes are on the Xenium 5K panel, so both
survive `FilterGroundTruth` and are computable. NICHES joins ligand and receptor with an
**em-dash `—` (U+2014)**, not a hyphen (`RunNeighborhoodToCell.R:61`) — every feature lookup
must use the em-dash.

## What NICHES actually computes

No statistical model, no null, no p-value at the interaction-scoring step. For an edge
(sending cell *i* → receiving cell *j*) and mechanism *L→R*, the score is the **product of the
normalized ligand expression in *i* and the normalized receptor expression in *j***, with
multi-subunit complexes multiplied across subunits (`Reduce('*', subunit.list)`,
`RunNeighborhoodToCell.R:38,55`). `NeighborhoodToCell` then averages that product over all
edges landing on *j* (`blend = "mean"`). The output is a **mechanism × cell matrix** that is
handed back as an ordinary Seurat assay; **all** inference (variable features, PCA, UMAP,
`FindAllMarkers`) is downstream Seurat run on that matrix. So "significance" in NICHES is
whatever test you run on the transformed object — here the vignette's ROC test — not a
property of the LR scoring itself.

## Spatial model

Mutual k-nearest-neighbour graph on the raw x/y coordinates, `k = 4` (default).
`ComputeEdgelist.R:46-56`: rank every cell's neighbours by euclidean distance, keep the
`k+1` nearest, then symmetrize with `adj & t(adj)` — so an edge survives only if **both**
cells list the other in their top-4. There is **no distance cutoff and no kernel**: `k` is
unitless, so the graph is invariant to whether coordinates are in µm, pixels or array
indices. `rad.set` offers a hard radius instead but is ignored whenever `k` is non-NULL.

## Gotchas found by reading the source

1. **`ComputeEdgelist` is dense O(N²).** `apply(df, 1, ...)` at `ComputeEdgelist.R:36`
   materialises the full N × N double distance matrix, then `adj_mat`, `t(adj_mat)` and
   `1*(...)` — roughly **4 × N² × 8 bytes**. At the full slide (100,197 cells) that is ~80 GB
   for the distance matrix alone. This is the hard reason the slide is split per TMA core.
2. **The `nn.method = 'aoz'` fast path is commented out** (`ComputeEdgelist.R:76-131`) and
   `RunNICHES` never passes `nn.method` anyway. The dense path is the only path.
3. **`as.matrix()` on the mechanism × edge matrices** (`RunNeighborhoodToCell.R:35,52,65`)
   densifies; memory scales with `n_mechanisms × n_edges`, not with sparsity.
4. **Em-dash, not hyphen**, in every feature name.
5. **`species` is silently ignored** when `LR.database = "custom"` (`LoadCustom.R:12`).
6. **`FilterGroundTruth` requires every subunit present** (or NA) in the object
   (`FilterGroundTruth.R:17-19`) — with the 5,119-gene Xenium panel this takes CellChatDB v2
   from 3,218 unique pairs down to **1,088**.
7. **Seurat v5 does not populate the `data` slot** on `CreateSeuratObject(counts=)`; NICHES
   patches it by hand (`RunNeighborhoodToCell.R:77-81`). The NICHES assay therefore holds the
   raw LR products in **both** `counts` and `data` — it is not re-normalized. Do not
   `NormalizeData` it again.
