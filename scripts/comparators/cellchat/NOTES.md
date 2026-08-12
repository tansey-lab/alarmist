# CellChat — tutorial-call contract

Source: **local clone** `/Users/jiayifan/tansey_lab/CellChat`, `DESCRIPTION` Version
**2.2.0.9001**, git `75253cd0c9e68410e6e721a6d3a0419a1d7e358f` (2026-03-04, "Update analysis.R").
No WebFetch — every call below was read out of the vendored `tutorial/*.Rmd` and verified
against the installed `R/*.R` source and `man/*.Rd`.

## Which vignettes govern

Our dataset (GBM / IDH-mut LGG TMA) is **spatial, single-cell resolution, multi-section,
two conditions**. No single vignette covers that, so the contract is assembled from four,
in this precedence order:

| # | Vignette | Governs |
| --- | --- | --- |
| V2 | `CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.Rmd` | **primary** — our data shape (many sections in one object via `meta$samples`, one `spatial.factors` row per section) |
| V1 | `CellChat_analysis_of_spatial_transcriptomics_data.Rmd` | the single-section spatial idioms V2 defers to ("check the vignette of applying CellChat to an individual spatially resolved dataset for detailed descriptions") |
| VF | `FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.Rmd` | **Xenium** `spatial.factors`, and the small-gene-panel escape hatch |
| VC | `Comparison_analysis_of_multiple_datasets.Rmd` | high-grade vs low-grade — the native cross-condition mode |
| VB | `CellChat-vignette.Rmd` | the downstream plot inventory both spatial vignettes defer to ("CellChat's various functionality can be used… check other functionalities in the basic tutorial") |

Not used: `CellChat_analysis_of_spatial_multiomics_data.Rmd` (multiome, not our data),
`Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.Rmd` (our two
conditions share all 9 cell types — see below), `Interface_with_other_single-cell_analysis_toolkits.Rmd`.

## Object layout for GBM — two objects, thirteen samples

VF:118 and V1:51 both state the rule explicitly: *"for comparison analysis across different
conditions, users still need to create a CellChat object separately for each condition"*, while
multiple sections **of the same condition** go into **one** object through a `samples` column.
So:

| Object | `meta$samples` levels | cells |
| --- | --- | --- |
| `high` | 7 cores — tma_id 1, 3, 5, 8, 10, 11, 13 | 79,998 |
| `low`  | 6 cores — tma_id 2, 4, 6, 9, 12, 14 | 20,199 |

13 cores total, 100,197 cells — matches `obs['tma_id']` (CLAUDE.md: the paper's "14 punch
cores, 7 low / 7 high" is wrong; id 7 is absent). All **9** cell types are present in all 13
cores, so `levels(meta$labels)` is identical across the two objects and both
`mergeCellChat` and the **functional**-similarity manifold analysis are applicable
(VC:179 — functional similarity requires the same cell-type composition).

## Call contract — one row per tutorial call

`run_cellchat.R` implements exactly this table, in this order. "Tutorial value" is the literal
argument in the vignette; where we differ the cell says so and the row is repeated under
**Deviations**.

### Stage A — input & object construction (per condition)

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| A1 | `data.input` | `Seurat::GetAssayData(seu, slot="data", assay="SCT")` — **normalized**, genes × cells | h5ad `X` (already log-normalized; raw counts live in `layers['counts']`), genes × cells `dgCMatrix` | V2:73-74 |
| A2 | `meta` | `data.frame(labels = Idents(seu), samples = "A1")` | `labels` ← `obs['cell_type']` (9), `samples` ← `obs['tma_id']`, both `factor` | V2:83-90 |
| A3 | `coordinates` | `Seurat::GetTissueCoordinates(seu, scale=NULL, cols=c("imagerow","imagecol"))` | `obsm['spatial']` (already µm), row-aligned to `data.input` columns | V2:96-99 |
| A4 | `spatial.factors` | `data.frame(ratio = 65/spot_diameter_fullres, tol = 65/2)` per section, `rownames = c("A1","A2")` | **Xenium**: `ratio = 1`, `tol = 10/2 = 5`, one row per core, `rownames = levels(meta$samples)` | VF:78-83 |
| A5 | `createCellChat(...)` | `createCellChat(object=data.input, meta=meta, group.by="labels", datatype="spatial", coordinates=spatial.locs, spatial.factors=spatial.factors)` | identical | V2:121-122 |

`spatial.factors` **must have one row per level of `meta$samples`, in `levels()` order** —
`computeRegionDistance` indexes them positionally as `ratio[k]`/`tol[k]` over
`samples.use <- levels(samples)` (`modeling.R:1165,1212`). Not documented in any vignette;
read off the source.

### Stage B — LR database

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| B1 | `CellChatDB <- CellChatDB.human` | `CellChatDB.human` (V2), `CellChatDB.mouse` (V1) | `CellChatDB.human` | V2:128 |
| B2 | `showDatabaseCategory(CellChatDB)` | same | same | V1:122 |
| B3 | `subsetDB(...)` | `subsetDB(CellChatDB, search="Secreted Signaling", key="annotation")` | **tier `default`**: identical. **tier `cellchatdb2`**: `subsetDB(CellChatDB)` — the commented-in alternative on the very next line, "use all CellChatDB except Non-protein Signaling" | V2:131, V1:127-133 |
| B4 | `cellchat@DB <- CellChatDB.use` | same | same | V2:133 |

**The bundled `CellChatDB.human` *is* CellChatDB v2** (V1:112 — "CellChatDB v2 contains ~3,300
validated molecular interactions"), which is the same resource as
`data/LRdatabase/CellChatDBv2.0.human.csv` (3,233 rows / 3,218 unique pairs; that CSV was
exported *from this R package*, CLAUDE.md). So the two tiers do **not** differ in resource —
only in **which annotation categories are used**, which is the only knob the vignettes expose.
The equivalence is audited quantitatively by `audit_db_equivalence.R`; see METHODS.md.
Re-importing the flat CSV through `updateCellChatDB` was rejected: it would *lose* the
`complex` / `cofactor` / `agonist` / `antagonist` columns that `computeCommunProb` multiplies
into `dataRavg` (`modeling.R:119-124`), i.e. it would degrade CellChat while adding nothing.

Both requested LRIs are `Secreted Signaling` and so are in **both** tiers:
`GRN → SORT1` (pathway GRN) and `ANXA1 → FPR1` (pathway ANNEXIN).

### Stage C — preprocessing

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| C1 | `subsetData(cellchat)` | same, no args | same | V2:141 |
| C2 | `future::plan("multisession", workers = 4)` | same | same | V2:142 |
| C3 | `identifyOverExpressedGenes(cellchat)` | all defaults (`do.DE=TRUE`, `do.fast=TRUE`, `thresh.pc=0`, `thresh.fc=0`, `thresh.p=0.05`, `min.cells=10`, `only.pos=TRUE`) | identical | V2:143 |
| C4 | `identifyOverExpressedInteractions(cellchat)` | V2 passes nothing → `variable.both = TRUE`; V1 passes `variable.both = F` | **`variable.both = TRUE`** (V2, the governing vignette; package default) | V2:144 |

VF:110-113 offers `identifyOverExpressedGenes(cellchat, do.DE = FALSE, min.cells = 10)` for
"a dataset with a small panel of genes". The Xenium 5K panel is **not** small in that sense —
it carries 5,119 genes and >1,000 CellChatDB pairs survive gene filtering — so the DE path
(the default) is used. The actual number of over-expressed interactions is recorded in
`run_manifest.json`; if it collapses to near zero, that is a reportable result, not a licence
to switch paths silently.

### Stage D — inference

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| D1 | `computeCommunProb(...)` | V2: `type="truncatedMean", trim=0.1, distance.use=FALSE, interaction.range=250, scale.distance=NULL, contact.dependent=TRUE, contact.range=100` | same, except **`contact.range = 10`** | V2:156-158 |
| D2 | `filterCommunication(cellchat, min.cells = 10)` | same | same | V2:164 |
| D3 | `subsetCommunication(cellchat)` | same; `slot.name="netP"` for pathway level | both, both persisted as CSV | V1:194 |
| D4 | `computeCommunProbPathway(cellchat)` | same | same | V2:173 |
| D5 | `aggregateNet(cellchat)` | same | same | V2:179 |
| D6 | `netAnalysis_computeCentrality(cellchat, slot.name = "netP")` | same | same | V2:224 |

Everything not listed is the package default: `raw.use = TRUE` (no `projectData` — both spatial
vignettes leave the PPI projection commented out), `population.size = FALSE`, `k.min = 10`,
`do.symmetric = TRUE`, `contact.dependent.forced = FALSE`, `nboot = 100`, `seed.use = 1L`,
`Kh = 0.5`, `n = 1`.

**`contact.range = 10`, not 100** — the vignettes' 100 is the 10X Visium spot centre-to-centre
distance. For single-cell resolution the authors pin the value themselves: *"Typically,
`contact.range = 10`, which is a typical human cell size… however, for low-resolution spatial
data such as 10X visium, it should be… `contact.range = 100`"* (VF:39, V1:60,
`man/computeCommunProb.Rd`). Xenium is single-cell, so 10 is the authors' value **for our
technology**. `computeCellDistance` is run to confirm the observed nearest-neighbour distance
and the number is recorded in the manifest.

**`distance.use = FALSE`** follows V2, the governing vignette. V1 (single-section) instead uses
`distance.use = TRUE, scale.distance = 0.01`; with `distance.use = TRUE` the communication
probability is scaled by `1/(d × scale.distance)` where `d` is the **cell-group-level** distance
(see Spatial model below), and `scale.distance` must be chosen so `min(d × scale.distance)`
lands in `[1,2]` or `computeCommunProb` hard-`stop()`s with the value to use
(`modeling.R:152-156`) — i.e. it is data-dependent and cannot be copied from a Visium tutorial.
A `distance.use = TRUE` variant is run as a **sensitivity check** into a separate
`sensitivity_distance/` tree with `scale.distance` picked by that rule, never mixed with the
tutorial-faithful run.

### Stage E — per-object visualization (VB + V1/V2)

| # | Call | Ours | Source |
| --- | --- | --- | --- |
| E1 | `netVisual_circle(cellchat@net$count / $weight, vertex.weight=rowSums(...), weight.scale=T, label.edge=F)` | both | V2:192-193 |
| E2 | `netVisual_heatmap(cellchat, measure="count" / "weight", color.heatmap="Blues")` | both | V2:196-197 |
| E3 | per-sender circle grid: loop rows of `cellchat@net$weight`, `edge.weight.max=max(mat)` | all 9 | VB:252-258 |
| E4 | `netVisual_aggregate(signaling=p, layout="circle")` | every pathway in `cellchat@netP$pathways` | V2:210 |
| E5 | `netVisual_aggregate(signaling=p, layout="chord")` | every pathway | VB:293 |
| E6 | `netVisual_aggregate(signaling=p, vertex.receiver=...)` — hierarchy | every pathway, `vertex.receiver` = index of the tumour cell groups | VB:287 |
| E7 | `netVisual_heatmap(signaling=p, color.heatmap="Reds")` | every pathway | VB:296 |
| E8 | `netVisual_aggregate(signaling=p, sample.use=<core>, layout="spatial", edge.width.max=2, vertex.size.max=1, alpha.image=0.2, vertex.label.cex=0)` | every pathway × **per core** (`sample.use` is mandatory in multi-sample mode) | V2:215 |
| E9 | same with `vertex.weight="incoming", vertex.size.max=6` | every pathway × per core | V2:234 |
| E10 | `netAnalysis_contribution(cellchat, signaling=p)` | every pathway | V2:240 |
| E11 | `extractEnrichedLR(cellchat, signaling=p, geneLR.return=FALSE)` → `netVisual_individual(..., pairLR.use=LR, layout="circle" / "chord")` | see "Per-LR plots" below | VB:315-323 |
| E12 | `netVisual_bubble(sources.use=…, targets.use=…, remove.isolate=FALSE)` | all-sources × all-targets | VB:351 |
| E13 | `netVisual_chord_gene(sources.use=…, targets.use=…)`, and `slot.name="netP"` | both | VB:381-387 |
| E14 | `netVisual_chord_cell(signaling=p, group=…)` | every pathway | VB:304 |
| E15 | `netAnalysis_signalingRole_network(signaling=p, width=8, height=2.5, font.size=10)` | every pathway | V2:227 |
| E16 | `netAnalysis_signalingRole_scatter(cellchat)` | aggregate + per pathway | VB:434 |
| E17 | `netAnalysis_signalingRole_heatmap(pattern="outgoing" / "incoming")` | both | VB:444-445 |
| E18 | `selectK(pattern="outgoing" / "incoming")` → `identifyCommunicationPatterns(k=…)` → `netAnalysis_river` + `netAnalysis_dot` | both patterns; `k` read off the selectK curve, recorded in the manifest | VB:474-501 |
| E19 | `computeNetSimilarity(type="functional" / "structural")` → `netEmbedding` → `netClustering` → `netVisual_embedding` | both types | VB:513-528 |
| E20 | `spatialFeaturePlot(features=…, sample.use=<core>, point.size=0.8, color.heatmap="Reds", direction=1)` | ligand+receptor genes of the requested LRIs | V2:247-248 |
| E21 | `spatialFeaturePlot(pairLR.use=…, sample.use=<core>, do.binary=FALSE/TRUE, cutoff=0.05, enriched.only=F, …)` | both continuous and binary, requested LRIs | V2:251-254 |
| E22 | `plotGeneExpression(cellchat, signaling=p, enriched.only=TRUE, type="violin")` | every pathway | VB:395 |

Not run, with reason: `runCellChatApp` (VB:542 — interactive Shiny, not a file artifact);
`netVisual_embeddingZoomIn` (VB:528 — only meaningful when a group has enough members;
attempted, skipped with a logged message when it errors); `projectData` (commented out in both
spatial vignettes).

### Stage F — cross-condition comparison, high vs low grade (VC)

| # | Call | Tutorial value | Ours | Source |
| --- | --- | --- | --- | --- |
| F1 | `mergeCellChat(object.list, add.names = names(object.list))` | `list(NL=…, LS=…)` | `list(low=…, high=…)` — **low first**, so "increased in the second dataset" reads as increased in high grade | VC:53-54 |
| F2 | `compareInteractions(cellchat, show.legend=F, group=c(1,2))` + `measure="weight"` | same | same | VC:78-79 |
| F3 | `netVisual_diffInteraction(weight.scale=T)` + `measure="weight"` | same | same | VC:91-92 |
| F4 | `netVisual_heatmap(cellchat)` + `measure="weight"` | same | same | VC:98-99 |
| F5 | `getMaxWeight(object.list, attribute=c("idents","count"))` → per-object `netVisual_circle` | same | same | VC:109-113 |
| F6 | `mergeInteractions(x, group.cellType)` → re-merge → `netVisual_circle(count.merged)` / `netVisual_diffInteraction(measure="count.merged")` | 12 clusters → 3 coarse types | 9 cell types → coarse grouping defined in the run script's `--coarse-map` (default: identity, i.e. **skipped**, since our 9 labels are already coarse) | VC:119-139 |
| F7 | `netAnalysis_signalingRole_scatter(object.list[[i]], weight.MinMax=…)` | same | same | VC:149-155 |
| F8 | `netAnalysis_signalingChanges_scatter(cellchat, idents.use=…)` | `"Inflam. DC"` | loop over **all 9** cell types | VC:163-164 |
| F9 | `computeNetSimilarityPairwise(type="functional")` → `netEmbedding` → `netClustering` → `netVisual_embeddingPairwise(label.size=3.5)` | same | same, **and** `type="structural"` (VC:199-204, `eval=FALSE` in the vignette but a documented step) | VC:189-193 |
| F10 | `rankSimilarity(cellchat, type="functional")` | same | same | VC:211 |
| F11 | `rankNet(mode="comparison", measure="weight", stacked=T/F, do.stat=TRUE)` | same | same | VC:223-224 |
| F12 | `netAnalysis_signalingRole_heatmap(object.list[[i]], pattern=…, signaling=pathway.union, width=5, height=6)` | outgoing / incoming / all | all three | VC:239-256 |
| F13 | `netVisual_bubble(comparison=c(1,2), angle.x=45)`, then `max.dataset=2` / `max.dataset=1`, `remove.isolate=T` | sources 4 → targets 5:11 | all sources → all targets | VC:267-275 |
| F14 | `identifyOverExpressedGenes(group.dataset="datasets", pos.dataset=…, features.name=…, only.pos=FALSE, thresh.pc=0.1, thresh.fc=0.05, thresh.p=0.05, group.DE.combined=FALSE)` | `pos.dataset="LS"` | `pos.dataset = "high"` | VC:293 |
| F15 | `netMappingDEG(features.name=…, variable.all=TRUE)` → `subsetCommunication(net=net, datasets=…, ligand.logFC=±0.05, receptor.logFC=NULL)` | same | same | VC:296-300 |
| F16 | `extractGeneSubsetFromPair(net.up / net.down, cellchat)` | same | same | VC:305-306 |
| F17 | `netVisual_bubble(pairLR.use=net.up[, "interaction_name", drop=F], comparison=c(1,2), angle.x=90, remove.isolate=T)` | same | same, up and down | VC:320-323 |
| F18 | `netVisual_chord_gene(object.list[[i]], slot.name='net', net=net.up/net.down, lab.cex=0.8, small.gap=3.5)` | same | same | VC:332-333 |
| F19 | `computeEnrichmentScore(net.down / net.up, species='human', variable.both=TRUE)` | same | same | VC:340-343 |
| F20 | `plotGeneExpression(cellchat, signaling=p, split.by="datasets", colors.ggplot=T, type="violin")` | `"CXCL"` | pathways of the requested LRIs (GRN, ANNEXIN) + top-ranked | VC:419 |
| F21 | `findEnrichedSignaling(object.list[[i]], features=…, idents=…, pattern="outgoing")` | `c("CCL19","CXCL12")` | the requested LRIs' ligands | VC:311 |

`thresh.fc = 0.05` at F14 is the vignette's own presto-adjusted value (VC:291) and is kept —
which means **presto must be installed**, otherwise `do.fast = TRUE` silently falls back to
`stats::wilcox.test` and returns systematically larger logFC against an unchanged threshold
(`utilities.R:434-445`). Presence of presto is asserted at startup and recorded in the manifest.

## Per-LR plots — two separate sets (skill invariant)

* `plots/top_lr/` — the LRs **CellChat itself** ranks highest: by communication probability in
  `subsetCommunication()` (`prob`, with `pval < 0.05`), and for the comparison by
  `netVisual_bubble(max.dataset=…)` / `net.up` / `net.down`.
* `plots/requested_lr/` — always **`GRN_SORT1`** and **`ANXA1_FPR1`** (the ALARMIST motif-1
  mGAM ⇄ MES-like loop), whatever their rank, at every level CellChat supports:
  `netVisual_individual` (circle + chord), `spatialFeaturePlot` (continuous + binary, per core),
  `netVisual_bubble(pairLR.use=…)`, and the pathway-level views for GRN and ANNEXIN.

CellChat joins ligand and receptor with an **underscore** in `interaction_name`
(`GRN_SORT1`, `ANXA1_FPR1`) and with a **hyphen + space** in `interaction_name_2`
(`GRN - SORT1`); `pairLR.use` for `netVisual_individual`/`netVisual_bubble` wants
`interaction_name`, `spatialFeaturePlot(pairLR.use=)` also wants `interaction_name`.

All four genes — GRN, SORT1, ANXA1, FPR1 — must be checked against the 5,119-gene Xenium panel
at runtime; if any is absent, or if the pair is dropped by
`identifyOverExpressedInteractions`, or its `pval >= 0.05`, **that is the result** and it is
written to `requested_lr/requested_lr_status.csv` rather than silently omitted.

## What CellChat actually computes

For each L-R pair and each ordered pair of **cell groups** (i → j), the communication
probability is a Hill function of the product of the group-average ligand expression in i and
group-average receptor expression in j:
`P = L·R / (Kh^n + L·R)` with `Kh = 0.5`, `n = 1` (`modeling.R`), where L and R are computed
over `data/max(data)` with a **10% truncated mean** per group (`type="truncatedMean", trim=0.1`),
multi-subunit complexes entering as the **geometric mean** across subunits
(`computeExpr_LR`), and the receptor term further multiplied by co-activation and divided by
co-inhibition receptor terms (`modeling.R:122-124`). Significance is a **permutation test**:
cell group labels are shuffled `nboot = 100` times, and `pval = #{P_boot >= P_obs}/nboot`
(`modeling.R:207,304`). So the unit of a CellChat result is a **(sender cell type, receiver
cell type, L-R pair)** triple — *not* a cell and *not* a spot. `netP` sums the probabilities
of all L-R pairs in a pathway; `aggregateNet` counts significant links (`net$count`) and sums
probabilities (`net$weight`) per cell-type pair.

## Spatial model

Spatial information enters at the **cell-group** level, not the cell level. For each sample k
and each ordered cell-group pair (i,j), `computeRegionDistance` (`modeling.R:1194-1228`) takes
every cell of group i, finds its **1-nearest neighbour** in group j (`BiocNeighbors::queryKNN`,
Annoy), converts that distance to µm by `× ratio[k]`, and takes the **10% trimmed mean** over
group i's cells → `d.spatial[i,j,k]`. Then:

* `adj.spatial[i,j,k] = 1` iff at least `k.min = 10` distinct group-j cells lie within
  `interaction.range + tol` (250 + 5 µm);
* `adj.contact[i,j,k] = 1` iff at least `k.min = 10` lie within `contact.range + tol` (10 + 5 µm);
* across samples: `d.spatial` is **averaged**, the adjacency matrices are `1` if **any** sample
  says 1 (`modeling.R:1231-1238`), then symmetrized (`adj * t(adj)`, `modeling.R:1242-1244`);
* group pairs with `adj.spatial == 0` are set to `NaN` and thereby **excluded entirely**.

With `distance.use = FALSE` (V2) `P.spatial` is an all-ones matrix zeroed at the excluded
pairs — i.e. distance acts as a **hard filter only**. With `distance.use = TRUE` (V1)
`P.spatial = 1/(d × scale.distance)`, a monotone down-weighting of distant cell-type pairs.
`contact.dependent = TRUE` restricts *only* the `Cell-Cell Contact` rows of the DB to
`adj.contact`; with the `default` (Secreted-Signaling-only) tier there are no such rows, so
CellChat prints *"Molecules of the input L-R pairs are diffusible"* and `contact.range` has no
effect — exactly as V1:174 says. It bites only in the `cellchatdb2` tier.

## Gotchas found by reading the source

1. **`spatial.factors` is indexed positionally by `levels(meta$samples)`** — one row per sample,
   in level order. Wrong order silently applies the wrong µm conversion.
2. **`contact.range` is mandatory when `contact.dependent = TRUE`** — `computeRegionDistance`
   `stop()`s if both `contact.range` and `contact.knn.k` are NULL (`modeling.R:1187-1189`).
3. **`scale.distance` is validated, not defaulted** — with `distance.use = TRUE`, if
   `min(d × scale.distance) < 1` CellChat aborts and prints the value to use
   (`modeling.R:152-156`). Never copy 0.01 from the Visium tutorial.
4. **Unused factor levels abort the run** — `computeCommunProb` `stop()`s unless
   `nlevels(idents) == length(unique(idents))` (`modeling.R:105-108`). Subsetting to one grade
   leaves stale `cell_type` levels; `droplevels` is required.
5. **presto changes the numbers, not just the speed.** `do.fast = TRUE` is the default and
   silently falls back to `stats` if presto is missing, giving larger logFC while VC's
   `thresh.fc = 0.05` was tuned for presto (`utilities.R:434-445`, VC:291).
6. **`sample.use` is mandatory for spatial plots in multi-sample mode** — `netVisual_aggregate(layout="spatial")`
   and `spatialFeaturePlot` plot one section at a time (V2:215, V2:243).
7. **`netEmbedding` prefers python `umap-learn` via reticulate** and only falls back to `uwot`
   when told (`analysis.R:652-653`). We pass `umap.method = "uwot"` explicitly so the run does
   not depend on a reticulate python.
8. **The permutation test densifies.** `aggregate(t(data.use), list(group), FUN=FunMean)` runs
   once per bootstrap over the signaling-gene × cell matrix; with `workers = 4` that is 4
   concurrent dense copies. Memory scales with `n_signaling_genes × n_cells`, so the `high`
   object (79,998 cells) is the memory high-water mark, not the number of LR pairs.
9. **`filterCommunication(min.cells = 10)` is applied per cell group across the whole object**,
   not per sample.

## Deviations from the tutorial

| Item | Tutorial | Ours | Why |
| --- | --- | --- | --- |
| `contact.range` | `100` (V1:182, V2:158) | **`10`** | 100 is the 10X Visium spot pitch. VF:39 / V1:60 / `man/computeCommunProb.Rd` pin `10` for **single-cell-resolution** platforms; Xenium is single-cell. `computeCellDistance` output recorded in the manifest. |
| `spatial.factors` | Visium `ratio = 65/spot_diameter_fullres`, `tol = 32.5` | `ratio = 1`, `tol = 5` | VF:78-83, the FAQ's **Xenium** row: coordinates already in µm, `spot.size = 10` (typical human cell). |
| normalization | `GetAssayData(slot="data", assay="SCT")` | h5ad `X`, used as-is | `X` is already log-normalized (CLAUDE.md); `layers['counts']` holds the raw counts. No re-normalization — `normalizeData` would log a second time. |
| `variable.both` | V1 passes `F`, V2 passes nothing (`TRUE`) | `TRUE` | V2 is the governing vignette and `TRUE` is the package default. |
| `distance.use` | V1 `TRUE`+`scale.distance=0.01`; V2 `FALSE`+`NULL` | `FALSE` (main run), `TRUE` as a separate sensitivity tree | V2 governs our data shape; V1's `0.01` is Visium-specific and would abort (gotcha 3). |
| `umap.method` | not passed (→ `umap-learn`) | `"uwot"` | avoids a reticulate python dependency; `uwot` is the package's own documented alternative. |
| cell groups in `sources.use`/`targets.use` | hard-coded indices (`4`, `5:11`) | all 9 cell types | the indices are specific to the skin dataset's 12 clusters. |
| `pathways.show` | one hand-picked pathway (`"IGF"`, `"EGF"`, `"CXCL"`) | **every** pathway in `netP$pathways` | skill invariant: generate every plot the standard workflow can produce. |
| `group.cellType` / `mergeInteractions` (F6) | 12 clusters → 3 coarse types | skipped by default | our 9 labels are already coarse; no biologically defensible 3-way grouping without asking. |
| tier `cellchatdb2` | — | `subsetDB(CellChatDB)` (all but Non-protein Signaling) | the bundled DB **is** CellChatDB v2, so the tier varies annotation scope rather than resource; re-importing our flat CSV would drop complex/cofactor/agonist/antagonist columns. Audited by `audit_db_equivalence.R`. |
