# Comparator methods — living document

One section per cell–cell communication method benchmarked against ALARMIST, written to the
template in `.claude/skills/comparator-benchmark/SKILL.md`. Purpose of this phase: record
**what each method computes, what it eats, and what it emits** — not to force a common unit.

Datasets: **GBM** (`data/xenium_mm_final_cell_id.h5ad`, LGG/GBM TMA, 13 cores, grade high/low)
and **LUAD** (`data/linghua/P{17,21}_{AIS,LUAD}_Xenium.h5ad`, AIS vs LUAD). See CLAUDE.md.
Tiers: `default` = the method's own tutorial LR resource; `cellchatdb2` = CellChatDB v2.0 human.

| Method | Language | Version | Env | Status |
|---|---|---|---|---|
| CytoSignal | R | 0.5.1 | `comp-cytosignal` | run (GBM full; LUAD crop only — full LUAD blocked) |
| stLearn | Python | 1.4.1 | `comp-stlearn` | ✅ GBM/`cellchatdb2` done (45.4 min) |
| SpatialDM | Python | 0.3.1 | `comp-spatialdm` | ✅ GBM/`cellchatdb2` done, 13 cores (7.1 min) |
| COMMOT | Python | 0.0.3 | `comp-commot` | ✅ GBM/`cellchatdb2` done, 11/13 cores (123.4 min) |
| LIANA+ | Python | 1.8.1 | `comp-liana` | ✅ GBM **both tiers**: `cellchatdb2` bivariate (2.0 min) + inflow (1.0 min) + NMF on both; `default` (LIANA consensus) bivariate (5.6 min) + inflow (1.3 min) + NMF on both. Plus 4 further branches of the authors' decision tree — Moran's-R local metric (3.6 min), **MOFA-Flex** on inflow (76.0 min), **LRIC / cross-PCF** per punch (2.6 min), **LR-MISTy** whole slide (2.9 min) — and (2026-08-06) a **reachability-normalised MOFA-Flex re-fit** (41.5 min) that keeps all 9 sender views |
| NICHES | R | 1.2.4 | `comp-niches` | ✅ GBM/`cellchatdb2` done, 13 cores × 2 imputation sub-runs (3.9 + 9.2 min) |
| CellChat | R | 2.2.0.9001 | `comp-cellchat` | not started — default DB is already CellChatDB v2, tiers may collapse |

## ⚠️ The methods do NOT share a spatial support — counts are not comparable

The GBM data is a **13-core TMA**, and the four methods run so far were applied at different
spatial granularity. This is deliberate and documented, but it means **any table placing their
LR-pair counts side by side is misleading unless normalised.**

| Method | Spatial support | Why |
|---|---|---|
| CytoSignal | **whole TMA** — one run over 100,197 cells | no constraint forcing a split; the *grade differential* separately builds 13 per-core objects and merges them via `mergeCytoSignal`, which is CytoSignal's own documented multi-sample mode |
| stLearn | **whole TMA** — one 12,562-spot grid | no constraint forcing a split |
| SpatialDM | **13 separate runs** | **Not** to prevent cross-core weights — measured, that is a non-issue here (see below). The binding reason is that `spatialdm_global` returns **one** Moran's R per LR pair for whatever object it is given, so a pooled run yields one tissue-wide number instead of 13 per-core ones, and `diff_utils.concat_obj` takes a **list of separately fitted objects** — the grade differential is not constructible from a pooled fit. Matches the authors' own multi-sample tutorial (`differential_test_intestine.ipynb` loops `for adata in samples:` running the full pipeline per sample) |
| COMMOT | **11 separate runs** (2 cores unanalysable) | **forced**: COMMOT materialises a dense N×N distance matrix — 80.3 GB at whole-slide scale |
| NICHES | **13 separate runs** | **forced**, same reason as COMMOT: `ComputeEdgelist` builds a dense N×N distance matrix plus three more N×N copies (~80 GB for the distance matrix alone whole-slide), and its fast `nn.method` path is commented-out dead code. Also matches the split-run-merge pattern its own multi-sample vignettes (04, 07) prescribe |

**Neither method's primary tutorial prescribes per-core splitting** — all of them are single
contiguous sections. Only SpatialDM's secondary multi-sample tutorial demonstrates the pattern.

### Cross-core contamination audit — measured per method, in each method's own units

The minimum **cell-to-cell** distance between different cores is **222.9 µm**. That figure bounds
*cell-level* graphs only — **it does not bound a gridded method**, because gridding relocates
cells to bin centres and spot centroids can end up closer than any two real cells. Each method
was therefore audited in the units it actually operates on:

| Method | Graph unit | Neighbourhood | Cross-core | Impact on results |
|---|---|---|---|---|
| CytoSignal — diffusion | cell | 200 µm ε-ball | **0 pairs** | none |
| CytoSignal — contact | cell | Delaunay, **pruned at `r.diffuse.scale` = 200 µm** | **0 edges** | none — but see below |
| stLearn | **51.3 µm spot** | `distance = 250` µm | **7 spot pairs** (0.0020% of 344,370); 10 of 12,562 spots affected | **40 of 45,065 significant spot-LR calls (0.0888%)** involve a contaminated spot |
| SpatialDM | cell | 135 µm, **per core** | n/a — split | none |
| COMMOT | cell | 365 µm, **per core** | n/a — split | none |
| LIANA+ — bivariate / inflow | cell | **28.2 µm** (Gaussian support; bandwidth 13.1454, cutoff 0.1) | **0 pairs** | none |
| LIANA+ — **LR-MISTy** | cell | **607.0 µm** nominal support (`misty_rbf`, bandwidth **200** = the tutorial's, cutoff 0.01) but **capped at 100 nearest neighbours, and the cap binds for 99.7% of cells** | **2,520 edges of 10,119,190 (0.0249%)**, touching **232 of 100,190 cells (0.23%)** | small but **not zero** — see the LR-MISTy section |
| LIANA+ — LRIC / cross-PCF | cell | annuli to 225 µm, **per punch** | n/a — split | none |
| NICHES | cell | **mutual kNN, k=4 — no radius at all**; measured *post hoc*: median **10.1 µm**, p95 26.6 µm, per-core median 8.1–22.2 µm | n/a — split | none. But note the neighbourhood **shrinks as density rises**, so cores are not compared at a common physical scale, and **23.9% of its edges are self-edges** (each cell is its own neighbour) |

**Two whole-slide runs have genuine cross-core edges: stLearn, and (since 2026-08-04) LIANA+'s
LR-MISTy.** stLearn's `distance = 250`
exceeds the 222.9 µm cell floor, and the gridding makes it worse: the shortest cross-core *spot*
pair is **205.4 µm**, below the cell-level minimum. Affected core pairs are (9,11), (10,12) and
(11,14); the worst-hit LR pairs are JAM3_JAM3 (2 of 457 significant spots), C3_C3AR1 (1 of 1,244)
and MPZ_MPZL1 (1 of 111). At **0.0888%** of all significant calls this does not change any
ranking, but it is **not zero** and should be stated rather than waved away.

**LR-MISTy is the second, and it is the *tutorial's own* parameter that does it.** `lrMistyData`'s
`bandwidth=200, cutoff=0.01` gives a nominal support radius of 200 × √(−2·ln 0.01) = **607.0 µm**,
2.7× the 222.9 µm inter-core floor. What keeps the damage small is `spatial_neighbors`'
`max_neighbours=100` default — which `lrMistyData` does **not expose** — so the graph is in
practice a 100-nearest-neighbour graph (degree median = max = 100, min 16, cap binding for 99.7%
of cells). Measured by rebuilding the graph and cross-tabulating `obs['tma_id']`:
**2,520 cross-core edges of 10,119,190 (0.0249%)**, involving **232 cells**, concentrated in the
core pairs (4,6) 787 edges, (12,11) 575, (5,6) 461, (12,10) 299, (11,14) 292. Same order as
stLearn's contamination and the same verdict: too small to move a ranking, too real to omit.
**Note this contradicts the "0 pairs" row that LIANA+ carries for its bivariate/inflow branches** —
that figure is specific to the 28.2 µm kernel and does not transfer to LIANA's other entry points.

**CytoSignal is clean only by a 23 µm margin.** Delaunay triangulation is *unbounded* — the raw
triangulation on these coordinates has **716 cross-core edges, the longest 12.3 mm**, spanning
the whole TMA. `findNNDT` removes them via `nn.valid <- dist <= max.r` with
`max.r = r.diffuse.scale = r.eps.real / scale.factor = 200 µm`, and the shortest cross-core DT
edge is 223 µm. **Had `r.eps.real` been set to 250 rather than 200, CytoSignal's contact-dependent
slot would silently have contained cross-patient "touching cell" edges.** Anyone re-running this
on TMA data with a larger radius must re-check.

**Consequences to respect when reading this document:**

1. Per-core methods report pairs *per core* (SpatialDM 1,133–1,661; COMMOT 51–217; NICHES a
   fixed 1,088 scored, of which 112–702 non-imputed / 529–1,088 imputed are actually detected),
   whole-TMA methods report a single number (stLearn 526; CytoSignal 1,088). **These are
   different quantities.** Compare within a method across cores, or compare *fractions* / ranks
   across methods — never raw counts.
2. Per-core analyses have smaller n and therefore lower power per test, and the LR-pair set
   differs per core because each core re-runs its own expression filter.
3. Harmonising would only be possible in the per-core direction, since COMMOT cannot run
   whole-slide at all. That was considered and deliberately **not** done; the difference is
   recorded here instead.

---

## CytoSignal — R, v0.5.1, env `comp-cytosignal`

Welch lab (Liu & Wang). Tutorials: `results/comparators/cytosignal/reference_notebook/*.html`
(4 vignettes: main workflow, differential multi-dataset, custom LR DB, container conversion).
There is **no local git clone** — the vendored HTML is the only source of truth.

### Core algorithm

CytoSignal is a **nonparametric, per-cell permutation test for ligand–receptor activity**. For
each interaction and *each individual cell* it asks: is the local ligand–receptor co-expression
higher than expected by chance? The pipeline is four steps:

1. **Neighborhood** — each cell gets a spatial neighborhood, defined differently for diffusible
   and contact-dependent interactions.
2. **Imputation** — how much ligand (`L`) and receptor (`R`) each cell *receives* from its
   neighbors, as a distance-weighted sum over that neighborhood.
3. **LRscore** — `L × R` per cell, then averaged over the cell's Delaunay neighborhood (this
   averaging is what the `_smooth` suffix on every score slot means).
4. **Spatial permutation test** — the null is built by randomizing **spatial locations**, giving
   a one-sided p-value per (cell × interaction), then a **spatial FDR** correction that accounts
   for cell-density differences.

Two properties matter for the comparison. First, **the unit of inference is the single cell**,
not a cell-type pair and not a spot — the output is a cell × interaction score matrix plus, for
each interaction, the set of cells where it is significant. Cell-type labels are used only for
colouring plots and (optionally) as NEBULA covariates; they do **not** enter the scoring. Second,
**diffusion and contact interactions are modelled separately end-to-end**, with different
kernels and separate result slots — an interaction is classified by the database, not inferred.
Multi-subunit complexes are handled natively (subunits enter as `L1..Ln` / `R1..Rn` columns).

Interactions are ranked either by the number of significant cells, or by **SPARK-X spatial
variability** of the LRscore (`rankIntrSpatialVar`, the `result.spx` tier).

### Spatial model

| Mode | Neighborhood | Weighting | Slot |
|---|---|---|---|
| Diffusion | Epsilon-ball, radius `r.eps.real` = **200 µm** (default) | Gaussian kernel on physical distance | `GauEps` → `diffusion-Raw` |
| Contact | Delaunay triangulation (immediate neighbors) | `dt.mode = "weight_sum_2"` | `DT` → `contact-Raw` |
| Same-spot | none (raw expression, no spatial smoothing) | — | `Raw-Raw` (for Visium-like multi-cell spots) |

`inferEpsParams(scale.factor, r.eps.real = 200, thresh = 0.001)` converts the physical radius into
kernel parameters (`eps`, `sigma`). **`scale.factor` = µm per coordinate unit and is external
knowledge** — the single most dangerous parameter in this method, because a wrong value silently
rescales every neighborhood. The tutorial uses 0.73 (Slide-tags) and 0.72 (Visium). Xenium
coordinates are already in µm, so **`scale.factor = 1`** and the epsilon ball is literally 200 µm.

The `Raw-Raw_smooth` slot is scored and tested regardless of platform; for single-cell-resolution
data like Xenium it is not the intended readout (it exists for multi-cell spots) but it is still
written out.

### LR database

**`default` — bundled CellPhoneDB v2** (re-sorted by the authors), loaded by
`addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)`. Verified in the installed package:
`g_to_u` 977 genes, `inter.index` 1,396 interactions, **754 diffusion + 109 contact** unique
interactions. The DB works in a **UniProt ID space**; `changeUniprot()` rewrites the gene-symbol
expression matrix into it and drops genes absent from the DB.

**`cellchatdb2` — CellChatDB v2.0 human**, built by `scripts/comparators/cytosignal/build_cellchat_db.R`
using the package's own exported `formatLRDB(interaction_type, ligands, receptors, gene_to_uniprot)`
(the documented custom-DB route). 3,233 CSV rows → **2,683 diffusion + 535 contact**, 1,383 LR
genes, **865 of them on the Xenium 5K panel**. Two decisions:

- **Gene symbols are used as the protein-ID space** (identity `g_to_u`: `gene_name = uniprot = symbol`).
  No UniProt mapping is needed, and it matches the already-gene-symbol expression matrix. This is
  legitimate — CytoSignal only requires the ID space to be internally consistent.
- **`signaling_type` → interaction type:** `Cell-Cell Contact` → contact-dependent; `Secreted
  Signaling` / `ECM-Receptor` / `Non-protein Signaling` → diffusion-dependent. This mirrors
  ALARMIST's juxtacrine-vs-secreted split.

### Input

| Argument | Form | Notes |
|---|---|---|
| `raw.data` | genes × cells `dgCMatrix` | **raw integer counts**, rownames = UPPERCASE gene symbols |
| `cells.loc` | cells × 2 numeric matrix | colnames must be lowercase `x`, `y`; rownames = barcodes |
| `clusters` | named `factor` | names must match `colnames(raw.data)` exactly |
| `scale.factor` | scalar | µm per coordinate unit — **not in the data**, supplied by the user |

Our on-disk exchange format (`results/comparators/cytosignal/<dataset>/input*/` — deliberately
*outside* the tier dirs, since the export does not depend on the LR database):
`counts.mtx` (genes × cells), `genes.tsv`, `barcodes.tsv`, `meta.csv` (`cell_id,x,y,celltype`),
`provenance.json`. GBM additionally has `meta_grade.csv` (`+tma_id,grade,patient`) for the
multi-sample test. Converters also exist in the other direction: `SeuratToCS`, `SCEToCS`.

### Workflow

| # | Call (tutorial argument values) | Produces |
|---|---|---|
| 1 | `createCytoSignal(raw.data, cells.loc, clusters)` | `CytoSignal` object |
| 2 | `addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)` | LR DB attached |
| 3 | `removeLowQuality(cs, counts.thresh = 300, gene.thresh = 50)` | QC-filtered cells/genes |
| 4 | `changeUniprot(cs)` | expression re-indexed into the DB's ID space |
| 5 | `inferEpsParams(cs, scale.factor = <µm/unit>, r.eps.real = 200)` | `eps`, `sigma` |
| 6 | `findNN(cs)` | `GauEps` + `DT` neighbor graphs |
| 7 | `imputeLR(cs)` | imputed L and R per cell |
| 8 | `inferIntrScore(cs, perm.size = 1e5)` (seed first) | LRscores + permutation null + p-values |
| 9 | `inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100)` | `result`, `result.hq` |
| 10 | `rankIntrSpatialVar(cs)` | `result.spx` (SPARK-X ranked) |
| 11 | `showIntr(cs, slot.use, signif.use, return.name = TRUE)` | significant interaction list |

`reads.thresh` = minimum reads for an interaction to be considered; `sig.thresh` = minimum number
of significant cells. Optional `recep.smooth = TRUE` in step 8 adds `diffusion-DT` / `contact-DT`
slots (receptor also DT-smoothed) — for sparse data; we did not use it.

### Data outputs

Object slots (`@lrscore[[slot]]`), for each of `diffusion-Raw_smooth`, `contact-Raw_smooth`,
`Raw-Raw_smooth` (and non-smooth `*-Raw` variants):

| Slot | Shape | Meaning |
|---|---|---|
| `@score` | cells × interactions | LRscore. **Near-dense** — the dominant memory/disk cost |
| `@res.list$result` | list per interaction | barcodes of cells with p < `p.thresh` |
| `@res.list$result.hq` | list per interaction | + passes `reads.thresh` / `sig.thresh` QC |
| `@res.list$result.spx` | list per interaction | + spatially variable by SPARK-X (the headline tier) |

Persisted by our `quant_io.R` into `<run>/quant/` (names alone are useless for comparison):

| File | Content |
|---|---|
| `score_<slot>.mtx` + `.cells.tsv` + `.intr.tsv` | cells × interactions LRscore (≤100k cells) |
| `score_<slot>.rds` | same as `dgCMatrix` with dimnames (>100k cells) |
| `reslist_<slot>.rds` | full `@res.list` (`result` / `result.hq` / `result.spx`) |
| `signif_summary_<slot>.csv` | `interaction_id, name, n_result, n_hq, n_spx` |

Multi-dataset: `@diff.results$<covariate_level>` → `interaction, logFC, se, p, padj`.
Signaling-associated DE: `inferIntrDEG(...)[[i]]$sign_genes`.
Export routes: `csToSeurat`, `csToSCE` (LRscores as `altExp`), `csToH5AD` (LRscores + significance
into `.obsm`, interaction metadata into `.uns['CytoSignal']`).

### Image outputs — full stock inventory

| Function | Shows | Used? |
|---|---|---|
| `plotCluster(cs)` | cell types in space (legend reference) | ✅ `cluster_map.png` both runs |
| `plotSignif(cs, intr, slot.use, signif.use, plot_dir, plot.fmt)` | per-interaction combination panel: imputed L & R, raw L & R, LRscore, cluster annotation, optional 3D edge. Writes `Rank_<n>_<L>-<R>.<fmt>` | ✅ LUAD crop (top 6 per mode); ✅ **GBM** (top 6 per mode in `plots/signif_<slot>/`, + the two requested LRIs in `plots/requested_diffusion_Raw_smooth/`) |
| `plotEdge(cs, intr, slot.use)` | 3D sender→receiver edge plot (senders bottom, receivers top, lines = edges) | ❌ not run |
| `plotNebulaVolcano(multics, covariate, intr.type, fdrThresh, logfcThresh)` | volcano of the multi-dataset test | ❌ unusable, see deviations |
| `plotNebulaAll(multics, intr.type)` | all-interaction summary of the multi-dataset test | ❌ unusable, see deviations |
| `heatmap_GO(intrDEG, goRes$result, intr, ...)` | GO enrichment of signaling-associated DEGs | ❌ not run |
| `plotREVIGO(revigoRes, labelSize)` | REVIGO semantic-similarity scatter of those GO terms | ❌ not run |
| via `csToSeurat` / `csToSCE` | `SpatialFeaturePlot`, `SpatialDimPlot`, `plotReducedDim` on LRscore / significance layers | ❌ not run |

Everything else in `GBM/cellchatdb2/run_full/plots/` (`ranking.png`, `top_{diffusion,contact}_grid.png`,
`spatial_scores.png`, `lr_panels_*.png`, `comparison_mGAM_vs_cytosignal.png`,
`reconstruct_motif1_from_cytosignal.png`, `motif1_top25_lris_cytosignal.png`,
`grade_comparison_2panel.png`) is **custom ALARMIST-comparison plotting, not stock CytoSignal.**

### Multi-sample / differential mode

Native and documented: build one CytoSignal object per sample (each prepped through
`findNN` + `imputeLR` — the merge needs the DT imputation), then

```r
multics <- mergeCytoSignal(objList, metadata = dataset.meta, name.by = "sample")
multics <- runNEBULA(multics, covariates = c('clusters', 'age'), cpc_thresh = 0.001, ncore = 4)
```

`runNEBULA` fits a **negative-binomial mixed model** per interaction (`nebula`), with dataset as
the subject/random effect and total counts as offset. Tutorial default in the signature is
`cpc_thresh = 0.005`, but the vignette body uses **0.001** — lower keeps *more* interactions.

For GBM we used **`tma_id` (13 cores) as the sample unit and `grade` as the covariate** —
correct per the sample-is-the-replicate rule, and grade is constant within a core.

### Gotchas

- **`removeLowQuality` default `counts.thresh = 300` drops ~half of a Xenium 5K panel.** Targeted
  panels have far lower per-cell totals than whole-transcriptome data.
- **`showIntr(return.name = TRUE)` errors on a custom DB.** `formatLRDB` writes only 3 columns to
  `intr.index`, but `getLigandNames`/`getReceptorNames` read columns 4–5 (`protein_name_a/b`).
  `build_cellchat_db.R` pads them with empty strings so naming falls back to `partner_a/b`.
- **`formatLRDB` returns `db.diff`/`db.cont` in a different column order than the bundled DB**
  (`ligands, receptors, combined` vs `combined, ligands, receptors`), and `getIntrValue` inside
  `plotSignif` reads them **by position** — so a raw custom DB swaps the ligand/receptor labels in
  plots. `build_cellchat_db.R` reorders to match. Scoring reads by name and is unaffected.
- **`@score` is effectively dense.** 89k cells × 919 interactions ≈ 1.9 GB as `.mtx`; whole objects
  are 2.7–4.3 GB per ~25k cells. `purgeBeforeSave()` keeps `@score`/`@res.list` and clears only the
  imputation/raw slots — it is not a small-file escape hatch.
- **`nosave` runs cannot be replotted — only recomputed.** Every plotting function needs
  `@score` + `@imputation` + `@counts` + `@cells.loc`, none of which survive a disk-safe run.
  Deciding to skip a plot at run time therefore costs a **full pipeline re-run** later
  (`plot_signif_rerun.R`, ~17 min for GBM), not a cheap replot. Plan the plot list up front.
- **`perm.size` is floored at `n_cells`.** At 498k cells the null matrix is enormous; this is what
  OOMs the full LUAD section (~57 GB peak on a 36 GB box).
- **`nebula` will not build in the conda env's R 4.3.3** (C++ vs Eigen 3.4/lgamma; CRAN needs R ≥ 4.4),
  while `cytosignal` will not link in system R 4.4.2 (Fortran/BLAS). See deviations.
- **`SPARK` is not auto-installed** despite the `Remotes` field — `rankIntrSpatialVar` fails without
  it, and `plotSignif(raster = TRUE)` needs `scattermore`.
- Degraded conda on this box: `conda activate`/`conda run` fall through to system R. Always
  `source scripts/comparators/cytosignal/activate_env.sh`.

### Deviations from the tutorial

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `removeLowQuality` | `counts.thresh = 300`, `gene.thresh = 50` | **100 / 20** | 300 drops ~half the cells on a 5K targeted panel. GBM: 100,197 → 89,035 cells. Forced by panel size, not preference. |
| `scale.factor` | 0.73 (Slide-tags) | **1** | Xenium coordinates are already µm. Required for a 200 µm ball to mean 200 µm. |
| `numCores` | 1 | 4 | Performance only; no effect on results. |
| `runNEBULA` | single call in the merged object | **reimplemented as a 2-stage script** (`run_nebula_grade.R`): stage 1 in conda R builds/merges the 13 core objects and extracts the exact model inputs via `cytosignal:::.setup.model(merged, "grade")`; stage 2 re-invokes itself in **system R 4.4.2** to call `nebula::nebula` on those inputs | `nebula` and `cytosignal` cannot be installed in the same R. Uses CytoSignal's own model-setup internals so the fitted model is the one `runNEBULA` would fit. **Uses a private `:::` function — the one place we depend on internals.** |
| NEBULA covariates | `c('clusters', 'age')` | **`c('grade')` only** | We test grade, and cell-type composition differs systematically by grade — including `clusters` would absorb the effect of interest. Flagging as a real modelling choice, not a technicality. |
| `plotNebulaVolcano` / `plotNebulaAll` | used | **not usable** — replaced by `build_grade_2panel.py` | The split-stage route never populates `multics@diff.results`, which both functions read. |
| p-value adjustment | internal to `runNEBULA` | **BH within each mode** (diff and cont separately) | Matches `runNEBULA` internals; pooling the two modes would be wrong. |
| `cpc_thresh` | 0.005 (signature) / 0.001 (vignette body) | **0.001** | Followed the vignette. |
| `plotSignif` for GBM | top 5 per mode | **top 6 per mode**, produced by a full pipeline re-run (`plot_signif_rerun.R`) | The GBM run was `nosave`, so no object survived to replot from. 6 (not the tutorial's 5) to match the LUAD crop run so the two datasets are comparable. The re-run was **verified to reproduce the original exactly** — see below. |
| requested LRIs | n/a | **GRN→SORT1 and ANXA1→FPR1 additionally plotted**, into `plots/requested_<slot>/` | User-specified LRIs of interest (ALARMIST motif-1 mGAM loop). Kept in a separate directory so they are never confused with CytoSignal's own ranking. |
| `Raw-Raw_smooth` plots | tutorial plots whichever slot | **skipped** | That slot is the multi-cell-spot (Visium-like) readout; it is not the intended output for single-cell-resolution Xenium. Scores are still persisted in `quant/`. |
| GO / REVIGO / `inferIntrDEG` | in the tutorial | not run | Downstream interpretation, out of scope for this phase. |

### Runs on our data

Migrated to the `<dataset>/<tier>/` convention on 2026-07-31 (was `GBM/run_full`,
`GBM/nebula_grade`, `LUAD/run_crop_2p5mm`). Layout under `results/comparators/cytosignal/`:

```
GBM/   input_full/                  input export — outside the tiers (DB-independent)
       default/     STATUS.md       not run (bundled CellPhoneDB v2)
       cellchatdb2/ run_full/       the GBM run
                    nebula_grade/   the high-vs-low grade test
LUAD/  input/  input_full/          2.5 mm crop and full-section exports
       default/     run_crop_2p5mm/ the LUAD run (crop only)
       cellchatdb2/ STATUS.md       blocked, needs >= 64 GB
bundle_bignode/                     portable bundle for the blocked LUAD cellchatdb2 run
reference_notebook/                 vendored tutorials      cellchat_db_human.rds  shared DB
```

The run-name subdirectory is kept inside the tier because it carries information the tier does
not (which section, full vs crop, which analysis) and leaves room for e.g. `LUAD/default/run_full/`.
`bundle_bignode/` holds `.R`/`.py` inside `results/` — a sanctioned exception to the
code/output separation, because it is a self-contained bundle meant to be copied to another machine.

| Dataset | Tier | Path (under `results/comparators/cytosignal/`) | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `GBM/cellchatdb2/run_full/` | ✅ done | 100,197 → **89,035** cells after QC; **919 diffusion / 169 contact** interactions scored; **895 / 166** with ≥1 significant cell (151 in `Raw-Raw_smooth`). Top diffusion by significant-cell count: WNT3–FZD*/LRP6 (~44k cells); top contact: OCLN–OCLN, CADM3–CADM3, PTPRM–PTPRM. **Stock `plotSignif` figures added 2026-07-31** (12 top-ranked + 2 requested) |
| GBM | `default` | `GBM/default/` (STATUS.md) | ❌ **not run** | bundled CellPhoneDB v2 tier missing; command in the STATUS file |
| GBM | grade test | `GBM/cellchatdb2/nebula_grade/` | ✅ done | 13 cores; **659 interactions** tested (535 diff / 124 cont); 147 raw p<0.05; **only 4 FDR-significant**, all generic junction/adhesion: CDH2–CDH2 (+0.53), NECTIN3–NECTIN1 (−0.63), F11R–F11R (−0.88), JAM3–F11R (−1.10). Of ALARMIST motif-1's top-100 LRIs (95 testable) only **1** (JAM3–F11R) is significant. ⚠️ **See the density–grade confound below** |

> ⚠️ **Confound affecting every cross-grade result in this document.** Core cell density
> correlates with grade: **r = 0.659, p = 0.014** (point-biserial, n = 13 cores; Mann-Whitney
> p = 0.005). Median cells within a 134.6 µm radius is **243 in high-grade cores (78–532) vs 60
> in low-grade (35–118)** — a **4× difference**. High cellularity is itself a WHO glioma grading
> criterion, so this is real biology rather than a pipeline artifact, but it means statistical
> power is not constant across grades for any method whose neighbourhood is distance-based —
> CytoSignal's 200 µm epsilon ball included. Differences in *counts* of significant interactions
> between grades are therefore partly a power effect and must not be read as pure biology.
> Also note the 13 cores come from only **7 patients** (patient 14007 contributes 4), so
> core-level tests are pseudoreplicated at the patient level.
| LUAD | `default` | `LUAD/default/run_crop_2p5mm/` | ⚠️ **crop only** | P21_LUAD, 2.5 mm window, 28,596 cells. **277 diffusion / 44 contact** significant (`result.spx`). Top diffusion: TGFB1/TGFB3–TGFβR1/R2, PDGFB–PDGFR; top contact: αMβ2/αLβ2/αXβ2–ICAM1, JAG1–NOTCH1, DLL4–NOTCH4 |
| LUAD | `cellchatdb2` | `LUAD/cellchatdb2/` (STATUS.md) → `bundle_bignode/` | ❌ **blocked** | Full P21 (560,183 → 498,422 cells) OOMs at ~57 GB. Self-contained bundle prepared for a ≥64 GB node |
| LUAD | AIS vs LUAD | — | ❌ **not run** | Only P21_LUAD was ever touched. The AIS↔LUAD comparison across `P{17,21}` has not been attempted |

### Requested LRIs — where CytoSignal puts the ALARMIST motif-1 loop (GBM, `cellchatdb2`)

Both are present in the DB, both are called significant, and both are **diffusion-only** —
consistent with GRN and ANXA1 being secreted ligands, and with CytoSignal classifying
interactions by database annotation rather than inferring the mode.

| LRI | ID | slot | rank (of 895 signif.) | significant cells | tiers |
|---|---|---|---|---|---|
| GRN → SORT1 | `CCI-01109` | `diffusion-Raw_smooth` | **66** | 27,630 / 89,035 (31%) | result = hq = spx |
| ANXA1 → FPR1 | `CCI-01088` | `diffusion-Raw_smooth` | **255** | 15,292 / 89,035 (17%) | result = hq = spx |
| both | — | `contact-Raw_smooth`, `Raw-Raw_smooth` | — | not significant | — |

`result = hq = spx` means each passes the p-value threshold, the read/bead QC, **and** SPARK-X
spatial variability — so CytoSignal considers both spatially structured, not diffuse background.
Neither is anywhere near CytoSignal's own top of the ranking, which is dominated by WNT3–FZD*/LRP6
(diffusion) and homotypic junction pairs (contact). Figures:
`plots/requested_diffusion_Raw_smooth/Rank_{66_GRN-SORT1,253_ANXA1-FPR1}.png`.

### Reproducing the GBM run

`plot_signif_rerun.R` rebuilds the object from `input_full/` with run_cytosignal.R's exact
parameters and **checks itself against the stored `quant/` before plotting**. Result
(`plots/reproduction_check.csv`) — every interaction in every slot matched exactly:

| slot | interactions | `n_hq` exact matches | corr | `n_spx` exact matches | corr |
|---|---|---|---|---|---|
| `diffusion-Raw_smooth` | 895 | 895 | 1.00 | 895 | 1.00 |
| `contact-Raw_smooth` | 166 | 166 | 1.00 | 166 | 1.00 |
| `Raw-Raw_smooth` | 151 | 151 | 1.00 | 151 | 1.00 |

So CytoSignal is **fully deterministic** here under `set.seed(42)` even with `numCores = 4` —
worth knowing, since it means a lost object is always recoverable at ~27 min (whole pipeline,
all slots ranked) or ~12 min (one slot). Cost: no checkpoint means recompute, never replot.

### Methods paragraph

> For CytoSignal (v0.5.1), we analyzed the data according to the default workflow, which consists
> of (1) defining spatial neighborhoods, (2) imputing ligand and receptor expression, (3)
> calculating LRscores with a spatial permutation test and (4) identifying significant
> interactions. We used the functions `inferEpsParams`, `findNN`, `imputeLR`, `inferIntrScore`,
> `inferSignif` and `rankIntrSpatialVar`, retaining the authors' default parameters
> (`r.eps.real = 200` µm, `perm.size = 1e5`, `p.thresh = 0.05`, `reads.thresh = 100`,
> `sig.thresh = 100`). Because Xenium coordinates are reported in microns, we set
> `scale.factor = 1`; because the 5,000-plex targeted panel yields lower per-cell counts than the
> whole-transcriptome data used in the tutorial, we relaxed `removeLowQuality` to
> `counts.thresh = 100` and `gene.thresh = 20`. Interactions were reported at the `result.spx`
> level, that is, significant, quality-controlled and spatially variable by SPARK-X. For the
> two-condition comparison, we built one CytoSignal object per TMA core, merged them with
> `mergeCytoSignal` and tested for differential interaction usage between high- and low-grade
> cores with a negative-binomial mixed model as implemented in `runNEBULA` (`cpc_thresh = 0.001`),
> treating each core as an independent sample and applying Benjamini–Hochberg correction
> separately within diffusion- and contact-dependent interactions.

---

## stLearn — Python, v1.4.1, env `comp-stlearn`

Tutorial: `stlearn.readthedocs.io/en/latest/tutorials/cell_cell_interaction_xenium.html` — the
**Xenium-specific** CCI vignette (the generic `cell_cell_interaction.html` targets Visium).
Call contract: `scripts/comparators/stlearn/NOTES.md`. Deviations: `DEVIATIONS.md`.

### Core algorithm

Two stages, and they answer different questions.

**Stage 1 — `run()`: where is an LR pair co-expressed more than chance?** Cells are first
aggregated into a **regular spatial grid of spots**; the LR score of a spot is the co-expression
of ligand and receptor across that spot and its neighbours within a physical radius. The null is
built by sampling `n_pairs` **random gene pairs matched on expression level** to the real ligand
and receptor, scoring them identically, and asking where the real pair exceeds its own
background. This yields a per-spot p-value per LR pair, BH-corrected. The unit of inference is
the **spot**, not the cell and not a cell-type pair.

**Stage 2 — `run_cci()`: which cell types sit in those hotspots?** Restricted to significant
spots, it counts cell-type→cell-type edges across the neighbourhood graph, then permutes the
**cell-type labels/proportions** (not expression) to test whether a given cell-type pair is
over-represented in that LR's hotspots. With `spot_mixtures=True` a spot may count as several
cell types at once, using the per-spot proportions that `grid()` stored.

Note what this implies: stLearn never scores an individual cell, and cell-type identity enters
only in stage 2 — the LR statistics themselves are cell-type agnostic. This is a different
target of inference from CytoSignal, which tests every cell individually.

### Spatial model

| | |
|---|---|
| aggregation | regular grid, `n_row` × `n_col` bins over the bounding box; **empty spots dropped** |
| neighbourhood | `distance=250` physical units (cKDTree radius), **independent of grid resolution** |
| ours | 12,562 occupied spots, 51.3 × 51.3 µm, 8.0 cells/spot, **61 median neighbours**, 0 isolated |

The grid is an aggregation knob (the tutorial's own text calls `n_=125` a resolution/compute
trade-off); `distance` is the actual signalling range.

### LR database

Default is **connectomeDB2020_lit** (from NATMI), 2,293 pairs, via
`st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')`.

**stLearn cannot represent multi-subunit complexes.** Its LR format is the string
`"LIGAND_RECEPTOR"`, so `_` is the L/R separator and cannot also mean "subunit of". Converting
CellChatDB v2 therefore drops every complex row: **3,233 → 1,371 pairs (57.5% lost)**, of which
**527 are on the GBM 5K panel** and **526 survived stLearn's own expression filter**. Compare
CytoSignal, which is complex-aware and tested **1,088** panel interactions on the same data.

**This is not a cosmetic loss.** All 40 `WNT3_*` rows in CellChatDB v2 have complex receptors
(`FZD*_LRP5/6`), so **stLearn cannot test a single WNT3 interaction** — while WNT3–FZD*/LRP6 was
CytoSignal's *entire top-6* diffusion ranking on this dataset. Part of the two methods'
disagreement is therefore structural, not biological.

### Input

| requirement | detail |
|---|---|
| `X` | **raw counts, never log1p'd**. The tutorial is explicit: the permutation null picks background genes of similar expression, so log-shrinking genes together breaks it. Ours comes from `layers['counts']`. |
| normalisation | `st.pp.normalize_total` only (library size), *after* QC |
| coordinates | `obs['imagecol']`/`obs['imagerow']`, plus a Visium-style `uns['spatial']` dict |
| cell labels | any `obs` column; the tutorial substitutes Leiden because its demo lacks annotation |
| gene names | **must not contain `_`** |

### Workflow

| # | Call | Produces |
|---|---|---|
| 1–2 | `st.pp.filter_genes/filter_cells(min_counts=10)` | QC |
| 3 | `adata.raw = adata` | raw store |
| 4 | `st.pp.normalize_total(adata)` | library-size norm, **no log** |
| 5 | `st.tl.cci.grid(n_row, n_col, use_label=...)` | gridded AnnData + per-spot cell-type proportions in `uns` |
| 6 | `st.tl.cci.load_lrs([...])` | LR pair array |
| 7 | `st.tl.cci.run(grid, lrs, min_spots, distance, n_pairs, random_state)` | per-spot LR scores + p-values |
| 8 | `st.tl.cci.run_cci(grid, label, spot_mixtures, cell_prop_cutoff, sig_spots, n_perms)` | cell-type interaction counts + p-values |

### Data outputs

| Object | Shape | Meaning |
|---|---|---|
| `uns['lr_summary']` | LR × **6** | stage 1: `n_spots`, `n_spots_sig`, `n_spots_sig_pval` — plus the three `*_<label>` columns `run_cci()` appends (`analysis.py:770-772`), which are the stage-2 ranking |
| `obsm['lr_scores']` | spots × LR | raw LR co-expression score |
| `obsm['p_vals']`, `['p_adjs']`, `['-log10(p_adjs)']` | spots × LR | significance |
| `obsm['lr_sig_scores']` | spots × LR | score masked to significant spots |
| `obsm['spot_neighbours']`, `['spot_neigh_bcs']` | spots | neighbourhood graph |
| `uns['<label>']` | spots × types | per-spot cell-type proportions |
| `uns['lr_cci_<label>']`, `['lr_cci_raw_<label>']` | types × types | interaction counts pooled over LRs, significant and raw |
| `uns['per_lr_cci_<label>']` | dict of types × types | **one count matrix per LR** (526) |
| `uns['per_lr_cci_pvals_<label>']` | dict of types × types | **one p-value matrix per LR** — the significance of stage 2 |
| `uns['per_lr_cci_raw_<label>']` | dict of types × types | uncorrected per-LR counts |

`run_stlearn.py` persists most of these under `<out>/data/` as gzipped CSVs + `grid.h5ad`
(809 MB). **Two gaps it left**, closed by `export_stlearn_quant.py` (read-only replay from
`grid.h5ad`, nothing recomputed):

- it looped only over `["lr_cci_<label>", "per_lr_cci_<label>"]`, so the stage-2 **p-values**
  (`per_lr_cci_pvals_<label>`) and both raw-count sets never reached disk — the per-LR CSVs
  carried interaction counts with no significance attached;
- it wrote `lr_summary.csv` *before* `run_cci()`, so the on-disk copy had 3 of the object's 6
  columns. The rewrite was verified to leave all 3 pre-existing columns byte-identical.

### Image outputs

stLearn ships two CCI vignettes and **they do not have the same figure set**. The Xenium one is
the correct vignette for our data; the generic one targets Visium and never calls
`st.tl.cci.grid()`. Provenance was verified by grepping both vignettes and then adversarially
re-checked, because several figures we had been treating as "the standard workflow" are in fact
generic-vignette-only, and three are in neither.

| Call | Shows | Xenium vignette | generic Visium | Ours |
|---|---|---|---|---|
| `st.pl.cluster_plot` | cell types, grid vs single-cell | ✅ | ✅ | ✅ |
| `st.pl.feat_plot` | per-spot cell-type proportion (`vmax=1`) | ✅ CELL 42 | ❌ | ✅ all 9 types |
| `st.pl.gene_plot` | gene expression, grid vs single-cell | ✅ CELL 44 | ❌ | ✅ |
| `st.pl.lr_summary` | LR ranking bar chart | ✅ ×2 | ✅ | ✅ |
| `st.pl.lr_result_plot` | **per-LR spatial map** (`lr_scores`, `-log10(p_adjs)`, `lr_sig_scores`) | ✅ CELL 67, top-1 | ✅ | ✅ |
| `st.pl.cci_check` | **diagnostic: LR significance vs cell-type frequency** — should show ~no dependence if the permutation controlled for abundance | ✅ | ✅ | ✅ |
| `st.pl.ccinet_plot` | cell-type interaction network, overall + per LR, **shared node layout** via `return_pos`/`pos` | ✅ CELL 81 | ✅ | ✅ |
| `st.pl.lr_chord_plot` | chord diagram | ✅ CELL 83 | ✅ | ✅ |
| `st.pl.lr_diagnostics` | expression-vs-significance diagnostic | ❌ | ✅ | `plots/` only |
| `st.pl.lr_n_spots` | spots per LR | ❌ | ✅ | `plots/` only |
| `st.pl.cci_map`, `st.pl.lr_cci_map` | cell-type interaction heatmaps | ❌ | ✅ | `plots/` only |
| `st.pl.lr_plot` | per-spot ligand/receptor detail | ❌ | ✅ ×8 | `plots_full/` only |
| `st.pl.lr_go` | GO enrichment of top LR genes | ❌ | ✅ | ❌ **unavailable** — `run_lr_go(r_path=...)` needs R with `clusterProfiler`; no comp env has it |
| `st.pl.het_plot`, `st.pl.grid_plot`, `st.pl.deconvolution_plot` | diversity / gridding / composition pies | ❌ | ❌ | `plots_full/` only — **in neither vignette** |

**Three plot directories, deliberately not merged:**

| Directory | Script | Scope | Count |
|---|---|---|---|
| `plots_tutorial/` | `plot_stlearn_tutorial.py` | **The Xenium vignette call-for-call — cite this one as "the authors' default workflow".** | 23 + 6 in `requested/` |
| `plots/` | `run_stlearn.py` | Xenium set + 4 generic-vignette calls; top-6 where the tutorial takes top-1/top-2 | 29 |
| `plots_full/` | `plot_stlearn_full.py` | adds `lr_plot` (generic) + `het_plot`/`grid_plot`/`deconvolution_plot` (neither) | 26 |

Requested LRIs are isolated in a `requested/` subdirectory in both `plots_tutorial/` and `plots/`.
`plot_stlearn_full.py`'s `deconvolution_plot` calls all failed (`KeyError: 'deconvolution'` — it
writes `uns`, the plot reads `obsm`) and its `feat_plot` covered only 4 of 9 cell types with a
per-type colour scale; both are superseded by `plots_tutorial/` and are recorded in
`DEVIATIONS.md` rather than repaired, since neither call is in the Xenium vignette.

### Multi-sample / differential mode

**None.** stLearn has no native multi-sample or between-condition test — no equivalent of
CytoSignal's `mergeCytoSignal` + `runNEBULA`. Any grade or AIS-vs-LUAD contrast would have to be
hand-rolled, which the project rules forbid. Report this as a capability gap.

### Gotchas

- **`run_cci` is broken on pandas 3 / numpy 2** — `obs[label].values.astype(str)` yields a
  `StringDtype`/object array that numba 0.66 cannot type (`TypingError: non-precise type
  array(pyobject, 1d, C)` at `het.py:227`). Needs object-dtype coercion or `pandas<3`.
- **Cell types are matched to deconvolution columns by SUBSTRING, first hit wins**
  (`get_data_for_counting`). Any label that is a substring of another (`mGAM` ⊂ `non-mGAM`) can
  silently bind to the wrong column and corrupt every interaction reported for it. Verified
  correct here only because alphabetical order favours the shorter name.
- `grid()` hard-requires `uns['spatial']`; a plain h5ad raises `KeyError: 'spatial'`.
- **`n_pairs` has a hard floor of 100** — below it `run()` exits with a message and stores nothing.
- **Cost scales with unique LR *genes*, not `n_pairs × n_LRs`** — backgrounds are cached per gene
  (`gene_bg_genes`). 100× the pairs cost only 7.6× the time, so there is no reason to run below
  the authors' recommended 10,000.
- **`n_cpus` does nothing measurable** — 3.5 min at 8 cores vs 3.4 min at 1.
- `lr_cci_map` raises `UnboundLocalError` on an empty `lrs` list.
- **`st.tl.cci.adj_pvals` is a no-op at the tutorial's own arguments.** The Xenium vignette calls
  it (CELL 62) with `correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh'`, which
  `run()` has already applied internally — its own docstring says so (`analysis.py:382`).
  Verified empirically here: all five `obsm` matrices and all six `lr_summary` columns come back
  identical. Its one real effect is that it re-sorts `lr_summary` with an **unstable** `argsort`,
  permuting **179 / 526** LRs *within ties on `n_spots_sig`*. So calling it is harmless but also
  pointless, and anyone diffing two runs on LR order should expect tie noise, not a result change.
- **`deconvolution_plot` reads `obsm['deconvolution']`, not `uns`** (`pl/deconvolution_plot.py:86`;
  `adds/add_deconvolution.py:35` shows the intended write). Putting the frame in `uns` raises
  `KeyError: 'deconvolution'` — and if the caller wraps plots in a try/except, the failure is
  visible only in the log.
- **`grid()` stores per-spot proportions under `uns[<label>]`, which collides with nothing but is
  easy to mistake for the categorical**; `obs[<label>]` is the dominant-type label, `uns[<label>]`
  is the spots × types proportion matrix. `feat_plot` needs the latter copied into `obs` first.

### Deviations from the tutorial

Full table in `scripts/comparators/stlearn/DEVIATIONS.md`. Headlines: real `cell_type` instead of
Leiden; grid `321 × 146` = 51.3 µm square (preserving the tutorial's **spot area**, since its own
spots are 60.2 × 43.8 µm — rectangular — on a 1.37:1 section, whereas our TMA is 2.20:1);
`n_pairs` 1,000 → **10,000** and `n_perms` 100 → **1,000** (both are tutorial-declared "example,
recommend higher"); 21 `Intergenic_Region_*` control probes dropped; CellChatDB v2 complexes
dropped as forced by the LR format.

**The ~51 µm grid is *not* aligned to ALARMIST's 50 µm patch** — it falls out of matching the
tutorial's spot area, and the proximity is a coincidence of the arithmetic.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `results/comparators/stlearn/GBM/cellchatdb2/` | ✅ 45.4 min | 100,197 cells → 12,562 spots (7.98 cells/spot); **526 LR pairs** tested; all 526 have ≥1 significant spot, **482 have ≥20**. Top: C3–C3AR1 (1,244 sig spots), CNTN2–CNTN2, GJA1–GJA1, APP–SORL1, C3–CR2, FGF1–FGFR2 |
| GBM | `default` | — | ❌ not run | connectomeDB2020_lit tier deferred (this pass is cellchatdb2-only) |
| LUAD | both | — | ❌ not run | deferred |

Post-hoc passes over that run, both replay-only (no recomputation, `grid.h5ad` opened read-only):

| Pass | Script | Wall | Produced |
|---|---|---|---|
| Xenium-vignette figure set | `plot_stlearn_tutorial.py` | 0.8 min | `plots_tutorial/` — 23 PNGs + 6 in `requested/`, all 29 succeeded |
| stage-2 quantitative export | `export_stlearn_quant.py` | 0.04 min | `per_lr_cci_pvals_cell_type/` (526), `per_lr_cci_raw_cell_type/` (526), `lr_cci_raw_cell_type.csv`, 6-column `lr_summary.csv` |

`plot_stlearn_tutorial.py` rebuilds the cell-level `adata` (QC + `normalize_total`, ~30 s,
deterministic) because CELL 40 / 42 / 44 are grid-vs-single-cell side-by-side panels that cannot
be drawn from `grid.h5ad` alone. It asserts the rebuild against `run_manifest.json`
(100,197 cells, 5,096 genes, 12,562 spots) before plotting anything, so a silent preprocessing
drift would abort rather than produce a mismatched figure.

### Requested LRIs — and the first cross-method comparison

| LRI | stLearn rank | percentile | sig. spots / expressing | CytoSignal rank | percentile |
|---|---|---|---|---|---|
| GRN → SORT1 | **21 / 526** | top 4.0% | 267 / 8,763 | 66 / 895 | top 7.4% |
| ANXA1 → FPR1 | **99 / 526** | top 18.8% | 122 / 5,721 | 255 / 895 | top 28.5% |

Two methods with different units (spot vs cell), different nulls (expression-matched background
genes vs spatial permutation), and different neighbourhoods (250 µm grid radius vs 200 µm
epsilon ball) **agree on the ordering and rough standing** of both LRIs: GRN→SORT1 comfortably
top-decile, ANXA1→FPR1 mid-table. That concordance is meaningful precisely because nothing about
the two pipelines is shared.

Their *top* rankings, by contrast, are disjoint — and partly by construction, since stLearn
cannot see the complex-receptor WNT3 interactions that dominate CytoSignal's.

### Methods paragraph

> For stLearn (v1.4.1), we followed the authors' Xenium cell–cell interaction workflow. Cells
> were filtered (`min_counts=10` for both genes and cells), library-size normalised with
> `normalize_total` without log-transformation as the authors require, and aggregated onto a
> regular spatial grid with `st.tl.cci.grid`, using a bin size chosen to preserve the spot area
> of the authors' own tutorial (51.3 µm, giving 12,562 occupied spots at 8.0 cells per spot). We
> performed the ligand–receptor permutation test with `st.tl.cci.run` (`distance=250`,
> `min_spots=20`, `n_pairs=10000`) and obtained adjusted P values, then identified significant
> cell-type interactions with `st.tl.cci.run_cci` (`spot_mixtures=True`, `cell_prop_cutoff=0.1`,
> `sig_spots=True`, `n_perms=1000`), in both cases using the authors' recommended rather than
> example permutation counts. Because stLearn encodes a ligand–receptor pair as a single
> underscore-delimited gene string, multi-subunit complexes cannot be represented, and the
> CellChatDB v2 resource was restricted to its 1,371 single-gene pairs (526 detectable on the
> Xenium panel).

---

## SpatialDM — Python, v0.3.1, env `comp-spatialdm`

Tutorials: `/Users/jiayifan/tansey_lab/SpatialDM/tutorial/{melanoma,differential_test_intestine}.ipynb`.
Call contract: `scripts/comparators/spatialdm/NOTES.md`. Deviations: `DEVIATIONS.md`.

### Core algorithm

SpatialDM tests **spatial co-expression** of a ligand and receptor using a **bivariate Moran's
statistic**, in two nested stages.

**Global (`spatialdm_global` → `sig_pairs`)** — for each LR pair, a bivariate Moran's R asks
whether ligand expression at a location is spatially associated with receptor expression at its
*neighbours*, given a spatial weight matrix `W`. The `z-score` method uses an **analytical
null**: the variance of R is derived in closed form from the structure of `W` itself
(`globle_st_compute`), not from permutation. BH-corrected, `fdr < 0.1` selects pairs.

**Local (`spatialdm_local` → `sig_spots`)** — for the globally selected pairs only, a local
Moran statistic is computed *per cell*, giving the subset of locations where that interaction
actually occurs. So the unit is the **pair** globally and the **cell** locally.

Two consequences matter for the comparison. First, this is a **correlation/autocorrelation**
framework, not a co-expression-product permutation like CytoSignal or stLearn — it asks whether
L and R are spatially *arranged* together, not whether their product exceeds a background.
Second, **because the null is analytical in `W`, statistical power depends directly on
neighbourhood structure** — which is what makes the density confound below so consequential.

### ⚠️ Overlaps with LIANA+ — these two are not independent votes

**LIANA+'s `morans` local metric is a reimplementation of this method**, and its source says so
(`liana/method/sp/_bivariate/_local_functions.py:458` cites the SpatialDM paper; the package even
carries a helper named `_spatialdm_weight_norm`). The formulae are identical —
`local_x = x*(W@y); local_y = y*(W@x); local_r = local_x + local_y` versus SpatialDM's
`local_I + local_I_R` — and for `local_name='morans'` LIANA+ applies SpatialDM's exact
preprocessing (`_norm_max` on both matrices, `W·N/ΣW`) and offers the same analytic z-score null
at `n_perms=0`.

Everything around the statistic differs, and on our data that dominates:

| | SpatialDM `cellchatdb2` | LIANA+ `cellchatdb2_morans` |
|---|---|---|
| support | 13 cores, separate fits | whole slide, one run |
| radius | **134.6 µm** | **28.2 µm** |
| null | analytic z-score | permutation, `n_perms=100` |
| filter | `min_cell=3` | `nz_prop=0.02` |
| pairs tested | **1,133–1,661 per core** | **131** |

So a pair agreeing between the two is a **sanity check, not corroboration**. SpatialDM's distinct
contributions to this benchmark are the analytic null (fast, and demonstrably liberal — see below)
and `diff_utils`, a native per-sample differential mode LIANA+ has no equivalent of. Note also
that `morans` is *not* LIANA+'s default local metric (`local_name='cosine'` is), and that Moran's
R is not NMF-admissible (33.4% negative entries), so LIANA+'s communication-program branch cannot
run on it.

### Spatial model

| | |
|---|---|
| kernel | RBF, `w(d) = exp(-d² / 2l²)`, then weights `< cutoff` zeroed |
| ours | `l = 75`, `cutoff = 0.2` → **effective radius 134.6 µm** |
| graph | `n_neighbors` per core (94–709), sized so `cutoff` truncates and the kNN cap never binds |
| adjacent graph | separate 6-NN graph → `obsp['nearest_neighbors']`, used for contact/ECM pairs |
| `single_cell=True` | zeroes the diagonal: a cell cannot signal to itself |

`W` (RBF) is used for `Secreted Signaling`; `nearest_neighbors` for `ECM-Receptor` and
`Cell-Cell Contact` — the split is by database annotation, and `geneInter` **must be sorted by
annotation** because `st` and `n_short_lri` index blocks positionally.

### LR database

Default is **CellChatDB v1** — **1,939 human interactions + 157 complexes**, an
`interaction_input` + `complex_input` pair. For the `cellchatdb2` tier we inject CellChatDB v2
directly, because `extract_lr` has no custom-DB hook. **Lossless for complexes** — SpatialDM
handles multi-subunit interactions natively via `Ligand0..N`/`Receptor0..N` columns, so unlike
stLearn (which must drop 57.5% of v2) nothing is lost for being a complex. Valid pairs per core:
**568–825 on v1**, **1,133–1,661 on v2**, versus stLearn's 526 whole-slide.

**`extract_lr`'s `datahost` argument is named backwards from what it does.** `datahost='package'`
reads the CSVs shipped inside the wheel (`spatialdm/datasets/LR_data/`); `datahost='builtin'` —
the function's own **default**, and therefore what the tutorial silently gets — **downloads from
figshare**, once per call. We pass `--datahost package`: identical CellChatDB v1, but offline,
pinned to the installed version, and not re-fetched 13 times. (The `package` branch needs
`pkg_resources`, which setuptools ≥ 81 no longer ships — see Gotchas.)

**v2's `Non-protein Signaling` category is remapped to `Secreted Signaling`** — see DEVIATIONS.md;
SpatialDM's v1-era code does not enumerate it and crashes otherwise. This is not cosmetic: core
13's entire top-6 (`SLC17A7_GLS2_GRIA2`, `SLC17A7_GLS_GRIN1_GRIN2A`, …) is glutamatergic and
would have been lost had the category been dropped instead.

### Input

Log-transformed expression in `.X` **and raw counts in `.raw`** — both are genuinely used
(`extract_lr`/global read `.X`; **local** spot selection reads `.raw`). Coordinates in
`obsm['spatial']`. Our GBM h5ad matches natively; `.raw` is wired from `layers['counts']`.

### Workflow

| # | Call | Produces |
|---|---|---|
| 1 | `weight_matrix(l, cutoff, n_neighbors, n_nearest_neighbors, single_cell)` | `obsp['weight']`, `obsp['nearest_neighbors']` |
| 2 | `extract_lr(species, min_cell)` *(or v2 injection)* | `uns['ligand'/'receptor'/'geneInter'/'num_pairs']` |
| 3 | `spatialdm_global(n_perm, method='z-score')` | `uns['global_I']`, `uns['global_stat']` |
| 4 | `sig_pairs(method='z-score', fdr=True, threshold=0.1)` | `uns['global_res']` with `.selected` |
| 5 | `spatialdm_local(n_perm, method='z-score')` | `uns['local_z']`, `uns['local_z_p']` (pairs × cells) |
| 6 | `sig_spots(method='z-score', fdr=False, threshold=0.1)` | `uns['selected_spots']`, `uns['local_stat']['n_spots']` |

### Data outputs

| File | Shape | Meaning |
|---|---|---|
| `global_res.csv` | pairs × 9 | Moran's R, p, **fdr**, `selected` — the ranking |
| `local_z.npz`, `local_z_p.npz` | selected-pairs × cells | **per-cell** local statistic and p-value (float32) |
| `selected_spots.csv.gz` | selected-pairs × cells | binary: where the interaction is called |
| `local_n_spots.csv` | pairs | number of selected cells per pair |
| `lr_{ligand,receptor}_subunits.csv` | pairs × subunits | resolved complex membership |
| `cell_meta.csv` | cells | coordinates + cell type |
| `per_split_summary.csv` | 13 rows | per-core diagnostics (see below) |

**Note `local_*` cover only globally-selected pairs**, by design — local selection runs on the
global survivors.

### Image outputs

| Call | Shows | Ours |
|---|---|---|
| `pl.global_plot` | global selection volcano | ✅ per core |
| `pl.plot_pairs` | **per-LR spatial maps** (L, R, local significance) | ✅ top 3 + requested, as `selected_pair_*.png`. The bundled `plot_pairs.pdf` is **truncated in cores 3/9/12** — `plot_selected_pair` indexes `spots`/`n_spots` by *globally selected* pairs only, so it raises on the first non-selected requested LR and the `PdfPages` context closes early. The per-pair PNGs are unaffected. |
| `pl.chord_celltype` | cell-type chord for a pair | ✅ top 2 per core (`obsm['celltypes']` built as one-hot from the single-cell labels — no deconvolution needed) |
| `pl.chord_LR` | self-self chord per cell type (senders/receivers are zipped, not crossed) | ⚠️ cellchatdb2 **12/13**, default **7/13**. For `default` the log gives both causes: 4 cores (2, 6, 10, 12) raise `zero-size array to reduction` — no edge survives `min_quantile=0.5` once `Links.value>0`, which the sparser v1 DB hits more often — and cores 1, 5 fail reproducibly inside the headless-browser render (`JavascriptException: root_view is undefined`); retried, identical. For `cellchatdb2` the single missing core (2) has **no log**, so its cause is unverified — core 2 is the smallest (819 cells, 686 pairs), so the `zero-size` path is the likely one, but that is inference. |
| `pl.dot_path` + `compute_pathway` | pathway enrichment | ⚠️ run per core, but **degenerate**: `dic={core: all selected}` instead of the tutorial's `dic={Pattern_i: …}`, so there is a single group and the background set is empty |
| `pl.ligand_ct` / `pl.receptor_ct` | per-cell-type ligand/receptor contribution | ✅ as `celltype_contrib_*.{csv,png}` — **these two return DataFrames, they are not plotters** |
| SparseAEH `plot_clusters` | spatial clustering of local patterns | ❌ optional extra package, not installed; this is what would make `dot_path` meaningful |
| weight-range scatter (melanoma cell 13) | visualises the RBF kernel extent | ❌ diagnostic, not produced |

**Two plotting passes were run per tier and both are on disk — see `plots/` vs `plots_full/`
below.** `plots_full/` holds **181 PNGs + 13 PDFs** (cellchatdb2) and **163 PNGs + 13 PDFs**
(default). A requested LR is plotted only where it is **globally selected** — SpatialDM computes
no local statistic otherwise, so there is nothing to draw; per-core status is in
`per_split_summary.csv`.

⚠️ **The `default` tier's `plots/` still contains 143 blank 7,544-byte PNGs** (130 in `plots/`,
13 in `plots/requested/`) from the first pass — the same bug described below. They were deleted
from the `cellchatdb2` tree at some point but have not been cleaned here. The only real files in
`default/*/plots/` are 13 `global_plot.png` and 13 `dot_path.png`.

#### `plots/` vs `plots_full/`

`plots/` was written by `run_spatialdm.py` during the compute run using
`fn(); plt.gcf().savefig(...)`, which is the **wrong mechanism for three classes of SpatialDM
function** and produced blank 7,544-byte PNGs: `plot_pairs` ends each iteration with
`plt.show(); plt.close()` so `gcf()` is a fresh empty figure; `ligand_ct`/`receptor_ct` are not
plotters at all; `chord_celltype`/`chord_LR` render through holoviews/bokeh and take their own
`save=`. Only `global_plot` survives there. `plots_full/` is the corrected pass
(`plot_spatialdm_full.py`), replotted **from the persisted `data/spatialdm.h5ad` with no
recomputation**, with a saver that refuses to write a figure that has no drawn content.
`plots/` is retained only for provenance; **read `plots_full/`**.

### Multi-sample / differential mode

**Yes — `spatialdm.diff_utils`**, demonstrated in `differential_test_intestine.ipynb`. **Run**
(2026-08-06, `run_diff_spatialdm.py` → `GBM/cellchatdb2/differential_grade/`): 7 high-grade vs
6 low-grade cores, `concat_obj` → `differential_test` → `group_differential_pairs`, then the
tutorial's clustermap / dendrogram / volcano / dot_path. Nothing is recomputed — it consumes the
persisted per-core objects. Two patches were required (`concat_db` hard-codes CellChatDB v1;
`dot_path`'s installed signature differs from the tutorial's) — see DEVIATIONS.md. Result below.

### Grade differential — the method finds nothing, and says why

| | |
|---|---|
| union LR pairs across 13 cores | **1,662** (581 testable in all 13; 12.2% of the pair × core grid is zero-filled) |
| pairs at differential FDR < 0.1 | **0** (minimum FDR **0.114**) |
| `high_specific` / `low_specific` | **0 / 0** |
| pairs selected in all 7 high cores and no low core | **0** |
| raw p < 0.05 | 162 of 1,662 (83 expected by chance) — signal exists, none survives BH at n=13 |

**The reason is the density confound, and three independent views agree on it.** (i) Per-core
mean \|z\| tracks neighbourhood size, **r = +0.752, p = 0.003** (core 13, `n_neighbors`=693 →
mean \|z\| 5.67; core 6, `n_neighbors`=94 → 0.60). The differential test regresses exactly that
quantity, so the high-grade arm carries most of the between-sample variance and the
likelihood-ratio test has little power left. (ii) `differential_dendrogram` clusters the 13 cores
almost perfectly by grade **except for cores 9 and 11** — the two cores where grade and density
disagree (9 is low-grade but dense, 11 is high-grade but sparse). (iii) Re-running the authors'
own test with a **median-density split**, which differs from the grade design at exactly those
two cores, yields **5 condition-specific pairs** (all five selected in all 7 dense cores and no
sparse core: `CSF1_CSF1R` FDR 0.0075, `BMP2_BMPR1B_ACVR2A`, `TNC_ITGAV_ITGB3`, `TNC_ITGAV_ITGB6`,
`IL16_CD4`) against **0** for grade.

**Conclusion for the benchmark: on this TMA, SpatialDM's native differential test separates
cellularity, not grade — and the two are not separable with 13 cores.** High cellularity is
itself a WHO grading criterion, so this is expected biology hitting a method whose null derives
from `W`; it is not a pipeline artifact.

**The design is also largely paired, and the test does not use it.** The 13 cores come from
**7 patients** (`scripts/research/gbm_supplemental.py`; the h5ad has no patient column, so that
hardcoded dict is the only provenance — and it corrects a "8 patients" figure that had propagated
through this file). **4 of the 7 contribute both a high- and a low-grade core, covering 10 of the
13 cores**, so the grade contrast is mostly within-patient. `differential_test` fits
`y ~ 1 + conditions` by OLS with no patient term, so it discards the pairing *and* treats patient
14007's four cores as four independent observations — pseudoreplication and lost power at once. A
paired/mixed model on the persisted `zscore_df` is the natural follow-up; it is not SpatialDM's
own workflow, so it was not run. See DEVIATIONS.md for the per-patient table.

Requested LRIs in the grade contrast: **GRN→SORT1** diff +4.39, FDR 0.62 (mean z 6.07 high vs
1.68 low, selected in 5/7 high and 2/6 low); **ANXA1→FPR1** diff +0.38, FDR 0.97 (2.95 vs 2.56,
4/7 and 2/6). Neither is condition-specific by this method.

### Gotchas

- **`eff_dist` is a squared distance despite its name** — `l = sqrt(-eff_dist/(2·ln cutoff))`
  against a kernel `exp(-d²/2l²)`. `eff_dist=135, cutoff=0.2` gives `l=6.48` and a weight of
  4e-95 at 135 µm instead of 0.2. Always set `l` directly.
- **`n_neighbors` and `n_nearest_neighbors` are independent** — the `n_neighbor_layers*31`
  derivation fires only when `n_neighbors is None`; `n_nearest_neighbors` sizes a *separate*
  small graph.
- **Once the kNN cap does not bind, `W` is numerically identical for any larger `n_neighbors`.**
  Only stored `nnz` differs, because `rbfweight` zeroes sub-cutoff entries without calling
  `eliminate_zeros()` — stored `nnz = n_cells × n_neighbors` regardless of `cutoff`. Compacting
  is safe (nothing downstream uses `.nonzero()`/`.nnz`) and here saved **17.4M** entries.
- **CellChatDB v2's `Non-protein Signaling` crashes `spatialdm_global`** (shape mismatch in
  `st`), and is inconsistently classified even beforehand.
- `sig_spots` **defaults to `fdr=True`** but the tutorial passes `False`.
- `local_z`/`local_z_p`/`local_perm_p` live at `uns` **top level**, not under `uns['local_stat']`.
- **`local_z_p` is NOT `norm.sf(local_z)`.** `spot_selection_matrix` ends with
  `np.where(pos.T == False, 1, local_z_p)` — the p-value is forced to 1 wherever *neither*
  standardised side is above its mean. On GBM core 1 that masks 86.8% of the matrix. The identity
  holds only on the unmasked entries (verified, max abs diff 4.4e-08).
- **`local_z.npz` carries placeholder row labels** (`'0','1','2',…`): `uns['local_z']` is a bare
  ndarray, so the runner's `getattr(v,'index',arange)` fell through. Row order **is** identical to
  `local_z_p.npz`, which does carry the pair names — verified via the identity above. Use
  `local_z_p`'s `pairs` array to label `local_z`.
- **Core 13 has one cell that yields NaN local statistics** for 753 of its 863 selected pairs —
  column 8370, cell id 68901, an `AC-like` cell. A single cell, always the same one; drop it
  before any correlation over `local_z`.
- **`extract_lr(datahost=…)` is named backwards**: `'builtin'` (the default) downloads from
  figshare; `'package'` reads the wheel's own CSVs. See LR database above.
- **`datahost='package'` needs `pkg_resources`, which setuptools ≥ 81 removed.** The env had
  setuptools 83 and the call died with `ModuleNotFoundError`. Pinned to `setuptools<81`
  (80.10.2) and frozen in `env.lock.yml`.
- **Chord plots need `geckodriver`/`firefox` on `PATH`, and they are not.** bokeh exports PNG by
  driving a headless browser. Both binaries *are* installed in `comp-spatialdm/bin/`, but this
  repo invokes the interpreter by absolute path (conda activate is broken here), which leaves
  that `bin/` off `PATH` — so every chord diagram fails with "Neither firefox and geckodriver …
  available on system PATH". Both runners now prepend `dirname(sys.executable)` to `PATH`.
  This is almost certainly why the `cellchatdb2` chords had to be produced by a separate,
  unpreserved pass.
- **`chord_LR(senders=…, receivers=…)` indexes `adata.obs[sender]`** — it wants one obs *column*
  per cell type (the intestine tutorial's deconvolution proportions), not category names, and the
  two lists are **zipped, not crossed** (so it draws the self-self diagonal). The replot script
  copies the one-hot `obsm['celltypes']` into `obs` first.
- `chord_LR` still hardcodes `title='Undifferentiated_Colonocytes'` from the authors' tutorial.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `results/comparators/spatialdm/GBM/cellchatdb2/<core>/` | ✅ 9.3 min, 13 cores | 100,197 cells; **1,133–1,661 valid pairs** per core; **54–863 significant** (FDR<0.1); kNN cap binds **0.000% in every core**; stored nnz 47.3M → **29.9M** after `eliminate_zeros` |
| GBM | `cellchatdb2` | `.../GBM/cellchatdb2/differential_grade/` | ✅ 0.4 min, peak RSS 3.6 GB | grade differential over the same 13 cores; **0 pairs at FDR < 0.1**; see the section above |
| GBM | `default` | `results/comparators/spatialdm/GBM/default/<core>/` | ✅ 12.9 min, 13 cores, peak RSS 6.0 GB | same 100,197 cells and same `W`; **568–825 valid pairs** per core (1,939-interaction v1 DB); **33–386 significant**; `--datahost package` (offline) |
| GBM | `default` | `.../GBM/default/differential_grade/` | ✅ 0.4 min | grade differential; **0 pairs at FDR < 0.1** — reproduces the v2 tier exactly |
| — | comparison | `.../GBM/tier_comparison/` | ✅ | `compare_tiers.py`: v1 ↔ v2 concordance, see below |
| LUAD | both | — | ❌ not run | deferred |

Per-core `n_neighbors`, stored nnz, cap-bind fraction, valid/significant pairs and runtime are
in `per_split_summary.csv`.

### `default` vs `cellchatdb2` — the LR database is not driving anything

Both tiers use the **same cells, same `W`, same `n_neighbors`** (verified identical per core), so
they differ only in DB membership. Pairs are matched across tiers by resolved subunit sets, not by
name (v1 names complexes `TGFB1_TGFBR1_TGFBR2`; our v2 names are `ligand_receptor`).

| | v1 (`default`) | v2 (`cellchatdb2`) |
|---|---|---|
| DB size | 1,939 interactions | 3,233 rows / 3,218 unique pairs |
| pooled valid pairs over 13 cores | 9,132 rows → **8,528 distinct** (6.6% redundant) | 18,971 rows → **16,057 distinct** (15.4% redundant) |
| pooled significant | **2,025** | **5,201** |
| of v2's pairs, share also in v1 | — | **52.6%** (median per core) |

Three findings, all pointing the same way:

1. **For any pair present in both DBs, SpatialDM computes the identical statistic — pooled
   Pearson r of the Moran z across all shared pairs = 1.000.** The database cannot change what
   the method says about a pair; it can only change *which* pairs are tested, and through that
   the BH denominator.
2. **Agreement on shared pairs is near-total: median Jaccard of the selected sets = 0.954**
   (0.78–1.00). The handful of disagreements per core are pairs sitting on the FDR boundary,
   moved by v2 testing roughly twice as many hypotheses.
3. **The two requested LRIs behave identically.** GRN→SORT1 testable in 13/13 and selected in
   exactly cores {1, 2, 5, 6, 8, 10, 13} under **both** DBs; ANXA1→FPR1 testable in the same
   11 cores and selected in exactly {1, 3, 8, 9, 12, 13} under both. Not a v2 artefact.

**So the `Non-protein Signaling` remap inflates the counts but does not change the conclusions.**
It supplies **48.6% of v2's distinct significant pairs** (2,065/4,248; 55.2% of raw selected rows
before collapsing duplicate resolved pairs — Non-protein rows are the more redundant ones, e.g.
four `SLC6A11/SLC32A1/SLC6A13/SLC6A8_GAD2_GABRA3_GABRB3_GABRQ` rows that reduce to one measurable
pair on a 5,119-gene panel). Removing it entirely — which is what the v1 tier effectively does —
leaves the requested-LRI verdict, the grade differential and the density confound unchanged.

**The density confound is not a v2 artefact either; it is stronger on v1:**

| | v1 | v2 |
|---|---|---|
| corr(`n_neighbors`, n_significant) | **r = +0.731, p = 0.005** | r = +0.709, p = 0.007 |
| corr(`n_neighbors`, **fraction** significant) | **r = +0.631, p = 0.021** | r = +0.442, p = 0.131 |
| high vs low grade, significant fraction | 0.275 vs 0.149, p = 0.138 | 0.308 vs 0.252, p = 0.628 |

And the grade differential reproduces on v1 to the digit: **0 pairs at FDR < 0.1** (min FDR 0.121
vs 0.114 on v2), `high_specific` = `low_specific` = 0, and the density-split sensitivity returns
**the identical five pairs** — `CSF1_CSF1R`, `BMP2_BMPR1B_ACVR2A`, `TNC_ITGAV_ITGB3`,
`TNC_ITGAV_ITGB6`, `IL16_CD4` — with the same `diff` for both requested LRIs (+4.39 / +0.38).

### The power confound, measured

Because the z-score null derives from `W`, power scales with neighbourhood size — and GBM core
density correlates with grade (r = 0.659, p = 0.014; 4× difference). The consequence is directly
visible in the results:

| relationship | statistic |
|---|---|
| corr(`n_neighbors`, **n_significant**) | **r = +0.709, p = 0.007** |
| corr(`n_neighbors`, **fraction** significant) | r = +0.442, p = 0.131 |
| high vs low grade, significant **fraction** | 0.308 vs 0.252, **p = 0.628** |

**Raw counts of significant pairs track neighbourhood size, and the apparent grade difference
disappears once normalised to the fraction of testable pairs.** Any cross-grade claim from this
method must therefore be made on fractions, never counts.

### Requested LRIs

Identical under both LR databases — same testable cores, same selected cores, not just the same
counts (`tier_comparison/requested_lr_by_tier.csv`):

| LRI | testable | selected (FDR<0.1) | which cores | high-grade | low-grade |
|---|---|---|---|---|---|
| GRN → SORT1 | **13/13 cores** | 7 | 1, 2, 5, 6, 8, 10, 13 | 5/7 | 2/6 |
| ANXA1 → FPR1 | 11/13 cores | 6 | 1, 3, 8, 9, 12, 13 | 4/6 → 4/7 | 2/6 |

Both lean high-grade, but with 13 cores from 7 patients this is not a significant association
(Fisher p = 0.29 and 0.59) and should not be read as one. GRN→SORT1 being testable in **every**
core is itself notable. Neither is condition-specific in the differential test (§ Grade
differential).

A requested LRI is **plotted only where it is globally selected** — SpatialDM computes no local
statistic for a non-selected pair, so `plot_selected_pair`, `ligand_ct` and `receptor_ct` have
nothing to index. The absence of `selected_pair_GRN_SORT1.png` in 6 cores is that, not a failure;
the per-core verdict is in `per_split_summary.csv` either way.

### Methods paragraph

> For SpatialDM (v0.3.1), we followed the authors' workflow: we computed the neighbour-weight
> graph, identified globally significant ligand–receptor pairs and then selected locally
> significant spots, using `weight_matrix`, `extract_lr`, `spatialdm_global`, `sig_pairs`,
> `spatialdm_local` and `sig_spots`. The radial-basis kernel was set to `l = 75` with
> `cutoff = 0.2`, reproducing the effective 135 µm signalling radius of the authors' intestine
> tutorial — the one example in which they state their reasoning for this parameter — and
> `single_cell=True` was used because Xenium resolves individual cells. The size of the
> neighbour graph was set per tissue core so that the distance cutoff, rather than the
> k-nearest-neighbour ceiling, determined the neighbourhood, matching the regime of both author
> tutorials. Each of the 13 tissue microarray cores was analysed independently, following the
> authors' multi-sample design, because the global bivariate Moran statistic is returned once per
> ligand–receptor pair per object and the differential test operates on separately fitted samples.
> The analysis was run twice, once with SpatialDM's own CellChatDB v1 resource and once with
> CellChatDB v2, whose multi-subunit complexes SpatialDM represents natively; for ligand–receptor
> pairs present in both resources the two runs gave identical Moran statistics (Pearson r = 1.000)
> and near-identical significance calls (median Jaccard 0.954). To compare high- and low-grade
> cores we used SpatialDM's native
> differential mode (`spatialdm.diff_utils.concat_obj`, `differential_test`,
> `group_differential_pairs`), which fits a likelihood-ratio test of each pair's per-core Moran
> z-score on the grade label; no ligand–receptor pair reached a differential false-discovery rate
> below 0.1.

---

## COMMOT — Python, v0.0.3, env `comp-commot`

Tutorials: `/Users/jiayifan/tansey_lab/COMMOT/docs/notebooks/{Basic_usage,visium-mouse_brain}.ipynb`.
Call contract: `scripts/comparators/commot/NOTES.md`. Deviations: `DEVIATIONS.md`.

### Core algorithm

COMMOT casts cell–cell communication as **collective optimal transport**. Ligand "mass" at
sender cells is transported to receptor "mass" at receiver cells so as to minimise a spatial
cost, subject to `dis_thr` forbidding any coupling beyond that distance. The defining property
— and the authors' main argument over pairwise scoring — is that transport is **competitive**:
multiple ligands and receptors compete for the same finite mass, so a receptor already saturated
by a nearby strong sender cannot also absorb signal from a distant one.

The output is a **cell × cell transport plan per LR pair**, summarised into per-cell amounts
sent and received. Crucially, **there is no significance test at the LR-pair level** — COMMOT
returns magnitudes, not p-values. Permutation p-values exist only at the **cell-type-pair**
level, per pathway, via `cluster_communication`. Any ranking of LR pairs here is therefore by
*total received signal*, a magnitude, and is **not comparable like-for-like** to CytoSignal's
significant-cell counts, stLearn's significant-spot counts, or SpatialDM's FDR.

### Spatial model

A single hard Euclidean cutoff, `dis_thr`, in the units of `obsm['spatial']`. No kernel, no
decay — coupling is either permitted or forbidden. **Ours is 365 µm**, derived by measuring the
tutorial dataset rather than copying its number (see below). Applied identically to every pair,
which is why only diffusible signalling types belong in the database subset.

### LR database

Built-in options are **CellChatDB (v1)** and **CellPhoneDB v4.0**, via
`ct.pp.ligand_receptor_database(database, species, signaling_type)`. `df_ligrec` is a plain
3-column frame — ligand, receptor, pathway — with heteromeric subunits joined by `_`, which is
**exactly the encoding our CellChatDB v2 CSV already uses**, so the v2 handover is direct and
lossless. Complexes are supported natively (`heteromeric=True`, `heteromeric_rule='min'`).

We use **Secreted Signaling + Non-protein Signaling** — 2,259 pairs, 859 of them on the Xenium
panel, **671 surviving `min_cell=100` evaluated once over all 100,197 cells** and reused for
every core (`--filter-scope global`).
The tutorial restricts to Secreted only; we keep that principle — a single 365 µm transport
radius is a diffusion model — but add v2's Non-protein category, which did not exist in the v1
the tutorial used and is equally diffusible. That inclusion is vindicated by the results:
**Glutamate is a top-3 pathway in cores 9, 13 and 14**, matching SpatialDM's independent finding
that glutamatergic pairs dominate core 13.

### Input

`adata.raw = adata`; `sc.pp.normalize_total`; `sc.pp.log1p` — applied to **raw counts**. COMMOT
requires "non-negative values that reasonably reflect the abundancy of signaling molecules".

### Workflow

| # | Call | Produces |
|---|---|---|
| 1 | `ct.pp.ligand_receptor_database(...)` *(or our v2 frame)* | `df_ligrec` |
| 2 | `ct.pp.filter_lr_database(df, adata, filter_criteria, min_cell)` | expression-filtered pairs |
| 3 | `ct.tl.spatial_communication(adata, database_name, df_ligrec, dis_thr, heteromeric=True, pathway_sum=True)` | transport plans + per-cell sums |
| 4 | `ct.tl.cluster_communication(adata, database_name, pathway_name=<top-5 by signal>, clustering='cell_type', n_permutations=100, random_seed=0)` | cell-type × cell-type matrix **+ permutation p-values** |
| 4b | same, with `lr_pair=('GRN','SORT1')` / `('ANXA1','FPR1')` | p-values for the **requested pairs** — the only significance test COMMOT gives a named pair |
| 5 | `ct.tl.communication_direction(..., k=5)` → `ct.pl.plot_cell_communication(background='summary')` | signalling vector fields |
| 6 | `ct.tl.cluster_position` → `ct.pl.plot_cluster_communication_network` | cell-type communication network |
| 7 | *(not run)* `communication_impact`, `communication_deg_detection` → `communication_deg_clustering` | signalling-DE genes |

### Data outputs

| Key / file | Shape | Meaning |
|---|---|---|
| `obsp['commot-<db>-<lig>-<rec>']` | cells × cells | transport plan per LR pair (sparse) |
| `obsp['commot-<db>-<pathway>']`, `-total-total` | cells × cells | pathway and total aggregates |
| `sum_sender.csv.gz`, `sum_receiver.csv.gz` | cells × (pairs+pathways) | amount sent / received per cell |
| `lr_total_received.csv` | pairs | **LR pairs only** — the ranking |
| `pathway_total_received.csv` | pathways | pathway aggregates, kept separate |
| `cluster_comm_<key>.csv`, `cluster_pval_<key>.csv` | types × types | cell-type communication + p-values. `<key>` = each of the top-5 pathways **and** each requested LR pair → 7 pairs of files per core |
| `lr_pairs_used.csv` | pairs × 3 | the post-filter database actually used (identical in all 13 cores under `--filter-scope global`) |
| `lr_pairs_global.csv` | 671 × 3 | *(run root)* the one globally-filtered pair set handed to every core |
| **`<core>/adata_commot.h5ad`** | — | **the full AnnData including `obsp`** — every transport plan, the sum-sender/receiver frames, the vector fields, `uns['commot-<db>-info']`, the 7 `cluster_communication` results and `cluster_pos`. `.raw` dropped (reconstructible from `--h5ad`). 0.01–0.60 GB/core, **2.13 GB** total |

**Why the h5ad matters:** COMMOT has no way to recover a transport plan except by re-solving the
OT. Before 2026-08-10 nothing persisted `obsp` (`--save-obsp` was declared at argparse and never
referenced — a dead flag), so `communication_direction`, `communication_impact`,
`deg_detection`, `spatial_autocorrelation` and the `group_*` family each implied a fresh ~5-hour
run. Verified by reloading `2/adata_commot.h5ad` and running `ct.tl.communication_direction`
on it with no OT re-run.

**Trap:** `pathway_sum=True` writes pathway columns *into the same* `sum_receiver` frame as the
per-pair columns. Ranking pairs without separating them inflates the denominator and shifts every
rank. In the 2026-08-01 artifacts this was **not** separated: `lr_total_received.csv` held
pairs + pathways (261 rows against 217 real pairs in core 1), `pathway_total_received.csv` was
never written, and the shipped ranks put GRN→SORT1 at median 8 instead of its true **median 2**.
Fixed and re-run; both files are now correct and disjoint.

### Image outputs

**41 figures per core × 13 cores = 533** (507 PNG + 26 PDF), after the 2026-08-10 re-run. The
2026-08-01 run produced **54 PNGs and not one native `ct.pl` figure**.

| Plot | What it shows | File | n |
|---|---|---|---|
| `sender_receiver_map` *(ours)* | per-cell sent/received, top **LR pairs** by total received | `signal_<lig>-<rec>.png` | 3/core |
| `sender_receiver_map` *(ours)* | same, **pathway** aggregates — labelled so they cannot be read as pairs | `pathway_<pw>.png` | 5/core |
| `sender_receiver_map` *(ours)* | the standing requested LRIs | `requested_rank<N>_<LR>.png` | 2/core |
| `ct.tl.communication_direction` → **`ct.pl.plot_cell_communication`** | signalling vector field, **tutorial arguments** (`plot_method='grid'`, `normalize_v=True`, `normalize_v_quantile=0.995`, `grid_density=0.4`, `ndsize=8`). Two backgrounds × sender/receiver × (5 pathways + 2 requested) | `native_vf_{sender,receiver}_{pathway_<pw>,requested_<LR>}_{sig,ct}.png` — `sig` = `background='summary'`/`Reds`, `ct` = `background='cluster'`/`Alphabet` | **28/core** |
| **`ct.pl.plot_cluster_communication_network`** | cell-type communication network + its colour legend | `native_cluster_network.pdf`, `native_cluster_network_cluster_legend.pdf` | 2/core |
| `ours_dotplot` *(ours — substitute)* | cell-type × cell-type dots, size = −log10(p), colour = strength | `ours_dotplot_top.png` | 1/core |
| `ct.pl.plot_cluster_communication_dotplot` | — | ❌ **unrunnable, version gap** (see Gotchas) | 0 |
| `ct.pl.plot_communication_impact` | downstream-impact heatmap | ⚠️ partially — see *Downstream impact* below | — |
| `ct.pl.plot_communication_dependent_genes` | signalling-DE gene heatmap | ❌ **blocked** — its input comes from `communication_deg_detection`, which needs R/tradeSeq | 0 |

⚠️ **The vector fields were drawn twice.** The first pass used the *function's* defaults
(`plot_method='cell'`, `scale=1.0`, `ndsize=1`, no `normalize_v`); they rendered without error but
were illegible at 3k–26k cells — invisible arrows, a colour scale swamped by near-zero values.
`scripts/comparators/commot/plot_commot_vf.py` re-drew all 364 with the **tutorial's** argument
set, reading the persisted `adata_commot.h5ad` so **no OT re-run was needed** (4 min for all 13
cores, 0 failures). `scale` is the one value that cannot be copied verbatim — it reaches
`quiver(scale_units='x')`, so arrow length is `|v|/scale` in *data* units, and the tutorial's
`0.00003` is tuned to ~9,000-unit Visium pixel coordinates against our ~1,800–2,400 µm cores.
Rescaled per core to hold arrow length at the tutorial's fraction of the field
(1.1e-4 – 1.5e-4). Same class of unit trap as `dis_thr`; see `DEVIATIONS.md`.

⚠️ **Correction to the 2026-08-01 write-up.** It stated the vector field "wants an H&E background
we do not have". That is wrong against the installed source: `plotting/_plotting.py:35` has
`background: str = "summary"` as the **default**, and no image is involved. Confirmed empirically
— 14 vector fields per core now render with no background image. Those plots were skipped for no
valid reason.

### Multi-sample / differential mode

**None native.** COMMOT has no cross-sample differential test; `communication_deg_detection`
relates signalling to gene expression *within* a sample. A grade contrast would have to be
hand-rolled from per-core outputs.

### Gotchas

- **Unusable on numpy ≥ 2** — `np.Inf` as a module-level default in
  `_optimal_transport/_usot.py`. Requires `numpy<2`.
- **A dense N×N distance matrix is materialised** → the full slide would need **80.3 GB**.
  Per-core is mandatory, not optional.
- **The tutorial's prose contradicts its own units** — `dis_thr=500` is described as "500 µm"
  but Visium coordinates are full-resolution pixels at 0.73 µm/unit, so the real constraint is
  **365 µm**. Copying 500 onto micron data gives a 37% wider neighbourhood.
- **`min_cell_pct=0.05` does not transfer to single cells** — 5% of multi-cell Visium spots vs 5%
  of individual cells retains 0.9–1.8% of pairs instead of the tutorial's 20.9%.
- **Runtime is not predictable from size** — the OT solver runs to convergence, and iteration
  count varies per core. This survives the re-run: at a uniform 671 pairs, core 10 (17,435 cells)
  took **91.1 min** while core 1 (26,456 cells) took **47.2 min**, and core 13 (9,126) took
  46.8 min against core 8's (9,300) 55.6 min. Do not project runtime from cell count.
- **`ct.pl.plot_cluster_communication_dotplot` cannot run on a current matplotlib/seaborn.**
  commot 0.0.3 was written against matplotlib <3.9 and seaborn <0.13; installed are **3.10.9**
  and **0.13.2**. Two independent breakages: `_plotting.py:788` reads `legend.legendHandles`
  (removed in matplotlib 3.9 → `legend_handles`), and past that seaborn 0.13 hands back `Line2D`
  rather than `PathCollection`, so `set_edgecolor` does not exist. Recorded, not patched —
  fixing the second needs the function rewritten. `ours_dotplot_top.png` substitutes.
- **`plot_cluster_communication_network` needs pygraphviz** (via
  `networkx.drawing.nx_agraph.to_agraph`), which is not a COMMOT dependency. Installed 1.14 into
  `comp-commot`; without it the call raises `ImportError` and the plot is silently lost to a
  try/except.
- **The local clone `/Users/jiayifan/tansey_lab/COMMOT/` is NOT the installed package** — it
  carries two post-0.0.3 upstream commits in `_optimal_transport/_cot.py` (`np.Inf` fix; sparse
  cost-matrix support). `tools/_spatial_communication.py` is identical. Verify signatures against
  the installed package, not the clone. Corollary: **the numpy≥2 blocker below is fixed upstream**,
  just not in the release we run.
- **The OT is normalised per run.** `_cot.py:269-271` divides both marginals by
  `max_amount = max(S.sum(), D.sum())` and `:335` multiplies the plan back, so units are restored
  — but `eps_p` and `rho` are fixed constants acting on the *normalised* masses, so the shape of
  the solution depends on which cells are in the run. Each core is normalised by its own
  constant: **ranks are comparable across cores, magnitudes are not.**
- No LR-pair-level significance; magnitudes only. Cell-type-level p-values exist for whatever
  `cluster_communication` was pointed at — and pointing it at the *alphabetically* first pathways
  instead of the strongest ones is an easy and silent mistake (it is what the 2026-08-01 run did).

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `results/comparators/commot/GBM/cellchatdb2/<core>/` | ✅ **296.1 min, 13/13 cores**, zero failures | 2,259 input pairs → **671 used in every core** (`--filter-scope global`); 771 obsp keys/core (671 pairs + 99 pathways + total); obsp up to **63.8M** nonzeros; peak RSS **11.5 GB**; **2.4 GB**, **533 figures** (41/core), 262 CSVs (20/core), 13 `adata_commot.h5ad` |
| GBM | `cellchatdb2` | `.../cellchatdb2_prefix_backup_20260810/` | 🗄 superseded | the 2026-08-01 run: 123.4 min, 11/13 cores, 51–217 pairs/core, 54 plots. Kept for provenance; **its `per_split_summary.csv` ranks are wrong** (see the defect table in `DEVIATIONS.md`) |
| GBM | `default` | — | ❌ not run | bundled CellChatDB v1 tier deferred (user scoped this pass to `cellchatdb2` only) |
| LUAD | both | — | ❌ not run | deferred |

**Cores 2 and 6 are no longer lost.** Under the old per-core filter they retained zero pairs and
were reported as a method limitation on sparse tissue. Both are **low-grade**, so dropping them
left 7 high + 4 low against the TMA's true **7 high + 6 low** — a grade-biased analysis set.
Evaluating the filter once on all 100,197 cells restores full coverage.

**The pair-count confound is gone.** `corr(n_cells, n_pairs_used) = 0.819, p = 0.002` described
the old run; with a single global pair set every core uses exactly 671, so **ranks are now
comparable across cores**. The density–grade confound documented for SpatialDM still applies to
*magnitudes*, and a second reason magnitudes are not cross-core comparable is that the OT is
normalised per run (`_cot.py:269-271`, `max_amount`) — see `DEVIATIONS.md`. **Compare ranks, not
magnitudes.**

### Requested LRIs — the strongest signal of the benchmark so far

Every core now ranks the same 671 pairs, so these ranks share one denominator.

| LRI | tested | median rank | range | |
|---|---|---|---|---|
| **GRN → SORT1** | **13/13 cores** | **2** (of 671) | **1–12** | #1 in 5 cores, top-3 in 9, **top-5 in 12** |
| ANXA1 → FPR1 | **13/13 cores** | 21 | 4–126 | never #1; top-5 in 2 |

COMMOT puts GRN→SORT1 in the top 5 of **12 of 13 cores** out of 671 candidates. The one
exception is core 2 (rank 12) — the smallest core, 819 cells. That is far higher than CytoSignal
(66/895) or stLearn (21/526) place it, but the statistic differs: COMMOT ranks by **transported
signal magnitude**, which favours abundantly expressed ligands, whereas the others rank by
significance against a null. The agreement worth noting is directional — all four methods place
GRN→SORT1 well above ANXA1→FPR1 — and COMMOT's margin should not be read as stronger evidence.

ANXA1→FPR1 is the clearer gain from the re-run: previously **untestable in core 14** and absent
from 3 cores' pair sets, it is now scored in all 13.

**New in this run — the only significance test COMMOT can give a named pair.**
`cluster_communication` was additionally run with `lr_pair=` for both requested LRIs, so each core
has `cluster_comm_{GRN_SORT1,ANXA1_FPR1}.csv` + `cluster_pval_*` — cell-type × cell-type
permutation p-values (`n_permutations=100`, `random_seed=0`). This is what the per-pair magnitude
ranking cannot provide.

Rank-1 LR pair per core: **FGF1–FGFR2** (7 cores), **GRN–SORT1** (5), WNT5A–FZD3 (1). Most
frequent members of a core's top 5: **GRN–SORT1 (12/13)**, FGF1–FGFR2 (7), PDGFA–PDGFRA (6),
ANGPTL2–TLR4 (5), PDGFC–PDGFRA (5), CSF1–CSF1R (4).

Top-5 pathways by frequency across the 13 cores: **FGF 13/13, PDGF 13/13**, Glutamate 8, GAS 6,
BMP 6, COMPLEMENT 5, ANGPTL 3, **GRN 3**, ncWNT 3, GALECTIN 2, IGF 2, PROS 1. Glutamate
persisting in 8 of 13 cores corroborates SpatialDM's independent finding that glutamatergic
signalling dominates core 13, and vindicates including v2's `Non-protein Signaling` class.

### What COMMOT actually found in this tumour — plain reading

COMMOT returns a **ranked list of how loudly each ligand–receptor pair is signalling**. It says
nothing about who is talking to whom (that needs `cluster_communication`) and nothing about
programs. Read at that level, four things come out.

**1. Growth-factor axes dominate, unambiguously.** **FGF and PDGF are top-5 pathways in 13 of 13
cores**, no exceptions — concretely `FGF1–FGFR2` (rank 1 in 7 cores), `PDGFA–PDGFRA` and
`PDGFC–PDGFRA` (top-5 in 6 and 5 cores). This is the least surprising and most reassuring result
in the benchmark: PDGFRA is a canonical glioma driver, and a method with no knowledge of this
tumour recovered both axes from 671 candidates. It is the sanity check that the run is not
mis-specified.

**2. Glutamatergic signalling is the most interesting hit.** **Glutamate is top-5 in 8 of 13
cores.** Three reasons it matters: it comes from CellChatDB v2's `Non-protein Signaling` class,
which we added on top of the tutorial's Secreted-only restriction — so the deviation paid off; it
**independently corroborates SpatialDM**, which ranked glutamatergic pairs as core 13's entire
top-6 under a completely different statistic; and it points at neuron–glioma synaptic signalling,
where glutamate drives proliferation and invasion.

**3. Myeloid cells are the communication hub; tumour cells are comparatively quiet.** Summing
COMMOT's transported signal over all 671 pairs (`total-total`, mean within-core percentile):

| Cell type | received | sent |
|---|---|---|
| **mGAM** | **0.804** | **0.777** |
| Vascular | 0.759 | 0.734 |
| non-mGAM | 0.725 | 0.694 |
| Lymphoid | 0.677 | 0.648 |
| AC-/MES-/NPC-/OPC-like | 0.455–0.483 | 0.400–0.474 |

mGAM tops **both** directions while every tumour state sits near the median. The molecular detail
agrees: `GRN–SORT1` (progranulin→sortilin, myeloid) is top-5 in 12/13 cores, `CSF1–CSF1R` in 4,
and `ANGPTL2–TLR4` and COMPLEMENT (`C3–C3AR1`) recur throughout. ⚠️ **Discount this.** COMMOT's
magnitudes grow with the number of cells inside the 365 µm radius, so "myeloid is the hub" is
partly a density artifact and COMMOT cannot separate the two. See the motif-1 section below,
where this confound is measured directly.

**4. The requested pair is near the top everywhere, but carries no grade signal.** `GRN→SORT1`
median rank **2 of 671**, top-5 in 12/13 cores, #1 in 5 — the sole exception is core 2, the
smallest (819 cells). A method sharing no mathematics with ALARMIST puts that pair at the front
on its own. `ANXA1→FPR1` median 21, range 4–126, never #1. **Neither pair's rank differs by
grade** (`GRN→SORT1` rank 1 in high vs 2 in low, p = 0.413), whereas ALARMIST's motif-1 loading
does (p = 0.022). COMMOT sees the interaction but not its clinical association.

### What the vector-field figures are, and whether they mean anything

**What they are.** For one LR pair (or pathway) COMMOT holds a cell × cell transport plan.
`communication_direction` collapses it to one 2-D arrow per cell — the k=5-weighted average
*direction* in which that cell's signal flows — and `plot_cell_communication(plot_method='grid')`
interpolates those onto a regular grid. **One arrow = "in this patch of tissue, this pathway's
signal is on average flowing that way."**

⚠️ **The single most common misreading: the arrow is a SPATIAL direction, not a cell-type
direction.** It says "towards the upper right of the image", **not** "from mGAM to MES-like". For
sender→receiver relationships between cell types you need `cluster_communication`'s matrix and
its permutation p-values, which is a different output entirely.

**Are they meaningful? On this dataset — presentational rather than analytical.** Four reasons:

1. **No test.** Neither an arrow's direction nor its length has a null, a p-value, or any notion
   of significance. Nothing here can be reported as "significantly directional".
2. **They largely trace tissue geometry.** Transported mass flows from dense to sparse regions.
   The `total-total` analysis shows COMMOT's signal distribution is dominated by abundance and
   density, so the arrows plausibly track a density gradient rather than biology.
3. **Edge artifacts.** In the figures inspected (core 13) the longest arrows sit at the tissue
   perimeter, where one side of the interpolation kernel is empty — geometry, not signal.
4. **Wrong setting.** The tutorial's mouse-brain Visium section has strong anatomical gradients
   for the arrows to trace. A 1–2 mm TMA punch of relatively homogeneous tumour has no comparable
   macro-gradient.

**What they are still good for:** a fast visual check of whether a pathway's signal is diffuse or
concentrated in a few sources, and — for the `_ct` variant, whose background is cell type — of
what the high-signal regions are made of. **If only one background is kept, keep `_ct`.** They
also document that the authors' standard workflow was run end to end, which has benchmark value.
**Do not** cite them as evidence that communication is directional; use `cluster_communication`
for that.

### Does COMMOT independently corroborate ALARMIST motif 1? — mostly **no**

`scripts/comparators/commot/compare_motif1_commot.py` → `results/comparators/commot/GBM/vs_alarmist/`.
Joins COMMOT's per-cell sent/received amounts to `results/GBM/single_cell/cell_loadings.npy`
(100,197 × 20, motif 1 = column 1). All 13 cores, all 100,197 cells; the join key is `obs_names`
and every COMMOT cell id resolves. **Alignment guard:** the script refuses to run unless motif 1
peaks on mGAM — it does (5.15e-3), and motif 1 is the only column among those checked that does
(others peak on AC-like / MES-like / NPC-like), which corroborates the row order *and* the column
index. *(Note motif 1's margin over Lymphoid, 5.12e-3, is razor-thin — motif 1 is high in mGAM,
Lymphoid and Vascular alike, so "the mGAM motif" is not a clean per-cell-type label.)*

Motif 1 is a bidirectional loop, which makes a falsifiable prediction: GRN→SORT1 should be sent
by mGAM and received by MES-like; ANXA1→FPR1 the reverse.

**The decisive control kills the naive reading.** COMMOT's grand totals over all 671 pairs
correlate with motif-1 loading **better than any motif-1-specific quantity does**:

| COMMOT quantity | median Spearman vs motif 1 (13 cores) | best-matching motif |
|---|---|---|
| **`s-total-total`** (all signalling sent) | **0.167** | motif 1 |
| **`r-total-total`** (all signalling received) | **0.156** | motif 1 |
| `composite-mGAM-side` (`s-GRN-SORT1` + `r-ANXA1-FPR1`) | 0.119 | motif 1 |
| `s-ANXA1-FPR1` | 0.089 | motif 9 |
| `r-ANXA1-FPR1` | 0.084 | motif 1 |
| `s-GRN-SORT1` | 0.079 | motif 1 |
| `r-GRN-SORT1` | 0.024 | motif 6 |
| `r-`/`s-FGF1-FGFR2` *(unrelated-pair control)* | −0.014 / 0.004 | motif 6 |

So motif-1 loading tracks **"how much signalling this cell does at all"** better than it tracks
COMMOT's estimate of the two LRIs that *define* motif 1. The FGF1–FGFR2 control does its job
(motif 1 is not correlated with everything — that pair matches motif 6), but the total-total
control shows the shared variance is a **generic communicativeness / local-density axis**, not
pair-specific agreement. All effect sizes are small (|ρ| ≤ 0.17).

**Directional test, before and after subtracting the matching `total-total` profile.** mGAM is
COMMOT's most communicative cell type overall (top of *both* grand totals: 0.804 received, 0.777
sent), so raw "mGAM is #1" is uninformative:

| Prediction | raw rank | baseline-corrected rank | residual |
|---|---|---|---|
| `s-GRN-SORT1` → mGAM | 1/9 ✅ | **7/9 ❌** | −0.121 |
| `r-GRN-SORT1` → MES-like | 4/9 ❌ | 2/9 (top: Glial-Neuronal) | +0.012 |
| `s-ANXA1-FPR1` → MES-like | 3/9 ❌ | **1/9 ✅** | +0.075 |
| `r-ANXA1-FPR1` → mGAM | 1/9 ✅ | **8/9 ❌** | −0.187 |

Correction **reverses** the raw reading: the mGAM sides, which looked confirmed, fall below what
mGAM's general communicativeness predicts; the MES-like sides, which looked refuted, rise. Only
one prediction lands first — **MES-like as the ANXA1 sender** — and even that is suggestive
rather than decisive, because MES-like sits mildly above baseline for *all four* quantities
(+0.012 … +0.075), so its positivity is not direction-specific. *(Minor conservatism: the
baseline includes the 2 tested pairs among its 671, so it is very slightly self-subtracting.)*

**Grade.** With the core as the replicate (7 high vs 6 low, Mann–Whitney), mean motif-1 loading
is higher in high-grade cores (9.1e-4 vs 8.9e-5, **p = 0.022**) — reproducing the grade
association CLAUDE.md records, and *not* explained by mGAM abundance (`frac_mGAM` p = 0.836).
**COMMOT does not recover it:** its GRN→SORT1 rank is 1 in high-grade vs 2 in low-grade cores,
p = 0.413.

**Bottom line — the two methods agree about *which LRIs matter*, not about *which cells are doing
it*.** COMMOT puts GRN→SORT1 in the top 5 of 12/13 cores out of 671 candidates, which is strong
independent support for the pair's importance. It provides **no independent per-cell support for
motif 1's cell-type attribution**, and what per-cell agreement exists is explained by a density /
communicativeness confound. This is a real limitation of using COMMOT as a validator: its
magnitudes scale with how many cells sit inside the 365 µm transport radius, so the most
abundant, most densely packed cell types dominate every pair. Treat the concordance as
LRI-level, and state the cell-type attribution as ALARMIST's claim, not a cross-validated one.

Outputs: `motif_vs_commot_spearman_{percore,median}.csv`,
`celltype_signal_percentile_{percore,mean,baseline_corrected}.csv`, `core_summary.csv`,
`summary.json`, and `{motif_vs_commot_heatmap,directional_celltype_test,spatial_concordance_core13}.{png,pdf,svg}`.

### Downstream impact — half of COMMOT's chain is unreachable on this machine

COMMOT does have a motif→gene analogue. The tutorial's chain (`visium-mouse_brain.ipynb`) is:

```python
df_deg, df_yhat  = ct.tl.communication_deg_detection(adata, database_name, pathway_name, summary='receiver')
df_deg_c, df_y_c = ct.tl.communication_deg_clustering(df_deg, df_yhat, deg_clustering_res=0.4)
top_de_genes     = ct.pl.plot_communication_dependent_genes(df_deg_c, df_y_c, top_ngene_per_cluster=5, return_genes=True)
df_impact        = ct.tl.communication_impact(adata, database_name, pathway_name=..., tree_combined=True,
                       method='treebased_score', tree_ntrees=100, tree_repeat=100, tree_method='rf',
                       ds_genes=top_de_genes, bg_genes=500, normalize=True)
ct.pl.plot_communication_impact(df_impact, summary='receiver')
```

**Step 1 needs R.** `communication_deg_detection` imports `rpy2` + `anndata2ri` **inside its own
body** and drives `tradeSeq::fitGAM` + `associationTest` and
`clusterExperiment::clusterExpressionPatterns`. It is the **only** function in this otherwise pure
Python package that touches R — `pip show commot` declares no R dependency, the imports are
function-local by design, and `import commot` plus every other call works without it. On this
machine `rpy2` and `anndata2ri` are absent and **none of the four R installs has tradeSeq**
(R 4.4.2 at `/usr/local/bin/R`; R 4.3.3 in `comp-cellchat`, `comp-niches`, `comp-cytosignal`).
Consequently `communication_deg_clustering` and `plot_communication_dependent_genes` — both pure
Python — are starved of input, and `ds_genes` cannot come from the tutorial's own source.

> The docstring's *"tradeSeq 1.0.1 with R 3.6.3"* is what the authors **tested**, not a declared
> minimum; current tradeSeq targets R ≥ 4.3. **`scripts/comparators/commot/SETUP_tradeSeq.md`** is
> the plan for enabling it in an isolated `comp-commot-r` env (nothing in it has been executed):
> install recipe, a pre-flight check, the **five specific places the modern stack can break** and
> what each failure looks like, the data-contract traps, and a runtime plan.
> **`run_commot_deg.py`** implements the authors' full chain against it — also never executed,
> written from the installed source, syntax-checked only.

**Steps 4–5 are runnable**, because `communication_impact`'s docstring sanctions a substitute:
*"A list of genes … for example, the highly variable genes."* `run_commot_impact.py` runs it with
every tutorial argument value intact and two `ds_genes` variants — `hvg` (the docstring's
substitute, method-faithful) and `alarmist` (top genes from
`results/GBM/impact/motif_1_celltype_mGAM_de_results.csv`, so both methods are asked about the
*same* genes and the answers are comparable).

> **Status 2026-08-11 — the substitute run was deliberately stopped.** 9 of 13 cores of the
> `GRN × alarmist` variant completed (~50 min) before the user elected to **install the R
> dependency and run the authors' actual chain** instead. The partial output is kept, clearly
> marked, at `results/comparators/commot/GBM/impact/` (`PARTIAL_RUN_DO_NOT_USE_AS_FINAL.md`); it
> has no `run_manifest.json`/`impact_summary.csv`, which is the marker of an incomplete run. The
> `hvg` control and the proposed motif-14 control never ran. **Nothing here is a finished result.**
>
> One observation from the completed cores, recorded but **uncontrolled**: on core 6 the
> receiver-side scores sat far above the 0.5 null (`r-GRN-SORT1` median 0.888, 28/30 genes > 0.8)
> while the sender side sat on it (`s-GRN-SORT1` median 0.519, 2/30 > 0.8) — i.e. how much
> GRN→SORT1 a cell *receives* tracked the motif-1 genes, how much it *sends* did not. Direction is
> what one would expect (the receiving cell is the one whose transcriptome should respond), but
> **the confound is not excluded**: motif-1's mGAM genes are myeloid genes and GRN→SORT1 receipt
> concentrates in myeloid-dense regions, so any mGAM-associated gene set might score the same. The
> control that would settle it — `ds_genes` from motif 14's mGAM impact table, the same method and
> cell type but only 6/30 gene overlap — was designed and not run. **Do not cite the agreement
> until it exists.**

**What `communication_impact` returns:** `df_impact`, rows = `s-<name>` and `r-<name>` for every
LR pair in the pathway plus the pathway aggregate, cols = the `ds_genes` plus `'average'`. With
`method='treebased_score'` and `tree_combined=True` each target gene gets one random forest whose
predictors are the communication features **plus ~500 background genes**, and the score is the
communication features' relative importance. Like ALARMIST's GLM, it **excludes the LR genes
themselves** (`exclude_lr_genes_list`) — the same anti-circularity guard.

Three implementation facts worth recording:

- **It requires `adata.raw`.** `run_commot.py` drops `.raw` before writing `adata_commot.h5ad`
  (to halve the file), so `run_commot_impact.py` rebuilds it from the source h5ad's
  `layers['counts']` — verified to be genuine integer counts (max 128), i.e. exactly what `.raw`
  held during the run. Still no OT re-run.
- **The stale `uns['log1p']` trap resurfaces:** `communication_impact` internally normalises and
  log1p's the background, and the key left by the main run makes scanpy warn "already
  log-transformed". Popped; behaviour-neutral.
- **`pathway_name='GRN'` is a degenerate, and convenient, case.** The GRN pathway contains exactly
  **one** LR pair — `GRN→SORT1` — so `df_impact` is 4 × (N+1) and the analysis is literally
  "what does GRN→SORT1 signal explain". The figure needs `cluster_knn` below the row count:
  `plot_communication_impact` clusters the `n_pairs + 1` rows of each summary, so the default
  `cluster_knn=5` raises `IndexError`. The tutorial hits the same wall and passes **`cluster_knn=2`**
  ("here we only have two LR pairs in PSAP" → 3 rows); GRN's 2 rows need **`cluster_knn=1`**,
  verified working. The runner now falls back 5 → 2 → 1. *(An earlier note here said the GRN
  figure was impossible — wrong; it just needed the tutorial's own knob turned further.)*
  `plot_communication_impact` takes `df_impact`, not the AnnData, so figures can be regenerated
  from the saved CSVs without re-fitting.

- **The score is a percentile against background genes, and its null is 0.5 — not 0.**
  `_utils/_similarity.py:85` fits `n_repeat` random forests on `[signal feature | ~500 background
  genes]`, takes the signal's rank among the feature importances, and returns
  `mean((n_bg - rank) / n_bg)`. So 0.97 means "more important than ~97 % of the background genes"
  (≈ rank 15 of 500), and a useless feature lands near **0.5**, not 0. This is why the tutorial's
  own `df_impact_PSAP` table sits entirely in 0.78–0.999 and looks uniformly high. **Always read
  these against 0.5.**

Cost: **3.9 min per core at 3,092 cells / 30 genes**, scaling with cells × genes;
`RandomForestRegressor` already runs at `n_jobs=-1`, so cores are saturated and running variants
in parallel buys nothing.

### What COMMOT's output objects actually are — measured, beside ALARMIST's

All shapes below were **measured on disk**, not read off documentation
(COMMOT: `results/comparators/commot/GBM/cellchatdb2/13/adata_commot.h5ad`, core 13, 9,126 cells —
one such object per core, ×13; ALARMIST: `results/GBM/`). 51 of them were independently re-measured
by a second pass; the three discrepancies it found are footnoted and none affect the picture.

#### The one structural fact that explains everything else

**ALARMIST has a latent axis; COMMOT does not.** Every ALARMIST matrix carries **motifs**
(K = 20) on one side — a learned axis that does not exist in the data. Every COMMOT object is
indexed by **observed LR pairs** (671) and by **cells**. There is no matrix anywhere in COMMOT
whose axis corresponds to ALARMIST's motif axis, and no way to construct one from its outputs.

A second, subtler difference: **ALARMIST's "LRI" axis already contains the cell-type direction.**
Its 25,271 columns are `celltype1|celltype2|ligand|receptor|signaling_type` (e.g.
`AC-like|AC-like|TGFB1|TGFBR2_TGFBR1|autocrine`). COMMOT's 671 LRIs are *only* ligand×receptor;
cell type enters exclusively at the `cluster_communication` step. That is why ALARMIST's factor
matrix can express "mGAM→MES-like via GRN→SORT1" as a single weighted feature and COMMOT cannot.

#### ALARMIST — `results/GBM/`

| Object | Shape (measured) | Axes |
|---|---|---|
| `bptf/lri_factors.npy` — **V** | **20 × 25,271** | motif × LRI-column *(note: transposed relative to the loadings)* |
| `bptf/patch_loadings.npy` — **W** | 13,113 × 20 | patch (50 µm tile) × motif |
| `single_cell/cell_loadings.npy` — **U** | **100,197 × 20** | cell × motif |
| `bptf/lri_motifs.csv` | 505,420 × 15 | long form of V (25,271 LRI × 20 motifs); carries `factor`, `factor_rescaled` and the prevalence-normalised **`score`** — rank on `score` |
| `impact/motif_<k>_celltype_<ct>_de_results.csv` | **180 files** (20 motifs × 9 cell types), 266,746 rows total, 4,018 distinct genes | one row = one **gene**; the (motif, cell type) key is in the *filename*. So the impact is really **gene × motif × cell type**, stored as 180 tables — not a flat gene × motif matrix |
| `impact/univariate_de_long.csv` | 266,746 × 6 | the same 180 tables stacked **with the cell-type label dropped** ⚠️ a (gene, motif) key repeats up to 9×; only 68,080 rows are unique on (gene, motif) |
| `markers/exclusion_matrix.csv` | 5,119 × 9 | gene × cell type, boolean |
| `patch_lri_matrix.npz` / `single_cell/cell_lri_matrix.npz` | 13,113 × 25,271 / 100,197 × 25,271 | the model **inputs**: integer co-occurrence counts (not binary) |

#### COMMOT — per core

| Object | Shape (measured, core 13) | Axes |
|---|---|---|
| `obsp['commot-cellchat-<lig>-<rec>']` | **671 matrices, each 9,126 × 9,126** (15.3 M nnz total) | a 3-way tensor: **LR pair × sender cell × receiver cell**. Entry = mass transported from cell *i* to cell *j* for that pair. **This is COMMOT's primary output.** |
| `obsp['commot-cellchat-<pathway>']` | 99 × (9,126 × 9,126), 11.4 M nnz | pathway × sender × receiver — **derived sums** of the member pairs, not independent estimates |
| `obsp['commot-cellchat-total-total']` | 9,126 × 9,126, 7.5 M nnz (9.0 % dense) | sender × receiver, all pairs collapsed |
| `obsm['commot-cellchat-sum-sender']` | **9,126 × 771** | cell × channel, where channel = 671 LR pairs + 99 pathways + 1 total. Row-sums of the plans above |
| `obsm['commot-cellchat-sum-receiver']` | 9,126 × 771 | same, column-sums (incoming) |
| `uns['commot_cluster-cell_type-cellchat-<key>']` | **7 keys × (9 × 9 matrix + 9 × 9 p-values)** | sender cell type × receiver cell type, per key. Computed only for the 5 top pathways + the 2 requested pairs — **the only place COMMOT produces a p-value** |
| `obsm['commot_{sender,receiver}_vf-…']` | 14 arrays, each 9,126 × 2 | cell × (dx, dy); plotting only, and only for those same 7 keys |
| `uns['commot-cellchat-info']['df_ligrec']` | 671 × 3 | LR pair × (ligand, receptor, pathway) — a flat vocabulary, never a loading matrix |

`obsp` holds 773 keys in total: 671 + 99 + 1 COMMOT matrices (34,279,269 nonzeros, matching
`per_split_summary.csv`) plus `connectivities` and `distances` carried in from the input h5ad.

#### Nearest-analogue map

| ALARMIST | COMMOT | Comment |
|---|---|---|
| `cell_loadings` **100,197 × 20** (cell × motif) | `sum-sender` / `sum-receiver` **9,126 × 771** per core (cell × LRI) | closest pair. But: latent motif vs observed LRI; one number vs a sent/received split; one object over all cells vs 13 objects with **no shared cell axis** |
| `lri_factors` **20 × 25,271** (motif × LRI) | — **nothing** | COMMOT never groups or weights LRIs. This is the missing piece that makes motifs impossible |
| impact **gene × motif × cell type** (180 tables) | — **nothing produced** | `communication_impact` / `deg_detection` exist but are per-pathway within one sample, and were not run |
| cell-type communication heatmap per motif (derived from V) | `cluster_communication` **9 × 9** + p-values, per key | ALARMIST's is per *motif* and has no p-value; COMMOT's is per *LR pair or pathway* and does. Complementary, not redundant |
| — **nothing** | `obsp` transport plans **671 × 9,126 × 9,126** | ALARMIST produces no cell-to-cell edges at all. This is COMMOT's unique contribution |
| patch loadings **13,113 × 20** | — **nothing** | no spatial-unit abstraction in COMMOT; its spatial model is the 365 µm cutoff on the transport cost |

**In one line:** ALARMIST factorises into *programs* and is fine-grained along the program axis;
COMMOT resolves individual *cell-to-cell edges* and is fine-grained along the cell axis. Neither
can be converted into the other.

*Three re-measurement corrections, none material:* `run_manifest.json`'s nested `versions` object
has **8** entries, not 9; `vf_replot_manifest.json` has **8** top-level keys, not 7; and 179 of the
180 impact files have 6 columns while `motif_6_celltype_Lymphoid_de_results.csv` is header-only
with 4 columns and 0 rows.

### COMMOT vs ALARMIST — what each can and cannot do

The two are not competing estimators of the same quantity; they answer different questions. The
gaps below are structural, not accuracy differences.

**Only ALARMIST can do these.**

1. **Discover sets of LRIs that switch on together.** COMMOT scores each pair and returns a
   671-row ranking. It has no representation in which `GRN→SORT1` and `ANXA1→FPR1` could belong
   to one program — co-activity is not a concept the model contains. ALARMIST's factorization is
   exactly that concept, which is why the motif-1 "mGAM ⇄ MES-like loop" hypothesis could only
   have come from ALARMIST.
2. **Give each cell a per-program score.** ALARMIST projects motifs onto cells (`U`, cells ×
   motifs) and binarizes to ON/OFF states. COMMOT's per-cell granularity is *per LR pair*; there
   is no "this cell is participating in program k".
3. **Compare across samples or conditions.** **COMMOT has no multi-sample mode at all.** We ran
   13 independent OT problems, and because each is normalised by its own constant only *ranks*
   are comparable across cores. ALARMIST fits one joint model so motif *k* means the same thing
   in every sample. The consequence is measurable: motif-1 loading separates high from low grade
   (p = 0.022, core as replicate) while COMMOT's rank for the same pair does not (p = 0.413).
   That is not COMMOT being wrong — it has no structure that can carry the question.
4. **Relate a signalling program to downstream transcription.** ALARMIST's Poisson GLM regresses
   raw counts of non-LR genes on the motif loading per cell type. COMMOT's `communication_impact`
   / `deg_detection` relate a *single pathway's* signal to genes within one sample; there is no
   program-level analogue. (We did not run them — downstream interpretation, and each needs its
   own fit.)
5. **Everything built on 1–4:** motif-derived gene signatures for external survival cohorts,
   combinatorial motif-state programs, motif-usage niches. None of these have a COMMOT entry point.
6. **Analyse the specimen as one object.** COMMOT materialises a dense N×N distance matrix —
   80.3 GB for 100,197 cells on a 36 GB machine. It *physically cannot* hold the slide at once,
   which is part of why (3) is impossible for it.

**Only COMMOT can do these.**

- **Competition for finite mass.** Collective optimal transport makes receptors saturable, so a
  nearby strong sender crowds out a distant one. ALARMIST's counting model has no competition.
- **An explicit physical distance constraint** (365 µm here) applied pair by pair, against
  ALARMIST's much coarser 50 µm patch binning.
- **Cell-to-cell edges.** COMMOT yields who sent how much to whom (`obsp`, persisted here).
  ALARMIST produces no such object.

**Net:** the concordance between the two stops at *which LRIs matter*. COMMOT independently puts
`GRN→SORT1` in the top 5 of 12/13 cores out of 671 — real, useful support for the pair. It cannot
confirm that motif 1 is a program, cannot confirm the program belongs to mGAM (the baseline-
corrected test argues against it), and cannot see the grade association. Report cross-method
agreement at the LRI level; report cell-type attribution and clinical association as ALARMIST's
claims.

### Methods paragraph

> For COMMOT (v0.0.3), we performed spatial communication inference with the
> `spatial_communication` function, following the authors' tutorial. Counts were normalised and
> log-transformed, and ligand–receptor pairs were taken from CellChatDB v2 restricted to the
> diffusible signalling classes, since COMMOT applies a single distance constraint to all pairs;
> heteromeric complexes were handled natively with `heteromeric=True` and the minimum rule.
> Pairs were filtered with `filter_lr_database` requiring at least 100 cells expressing each
> side, an absolute criterion chosen because the tutorial's 5%-of-spots threshold does not
> transfer to single-cell resolution; the filter was evaluated once across the whole specimen and
> the resulting 671 pairs applied to every core, so that pair sets and hence rankings are
> comparable between cores. The spatial distance constraint was set to 365 µm, reproducing the
> physical constraint of the authors' own example after converting their Visium pixel coordinates
> to microns. Cell-type-level communication and its permutation p-values were obtained with
> `cluster_communication` (`n_permutations=100`, fixed random seed) for the five pathways carrying
> the most received signal in each core and, separately, for each ligand–receptor pair of interest.
> Signalling directions were interpolated with `communication_direction` and visualised with
> `plot_cell_communication`. Each of the 13 tissue microarray cores was analysed separately, as
> COMMOT materialises a dense pairwise distance matrix that is intractable at whole-slide scale.

---

## LIANA+ — Python, v1.8.1, env `comp-liana`

Tutorial: `/Users/jiayifan/tansey_lab/liana-py/docs/notebooks/bivariate.ipynb`.
Call contract: `scripts/comparators/liana/NOTES.md`. Deviations: `DEVIATIONS.md`.

### Core algorithm

LIANA+'s spatial mode computes a **spatially-weighted bivariate similarity** between a ligand
and a receptor — by default a weighted cosine, one of six available metrics
(`li.mt.bivariate.show_functions()`). One call returns two levels simultaneously:

- **Local** — a score *per cell per LR pair* (`lrdata.X`), with permutation p-values from
  shuffling cell labels (`layers['pvals']`) and categorical labels (`high-high`, `low-low`, …)
  in `layers['cats']`.
- **Global** — per-pair summaries in `lrdata.var`: `mean`, `std`, and **bivariate Moran's R**
  (Lee's statistic, the two-variable extension of Moran's I) with `morans_pvals`.

**LIANA+ is the only method in this benchmark that natively returns both a per-cell score and a
per-pair global statistic with p-values from a single call.** Note what Moran's R measures: not
"is co-expression high" but "are ligand and receptor *spatially arranged* together" — a pair
that is ubiquitous and uniform scores low (the tutorial makes this point with TIMP1–CD63).

### Spatial model

Gaussian kernel `exp(-d²/(2·bandwidth²))`, weights below `cutoff` zeroed, on a
`max_neighbours`-capped KNN graph. Effective radius = `bandwidth × sqrt(-2·ln cutoff)`.

**Ours: bandwidth 13.1454 µm, cutoff 0.1 → effective (support) radius 28.2 µm**, **median 14**
neighbours per cell, **max 52** (`cellchatdb2/run_manifest.json`, verified over all 100,197
cells). The value comes from an **equal-area correspondence to a 50 µm square patch**; the full
derivation, and the contract deviation it entails, are in *Spatial bandwidth — final parameter
table* and *Contract deviations* in the second LIANA+ section below. Both branches (bivariate
and inflow) use the identical graph.

> **⚠️ Superseded — do not quote.** An earlier configuration used **bandwidth 18.75 µm →
> effective radius 40.2 µm, median 27 neighbours**, derived with the authors' own
> `li.ut.query_bandwidth` from the *bivariate* tutorial's criterion (*"a bandwidth of 150-200
> roughly includes 6 neighbours i.e. the first ring of neighbours"*). **That value is VOID** —
> "first ring" is a property of Visium's hexagonal lattice and has no referent on irregularly
> packed single cells (full argument below). **All five `run_manifest.json` under
> `results/comparators/liana/GBM/` that carry a `bandwidth` key at all — `cellchatdb2`,
> `cellchatdb2_inflow`, `cellchatdb2_morans`, `default`, `default_inflow` — record
> `bandwidth: 13.1454`, so nothing currently on disk was produced at 18.75**
> *(this read "all four" until 2026-08-04; there are now 16 `run_manifest.json` in that tree, and
> the count of those carrying the key is five — see the CD-1 row for the full audit)* — the only surviving trace of that run is the sensitivity table below,
> which is kept deliberately. The record is kept for provenance only; 18.75 / 40.2 µm / 27
> neighbours must not be reported as this method's spatial scale.

This is the **smallest neighbourhood of any method here** (SpatialDM 135 µm, CytoSignal 200,
stLearn 250, COMMOT 365), which is appropriate to single-cell resolution — the relevant range is
tens of µm, whereas on Visium it is ~100 µm because each unit is a 55 µm multi-cell spot. Note,
however, that unlike the other methods' figures **this one is not LIANA's own default**: LIANA
pins no default bandwidth at all, and ours was set by geometric reference to ALARMIST's patch
size. See *Contract deviations* below.

### LR database

Default is LIANA's **consensus** resource. `li.mt.bivariate` also accepts
`resource: pd.DataFrame` with `['ligand','receptor']`, and LIANA joins heteromeric subunits with
`_` — exactly our CellChatDB v2 encoding, so the v2 handover is **direct and lossless**, with
native complex support. Pairs are named `<ligand>^<receptor>`; `^` and `_` mean different things.

3,218 unique v2 pairs → 27 self-interactions removed → **131 pairs** survive `nz_prop = 0.02`.

#### ⚠️ Which CellChatDB export — a silent join trap against the ALARMIST GBM run

**The LIANA runs consumed the re-exported (post-2026-07-28)
`data/LRdatabase/CellChatDBv2.0.human.csv`.** Proven by complex-subunit ordering: the LIANA
outputs carry **450 × `TGFBR1_TGFBR2` and 0 × `TGFBR2_TGFBR1`**, and **633 / 633** of the run's
distinct LR keys match the current file (only **496 / 633** match
`CellChatDBv2.0.human.old.csv`).

**But `results/GBM/patch_lri_columns.csv` — the ALARMIST run this is compared against — has
210 × `TGFBR2_TGFBR1` and 0 reversed**, i.e. the **old** export. Neither run is wrong (CLAUDE.md
directs fresh runs at the current file, and the ALARMIST result predates the re-export), but the
two label heteromeric complexes in **different subunit order**, so a raw-string LR join between
them silently drops keys:

| join of LIANA's 633 LR keys against ALARMIST's 712 | shared |
|---|---|
| raw string | **495** |
| after canonicalisation | **545** |
| of ALARMIST's 205 heteromeric keys: raw / canonical | **94** / **144** |

So **50 heteromeric keys (24% of them) vanish from a naive join** — silently, with no error.
**Remedy: canonicalise each side before joining**, applying
`'_'.join(sorted(s.split('_')))` to the ligand and the receptor *separately*:

```python
def canon(lr):                      # 'TGFB1^TGFBR2_TGFBR1' -> 'TGFB1^TGFBR1_TGFBR2'
    l, r = lr.split('^')
    return '_'.join(sorted(l.split('_'))) + '^' + '_'.join(sorted(r.split('_')))
```

**`GRN^SORT1` and `ANXA1^FPR1` are identical in both files**, so the two requested LRs are
unaffected by this trap.

### Input

`layers['counts'] = X.copy()`; `sc.pp.normalize_total(target_sum=1e4)`; `sc.pp.log1p` — the
tutorial's exact recipe, applied to genuine raw counts.

### Workflow

| # | Call | Produces |
|---|---|---|
| 1 | `li.ut.query_bandwidth(coordinates, ...)` | bandwidth calibration curve |
| 2 | `li.ut.spatial_neighbors(adata, bandwidth, cutoff, kernel, set_diag)` | `obsp['spatial_connectivities']` |
| 3 | `li.mt.bivariate(adata, resource, local_name='cosine', global_name='morans', n_perms=100, mask_negatives=False, add_categories=True, nz_prop, use_raw=False)` | `lrdata` (cells × pairs) with local scores, p-values, categories, and global stats |

### Data outputs

| File | Shape | Meaning |
|---|---|---|
| `global_scores.csv` | 131 × 10 | ligand/receptor, `*_means`, `*_props`, **`morans`**, `morans_pvals`, `mean`, `std` |
| `local_scores.npz` | cells × pairs | per-cell local bivariate score (float32) |
| `local_pvals.npz` | cells × pairs | per-cell permutation p-value |
| `local_categories.npz` | cells × pairs | local category labels |
| `cell_meta.csv` | cells | coordinates + cell type |

### Image outputs

`run_liana.py` writes **35 PNGs** into `cellchatdb2/plots/` (27 at top level + 8 in `requested/`),
plus `bandwidth_query.csv`. Requested pairs are segregated into `plots/requested/` so they cannot
be mistaken for the method's own ranking (`SKILL.md:56-61`).

| plot function | what it shows | file(s) written |
|---|---|---|
| `sc.pl.spatial(lrdata, color=<pair>)` | local bivariate cosine map, top-6 by Moran's R (`--n-top-plots 6`) | `plots/local_{COL4A1-CD44, COL4A2-CD44, APP-TNFRSF21, DLL3-NOTCH1, C3-C3AR1, TNC-SDC4}.png` (6) |
| `sc.pl.spatial(layer='pvals')` | per-cell permutation p-value map for each of the same top-6 | `plots/local_*_pvals.png` (6) |
| `sc.pl.spatial(layer='cats')` | high-high / low-low local category map for each of the same top-6 | `plots/local_*_cats.png` (6) |
| `sc.pl.spatial(adata, color=[lig, rec])` | constituent ligand + receptor gene maps beside each LR map | `plots/genes_*.png` (6) |
| `sc.pl.spatial`, ranked bar | the global Moran's R ranking | `plots/top_morans.png` |
| `li.ut.query_bandwidth` | neighbours-vs-radius calibration curve | `plots/bandwidth_query.png` (+ `plots/bandwidth_query.csv`) |
| `li.pl.connectivity` | the spatial kernel around one representative cell | `plots/connectivity_idx57404.png` |
| requested pairs — score, p-value, categories **and** gene maps | the two ALARMIST motif-1 arms, whatever their rank | `plots/requested/rank33_GRN-SORT1{,_pvals,_cats,_genes}.png`, `plots/requested/rank60_ANXA1-FPR1{,_pvals,_cats,_genes}.png` (8) |

> **Corrected 2026-08-04 (evening).** This section previously read "**11 PNGs**" and listed the
> `cats` maps, the top-6 p-value maps, the gene maps, the `query_bandwidth` curve and
> `li.pl.connectivity` as **not produced** for this tier. All five have since been produced —
> `run_liana.py` was re-run with `--plots-only` and the directory now holds 35 PNGs, verified by
> `find results/comparators/liana/GBM/cellchatdb2/plots -name '*.png' | wc -l`. The record of the
> gap is kept here for provenance; the gap itself is closed.

**Plots the bivariate workflow can produce and we did not** — recorded rather than produced:

| plot | tutorial cell | status / reason |
|---|---|---|
| `mask_negatives=True` re-plot | `bivariate.ipynb` 42-45 | ❌ not run; reproducible from the persisted score + category matrices with no re-fit. |
| cell2location compositions, decoupler TF activities (MuData route) | `bivariate.ipynb` 62-84 | **N/A** — single-cell Xenium has neither modality. *(Corrected 2026-08-04: this row previously also said "`decoupler` is not installed". It is — **decoupler 2.2.0**, and `annotate_factors.py` uses `dc.mt.mlm` / `dc.mt.ulm`. The remaining reason stands.)* |
| MISTy / multi-view extensions | `misty.ipynb`, `sma.ipynb` | ✅ **run 2026-08-04** — `misty.ipynb`'s LR-MISTy configuration on the full slide; see the *LIANA+ — LR-MISTy* section. `sma.ipynb` still not run. *(This row previously read "N/A here".)* |

### Multi-sample / differential mode

LIANA+ has extensive multi-sample machinery (`liana_c2c`, MOFA-based `mofatalk`,
`liana_pyCrossTalkeR`) — but those operate on **cell-type-level** LR results, not on the spatial
bivariate scores used here. For a spatial grade contrast the bivariate route offers no native
test; it would be hand-rolled from per-cell scores.

**None of it was used.** `grep -n by_sample scripts/comparators/liana/*.py` returns nothing;
neither `by_sample`, `dotplot_by_sample`, `lrs_to_views` nor `to_tensor_c2c` was ever called, and
the 13 `obs['tma_id']` cores never enter any test. This is an **open contract deviation** —
see *Contract deviations* in the second LIANA+ section.

### Gotchas

- **`query_bandwidth` returns a column spelled `bandwith`** (missing the second `d`). Indexing
  `'bandwidth'` raises `KeyError`.
- **`max_neighbours=100`** is a KNN ceiling on the graph — the same class of cap that binds
  destructively in SpatialDM. Verified non-binding here (**max 52** neighbours, 0.0000%).
- **`nz_prop=0.2` in the tutorial is 20% of Visium *spots*.** On our single-cell data the median
  gene is detected in **4.3%** of cells, so 0.2 is above the 90th percentile of genes and keeps
  only **13** pairs. See DEVIATIONS.md for the binomial spot→cell conversion.
- **Gene names containing `_` collide with the complex-subunit separator.** LIANA only *warns*
  (stLearn hard-errors); the 21 `Intergenic_Region_*` control probes are dropped anyway.
- **With `n_perms=100` the p-value grid is coarse** — the floor is 0, and **106/131** pairs report
  `morans_pvals = 0` (**113/131** at p < 0.05). p-values here separate "spatially structured"
  from "not"; they do not rank.
- Moran's R rewards spatial *structure*, not abundance. A ubiquitous, uniformly expressed pair
  scores near zero by design.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `results/comparators/liana/GBM/cellchatdb2/` | ✅ **2.0 min** | 100,197 cells **whole slide**; 3,218 → **131 pairs**; **113 with `morans_pvals` < 0.05** (106 at exactly 0); 10,119,897 connectivities; 35 plots |
| GBM | `cellchatdb2` (inflow) | `results/comparators/liana/GBM/cellchatdb2_inflow/` | ✅ **1.0 min** + 1.8 min downstream | 100,190 cells × **4,608** features; 80 figures. See the inflow section. |
| GBM | `cellchatdb2` (NMF) | `results/comparators/liana/GBM/{nmf_bivariate,nmf_inflow}/` | ✅ 0.9 min each | rank **6** / rank **7**, `k_range` **1..20**. See the inflow section. |
| GBM | `cellchatdb2` (local = Moran's R) | `results/comparators/liana/GBM/cellchatdb2_morans/` | ✅ **3.6 min** | same 131 pairs, `local_name='morans'` the only change; 35 plots. Global Moran's R **bit-identical** to the cosine run. See *Second local metric*. |
| GBM | **`default`** (bivariate) | `results/comparators/liana/GBM/default/` | ✅ **5.6 min** | LIANA **consensus**, 4,624 pairs → **388** tested; 35 plots. GRN^SORT1 **77/388**, ANXA1^FPR1 **187/388**. |
| GBM | **`default`** (inflow) | `results/comparators/liana/GBM/default_inflow/` | ✅ **1.3 min** + 3.4 min downstream | 100,190 cells × **9,448** features = 9 senders × **1,217** distinct LR pairs (650–1,173 per sender); **99.3555%** zeros; 80 figures. |
| GBM | **`default`** (NMF) | `results/comparators/liana/GBM/{nmf_bivariate_default,nmf_inflow_default}/` | ✅ refit 2026-08-04 | rank **6** / rank **5** on `k_range` **1..20**, matching the `cellchatdb2` tier. Inflow 9,448 → **6,178** after the ≥5/13 punch filter. *(The first pass reported 4 / 4; that was a `k_range` 1..10 artifact — `run_default_tier.sh` omitted `--k-max` and `run_nmf.py` defaults it to 11. Refitted on the matched window; the script now passes `--k-max 21` explicitly.)* |
| GBM | `cellchatdb2` (MOFA-Flex on inflow) | `results/comparators/liana/GBM/mofaflex_inflow/` | ✅ **76.0 min** (fit 70.5) | K=20, 17 active factors, 6 of 9 sender views survive QC. **The authors' QC deletes both arms of motif 1.** See the MOFA-Flex section. |
| GBM | `cellchatdb2` (MOFA-Flex, **reachability-normalised** QC) | `results/comparators/liana/GBM/mofaflex_inflow_reachnorm/` | ✅ **41.5 min** wall (fit 40.1), 2026-08-06 | K=20 requested, **19 active**; 4,608 → **779** features; **9 of 9** sender views kept, mGAM and Lymphoid restored. Grade test **0 / 19**. 32 PNGs. See the *reachability-normalised QC* section. |
| GBM | `cellchatdb2` (LRIC + cross-PCF) | `results/comparators/liana/GBM/lric_percore/` | ✅ **2.6 min** | **13 separate per-punch runs** + a pooled whole-slide control; 12 figures × png/pdf/svg. See the LRIC section. |
| GBM | `cellchatdb2` (LR-MISTy) | `results/comparators/liana/GBM/misty/linear_fullslide/` | ✅ **2.9 min** | whole slide, 382 receptor targets × 37 ligand predictors, 7 figures. RF secondary measured and **not** run (9.31 h). See the LR-MISTy section. |
| GBM | (NMF factor annotation) | `results/comparators/liana/GBM/factor_annotation/` | ✅ **0.22 min** | PROGENy + CellChatDB pathway annotation of the existing NMF factors; no re-fit. See *Factor annotation*. |
| LUAD | all | — | ❌ not run | deferred |

> **The `default` tier gap is closed (2026-08-04).** This table previously carried a
> `❌ still not run` row reading *"there is no `default` tree"*. There now is:
> `results/comparators/liana/GBM/` holds `default, default_inflow, nmf_bivariate_default,
> nmf_inflow_default`, produced by `scripts/comparators/liana/run_default_tier.sh` with every
> parameter except `--db` identical to the `cellchatdb2` tier. `SKILL.md:51-54` is satisfied.

**Whole-slide was verified safe rather than assumed**: LIANA's **28.2 µm** support radius against
a measured **222.9 µm minimum inter-core cell–cell distance** gives **zero** cross-core pairs
(still zero at 200 µm). A **7.9× margin** (222.9 / 28.2095). This is the only spatially-aware
method here run on the merged TMA by choice rather than by default.

### Top pairs and the requested LRIs

Top 8 by Moran's R (recomputed from `cellchatdb2/data/global_scores.csv`): **COL4A1^CD44
(0.1155), COL4A2^CD44 (0.1051), APP^TNFRSF21 (0.0900), DLL3^NOTCH1 (0.0895), C3^C3AR1 (0.0870),
TNC^SDC4 (0.0752), HLA-DQA1^CD4 (0.0739), IGF2^IGF2R (0.0694)** — ECM–CD44 adhesion dominates the
head of the list. (`LAMC1^CD44` is rank 11 at 0.0617, not top-8.)

| LRI | rank | Moran's R | p | ligand prop | receptor prop |
|---|---|---|---|---|---|
| GRN^SORT1 | **33 / 131** | 0.035702 | 0.0 | 0.109 | 0.188 |
| ANXA1^FPR1 | 60 / 131 | 0.013804 | 0.0 | 0.121 | 0.029 |

These ranks agree with `cellchatdb2/run_manifest.json` (`requested`) and with the filenames
`plots/requested/rank33_GRN-SORT1.png` / `rank60_ANXA1-FPR1.png`.

Both are **spatially significant** but with modest effect sizes. GRN→SORT1 sits at the **25th
percentile** (33 / 131) — better than chance but firmly mid-table, and well below COMMOT's rank
1–4. *(An earlier version of this document claimed rank 26 and "top 20%"; both were wrong.)*

### ⚠️ Merging the cores dilutes core-specific signal — measured

The smoke run on cores 13+14 alone ranked **`SLC17A7_GLS^GRIA1` (glutamatergic) first** by
Moran's R. On the full 13-core slide that pair drops out of the top entirely, replaced by
ECM–CD44 adhesion. Nothing is wrong: Moran's R is computed over all cores jointly, so a program
that is strong in two cores is averaged against eleven where it is not.

This is the concrete cost of the merged design — and it is worth weighing against the fact that
**SpatialDM and COMMOT, both run per-core, independently identified glutamatergic signalling as
top in core 13**. A whole-slide LIANA+ run cannot see that. If core-level heterogeneity is of
interest, LIANA+ should be re-run per core; the whole-slide result answers a different question
("what is spatially structured across this TMA overall").

### Methods paragraph

> For LIANA+ (v1.8.1), we identified spatial neighbours using the `spatial_neighbors` function
> and computed bivariate scores using the `bivariate` function, following the authors' tutorial.
> The Gaussian kernel bandwidth was set to 13.1454 µm with a cutoff of 0.1, corresponding to a
> support radius of 28.2 µm; this value was fixed by an area-preserving correspondence between
> the kernel's support disk and a 50 µm square patch (support radius = *s*/√π, bandwidth =
> support / √(−2 ln cutoff) = *s*/3.804), because the inflow branch of LIANA+ that applies to
> single-cell data specifies no numeric bandwidth rule. The resulting graph gave a median of 14
> neighbours per cell (maximum 52, below the `max_neighbours` ceiling of 100). Local scores were
> computed with the weighted cosine metric and global scores with bivariate Moran's R, using 100
> permutations and a minimum non-zero expression proportion (`nz_prop`) of 0.02; the latter
> replaces the tutorial's 0.2, which refers to multi-cell Visium spots and is not transferable to
> single-cell resolution. Ligand–receptor pairs were taken from CellChatDB v2, whose heteromeric
> complexes LIANA+ represents natively; 131 of 3,218 pairs passed the expression filter. All
> thirteen tissue microarray cores were analysed together, after verifying that the minimum
> distance between cells of different cores (222.9 µm) exceeds the kernel's support radius by
> 7.9-fold, so that no cross-core interaction can contribute. **No differential test between
> tumour grades was performed on the bivariate scores**, as LIANA+ provides no native spatial
> differential mode. (A grade contrast was run separately on the *inflow* NMF factors, aggregated
> to the thirteen cores; see the punch-level section.)

---

## LIANA+ — Inflow score and NMF communication programs

**This is the only comparator output that shares ALARMIST's shape** — a locations × factors
matrix (`NMF_W`) and a features × factors matrix (`NMF_H`) — so it is the only place where a
structural comparison against BPTF motifs is even possible. It therefore gets its own section.

### We used the wrong branch first, and it mattered

LIANA's README decision tree routes **Spatially-resolved → Single-cell → Interaction scoring →
Inflow Score**, while `li.mt.bivariate` (used above) sits under the **Spot-based** branch. Xenium
is single-cell. The package encodes the same judgement in its defaults: `_inflow.py` sets
`nz_prop = 0.001` against `bivariate`'s `0.05` — **the authors addressed single-cell sparsity by
writing a different method, not by asking users to lower a threshold.**

| | `li.mt.bivariate` | `li.mt.inflow` |
|---|---|---|
| branch | spot-based | **single-cell** |
| `nz_prop` | 0.05 (we used 0.02) | **0.001** |
| features | **131** | **4,608** (35.2×) |
| unique LR interactions | 131 | **633** (4.83×) |
| feature identity | LR pair | **sender cell type × LR pair** |
| non-negative | yes (0 of **13,125,807** = 100,197 × 131) | yes (0 of **461,675,520** = 100,190 × 4,608) |
| sparsity | 83.6214% zeros | **99.4538%** zeros |
| runtime | 2.0 min | 1.0 min |

**The inflow feature space is ragged, not a grid.** Splitting the 4,608 feature names on the
first `^` gives **9 senders × 633 distinct LR pairs**, with **250 to 602 LR pairs per sender**
(Lymphoid 250, non-mGAM 382, mGAM 508, Vascular 520, Glial-Neuronal 559, AC-like 590, OPC-like
596, NPC-like 601, MES-like 602). It is **not** a complete 512 × 9 grid — each sender keeps only
the LR pairs that clear `nz_prop` for *that* sender. *(Earlier statements of "535 unique
interactions / 4.1×" and "4,608 = 512 interactions × 9 senders" were both wrong; the 535 appears
to have been carried over from the CellChatDB Cell-Cell-Contact row count elsewhere in this
document.)*

Inflow carries **sender identity inside the feature** — `C_{j,s}` is a hard cell-type indicator
in `Inflow_{i,s,l,r} = (Σ_j W_ij L_{j,l} C_{j,s} / Σ_j W_ij) · R_{i,r}` — so a feature reads as
"cell type *s* sending ligand *l* to receptor *r* here". `bivariate` has no such axis.

Incidentally, the tutorial's SVG gene pre-filter was a **no-op** on our data: 5,097 → 5,097
genes, i.e. every gene on the Xenium panel is spatially autocorrelated at FDR<0.05 & I>0.01.
**The no-op is carried entirely by the `I > 0.01` half of the criterion**: min I = **0.159**, so
every gene clears it by an order of magnitude — while all 5,097 `pval_norm` values in
`data/gene_moranI.csv` are **exactly 0.0**, which at n = 100,190 is floating-point underflow of
the normal approximation, not evidence. Do not quote the FDR half as if it discriminated.

**The tutorial's *optional* SVI interaction filter is also a no-op — measured 2026-08-04, so the
deviation is now a number rather than an argument.** It was left off deliberately (pre-selecting the
feature space for spatial structure would confound the NMF comparison), and it is exposed as
`--svi-filter`. Applying it post-hoc: **4,608 → 4,608 interactions, none removed**, reproduced
twice, **19.4 s** at full slide over all **100,190** cells — no extrapolation was needed.
Same structure as the gene filter: **`I > 0.01` carries it** (min I **0.0578**,
`NPC-like^WNT4^FZD6_LRP6`; max **0.9962**, `Glial-Neuronal^CXCL2^CXCR1`; median **0.320**), and
**the FDR half is vacuous** — all 4,608 `pval_norm` and `pval_norm_fdr_bh` are **exactly 0.0**,
underflow at n = 100,190, exactly as for the genes. Table:
`cellchatdb2_inflow/data/interaction_moranI.csv`, with the measurement's own provenance record
(both replicates, timing, peak RSS) at `_benchmarks/svi_filter_measurement.json` — the third and
last entry in `_benchmarks/`, alongside `lric_punch3/` and `misty_punch4/`, and the only one this
document had not named.

### NMF: provenance of the composition

**NMF is demonstrated only in `bivariate.ipynb`; `inflow_score.ipynb` has no NMF section.** So
NMF-on-inflow is **our composition**, not an author-demonstrated path. Furthermore the decision
tree's own answer for unsupervised decomposition at single-cell resolution is
**"Communication Programs — Inflow + MOFA-Flex"** (`inflow_mofaflex.ipynb`) — a *different*
factorisation. ~~**MOFA-Flex on inflow remains the untested, author-sanctioned route.**~~
**Superseded 2026-08-04: MOFA-Flex on inflow has now been run** —
`results/comparators/liana/GBM/mofaflex_inflow/`, 76.0 min, K = 20, 17 active factors. See the
*MOFA-Flex on inflow* section, and note its headline: the tutorial's own QC removes every
`ANXA1^FPR1` feature and the whole mGAM view. The NMF runs are still **our** composition and that
does not change.

Both NMFs are kept: on `bivariate` as the tutorial-sanctioned composition, on `inflow` as the
resolution-appropriate input.

### Choosing `k_range` — the default is too narrow, but a wide range is not comparable either

`li.multi.nmf` defaults to `k_range = range(1, 11)`. Kneedle **normalises the curve over the
range it is given**, so a rank obtained on 1..11 and one obtained on 1..41 are not comparable
numbers. Three configurations were run; the **1..20** row is the one in force, because both
branches were fitted on the *same* window:

| `k_range` | bivariate rank | inflow rank | note |
|---|---|---|---|
| `range(1,11)` (default) | 3 | — | inflow's elbow lies outside this window |
| `range(1,41)` | 3 | 11 | historical; ranks not comparable to the row below |
| **`range(1,21)` — in force** | **6** | **7** | both branches, same window, neither at a boundary |

The runner warns when a rank lands on a `k_range` endpoint. NMF cost is not the constraint: fits
get *faster* with k (18.7 s at k=5, 10.8 s at k=40).

**Never quote "rank 6" / "rank 7" bare — the elbow is weakly supported.** `li.multi.nmf`'s elbow
metric is **MAE**, and on a matrix that is 83.6% (bivariate) / 99.45% (inflow) zeros the
zero-predictor is a formidable baseline. Recomputed directly from the stored `W`, `H` and the
input matrices:

| branch | rank | relative Frobenius error | fraction of SS captured | achieved MAE | zero-predictor MAE |
|---|---|---|---|---|---|
| bivariate | 6 | **0.7607** | **42%** | 0.080583 | 0.072522 |
| inflow | 7 | **0.7771** | **40%** | 0.015287 | 0.011817 |

**Neither fit beats the trivial all-zero predictor on the metric the elbow is computed with.**
That is a property of an L1 elbow on an extremely sparse non-negative matrix rather than a defect
in the run — but the rank should always be reported with `rel_frob` and the zero baseline
alongside it, not as a count of programs the tissue "has".

### How much does the program structure change? — a lot

Spatial correlation of the two `NMF_W` matrices over the 100,190 shared locations
(current configuration: bandwidth 13.1454, `k_range` 1..20, cross-punch filter applied;
full matrix in `nmf_factor_correlation.csv`):

| bivariate factor | best-matching inflow factor | r |
|---|---|---|
| F1 | inflow F7 | +0.230 |
| F2 | inflow F1 | **−0.041** |
| F3 | inflow F3 | +0.378 |
| F4 | inflow F3 | +0.232 |
| F5 | inflow F3 | +0.225 |
| F6 | inflow F7 | +0.430 |

**No bivariate program is well recovered by any inflow program** (best r = 0.43), and **5 of 7
inflow factors have max |r| < 0.3 to any bivariate factor**. Three separate bivariate factors
map onto the same inflow F3, so the mapping is not even one-to-one. The decompositions are
**not nested** — bivariate is not a coarse summary of inflow, it is a different partition.

The organising axis differs qualitatively, and this is the clearer result:

| | top features per factor |
|---|---|
| **bivariate — by LR family** | F1 JAG1^NOTCH1/2 · F2 immune/complement (C3^C3AR1, GAS6^AXL/MERTK) · F3 ECM (COL4A1/A2^CD44, LAMC1^CD44) · F4 TNC^SDC4, APP^SORL1 · F5 neurexin (NRXN2^CLSTN1/ADGRL1/NLGN3) · F6 DLL3/DLL1^NOTCH1/2 |
| **inflow — by sender cell type** | F1 Glial-Neuronal · F2 OPC-like · F3 MES-like · F4 NPC-like · F5 AC-like · F6 Glial-Neuronal · F7 NPC-like |

**Every one of inflow's 7 factors is anchored to a single sender cell type**; not one bivariate
factor is. Inflow also splits the NOTCH ligands by sender — OPC-like^DLL3^NOTCH1/2 (F2) vs
NPC-like^DLL3^NOTCH1/2 (F7) — where bivariate collapses all DLL3/DLL1→NOTCH signalling into a
single F6 regardless of who is sending. That is a genuine resolution gain attributable to the
finer feature space, not a re-labelling.

**The finding to carry forward: the recovered "communication programs" are strongly determined by
the upstream expression filter and feature construction, not only by the tissue.** A reviewer
shown only the bivariate NMF would conclude the GBM TMA contains six communication programs
organised by ligand family; shown only the inflow NMF, seven organised by sender cell type.
Both are defensible runs of the same package. Any claim about "communication programs" from this
class of method must state the upstream configuration.

### Spatial bandwidth — final parameter table

**An earlier run used bandwidth 18.75 µm; that value is VOID.** It was derived from the
*bivariate* tutorial's criterion — *"roughly includes 6 neighbours i.e. the first ring of
neighbours in the hexagonal grid"* — via `li.ut.query_bandwidth`. **"First ring" is a topological
property of Visium's hexagonal lattice. Xenium cells are irregularly packed and have no lattice,
so the criterion has no referent on this data and the number it produced is not defensible.**
The *inflow* tutorial, which is the branch that actually applies here, states no numeric rule at
all — only that bandwidth "should reflect the typical range of molecular signaling in the tissue",
traded off against resolution.

**Replacement: equal-area correspondence to a 50 µm square patch.**

| step | value |
|---|---|
| `k = sqrt(-2·ln cutoff)`, cutoff = 0.1 | 2.14597 |
| equal-area disk radius of an s×s patch, `s/√π` | **28.2095 µm** ← support radius |
| `bandwidth = R / k = s / 3.804` | **13.1454 µm** |

| parameter | bivariate | inflow |
|---|---|---|
| bandwidth | **13.1454 µm** | **13.1454 µm** |
| cutoff | 0.1 | 0.1 |
| **support radius** | **28.2 µm** | **28.2 µm** |
| median neighbours / cell | **14** | 14 |
| **max neighbours / cell** | **52** | 52 |
| `max_neighbours` cap | 100 (LIANA default) | 100 |
| **cap binding?** | **NO — 0.0000%** (52 < 100) | **NO — 0.0000%** |
| connectivity nonzeros | 10,119,897 | 10,119,897 |
| `nz_prop` | 0.02 | 0.001 (inflow default) |
| features tested | 131 LR pairs | **4,608** = 9 senders × 633 distinct LR pairs, **ragged** (250–602 per sender) |
| after cross-punch filter (≥5/13) | 131 (100% kept) | **2,704** (58.7% kept) |
| NMF `k_range` | **`range(1,21)`** | **`range(1,21)`** |
| **NMF rank** | **6** | **7** |
| rank at k_range boundary? | no | no |

The cap was left at LIANA's default because it is **already non-binding** (max 52 < 100); raising
it would have been a gratuitous deviation. Verified over all 100,197 cells, not a sample.

**Provenance note, stated plainly because a reviewer will probe it.** The 50 µm in this
derivation is ALARMIST's patch edge length. Only the *geometry* is borrowed — an area-preserving
square→disk conversion — and no ALARMIST output (motifs, loadings, LRI tables, factor counts)
enters at any point. It nevertheless means LIANA's spatial scale is **set by reference to the
method it is being compared against**, which is not the same as being set by LIANA's own defaults
or by the tissue. Anyone reporting these results must say so. The alternative — the authors'
qualitative "typical range of molecular signaling" — pins no number at all, which is why a
geometric convention was used instead.

#### The bandwidth was re-examined properly (2026-08-06) — and the number does not move

`scripts/comparators/liana/choose_bandwidth.py` →
`results/comparators/liana/GBM/bandwidth_choice/` (2 figures, 4 CSVs). It runs the tutorial's own
exploration step the way the tutorial runs it, then asks what would happen if the answer were
taken seriously. **Outcome: the bandwidth STAYS at σ = 13.1454 µm / support R = 28.2096 µm — a
user decision, not a derived one. CD-1 below therefore remains ❌ OPEN and must not be read as
resolved.** What the exercise did produce is the evidence that the exploration step *cannot*
settle it:

- **`li.ut.query_bandwidth` returns `ceil(MEDIAN) − 1`, not the mean**, despite the variable being
  named `avg_nn` (`liana/utils/query_bandwidth.py:71-72`:
  `avg_nn = np.ceil(np.median(num_neighbors))`, then `− 1`). This reconciles three numbers that
  otherwise disagree: our BallTree **mean** of 14.65 neighbours at R = 28.21, the **13** you read
  off the curve there (`query_bandwidth_tutorial_5_35.csv`: 13 at R = 28.08, 14 at R = 28.85), and
  the manifests' **median 14**. `ceil(14) − 1 = 13`.
- **⚠️ Its x-axis is a HARD query radius, not a σ.** `query_bandwidth` uses
  `BallTree.query_radius`, whereas `spatial_neighbors(bandwidth=…)` takes a Gaussian **σ**
  truncated at `cutoff`, with reach = σ × 2.145966. **A value read off that curve must be divided
  by 2.146 before being passed as `bandwidth=`**; passing it straight through inflates the
  neighbourhood **area by 4.6×**. This is an inconsistency in the tutorial itself, not in our
  runs. **Our figures do not inherit it** — see the retraction below.
- **✏️ Retracted 2026-08-07 — `bandwidth_query.png` was never on the wrong scale.** The
  2026-08-06 pass wrote here (and in two other places in this file) that
  `cellchatdb2_inflow/plots/global/bandwidth_query.png` draws its guide at the σ on a radius axis
  and should be regenerated. **False, and never checked against the code or the image.**
  `run_inflow_downstream.py:163-167` computes `R = a.bandwidth * np.sqrt(-2 * np.log(a.cutoff))`
  and passes it to `geom_vline(xintercept=float(R))`, i.e. the guide sits at the **support radius
  28.2096 µm** — the correct scale — labelled `support radius = 28.2 um (gaussian sigma =
  13.1454)`. The PNG on disk agrees: dashed line at ≈28 µm, crossing the curve at ≈14 neighbours.
  The only defect is cosmetic — the rotated annotation is anchored at `y = neighbours.max()` with
  `va="bottom"`, so it is clipped above the panel and reads as `su`. **Nothing needs regenerating
  for correctness.** Recorded rather than deleted because the wrong claim survived three documents.
- **The curve has no plateau, elbow or inflection anywhere between 5 and 120 µm**, so the
  exploration step does not select a value on this tissue. Mean neighbours per cell
  (`neighbours_vs_radius.csv`, interpolated at the quoted radii, all 100,190 cells):

  | R (µm) | 10 | 20 | 28.21 | 40 | 57.94 | 70 | 120 |
  |---|---|---|---|---|---|---|---|
  | mean neighbours | 1.75 | 7.49 | **14.65** | 28.88 | 59.16 | 85.21 | 240.06 |

- **Cell spacing gives a floor, not a choice.** Median nearest-neighbour distance **7.86 µm**
  (IQR 6.30–10.47, p95 19.97), so a strictly juxtacrine reach is ~8–10 µm (σ ≈ 3.7–4.9).
- **There is no characteristic signalling length scale in this tissue**, so the tutorial's
  biological criterion ("reflect the typical range of molecular signaling") has nothing to bind
  to. Pooling the LRIC radial profiles over **all 1,088 resolvable LR pairs × 13 punches =
  11,795 pair×punch observations** (`lric_percore/punches/punch_*/lric_agnostic_matrix.csv.gz`),
  the median *g(r)* by annulus is **1.459** (0–50 µm), 1.419, 1.403, 1.394, 1.401, 1.393, 1.381,
  **1.395** (200–225 µm) — a **4.4% decline across the entire range**. *(An earlier version of
  this argument rested on the **2** required LR pairs only; it is now 1,088 pairs and should be
  stated with that denominator.)*
- **LIANA's own `max_neighbours=100` default supplies the upper bound**: at R = 70 the mean is
  already 85.2. The defensible window is therefore R ∈ [20, 58], σ ∈ [9.3, 27].
- **The alternative not taken:** `inflow_mofaflex.ipynb`'s own value is **σ = 27 µm → R = 57.94 µm**
  (59.2 neighbours). Ours is **half that spatial scale, a quarter of the area**.
- **The bandwidth and the QC attrition are coupled**, which is new. A wider kernel would by itself
  have prevented much of the view loss described in the reachability section below: mGAM
  reachability rises **0.319 → 0.712** going from R = 28.21 to R = 57.94, and Vascular
  0.157 → 0.417 (`bandwidth_choice.json`).

### ⚠️ Contract deviations — OPEN, awaiting sign-off

These are departures from `.claude/skills/comparator-benchmark/SKILL.md`, **not** from the LIANA
tutorials. They were previously recorded only as tutorial deviations (bandwidth) or not at all
(differential mode). Recording them here as **open**, not settled. **Neither parameter is being
changed by this correction pass.**

| # | Invariant | What SKILL.md says | What was actually done | Status |
|---|---|---|---|---|
| **CD-1** | kernel scale | `SKILL.md:45-46`, verbatim: *"Keep each method's own neighborhood/kernel definition at its default. Do NOT harmonize spatial scale across methods, and do not match it to ALARMIST's patch size."* `SKILL.md:105` additionally lists an unpinned kernel scale under **STOP-and-ask**. | bandwidth **13.1454 µm** derived from ALARMIST's 50 µm patch edge (`run_liana.py:30-31`; derivation above). Every bivariate / inflow `run_manifest.json` carries it — **five of them across both tiers** (`cellchatdb2`, `cellchatdb2_inflow`, `cellchatdb2_morans`, `default`, `default_inflow`). *(Corrected 2026-08-04: this read "**eight** of them" and included the NMF trees. Parsing all 16 manifests: the four `nmf_*` ones carry **no** `bandwidth` key. Eight `run_manifest.json` do contain the literal string `13.1454`, but three of those hits — `lric_percore`, `misty/linear_fullslide`, `misty/rf_rate_probe` — are `params_provenance` strings stating the opposite, so the old count was a grep counting the negative statements.)* ⚠️ **Scope, 2026-08-04: CD-1 does NOT apply to the LRIC/cross-PCF or LR-MISTy branches.** LRIC consumes no connectivity graph (own `cKDTree`), and LR-MISTy was run at the *tutorial's* `bandwidth=200`. | ❌ **OPEN — violates the invariant.** Not a defensible "default", because LIANA pins none. **It is load-bearing:** the sensitivity table below measures this one parameter moving the inflow NMF rank **11 → 7**, so the headline "LIANA finds 7 programs" is a function of it. Needs sign-off, or a re-run at an independently justified scale. |
| **CD-2** | native multi-sample / differential mode | `SKILL.md:47-49`: GBM → split by `obs['grade']` (high vs low), with the **13 `obs['tma_id']` cores as the units**. | `li.mt.compute_global_specificity(groupby='grade')` only. Reading the installed `liana/method/sp/_compute_global_specificity.py`, that is a one-sided **per-group specificity** test that permutes labels across **cells** — not a contrast. `region_global_interactions.csv` is **9,216 rows = 4,608 × {high, low}** with **no contrast column**. `grep -n tma_id scripts/comparators/liana/run_inflow_downstream.py` returns **nothing** — the cores never enter the grade analysis at all (they are used only by `run_nmf.py`'s ≥5/13 presence filter). | ❌ **OPEN — requirement NOT satisfied, and in fact inapplicable.** LIANA has no native spatial differential mode, and `SKILL.md:49` says *"If the method has no multi-sample mode, say so — do not hand-roll one."* The 5,417 rows at p<0.05 are **cell-level permutation p-values (n = 100,190)**, pseudoreplicated by ~4 orders of magnitude relative to the 13 cores that are the real replicate unit. **A punch-level test has since been run as an *additional* analysis** — see below — and is null. |

**Wherever a grade-associated p-value from this method is quoted, carry the CD-2 caveat**: it is
a cell-level p-value, not a per-core one, and its magnitude is not interpretable. Treat the
*direction* of an effect as the result, not the *p*.

#### Punch-level grade test — the correct replicate unit, and it is null

`scripts/comparators/liana/analyse_existing.py` →
`results/comparators/liana/GBM/nmf_inflow/punch_level/` (2026-08-04, **0.12 wall-min**, no re-fit —
it reads `nmf_WH.npz` and `inflow_lrdata.h5ad` only). It aggregates the inflow `NMF_W` by
`obs['tma_id']` and runs a two-sided Mann-Whitney over the **7 high vs 6 low punches**, BH-corrected.

| test | signif. at BH q<0.05 | smallest raw p | smallest q |
|---|---|---|---|
| 7 inflow NMF factors vs grade | **0 / 7** | 0.013986 (Factor4, log2FC **+1.51**) | 0.0816 |
| 20 required-LR features vs grade | **0 / 20** | 0.013986 (`MES-like^GRN^SORT1`, log2FC **+2.62**) | 0.1399 |

**Nothing survives correction at the correct replicate unit**, against 5,417 "significant" rows at
the cell level — that gap is CD-2's pseudoreplication made concrete. But read the null carefully:
a two-sided rank test on 7 vs 6 units has a hard floor of **p = 0.001166**, so with BH over 7 (or
20) tests this design is underpowered and the result is *not* evidence of no grade effect. The
uncorrected direction is consistent across both arms of motif 1 (`GRN^SORT1` up in high grade in
the MES-like and AC-like senders), which is suggestive only.

**Outputs** — `nmf_inflow/punch_level/`, recorded 2026-08-04 because only two of its eleven files
had ever been cited:

| | |
|---|---|
| `data/` (4) | `punch_factor_tests.csv`, `punch_requiredLR_tests.csv` (the Mann-Whitney results above, with log2FC and BH q), and the underlying `punch_factor_means.csv`, `punch_requiredLR_means.csv` — 13 punches × factors / required-LR features, which is what the tests consume and what any re-analysis at a different replicate unit should start from |
| `plots/` (7) | `punch_factor_by_grade.png`, `punch_requiredLR_by_grade.png` (per-punch means split high vs low); `dotplot_top25.png`, `tileplot_top25.png`; `tileplot_GRN-SORT1.png`, `tileplot_ANXA1-FPR1.png`; `circle_plot_top50.png` |
| | `manifest.json` |

Those 7 PNGs are part of `nmf_inflow`'s 64-file total (2 in `plots/` + 55 in `plots_full/` + 7
here), which is how that count reconciles. They are **our** figures, not tutorial ones.

Note this test is **hand-rolled**, which `SKILL.md:49` discourages; it is recorded as an additional
analysis and does **not** convert CD-2 to satisfied.

**Updated 2026-08-04 — three more punch-level grade results, all null or untestable, all at the
7-high / 6-low replicate unit with the same `p = 0.0011655` floor:**

| analysis | grade result |
|---|---|
| MOFA-Flex primary (**17** active factors) | **0 significant**; smallest raw p **0.013986** (Factor 18) → q 0.217949 |
| MOFA-Flex sensitivity `nzf>0.001` (20 factors) | **0 significant**; smallest raw p **0.008159** (Factor 11) → q 0.163170 |
| LRIC / cross-PCF, cell-type-resolved | ⚠️ **not testable** — `min_cells` + `min_expressing` leave 8 of 13 punches informative and **7 of those 8 are high grade**, and the loss is *systematic* (the low-grade cores are the ones short of mGAM or of SORT1+/GRN+/FPR1+ MES-like cells). No test was run. The cell-type-**agnostic** fallback is null: `GRN^SORT1` p = 0.234, `ANXA1^FPR1` p = 0.788 |
| **MOFA-Flex reachability-normalised (19 active factors)** — *added 2026-08-06* | **0 significant**; smallest raw p **0.022145** (Factor 18) → q **0.332168** |

*(The first row read "(20 factors)". `mofaflex_inflow/data/factor_grade_punch_mannwhitney.csv` has
**17** rows, one per active factor — the test is run on the active set, not on the requested K.
The sensitivity row's 20 is right because all 20 are active there.)*

**Five** independent punch-level tests, five nulls (one untestable). That is consistent, but it is
**not** evidence of no grade effect — every one of them is bounded by the same 7-vs-6 floor.
CD-2 remains ❌.

### Sensitivity: the number of communication programs moves with the bandwidth

Both rows below were fitted at `k_range = range(1,41)`, so the bandwidth is the only thing that
differs between them — the rank change is attributable to bandwidth alone:

| bandwidth | support | inflow features | **inflow NMF rank** | bivariate NMF rank |
|---|---|---|---|---|
| 18.75 µm (void) | 40.2 µm | 4,815 | **11** | 3 |
| **13.1454 µm** (equal-area) | **28.2 µm** | **4,608** | **7** | **3** |

A 30% reduction in support radius took inflow's recovered program count from **11 to 7** — a
36% change — while bivariate stayed at 3 throughout. So the inflow decomposition's rank is
**sensitive to the neighbourhood scale**, and the bivariate one is not (it is pinned by its
131-feature ceiling, itself an artifact of `nz_prop`). Neither rank should be quoted as a
property of the tissue without stating the bandwidth that produced it.

*(The bivariate "3" in this table is the `range(1,41)` fit. On the in-force `range(1,21)` window
bivariate gives **6** — a `k_range` effect, not a bandwidth effect. Never compare a rank across
rows that used different windows.)*

### ⚠️ The inflow run was truncated — half of `inflow_score.ipynb` was never executed

`run_inflow.py` stopped immediately after `li.mt.inflow`, wrote its matrices, and exited. It
also created `plots/` and left it **empty** — it contained no plotting code at all, unlike its
bivariate sibling `run_liana.py`. Discovered 2026-08-04 because the empty directory was noticed,
not because anything failed: **every step below returned no error, they were simply never
called.** The tutorial continues well past `li.mt.inflow`:

| tutorial cell | call | was run? |
|---|---|---|
| 54 | `li.mt.compute_global_specificity(groupby='cell_type')` → `uns['global_interactions']` | ❌ → ✅ |
| 62 | same, `groupby=<region>` (here `grade`) → `uns['region_global_interactions']` | ❌ → ✅ **as a tutorial cell** — but this does **NOT** satisfy the benchmark's differential requirement. It is a per-group specificity test permuting labels across **cells**, not a high-vs-low contrast over the 13 cores. See **CD-2** above. |
| 67 | same, on `cell_type::region` composite labels | ❌ → ✅ |
| 75 | `li.ut.spatial_pair_proximity` | ❌ → ✅ |
| 83 | `li.mt.rank_aggregate` (spatially-constrained standard LR) | ❌ → ✅ |

Without cell 54 there is **no source→target significance for the inflow branch at all**, so the
branch had no per-LR result to report even if plotting code had existed. `scripts/comparators/
liana/run_inflow_downstream.py` now performs all five steps and every figure the tutorial draws
(**80 figures, 0 blank**, in `cellchatdb2_inflow/plots/{global,interactions}/` — that count and
that "0 blank" are scoped to *this* directory. **Four** blank PNGs have been found across the wider
LIANA tree, all the same plotnine return-type trap and all since regenerated, in **three** scripts:
two `connectivity.png` under `nmf_*/plots_full/global/` (`plot_liana_full.py`, described under
*Outputs*) and two `dotplot_focus_factors.png` under `mofaflex_inflow/plots/` and
`mofaflex_inflow/sensitivity_nzf0.001/plots/` (`run_mofaflex.py`, regenerated 2026-08-04 22:36/22:37;
see the MOFA-Flex section). Final sweep over the whole tree: **420 PNGs, 0 blank**.)
*(This passage read "the two blank PNGs found by audit were in `nmf_*/plots_full/global/`, written
by a different script" — correct when written, but it implied the audit was exhaustive at two
files, which it was not.)*

`compute_global_specificity` on 100,190 × 4,608 at `n_perms=1000` takes **8 s**, and the whole
downstream script **1.8 min** — the omission cost nothing in compute, it was simply not written.

### Requested LRIs — inflow corroborates ONE arm of the ALARMIST motif-1 loop, not both

`global_interactions`: 41,472 source×target×LR rows, 5,417 at p<0.05; 616 LRs significant
somewhere (of **633** distinct LRs — 17 are significant nowhere). p is permutation-based with
`n_perms=1000`, so **p = 0.000999 is the floor**, not a point estimate. **All of these p-values
are cell-level (n = 100,190), not per-core — see CD-2.**

| ALARMIST motif-1 arm | `lr_mean` | p | rank within that LR (of 81 pairs) |
|---|---|---|---|
| **ANXA1→FPR1, MES-like → mGAM** | 0.1263 | **0.000999** (floor) | **2 / 81** |
| **GRN→SORT1, mGAM → MES-like** | 0.0202 | **1.0** | 20 / 81 |

**The MES-like→mGAM arm is strongly corroborated; the mGAM→MES-like arm is not corroborated at
all** (p = 1.0). As whole LRs both rank mid-table among the 616 — GRN^SORT1 54th, ANXA1^FPR1
72nd in `lr_ranking_by_lr_mean.csv` — so neither is a headline interaction for the method; the
corroboration is specific to the *directed cell-type pair*, which is exactly what ALARMIST claims.

> **`lr_ranking_by_lr_mean.csv` is a *significance-conditioned maximum*, not a ranking by
> `lr_mean`.** `run_inflow_downstream.py` computes `gi[gi.pval < 0.05].groupby('lr')['lr_mean'].max()`
> — a maximum over the 81 source→target pairs, not a mean and not an unconditional ranking. Hence
> **616 rows against 633 distinct LRs**: 17 LRs are absent *by construction* because they are
> significant nowhere, and every rank in it is conditional on p < 0.05 somewhere. The six "top"
> LRs that were plotted were selected by this conditional maximum.
>
> *Since 2026-08-04 the file is self-documenting* — its columns are now
> `lr, max_lr_mean_over_signif_pairs, n_signif_pairs` rather than a bare `lr_mean` — and an
> unconditional companion **`lr_ranking_all_pairs.csv`** (633 rows, column `max_lr_mean_all_pairs`)
> sits beside it. Both put `GRN^SORT1` 54th and `ANXA1^FPR1` 72nd, so the ranks quoted above are
> unaffected by the conditioning.

**Caveat that limits how far this can be pushed.** For both LRs the ranking is dominated by
**self-self pairs** — GRN→SORT1's top four are Vascular→Vascular, mGAM→mGAM,
Glial-Neuronal→Glial-Neuronal, non-mGAM→non-mGAM; ANXA1→FPR1's top is mGAM→mGAM. The inflow
score multiplies a neighbourhood-averaged ligand term by the receiver's own receptor expression,
so a cell type that both expresses the ligand and clusters with itself scores highly without any
cross-type signalling. Read the off-diagonal cells, and treat the diagonal as a co-localisation
baseline rather than a result.

CellChat recovers **both** arms on this dataset (see its section) — so the disagreement is
between comparators, not between ALARMIST and the field.

### Relation to ALARMIST

ALARMIST used **K = 20** BPTF motifs on this dataset. Inflow+NMF gives **7** factors and
bivariate+NMF gives **6** (`k_range` 1..20, cross-punch filter applied), so neither is a
like-for-like reproduction, and the comparison is
structural rather than numerical. The meaningful similarity is that **inflow's factors are
cell-type-anchored signalling programs** — the same kind of object as an ALARMIST motif — whereas
bivariate's are ligand-family groupings. Any future motif-to-factor matching should therefore use
the **inflow** decomposition; matching against bivariate would compare objects of different kinds.

**Updated 2026-08-04.** There are now three more factor counts on the same input. **MOFA-Flex** on
inflow gives **17 active factors of 20 requested** (K is a user parameter there, not an elbow), and
the `default`-tier NMFs, refitted on the matched `k_range` 1..20, give **6 / 5** against
`cellchatdb2`'s **6 / 7**.

So the LR-database effect on the factor count is **much smaller than the first pass claimed**:
**none at all for bivariate** (6 → 6, despite 131 → 388 features) and **7 → 5 for inflow**. The
initially reported "6/7 → 4/4" was mostly the elbow window, not the database — see the confound
note. The spread that remains is still real (5, 6, 7, or 17 on the same tissue, depending on
resource, window and factorisation) but it is driven by the **factorisation and the window**, not
by the resource. **The factor count is not a property of the tissue and must never be compared to
ALARMIST's K = 20 as if it were.** The
*structural* comparison stands; the numerical one does not. And per *Factor annotation*, six of the
seven inflow NMF factors are ≥ 75% one sender, so even structurally they are closer to cell-type
identity than to ALARMIST motifs.

**Strengthened 2026-08-06.** The MOFA-Flex counts (17 / 19 / 20) are not merely "a user parameter"
— **`n_factors` is a ceiling and it binds in every fit**; see *`n_factors = 20` is a BINDING
CEILING* in the MOFA-Flex section. So there is no sense in which "MOFA-Flex found 19 programs" can
be set beside "ALARMIST K = 20": only one of the two numbers is even nominally an output.

#### Factor-vs-motif cosine — they agree on vocabulary, not on programs

`scripts/comparators/liana/cosine_factors_vs_motifs.py` →
`results/comparators/liana/GBM/vs_alarmist/` (2026-08-06), following the matching procedure in
`.claude/skills/alarmist`. It compares the **reachability-normalised** MOFA-Flex loadings against
`results/GBM/bptf/lri_motifs.csv`. Three confounds had to be handled explicitly, and each one
costs something:

1. **The feature spaces are not the same object.** ALARMIST is
   (sender, receiver, ligand, receptor, contact mode) — **25,271** features per motif; MOFA-Flex is
   (sender, ligand, receptor) — **779**. The only common space is the latter, so **ALARMIST must be
   collapsed by SUMMING over receiver and contact mode**: 25,271 rows → **4,756** keys, a median of
   **5** rows merged per key (mean 5.31, max 10). *Worked example — `mGAM|GRN|SORT1` on motif 1:*
   the MES-like receiver arm scores **3.091** of a summed **12.369**, i.e. **25%**. The
   biologically meaningful direction is a quarter of the number that enters the cosine, and
   mGAM→mGAM autocrine is summed in as if equivalent. **The comparison is therefore deliberately
   biased in LIANA's favour — it is an UPPER BOUND on agreement, not a neutral measurement.**
   ALARMIST mass retained after collapse + join: **90.9%** on raw `V`, **74.0%** on `V*`.
2. **DB export mismatch** (the subunit-ordering trap documented above). Keys are canonicalised by
   sorting subunits: raw overlap **713** keys → **742** after canonicalisation.
3. **⚠️ Signed vs non-negative — this one silently manufactured a result.** MOFA-Flex weights are
   signed; BPTF factors are not. Three modes are now implemented, `--sign-mode {poles,abs,signed}`:

| mode | what it does | verdict |
|---|---|---|
| `abs` | take `abs(weight)` | **WRONG, and it is what was originally used.** **57.2%** of weights are negative and the minor pole holds a median **38.5%** of a factor's mass (range 13.9–49.3%), so `abs` **merges two anti-correlated poles** and manufactures similarity. Kept only for provenance. |
| `poles` | split each factor into `max(w,0)` and `max(−w,0)` → 40 non-negative vectors | **default.** Each pole is a genuine non-negative feature set, directly comparable to a motif. |
| `signed` | signed weights as-is; cosine may be negative | plotted on a diverging `RdBu_r` centred at 0. |

Results, with ALARMIST scored on **`V*` = V/(mean_LR + 1)** (the prevalence-normalised column the
skill says to rank on), 742 shared keys, 200-permutation null preserving sparsity and magnitude:

| sign mode | max cosine | median best-match | motifs with a match > 0.5 | motifs beating the null at p < 0.05 | null median max-cosine |
|---|---|---|---|---|---|
| `abs` | 0.671 | 0.432 | **3 / 20** | 18 / 20 † | 0.241 † |
| **`poles`** | **0.743** | **0.517** | **13 / 20** | **20 / 20** | 0.221 |
| `signed` | +0.478 (max \|cos\| 0.643; matrix spans −0.643…+0.478) | 0.345 | **0 / 20** | 20 / 20 | 0.088 |

**⚠️ CORRECTION, 2026-08-06.** An earlier reading of these numbers — *"3 / 20; the methods are not
recovering the same programs"* — **was an artefact of the `abs()` handling**, not a finding.
Corrected reading: the two methods **agree substantially on LR VOCABULARY** (13/20 above 0.5, 20/20
above the null under `poles`) but **not on how that vocabulary groups into programs** (the
cell-space Spearman ceiling below is |ρ| ≈ 0.46, with 9 of 20 motifs collapsing onto one hub
factor) and **not on receiver or direction**, which are absent from LIANA by construction.

Three guards on that conclusion:

- **It is not a multiple-comparisons artefact.** Going from 20 to 40 candidate vectors *lowered*
  the permutation null — median max-cosine **0.241 → 0.221**, because pole vectors are sparser — so
  the `poles` numbers are not bought with extra candidates.
  **† The `abs` null is a RECOMPUTATION, not a file on disk.**
  `vs_alarmist/cosine_mofaflex_reachnorm_vs_alarmist.json` predates the permutation block and
  carries no `null_*` key (its filename also lacks the `_abs` suffix the current script emits).
  Re-running the script's own null procedure at its own seed returns **0.2413 / 18-of-20** for
  `abs`, and returns the stored `poles` values **0.2211 / 20-of-20 exactly**, which is the check
  that the recomputation is faithful. *(Corrected 2026-08-07: this bullet previously said the
  `abs` null "cannot be reproduced", which contradicted `DEVIATIONS.md` C2, where it had been
  reproduced. Only the `poles` (0.221) and `signed` (0.088) nulls are **recorded on disk**; the
  `abs` figures are reproducible on demand.)*
- **`aggfunc` is not load-bearing** (checked by recomputation): summing / max / mean when
  collapsing ALARMIST gives max |cos| **0.643 / 0.678 / 0.660** with **1 / 20** motifs above 0.5 in
  all three under `signed`, and 0.743 / 0.783 / 0.760 with 13 / 15 / 13 under `poles`. The choice
  moves the third decimal, not the conclusion.
- **On raw `V` instead of `V*`, 19 / 20 motifs clear 0.5 under `abs`.** That is a shared-prevalence
  artefact — both methods are dominated by the same ubiquitous adhesion pairs — not agreement.

Figures: `cosine_*` (sorted) and `clustermap_*` (clustered on both axes), `Reds` for the
non-negative modes and `RdBu_r` for `signed`, each as **PNG + PDF + SVG** and suffixed by sign
mode — 12 PNGs in `vs_alarmist/figures/`.

#### Per-cell Spearman — one hub factor, and it is an activity axis

`scripts/comparators/liana/compare_programs_to_alarmist.py` →
`vs_alarmist/comparison_summary.json` + `data/rho_*.csv`. Both methods emit a per-cell × per-program
matrix over the **same 100,190 cells**, so Spearman is directly computable — no collapsing, no key
repair. Against ALARMIST's 20 motifs:

| LIANA program set | programs | max \|ρ\| | motifs whose best match is the SAME factor | hub vs total ALARMIST loading | motif 1's best match |
|---|---|---|---|---|---|
| MOFA-Flex, tutorial QC | 17 | **0.518** | **13 / 20** → Factor 18 | **+0.506** (and +0.476 vs total inflow) | Factor 11, ρ **−0.222** |
| MOFA-Flex, sensitivity `nzf>0.001` | 20 | 0.458 | 6 / 20 → Factor 18 | +0.400 | Factor 7, ρ −0.266 |
| MOFA-Flex, reachability-normalised | 19 | 0.456 | **9 / 20** → Factor 18 | **+0.422** (+0.362 vs total inflow) | Factor 18, ρ **+0.214** |

**The hub is a general-activity axis, not a program.** It correlates +0.42 to +0.51 with a cell's
*total* ALARMIST loading and +0.30 to +0.48 with its *total* inflow — so most of what looks like
motif↔factor agreement in cell space is "this cell has a lot of signalling in it". By contrast
**ALARMIST motif 1 vs total inflow is only +0.146**, i.e. motif 1 is **specific**, not an activity
axis, which is precisely why it has no good match.

**Cell-type placement disagrees at the top.** ALARMIST motif 1 is highest in **mGAM (0.794)** and
lowest in Glial-Neuronal (0.307). Its best-matching factor (reach-norm Factor 18) peaks instead in
**NPC-like (0.654)** and puts mGAM **fifth** (0.504); the two agree only at the bottom, both
lowest in Glial-Neuronal (`data/celltype_placement_mofa-flex.csv`).

⚠️ **The NMF row of that script is UNVERIFIED — do not quote it.** For `kind == "csv_nmf"` the
script re-reads `nmf_inflow/data/NMF_W_factor_scores.csv` **without `index_col`**
(`compare_programs_to_alarmist.py:80-84`), so the numeric `cell_id` column survives
`select_dtypes("number")` and is scored as if it were a factor: the summary reports **8 factors for
a rank-7 NMF**, and `data/rho_nmf_8f.csv` has 8 columns. The header on disk is
`cell_id,Factor1,…,Factor7`. Open issue; the three MOFA-Flex rows are unaffected.

### Outputs

`results/comparators/liana/GBM/{nmf_bivariate,nmf_inflow}/`: `elbow.png`, `factor_maps.png`,
`NMF_W_factor_scores.csv`, `NMF_H_loadings.csv`, `nmf_WH.npz`, `nmf_errors.csv`,
`top10_loadings_per_factor.csv`, `punch_presence.csv`, plus `plots_full/{factors,global,
interactions}/` and `nmf_factor_correlation.csv` + `nmf_error_vs_k.png` at the GBM level.

#### ⚠️ `plots/` vs `plots_full/` — DISJOINT, and `plots_full` does **not** mean "the full set"

Four different scripts write a directory called `plots`, and the NMF dirs have two of them. They
must not be confused, and **neither is a superset of the other — not one filename appears in
both.**

| directory | written by | contents |
|---|---|---|
| `nmf_{bivariate,inflow}/plots/` | `run_nmf.py` | exactly **2** files: `elbow.png` (error vs k over 1..20 with the chosen rank marked) and `factor_maps.png` (one spatial `NMF_W` panel per factor) |
| `nmf_{bivariate,inflow}/plots_full/` | `plot_liana_full.py` | **30** real files (bivariate: 25 PNG + 4 CSV + `plot_manifest.json`) / **62** (inflow: 55 PNG + 6 CSV + manifest), under `factors/`, `global/`, `interactions/` (incl. `interactions/requested/`). *(Before the 2026-08-04 regeneration these were 22 / 31; the increase is the `requested/` dirs and the `feature_by_group` loop — see the defect write-up below.)* ✅ **`plot_manifest.json` now agrees with disk exactly — `n_files` 30 / 62 and a separate `n_png` 25 / 55.** *(It previously read 31 / 63, one higher in each tree, because the end-of-run glob counted the macOS `.DS_Store`; `plot_liana_full.py:479-480` now excludes dotfiles and counts PNGs separately. Verified with `find`.)* |

**`plots_full` means "the output of `plot_liana_full.py`" — a second, *additive* pass — not "the
complete set of plots."** The name is misleading and is also overloaded across this document
(stLearn uses `plots_full` to mean "went beyond the tutorial"). Concretely:

- **`elbow.png` and `factor_maps.png` exist only under `plots/`.** They have no counterpart
  anywhere under `plots_full/`. **Deleting `plots/` as "redundant" would destroy the entire
  rank-selection evidence.**
- Conversely every `plots_full/` file is absent from `plots/`.
- The only *conceptual* overlap: `plots/factor_maps.png` and the last panel of each
  `plots_full/factors/identity_Factor*.png` render the same `NMF_W` column spatially, magma at
  the 99th-percentile vmax, differing only in dot size.
- **Both directories visualise the same fit** — `plot_liana_full.py` reads the same
  `data/nmf_WH.npz` that `run_nmf.py` wrote, factor counts agree across `run_manifest.json`,
  `plot_manifest.json`, the elbow titles and the panel counts, and mtimes place `plots_full`
  strictly after each `nmf_WH.npz` with no refit between. **There is no stale-figure problem.**

#### Image outputs — `plot_liana_full.py` → `nmf_*/plots_full/`

| plot function | what it shows | file(s) written |
|---|---|---|
| `sc.pl.spatial` | all 100,197 cells coloured by `cell_type`, per-type counts in the legend | `global/cell_types.png` |
| `li.pl.connectivity` | the spatial kernel around one central cell | `global/connectivity.png` (both branches) — **real since the 2026-08-04 fix** (1080 × 1080, ~2,485 unique colours, mean 217.1; previously blank, see below) |
| — (CSV) | per-feature column sum of the score matrix, descending | `global/interaction_total_inflow.csv` |
| `li.mt.compute_global_specificity` → heatmap | source × target for the single top interaction | `global/{global_interactions.SUPERSEDED_nperms100_unseeded.csv, sender_receiver_heatmap.png, sender_receiver_<top interaction>.png}` — **inflow only** (`--global-specificity`). The CSV was **renamed** in the 2026-08-04 pass so it can no longer be mistaken for the canonical `cellchatdb2_inflow/data/global_interactions.csv` — see the duplicate-p-value note below. |
| — (CSV + heatmap) | factor × sender cell type, summed `H` loadings | `factors/factor_by_SENDER_celltype.csv`, `factors/heatmap_factor_by_SENDER.png` — **inflow only** since 2026-08-04; the bivariate branch has no sender axis and now writes neither (see below) |
| — (CSV + heatmap) | factor × receiver cell type, mean `NMF_W` | `factors/factor_by_RECEIVER_celltype.csv`, `factors/heatmap_factor_by_RECEIVER.png` |
| — (CSV + heatmap) | factor × LR pair | `factors/factor_by_LRpair.csv`, `factors/heatmap_factor_by_LRpair_top30.png` |
| factor identity panel | per factor: top-10 `NMF_H` features, sender composition (inflow only), receiver composition, spatial `NMF_W` map | `factors/identity_Factor{1..6}.png` (bivariate) / `{1..7}.png` (inflow) |
| lollipop grid | top-20 LR pairs per factor | `factors/top_lri_dot_by_factor.png` |
| receiver × sender heatmap per factor | rank-1 outer product of the receiver and sender profiles | `factors/celltype_communication_by_factor.png` — **inflow only** (correctly raises for bivariate, which has no sender axis) |
| 4-panel interaction map | inflow score map + ligand gene map + receptor gene map + sender-cell-type mask | `interactions/panel_<sender>-<L>-<R>.png` ×5 — **inflow only** (gated on `HAS_SENDER`) |
| `sc.pl.violin` | score by receiver cell type, top-5 features by total score | `interactions/violin_*.png` ×5 (both branches) |
| `li.pl.feature_by_group` | score × group dotplot | `interactions/feature_by_group_*.png` **×5**, one per top interaction. *(Before 2026-08-04 this was ×1 only — the call was a single `guard(...)`, not a loop, unlike the panel and violin loops.)* |
| **requested LRs (GRN→SORT1, ANXA1→FPR1)** | violin + `feature_by_group` (+ `panel_` and `sender_receiver_` on inflow) for every feature carrying a requested LR | ✅ `interactions/requested/` — **5 files** (bivariate: `violin_`/`feature_by_group_` × 2 LRs + `requested_lr_ranks.csv`), **25 files** (inflow: 6 plotted sender-resolved features × 4 figure types + the CSV). *(Before 2026-08-04 this was **absent entirely** from both trees — `plot_liana_full.py` had no required-LR mechanism, selection being purely `top_int = strength.index[:n_top_interactions]`, in violation of `SKILL.md:56-61`.)* |

#### Image outputs — `run_inflow_downstream.py` → `cellchatdb2_inflow/plots/` (80 figures)

`downstream_manifest.json`: `n_figures_saved: 80, n_skipped: 0`. Counts reconcile exactly:
16 in `global/` + 64 in `interactions/` (8 LRs × 8 figure types).

| plot function | what it shows | file(s) written |
|---|---|---|
| `sc.pl.spatial` | cell types, and grade, over the whole TMA | `global/spatial_{cell_type,grade}.png` |
| `li.ut.query_bandwidth` | neighbours-vs-radius calibration curve | `global/bandwidth_query.png` (+ `.csv` in `data/`) |
| `li.pl.connectivity` | kernel around two representative cells | `global/connectivity_idx{25047,50095}.png` |
| `li.pl.dotplot` | global source→target specificity dotplot | `global/dotplot_global.png` |
| `li.ut.spatial_pair_proximity` | source×target proximity heatmap | `global/pair_proximity.png` (+ `data/pair_proximity.csv`) |
| `li.mt.rank_aggregate` | spatially-constrained standard LR ranking | `global/rank_aggregate_dotplot.png`, `global/rank_aggregate_<LR>.png` ×8 |
| per-LR suite ×8 | `_inflow_by_sender`, `_genes`, `_sender_<type>`, `_violin`, `_feature_by_group`, `_source_target`, `_region`, `_niche_dotplot` | `interactions/<LR>_*.png` = **64** files for `CNTN2-CNTN2, NCAM1-NCAM1, C3-C3AR1, PECAM1-PECAM1, JAM3-JAM3, CD99-CD99, GRN-SORT1, ANXA1-FPR1` |

⚠️ **The two requested LRs are present but not segregated.** `GRN-SORT1_*` and `ANXA1-FPR1_*` sit
in the same 64-file directory as the method's own top-6, and `rank_aggregate_GRN-SORT1.png` sits
beside `rank_aggregate_PECAM1-PECAM1.png` — `run_inflow_downstream.py:189-193` appends the
required LRs to a single `top_lrs` list. `SKILL.md:56-61` asks for **separate output dirs** so a
requested pair is never mistaken for the method's own ranking (as `cellchatdb2/plots/requested/`
correctly does). Partial compliance; the ranking of each requested LR is recorded in this
document and in the run log, so nothing is actually misattributed today.

#### Three defects in `plot_liana_full.py` — all found by audit and since fixed

*(The third — the missing required-LR mechanism, in violation of `SKILL.md:56-61` — is described in
the image-outputs table above rather than repeated here; `liana/DEVIATIONS.md` heads the same
material "three defects" and writes all three up in full. This heading read "Two figure defects"
until 2026-08-04, so the two documents disagreed on the count for the same script.)*

Recorded because both silently produced *plausible-looking* output. Neither raised an exception,
and `guard()` logged "ok" in both cases.

1. **Blank `connectivity.png` in both `plots_full/global/`.** `li.pl.connectivity` draws with
   **plotnine**, but the saver fell back to `plt.gcf()` — an empty matplotlib canvas — so both
   files were written as a pure-white PNG (a single unique colour, mean pixel value 255.0). The
   correct idiom already existed 200 lines away in `run_inflow_downstream.py` (`return_fig=True`
   + `gg.save(...)`), whose `connectivity_idx*.png` outputs are real. **Fixed**:
   `plot_liana_full.py` now passes `return_fig=True` and routes plotnine objects through a
   dedicated `save_gg` helper, and `_blank_figure()` guards against a recurrence. **Verified after
   regeneration:** both files are now 1080 × 1080 with ~2,485 unique colours and mean pixel value
   217.1 — real figures, not white canvases. Impact was cosmetic — connectivity is illustrative
   and nothing downstream read it.
2. **Row-order corruption of the SENDER and LR-pair axes.**
   `li.ut.get_variable_loadings` **re-sorts its rows by |Factor1| descending**, so
   `data/NMF_H_loadings.csv` is *not* in `nmf_WH.npz['features']` order. The script built its
   annotation vectors (`sender_of`, `lr_of`) from the npz ordering and its value vectors from the
   CSV ordering, then masked one by the other — a pure permutation, so every figure rendered
   normally with the wrong labels. Affected:
   `plots_full/factors/{factor_by_SENDER_celltype.csv, factor_by_LRpair.csv,
   heatmap_factor_by_SENDER.png, heatmap_factor_by_LRpair_top30.png,
   celltype_communication_by_factor.png, identity_Factor*.png}`. Diagnostic symptom:
   `nmf_inflow/.../identity_Factor1.png` showed ten `Glial-Neuronal^…` bars under a suptitle
   reading "top sender: AC-like". **Fixed**: `plot_liana_full.py:116-125` now hard-fails on any
   feature missing from the CSV, reindexes `loadings` onto the npz feature order, and asserts
   alignment. **Not affected** (verified): `run_nmf.py`'s own `NMF_H_loadings.csv`,
   `top10_loadings_per_factor.csv`, `nmf_WH.npz`, `plots/elbow.png`, `plots/factor_maps.png`, the
   RECEIVER matrices (built from `W`, which is npz-ordered), the top-10 bar panels, and
   `top_lri_dot_by_factor.png` (label-indexed).

**The prose sender list in *How much does the program structure change?* above is correct** — it
was derived from `top10_loadings_per_factor.csv`, which was never affected. Where that list and
an older `plots_full` figure disagreed, the prose was right.

⚠️ **A degenerate bivariate SENDER table used to be written, and no longer is.** Before
2026-08-04, `nmf_bivariate/plots_full/factors/factor_by_SENDER_celltype.csv` was a single row
labelled `(none)` — bivariate features carry no sender — and `heatmap_factor_by_SENDER.png` was
drawn from it anyway. Misleading, never numerically wrong (a full-column sum is order-invariant).
**Both files are now absent from the bivariate tree**, verified by `ls
nmf_bivariate/plots_full/factors/`; the sender pair is written for `nmf_inflow` only, which is the
same guard `celltype_communication_by_factor` already had.

⚠️ **Two `global_interactions.csv` with different p-values — resolved by renaming.** The copy under
`nmf_inflow/plots_full/global/` was written by `plot_liana_full.py` at its default `--n-perms 100`,
**unseeded**; `cellchatdb2_inflow/data/` was written by `run_inflow_downstream.py` at
`n_perms=1000, seed=1337`. Both are 41,472 rows with identical `lr_mean`, but **393 rows flip
across p < 0.05** (5,290 vs 5,417). **`cellchatdb2_inflow/data/global_interactions.csv` is the
canonical copy.** Since 2026-08-04 the `plots_full` copy is named
`global_interactions.SUPERSEDED_nperms100_unseeded.csv`, so the two can no longer be confused by
filename — but note the manifest still does not record `n_perms`.

⚠️ **`plots_full/global/interaction_total_inflow.csv` has 4,608 rows** inside a tree whose NMF
used only the 2,704 post-punch-filter features (`plot_liana_full.py` loads the unfiltered h5ad).
Verified that no co-indexed array actually mixes the two spaces and all five plotted top
interactions are in the kept set — a latent hazard, not a current error. The bivariate copy of
this file is also mis-named: the column `total_inflow` holds summed local **cosine**.

⚠️ **`nmf_factor_correlation.csv` has no producing script in the repo — but its GBM-level neighbour
does.** `nmf_error_vs_k.png` is written by `scripts/comparators/liana/plot_nmf_errors.py`, which
reads only each branch's saved `data/nmf_errors.csv` and `run_manifest.json['nmf_rank']` and draws
three panels: error vs k per branch with the chosen rank dashed, plus both branches min–max scaled
onto a common k axis. Its docstring is explicit that it performs **no** stability selection,
resampling or rank test that LIANA does not itself perform. So the missing-producer warning below
applies to `nmf_factor_correlation.csv` **alone**. *(Added 2026-08-04: `plot_nmf_errors.py` was
named nowhere in either document, so this warning read as if both GBM-level artifacts were
orphaned.)*

`grep -rn nmf_factor_correlation scripts/` matches only this document. The **values are correct**
— the stored 6×7 matrix reproduces the cell-ID-aligned correlation and is inconsistent with
positional truncation; max |r| = 0.4299 matches the 0.43 quoted above. Only the code is missing,
so the table cannot currently be regenerated.

`results/comparators/liana/GBM/cellchatdb2_inflow/`:

| | |
|---|---|
| `data/` | `inflow_lrdata.h5ad` (now carries `global_interactions`, `region_global_interactions`, `niche_interactions` in `uns`), `inflow_scores.npz`, `global_interactions.csv`, `region_global_interactions.csv`, `niche_interactions.csv`, `pair_proximity.csv`, `rank_aggregate_liana_res.csv`, `lr_ranking_by_lr_mean.csv`, `lr_ranking_all_pairs.csv`, `bandwidth_query.csv`, `cell_meta.csv`, `feature_var.csv`, `gene_moranI.csv` |
| `plots/global/` (16) | `spatial_{cell_type,grade}`, `bandwidth_query`, `connectivity_idx*` ×2, `dotplot_global`, `pair_proximity`, `rank_aggregate_dotplot`, `rank_aggregate_<LR>` ×8 |
| `plots/interactions/` (64) | 8 LRs × 8 figures: `_inflow_by_sender`, `_genes`, `_sender_<type>`, `_violin`, `_feature_by_group`, `_source_target`, `_region`, `_niche_dotplot` |

**Two plotting traps, both hit here — and the return-type one has now bitten three separate
scripts.** `li.multi.nmf`'s internal `_plot_elbow`, `li.pl.dotplot`
and `li.pl.connectivity` return **plotnine** objects (save with `.save()`, and `plt.gcf()`
captures *nothing* after them), while `li.pl.feature_by_group` returns a **matplotlib `(fig, ax)`
tuple**. The three entry points disagree on return type, so a single save helper does not work.
`run_inflow_downstream.py` handled it correctly from the start; `plot_liana_full.py` did not (two
blank `connectivity.png`, fixed 2026-08-04); and `run_mofaflex.py` did not either (two blank
`dotplot_focus_factors.png`, found and fixed later the same day — see the MOFA-Flex section).
**Four blank PNGs across three scripts, all now regenerated; the tree sweeps clean at 420 PNGs,
0 blank.**
The elbow plots are regenerated inside `run_nmf.py` from the saved `nmf_errors.csv`, with the
same content as `_plot_elbow` (error vs k, dashed line at the chosen rank) — done in the runner
so it survives a re-run rather than as a one-off.

**Colour scaling on the inflow score maps.** The inflow matrix is **99.45% zeros** with a long
tail, so full-range scaling renders an all-black map. The maps use `vmin=0, vmax='p99.5',
sort_order=True`; percentile scaling is the tutorial's own idiom (`percentile_scaling=(1,97)` in
`li.pl.feature_by_group`).

### Reading traps in the persisted outputs

Small items, none of which changes a result, but each of which will mislead someone who reads the
file at face value.

- **`pair_proximity.csv`'s `interacting` column is `1` for all 81 source→target pairs** — it
  carries no information. That is expected at the deliberately coarse **100 µm** proximity
  bandwidth used for that call (distinct from the 13.1454 µm scoring kernel). **The usable output
  is the continuous `proximity` column**, which ranges 0.0435 – 0.9933.
- ~~**`bandwidth_query.png` draws its vertical line on the wrong scale.**~~ **Retracted
  2026-08-07 — this bullet was wrong on both counts.** The vline is *not* at the σ:
  `run_inflow_downstream.py:163-167` draws it at
  `R = bandwidth * sqrt(-2 ln cutoff)` = the **support radius 28.2096 µm**, correctly, and labels
  it as such; the PNG on disk shows the guide at ≈28 µm crossing the curve at ≈14 neighbours. And
  the quoted **13.46** at that radius reproduces from no file — the saved 2.95 µm grid
  (`cellchatdb2_inflow/data/bandwidth_query.csv`, `start=5, end=120, interval_n=40`) interpolates
  to **13.61**, the fine 0.77 µm sweep (`bandwidth_choice/data/query_bandwidth_tutorial_5_35.csv`)
  to **13.17**, reading exactly 13 at r = 28.077. The substantive reading is unchanged and simpler
  than the bullet claimed: **read the curve at 28.2 µm, where it gives 13, which is the realised
  median of 14 minus the self that `query_bandwidth` subtracts.** The figure's only real blemish is
  a clipped annotation. See `NOTES.md`, *`query_bandwidth`'s x-axis is a RADIUS*.
- **`inflow_lrdata.h5ad` still carries ALARMIST outputs**: `obs['motif']` (15 categories),
  `obs['patch_id']`, `uns['motif_colors']`, plus orphaned `obsp['connectivities']` /
  `obsp['distances']`. **Nothing in the LIANA pipeline ever read them** — `run_inflow_downstream.py`
  hard-refuses them as `--region-col`, and `downstream_manifest.json` records `region_col: grade`
  — so no result is contaminated. But it is a **leakage trap** for any future script that iterates
  `adata.obs` keys generically: drop `motif` and `patch_id` at load.
---

---

## LIANA+ — the `default` tier, and what the LR database does and does not confound

Driver: `scripts/comparators/liana/run_default_tier.sh` (5 steps, one log each under
`logs/comparators/liana-GBM-default-*.log`). Every parameter except `--db` is identical to the
`cellchatdb2` tier. This closes the `SKILL.md:51-54` two-tier requirement, previously recorded
in the runs table as *"still not run"*.

**The two resources are a real comparison, not a formality.** LIANA's own `consensus` resource has
**4,624** unique `(ligand, receptor)` pairs against CellChatDB v2's **3,218**, and they share only
**1,663 — 36.0% of consensus, 51.7% of CellChatDB**. Both required LRs (`GRN^SORT1`, `ANXA1^FPR1`)
are present in both.

| | `cellchatdb2` | `default` (consensus) |
|---|---|---|
| resource pairs | 3,218 | **4,624** |
| bivariate pairs after `nz_prop=0.02` | 131 | **388** |
| inflow features | 4,608 | **9,448** |
| inflow distinct LR pairs × senders | 633 × 9, ragged (250–602 per sender) | **1,217 × 9, ragged (650–1,173 per sender)** |
| inflow % zeros | 99.4538% | **99.3555%** |
| inflow features after ≥5/13 punch filter | 2,704 (58.7%) | **6,178 (65.4%)** |
| global-specificity rows (features × 9 targets) | 41,472 → 5,417 at p<0.05 | **85,032 → 12,902** |
| NMF rank (bivariate / inflow), `k_range` 1..20 both | 6 / 7 | **6 / 5** |
| wall time | 2.0 / 1.0 + 1.8 / 0.9 + 0.9 min | **5.6 / 1.3 + 3.4 / 3.3 + 2.2 min** |
| figures (bivariate / inflow) | 35 / 80 | **35 / 80** |
| NMF figures (bivariate / inflow) | 27 / 64 | **2 / 2** — see below |

⚠️ **`nmf_bivariate_default/` and `nmf_inflow_default/` hold 6 data files and exactly 2 PNGs each**
(`plots/{elbow,factor_maps}.png`, from `run_nmf.py`) and **no `plots_full/` tree at all** —
`plot_liana_full.py` was never run on the default tier. Against the `cellchatdb2` NMFs' 27 and 64
PNGs that is a 25 / 62 figure deficit, and it is a **deliberate scope limit, not a loss**: the
default tier exists to answer whether the resource confounds the statistics, and the data files
answer that. Data files in both dirs: `NMF_W_factor_scores.csv`, `NMF_H_loadings.csv`,
`nmf_WH.npz`, `nmf_errors.csv`, `top10_loadings_per_factor.csv`, `punch_presence.csv`, plus
`run_manifest.json`. *(Recorded 2026-08-04 — the asymmetry is on disk and was documented nowhere,
so a reader comparing the two tiers' figure trees would have assumed figures had been lost.)*

### ✅ The key result: the LR database is NOT a confounder for any per-interaction statistic

Both tiers were joined on their shared entries and compared value-by-value:

| comparison | shared | max &#124;difference&#124; |
|---|---|---|
| bivariate global Moran's R | **79 pairs** | **0.000e+00** (means also identical: 0.02256821 both) |
| bivariate `morans_pvals`, `mean`, `std` | 79 pairs | **0.0** each |
| inflow `lr_mean`, source×target×LR | **23,787 rows** | **0.0** |
| inflow `pval`, source×target×LR | 23,787 rows | **0 rows differ** |

This is not luck — it follows from LIANA scoring each pair independently, with no cross-pair
normalisation anywhere in `bivariate` or `inflow`. Swapping the resource therefore changes exactly
three things:

1. **which pairs are tested** — 388 vs 131, 9,448 vs 4,608 features;
2. **the multiple-testing / ranking denominator** — `GRN^SORT1` moves from rank **33 / 131** to
   **77 / 388** with an *unchanged* Moran's R of 0.035702, purely because the denominator grew
   (`ANXA1^FPR1`: 60/131 → 187/388, R unchanged at 0.013804);
3. **the feature space the factorisation sees** — and this is where it bites.

**Consequence for every rank quoted from this method: a LIANA rank is a statement about the
resource, not about the interaction.** The score is not.

### ✅ The factor count was confounded by `k_range` — caught, refitted, and the effect shrank

The first reading of "6 / 7 → 4 / 4" as a **database** effect did **not** survive checking the
manifests. The initial `nmf_*_default/run_manifest.json` recorded `"k_range": [1, 11, 1]` against
`[1, 21, 1]` for the `cellchatdb2` pair, because `run_default_tier.sh` called `run_nmf.py` without
`--k-max` and `run_nmf.py:31` defaults it to **11** — walking straight into the trap documented in
*Choosing `k_range`* above: **Kneedle normalises the elbow curve over the window it is given, so
ranks from different windows are not comparable numbers.**

**Refitted 2026-08-04 with `--k-max 21`** (and the driver now passes it explicitly, with a comment
saying it must never be omitted). On the matched 1..20 window:

| branch | `cellchatdb2` | `default` | features | rel-Frobenius (default) |
|---|---|---|---|---|
| bivariate | rank **6** | rank **6** | 131 → 388 | 0.7582, 42.5% SS |
| inflow | rank **7** | rank **5** | 2,704 → 6,178 | 0.8294, 31.2% SS |

So the LR database's effect on the factor count is **far smaller than first reported**: **none at
all for bivariate** — rank 6 either way, despite the feature count tripling — and **7 → 5** for
inflow. The database changes *which* features are factorised, not *how many* programs the elbow
finds. Both default fits still fail the zero-predictor check (MAE 0.085256 and 0.016617 vs
0.077797 and 0.012344), exactly as the `cellchatdb2` fits do.

What *is* safely established, and should sit next to any factor-count claim:

| upstream choice | measured effect on the inflow NMF rank | status |
|---|---|---|
| bandwidth 18.75 → 13.1454 µm (`k_range` 1..40 both) | **11 → 7** | clean — only bandwidth differs |
| `k_range` 1..10 vs 1..40 vs 1..20 (bivariate) | 3 / 3 / 6 | clean — only the window differs |
| LR database, `cellchatdb2` → `default` (`k_range` 1..20 both) | **7 → 5** (inflow); **6 → 6** (bivariate) | clean — refitted 2026-08-04 on the matched window |

So the rank is sensitive to the **bandwidth** and to the **elbow window**, both strongly. The
**database** effect, once measured cleanly, is the weakest of the three: it does not move the
bivariate rank at all (6 → 6 across a 3× change in feature count) and moves inflow only 7 → 5.
The earlier "7 → 4" was the window, not the resource.

### The zero-predictor warning fires on both default NMF runs too

Same pattern as the `cellchatdb2` fits, recomputed by `run_nmf.py` itself and stored in each
manifest:

| branch | rank | rel. Frobenius error | SS captured | achieved MAE | zero-predictor MAE | beats zero? |
|---|---|---|---|---|---|---|
| `nmf_bivariate_default` | 6 | **0.758200** | **42.5%** | 0.08525600 | 0.07779735 | **no** |
| `nmf_inflow_default` | 5 | **0.829400** | **31.2%** | 0.01661700 | 0.01234411 | **no** |

Both are *worse* than the `cellchatdb2` pair (0.7607 / 42% and 0.7771 / 40%). The caveat already
attached to "LIANA finds 6 / 7 programs" applies unchanged to the default tier's 6 / 5.

### Small provenance blemishes in the default tier — recorded, none affects a number

- `default/run_manifest.json` records `"dataset": "default"` (the dataset is GBM; the string is
  the output-dir name).
- All four default-tier manifests, and `default_inflow/downstream_manifest.json`, record
  `"tier": "cellchatdb2"` — the tier string is hardcoded in `run_nmf.py` /
  `run_inflow_downstream.py` and was not parameterised. **The `resource` fields are correct**
  (`"LIANA consensus"`, `resource_n_pairs: 4624`), so the runs are not misidentifiable; the `tier`
  key alone is wrong.
- `default_inflow/downstream_manifest.json` records
  `"db": "/Users/jiayifan/tansey_lab/alarmist/consensus"` — a path that does not exist. The literal
  `consensus` was resolved against the repo root by the manifest writer, not by the loader;
  `run_inflow_downstream.py:387-388` correctly branches on `a.db == "consensus"` and calls
  `li.rs.select_resource("consensus")`. Cosmetic, but do not read that key as a file.
- `default/run_manifest.json` has `resource_fingerprint: null` (the fingerprint block only runs for
  file-backed resources).

---

## LIANA+ — second local metric: `cosine` vs `morans`

`results/comparators/liana/GBM/cellchatdb2_morans/`, produced by
`run_liana.py --local-name morans` (**3.6 min**, 35 figures, same 131 pairs, same resource, same
bandwidth). `local_name` is the **only** thing that changed. Both metrics are in
`li.mt.bivariate.show_functions()`, and that table cites the SpatialDM paper for `morans`.

### The global statistic is unchanged; the local one is a different object entirely

`global_name` was already `morans` in the original run, so global Moran's R, its p-values, and the
whole pair ranking are **bit-identical** across the two runs (max |Δ| = 0.0 over all 131 pairs).
Only `lrdata.X` differs — and it differs qualitatively:

| | `cosine` | `morans` |
|---|---|---|
| range over 100,197 × 131 | **[0, 1.0000002]** | **[−14.084, 73.273]** |
| zeros | **83.6214%** | **0.0000%** |
| negatives | **0** | **4,378,791 (33.4%)** |

### Agreement between the two: global yes, local barely

Per-pair Pearson correlation between the two local score matrices, over all 100,197 cells:

| | r |
|---|---|
| median over 131 pairs | **0.195** |
| range | **0.095 – 0.337** |
| `GRN^SORT1` | 0.189 |
| `ANXA1^FPR1` | 0.253 |

**"Is this interaction spatially structured overall" is robust to the local metric; "which cells
carry it" is almost entirely determined by it.** Any figure, niche assignment or downstream
factorisation built on the *local* scores is a statement about the chosen metric as much as about
the tissue.

### This retroactively justifies `cosine` for the NMF branch — and corrects a filter claim

Moran's R is **not NMF-admissible**: 33.4% of its entries are negative and none are zero.
`li.multi.nmf` requires a non-negative matrix, so the NMF branch could not have been run on it.
That is a *post hoc* justification, not the reason `cosine` was originally chosen, and it is
recorded as such.

It also **corrects a latent assumption in the cross-punch reproducibility filter.** That filter
counts a feature as present in a punch if it is non-zero there, and `run_nmf.py` implements it as
`X > 0`. For `cosine` (0 negatives) `X > 0` and `X != 0` coincide; **for `morans` they do not**, and
a third of the matrix would be silently treated as absent. The filter is only correct for
non-negative local metrics — it happened to be applied to one.

### Cross-method: LIANA whole-slide Moran's R vs SpatialDM per-core z

All **131** LIANA pairs appear in SpatialDM's union of **1,662** interactions across the 13 cores,
so the comparison needs no key repair. Taking LIANA's whole-slide Moran's R against SpatialDM's
**median per-core z**:

| | value |
|---|---|
| Spearman | **0.616** |
| Pearson | **0.588** |
| `GRN_SORT1` | selected in **7 / 13** cores, median z **2.266** |
| `ANXA1_FPR1` | selected in **6 / 11** cores present, median z **2.390** |

⚠️ **The caveat is not optional.** This is agreement between two implementations **at their own
default spatial scales** — SpatialDM per-core at ~709 neighbours/cell, LIANA whole-slide at a
28.2 µm support with a median of 14 neighbours — because `SKILL.md:45-46` forbids harmonising them.
The residual disagreement therefore mixes *implementation* and *spatial scale*, and **cannot be
attributed to either.** An r of 0.62 across a ~50× difference in neighbourhood size is the
headline; do not read it as an implementation-agreement coefficient.

### Outputs

Recorded 2026-08-04 — this branch had no outputs inventory. `cellchatdb2_morans/` mirrors the
`cellchatdb2` tree exactly (same script, same filenames, same figure inventory; only `local_name`
differs).

| | |
|---|---|
| `data/` (6) | `global_scores.csv` (131 × 10 — **bit-identical to the cosine run**), `local_scores.npz`, `local_pvals.npz`, `local_categories.npz` (cells × 131, all three *different* from the cosine run), `cell_meta.csv`, `bandwidth_query.csv` |
| `plots/` (27) | `local_*{,_pvals,_cats}.png` and `genes_*.png` for the same top-6 by Moran's R, plus `top_morans.png`, `bandwidth_query.png`, `connectivity_idx57404.png` |
| `plots/requested/` (8) | `rank33_GRN-SORT1{,_pvals,_cats,_genes}.png`, `rank60_ANXA1-FPR1{,_pvals,_cats,_genes}.png` |
| | `run_manifest.json` |

⚠️ **The `requested/` filenames still carry rank 33 / 60** — global Moran's R sets the ranking and is
unchanged; only `lrdata.X` differs. Do not read those ranks as a `morans`-*local* result.

---

## LIANA+ — LRIC and cross-PCF (spatial co-occurrence branch)

`scripts/comparators/liana/run_lric.py` →
`results/comparators/liana/GBM/lric_percore/`. This is the decision tree's
**Single-cell → Spatial co-occurrence → LRIC** branch, and the one place in this benchmark where a
method resolves **direction** (sender → receiver) as an explicit argument rather than as a feature
name.

**Core algorithm.** `li.mt.cross_pcf` is a cross-type pair-correlation function: for an annulus
[r, r+w) it counts observed *A*-near-*B* pairs and divides by the count expected under complete
spatial randomness, giving g(r) with g = 1 meaning "no more co-located than chance".
`li.mt.lric` is the same construction restricted to the *ligand-expressing* members of *A* and the
*receptor-expressing* members of *B*. **The LRIC / cross-PCF ratio is therefore the whole point**:
it asks whether expression adds anything over cell-type co-location.

**Spatial model.** Annuli, no kernel and no connectivity graph: both functions build their own
`cKDTree`. `li.ut.spatial_neighbors` is deliberately **not** called, so the repo's 13.1454 µm
bandwidth does not enter this branch at all. Realised bins:
`[0,50) [50,75) [75,100) [100,125) [125,150) [150,175) [175,200) [200,225)` — 8 bins, the first
double-width because `extend_first_annulus` merges the contact band into it, and `max_radius=200`
is the **inner** edge of the last bin so the true reach is **225 µm**.

**Gates.** `min_cells=50` (a cell type is dropped from a punch entirely) and `min_expressing=20`
(an LR pair is NaN'd for that directed type pair). Both bite here — see *Support* below.

### Why per punch — measured, not asserted

Both functions normalise by density computed as `n_points / bounding-box area`
(`_LRIC.py::_corrected_areas`, `CrossPCF._compute_pair`). On a TMA the global bounding box is
mostly empty:

| | value |
|---|---|
| global bbox | **123,358,214.8 µm²** |
| Σ of the 13 punch bboxes | **52,094,949.7 µm²** |
| occupancy | **42.2306%** |
| ⇒ whole-slide density understated by | **2.3679×** |

A pooled whole-slide control was run anyway (`--whole-slide-check`,
`lric_percore/whole_slide_check/`), and the **measured** inflation exceeds that lower bound,
because each punch bbox itself still has empty corners:

| | median inflation, whole slide over the per-punch median |
|---|---|
| LRIC g(r) | **3.4284×** |
| cross-PCF g(r) | **2.9445×** |

**And it distorts the ratio, which is the quantity under test.** Whole-slide gives
`GRN^SORT1` mGAM→MES-like ratios of **1.212–1.284** per bin and `ANXA1^FPR1` mGAM→MES-like
**1.737–2.261**, against per-punch medians of **0.951–1.073** and **0.926–1.134**. A whole-slide run
would have manufactured a **false positive** — "expression adds ~2× over co-location" — for exactly
the claim being tested. This is recorded in `run_manifest.json` under `why_per_punch`.

**Note the contrast with bivariate/inflow, which *were* safe whole-slide.** Their 28.2 µm Gaussian
support against a 222.9 µm inter-punch floor gives zero cross-punch pairs — but **density
normalisation is a different exposure and is not protected by that argument.** "LIANA was fine
whole-slide" is a per-branch statement.

### Cost and support

13 punches at 816–26,456 cells, 3,418–5,092 genes, 299–1,087 LR pairs, 6–9 cell types kept at
`min_cells=50`. Per-punch LIANA call time **65.3 s** (cross-PCF 8.1 s, LRIC agnostic 24.3 s, LRIC
cell-type 32.9 s); pooled control a further **61.8 s**; whole script **155.7 s wall**, peak process
RSS **17.017 GB** (the pooled 100k-cell control dominates it — the largest single punch, punch 1 at
26,456 cells, peaked at 11.775 GB and took 21.67 s).

**Can the test even run?** mGAM clears `min_cells=50` in **11 / 13** punches (failures: punch 2 with
13 mGAM, punch 6 with 2 — both low grade). `min_expressing=20` removes more:

| target | punches with both directions available |
|---|---|
| `GRN^SORT1` | **8 / 13**, both directions |
| `ANXA1^FPR1` | **8 / 13** (MES-like→mGAM), **6 / 13** (mGAM→MES-like) |

The extra losses are punches 4, 12, 14 (MES-like SORT1+ = 7 / 15 / 7 and MES-like GRN+ = 7 / 5 / 3),
plus punches 5 and 11 for `ANXA1^FPR1` mGAM→MES-like (MES-like FPR1+ = 10 and 8). **Every gate is
attributed per punch** in `combined/target_availability.csv` and
`combined/target_expression_support.csv` — this is a design feature, not a footnote: the reason a
number is missing is on disk next to the numbers that are present.

### Result 1 — direction asymmetry: holds for GRN→SORT1, does not hold for ANXA1→FPR1

Punches, not cells, are the replicate unit throughout.

| LR | direction | n punches | median mean-g | IQR | range |
|---|---|---|---|---|---|
| `GRN^SORT1` | **mGAM→MES-like** | 8 | **1.566** | 1.424–1.668 | 1.329–2.961 |
| `GRN^SORT1` | MES-like→mGAM | 8 | 1.410 | 1.312–1.475 | 1.198–2.227 |
| `ANXA1^FPR1` | mGAM→MES-like | 6 | 1.524 | 1.459–1.626 | 1.355–1.978 |
| `ANXA1^FPR1` | MES-like→mGAM | 8 | 1.560 | 1.438–1.685 | 1.385–2.665 |

- `GRN^SORT1`: forward > reverse in **7 of 8** punches, paired Wilcoxon **p = 0.0156** (the minimum
  attainable at n = 8), median paired difference **+0.122**. Within each punch's own valid LR pairs
  the forward direction has median rank **41.5 of 70.5** (rank-pct 0.41) and the reverse **69 of
  87.5** (0.72) — the reverse direction is consistently the worse-ranked one.
- `ANXA1^FPR1`: paired on the 6 punches carrying both, median difference **+0.013 in favour of
  MES-like→mGAM**, **p = 0.844**, 4/6 punches. **The single-punch asymmetry does not replicate.**

### Result 2 — the LRIC / cross-PCF ratio is ~1 everywhere. This is the central finding, and it is negative

| LR | direction | per-punch median ratio | IQR | Wilcoxon vs 1 |
|---|---|---|---|---|
| `GRN^SORT1` | mGAM→MES-like | **1.043** | 0.965–1.077 | p = 0.383 |
| `GRN^SORT1` | MES-like→mGAM | **0.963** | 0.900–1.004 | p = 0.109 |
| `ANXA1^FPR1` | mGAM→MES-like | **1.040** | 0.974–1.097 | p = 0.438 |
| `ANXA1^FPR1` | MES-like→mGAM | **1.033** | 0.999–1.058 | p = 0.313 |

Per bin the median ratio stays in **0.93–1.13** everywhere except one dip (`GRN^SORT1`
MES-like→mGAM in [50,75) at **0.785**). **So across punches the spatial co-occurrence LIANA reports
for these pairs is fully explained by cell-type co-location plus the underlying point pattern;
ligand/receptor expression contributes nothing detectable.** Reported as a negative result.

**This also refutes the "the 28.2 µm kernel missed a longer-range interaction" hypothesis.** The
ratio is ~1 in *every* annulus out to **225 µm** — 8× the bivariate/inflow support radius. There is
no longer-range LR-specific structure for a wider kernel to have found.

### Result 3 — the grade contrast is not testable here, and the loss is systematic

Of the 8 punches supporting `GRN^SORT1`, **7 are high grade (1, 3, 5, 8, 10, 11, 13) and exactly 1
is low grade (punch 9)**. The loss is not random: the low-grade punches are precisely the ones where
mGAM falls below `min_cells` (2, 6) or MES-like SORT1+/GRN+/FPR1+/ANXA1+ falls below
`min_expressing` (4, 12, 14). **n = 7 vs 1 supports no test, and none was run.**

Fallback on the cell-type-**agnostic** target curves (13/13 punches for `GRN^SORT1`, 11/13 for
`ANXA1^FPR1`), Mann-Whitney high vs low:

| LR | n high / low | median high | median low | p |
|---|---|---|---|---|
| `GRN^SORT1` | 7 / 6 | 1.556 | 2.053 | **0.234** |
| `ANXA1^FPR1` | 7 / 4 | 1.580 | 1.775 | **0.788** |

No grade difference. **And do not read the direction**: the low-grade punches pushing the agnostic
median up (2 at 4.14, 12 at 2.43, 6 at 2.21, 4 at 2.08) are the small, sparse ones, so this is more
likely residual within-punch density structure than biology.

### Reproducibility, and one off-by-one in the earlier benchmark

Punch 3 was re-run inside the 13-punch set so the aggregate is homogeneous;
`_benchmarks/lric_punch3/` was not touched. The re-run is **bit-identical** to the benchmark: max
absolute difference **0.0** across `lric_celltype_summary` (58,744 rows × mean/max/peak-radius, NaN
pattern identical) and `cross_pcf_long` (448 rows).

⚠️ **One discrepancy, in the benchmark's derived CSV rather than in its data.**
`_benchmarks/lric_punch3/TARGETS_lric_vs_crosspcf.csv` lists `rank_by_mean = 42` for `GRN^SORT1`
mGAM→MES-like, but the underlying `lric_celltype_summary.csv.gz` ranks it **41 of 60** valid pairs —
verified with no ties, its neighbours being `CDH2^CDH2` at 1.606021 and `LAIR1^LILRB4` at 1.603856
around `GRN^SORT1` at 1.605142. The other three target ranks (13/60, 40/75, 75/75) match exactly.
**The g values are unaffected; the 42 is an off-by-one in that summary file only.**

Also note that punch 3's "75 / 75, dead last" for the reverse direction is the **extreme of the
distribution, not the typical case** (median reverse rank-pct 0.72), and that punch 3's
`ANXA1^FPR1` preference pointed the *opposite* way to the biology — both are reasons the
single-punch benchmark was repeated across all 13.

### Outputs

`combined/`: `targets_long.csv`, `per_punch_feature_summary.csv`, `aggregate_per_bin.csv`,
`aggregate_per_feature.csv`, `statistics.csv`, `target_availability.csv`,
`target_expression_support.csv`, `agnostic_per_punch_summary.csv`, `targets_agnostic_long.csv`,
`whole_slide_vs_perpunch.csv`, `punch_summary.csv`, `celltype_kept_skipped.csv`, `timings.csv`.
`punches/punch_<id>/` × 13: `cross_pcf_long.csv`, `lric_agnostic_matrix.csv.gz`,
`lric_agnostic_pair_scores.csv`, `lric_celltype.npz`, `lric_celltype_summary.csv.gz`,
`punch_report.json`. `whole_slide_check/` holds the pooled control (labelled a **control**, not a
result, in the code, the manifest and the figure legend). `figures/`: **12 figures × png + pdf +
svg = 36 files** — the only LIANA branch that exports all three formats.

---

## LIANA+ — LR-MISTy (multi-view modelling)

`scripts/comparators/liana/run_misty.py` →
`results/comparators/liana/GBM/misty/linear_fullslide/`. Tutorial:
`misty.ipynb`'s own LR-MISTy configuration. This replaces the earlier punch-4 benchmark, which had
two defects: it used the wrong model for the headline, and it was run on a **low-grade punch where
neither GRN nor ANXA1 entered the extra view at all** (see below).

**Core algorithm.** `lrMistyData` builds two views over the same cells: an **intra** view of
receptor expression (the targets) and an **extra** view of *spatially-weighted neighbourhood ligand*
expression (the predictors). Each target is regressed on the extra view under 10-fold CV; the
reported statistic is the CV R². It is a *predictive* model, not a test — there is no null and no
p-value.

**Configuration**, all tutorial or package values, recorded in the manifest's `params_provenance`:
`lrMistyData(bandwidth=200, set_diag=False, cutoff=0.01, nz_threshold=0.1, kernel='misty_rbf')`,
then `misty(bypass_intra=True, model=LinearModel, k_cv=10, seed=1337)`. **Bandwidth 200 is the
tutorial's, not this repo's 13.1454** — stated explicitly, because this is the one LIANA branch
here that does *not* carry CD-1.

### Views, cost, and the two structural facts that determine how to read the output

| | value |
|---|---|
| cells | 100,197 loaded → **100,190** after `filter_cells(min_genes=10)` (7 dropped) |
| genes | 5,119 → **5,097** (21 `Intergenic_Region_*` control probes dropped, then `filter_genes(min_cells=3)`) |
| **intra view** | 100,190 × **382** receptor targets |
| **extra view** | 100,190 × **37** ligand predictors |
| connectivity | 10,119,190 nnz, degree **median = max = 100**, min 16, **cap binds for 99.7% of cells** |
| fit | **2.69 min** = 0.423 s/target × 382; construct 8.6 s; **2.91 min** total; peak RSS **6.12 GB** |

That is ~7× faster than the punch-4 extrapolation predicted, and it reproduced bit-identically
across two runs.

⚠️ **Fact 1 — `bypass_intra=True` means `gain_R2` is not a gain.** `_Misty.py` sets
`intra_r2 = ... if not bypass_intra else 0`, so `intra_R2` is **0.0 for all 379 non-NaN rows** and
`gain_R2 == multi_R2` **exactly** — it is the CV R² of the ligand extra view *alone*, with nothing
subtracted. The question "does the spatial extra view add anything over an intrinsic model" **is
never computed by the tutorial configuration.** The only place that comparison exists is the
benchmark's non-tutorial `tutorial_linear_with_intra` variant in `_benchmarks/misty_punch4/`.

⚠️ **Fact 2 — with one view the contribution column is 1.0 by construction.** The meta-model
coefficient is normalised to itself, so `contributions.csv` reads `extra = 1.0` for every target and
`plots/contributions_top.png` is uninformative **by design, not broken**.

### Target metrics — the hypothesised structural null does NOT hold

| statistic over 382 targets | value |
|---|---|
| median `gain_R2` | **0.00414** |
| q25 / q75 / q90 / q95 / q99 | 0.00033 / 0.01539 / 0.03382 / 0.05109 / 0.11107 |
| max | **0.57675** (MET) |
| > 0.01 | **133 / 382 (34.8%)** |
| > 0.05 / > 0.1 | 20 / 6 |
| exactly 0 | **57** |
| NaN | **3** |

Top targets: MET 0.5767, CNTN2 0.1779, CADM4 0.1212, GRIK3 0.1155, ADCYAP1R1 0.1098, CD44 0.1061,
CALCRL 0.0914, APLNR 0.0856, PDGFRA 0.0793, NOTCH1 0.0790.

The **57 exact zeros** are `clip(min=0)` applied to a *negative* mean CV R² — the extra view predicts
worse than the intercept — and they are the ultra-sparse receptors (SSTR3 at 0.22% non-zero, SELE
0.09%, SSTR4, HCRTR1, TREML2, CEACAM8, MC3R/MC4R…). The **3 NaNs** are heteromeric complexes
(`IL15RA_IL2RB_IL2RG`, `GABRA6_GABRB3_GABRG2`, `GABRA6_GABRB3_GABRD`): a complex is the elementwise
**min** over subunits, no cell co-expresses all three, the variance is 0, and liana's own log says
*"Variance of '<x>' is 0.0, metrics set to NaN"*.

### But the signal is spatial context, not ligand–receptor mechanism

Three concrete pieces of evidence, and the third is the important one:

1. **MET, the best target at R² = 0.577, is predicted without its own ligand.** Its top predictors
   are THY1 (t = −38.1), TNR (+32.0), TNC (+28.8), NTN1 (−27.0), CD99 (−25.0). **HGF is absent from
   the model entirely** — it is non-zero in only 3.77% of cells and fails `nz_threshold=0.1`. The
   model reconstructs MET from a mesenchymal/glial neighbourhood signature.
2. **Nothing is subtracted** (Fact 1 above), so `gain_R2` cannot separate "spatial ligand
   information" from "anything correlated with position".
3. **The extra view is 37 ligands wide out of 3,218 CellChatDB pairs**, and the 100-NN cap binds for
   99.7% of cells. What is being regressed is a heavily smoothed 100-cell neighbourhood average of
   37 abundant genes — a **niche-composition proxy**.

**Net: LR-MISTy on this slide is a spatial-context predictor of receptor abundance, not an LRI
detector, and its `gain_R2` ranking must not be read as an LRI ranking.**

### The requested LRIs — a genuine partial positive, and why the punch-4 negative was uninformative

`nz_fraction` is the gate that decides whether a ligand is even in the model:

| gene | full slide | punch 4 (the old benchmark) |
|---|---|---|
| GRN | **10.88% ✅** | 7.27% ❌ |
| ANXA1 | **12.06% ✅** | 2.52% ❌ |
| SORT1 | 18.75% | 24.82% |
| FPR1 | 2.94% | 4.36% |

**On punch 4 neither ligand cleared `nz_threshold=0.1`, so neither was in the extra view and the
model could not have found them.** The punch-4 negative was uninformative, exactly as suspected.
On the full slide both clear it, so this is a fair test:

Predictor ranks below are by **`|importance|`** — MISTy importances are signed, and this table asks
"how much does this predictor matter, in either direction". *(Convention stated explicitly
2026-08-07: the report README ranks the same rows by **signed** importance and therefore gives
different numbers — GRN 2nd and ANXA1 15th for FPR1, GRN 6th for SORT1. Neither is a correction of
the other; both are recomputed from `misty/linear_fullslide/data/interactions.csv`, `view ==
'extra'`, 37 predictors.)*

| target | `gain_R2` | rank | GRN as predictor | ANXA1 as predictor |
|---|---|---|---|---|
| **SORT1** | 0.06210 | **17 / 382** (top 4.5%) | **t = +5.46, rank 9 / 37** by \|imp\| (6 / 37 signed) | t = **−3.50**, rank 11 / 37 by \|imp\| (33 / 37 signed) |
| **FPR1** | 0.02615 | **54 / 382** (top 14%) | **t = +7.08, rank 3 / 37** by \|imp\| (2 / 37 signed) | t = **+1.34**, rank 22 / 37 by \|imp\| (15 / 37 signed) |

Top-5 predictors of SORT1: SEMA4D +8.8, NTN4 −6.7, SEMA4C +6.5, PRNP +6.2, DHCR24 +6.2.
Of FPR1: C3 +14.1, APP −8.4, **GRN +7.1**, TGFB2 +6.5, DHCR24 +4.1.

**GRN→SORT1 and GRN→FPR1 both come out with strong positive importance and both targets sit in the
upper half of the ranking. ANXA1→FPR1 does not replicate (t = +1.34, unranked), and ANXA1 enters
SORT1 with the wrong sign.** Read against the caveat above: this is evidence that neighbourhood GRN
tracks these receptors, not that a receptor–ligand mechanism was detected.

### RandomForest secondary — measured, not run, with the number that forced it

`RandomForestModel` is the tutorial's other model. A 4-target probe on the **full slide**
(`misty/rf_rate_probe/`) took **5.85 min = 87.7 s/target**, so 382 targets is **9.31 h** — ~9× the
~1 h budget, and **2.8× worse than the punch-4 extrapolation's 3.3 h central estimate**, so the
punch-4 scaling law under-predicts. The probe agrees with the linear model on ordering: MET 0.5950,
CD44 0.1169, SORT1 0.0444, FPR1 0.0190.

### Multi-sample / differential mode

**None was run, deliberately** — recorded 2026-08-04 because this is the one LIANA branch whose
silence on grade could otherwise be mistaken for an oversight. Every other branch states its grade
position explicitly (inflow via CD-2, NMF and MOFA-Flex via punch-level Mann-Whitney, LRIC via "not
testable, 7 of 8 informative punches are high grade"). LR-MISTy returns a cross-validated R² **per
target**, not a per-cell or per-core score, so there is no quantity to aggregate to the 13
`obs['tma_id']` cores and no native contrast in the package. Verified on disk: the branch writes no
grade or per-punch file at all, and `run_misty.py` reads `tma_id` only to report per-group cell
counts and non-zero fractions. Fitting per core would give 13 incomparable models — each with its
own `nz_threshold`-determined predictor set — which is why the punch-4 benchmark was abandoned for a
single full-slide fit. **So this branch contributes no grade result by design.**

### Deviations from `misty.ipynb`

| # | Item | Tutorial | Ours | Number that forced it |
|---|---|---|---|---|
| **DEV-1** | HVG pre-step `adata[:, hvg]` | applied | **omitted** | Measured with `--hvg 1 --n-top-genes 1521 --construct-only 1`: full slide intra **100,190 × 82** / extra **100,190 × 13**, against 382 / 37 without it — **78% of receptor targets and 65% of ligand predictors discarded**. On a single 5,363-cell punch it collapses the extra view to **one** predictor (intra 5363 × 54, extra 5363 × 1), i.e. vacuous. The tutorial frames HVG as *"for the sake of computational speed"* on a genome-wide Visium slide; this is a 5,119-gene targeted panel and the full-panel fit takes 2.7 min, so the speed argument does not apply. Kept as `misty/hvg_probe_construct_only/`. |
| DEV-2 | *(reading note, not a departure)* | `bypass_intra=True` | same | see Fact 1 above |
| DEV-3 | bandwidth | 200 | **200 — the tutorial's** | this branch deliberately does **not** inherit the repo's 13.1454 µm, so **CD-1 does not apply to it** |
| DEV-4 | `max_neighbours` | not exposed by `lrMistyData` | package default 100 | cap binds for **99.7%** of cells, so the 200 µm bandwidth is largely inert — the extra view is a 100-NN neighbourhood, not a 200 µm kernel. Tutorial behaviour, not a change. |
| DEV-5 | preprocessing | tutorial QC | `layers['counts']` → drop 21 `_`-genes → `filter_cells(min_genes=10)` + `filter_genes(min_cells=3)` → `normalize_total(1e4)` + `log1p` | matches `run_inflow.py`. ⚠️ `run_liana.py`'s bivariate branch **skips the QC**, which is why that branch has 100,197 cells against this run's 100,190 — **cross-branch joins must be on cell ID, not position.** |
| DEV-6 | resource | `consensus` | CellChatDB v2 (3,218 pairs) | consistency with the rest of `results/comparators/liana/GBM/` |
| DEV-7 | RandomForest secondary | run | **not run** | 87.7 s/target measured → 9.31 h; see above |

### Outputs

`linear_fullslide/data/`: `target_metrics.csv`, `target_metrics_ranked.csv`, `interactions.csv`
(14,014 rows), `interactions_SORT1.csv`, `interactions_FPR1.csv`, `contributions.csv`,
`view_features.csv`, `nz_fractions.csv`; plus `run_manifest.json` and `run.log`.
`linear_fullslide/plots/`: **7 PNGs** — `target_metrics_{gain_R2,multi_R2}_top`, `contributions_top`
(uninformative by construction, see Fact 2), `interactions_extra_top200`, `gain_R2_hist`,
`importances_{SORT1,FPR1}`. All three previously-unexercised `li.pl.*` MISTy plotting functions
(`target_metrics`, `contributions`, `interactions`) ran clean against liana 1.8.1.
⚠️ **PNG only — no PDF/SVG**, unlike the LRIC branch.

Two benign runtime messages, neither an error: `RuntimeWarning: invalid value encountered in divide`
at `_Misty.py:354` (that is `coefs / coefs.sum()` when the single view's mean coefficient clips to
0 — it is exactly the 57 zero-gain targets), and a mudata `pull_on_update` FutureWarning.

---

## LIANA+ — MOFA-Flex on inflow (the authors' prescribed unsupervised route)

`scripts/comparators/liana/run_mofaflex.py` →
`results/comparators/liana/GBM/mofaflex_inflow/`, following `inflow_mofaflex.ipynb` cell by cell.
**This closes the gap this document and `DEVIATIONS.md` both called "the untested,
author-sanctioned alternative"** — the decision tree's answer for *Single-cell → Unsupervised* is
**Communication Programs = Inflow + MOFA-Flex**, not NMF.

⚠️ **CD-1 scope — it DOES apply to this branch**, and that was never stated. MOFA-Flex re-specifies
no spatial model of its own here: it factorises `cellchatdb2_inflow/data/inflow_lrdata.h5ad`
(confirmed from the manifest's `argv`), so all spatial weighting already happened inside
`li.mt.inflow` at the ALARMIST-derived **13.1454 µm** bandwidth. Unlike LRIC (own `cKDTree`) and
LR-MISTy (tutorial `bandwidth=200`), which the CD-1 row explicitly exempts, this branch inherits the
contract deviation — record it wherever a MOFA-Flex factor count is quoted. *(Recorded 2026-08-04;
CD-1's scope note named only bivariate / inflow / NMF, LRIC and LR-MISTy, leaving MOFA-Flex
unclassified.)*

### ⚠️ The single most reportable result: the authors' own QC deletes both arms of ALARMIST motif 1

Two consecutive tutorial cells do it, and neither errors:

| tutorial cell | filter | effect on this dataset |
|---|---|---|
| 19 | `nonzero_fraction > 0.01` | **4,608 → 447 features.** Removes **every** `ANXA1^FPR1` feature: its maximum non-zero fraction across all 9 senders is **0.009422** (MES-like); mGAM is 0.004911. All 9 fail. |
| 23 | `lrdata_to_mudata(min_features=25)` | drops the **entire mGAM view**, which retained exactly **24** features — **one short**. (`non-mGAM` also drops, at 4.) So `mGAM^GRN^SORT1` is gone even though its non-zero fraction, 0.031770, comfortably cleared cell 19. |

**Run exactly as its authors demonstrate, LIANA+ structurally cannot see either arm of motif 1 on
this dataset.** This is a finding about the comparator, not about the biology, and it is the reason
a sensitivity fit was added.

**Updated 2026-08-06 — and the cell-19 cut is worse than "a threshold that happens to be strict".**
It is an **abundance filter in disguise**: an inflow feature can only be non-zero where its sender
is reachable, so a single global threshold imposes a different effective stringency per sender
(surviving-feature count vs sender abundance, **Spearman ρ = 0.917**). A reachability-normalised
re-fit that keeps **all nine** views is documented in its own section below; it restores mGAM
(24 → 66 features) and Lymphoid (0 → 39), **and the two arms still do not co-load**.

### Primary fit

CPU, K = 20, batch 2048, lr 0.005, patience 50, seed 0; early-stopped at **632 / 1000** epochs.
Fit **70.5 min**, total wall **76.0 min**, peak RSS **4.11 GB**. A 20-epoch determinism probe is
bit-identical (`max_abs_weight_diff: 0.0`). ⚠️ **These numbers now come from
`logs/mofaflex_primary.log`, not from the manifest** — the 2026-08-04 22:36 replot pass overwrote
`run_manifest.json`, which today records `fit_seconds: null`, `determinism_probe: null` and a
`wall_seconds` / `peak_rss_gb` describing the replot. See the provenance note under *Outputs*.

| | value |
|---|---|
| views surviving | **6 of 9 senders**, **419 features** (NPC-like 98, MES-like 97, AC-like 92, OPC-like 78, Glial-Neuronal 28, Vascular 26) |
| dropped views | mGAM (24 features), non-mGAM (4); Lymphoid absent already |
| active factors (R² ≥ 0.02 in ≥1 view) | **17 / 20** |
| max single-view R² | **0.1681** |
| max &#124;r&#124; between factor scores | **0.243** — **zero** pairs exceed the tutorial's 0.6 redundancy flag |

R² per sender view: Glial-Neuronal **0.602**, Vascular 0.399, MES-like 0.398, NPC-like 0.394,
OPC-like 0.338, AC-like 0.302.

### ⚠️ `n_factors = 20` is a BINDING CEILING, not a selection (2026-08-06)

**MOFA-Flex does not choose K.** `n_factors` is a hard ceiling: the model fits exactly that many
and inactive ones are pruned *afterwards* by the 2%-R²-in-one-view floor. The 20 was taken from
`inflow_mofaflex.ipynb` cell 25, which uses `n_factors=20`. **It coincidentally equals ALARMIST's
K = 20; that coincidence is not the justification and must never be presented as one.**

The ceiling **binds** in all three fits, and there is no taper at the bottom of the list:

| fit | active / requested | share of total R² held by the 3 weakest active factors | weakest factor's best single-view R² (floor 0.02) |
|---|---|---|---|
| primary (tutorial QC) | **17 / 20** | 6.7% | — |
| sensitivity `nzf>0.001` | **20 / 20** | 8.3% | — |
| reachability-normalised | **19 / 20** | 7.7% | **0.0287** (Factor 6) — clears the floor by 44% |

So **17 / 19 / 20 are artefacts of the argument, not discovered numbers**, and none of them may be
compared to ALARMIST's K = 20 as though both were fitted. This is the same failure class as the
`k_range` confound already recorded for `run_nmf.py`. `run_mofaflex.py` now emits a warning when
≥ 90% of factors are active and records `n_factors_requested` / `n_factors_active` /
`ceiling_binding` in the manifest (`run_mofaflex.py:505-519`) — **but no manifest currently on
disk carries those three keys**, because the code postdates the last write of each; they will
appear on the next fit. **A K-sweep (20 / 30 / 40 / 60) has NOT been run**, so where the active
count saturates is unknown.

Note for symmetry: `results/GBM/analysis_parameters.csv` records no selection criterion for
ALARMIST's K = 20 either. That is a statement about what is recorded, not a demonstration that
the two choices are equally unprincipled — nothing further was checked.

### Sensitivity fit at `nzf > 0.001` — where all four loop features survive

`sensitivity_nzf0.001/`: 4,608 → **1,550** features → **1,541** in **8** views (Lymphoid's 9 dropped
by `min_features`); early-stopped at **199 / 1000**; 20/20 factors active.

| feature | peak factor | weight | rank |
|---|---|---|---|
| `mGAM^GRN^SORT1` | **Factor 19** | **+0.651** | 73 / 1541 |
| `MES-like^ANXA1^FPR1` | **Factor 7** | **−0.278** | 136 / 1541 |

**Different factors — the GP-prior factorisation does not put the two arms of the loop on one
program.** Factor 7 does carry three of the four loop features with concordant sign
(`mGAM^GRN^SORT1` −0.398, `mGAM^ANXA1^FPR1` −0.308, `MES-like^ANXA1^FPR1` −0.278) and is
myeloid-anchored (its largest weight is `mGAM^C3^C3AR1` at **−2.892**).

⚠️ **Correction from verification — soften the "coherent mGAM-centred axis" reading.** Factor 7 is
the **eighth** strongest factor for `mGAM^GRN^SORT1`, not its second. Ordered by |weight|:
F19 0.651, F15 0.590, F1 0.558, F16 0.556, F17 0.531, F9 0.451, F6 0.413, **F7 0.398**. The
concordant signs are real; the *preference* is not.

### Punch-level grade test — null, at the correct replicate unit

**Nothing significant at BH q < 0.05 in either fit.** The test runs over the **active** factors,
not the requested K, so the denominators differ.

| fit | factors tested | smallest raw p | q |
|---|---|---|---|
| primary | **17** | Factor 18, **0.013986** | 0.217949 |
| sensitivity | **20** | Factor 11, **0.008159** | 0.163170 |

*(This read "0 of 20 factors … in both fits". `mofaflex_inflow/data/factor_grade_punch_mannwhitney.csv`
has 17 rows; corrected 2026-08-06. The conclusion is unchanged.)*

Same hard floor as everywhere else in this comparator: a two-sided rank test on **7 high vs 6 low**
punches cannot go below **p = 0.0011655**. Underpowered; not evidence of absence.

### Version gap — recorded, not guessed

`inflow_mofaflex.ipynb` cell 5 states verbatim that it *"uses the MOFA-Flex `0.2.0` API, which is not
yet released on PyPI"* and instructs installing from git main. The installed build is
**`0.1.0.post2.dev179+g9792b435f`** from git main. **Every symbol the notebook uses exists in the
installed build with matching argument names**, so nothing had to be guessed — but the version
strings do not match and that is on record.

### Outputs

`data/`: `mofaflex_loadings.csv`, `factor_scores.csv.gz`, `lr_of_interest_loadings.csv`,
`r2_per_view.csv`, `r2_per_factor_view.csv`, `active_factors.csv`,
`factor_grade_punch_mannwhitney.csv`, `factor_means_per_punch.csv`, `view_feature_counts.csv`,
`moranI_features.csv`, `inflow_means_by_receiver.csv.gz`, `factor_interactions_Factor*.csv`.
`models/`: the HDF5 fit plus a two-run `determinism_probe/`. `plots/`: **15 PNGs** — `qc_inflow_distributions`,
`data_overview`, `variance_explained`, `variance_explained_by_view`, `factor_correlation`,
`top_weights`, `dotplot_focus_factors`, `circle_plot_Factor{18,5}`, `spatial_focus_factors`,
`umap_focus_factors`, `spatial_leiden_annotations`, `umap_leiden_annotations`,
`focus_factors_by_region`, and `factor_by_punch_grade` (**ours**, not a tutorial figure).
`sensitivity_nzf0.001/` mirrors the whole tree, 15 PNGs again (its circle plots are
`Factor{16,3}`) — **30 across these two fits**. *(This said 14 PNGs; recounted 2026-08-04,
`dotplot_focus_factors.png` is the 15th, and 15 + 15 = 30 is what the whole-tree sweep reports.)*
*(2026-08-06: the sibling `mofaflex_inflow_reachnorm/` adds **32** more — 13 of the same figures
plus a `circle_plot_Factor*.png` for every active factor — so the MOFA-Flex branch as a whole is
now 62. The "30" above is the count for `mofaflex_inflow/` + its `sensitivity_nzf0.001/` only.)*

⚠️ **Two of those 15 were blank until 2026-08-04.** `dotplot_focus_factors.png`, here and in the
sensitivity tree, was a 4,377-byte all-white canvas: `li.pl.dotplot` returns a **plotnine** ggplot
and `run_mofaflex.py` saved `plt.gcf()` — the same return-type trap as `plot_liana_full.py`'s
`connectivity.png`, and the third script in this codebase to hit it. `run_mofaflex.py` now routes
plotnine returns through `save_gg()` and its `save_current()` refuses to write a blank canvas. Both
were regenerated from the cached model with **no refit**
(`logs/comparators/liana-mofaflex-replot{,-sens}.log`).

⚠️ **Provenance gap, and the replot widened it.** It was recorded here as "`run_manifest.json`
records `fit_seconds` but **not the epoch reached**". Both are now missing: the regeneration pass
above **overwrote both manifests** (`refit: false`), so `mofaflex_inflow/run_manifest.json` records
`fit_seconds: null`, `determinism_probe: null`, `wall_seconds: 73.6` and `peak_rss_gb: 3.24` — the
*replot's* cost, not the fit's — and the sensitivity manifest `null` / 81.5 s / 4.24 GB. The fit's
70.5 min, the 76.0 min wall, the 4.1 GB peak, the stopping epochs **632** (primary) and **199**
(sensitivity), and the determinism probe now all live only in
`logs/mofaflex_{primary,sensitivity}.log`. The manifest should carry the stopping epoch, and a
replot must not clobber fit provenance.

**RESOLVED 2026-08-04.** The clobbered keys were restored into both manifests from `logs/mofaflex_{primary,sensitivity}.log`, the only surviving source: `fit_seconds` 4,230 s / 2,022 s (70.5 / 33.7 min), `wall_seconds` 76.0 / 35.3 min, `peak_rss_gb` 4.1 / 8.3, `n_epochs` **632** / **199** against a 1000 cap (the log says *Training converged after N epochs*, so this is convergence, not a cap hit), and the primary's `determinism_probe` `{'epochs': 20, 'max_abs_weight_diff': 0.0, 'bitwise_identical': True}`. `run_mofaflex.py` now carries `fit_seconds` / `n_epochs` / `determinism_probe` / `peak_rss_gb` / `wall_seconds` forward when it reuses a cached model, and records the replot separately under `last_replot`, so this cannot recur. Each manifest also carries a `provenance_note` saying the values were restored rather than remeasured.

### ⚠️ What the three MOFA-Flex figures actually show — read before quoting one (2026-08-06)

All three were read out of the installed `mofaflex/pl/_plotting.py` rather than inferred from the
picture, because two of them are easy to describe wrongly.

**`top_weights.png`** (`mfl.pl.top_weights`, `_plotting.py:1118-1170`):

- the x axis is **`| Weight |`** — the *absolute* loading. **Sign is carried only by the glyph**,
  ⊕ for w ≥ 0 and ⊖ for w < 0, and `scale_shape_manual(..., guide=None)` means **no legend for it
  is ever drawn**. A reader who has not been told this will read every bar as positive.
- within each facet the top *n* are taken by `|weight|` and then sorted **ascending**, so the
  **largest sits at the TOP** of the panel.
- `facet_wrap("factor", scales="free")` — **every panel has its own x scale**, so bar lengths are
  **not comparable across factors**.
- these are **raw** weights: no prevalence normalisation of any kind.
- now regenerated at **top-10** (was top-5) via `--top-weights-n`.

**`variance_explained.png`** (`mfl.pl.variance_explained`, `_plotting.py:477-517`): **it is a
`geom_tile` HEATMAP, not a bar chart.** Rows = **factors**, columns = **views = sender cell
types**, fill = R², wrapped into one facet per *group* — and there is a single group here
(`group_1`), so there is one panel. Rows are ordered by total R² descending with the **largest at
the BOTTOM**. `r2_per_view.csv` is exactly the **column sum over factors** (verified to 4 d.p.
against `r2_per_factor_view.csv`). On the reachability-normalised fit the single darkest cell is
**Factor 5 × Glial-Neuronal = 0.1641**, and the **Lymphoid column is effectively blank** (largest
cell 0.00021). `variance_explained_by_view.png` is the same data at `group_by="view"` (x = group,
faceted by view) — with one group that is nine one-column facets, i.e. the less useful layout.

**`circle_plot_Factor*.png`** — **the edges are factor-independent.** The factorisation never sees
the receiver. `source` comes from the feature name (model), and *which* 10 edges are drawn comes
from `|loading|` (model) — but **`target` and the edge weight come from
`inflow_means = lrdata.to_df().groupby(obs[cell_type]).mean()`** (`run_mofaflex.py:624-638`), i.e.
the *receiving* cell's own annotation and the **raw** mean inflow, not anything factor-weighted.
Verified on the reachability-normalised run: **Factor 1 and Factor 6 share 108 of their 225
(source, LR, target) rows with edge weights identical to 0.000e+00**, while their loadings on
those same rows differ by up to **1.9001**. Consequences: two factors that share a feature display
the *identical* sub-network; the edges ignore the factor's sign; and the apparent "per-factor
network" is mostly not the factor. **Read it as "the top-10 interactions this factor selects, and
where those interactions generally go", NOT as "this factor's sender→receiver structure".** Now
emitted for **all** active factors (was 2 focus factors) via `--circle-all-factors` /
`--circle-top-n`.

---

## LIANA+ — MOFA-Flex with reachability-normalised QC (our deviation; all nine views survive)

`scripts/comparators/liana/run_mofaflex.py --nzf-mode reachability` →
`results/comparators/liana/GBM/mofaflex_inflow_reachnorm/` (2026-08-06). This is the direct answer
to the section above: *the authors' QC deletes both arms of ALARMIST motif 1*. **The default stays
`--nzf-mode global`, so the tutorial-faithful run remains reproducible; this is an additional fit,
not a replacement.**

⚠️ **CD-1 applies here too, twice over.** It factorises the same
`cellchatdb2_inflow/data/inflow_lrdata.h5ad`, so all spatial weighting already happened inside
`li.mt.inflow` at the ALARMIST-derived 13.1454 µm bandwidth — *and* the reachability denominator
is itself computed at that same support radius (`--bandwidth 13.1454 --cutoff 0.1` in the
manifest's `argv`). Both the feature set and the correction to it are conditioned on the deviating
kernel scale.

### The tutorial's global cut is an abundance filter, not an expression filter

`inflow_mofaflex.ipynb` cell 19 applies **one global `nonzero_fraction > 0.01`** to every feature.
But an inflow feature `<sender>^<lig>^<rec>` can only be non-zero for a cell that has a cell of
that **sender** inside the kernel support, so its `nonzero_fraction` is bounded above by the
sender's **reachability** — P(≥1 neighbour of type *s* within R). Measured on this dataset:

- surviving-feature count vs sender abundance is **Spearman ρ = 0.917** (p = 5.1e-4, n = 9 senders)
  — the filter is tracking how common the sender is, not how informative the feature is;
- reachability at R = 28.2096 µm: NPC-like **0.761**, AC-like 0.750, OPC-like 0.748,
  MES-like 0.727, Glial-Neuronal 0.379, mGAM **0.319**, non-mGAM 0.205, Vascular **0.157**,
  Lymphoid **0.025**;
- **Lymphoid cannot pass the cut for arithmetic rather than biological reasons.** The largest
  `nonzero_fraction` of any Lymphoid feature is **0.004332** — below 0.01 — so **zero** Lymphoid
  features can survive whatever the biology is;
- expressed *as a fraction of the cells that could receive it*, the same 0.01 cut is **1.3%** for
  NPC-like and **39%** for Lymphoid;
- **Vascular is penalised for clustering, not for rarity** — 3.2% of cells but only 15.7%
  reachable, because vessels are spatially aggregated.

Views the tutorial recipe consequently loses: **mGAM (24 features, exactly ONE short of
`min_features=25`), non-mGAM (4), Lymphoid (0)**.

### The deviation, and what it recovers

`nzf_norm = nonzero_fraction / reach[sender]`, keep `> 0.01` — "is this feature non-zero in an
appreciable share of the cells that *could* receive it", which is scale-free in sender abundance.

| | tutorial QC (`global`) | **reachability-normalised** |
|---|---|---|
| features after cut | 4,608 → **447** | 4,608 → **779** |
| SVI filter | no-op | **no-op again** (779 → 779, measured) |
| views ≥ `min_features=25` | **6 of 9** | **9 of 9** |
| NPC-like / MES-like / AC-like / OPC-like | 98 / 97 / 92 / 78 | 114 / 122 / 111 / 92 |
| Glial-Neuronal | 28 | 63 |
| **Vascular** | 26 | **138** (5.3× recovery) |
| **mGAM** | **24 → view dropped** | **66** |
| **non-mGAM** | **4 → view dropped** | **34** |
| **Lymphoid** | **0 → absent** | **39** |

**Fit.** K = 20 requested, CPU, batch 2048, lr 0.005, patience 50, seed 0; **converged after 294
epochs** against a 1000 cap (`logs/comparators/liana-mofaflex-reachnorm.log`, *"Training converged
after 294 epochs"* — convergence, not a cap hit). Fit **40.1 min**, total wall **41.5 min**, peak
RSS **6.78 GB**. **19 / 20 factors active** — see the ceiling warning above; this is not a
discovered number. No determinism probe was run for this fit (`determinism_probe: null`).

R² per sender view — **and the view the tutorial deleted is the second best**:

| view | R² | | view | R² |
|---|---|---|---|---|
| Glial-Neuronal | **0.5058** | | Vascular | 0.2708 |
| **mGAM** | **0.4453** | | AC-like | 0.2629 |
| NPC-like | 0.3931 | | non-mGAM | 0.1397 |
| MES-like | 0.3542 | | **Lymphoid** | **0.0013** |
| OPC-like | 0.3254 | | | |

⚠️ **Caveat that must travel with this run: Lymphoid is admitted but uninformative.** Its R² is
**0.0013**, its 39 features rest on **2.5% reachability**, and its largest single factor×view cell
is 0.00021. The criterion lets it in; the model explains essentially nothing there. **Do not
interpret Lymphoid factors.** No absolute cell-count floor was applied — the user asked for the
normalisation alone — so nothing else guards against a thin cell base.

**Punch-level grade test — null again, same replicate unit, same floor.** **0 of 19** factors at
BH q < 0.05; smallest raw p **0.022145** (Factor 18) → q **0.332168**
(`data/factor_grade_punch_mannwhitney.csv`). Underpowered as everywhere else: 7 high vs 6 low
punches floors a two-sided rank test at p = 0.0011655. **Restoring the mGAM view did not change
the grade result, and it did not put the two motif-1 arms on one factor** — see next.

### Where the two motif-1 arms actually rank in this fit

Ranks are within-factor by `|weight|` over all **779** features
(`data/mofaflex_loadings.csv`):

| feature | peak factor | weight | rank | next strongest |
|---|---|---|---|---|
| `mGAM^GRN^SORT1` | **Factor 19** | **+0.329** | **67 / 779** (top 8.6%) | F7 −0.327 (80), F9 −0.300 (62), F16 −0.286 (119), F1 +0.279 (99) |
| `MES-like^ANXA1^FPR1` | **Factor 1** | **+0.078** | **276 / 779** (top 35.4%) | F7 −0.071 (305), F15 −0.054 (308) |
| `MES-like^GRN^SORT1` (reverse) | Factor 4 | +0.258 | 121 / 779 | F12 +0.148 (169) |
| `mGAM^ANXA1^FPR1` (reverse) | Factor 1 | +0.165 | 164 / 779 | F7 −0.150 (176) |

Three things follow, and none of them is favourable:

1. **`mGAM^GRN^SORT1` is spread flat across factors with mixed signs** — its top **six** factors
   lie within **0.057** of each other, and their signs alternate — so **no factor claims it**.
2. **`MES-like^ANXA1^FPR1` is nowhere**: roughly two-thirds of all features load more strongly on
   its own best factor than it does.
3. **The autocrine-ish myeloid direction `mGAM^ANXA1^FPR1` outranks the biologically meaningful
   `MES-like^ANXA1^FPR1` on 17 of the 20 factors**, and its peak rank (164) is far ahead of the
   MES-like arm's (276).

**Neither arm is a top-10 feature of any factor**, so neither ever appears in `top_weights.png`.
The best joint factor by worst-of-the-two rank is **Factor 1** (ranks 99 and 276), then F6, F7,
F15 — all four carry the **same sign** for both arms (in fact 15 of 20 factors do). So the model
is **not** placing the two arms at opposite poles; it simply places neither anywhere prominent.

### Outputs

| | |
|---|---|
| `data/` (**30**) | `mofaflex_loadings.csv`, `factor_scores.csv.gz`, `lr_of_interest_loadings.csv`, `r2_per_view.csv`, `r2_per_factor_view.csv`, `active_factors.csv`, `factor_grade_punch_mannwhitney.csv`, `factor_means_per_punch.csv`, `view_feature_counts.csv`, `moranI_features.csv`, `inflow_means_by_receiver.csv.gz`, plus **19** `factor_interactions_Factor*.csv` (one per active factor) |
| `models/` (1) | `mofaflex_inflow_reachnorm0.01.hdf5`. **No `determinism_probe/`** — not run for this fit |
| `plots/` (**32** PNG) | `qc_inflow_distributions`, `data_overview`, `variance_explained`, `variance_explained_by_view`, `factor_correlation`, `top_weights` (now top-10), `dotplot_focus_factors`, `spatial_focus_factors`, `umap_focus_factors`, `spatial_leiden_annotations`, `umap_leiden_annotations`, `focus_factors_by_region`, `factor_by_punch_grade` (**ours**, not a tutorial figure), plus **19** `circle_plot_Factor*.png` — read them with the caveat above |
| | `run_manifest.json` (carries `fit_seconds` / `peak_rss_gb` / `wall_seconds` and a separate `last_replot` block, so the 2026-08-04 clobbering defect cannot recur) |

The whole `results/comparators/liana/` tree now holds **466 PNGs** (was 420): +32 here,
+12 under `vs_alarmist/`, +2 under `bandwidth_choice/`. None of the 46 new files is under 20 kB,
so none is a repeat of the blank-canvas trap.

---

## LIANA+ — annotating the NMF factors (PROGENy + CellChatDB pathways)

`scripts/comparators/liana/annotate_factors.py` →
`results/comparators/liana/GBM/factor_annotation/` (**0.22 min**, 0.83 GB, **no re-fit** — it reads
`nmf_{bivariate,inflow}/data/NMF_H_loadings.csv` only). Methods: `dc.mt.mlm` (the authors' choice in
`mofatalk`) and `dc.mt.ulm` for robustness, decoupler **2.2.0**, BH across pathways within each
factor. Deviation from `mofatalk`: the LR universe is CellChatDB v2 rather than
`select_resource('consensus')`.

### ⚠️ PROGENy is a poor fit here and must be reported as one

| branch | LR pairs covered by the PROGENy LR set |
|---|---|
| bivariate | **47 / 131 (35.9%)** |
| inflow | **160 / 524 (30.5%)** |

And most of what it returns is an artifact of a single adhesion pair. Of the inflow MLM hits at
FDR < 0.05, **6 are WNT and 5 are p53, every one with negative t** — and for F2/F4/F6 these are
restatements of *"this factor loads `NCAM1^NCAM1`"* (PROGENy weight `WNT ← NCAM1^NCAM1` =
**−0.877**; those factors' NCAM1 shares are 23.5% / 39.9% / 30.5%). For F7 the driver is instead
`DLL3^NOTCH1` (`WNT` weight **−0.626**), since **F7's NCAM1 share is exactly 0%**.

**Only the Hypoxia hits are defensible:**

| branch | factor | t | FDR |
|---|---|---|---|
| bivariate | F4 | **3.55** | 3.2e-3 |
| inflow | F5 | **3.42** | 3.9e-3 |
| inflow | F3 | **2.74** | 2.7e-2 |

### CellChatDB pathway composition carries the annotation instead — 100% coverage

| branch | factor | call | evidence |
|---|---|---|---|
| bivariate | **F1** | NOTCH / JAG1 lateral induction | NOTCH 28.7% of loading; ULM t = **4.01**, FDR 3.7e-3 |
| bivariate | **F2** | complement + galectin, myeloid-immune | GAS 18.2%, COMPLEMENT 14.4%, GALECTIN 11.9% — **composition only, no FDR < 0.05 hit** |
| bivariate | **F3** | basement-membrane COLLAGEN | 33.2%; t = **9.65**, FDR **2.3e-15** — the strongest call anywhere |
| bivariate | **F4** | FGF–tenascin | FGF 24.4%, TENASCIN 11.8%; t = **5.03**, FDR 5.8e-5. ⚠️ the third-largest component is **APP at 10.5%, NOT THBS** |
| bivariate | **F5** | synaptic NRXN | 30.7% |
| bivariate | **F6** | DLL3–NOTCH | NOTCH **55.5%**; t = 6.90, FDR 7.5e-9 |
| inflow | **F7** | NOTCH | 50.6%; t = **9.73**, FDR **1.1e-17** |

⚠️ **The CellChatDB column is named `pathway`, not `pathway_name`.**

### ⚠️ Structural finding — the inflow NMF factors are mostly sender identity plus one adhesion pair

**Six of seven inflow factors are ≥ 75% a single sender:** F1 Glial-Neuronal **92.9%**, F2 OPC-like
**92.2%**, F3 MES-like **82.5%**, F4 NPC-like **86.4%**, F5 AC-like **75.6%**, F6 Glial-Neuronal
**93.5%**. Only F7 is mixed (NPC-like 63.2%).

**And five of seven share the same top feature, `X^NCAM1^NCAM1`** — F2, F3, F4, F5, F6. (F1's top is
`Glial-Neuronal^CNTN2^CNTN2`; F7's is `NPC-like^DLL3^NOTCH1`.) *(Verification corrected this from a
reported 4.)*

**So these are largely sender identity plus one dominant adhesion pair, not distinct communication
programs.** That is a stronger statement of the same point the *How much does the program structure
change?* section makes from `top10_loadings_per_factor.csv`.

### Neither branch co-loads the two arms of motif 1

| arm | inflow | bivariate |
|---|---|---|
| GRN→SORT1 | **F1** (92.9% Glial-Neuronal) and **F6** (93.5% Glial-Neuronal) | **F1** (loading 3.042, rank 7/131) |
| ANXA1→FPR1 | **F3** (82.5% MES-like) | **F2** (loading 0.851, rank 19/131) |

**And the strongest `GRN^SORT1` feature anywhere in the inflow decomposition is
`Glial-Neuronal^GRN^SORT1` — rank 9 of 2,704 on F1 — not the mGAM one** (`mGAM^GRN^SORT1` peaks on
F6). CellChatDB pathways: `GRN^SORT1` → **GRN** (Secreted Signaling, a one-member pathway);
`ANXA1^FPR1` → **ANNEXIN**.

### Outputs

`data/` (**23** CSVs — recounted on disk 2026-08-04; this said 24) and `plots/` (10 PNGs), covering PROGENy MLM/ULM, CellChatDB ULM,
pathway-fraction and loading-sum matrices, the top-25 annotated features per branch, the inflow
sender-fraction matrix, and the focus-LR loading tables. `factor_calls_summary.csv` is the digest.

---

## LIANA+ — the converging conclusion, stated once

**Four independent LIANA analyses now agree that LIANA does not reconstruct the ALARMIST motif-1
mGAM ⇄ MES-like loop as a single program.**

| analysis | what it says about the loop |
|---|---|
| inflow global specificity | **ANXA1→FPR1 (MES-like→mGAM) rank 2 / 81, p = 0.000999 (floor); GRN→SORT1 (mGAM→MES-like) p = 1.00000.** One arm only. |
| NMF (bivariate and inflow) | the arms land on **different factors with different dominant senders**; the strongest `GRN^SORT1` feature is Glial-Neuronal, not mGAM |
| MOFA-Flex | **different peak factors** (F19 vs F7) at `nzf > 0.001`; and at the tutorial's own QC **both arms are deleted outright**. With the reachability-normalised QC that restores the mGAM view, they *still* land differently (F19 vs F1) and **neither is a top-10 feature of any factor** |
| LRIC / cross-PCF | GRN→SORT1 is genuinely directional (7/8 punches, p = 0.0156) but the **LRIC/cross-PCF ratio is ~1 in every bin**, so co-occurrence is fully explained by cell-type co-location |

### ⚠️ The dominant reason is the UNIT OF ANALYSIS — feature indexing is secondary (2026-08-06)

**This corrects the emphasis carried by the four rows above and by `DEVIATIONS.md`.** Those rows
describe *symptoms* of the factorisations; the cause sits one level upstream, in what a **row** of
LIANA's matrix is. `scripts/comparators/liana/why_no_mgam_motif.py` →
`vs_alarmist/why_no_mgam_motif.json`.

**Inflow scores the RECEIVING cell**, so the two arms of a bidirectional loop land on **disjoint
populations** by construction:

| feature | non-zero cells | top receivers |
|---|---|---|
| `mGAM^GRN^SORT1` | **3,183 (3.18%)** | MES-like 19.2%, AC-like 17.0%, NPC-like 16.7%, Glial-Neuronal 16.0% |
| `MES-like^ANXA1^FPR1` | **944 (0.94%)** | **mGAM 45.2%**, MES-like 22.7% |
| **both** | **95 (0.095%)** | Pearson **+0.0177**, Spearman **+0.0385** |

Now hold the two interactions, the tissue and the LR database fixed and **change only the unit**:

| | cells (LIANA inflow) | 50 µm patches (ALARMIST) |
|---|---|---|
| rows | 100,190 | 13,113 |
| rows carrying **both** arms | **0.095%** | **1.235%** |
| Pearson *r* | **+0.0177** | **+0.4562** |
| Spearman ρ | +0.0385 | +0.4044 |
| P(arm 2 \| arm 1) ÷ marginal | **3.2×** | **14.2×** |

**Pearson rises 26-fold purely by aggregating cells into patches.** The loop is a property of a
*neighbourhood*, and **a neighbourhood is not a row in LIANA's matrix** — so no factorisation
applied to that matrix can recover structure the input does not contain. Feature indexing (sender
inside the feature name, receiver nowhere) is real and compounds this, but it is the **secondary**
effect: correcting the QC to restore every sender view (see the reachability-normalised section)
changed neither the co-loading nor the grade result.

Two further points that only fall out of the whole set:

- **LRIC is the only branch that resolves direction**, and it is the only one that supports the
  forward arm as directional at the correct replicate unit.
- **LRIC also refutes the "the 28.2 µm kernel missed a longer-range interaction" explanation**: its
  ratio stays at ~1 out to **225 µm**, 8× the bivariate/inflow support radius.

The one dissent is **LR-MISTy**, which finds GRN a top-10 predictor of both SORT1 and FPR1 — but
that branch is a niche-composition predictor rather than an LRI detector (see its section), so it
does not overturn the above. **CellChat recovers both arms** (see its section), so the disagreement
remains between comparators, not between ALARMIST and the field.

---

## LIANA+ — housekeeping and open issues (2026-08-04)

- **`scripts/comparators/` is STILL UNTRACKED IN GIT.** `git ls-files scripts/comparators` returns
  nothing, and `git check-ignore` says it is not ignored either — it is simply never added.
  **This has now cost verifiability three separate times today**: an edit can be confirmed to have
  the intended content, but there is no way to prove no other line moved. **Open issue — add
  `scripts/comparators/` to git.** Every run manifest records `git_sha: 95208de`, which pins the
  *package*, not the comparator scripts that produced the results.
- `env.lock.yml` re-frozen after installing torch / decoupler / mofaflex — **167 lines**, carrying
  `decoupler==2.2.0`, `torch==2.13.0`, `pyro-ppl==1.9.1`, `gpytorch==1.15.2`,
  `mofaflex==0.1.0.post2.dev179+g9792b435f`.
- `run_inflow_downstream.py --db` now accepts the literal **`consensus`** (it was file-path only),
  which is what let the default tier reach `rank_aggregate`.
- `plot_liana_full.py` previously counted `.DS_Store` in `n_files`; it now excludes dotfiles and
  reports `n_png` separately (`plot_liana_full.py:479-480`). **Both `plots_full` manifests now match
  disk exactly: 30 files / 25 PNG (bivariate) and 62 / 55 (inflow)** — verified by `find`. The
  earlier warning that `n_files` reads "31 / 63" is now **historical**; see the corrected note in
  *`plots/` vs `plots_full/`*.
- ⚠️ A stray **`.DS_Store` was written into `cellchatdb2_inflow/`** during the SVI measurement
  (2026-08-04 19:40), and one sits at `results/comparators/liana/GBM/` too. Harmless, but it is what
  inflated the old counts — delete rather than count around them.
- **Two more blank PNGs were found and fixed (2026-08-04, later than everything above).**
  `mofaflex_inflow/plots/dotplot_focus_factors.png` and its `sensitivity_nzf0.001/` twin were
  4,377-byte all-white canvases — `li.pl.dotplot` returns a plotnine ggplot and `run_mofaflex.py`
  saved `plt.gcf()`, the **third** script in this codebase to hit that trap. `run_mofaflex.py` now
  has `save_gg()` for plotnine returns and a `save_current()` that **refuses** to write a blank
  canvas; both figures were regenerated from the cached model with no refit. **Four blanks total
  across three scripts; the final whole-tree sweep is 420 PNGs, 0 blank.** *(That sweep is dated
  2026-08-04. The tree now holds **466** PNGs — the 46 added on 2026-08-06 are the 32 under
  `mofaflex_inflow_reachnorm/plots/`, 12 under `vs_alarmist/figures/` and 2 under
  `bandwidth_choice/figures/`; none is under 20 kB, so none is a blank canvas, but they were not
  re-swept by the original blank-detection pass.)*
- MOFA-Flex `run_manifest.json` does not record the stopping epoch — and since the replot pass
  above overwrote it, it no longer records `fit_seconds` or the determinism probe either, and its
  `wall_seconds` / `peak_rss_gb` now describe the replot. The fit's provenance survives only in
  `logs/mofaflex_{primary,sensitivity}.log` (see that section). **A replot must not clobber fit

**RESOLVED 2026-08-04.** The clobbered keys were restored into both manifests from `logs/mofaflex_{primary,sensitivity}.log`, the only surviving source: `fit_seconds` 4,230 s / 2,022 s (70.5 / 33.7 min), `wall_seconds` 76.0 / 35.3 min, `peak_rss_gb` 4.1 / 8.3, `n_epochs` **632** / **199** against a 1000 cap (the log says *Training converged after N epochs*, so this is convergence, not a cap hit), and the primary's `determinism_probe` `{'epochs': 20, 'max_abs_weight_diff': 0.0, 'bitwise_identical': True}`. `run_mofaflex.py` now carries `fit_seconds` / `n_epochs` / `determinism_probe` / `peak_rss_gb` / `wall_seconds` forward when it reuses a cached model, and records the replot separately under `last_replot`, so this cannot recur. Each manifest also carries a `provenance_note` saying the values were restored rather than remeasured.
  provenance.**

### Code added in the 2026-08-06 pass — `scripts/comparators/liana/`

| script | what it produces |
|---|---|
| `choose_bandwidth.py` | the bandwidth evidence + the two `query_bandwidth` figures → `GBM/bandwidth_choice/` (2 PNG, 4 CSV) |
| `why_no_mgam_motif.py` | the cell-vs-patch unit demonstration → `vs_alarmist/why_no_mgam_motif.json` |
| `compare_programs_to_alarmist.py` | per-cell Spearman, LIANA programs vs ALARMIST motifs → `vs_alarmist/comparison_summary.json` + `data/rho_*.csv` |
| `cosine_factors_vs_motifs.py` | factor-vs-motif cosine, 3 sign modes, permutation null → `vs_alarmist/` (12 PNG + PDF + SVG) |
| `build_report_figures.py` | assembles `reports/liana_plus_GBM_cellchatdb2/` — **159 PNGs** across **5** tutorial sub-trees (`bivariate`, `inflow`, `lric`, `misty`, `mofaflex`) + `figure_manifest.json` |
| `build_liana_report.py` + `_liana_report_sections.json` | builds `reports/liana_plus_GBM_cellchatdb2/liana_plus_GBM_methods.html` |

`run_mofaflex.py` gained `--nzf-mode {global,reachability}` (default **`global`**, so the
tutorial-faithful path is untouched), `--bandwidth`, `--cutoff`, `--bandwidth-support-radius`,
`--xy-sep`, `--top-weights-n`, `--circle-all-factors`, `--circle-top-n`, the `n_factors`-ceiling
warning, and a fit-provenance carry-forward so a cached-model replot no longer clobbers
`fit_seconds` / `n_epochs` / `determinism_probe`.

### Open issues after the 2026-08-06 pass

- ❌ **CD-1 stands.** Bandwidth **13.1454 µm** remains a `SKILL.md:45-46` violation, kept by user
  decision after the evidence review — **not** resolved by it.
- ❌ **CD-2 stands.** No native multi-sample / differential mode; the punch-level test is a
  documented hand-rolled substitute, and it is null in **five** places now (4 previously + the
  reachability-normalised fit's 0 / 19).
- ❌ **No K-sweep for MOFA-Flex.** The ceiling binds in all three fits and nobody has looked for
  where the active count saturates.
- ❌ **No prevalence normalisation of the MOFA-Flex loadings.** ALARMIST has `V*`; LIANA has no
  equivalent, so the cosine comparison normalises one side only.
- ❌ **`compare_programs_to_alarmist.py` NMF index bug** — reads 8 numeric columns from a rank-7
  `NMF_W_factor_scores.csv` (written `index=False`). That row of the summary is **unverified**.
- ✅ ~~`cellchatdb2_inflow/plots/global/bandwidth_query.png` draws its guide line at the σ on a
  radius axis.~~ **NOT AN ISSUE — retracted 2026-08-07.** The guide is at the support radius
  28.2096 µm and is correct; only the annotation is clipped. See the retraction in the inflow
  bandwidth section above. This was never verified before being listed.
- ❌ **`scripts/comparators/` is still UNTRACKED IN GIT** (unchanged from the entry above).
- ❌ **`reports/liana_plus_GBM_cellchatdb2/liana_plus_GBM_methods.html` was NOT updated by this
  pass.** It is built from `_liana_report_sections.json` and is now **stale** relative to this
  document.

## NICHES — R, v1.2.4, env `comp-niches`

Local source `/Users/jiayifan/tansey_lab/NICHES` @ `d698e37b8c38ebd103c34acd0b35b03b48d3c5a3`
(2026-01-29). R 4.3.3, Seurat 5.3.0, SeuratObject 5.2.0. Code
`scripts/comparators/niches/`, contract `NOTES.md`, departures `DEVIATIONS.md`.

**Core algorithm** — NICHES is a *transformation*, not a test. For a directed edge
(sending cell *i* → receiving cell *j*) and mechanism *L→R* it computes the **product of the
normalized ligand expression in *i* and the normalized receptor expression in *j***, with
multi-subunit complexes multiplied across subunits (`Reduce('*', subunit.list)`,
`RunNeighborhoodToCell.R:38,55`). `NeighborhoodToCell` averages that product over every edge
landing on *j* (`blend = "mean"`), giving one **niche vector per cell**; `CellToCellSpatial`
keeps one vector **per edge**. The result is handed back as an ordinary Seurat assay
(mechanisms as "genes", cells or edges as "cells"). **There is no null model, no permutation
and no p-value anywhere in the scoring step.** Every inferential statement comes from whatever
Seurat test you then run on that matrix — here the vignettes' `FindAllMarkers(test.use = "roc")`.
So a "significant interaction" in NICHES is a *downstream* claim about a group contrast, and
the unit it attaches to is the receiving cell (`NeighborhoodToCell`) or the cell–cell edge
(`CellToCellSpatial`), never the interaction itself.

**Spatial model** — **mutual k-nearest-neighbour graph, `k = 4` (package default), no kernel
and no distance cutoff.** `ComputeEdgelist.R:46-56` ranks every cell's neighbours by euclidean
distance, keeps the `k+1` nearest, and symmetrises with `adj & t(adj)`, so an edge survives
only if each cell is in the other's top 4. `rad.set` offers a hard radius instead but is
ignored whenever `k` is non-NULL, and `RunNICHES` always passes `k`. Two consequences specific
to NICHES, both measured on our data (`summary_edge_audit.csv`, `audit_edges.py`):

* **`k` is unitless, so NICHES has no radius to quote** — the neighbourhood must be measured
  after the fact. Pooled over all 13 cores the real (non-self) edges have **median 10.1 µm,
  p95 26.6 µm, p99 37.8 µm, max 243.1 µm**. That is far and away the **tightest** neighbourhood
  of any method here (LIANA+ 28.2 µm, SpatialDM 135 µm, CytoSignal 200 µm, stLearn 250 µm,
  COMMOT 365 µm) — NICHES is effectively scoring *directly abutting* cells.
* **The radius is density-dependent, and density tracks grade.** Because `k` is fixed, a denser
  core gets a physically *smaller* neighbourhood: per-core median edge length runs from
  **8.1 µm (core1, 26,456 cells) to 22.2 µm (core6, 3,092 cells)**, a 2.7× spread. High-grade
  cores average 3.39× more cells than low-grade cores, so the high- and low-grade arms of this
  dataset are **not being measured at the same physical scale**. See the confound note below.
* **Every cell is its own neighbour.** `order(dis_vec)[1:(k+1)]` includes the cell itself at
  distance 0 and the mutual-NN step always keeps it, so the edge list contains one self-edge
  per cell: **100,197 of 419,299 edges (23.9%)** are `cell—itself`. `NeighborhoodToCell`
  therefore mixes each cell's own **autocrine** ligand×receptor product into its "neighbourhood"
  average. This is package behaviour, not a configuration choice, and it is not mentioned in
  the vignettes.

**LR database** — default is **FANTOM5** (`ncomms8866_human`, bundled `.rda`); OmniPath is the
other built-in. For the `cellchatdb2` tier we pass `LR.database = "custom"` with
`custom_LR_database` = `data/LRdatabase/CellChatDBv2.0.human.csv`. **The conversion is trivial** —
`LoadCustom` wants exactly a 2-column data.frame whose first column is the ligand subunits and
second the receptor subunits, each `_`-separated, which is already the format of our export, so
the "conversion" is a straight `db[["ligand","receptor"]]` selection with no complex
re-encoding. Complexes are handled by **multiplying subunit expression**, and
`FilterGroundTruth` keeps a mechanism only if **every** subunit is present in the object. On the
5,119-gene Xenium panel that takes CellChatDB v2 from **3,218 unique pairs → 1,088 mechanisms**,
identical in every core. `species` is silently ignored in custom mode (`LoadCustom.R:12`).
The FANTOM5 `default` tier was **not run** (this pass is `cellchatdb2`-only, matching stLearn /
SpatialDM / COMMOT / LIANA+).

**Input** — a Seurat object with (a) **raw counts**, normalized by `NormalizeData` inside the
runner (the vignettes normalize themselves, so counts must be handed in un-normalized — we read
`layers['counts']`, since this h5ad's `X` is already log-normalized); (b) **`x` and `y` as
ordinary `meta.data` columns**, not an `@images` object — units are irrelevant because `k` is
unitless, but ours are µm from `obsm['spatial']`; (c) a **cell-type column** (`cell_type`,
9 types), which `RunNICHES.Seurat` copies into `Idents`. No QC filter is applied:
`min.cells.per.ident` and `min.cells.per.gene` are left at their `NULL` defaults, so all
100,197 cells and all 5,119 genes enter.

**Workflow**

| Step | Call with argument values | Produces |
|---|---|---|
| 1 | `CreateSeuratObject(counts, meta.data)` — per TMA core | 5,119 × n_cells object |
| 2 | `NormalizeData(obj)` — `LogNormalize`, `scale.factor = 1e4` | `RNA` `data` layer |
| 3 | *(alra sub-run only)* `RunALRA(obj)` — `k = NULL` (auto), `q = 10`, `use.mkl = FALSE` | `alra` assay |
| 4 | `RunNICHES(assay = "RNA"\|"alra", LR.database = "custom", custom_LR_database = CellChatDBv2, species = "human", cell_types = "cell_type", position.x = "x", position.y = "y", k = 4, rad.set = NULL, blend = "mean", min.cells.per.ident = NULL, min.cells.per.gene = NULL, CellToCellSpatial = TRUE, NeighborhoodToCell = TRUE, CellToCell = FALSE, output_format = "seurat")` | 2 Seurat objects per core |
| 5 | tag `Condition <- grade`, `Core <- tma_id`; `merge()`; `JoinLayers()` | one 13-core object |
| 6 | *(CellToCellSpatial only)* `subset(nFeature_CellToCellSpatial > 5)` | vignette 04's low-information filter |
| 7 | `ScaleData` → `FindVariableFeatures(selection.method = "disp")` → `RunPCA(npcs = 100)` → `RunUMAP(dims = 1:50)` | embedding |
| 8 | `FindAllMarkers(min.pct = 0.25, only.pos = TRUE, test.use = "roc")` with `Idents` = ReceivingType, then = grade | marker tables |
| 9 | per ReceivingType: `subset` → `ScaleData` → `FindVariableFeatures` → `RunPCA(npcs = 50)` → `RunUMAP(dims = 1:40)` → `FindAllMarkers` → `DoHeatmap(top_n(20, myAUC))` | per-population grade contrast |

**Data outputs** (per core, per imputation sub-run)

| File | Shape / schema | Meaning |
|---|---|---|
| `objects/{NeighborhoodToCell,CellToCellSpatial}.rds` | Seurat object | full NICHES output; every downstream step re-reads this, so the O(N²) edgelist is never recomputed |
| `quant/<org>_scores.mtx` | 1,088 × n_columns sparse | the LR-product score matrix itself |
| `quant/<org>_features.tsv` | 1,088 lines | mechanism names, em-dash joined |
| `quant/<org>_columns.tsv` | n_columns lines | receiving cell (N2C) or `sender—receiver` edge (C2CS) |
| `quant/<org>_metadata.csv` | n_columns × meta | `ReceivingType`/`SendingType`/`VectorType`, grade, tma_id, x, y |
| `quant/<org>_mechanism_summary.csv` | 1,088 rows | mean/max score, n and frac nonzero per mechanism |
| `run_manifest.json` | — | every parameter, versions, seed, wall time, peak heap, git SHA |
| `DONE` | sentinel | makes the whole pipeline resumable |

Tier-level: `summary_per_core.csv`, `summary_per_tier.csv`, `summary_confound.csv`,
`summary_requested_lr.csv`, `summary_edge_audit.csv`, `run_timings.csv`.
Analysis-level (`_analysis/<org>/`): `merged.rds`, `mechanism_detection.csv`,
`markers/markers_*.csv`, `differential_summary.csv`, `composition_*.csv`,
`requested_*_by_ReceivingType_grade.csv`, `analysis_manifest.json`.

**Image outputs** — all as `.png` + `.pdf` (+ `.svg` under 50k points; skipped above that and
logged, since an SVG of a 100k-point UMAP is hundreds of MB).

| Plot | Shows | File |
|---|---|---|
| `ElbowPlot` | PC variance, 100 PCs | `plots/elbow.*` |
| `PCHeatmap` | loadings, PCs 40–48 | `plots/pc_heatmap.*` |
| `DimPlot` × ReceivingType / Condition / Core | niche UMAP | `plots/umap_<var>.*` |
| `DimPlot` × VectorType / SendingType (C2CS) | edge UMAP | `plots/umap_<var>.*` |
| `DoHeatmap` per ReceivingType (vignette 01) | niche mechanisms per cell type | `plots/heatmap_ALL_by_ReceivingType.*` |
| `DoHeatmap` grade contrast | top-20 by `myAUC`, high vs low | `plots/heatmap_ALL_by_grade.*`, `plots/heatmap_<type>.*` |
| `DimPlot` per population | that population's niche, split by grade | `plots/umap_<type>_by_grade.*` |
| `FeaturePlot`, method's own top LRs | top `myAUC` mechanisms per population | `plots/top_lr/` |
| `FeaturePlot` + `VlnPlot` ×2, requested LRs | GRN—SORT1, ANXA1—FPR1 | `plots/requested_lr/` |
| spatial scatter, requested LRs | niche score in situ, faceted by core | `plots/requested_lr/*_spatial.*` |

Not produced: `SpatialFeaturePlot` (needs a Seurat `@images` object the Xenium h5ad has none of
— replaced by an equivalent ggplot scatter, D7), and `ALRAChooseKPlot` (diagnostic for ALRA's
rank choice; the chosen `k` is recorded per core instead).

**Multi-sample / differential mode** — **native, and it is the split-run-merge pattern**, not a
joint model. Vignettes 04 and 07 both say to *split the object by condition first*, run NICHES
separately on each, tag the outputs, `merge()`, and then do ordinary Seurat differential testing
on the merged matrix. We apply it at the TMA-core level (13 runs) and contrast `grade` on the
merged object. There is no shared latent space, no batch term and no random effect for core —
**every cell is an independent observation in the marker test**, which is pseudoreplication with
respect to the 13 cores. That is the method's own design, not our choice.

**Gotchas**

1. **`ComputeEdgelist` is dense O(N²) and there is no fast path.** `apply(df, 1, ...)`
   (`ComputeEdgelist.R:36`) materialises the full N × N double distance matrix, then `adj_mat`,
   `t(adj_mat)` and `1*(...)` — roughly 4 × N² × 8 bytes. At the whole slide (100,197 cells)
   that is **~80 GB for the distance matrix alone**. The `nn.method = 'aoz'` spNNGP fast path
   exists in the signature but its body is **commented out** (`ComputeEdgelist.R:76-131`) and
   `RunNICHES` never passes the argument. This is why the slide is split per core.
2. **Every cell is its own spatial neighbour** — 23.9% of edges are self-edges, so
   `NeighborhoodToCell` blends autocrine signal into the niche. Undocumented in the vignettes.
3. **`k` is unitless, so the physical neighbourhood shrinks as density rises** — 8.1 µm in the
   densest core vs 22.2 µm in the sparsest. High- and low-grade arms are not measured at the
   same scale.
4. **NICHES' `data` slot is not log-space.** It holds raw LR products in both `counts` and
   `data` (`RunNeighborhoodToCell.R:77-81`). Any downstream Seurat call that assumes log1p —
   including the vignette's own `FindVariableFeatures(selection.method = "disp")` — is being fed
   the wrong space. On the imputed tier the products reach 4,933.5 and `exp()` overflows,
   crashing `CalcDispersion` outright (D10). **Never `NormalizeData` a NICHES assay again.**
5. **Em-dash, not hyphen.** Mechanisms are `LIGAND—RECEPTOR` with U+2014. Every `FetchData` /
   `FeaturePlot` lookup must use it.
6. **Seurat rewrites underscores in feature names.** `CreateSeuratObject` warns
   `Feature names cannot have underscores ('_'), replacing with dashes ('-')`, so a multi-subunit
   mechanism is stored as `TGFB1—TGFBR1-TGFBR2`. Both requested LRs are single-subunit and unaffected.
7. **`species` is silently ignored** when `LR.database = "custom"` (`LoadCustom.R:12`).
8. **`BisRNA` is declared in `Imports` but never called** anywhere in `R/`, and it is archived on
   CRAN — install from the archive purely to satisfy `R CMD INSTALL`.
9. **`SeuratWrappers` 0.4.0 will not install** because of its `Banksy` dependency chain (D3);
   `RunALRA` itself needs none of it.
10. **ALRA picks its rank per object.** Across our 13 cores `RunALRA` chose **k = 20–56**, and
    higher rank means less smoothing: core6 (k = 56) ends up with **611 / 1,088** mechanisms
    detected while every comparable core reaches ~1,000. Since cores are fit independently,
    **imputation strength is not comparable across cores** — a confound layered on top of the
    density confound below.
11. **No memory cliff, but a wall-clock one.** Peak RSS was well below the naive estimate
    (14.5 GB noimpute / 19.8 GB alra on the 26,456-cell core, vs ~26 GB projected) because R
    reclaims the distance matrix before the scoring stage.

**⚠️ Two confounds specific to running NICHES on this TMA**

Both are quantified in `summary_confound.csv` / `summary_edge_audit.csv`, and both push in the
**same direction as grade**, so a naive reading of the high-vs-low contrast will overstate it.

1. **Detection density tracks core cellularity, and high-grade cores are 3.39× larger.**
   High-grade cores average **11,428** cells vs **3,366** for low-grade. Interaction density
   correlates with core size at **r = 0.881** (alra) / **0.410** (noimpute). The consequence is
   visible in `mechanism_detection.csv`: for essentially *every* mechanism, `frac_high` is about
   twice `frac_low` (NCAM1—NCAM1 55.3% vs 39.0%, DLL3—NOTCH1 32.5% vs 14.7%, …). Some of the
   "high grade signals more" signal is simply "high-grade cores have more cells".
2. **The neighbourhood is physically smaller in the denser (high-grade) cores** — 8.1 µm vs
   22.2 µm median edge length — because `k` is fixed and unitless.

This is the same class of density–grade confound already recorded for CytoSignal. The
**cell-type** contrasts (`markers_ALL_by_ReceivingType.csv`) are far less exposed to it, since
they compare populations *within* the same cores.

**Deviations from the tutorial** — full table with justifications in
`scripts/comparators/niches/DEVIATIONS.md`.

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| Unit of a run | one section | one **TMA core**, 13 runs, merged | dense O(N²) edgelist = ~80 GB whole-slide; cores are physically disjoint anyway |
| Imputation | vignettes 01/07 use ALRA; 04 does not | **both**, as parallel sub-runs | requested; and it is decisive — mean mechanisms detected per core **526 vs 966** of 1,088, and the whole grade differential is **null without imputation** (6 vs 229 mechanisms reach `min.pct = 0.25`) |
| Source of `RunALRA` | `library(SeuratWrappers)` | authors' **unmodified** `alra.R` + `internal.R` @ `8df8343`, vendored | SeuratWrappers 0.4.0 won't install (Banksy's Bioc chain fails, 16 pkgs); `RunALRA` needs none of it |
| LR database | `"fantom5"` | `"custom"` = CellChatDB v2.0 human | the `cellchatdb2` tier; format already matches `LoadCustom`, no re-encoding |
| `meta.data.to.map` | vignette 07 omits (all cols) | names 6 columns | guarantees x/y/grade carry through for stages C8–C12 |
| ReceivingType loop | one hand-picked population | **all 9** | no a-priori population; hand-picking would be a silent choice |
| `SpatialFeaturePlot` | Visium `@images` | equivalent `ggplot` scatter | Xenium h5ad has bare coordinates, no image object to dispatch on |
| `nFeature > 5` filter | vignette 04, on CellToCell | CellToCellSpatial only | vignette 07 applies none to NeighborhoodToCell |
| `selection.method = "disp"` | vignettes 01/04/07 | `"disp"` on noimpute; **`"vst"` on alra** | the vignette's own call **crashes** on imputed data — NICHES stores raw LR products in `data`, Seurat's `disp` exponentiates them, and max product 4,933.5 overflows `exp()` |
| `JoinLayers` after `merge` | n/a (pre-v5 vignettes) | added | Seurat v5 keeps one layer per merged object |

**Runs on our data**

| Dataset | Tier | Sub-run | Status | Key numbers | Output path |
|---|---|---|---|---|---|
| GBM | `cellchatdb2` | `noimpute` | ✅ 13/13 cores, 3.9 min | 100,197 cells → **419,299 edges** (319,102 real + 100,197 self); **1,088** mechanisms; **526 detected/core on average** (112–702); mean niche density **0.73%**; peak RSS 14.5 GB | `results/comparators/niches/GBM/cellchatdb2/noimpute/` |
| GBM | `cellchatdb2` | `alra` | ✅ 13/13 cores, 9.2 min | same 419,299 edges; **966 detected/core on average** (529–1,088); mean niche density **9.81%**; ALRA rank k = 20–56 per core; peak RSS 19.8 GB | `results/comparators/niches/GBM/cellchatdb2/alra/` |
| GBM | `cellchatdb2` | `noimpute` differential | ✅ **grade contrast null** | only **6 / 1,088** mechanisms reach `min.pct = 0.25`; **zero grade markers** — globally and in all 9 niches. The cell-type contrast is not quite empty: **1** marker survives (`CNTN2—CNTN2` for Glial-Neuronal, `myAUC` 0.795) | `.../noimpute/_analysis/NeighborhoodToCell/` |
| GBM | `cellchatdb2` | `alra` differential | ✅ | **229 / 1,088** mechanisms reach `min.pct = 0.25`; **180** cell-type markers, **28** global grade markers, 3–98 per-niche grade markers (Glial-Neuronal: 0) | `.../alra/_analysis/NeighborhoodToCell/` |
| GBM | `cellchatdb2` | `noimpute` CellToCellSpatial | ✅ | vignette-04 `nFeature > 5` filter keeps **89,648 / 419,299** edges (21.4%); **13 / 1,088** mechanisms reach `min.pct = 0.25`; **53** VectorType markers, **0** grade markers; top-10 VectorTypes analysed | `.../noimpute/_analysis/CellToCellSpatial/` |
| GBM | `cellchatdb2` | `alra` CellToCellSpatial | ✅ | filter keeps **381,690 / 419,299** edges (91.0%); **116 / 1,088** mechanisms reach `min.pct = 0.25`; **215** VectorType markers, **1** global grade marker; `disp` overflowed again (max product **9,693.1**) → `vst` | `.../alra/_analysis/CellToCellSpatial/` |
| GBM | `cellchatdb2` | `default` (FANTOM5) | ❌ not run | this pass is `cellchatdb2`-only, matching stLearn / SpatialDM / COMMOT / LIANA+ | — |
| LUAD | — | — | ❌ not run | not attempted this pass | — |

**Requested LRIs — where NICHES puts the ALARMIST motif-1 mGAM loop (GBM, `cellchatdb2`)**

Both are in CellChatDB v2.0 human and all four genes are on the Xenium 5K panel, so both survive
`FilterGroundTruth` and are scored **in all 13 cores in both sub-runs** — nothing is missing.

| LR | Sub-run | Detection | Rank by detection | Verdict |
|---|---|---|---|---|
| **GRN—SORT1** | `alra` | **69.2%** of cells (mean 58.1% per core) | mean rank **22** of 1,088; top-10 in 6/13 cores | **strongly recovered** |
| **GRN—SORT1** | `noimpute` | 8.5% of cells | mean rank 28 | recovered but sparse |
| **ANXA1—FPR1** | `alra` | 12.1% of cells | mean rank 307 | present, low-ranked overall — **but see below** |
| **ANXA1—FPR1** | `noimpute` | 1.2% of cells | mean rank 127 | barely detectable |

Two results worth stating plainly, both from the `alra` sub-run:

* **`ANXA1—FPR1` is the single best marker of the mGAM niche** in the vignette-01 cell-type
  contrast — `myAUC = 0.866`, `avg_log2FC = +7.60` for `ReceivingType == mGAM`
  (`markers_ALL_by_ReceivingType.csv`). Its globally low rank is a *prevalence* statement
  (FPR1 is a rare myeloid receptor); its **specificity** to mGAM is the strongest of any
  mechanism NICHES scored. This is an independent corroboration of the ANXA1→FPR1 arm.
* **`GRN—SORT1` is a high-grade marker in 5 of the 9 niches** — AC-like (`myAUC` 0.787),
  MES-like (0.806), NPC-like (0.763), OPC-like (0.747), Vascular (0.708), all with positive
  `avg_log2FC`, i.e. up in high grade. It is also a cell-type marker of the Vascular niche
  (0.739). Read against confound #1 above, the *direction* is consistent but the effect size
  is inflated by core cellularity.

One cross-method concordance: **`JAM3—F11R` is the top low-grade mGAM marker** here
(`myAUC = 0.291`, i.e. down in high grade), matching CytoSignal's `JAM3–F11R` (log-FC −1.10),
one of only four FDR-significant interactions in its grade test.

**Methods paragraph** — Cell–cell signalling was inferred with NICHES v1.2.4 (Raredon et al.)
in R 4.3.3 / Seurat 5.3.0. Each of the 13 TMA cores was processed independently: raw counts were
normalized with `Seurat::NormalizeData` (`LogNormalize`, scale factor 1e4) and, in a parallel
sub-run, imputed with `RunALRA` at its default automatic rank. `NICHES::RunNICHES` was then
applied per core with `LR.database = "custom"` supplying CellChatDB v2.0 human (3,218 ligand–
receptor pairs, of which 1,088 had every subunit represented on the 5,119-gene Xenium panel),
`cell_types = "cell_type"`, spatial coordinates in microns, and the package-default mutual
k-nearest-neighbour graph (`k = 4`, `blend = "mean"`), requesting the `CellToCellSpatial` and
`NeighborhoodToCell` organizations. Per-core outputs were tagged with tumour grade, merged, and
embedded with `ScaleData`, `FindVariableFeatures`, `RunPCA` (100 PCs) and `RunUMAP` (50 PCs).
Differential signalling was tested with `Seurat::FindAllMarkers` (`test.use = "roc"`,
`min.pct = 0.25`, `only.pos = TRUE`) both across the nine receiving cell types and, within each,
between high- and low-grade cores.

---

## CellChat — R, v2.2.0.9001, env `comp-cellchat`

Local source `/Users/jiayifan/tansey_lab/CellChat` @ `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`
(2026-03-04). R 4.3.3, Seurat 5.1.0, presto 1.0.0, NMF 0.27. Code
`scripts/comparators/cellchat/`, contract `NOTES.md`.

Our data is spatial, single-cell resolution, multi-section, two conditions, and no single
vignette covers that, so the contract is assembled from five in a stated precedence order:
`CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.Rmd` (**V2, primary** — it is the
one that matches our shape: many sections in one object via `meta$samples`),
`CellChat_analysis_of_spatial_transcriptomics_data.Rmd` (**V1**),
`FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.Rmd` (**VF** — the Xenium
`spatial.factors`), `Comparison_analysis_of_multiple_datasets.Rmd` (**VC** — high vs low grade)
and `CellChat-vignette.Rmd` (**VB** — the downstream plot inventory both spatial vignettes defer
to).

### Core algorithm

For each L-R pair and each ordered pair of **cell groups** (i → j), CellChat computes a Hill
function of the product of group-average ligand expression in i and group-average receptor
expression in j: `P = L·R / (Kh^n + L·R)` with `Kh = 0.5`, `n = 1`. L and R are computed over
`data/max(data)` with a **10% truncated mean** per group; multi-subunit complexes enter as the
**geometric mean** across subunits (`computeExpr_LR`); the receptor term is further multiplied by
co-activation and divided by co-inhibition receptor terms (`modeling.R:122-124`). Significance is
a **permutation test**: cell-group labels are shuffled `nboot = 100` times and
`pval = #{P_boot >= P_obs}/nboot` (`modeling.R:207,304`).

The unit of a CellChat result is therefore a **(sender cell type, receiver cell type, L-R pair)**
triple — not a cell, not a spot, not an edge. `netP` sums the probabilities of all L-R pairs in a
pathway; `aggregateNet` counts significant links (`net$count`) and sums probabilities
(`net$weight`) per cell-type pair.

### Spatial model

Space enters at the **cell-group** level, not the cell level. For each sample k and ordered group
pair (i,j), `computeRegionDistance` (`modeling.R:1194-1228`) takes every cell of group i, finds
its **1-nearest neighbour** in group j (`BiocNeighbors::queryKNN`, Annoy), converts to µm by
`× ratio[k]`, and takes the **10% trimmed mean** → `d.spatial[i,j,k]`. Then

* `adj.spatial[i,j,k] = 1` iff ≥ `k.min = 10` distinct group-j cells lie within
  `interaction.range + tol` (250 + 5 µm);
* `adj.contact[i,j,k] = 1` iff ≥ 10 lie within `contact.range + tol` (10 + 5 µm);
* across samples `d.spatial` is **averaged** and the adjacency matrices are 1 if **any** sample
  says 1 (`modeling.R:1231-1238`), then symmetrised (`adj * t(adj)`);
* group pairs with `adj.spatial == 0` become `NaN` and are **excluded entirely**.

With `distance.use = FALSE` (V2's value, ours) distance is a **hard filter only**; with
`distance.use = TRUE` (V1) it additionally down-weights by `1/(d × scale.distance)`.

**Consequence for this benchmark: CellChat produces no per-cell and no per-spot output at all.**
Its finest spatial statement is "these two cell types have ≥10 mutual neighbours within 250 µm
somewhere in this section" — coarser than every other method here, and the reason its counts are
not comparable to CytoSignal's per-cell or stLearn's per-spot tallies (see the warning at the top
of this file).

### LR database — and why the two tiers cannot differ by resource

The bundled **`CellChatDB.human` *is* CellChatDB v2** (3,233 rows, 338 complexes, 32 cofactors).
`audit_db_equivalence.R` re-derives the repo CSV's flattening from it using the mapping CLAUDE.md
documents (ligand/receptor ← `complex` subunits joined with `_`) and diffs the key sets:

| | bundled `CellChatDB.human` | `data/LRdatabase/CellChatDBv2.0.human.csv` |
|---|---|---|
| rows | 3,233 | 3,233 |
| unique ligand\|receptor keys | 3,218 | 3,218 |
| Secreted / Non-protein / Cell-Cell Contact / ECM-Receptor | 1280 / 994 / 535 / 424 | 1280 / 994 / 535 / 424 |

**Shared keys 3,218 — only-bundled 0, only-repo 0, Jaccard 1.0000.** The 3,233 vs 3,218 gap is 15
keys listed twice under two annotations; the three `POMC|OPR*` pairs appear as both Secreted and
Non-protein Signaling, and the merge cross-products them, which is the entire content of the "6
signaling_type disagreements" the audit reports. Both sides carry both rows.

So the tiers differ by **annotation scope**, the only knob the vignettes expose:

* **`default`** — `subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")`, the
  literal call in V1:127 and V2:131 → **1,280 interactions, 158 pathways**.
* **`cellchatdb2`** — `subsetDB(CellChatDB)`, the alternative commented in on the next line of the
  same vignette ("use all CellChatDB except for Non-protein Signaling") → **2,239 interactions,
  252 pathways**, i.e. Secreted + ECM-Receptor + Cell-Cell Contact. Verified against
  `database.R:137-141`: with no `search`, `subsetDB` defaults to exactly those three categories.
  This is the scope closest to the resource ALARMIST uses.

Re-importing the flat repo CSV through `updateCellChatDB` was **rejected, not overlooked**: it
would discard 22 columns including `agonist`, `antagonist`, `co_A_receptor` and `co_I_receptor` —
precisely the terms `computeCommunProb` multiplies into `dataRavg` — while adding zero
interactions. Degrading the method to re-supply a database it already ships is not a second tier.

**Widening the DB does not change any individual pair's probability.** The requested LRIs come out
bit-identical between tiers (`max_prob` equal to 16 significant figures, see below), because the
Hill function for a pair depends only on that pair's ligand/receptor expression. The tier changes
*coverage* and *pathway-level aggregation*, nothing else.

### Input

| requirement | detail |
|---|---|
| `data.input` | **normalized** data, genes × cells (V1:49). `adata.X` is already `log1p(CP10K)` — `prepare_gbm_input.py` asserts it and records the implied library size (median exactly 10,000.0). `layers['counts']` is **not** used. |
| `meta` | `labels` ← `obs['cell_type']` (9 types), and a `samples` column ← `obs['tma_id']`. Both factors; unused levels abort `computeCommunProb` (`modeling.R:105-108`). |
| `coordinates` | `obsm['spatial']`, already µm |
| `spatial.factors` | `ratio = 1`, `tol = 5`, **one row per sample in `levels(samples)` order** |

**Object layout.** The vignettes are explicit that several sections of *one* condition go into
**one** object via `samples`, but comparison across conditions needs a **separate object per
condition** (V2:82, VF:118). GBM is therefore two objects: `high` (7 cores — tma_id
1,3,5,8,10,11,13; 79,998 cells) and `low` (6 cores — 2,4,6,9,12,14; 20,199 cells). All 9 cell
types occur in all 13 cores, so `mergeCellChat` and the **functional**-similarity manifold
analysis are both applicable (VC:179).

### Workflow

| # | Call with argument values | Produces |
|---|---|---|
| A5 | `createCellChat(object=data.input, meta, group.by="labels", datatype="spatial", coordinates, spatial.factors)` | CellChat object |
| B3 | `subsetDB(...)` per tier | `cellchat@DB` |
| C1 | `subsetData(cellchat)` | `data.signaling` |
| C2 | `future::plan("multisession", workers=4)` | parallel backend |
| C3 | `identifyOverExpressedGenes(cellchat)` — all defaults | over-expressed genes |
| C4 | `identifyOverExpressedInteractions(cellchat)` — `variable.both=TRUE` (V2 default) | `LR$LRsig` |
| D1 | `computeCommunProb(type="truncatedMean", trim=0.1, distance.use=FALSE, interaction.range=250, scale.distance=NULL, contact.dependent=TRUE, contact.range=10, nboot=100, seed.use=1)` | `net$prob`, `net$pval` |
| D2 | `filterCommunication(min.cells=10)` | filtered net |
| D4 | `computeCommunProbPathway(cellchat)` | `netP$prob` |
| D5 | `aggregateNet(cellchat)` | `net$count`, `net$weight` |
| D6 | `netAnalysis_computeCentrality(slot.name="netP")` | `netP$centr` |
| F1 | `mergeCellChat(list(low=…, high=…), add.names=…)` | merged object |
| F14 | `identifyOverExpressedGenes(group.dataset="datasets", pos.dataset="high", only.pos=FALSE, thresh.pc=0.1, thresh.fc=0.05, thresh.p=0.05)` | cross-condition DEA |
| F15 | `netMappingDEG` → `subsetCommunication(ligand.logFC=±0.05)` | `net.up` / `net.down` |

Everything not listed is the package default: `raw.use = TRUE` (no `projectData` — both spatial
vignettes leave the PPI projection commented out), `population.size = FALSE`, `k.min = 10`,
`do.symmetric = TRUE`, `contact.dependent.forced = FALSE`, `Kh = 0.5`, `n = 1`.

### Data outputs

Per tier dir `results/comparators/cellchat/GBM/<tier>/`:

| file | shape / schema | meaning |
|---|---|---|
| `objects/{low,high}.rds` | CellChat object | full state; replotting re-reads this |
| `objects/merged.rds`, `objects/object_list.rds` | merged + list | the VC comparison object |
| `quant/<c>_net_full.csv.gz` | source, target, interaction_name, prob, pval, pathway, ligand, receptor, annotation | **every** non-zero cell of the prob array, with its p-value — not only the significant subset |
| `quant/<c>_net_significant.csv` | `subsetCommunication()` | significant LR-level links |
| `quant/<c>_netP_significant.csv` | `subsetCommunication(slot.name="netP")` | significant pathway-level links |
| `quant/<c>_netP_prob.csv.gz` | source × target × pathway | pathway probabilities |
| `quant/<c>_net_{count,weight}.csv` | 9 × 9 | aggregated links / summed probability |
| `quant/<c>_centrality.csv` | pathway × cell type | outdeg, indeg, flowbet, info |
| `quant/<c>_d_spatial.csv` | 9 × 9 | mean cell-group centroid distance (µm) actually used |
| `quant/<c>_group_sizes{,_by_sample}.csv` | cell type (× sample) | group sizes |
| `quant/<c>_db_used.csv` | DB rows | the exact interaction table used |
| `quant/<c>_options.json` | — | `object@options$parameter` + run time |
| `quant/<c>_lr_ranked.csv` | interaction_name, summed prob | CellChat's own LR ranking |
| `quant/<c>_selectK_<pattern>_measures.csv` | rank, cophenetic, silhouette | the curve `k` was chosen from |
| `quant/<c>_pattern_<pattern>_{cell,signaling}.csv` | NMF W / H | communication patterns |
| `quant/<c>_embedding_{functional,structural}.csv` | pathway, UMAP1/2, group | manifold learning |
| `quant/summary_by_condition.csv` | condition × 3 | links, strength, pathways |
| `quant/diff_{count,weight}_high_minus_low.csv` | 9 × 9 | differential networks |
| `quant/information_flow_by_pathway.csv` | condition × pathway | what `rankNet` plots |
| `quant/netMappingDEG.csv`, `net_{up,down}_in_high.csv`, `genes_{up,down}.txt` | — | VC Part III DEA |
| `quant/requested_lr_status.csv`, `requested_lr_in_DEA.csv`, `<c>_requested_lr_plot_status.csv` | — | where the requested LRIs landed, including when they landed nowhere |
| `run_manifest.json` | — | DB, every parameter, per-condition counts, wall time, versions, git SHA |

### Image outputs

Every plot goes through one `save_all_formats()` → **png + pdf + svg** (CLAUDE.md rule). Per
condition: `plots/<cond>/aggregate/` (circle count/weight, heatmap count/weight, one circle per
sender), `plots/<cond>/pathways/<pathway>/` (circle, chord, hierarchy, heatmap, LR contribution,
signalingRole network, chord_cell, gene-expression violin, plus circle+chord for every enriched
L-R pair), `plots/<cond>/spatial/core<id>/<pathway>/` (spatial network, and the incoming-weighted
variant — **every pathway × every core**), `plots/<cond>/systems/` (role scatter, outgoing /
incoming / all role heatmaps, selectK curve, river + dot pattern plots, functional and structural
embeddings). Cross-condition: `plots/comparison/` (compareInteractions, diffInteraction
count/weight, differential heatmaps, per-condition circles, role scatter, `signalingChanges`
scatter **for all 9 cell types**, pairwise functional+structural embeddings, rankSimilarity,
rankNet stacked/unstacked, role heatmaps outgoing/incoming/all side by side, bubble
comparison/increased/decreased, bubble+chord for up- and down-regulated pairs, enrichment
wordclouds, split gene-expression violins).

Per-LR plots are kept in **two separate trees**, as the benchmark requires:
`plots/top_lr/<cond>/<LR>/` for the pairs CellChat itself ranks highest by summed communication
probability, and `plots/requested_lr/<cond>/<LR>/` for **GRN_SORT1 and ANXA1_FPR1 regardless of
rank** (circle, chord, bubble, pathway contribution, and per core: gene spatial maps, continuous
LR spatial map, binary LR spatial map).

**Totals: `default` 1,632 plots and `cellchatdb2` 3,666 plots, each written in all three formats**
(4,896 and 10,998 files; 6.3 GB and 13 GB).

**Not produced, with reasons** — written to `plots_not_produced.txt` on every run, never dropped
silently. After the final pass exactly **one** entry remains in each tier:
`netVisual_chord_gene` at *all* sources × *all* targets, which is not renderable — the high-grade
network has **161 distinct ligand/receptor sectors** and circlize fails to allocate them at every
`small.gap` tried (1, 0.5, 0.2, 0.1). The tutorial never draws it at that scope either (VB:381 and
VC:397 both use a single sender), so the covering set is **one chord per sender**, which renders
for all 9 cell types in both conditions, plus per-sender up- and down-regulated chords.
Also deliberately not run: `runCellChatApp` (interactive Shiny, not a file artifact);
`netVisual_embeddingZoomIn` (attempted, skipped when a group is too small); `projectData`/PPI
smoothing (commented out in both spatial vignettes); `mergeInteractions` coarse-cell-type
regrouping (VC option D — our 9 labels are already coarse and no defensible 3-way grouping exists
without asking); and any network view of a requested LR with **no significant link**, which is
itself recorded as a result.

### Multi-sample / differential mode

Native, and used in both forms. **Within** a condition, all cores live in one object as
`meta$samples` levels and CellChat computes cell-group distances per sample then averages them
(`modeling.R:1194-1231`) — so physically separate punches are never treated as neighbours.
**Across** conditions, `mergeCellChat` + the whole of VC: `compareInteractions`,
`netVisual_diffInteraction`, `rankNet(do.stat=TRUE)`, `netAnalysis_signalingChanges_scatter`,
pairwise manifold learning, and the presto-backed cross-condition DEA.

### Gotchas

1. **`spatial.factors` is indexed positionally by `levels(meta$samples)`** — one row per sample in
   level order. Wrong order silently applies the wrong µm conversion. No vignette says so.
2. **`contact.range` is mandatory when `contact.dependent = TRUE`** (`modeling.R:1187-1189`).
3. **`scale.distance` is validated, not defaulted** — with `distance.use = TRUE`, if
   `min(d × scale.distance) < 1` CellChat aborts and prints the value to use
   (`modeling.R:152-156`). Never copy `0.01` from the Visium tutorial.
4. **Unused factor levels abort the run** (`modeling.R:105-108`); `droplevels` is required after
   subsetting to one grade.
5. **presto changes the numbers, not just the speed.** `do.fast = TRUE` is the default and falls
   back to `stats::wilcox.test` **silently** if presto is missing, giving larger logFC while VC's
   `thresh.fc = 0.05` was tuned *for* presto (`utilities.R:434-445`, VC:291). presto 1.0.0 is
   installed and its presence is asserted at startup and recorded in the manifest.
6. **`computeCellDistance` is dense O(N²)** — it needs ~51 GB at 79,998 cells and returns NA.
   Replaced by `BiocNeighbors::findKNN(k = 1)`, which gives the identical quantity in O(N log N).
7. **`sample.use` is mandatory** for `netVisual_aggregate(layout="spatial")` and
   `spatialFeaturePlot` in multi-sample mode.
8. **`netEmbedding` prefers python `umap-learn` via reticulate**; we pass `umap.method = "uwot"`
   so the run does not depend on a reticulate python.
9. **`selectK` cannot be made to work inside a long plotting session.** It calls
   `NMF::nmfEstimateRank` with NMF's default parallel foreach backend and exposes no override; in
   a session that has already run Seurat/ComplexHeatmap work, every run dies with "All the runs
   produced an error", while the *identical call succeeds in a clean session*. We compute the
   measures with `.pbackend = "seq"` and redraw the curve from them.
10. **`identifyCommunicationPatterns` requires `k`** (`analysis.R:385-387` stops on NULL) and the
    vignette picks it by eye. We apply the vignette's stated rule — "the one at which Cophenetic
    and Silhouette begin to drop suddenly" — programmatically to the same measures, and persist
    the measures so the choice is auditable. NMF's rank must also stay below both matrix
    dimensions, so `k.range` is capped at `min(10, n_celltypes-1, n_pathways-1)`.
11. **The `future.rng.onMisuse` warnings are benign.** CellChat draws its permutation matrix once
    under `set.seed(seed.use)` *before* the parallel loop (`modeling.R:206-207`). Two runs with
    identical arguments produced **bit-identical** `net_full`, `net_count` and `net_weight` —
    verified, not assumed.
12. **`netAnalysis_river` needs `ggalluvial` ATTACHED, not just installed.** It builds a
    ggalluvial plot with `stat = "stratum"`, and ggplot2 resolves stats by name off the search
    path, so `ggalluvial::` is not enough — it dies with `Can't find stat called "stratum"`.
    VB:469 loads it explicitly and so must any script calling this. Cost us all 8 river plots on
    the first pass.
13. **`netVisual_chord_gene` has a hard scale ceiling.** At 161 ligand/receptor sectors circlize
    cannot allocate the layout at any `small.gap`. Restricting to one sender — the tutorial's own
    scope — always works.

### Deviations from the tutorial

| item | tutorial | ours | why |
|---|---|---|---|
| `contact.range` | `100` (V1:182, V2:158) | **`10`** | 100 is the 10X Visium spot pitch. VF:39 / V1:60 / `man/computeCommunProb.Rd` pin `10` for single-cell-resolution platforms. **Measured on our data: median nearest-neighbour distance 11.5 µm (low) and 7.45 µm (high)** — so 10 µm is the right scale and 100 µm would have made "contact" meaningless. |
| `spatial.factors` | Visium `ratio = 65/spot_diameter_fullres`, `tol = 32.5` | `ratio = 1`, `tol = 5` | VF:78-83, the FAQ's own **Xenium** row: coordinates already µm, `spot.size = 10` = typical human cell |
| normalization | `GetAssayData(slot="data", assay="SCT")` | `adata.X` as-is | already `log1p(CP10K)`, asserted at export; re-normalizing would log twice |
| `distance.use` | V1 `TRUE` + `scale.distance=0.01`; V2 `FALSE` + `NULL` | **`FALSE`** | V2 governs our data shape; V1's `0.01` is Visium-specific and would abort (gotcha 3) |
| `variable.both` | V1 passes `F`, V2 passes nothing (`TRUE`) | `TRUE` | V2 is the governing vignette and `TRUE` is the package default |
| `umap.method` | not passed (→ `umap-learn`) | `"uwot"` | avoids a reticulate python dependency; the package's own documented alternative |
| `sources.use`/`targets.use` | hardcoded indices (`4`, `5:11`) | all 9 cell types | the indices are specific to VC's 12-cluster skin dataset |
| `pathways.show` | one hand-picked pathway (`"IGF"`, `"EGF"`, `"CXCL"`) | **every** pathway in `netP$pathways` | benchmark invariant: produce every plot the standard workflow can produce |
| `selectK` | called directly | measures recomputed with `.pbackend="seq"`, curve redrawn | gotcha 9 — the stock call cannot complete in this session |
| `k` for patterns | read off the curve by eye | same rule applied programmatically, measures persisted | reproducibility |
| `computeCellDistance` | used as-is | `BiocNeighbors::findKNN(k=1)` | gotcha 6 — dense O(N²) OOMs at 79,998 cells |
| `mergeInteractions` (VC option D) | 12 clusters → 3 coarse types | skipped | our 9 labels are already coarse; no defensible 3-way grouping without asking |
| `netVisual_chord_gene` scope | one sender (`sources.use = 4`) | all-sources attempted, **plus one chord per sender** | all-vs-all is 161 sectors and circlize cannot lay it out (gotcha 13); per-sender is the tutorial's own scope and covers all 9 cell types |
| env: CRAN source | CRAN HEAD | dated snapshot `2024-06-01` | the env is R 4.3.3 and current CRAN sources need ≥ 4.4 (e.g. Deriv 4.2.0 uses `Rf_allocLang`) |
| env: igraph | stock build | `--disable-graphml` + `xml2-config` shadowed | base anaconda's `xml2-config` leaks in and igraph links `libxml2.2.dylib`, absent from the env |
| smoke test | on the tutorial's demo data | on a 819-cell GBM core | **CellChat ships no demo expression data** — only `CellChatDB.*.rda` and `PPI.*.rda`; the vignettes load from the author's local OneDrive paths |

### Runs on our data

| dataset | tier | status | key numbers | output |
|---|---|---|---|---|
| GBM | `default` | ✅ | DB 1,280 interactions / 158 pathways (Secreted only); 525 signaling genes on the 5,119-gene panel. **low**: 20,199 cells / 6 cores, 355 LR tested, **270 significant links**, 23 pathways, 511 s. **high**: 79,998 cells / 7 cores, 382 LR tested, **398 significant links**, 25 pathways, 778 s. | `results/comparators/cellchat/GBM/default/` |
| GBM | `cellchatdb2` | ✅ | DB 2,239 interactions / 252 pathways (Secreted + ECM + Contact); 722 signaling genes on the panel. **low**: 568 LR tested, **927 significant links**, 61 pathways, 1,257 s. **high**: 601 LR tested, **1,191 significant links**, 67 pathways, 1,854 s. | `results/comparators/cellchat/GBM/cellchatdb2/` |
| GBM | db audit | ✅ | Jaccard 1.0000 vs the repo CSV | `results/comparators/cellchat/db_audit/` |

Aggregate comparison (`default`): low 270 links / total strength 3.074 / 23 pathways versus high
398 links / 3.227 / 25 pathways — i.e. high grade has **47% more significant links but only 5%
more total interaction strength**, so the difference is breadth, not intensity. The presto DEA
found **170 up** and **65 down** L-R pairs in high grade; the up list is led by BMP (32 pairs),
COMPLEMENT (26), PDGF (17), TGFb (15) and IGF (15), the down list by COMPLEMENT (12), GAS (12),
GRN (10), FGF (9) and PDGF (7). COMPLEMENT, PDGF and GRN appearing on both sides is the effect VC
warns about at line 282: the DEA is run per cell group, so one pair can be up in one group and
down in another.

### Requested LRIs — CellChat recovers both arms of the ALARMIST motif-1 loop, and only in high grade

Both are `Secreted Signaling`, so both are in **both** tiers, and both survive
`identifyOverExpressedInteractions` in both conditions. Tier makes no difference to the numbers.

| LR | condition | significant cell-type pairs | max prob | ALARMIST direction | that direction's prob |
|---|---|---|---|---|---|
| `GRN_SORT1` (pathway GRN) | low | 15 | 0.04775 | mGAM → MES-like | **absent (n.s.)** |
| `GRN_SORT1` | high | **25** | 0.04434 | mGAM → MES-like | **0.01128, p < 0.001** |
| `ANXA1_FPR1` (pathway ANNEXIN) | low | 8 | 0.02521 | MES-like → mGAM | 0.00140 |
| `ANXA1_FPR1` | high | **10** | 0.02493 | MES-like → mGAM | **0.01624, p < 0.001** (3rd-ranked pair of that LR) |

So the **complete bidirectional mGAM ⇄ MES-like loop is significant only in high grade**: in low
grade the GRN→SORT1 arm is not called at all and the ANXA1→FPR1 arm is ~12× weaker. This is an
independent corroboration of ALARMIST motif 1 being grade-associated, from a method with a
completely different inference target (cell-type pairs, permutation null) and no knowledge of the
motif decomposition.

Two caveats worth stating. First, CellChat's **own** ranking does not put the ALARMIST direction
on top: the strongest GRN_SORT1 pair in both conditions is mGAM → Glial-Neuronal, and the
strongest ANXA1_FPR1 pairs are Vascular → mGAM and Lymphoid → mGAM. The ALARMIST direction is
present and significant but not the maximum. Second, in the cross-condition DEA `GRN_SORT1`
appears in **both** `net.up` and `net.down` (the per-cell-group artefact above), while
`ANXA1_FPR1` appears only in `net.up`.

Also of note: ANXA1_FPR1's receivers in high grade are **exactly `mGAM` and `non-mGAM`** and
nothing else, which is consistent with FPR1 being mGAM-restricted on this panel.

### Reproducing the GBM run

```bash
source scripts/comparators/cellchat/activate_env.sh
bash scripts/comparators/cellchat/run_all_gbm.sh all
```

`prepare_gbm_input.py` (env `bptf`) regenerates the input tree; `run_all_gbm.sh` runs the DB
audit, then both tiers (inference + plots). Env: `env.lock.yml` (conda) +
`r_packages.lock.csv` (249 R packages with versions); `install_env.R` rebuilds the R library.

### Methods paragraph

Cell–cell communication was inferred with CellChat v2.2.0.9001 (Jin et al.) in R 4.3.3. For each
tumour grade a separate CellChat object was created with `createCellChat(datatype = "spatial")`
from log-normalized Xenium expression, the nine annotated cell types as `labels`, and the
constituent TMA cores as `samples`, so that cell-group distances are computed within each core
and averaged across cores; spatial factors were set to the authors' recommended values for
single-cell-resolution platforms (`ratio = 1`, `tol = 5` µm). The ligand–receptor database was
CellChatDB v2 (human), used both as the tutorial's `Secreted Signaling` subset (1,280
interactions) and as the full protein-coding database (`subsetDB`, 2,239 interactions).
Over-expressed genes and interactions were identified with `identifyOverExpressedGenes` and
`identifyOverExpressedInteractions` at package defaults (Wilcoxon test via presto), and
communication probabilities were computed with `computeCommunProb` using a 10% truncated mean
(`type = "truncatedMean"`, `trim = 0.1`), a 250 µm diffusion range as a hard spatial filter
(`distance.use = FALSE`), a 10 µm contact range, and a 100-permutation label-shuffling test
(`seed.use = 1`). Communications were filtered with `filterCommunication(min.cells = 10)`,
aggregated to pathway level with `computeCommunProbPathway` and to cell-type level with
`aggregateNet`, and network centralities computed with `netAnalysis_computeCentrality`. High- and
low-grade objects were combined with `mergeCellChat` and contrasted with `compareInteractions`,
`netVisual_diffInteraction`, `rankNet` (paired Wilcoxon), `netAnalysis_signalingChanges_scatter`
and a cross-condition differential expression analysis (`identifyOverExpressedGenes` with
`group.dataset = "datasets"`, `thresh.fc = 0.05`) mapped back onto the inferred communications
with `netMappingDEG`.
