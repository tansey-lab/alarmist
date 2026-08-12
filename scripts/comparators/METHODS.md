# Comparator methods — living document

One section per cell–cell communication method benchmarked against ALARMIST, written to the
template in `.claude/skills/comparator-benchmark/SKILL.md`. Purpose: record **what each method
computes, what it eats, and what it emits** — not to force a common unit.

Datasets: **GBM** (`data/xenium_mm_final_cell_id.h5ad`, LGG/GBM TMA, 13 cores, grade high/low)
and **LUAD** (`data/linghua/P{17,21}_{AIS,LUAD}_Xenium.h5ad`, AIS vs LUAD). See CLAUDE.md.
Tiers: `default` = the method's own tutorial LR resource; `cellchatdb2` = CellChatDB v2.0 human.

## How to read this document

**This file states what is true now.** It is not the record of how we got here.

| you want | read |
|---|---|
| what a method computes, its parameters, its numbers, its traps | **here** |
| why we departed from a tutorial, and the number that forced it | `<method>/DEVIATIONS.md` |
| the call-for-call contract against the tutorial | `<method>/NOTES.md` |
| what a claim used to say before it was corrected | `_archive/METHODS.2026-08-12.md` (verbatim 4,370-line predecessor) |

Restructured 2026-08-12. The predecessor had grown to 4,370 lines, roughly half of it the
history of its own corrections — 85 `corrected` / `retracted` / `superseded` markers, one
paragraph duplicated verbatim, and LIANA+ split across eight top-level sections. That history
is preserved in the archive above and, where it explains a methodological choice, in each
method's `DEVIATIONS.md`. **Nothing was deleted without a home.**

Convention for this file going forward: a correction **replaces** the wrong statement rather
than being appended beside it. If knowing the old value changes how someone reads a number on
disk, say so in one clause; otherwise let the archive carry it.

## Status

| Method | Language | Version | Env | GBM `cellchatdb2` | GBM `default` | LUAD |
|---|---|---|---|---|---|---|
| CytoSignal | R | 0.5.1 | `comp-cytosignal` | ✅ full slide, 27 min | ❌ | ⚠️ 2.5 mm crop only |
| stLearn | Python | 1.4.1 | `comp-stlearn` | ✅ 45.4 min | ❌ | ❌ |
| SpatialDM | Python | 0.3.1 | `comp-spatialdm` | ✅ 13 cores, 9.3 min | ✅ 12.9 min | ❌ |
| COMMOT | Python | 0.0.3 | `comp-commot` | ✅ 13/13 cores, 296.1 min | ❌ | ❌ |
| LIANA+ | Python | 1.8.1 | `comp-liana` | ✅ 8 branches | ✅ bivariate + inflow + NMF | ❌ |
| NICHES | R | 1.2.4 | `comp-niches` | ✅ 13 cores × 2 imputation | ❌ | ❌ |
| CellChat | R | 2.2.0.9001 | `comp-cellchat` | ✅ | ✅ | ❌ |

**GBM is complete for every method. LUAD is essentially untouched** — only CytoSignal ever ran
on it, on a crop. The `default` tier exists for 3 of 7. Both are open scope, not defects; see
*Open issues*.

---

## Reading the numbers

### The methods do NOT share a spatial support — counts are not comparable

The GBM data is a **13-core TMA**, and the methods were applied at different spatial
granularity. This is deliberate — `SKILL.md:45-46` forbids harmonising kernel scale — but it
means **any table placing their LR-pair counts side by side is misleading unless normalised.**

| Method | Spatial support | Why |
|---|---|---|
| CytoSignal | **whole TMA**, one run over 100,197 cells | no constraint forcing a split; the grade differential separately builds 13 per-core objects and merges via `mergeCytoSignal` |
| stLearn | **whole TMA**, one 12,562-spot grid | no constraint forcing a split |
| SpatialDM | **13 separate runs** | `spatialdm_global` returns **one** Moran's R per pair per object, so a pooled run yields one tissue-wide number instead of 13; `diff_utils.concat_obj` takes a list of separately fitted objects. Matches the authors' own multi-sample tutorial |
| COMMOT | **13 separate runs** | **forced** — dense N×N distance matrix, 80.3 GB whole-slide |
| NICHES | **13 separate runs** | **forced** — `ComputeEdgelist` builds a dense N×N matrix plus three more N×N copies; its fast `nn.method` path is commented-out dead code |
| LIANA+ | **whole slide** (bivariate/inflow/MISTy) or **13 punches** (LRIC) | verified safe per branch, not assumed — see the LIANA+ section |
| CellChat | **two objects** (one per grade), cores as `meta$samples` | its own multi-sample design; distances computed per core then averaged |

**No method's primary tutorial prescribes per-core splitting** — all are single contiguous
sections. Only SpatialDM's secondary multi-sample tutorial demonstrates the pattern.

Consequences:

1. Per-core methods report pairs *per core* (SpatialDM 1,133–1,661; COMMOT a uniform 671;
   NICHES 1,088 scored, 112–702 non-imputed / 529–1,088 imputed detected); whole-TMA methods
   report a single number (stLearn 526; CytoSignal 1,088). **These are different quantities.**
   Compare within a method across cores, or compare *fractions* / ranks across methods — never
   raw counts.
2. Per-core analyses have smaller n and lower power per test, and the LR-pair set differs per
   core because each core re-runs its own expression filter. COMMOT is the exception: its
   `--filter-scope global` evaluates the filter once over all 100,197 cells, so all 13 cores
   share one 671-pair denominator and its ranks *are* cross-core comparable.
3. Harmonising would only be possible in the per-core direction, since COMMOT cannot run
   whole-slide at all. Considered and deliberately not done.

### Cross-core contamination audit — measured in each method's own units

The minimum **cell-to-cell** distance between different cores is **222.9 µm**. That bounds
*cell-level* graphs only — **it does not bound a gridded method**, because gridding relocates
cells to bin centres and spot centroids can end up closer than any two real cells.

| Method | Graph unit | Neighbourhood | Cross-core | Impact |
|---|---|---|---|---|
| CytoSignal — diffusion | cell | 200 µm ε-ball | **0 pairs** | none |
| CytoSignal — contact | cell | Delaunay, pruned at 200 µm | **0 edges** | none — but only by a 23 µm margin, see below |
| stLearn | **51.3 µm spot** | `distance = 250` µm | **7 spot pairs** of 344,370 (0.0020%); 10 of 12,562 spots | **40 of 45,065 significant spot-LR calls (0.0888%)** |
| SpatialDM | cell | 135 µm, per core | n/a — split | none |
| COMMOT | cell | 365 µm, per core | n/a — split | none |
| LIANA+ — bivariate / inflow | cell | **28.2 µm** support | **0 pairs** | none (7.9× margin) |
| LIANA+ — **LR-MISTy** | cell | 607 µm nominal, **capped at 100-NN, cap binds for 99.7% of cells** | **2,520 of 10,119,190 edges (0.0249%)**, 232 of 100,190 cells | small but **not zero** |
| LIANA+ — LRIC / cross-PCF | cell | annuli to 225 µm, per punch | n/a — split | none |
| NICHES | cell | **mutual kNN, k=4 — no radius**; measured median **10.1 µm**, p95 26.6 | n/a — split | none |
| CellChat | **cell group** | 250 µm within a core, per `meta$samples` | n/a — per sample | none |

**Two whole-slide runs have genuine cross-core edges: stLearn and LIANA+'s LR-MISTy.** Both are
the *tutorial's own* parameter doing it — stLearn's `distance = 250` exceeds the 222.9 µm floor
(shortest cross-core *spot* pair 205.4 µm; affected core pairs (9,11), (10,12), (11,14); worst-hit
LRs JAM3_JAM3 2/457, C3_C3AR1 1/1,244, MPZ_MPZL1 1/111), and `lrMistyData`'s `bandwidth=200,
cutoff=0.01` gives 607 µm nominal reach, 2.7× the floor (concentrated in core pairs (4,6) 787
edges, (12,11) 575, (5,6) 461). Both are too small to move a ranking and too real to omit.
**LIANA's "0 pairs" row is specific to the 28.2 µm kernel and does not transfer to its other
entry points.**

**CytoSignal is clean only by a 23 µm margin.** Delaunay triangulation is *unbounded* — the raw
triangulation has **716 cross-core edges, the longest 12.3 mm**. `findNNDT` removes them via
`dist <= max.r` with `max.r = r.eps.real / scale.factor = 200 µm`, and the shortest cross-core
DT edge is 223 µm. **Had `r.eps.real` been 250 rather than 200, the contact-dependent slot would
silently have contained cross-patient "touching cell" edges.** Re-check before reusing this on
TMA data with a larger radius.

### ⚠️ The density–grade confound — applies to every cross-grade result in this document

Stated once here; every method's grade result refers back to it.

Core cell density correlates with grade: **r = 0.659, p = 0.014** (point-biserial, n = 13;
Mann-Whitney p = 0.005). Median cells within a 134.6 µm radius is **243 in high-grade cores
(78–532) vs 60 in low-grade (35–118)** — a **4× difference**; high-grade cores average
**11,428** cells against **3,366**. High cellularity is itself a WHO glioma grading criterion,
so this is real biology, not a pipeline artifact — but it means **statistical power is not
constant across grades for any method whose neighbourhood is distance-based**, and for `k`-based
neighbourhoods (NICHES) the *physical* scale differs between the arms as well.

Measured consequences, per method:

| method | what tracks density |
|---|---|
| SpatialDM | corr(`n_neighbors`, n_significant) **r = +0.709, p = 0.007**; the grade difference **disappears** once normalised to the fraction of testable pairs (0.308 vs 0.252, p = 0.628) |
| SpatialDM | per-core mean \|z\| vs neighbourhood size **r = +0.752, p = 0.003** |
| NICHES | interaction density vs core size **r = 0.881** (alra) / 0.410 (noimpute); `frac_high` ≈ 2 × `frac_low` for essentially every mechanism |
| NICHES | median edge length **8.1 µm** (densest core) vs **22.2 µm** (sparsest) — a 2.7× spread, because `k` is unitless |
| COMMOT | magnitudes grow with cells inside the 365 µm radius, so "myeloid is the hub" is partly density |
| CytoSignal | 200 µm ε-ball, same exposure |

**Read cross-grade claims as fractions or ranks, never counts.**

**A second, independent limit: 13 cores from 7 patients**, patient 14007 contributing 4 — so
core-level tests are pseudoreplicated at the patient level. And **4 of the 7 patients contribute
both a high- and a low-grade core, covering 10 of the 13**, so the grade contrast is mostly
*within*-patient; no method's native test uses that pairing.

**The power floor.** A two-sided rank test on **7 high vs 6 low** punches cannot go below
**p = 0.0011655**. Every punch-level null in this document is bounded by it. **Those nulls are
underpowered designs, not evidence of absence.**

---

## CytoSignal — R, v0.5.1, env `comp-cytosignal`

Welch lab (Liu & Wang). Tutorials: `results/comparators/cytosignal/reference_notebook/*.html`
(4 vignettes: main workflow, differential multi-dataset, custom LR DB, container conversion).
There is **no local git clone** — the vendored HTML is the only source of truth.

Code: `scripts/comparators/cytosignal/` — `run_cytosignal.R` (the run), `build_cellchat_db.R`
(the `cellchatdb2` DB), `quant_io.R` (persistence), `run_nebula_grade.R` (the grade test),
`plot_signif_rerun.R` (stock figures, by re-run), `build_grade_2panel.py` (the grade figure),
`activate_env.sh`.

> ⚠️ **This method has no `NOTES.md`, no `DEVIATIONS.md` and no `env.lock.yml`**, although
> `SKILL.md:71` names it the reference implementation for the others. The deviations table
> below is its only contract record. Open issue.

### Core algorithm

A **nonparametric, per-cell permutation test for ligand–receptor activity**. For each
interaction and *each individual cell*: is local L–R co-expression higher than chance?

1. **Neighborhood** — separate definitions for diffusible and contact-dependent interactions.
2. **Imputation** — how much L and R each cell *receives* from neighbors, distance-weighted.
3. **LRscore** — `L × R` per cell, averaged over the Delaunay neighborhood (this is what the
   `_smooth` suffix on every score slot means).
4. **Spatial permutation test** — null built by randomizing **locations**, giving a one-sided
   p-value per (cell × interaction), then a **spatial FDR** accounting for density differences.

Two properties matter. **The unit of inference is the single cell** — output is a cell ×
interaction score matrix plus, per interaction, the set of cells where it is significant.
Cell-type labels are used only for plot colouring and optional NEBULA covariates; they do **not**
enter the scoring. And **diffusion and contact are modelled separately end to end**, with
different kernels and separate result slots — an interaction is classified by the database, not
inferred. Multi-subunit complexes are native (`L1..Ln` / `R1..Rn` columns).

Ranking is either by number of significant cells or by **SPARK-X spatial variability**
(`rankIntrSpatialVar`, the `result.spx` tier).

### Spatial model

| Mode | Neighborhood | Weighting | Slot |
|---|---|---|---|
| Diffusion | ε-ball, `r.eps.real` = **200 µm** (default) | Gaussian on physical distance | `GauEps` → `diffusion-Raw` |
| Contact | Delaunay (immediate neighbors) | `dt.mode = "weight_sum_2"` | `DT` → `contact-Raw` |
| Same-spot | none (raw expression) | — | `Raw-Raw` (for Visium-like multi-cell spots) |

`inferEpsParams(scale.factor, r.eps.real = 200, thresh = 0.001)` converts the physical radius
into kernel parameters. **`scale.factor` = µm per coordinate unit and is external knowledge** —
the single most dangerous parameter here, because a wrong value silently rescales every
neighborhood. The tutorial uses 0.73 (Slide-tags) / 0.72 (Visium). Xenium coordinates are
already µm, so **`scale.factor = 1`** and the ε-ball is literally 200 µm.

`Raw-Raw_smooth` is scored regardless of platform; for single-cell data it is not the intended
readout but is still written out.

### LR database

**`default` — bundled CellPhoneDB v2** (re-sorted by the authors), via
`addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)`. Verified in the installed package:
`g_to_u` 977 genes, `inter.index` 1,396 interactions, **754 diffusion + 109 contact** unique.
The DB works in a **UniProt ID space**; `changeUniprot()` rewrites the gene-symbol matrix into
it and drops genes absent from the DB.

**`cellchatdb2` — CellChatDB v2.0 human**, built by `build_cellchat_db.R` using the package's
own exported `formatLRDB(...)` (the documented custom-DB route). 3,233 CSV rows → **2,683
diffusion + 535 contact**, 1,383 LR genes, **865 on the Xenium 5K panel**. Two decisions:

- **Gene symbols are the protein-ID space** (identity `g_to_u`). Legitimate — CytoSignal only
  requires the ID space to be internally consistent — and it matches the expression matrix.
- **`signaling_type` → interaction type:** `Cell-Cell Contact` → contact; `Secreted Signaling` /
  `ECM-Receptor` / `Non-protein Signaling` → diffusion. Mirrors ALARMIST's juxtacrine/secreted split.

### Input

| Argument | Form | Notes |
|---|---|---|
| `raw.data` | genes × cells `dgCMatrix` | **raw integer counts**, rownames = UPPERCASE symbols |
| `cells.loc` | cells × 2 numeric | colnames lowercase `x`, `y`; rownames = barcodes |
| `clusters` | named `factor` | names must match `colnames(raw.data)` exactly |
| `scale.factor` | scalar | µm per coordinate unit — **not in the data**, supplied by the user |

On-disk exchange format (`<dataset>/input*/`, deliberately *outside* the tier dirs since the
export does not depend on the LR database): `counts.mtx`, `genes.tsv`, `barcodes.tsv`,
`meta.csv` (`cell_id,x,y,celltype`), `provenance.json`. GBM adds `meta_grade.csv`
(`+tma_id,grade,patient`). Converters exist both ways: `SeuratToCS`, `SCEToCS`.

### Workflow

| # | Call (tutorial argument values) | Produces |
|---|---|---|
| 1 | `createCytoSignal(raw.data, cells.loc, clusters)` | `CytoSignal` object |
| 2 | `addIntrDB(cs, g_to_u, db.diff, db.cont, inter.index)` | LR DB attached |
| 3 | `removeLowQuality(cs, counts.thresh = 300, gene.thresh = 50)` | QC |
| 4 | `changeUniprot(cs)` | expression re-indexed into the DB's ID space |
| 5 | `inferEpsParams(cs, scale.factor = <µm/unit>, r.eps.real = 200)` | `eps`, `sigma` |
| 6 | `findNN(cs)` | `GauEps` + `DT` graphs |
| 7 | `imputeLR(cs)` | imputed L and R per cell |
| 8 | `inferIntrScore(cs, perm.size = 1e5)` (seed first) | LRscores + null + p-values |
| 9 | `inferSignif(cs, p.thresh = 0.05, reads.thresh = 100, sig.thresh = 100)` | `result`, `result.hq` |
| 10 | `rankIntrSpatialVar(cs)` | `result.spx` (SPARK-X ranked) |
| 11 | `showIntr(cs, slot.use, signif.use, return.name = TRUE)` | significant interaction list |

Optional `recep.smooth = TRUE` in step 8 adds `diffusion-DT` / `contact-DT` slots for sparse
data; not used.

### Data outputs

Object slots (`@lrscore[[slot]]`), for each of `diffusion-Raw_smooth`, `contact-Raw_smooth`,
`Raw-Raw_smooth` (and non-smooth `*-Raw` variants):

| Slot | Shape | Meaning |
|---|---|---|
| `@score` | cells × interactions | LRscore. **Near-dense** — the dominant memory/disk cost |
| `@res.list$result` | list per interaction | barcodes with p < `p.thresh` |
| `@res.list$result.hq` | list per interaction | + passes `reads.thresh` / `sig.thresh` |
| `@res.list$result.spx` | list per interaction | + spatially variable by SPARK-X (the headline tier) |

Persisted by `quant_io.R` into `<run>/quant/` — names alone are useless for comparison:

| File | Content |
|---|---|
| `score_<slot>.mtx` + `.cells.tsv` + `.intr.tsv` | cells × interactions LRscore (≤100k cells) |
| `score_<slot>.rds` | same as `dgCMatrix` with dimnames (>100k cells) |
| `reslist_<slot>.rds` | full `@res.list` |
| `signif_summary_<slot>.csv` | `interaction_id, name, n_result, n_hq, n_spx` |

Multi-dataset: `@diff.results$<level>` → `interaction, logFC, se, p, padj`. Export routes:
`csToSeurat`, `csToSCE`, `csToH5AD`.

### Image outputs

| Function | Shows | Used? |
|---|---|---|
| `plotCluster` | cell types in space | ✅ `cluster_map.png` both runs |
| `plotSignif` | per-interaction panel: imputed L & R, raw L & R, LRscore, clusters | ✅ GBM top 6 per mode in `plots/signif_<slot>/`, + the two requested LRIs in `plots/requested_diffusion_Raw_smooth/`; ✅ LUAD crop |
| `plotEdge` | 3D sender→receiver edge plot | ❌ not run |
| `plotNebulaVolcano` / `plotNebulaAll` | volcano / summary of the multi-dataset test | ❌ unusable — the split-stage route never populates `multics@diff.results`, which both read |
| `heatmap_GO`, `plotREVIGO` | GO enrichment of signalling-associated DEGs | ❌ not run — downstream interpretation, out of scope |
| via `csToSeurat` / `csToSCE` | `SpatialFeaturePlot` etc. on LRscore layers | ❌ not run |

Everything else in `run_full/plots/` (`ranking.png`, `top_{diffusion,contact}_grid.png`,
`spatial_scores.png`, `lr_panels_*.png`, `comparison_mGAM_vs_cytosignal.png`,
`reconstruct_motif1_from_cytosignal.png`, `motif1_top25_lris_cytosignal.png`,
`grade_comparison_2panel.png`) is **custom ALARMIST-comparison plotting, not stock CytoSignal.**
The six one-off scripts that produced them were archived 2026-08-12 to
`scripts/comparators/_archive/cytosignal/` — superseded by `figure6_*` / `figure7_*`, but kept
because their outputs are still on disk.

### Multi-sample / differential mode

Native and documented: build one object per sample (each through `findNN` + `imputeLR` — the
merge needs the DT imputation), then

```r
multics <- mergeCytoSignal(objList, metadata = dataset.meta, name.by = "sample")
multics <- runNEBULA(multics, covariates = c('clusters', 'age'), cpc_thresh = 0.001, ncore = 4)
```

`runNEBULA` fits a **negative-binomial mixed model** per interaction, dataset as the random
effect, total counts as offset. Tutorial signature default `cpc_thresh = 0.005`; the vignette
body uses **0.001** (lower keeps *more* interactions) — we followed the vignette.

For GBM we used **`tma_id` (13 cores) as the sample unit and `grade` as the covariate**.

### Gotchas

- **`removeLowQuality` default `counts.thresh = 300` drops ~half of a Xenium 5K panel.**
  Targeted panels have far lower per-cell totals than whole-transcriptome data.
- **`showIntr(return.name = TRUE)` errors on a custom DB.** `formatLRDB` writes only 3 columns
  to `intr.index`; `getLigandNames`/`getReceptorNames` read columns 4–5. `build_cellchat_db.R`
  pads them so naming falls back to `partner_a/b`.
- **`formatLRDB` returns `db.diff`/`db.cont` in a different column order than the bundled DB**,
  and `getIntrValue` inside `plotSignif` reads **by position** — a raw custom DB swaps the
  ligand/receptor labels in plots. `build_cellchat_db.R` reorders. Scoring reads by name and is
  unaffected.
- **`@score` is effectively dense.** 89k cells × 919 interactions ≈ 1.9 GB as `.mtx`.
  `purgeBeforeSave()` keeps `@score`/`@res.list` and clears only imputation/raw slots — it is
  not a small-file escape hatch.
- **`nosave` runs cannot be replotted, only recomputed.** Every plotting function needs
  `@score` + `@imputation` + `@counts` + `@cells.loc`, none of which survive a disk-safe run.
  Skipping a plot at run time costs a **full pipeline re-run** later (~17 min for GBM), not a
  cheap replot. **Plan the plot list up front.**
- **`perm.size` is floored at `n_cells`.** At 498k cells the null matrix is enormous; this is
  what OOMs the full LUAD section (~57 GB peak on a 36 GB box).
- **`nebula` will not build in the conda env's R 4.3.3** (C++ vs Eigen 3.4/lgamma; CRAN needs
  R ≥ 4.4), while `cytosignal` will not link in system R 4.4.2 (Fortran/BLAS).
- **`SPARK` is not auto-installed** despite the `Remotes` field — `rankIntrSpatialVar` fails
  without it, and `plotSignif(raster = TRUE)` needs `scattermore`.
- Degraded conda: `conda activate` / `conda run` fall through to system R. Always
  `source scripts/comparators/cytosignal/activate_env.sh`.

### Deviations from the tutorial

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `removeLowQuality` | 300 / 50 | **100 / 20** | 300 drops ~half the cells on a 5K targeted panel. GBM 100,197 → 89,035. Forced by panel size |
| `scale.factor` | 0.73 (Slide-tags) | **1** | Xenium coordinates are already µm. Required for 200 µm to mean 200 µm |
| `numCores` | 1 | 4 | performance only |
| `runNEBULA` | single call on the merged object | **2-stage script** (`run_nebula_grade.R`): stage 1 in conda R builds/merges the 13 core objects and extracts the model inputs via `cytosignal:::.setup.model(merged, "grade")`; stage 2 re-invokes itself in **system R 4.4.2** to call `nebula::nebula` | `nebula` and `cytosignal` cannot live in the same R. Uses CytoSignal's own model-setup internals so the fitted model is the one `runNEBULA` would fit. **Depends on a private `:::` function — the one place we touch internals** |
| NEBULA covariates | `c('clusters','age')` | **`c('grade')` only** | cell-type composition differs systematically by grade; including `clusters` would absorb the effect of interest. A real modelling choice, not a technicality |
| `plotNebulaVolcano` / `plotNebulaAll` | used | replaced by `build_grade_2panel.py` | the split-stage route never populates `multics@diff.results` |
| p-value adjustment | internal to `runNEBULA` | **BH within each mode** (diff and cont separately) | matches `runNEBULA` internals; pooling the modes would be wrong |
| `cpc_thresh` | 0.005 (signature) / 0.001 (vignette) | **0.001** | followed the vignette |
| `plotSignif` for GBM | top 5 per mode | **top 6 per mode**, via a full re-run | the run was `nosave`. 6 to match the LUAD crop. The re-run reproduces the original exactly — see below |
| requested LRIs | n/a | **GRN→SORT1 and ANXA1→FPR1 additionally plotted** into `plots/requested_<slot>/` | `SKILL.md:56-61`; separate directory so they are never confused with the method's own ranking |
| `Raw-Raw_smooth` plots | tutorial plots whichever slot | **skipped** | that slot is the multi-cell-spot readout, not the intended output for single-cell Xenium. Scores still persisted |
| GO / REVIGO / `inferIntrDEG` | in the tutorial | not run | downstream interpretation, out of scope |

### Runs on our data

Layout under `results/comparators/cytosignal/`:

```
GBM/   input_full/                  input export — outside the tiers (DB-independent)
       default/     STATUS.md       not run (bundled CellPhoneDB v2)
       cellchatdb2/ run_full/       the GBM run
                    nebula_grade/   the high-vs-low grade test
LUAD/  input/  input_full/          2.5 mm crop and full-section exports
       default/     run_crop_2p5mm/ the LUAD run (crop only)
       cellchatdb2/ STATUS.md       blocked, needs >= 64 GB
bundle_bignode/                     portable bundle for the blocked LUAD run
reference_notebook/                 vendored tutorials    cellchat_db_human.rds   shared DB
```

The run-name subdirectory sits inside the tier because it carries information the tier does not
(which section, full vs crop, which analysis). `bundle_bignode/` holds `.R`/`.py` inside
`results/` — a sanctioned exception to the code/output separation, being a self-contained bundle
meant to be copied to another machine.

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `GBM/cellchatdb2/run_full/` | ✅ | 100,197 → **89,035** cells after QC; **919 diffusion / 169 contact** scored; **895 / 166** with ≥1 significant cell (151 in `Raw-Raw_smooth`). Top diffusion by significant-cell count: WNT3–FZD*/LRP6 (~44k cells); top contact: OCLN–OCLN, CADM3–CADM3, PTPRM–PTPRM |
| GBM | `default` | `GBM/default/STATUS.md` | ❌ not run | bundled CellPhoneDB v2 tier; command in the STATUS file |
| GBM | grade test | `GBM/cellchatdb2/nebula_grade/` | ✅ | 13 cores; **659 interactions** tested (535 diff / 124 cont); 147 raw p<0.05; **only 4 FDR-significant**, all generic junction/adhesion: CDH2–CDH2 (+0.53), NECTIN3–NECTIN1 (−0.63), F11R–F11R (−0.88), JAM3–F11R (−1.10). Of motif 1's top-100 LRIs (95 testable) only **1** (JAM3–F11R) is significant. Read against the density confound |
| LUAD | `default` | `LUAD/default/run_crop_2p5mm/` | ⚠️ crop only | P21_LUAD, 2.5 mm window, 28,596 cells. **277 diffusion / 44 contact** significant (`result.spx`). Top diffusion: TGFB1/TGFB3–TGFβR1/R2, PDGFB–PDGFR; top contact: αMβ2/αLβ2/αXβ2–ICAM1, JAG1–NOTCH1, DLL4–NOTCH4 |
| LUAD | `cellchatdb2` | `LUAD/cellchatdb2/STATUS.md` → `bundle_bignode/` | ❌ blocked | full P21 (560,183 → 498,422 cells) OOMs at ~57 GB. Bundle prepared for a ≥64 GB node |
| LUAD | AIS vs LUAD | — | ❌ not run | only P21_LUAD was ever touched |

### Requested LRIs

Both are in the DB, both are called significant, and both are **diffusion-only** — consistent
with GRN and ANXA1 being secreted, and with CytoSignal classifying by database annotation rather
than inferring the mode.

| LRI | ID | slot | rank (of 895 signif.) | significant cells | tiers |
|---|---|---|---|---|---|
| GRN → SORT1 | `CCI-01109` | `diffusion-Raw_smooth` | **66** | 27,630 / 89,035 (31%) | result = hq = spx |
| ANXA1 → FPR1 | `CCI-01088` | `diffusion-Raw_smooth` | **255** | 15,292 / 89,035 (17%) | result = hq = spx |
| both | — | `contact-`, `Raw-Raw_smooth` | — | not significant | — |

`result = hq = spx` means each passes the p-value threshold, the read/bead QC **and** SPARK-X
spatial variability — CytoSignal considers both spatially structured, not diffuse background.
Neither is near the top of its own ranking, which is dominated by WNT3–FZD*/LRP6 (diffusion) and
homotypic junction pairs (contact). Figures:
`plots/requested_diffusion_Raw_smooth/Rank_{66_GRN-SORT1,253_ANXA1-FPR1}.png`.

### Reproducing the GBM run

`plot_signif_rerun.R` rebuilds the object from `input_full/` with `run_cytosignal.R`'s exact
parameters and **checks itself against the stored `quant/` before plotting**
(`plots/reproduction_check.csv`) — every interaction in every slot matched exactly:

| slot | interactions | `n_hq` exact matches | corr | `n_spx` exact matches | corr |
|---|---|---|---|---|---|
| `diffusion-Raw_smooth` | 895 | 895 | 1.00 | 895 | 1.00 |
| `contact-Raw_smooth` | 166 | 166 | 1.00 | 166 | 1.00 |
| `Raw-Raw_smooth` | 151 | 151 | 1.00 | 151 | 1.00 |

So CytoSignal is **fully deterministic** here under `set.seed(42)` even at `numCores = 4` — a
lost object is always recoverable at ~27 min (whole pipeline) or ~12 min (one slot). Cost: no
checkpoint means recompute, never replot.

### Methods paragraph

> For CytoSignal (v0.5.1), we analyzed the data according to the default workflow, which
> consists of (1) defining spatial neighborhoods, (2) imputing ligand and receptor expression,
> (3) calculating LRscores with a spatial permutation test and (4) identifying significant
> interactions. We used the functions `inferEpsParams`, `findNN`, `imputeLR`, `inferIntrScore`,
> `inferSignif` and `rankIntrSpatialVar`, retaining the authors' default parameters
> (`r.eps.real = 200` µm, `perm.size = 1e5`, `p.thresh = 0.05`, `reads.thresh = 100`,
> `sig.thresh = 100`). Because Xenium coordinates are reported in microns, we set
> `scale.factor = 1`; because the 5,000-plex targeted panel yields lower per-cell counts than
> the whole-transcriptome data used in the tutorial, we relaxed `removeLowQuality` to
> `counts.thresh = 100` and `gene.thresh = 20`. Interactions were reported at the `result.spx`
> level, that is, significant, quality-controlled and spatially variable by SPARK-X. For the
> two-condition comparison, we built one CytoSignal object per TMA core, merged them with
> `mergeCytoSignal` and tested for differential interaction usage between high- and low-grade
> cores with a negative-binomial mixed model as implemented in `runNEBULA` (`cpc_thresh = 0.001`),
> treating each core as an independent sample and applying Benjamini–Hochberg correction
> separately within diffusion- and contact-dependent interactions.

---

## stLearn — Python, v1.4.1, env `comp-stlearn`

Tutorial: the **Xenium-specific** CCI vignette (`cell_cell_interaction_xenium`); the generic
`cell_cell_interaction` one targets Visium. Contract: `stlearn/NOTES.md`. Deviations:
`stlearn/DEVIATIONS.md`. Code: `run_stlearn.py` (the run), `build_cellchat_lrs.py` (DB
conversion), `plot_stlearn_tutorial.py` (the vignette figure set), `plot_stlearn_full.py`
(beyond-vignette figures), `export_stlearn_quant.py` (stage-2 export), `build_stlearn_report.py`.

### Core algorithm

Two stages answering different questions.

**Stage 1 — `run()`: where is an LR pair co-expressed more than chance?** Cells are aggregated
into a **regular grid of spots**; a spot's LR score is L–R co-expression across that spot and its
neighbours within a physical radius. The null samples `n_pairs` **random gene pairs matched on
expression level**, scores them identically, and asks where the real pair exceeds its own
background. Per-spot p-value per LR pair, BH-corrected. **The unit is the spot.**

**Stage 2 — `run_cci()`: which cell types sit in those hotspots?** Restricted to significant
spots, it counts cell-type→cell-type edges across the neighbourhood graph, then permutes the
**cell-type labels/proportions** (not expression) to test over-representation. With
`spot_mixtures=True` a spot may count as several cell types, using the per-spot proportions
`grid()` stored.

So stLearn never scores an individual cell, and cell-type identity enters only in stage 2 — the
LR statistics themselves are cell-type agnostic. A different target of inference from
CytoSignal, which tests every cell.

### Spatial model

| | |
|---|---|
| aggregation | regular grid, `n_row` × `n_col` bins over the bounding box; **empty spots dropped** |
| neighbourhood | `distance=250` physical units (cKDTree radius), **independent of grid resolution** |
| ours | 12,562 occupied spots, 51.3 × 51.3 µm, 8.0 cells/spot, **61 median neighbours**, 0 isolated |

The grid is an aggregation knob (the tutorial calls `n_=125` a resolution/compute trade-off);
`distance` is the actual signalling range.

### LR database

Default **connectomeDB2020_lit** (from NATMI), 2,293 pairs, via `st.tl.cci.load_lrs`.

**stLearn cannot represent multi-subunit complexes.** Its LR format is the string
`"LIGAND_RECEPTOR"`, so `_` is the L/R separator and cannot also mean "subunit of". Converting
CellChatDB v2 drops every complex row: **3,233 → 1,371 pairs (57.5% lost)**, of which **527 are
on the GBM 5K panel** and **526 survived stLearn's own expression filter**. Compare CytoSignal,
complex-aware, which tested **1,088** panel interactions on the same data.

**This is not cosmetic.** All 40 `WNT3_*` rows in CellChatDB v2 have complex receptors
(`FZD*_LRP5/6`), so **stLearn cannot test a single WNT3 interaction** — while WNT3–FZD*/LRP6 was
CytoSignal's *entire top-6* diffusion ranking on this dataset. Part of the two methods'
disagreement is structural, not biological.

### Input

| requirement | detail |
|---|---|
| `X` | **raw counts, never log1p'd**. The permutation null picks background genes of similar expression, so log-shrinking genes together breaks it. Ours from `layers['counts']` |
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
| `uns['lr_summary']` | LR × **6** | stage 1: `n_spots`, `n_spots_sig`, `n_spots_sig_pval` — plus the three `*_<label>` columns `run_cci()` appends, which are the stage-2 ranking |
| `obsm['lr_scores']` | spots × LR | raw LR co-expression score |
| `obsm['p_vals']`, `['p_adjs']`, `['-log10(p_adjs)']` | spots × LR | significance |
| `obsm['lr_sig_scores']` | spots × LR | score masked to significant spots |
| `obsm['spot_neighbours']`, `['spot_neigh_bcs']` | spots | neighbourhood graph |
| `uns['<label>']` | spots × types | per-spot cell-type proportions |
| `uns['lr_cci_<label>']`, `['lr_cci_raw_<label>']` | types × types | interaction counts pooled over LRs |
| `uns['per_lr_cci_<label>']` | dict of types × types | **one count matrix per LR** (526) |
| `uns['per_lr_cci_pvals_<label>']` | dict of types × types | **one p-value matrix per LR** — the significance of stage 2 |
| `uns['per_lr_cci_raw_<label>']` | dict of types × types | uncorrected per-LR counts |

`run_stlearn.py` persists most of these under `<out>/data/` as gzipped CSVs + `grid.h5ad`
(771 MB). **Two gaps it left**, closed by `export_stlearn_quant.py` (read-only replay from
`grid.h5ad`, nothing recomputed): it looped only over `["lr_cci_<label>", "per_lr_cci_<label>"]`,
so the stage-2 **p-values** and both raw-count sets never reached disk; and it wrote
`lr_summary.csv` *before* `run_cci()`, so the on-disk copy had 3 of the object's 6 columns. The
rewrite was verified to leave all 3 pre-existing columns byte-identical.

### Image outputs

stLearn ships two CCI vignettes and **they do not have the same figure set**. Provenance was
verified by grepping both, because several figures we had treated as "the standard workflow" are
generic-vignette-only, and three are in neither.

| Call | Shows | Xenium vignette | generic Visium | Ours |
|---|---|---|---|---|
| `st.pl.cluster_plot` | cell types, grid vs single-cell | ✅ | ✅ | ✅ |
| `st.pl.feat_plot` | per-spot cell-type proportion | ✅ | ❌ | ✅ all 9 types |
| `st.pl.gene_plot` | gene expression, grid vs single-cell | ✅ | ❌ | ✅ |
| `st.pl.lr_summary` | LR ranking bar chart | ✅ ×2 | ✅ | ✅ |
| `st.pl.lr_result_plot` | **per-LR spatial map** | ✅ top-1 | ✅ | ✅ |
| `st.pl.cci_check` | **diagnostic: LR significance vs cell-type frequency** | ✅ | ✅ | ✅ |
| `st.pl.ccinet_plot` | cell-type interaction network, shared node layout | ✅ | ✅ | ✅ |
| `st.pl.lr_chord_plot` | chord diagram | ✅ | ✅ | ✅ |
| `st.pl.lr_diagnostics` | expression-vs-significance diagnostic | ❌ | ✅ | `plots/` only |
| `st.pl.lr_n_spots` | spots per LR | ❌ | ✅ | `plots/` only |
| `st.pl.cci_map`, `st.pl.lr_cci_map` | cell-type interaction heatmaps | ❌ | ✅ | `plots/` only |
| `st.pl.lr_plot` | per-spot ligand/receptor detail | ❌ | ✅ ×8 | `plots_full/` only |
| `st.pl.lr_go` | GO enrichment of top LR genes | ❌ | ✅ | ❌ **unavailable** — `run_lr_go(r_path=...)` needs R with `clusterProfiler`; no comp env has it |
| `st.pl.het_plot`, `grid_plot`, `deconvolution_plot` | diversity / gridding / composition pies | ❌ | ❌ | `plots_full/` only — **in neither vignette** |

**Three plot directories, deliberately not merged:**

| Directory | Script | Scope | Count |
|---|---|---|---|
| `plots_tutorial/` | `plot_stlearn_tutorial.py` | **the Xenium vignette call-for-call — cite this one as "the authors' default workflow"** | 23 + 6 in `requested/` |
| `plots/` | `run_stlearn.py` | Xenium set + 4 generic-vignette calls; top-6 where the tutorial takes top-1/top-2 | 29 |
| `plots_full/` | `plot_stlearn_full.py` | adds `lr_plot` (generic) + `het_plot`/`grid_plot`/`deconvolution_plot` (neither) | 26 |

Requested LRIs are isolated in `requested/` in both `plots_tutorial/` and `plots/`.
`plot_stlearn_full.py`'s `deconvolution_plot` calls all failed (`KeyError: 'deconvolution'` — it
writes `uns`, the plot reads `obsm`) and its `feat_plot` covered only 4 of 9 cell types with a
per-type colour scale; both are superseded by `plots_tutorial/` and are recorded in
`DEVIATIONS.md` rather than repaired, since neither call is in the Xenium vignette.

`plot_stlearn_tutorial.py` rebuilds the cell-level `adata` (QC + `normalize_total`, ~30 s,
deterministic) because three vignette cells are grid-vs-single-cell side-by-side panels that
cannot be drawn from `grid.h5ad` alone. It asserts the rebuild against `run_manifest.json`
(100,197 cells, 5,096 genes, 12,562 spots) before plotting, so a silent preprocessing drift
aborts rather than producing a mismatched figure.

### Multi-sample / differential mode

**None.** No native multi-sample or between-condition test — no equivalent of CytoSignal's
`mergeCytoSignal` + `runNEBULA`. Any grade or AIS-vs-LUAD contrast would have to be hand-rolled,
which `SKILL.md:49` forbids. **Report this as a capability gap.**

### Gotchas

- **`run_cci` is broken on pandas 3 / numpy 2** — `obs[label].values.astype(str)` yields a
  `StringDtype`/object array numba 0.66 cannot type (`TypingError` at `het.py:227`). Needs
  object-dtype coercion or `pandas<3`.
- **Cell types are matched to deconvolution columns by SUBSTRING, first hit wins**
  (`get_data_for_counting`). Any label that is a substring of another (`mGAM` ⊂ `non-mGAM`) can
  silently bind to the wrong column and corrupt every interaction reported for it. **Verified
  correct here only because alphabetical order favours the shorter name.**
- `grid()` hard-requires `uns['spatial']`; a plain h5ad raises `KeyError: 'spatial'`.
- **`n_pairs` has a hard floor of 100** — below it `run()` exits with a message and stores nothing.
- **Cost scales with unique LR *genes*, not `n_pairs × n_LRs`** — backgrounds are cached per gene.
  100× the pairs cost only 7.6× the time, so there is no reason to run below the authors'
  recommended 10,000.
- **`n_cpus` does nothing measurable** — 3.5 min at 8 cores vs 3.4 min at 1.
- `lr_cci_map` raises `UnboundLocalError` on an empty `lrs` list.
- **`st.tl.cci.adj_pvals` is a no-op at the tutorial's own arguments.** The Xenium vignette calls
  it with `correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh'`, which `run()` has
  already applied internally — its own docstring says so. Verified empirically: all five `obsm`
  matrices and all six `lr_summary` columns come back identical. Its one real effect is
  re-sorting `lr_summary` with an **unstable** `argsort`, permuting **179 / 526** LRs *within
  ties on `n_spots_sig`*. Harmless but pointless; anyone diffing two runs on LR order should
  expect tie noise, not a result change.
- **`deconvolution_plot` reads `obsm['deconvolution']`, not `uns`.** Putting the frame in `uns`
  raises `KeyError` — and if the caller wraps plots in try/except, the failure is visible only
  in the log.
- **`grid()` stores per-spot proportions under `uns[<label>]`** while `obs[<label>]` is the
  dominant-type label. `feat_plot` needs the former copied into `obs` first.

### Deviations from the tutorial

Full table in `stlearn/DEVIATIONS.md`. Headlines: real `cell_type` instead of Leiden; grid
`321 × 146` = 51.3 µm square (preserving the tutorial's **spot area**, since its own spots are
60.2 × 43.8 µm — rectangular — on a 1.37:1 section, whereas our TMA is 2.20:1); `n_pairs`
1,000 → **10,000** and `n_perms` 100 → **1,000** (both tutorial-declared "example, recommend
higher"); 21 `Intergenic_Region_*` control probes dropped; CellChatDB v2 complexes dropped as
forced by the LR format.

**The ~51 µm grid is *not* aligned to ALARMIST's 50 µm patch** — it falls out of matching the
tutorial's spot area, and the proximity is a coincidence of the arithmetic.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `stlearn/GBM/cellchatdb2/` | ✅ 45.4 min | 100,197 cells → 12,562 spots (7.98 cells/spot); **526 LR pairs** tested; all 526 have ≥1 significant spot, **482 have ≥20**. Top: C3–C3AR1 (1,244 sig spots), CNTN2–CNTN2, GJA1–GJA1, APP–SORL1, C3–CR2, FGF1–FGFR2 |
| GBM | `default` | — | ❌ not run | connectomeDB2020_lit tier deferred |
| LUAD | both | — | ❌ not run | deferred |

Post-hoc passes, both replay-only (`grid.h5ad` opened read-only): the Xenium-vignette figure set
(`plot_stlearn_tutorial.py`, 0.8 min, 23 PNGs + 6 in `requested/`, all 29 succeeded) and the
stage-2 quantitative export (`export_stlearn_quant.py`, 0.04 min).

### Requested LRIs — and the first cross-method comparison

| LRI | stLearn rank | percentile | sig. spots / expressing | CytoSignal rank | percentile |
|---|---|---|---|---|---|
| GRN → SORT1 | **21 / 526** | top 4.0% | 267 / 8,763 | 66 / 895 | top 7.4% |
| ANXA1 → FPR1 | **99 / 526** | top 18.8% | 122 / 5,721 | 255 / 895 | top 28.5% |

Two methods with different units (spot vs cell), different nulls (expression-matched background
genes vs spatial permutation) and different neighbourhoods (250 µm grid radius vs 200 µm ε-ball)
**agree on the ordering and rough standing** of both LRIs. That concordance is meaningful
precisely because nothing about the two pipelines is shared. Their *top* rankings, by contrast,
are disjoint — partly by construction, since stLearn cannot see the complex-receptor WNT3
interactions that dominate CytoSignal's.

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
Contract: `spatialdm/NOTES.md`. Deviations: `spatialdm/DEVIATIONS.md`. Code: `run_spatialdm.py`,
`run_diff_spatialdm.py` (grade differential), `plot_spatialdm_full.py` (the corrected replot),
`compare_tiers.py`.

### Core algorithm

Tests **spatial co-expression** of a ligand and receptor with a **bivariate Moran's statistic**,
in two nested stages.

**Global (`spatialdm_global` → `sig_pairs`)** — per LR pair, a bivariate Moran's R asks whether
ligand expression at a location is spatially associated with receptor expression at its
*neighbours*, given a weight matrix `W`. The `z-score` method uses an **analytical null**: the
variance of R is derived in closed form from the structure of `W` itself, not from permutation.
BH-corrected, `fdr < 0.1` selects pairs.

**Local (`spatialdm_local` → `sig_spots`)** — for the globally selected pairs only, a local Moran
statistic per cell, giving where that interaction actually occurs. **The unit is the pair
globally and the cell locally.**

Two consequences. This is a **correlation/autocorrelation** framework, not a co-expression-product
permutation like CytoSignal or stLearn — it asks whether L and R are spatially *arranged*
together, not whether their product exceeds a background. And **because the null is analytical in
`W`, power depends directly on neighbourhood structure** — which is what makes the density
confound so consequential here.

### ⚠️ Overlaps with LIANA+ — these two are not independent votes

**LIANA+'s `morans` local metric is a reimplementation of this method**, and its source says so
(`liana/method/sp/_bivariate/_local_functions.py:458` cites the SpatialDM paper; the package even
carries a helper named `_spatialdm_weight_norm`). The formulae are identical, and for
`local_name='morans'` LIANA+ applies SpatialDM's exact preprocessing and offers the same analytic
z-score null at `n_perms=0`.

Everything around the statistic differs, and on our data that dominates:

| | SpatialDM `cellchatdb2` | LIANA+ `cellchatdb2_morans` |
|---|---|---|
| support | 13 cores, separate fits | whole slide, one run |
| radius | **134.6 µm** | **28.2 µm** |
| null | analytic z-score | permutation, `n_perms=100` |
| filter | `min_cell=3` | `nz_prop=0.02` |
| pairs tested | **1,133–1,661 per core** | **131** |

**A pair agreeing between the two is a sanity check, not corroboration.** SpatialDM's distinct
contributions are the analytic null (fast, demonstrably liberal) and `diff_utils`, a native
per-sample differential mode LIANA+ has no equivalent of. Note also that `morans` is *not*
LIANA+'s default local metric, and that Moran's R is not NMF-admissible (33.4% negative entries),
so LIANA+'s communication-program branch cannot run on it.

### Spatial model

| | |
|---|---|
| kernel | RBF, `w(d) = exp(-d² / 2l²)`, weights `< cutoff` zeroed |
| ours | `l = 75`, `cutoff = 0.2` → **effective radius 134.6 µm** |
| graph | `n_neighbors` per core (94–709), sized so `cutoff` truncates and the kNN cap never binds |
| adjacent graph | separate 6-NN graph → `obsp['nearest_neighbors']`, for contact/ECM pairs |
| `single_cell=True` | zeroes the diagonal: a cell cannot signal to itself |

`W` (RBF) is used for `Secreted Signaling`; `nearest_neighbors` for `ECM-Receptor` and
`Cell-Cell Contact` — the split is by database annotation, and `geneInter` **must be sorted by
annotation** because `st` and `n_short_lri` index blocks positionally.

### LR database

Default **CellChatDB v1** — 1,939 human interactions + 157 complexes. For `cellchatdb2` we inject
CellChatDB v2 directly, because `extract_lr` has no custom-DB hook. **Lossless for complexes** —
SpatialDM handles multi-subunit interactions natively via `Ligand0..N`/`Receptor0..N` columns, so
unlike stLearn (which must drop 57.5% of v2) nothing is lost for being a complex. Valid pairs per
core: **568–825 on v1**, **1,133–1,661 on v2**, versus stLearn's 526 whole-slide.

**`extract_lr`'s `datahost` argument is named backwards.** `datahost='package'` reads the CSVs
shipped inside the wheel; `datahost='builtin'` — the function's own **default**, and therefore
what the tutorial silently gets — **downloads from figshare**, once per call. We pass
`--datahost package`: identical CellChatDB v1, but offline, pinned, and not re-fetched 13 times.

**v2's `Non-protein Signaling` category is remapped to `Secreted Signaling`** — SpatialDM's
v1-era code does not enumerate it and crashes otherwise. Not cosmetic: core 13's entire top-6
(`SLC17A7_GLS2_GRIA2`, `SLC17A7_GLS_GRIN1_GRIN2A`, …) is glutamatergic and would have been lost
had the category been dropped instead.

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
| 5 | `spatialdm_local(n_perm, method='z-score')` | `uns['local_z']`, `uns['local_z_p']` |
| 6 | `sig_spots(method='z-score', fdr=False, threshold=0.1)` | `uns['selected_spots']`, `n_spots` |

### Data outputs

| File | Shape | Meaning |
|---|---|---|
| `global_res.csv` | pairs × 9 | Moran's R, p, **fdr**, `selected` — the ranking |
| `local_z.npz`, `local_z_p.npz` | selected-pairs × cells | **per-cell** local statistic and p-value |
| `selected_spots.csv.gz` | selected-pairs × cells | binary: where the interaction is called |
| `local_n_spots.csv` | pairs | number of selected cells per pair |
| `lr_{ligand,receptor}_subunits.csv` | pairs × subunits | resolved complex membership |
| `cell_meta.csv` | cells | coordinates + cell type |
| `per_split_summary.csv` | 13 rows | per-core diagnostics |

**`local_*` cover only globally-selected pairs**, by design.

### Image outputs

| Call | Shows | Ours |
|---|---|---|
| `pl.global_plot` | global selection volcano | ✅ per core |
| `pl.plot_pairs` | **per-LR spatial maps** (L, R, local significance) | ✅ top 3 + requested. The bundled `plot_pairs.pdf` is **truncated in cores 3/9/12** — `plot_selected_pair` indexes by *globally selected* pairs only, so it raises on the first non-selected requested LR and the `PdfPages` context closes early. The per-pair PNGs are unaffected |
| `pl.chord_celltype` | cell-type chord for a pair | ✅ top 2 per core (`obsm['celltypes']` built as one-hot — no deconvolution needed) |
| `pl.chord_LR` | self-self chord per cell type | ⚠️ cellchatdb2 **12/13**, default **7/13**. For `default` the log gives both causes: 4 cores raise `zero-size array to reduction` (no edge survives `min_quantile=0.5`, which the sparser v1 DB hits more often) and 2 fail reproducibly inside the headless-browser render. For `cellchatdb2` the single missing core (2) has **no log**, so its cause is unverified |
| `pl.dot_path` + `compute_pathway` | pathway enrichment | ⚠️ run per core but **degenerate**: `dic={core: all selected}` instead of the tutorial's `dic={Pattern_i: …}`, so there is one group and the background set is empty |
| `pl.ligand_ct` / `pl.receptor_ct` | per-cell-type contribution | ✅ as `celltype_contrib_*.{csv,png}` — **these two return DataFrames, they are not plotters** |
| SparseAEH `plot_clusters` | spatial clustering of local patterns | ❌ optional extra package, not installed; this is what would make `dot_path` meaningful |
| weight-range scatter | visualises the RBF kernel extent | ❌ diagnostic, not produced |

**Read `plots_full/`, not `plots/`.** `plots/` was written during the compute run using
`fn(); plt.gcf().savefig(...)`, the **wrong mechanism for three classes of SpatialDM function**,
and produced blank 7,544-byte PNGs: `plot_pairs` ends each iteration with `plt.show(); plt.close()`
so `gcf()` is a fresh empty figure; `ligand_ct`/`receptor_ct` are not plotters; `chord_*` render
through holoviews/bokeh and take their own `save=`. Only `global_plot` survives there.
`plots_full/` is the corrected pass (`plot_spatialdm_full.py`), replotted **from the persisted
`data/spatialdm.h5ad` with no recomputation**, with a saver that refuses to write a figure with
no drawn content. It holds **181 PNGs + 13 PDFs** (cellchatdb2) and **163 + 13** (default).

⚠️ **The `default` tier's `plots/` still contains 143 blank 7,544-byte PNGs** from the first pass
(130 in `plots/`, 13 in `plots/requested/`). They were cleaned from the `cellchatdb2` tree but not
here. The only real files in `default/*/plots/` are 13 `global_plot.png` and 13 `dot_path.png`.

A requested LR is plotted only where it is **globally selected** — SpatialDM computes no local
statistic otherwise, so there is nothing to draw. Per-core status is in `per_split_summary.csv`.

### Multi-sample / differential mode

**Yes — `spatialdm.diff_utils`**, demonstrated in `differential_test_intestine.ipynb`. Run via
`run_diff_spatialdm.py`: 7 high-grade vs 6 low-grade cores, `concat_obj` → `differential_test` →
`group_differential_pairs`, then the tutorial's clustermap / dendrogram / volcano / dot_path.
Nothing is recomputed — it consumes the persisted per-core objects. Two patches were required
(`concat_db` hard-codes CellChatDB v1; `dot_path`'s installed signature differs from the
tutorial's) — see `DEVIATIONS.md`.

### Grade differential — the method finds nothing, and says why

| | |
|---|---|
| union LR pairs across 13 cores | **1,662** (581 testable in all 13; 12.2% of the pair × core grid zero-filled) |
| pairs at differential FDR < 0.1 | **0** (minimum FDR **0.114**) |
| `high_specific` / `low_specific` | **0 / 0** |
| pairs selected in all 7 high cores and no low core | **0** |
| raw p < 0.05 | 162 of 1,662 (83 expected by chance) — signal exists, none survives BH at n=13 |

**The reason is the density confound, and three independent views agree.** (i) Per-core mean |z|
tracks neighbourhood size, **r = +0.752, p = 0.003** (core 13, `n_neighbors`=693 → mean |z| 5.67;
core 6, `n_neighbors`=94 → 0.60). The differential test regresses exactly that quantity.
(ii) `differential_dendrogram` clusters the 13 cores almost perfectly by grade **except cores 9
and 11** — the two where grade and density disagree (9 is low-grade but dense, 11 is high-grade
but sparse). (iii) Re-running the authors' own test with a **median-density split**, which differs
from the grade design at exactly those two cores, yields **5 condition-specific pairs** (all five
selected in all 7 dense cores and no sparse core: `CSF1_CSF1R` FDR 0.0075,
`BMP2_BMPR1B_ACVR2A`, `TNC_ITGAV_ITGB3`, `TNC_ITGAV_ITGB6`, `IL16_CD4`) against **0** for grade.

**Conclusion: on this TMA, SpatialDM's native differential test separates cellularity, not grade
— and the two are not separable with 13 cores.**

**The design is also largely paired, and the test does not use it.** `differential_test` fits
`y ~ 1 + conditions` by OLS with no patient term, so it discards the within-patient pairing *and*
treats patient 14007's four cores as four independent observations — pseudoreplication and lost
power at once. A paired/mixed model on the persisted `zscore_df` is the natural follow-up; it is
not SpatialDM's own workflow, so it was not run.

Requested LRIs in the grade contrast: **GRN→SORT1** diff +4.39, FDR 0.62 (mean z 6.07 high vs
1.68 low, selected in 5/7 high and 2/6 low); **ANXA1→FPR1** diff +0.38, FDR 0.97 (2.95 vs 2.56,
4/7 and 2/6). Neither is condition-specific by this method.

### Gotchas

- **`eff_dist` is a squared distance despite its name** — `l = sqrt(-eff_dist/(2·ln cutoff))`
  against a kernel `exp(-d²/2l²)`. `eff_dist=135, cutoff=0.2` gives `l=6.48` and a weight of
  4e-95 at 135 µm instead of 0.2. **Always set `l` directly.**
- **`n_neighbors` and `n_nearest_neighbors` are independent** — the `n_neighbor_layers*31`
  derivation fires only when `n_neighbors is None`; `n_nearest_neighbors` sizes a *separate* graph.
- **Once the kNN cap does not bind, `W` is numerically identical for any larger `n_neighbors`.**
  Only stored `nnz` differs, because `rbfweight` zeroes sub-cutoff entries without calling
  `eliminate_zeros()`. Compacting is safe and here saved **17.4M** entries.
- **CellChatDB v2's `Non-protein Signaling` crashes `spatialdm_global`** (shape mismatch in `st`).
- `sig_spots` **defaults to `fdr=True`** but the tutorial passes `False`.
- `local_z`/`local_z_p`/`local_perm_p` live at `uns` **top level**, not under `uns['local_stat']`.
- **`local_z_p` is NOT `norm.sf(local_z)`.** `spot_selection_matrix` ends with
  `np.where(pos.T == False, 1, local_z_p)` — the p-value is forced to 1 wherever *neither*
  standardised side is above its mean. On core 1 that masks 86.8% of the matrix. The identity
  holds only on the unmasked entries (verified, max abs diff 4.4e-08).
- **`local_z.npz` carries placeholder row labels** (`'0','1','2',…`): `uns['local_z']` is a bare
  ndarray, so the runner's `getattr(v,'index',arange)` fell through. Row order **is** identical to
  `local_z_p.npz`, which does carry the pair names. **Use `local_z_p`'s `pairs` array to label
  `local_z`.**
- **Core 13 has one cell that yields NaN local statistics** for 753 of its 863 selected pairs —
  column 8370, cell id 68901, an `AC-like` cell. Always the same one; drop it before any
  correlation over `local_z`.
- **`datahost='package'` needs `pkg_resources`, which setuptools ≥ 81 removed.** Pinned to
  `setuptools<81` (80.10.2) and frozen in `env.lock.yml`.
- **Chord plots need `geckodriver`/`firefox` on `PATH`, and they are not.** Both binaries *are*
  installed in `comp-spatialdm/bin/`, but this repo invokes the interpreter by absolute path
  (conda activate is broken here), which leaves that `bin/` off `PATH`. Both runners now prepend
  `dirname(sys.executable)` to `PATH`.
- **`chord_LR(senders=…, receivers=…)` indexes `adata.obs[sender]`** — it wants one obs *column*
  per cell type, not category names, and the two lists are **zipped, not crossed** (so it draws
  the self-self diagonal). It also still hardcodes `title='Undifferentiated_Colonocytes'` from the
  authors' tutorial.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `spatialdm/GBM/cellchatdb2/<core>/` | ✅ 9.3 min, 13 cores | 100,197 cells; **1,133–1,661 valid pairs** per core; **54–863 significant** (FDR<0.1); kNN cap binds **0.000% in every core**; stored nnz 47.3M → **29.9M** after `eliminate_zeros` |
| GBM | `cellchatdb2` | `.../differential_grade/` | ✅ 0.4 min, peak RSS 3.6 GB | **0 pairs at FDR < 0.1** |
| GBM | `default` | `spatialdm/GBM/default/<core>/` | ✅ 12.9 min, peak RSS 6.0 GB | same cells and same `W`; **568–825 valid pairs** per core; **33–386 significant**; `--datahost package` |
| GBM | `default` | `.../default/differential_grade/` | ✅ 0.4 min | **0 pairs at FDR < 0.1** — reproduces the v2 tier exactly |
| — | comparison | `.../GBM/tier_comparison/` | ✅ | `compare_tiers.py`: v1 ↔ v2 concordance |
| LUAD | both | — | ❌ not run | deferred |

### `default` vs `cellchatdb2` — the LR database is not driving anything

Both tiers use the **same cells, same `W`, same `n_neighbors`** (verified identical per core), so
they differ only in DB membership. Pairs are matched across tiers by resolved subunit sets, not by
name (v1 names complexes `TGFB1_TGFBR1_TGFBR2`; our v2 names are `ligand_receptor`).

| | v1 (`default`) | v2 (`cellchatdb2`) |
|---|---|---|
| DB size | 1,939 interactions | 3,233 rows / 3,218 unique pairs |
| pooled valid pairs over 13 cores | 9,132 rows → **8,528 distinct** | 18,971 rows → **16,057 distinct** |
| pooled significant | **2,025** | **5,201** |
| of v2's pairs, share also in v1 | — | **52.6%** (median per core) |

1. **For any pair present in both DBs, SpatialDM computes the identical statistic — pooled
   Pearson r of the Moran z across all shared pairs = 1.000.** The database cannot change what
   the method says about a pair; only *which* pairs are tested, and through that the BH denominator.
2. **Agreement on shared pairs is near-total: median Jaccard of the selected sets = 0.954**
   (0.78–1.00). Disagreements sit on the FDR boundary, moved by v2 testing twice as many hypotheses.
3. **The two requested LRIs behave identically.** GRN→SORT1 testable in 13/13 and selected in
   exactly cores {1, 2, 5, 6, 8, 10, 13} under **both** DBs; ANXA1→FPR1 testable in the same 11
   cores and selected in exactly {1, 3, 8, 9, 12, 13} under both. Not a v2 artefact.

**So the `Non-protein Signaling` remap inflates the counts but does not change the conclusions.**
It supplies **48.6% of v2's distinct significant pairs** (2,065/4,248). Removing it entirely —
which is what the v1 tier effectively does — leaves the requested-LRI verdict, the grade
differential and the density confound unchanged. **The density confound is not a v2 artefact
either; it is stronger on v1** (corr(`n_neighbors`, fraction significant) r = +0.631, p = 0.021 on
v1 vs +0.442, p = 0.131 on v2), and the grade differential reproduces to the digit, returning
**the identical five density-split pairs**.

### Requested LRIs

Identical under both LR databases — same testable cores, same selected cores
(`tier_comparison/requested_lr_by_tier.csv`):

| LRI | testable | selected (FDR<0.1) | which cores | high-grade | low-grade |
|---|---|---|---|---|---|
| GRN → SORT1 | **13/13 cores** | 7 | 1, 2, 5, 6, 8, 10, 13 | 5/7 | 2/6 |
| ANXA1 → FPR1 | 11/13 cores | 6 | 1, 3, 8, 9, 12, 13 | 4/7 | 2/6 |

Both lean high-grade, but with 13 cores from 7 patients this is **not** a significant association
(Fisher p = 0.29 and 0.59) and must not be read as one. GRN→SORT1 being testable in **every** core
is itself notable.

### Methods paragraph

> For SpatialDM (v0.3.1), we followed the authors' workflow: we computed the neighbour-weight
> graph, identified globally significant ligand–receptor pairs and then selected locally
> significant spots, using `weight_matrix`, `extract_lr`, `spatialdm_global`, `sig_pairs`,
> `spatialdm_local` and `sig_spots`. The radial-basis kernel was set to `l = 75` with
> `cutoff = 0.2`, reproducing the effective 135 µm signalling radius of the authors' intestine
> tutorial — the one example in which they state their reasoning for this parameter — and
> `single_cell=True` was used because Xenium resolves individual cells. The size of the neighbour
> graph was set per tissue core so that the distance cutoff, rather than the k-nearest-neighbour
> ceiling, determined the neighbourhood, matching the regime of both author tutorials. Each of the
> 13 tissue microarray cores was analysed independently, following the authors' multi-sample
> design, because the global bivariate Moran statistic is returned once per ligand–receptor pair
> per object and the differential test operates on separately fitted samples. The analysis was run
> twice, once with SpatialDM's own CellChatDB v1 resource and once with CellChatDB v2, whose
> multi-subunit complexes SpatialDM represents natively; for ligand–receptor pairs present in both
> resources the two runs gave identical Moran statistics (Pearson r = 1.000) and near-identical
> significance calls (median Jaccard 0.954). To compare high- and low-grade cores we used
> SpatialDM's native differential mode (`spatialdm.diff_utils.concat_obj`, `differential_test`,
> `group_differential_pairs`), which fits a likelihood-ratio test of each pair's per-core Moran
> z-score on the grade label; no ligand–receptor pair reached a differential false-discovery rate
> below 0.1.

---

## COMMOT — Python, v0.0.3, env `comp-commot`

Tutorials: `/Users/jiayifan/tansey_lab/COMMOT/docs/notebooks/{Basic_usage,visium-mouse_brain}.ipynb`.
Contract: `commot/NOTES.md`. Deviations: `commot/DEVIATIONS.md`. Code: `run_commot.py`,
`plot_commot_vf.py` (vector-field replot), `compare_motif1_commot.py` (vs ALARMIST),
`run_commot_impact.py`, `run_commot_deg.py` (**never executed** — see *Downstream impact*),
`SETUP_tradeSeq.md`.

### Core algorithm

Casts communication as **collective optimal transport**. Ligand "mass" at senders is transported
to receptor "mass" at receivers so as to minimise a spatial cost, subject to `dis_thr` forbidding
any coupling beyond that distance. The defining property — the authors' main argument over
pairwise scoring — is that transport is **competitive**: multiple ligands and receptors compete
for the same finite mass, so a receptor already saturated by a nearby strong sender cannot also
absorb signal from a distant one.

Output is a **cell × cell transport plan per LR pair**, summarised into per-cell amounts sent and
received. Crucially, **there is no significance test at the LR-pair level** — COMMOT returns
magnitudes, not p-values. Permutation p-values exist only at the **cell-type-pair** level, per
pathway, via `cluster_communication`. Any ranking of LR pairs here is by *total received signal*,
a magnitude, and is **not comparable like-for-like** to CytoSignal's significant-cell counts,
stLearn's significant-spot counts, or SpatialDM's FDR.

### Spatial model

A single hard Euclidean cutoff, `dis_thr`, in the units of `obsm['spatial']`. No kernel, no decay
— coupling is either permitted or forbidden. **Ours is 365 µm**, derived by measuring the tutorial
dataset rather than copying its number (see Gotchas). Applied identically to every pair, which is
why only diffusible signalling types belong in the database subset.

### LR database

Built-in options are **CellChatDB (v1)** and **CellPhoneDB v4.0**, via
`ct.pp.ligand_receptor_database`. `df_ligrec` is a plain 3-column frame — ligand, receptor,
pathway — with heteromeric subunits joined by `_`, which is **exactly the encoding our CellChatDB
v2 CSV already uses**, so the v2 handover is direct and lossless. Complexes are supported natively
(`heteromeric=True`, `heteromeric_rule='min'`).

We use **Secreted Signaling + Non-protein Signaling** — 2,259 pairs, 859 on the Xenium panel,
**671 surviving `min_cell=100` evaluated once over all 100,197 cells** and reused for every core
(`--filter-scope global`). The tutorial restricts to Secreted only; we keep that principle — a
single 365 µm transport radius is a diffusion model — but add v2's Non-protein category, which did
not exist in the v1 the tutorial used and is equally diffusible. **That inclusion is vindicated by
the results: Glutamate is a top-3 pathway in cores 9, 13 and 14**, matching SpatialDM's
independent finding that glutamatergic pairs dominate core 13.

### Input

`adata.raw = adata`; `sc.pp.normalize_total`; `sc.pp.log1p` — applied to **raw counts**. COMMOT
requires "non-negative values that reasonably reflect the abundancy of signaling molecules".

### Workflow

| # | Call | Produces |
|---|---|---|
| 1 | `ct.pp.ligand_receptor_database(...)` *(or our v2 frame)* | `df_ligrec` |
| 2 | `ct.pp.filter_lr_database(df, adata, filter_criteria, min_cell)` | expression-filtered pairs |
| 3 | `ct.tl.spatial_communication(adata, database_name, df_ligrec, dis_thr, heteromeric=True, pathway_sum=True)` | transport plans + per-cell sums |
| 4 | `ct.tl.cluster_communication(..., pathway_name=<top-5 by signal>, clustering='cell_type', n_permutations=100, random_seed=0)` | cell-type × cell-type matrix **+ permutation p-values** |
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
| `cluster_comm_<key>.csv`, `cluster_pval_<key>.csv` | types × types | cell-type communication + p-values. `<key>` = each top-5 pathway **and** each requested LR pair → 7 pairs of files per core |
| `lr_pairs_used.csv` | pairs × 3 | the post-filter database actually used (identical in all 13 cores) |
| `lr_pairs_global.csv` | 671 × 3 | *(run root)* the one globally-filtered pair set |
| **`<core>/adata_commot.h5ad`** | — | **the full AnnData including `obsp`** — every transport plan, the sum frames, the vector fields, the 7 `cluster_communication` results and `cluster_pos`. `.raw` dropped (reconstructible). 0.01–0.60 GB/core, **2.13 GB** total |

**Why the h5ad matters:** COMMOT has no way to recover a transport plan except by re-solving the
OT. Persisting `obsp` is what makes `communication_direction`, `communication_impact`,
`deg_detection`, `spatial_autocorrelation` and the `group_*` family reachable without a fresh
~5-hour run. Verified by reloading `2/adata_commot.h5ad` and running `communication_direction`
on it with no OT re-run.

**Trap:** `pathway_sum=True` writes pathway columns *into the same* `sum_receiver` frame as the
per-pair columns. **Ranking pairs without separating them inflates the denominator and shifts
every rank.** `lr_total_received.csv` and `pathway_total_received.csv` are now correct and disjoint.

### Image outputs

**41 figures per core × 13 cores = 533** (507 PNG + 26 PDF).

| Plot | What it shows | File | n |
|---|---|---|---|
| `sender_receiver_map` *(ours)* | per-cell sent/received, top **LR pairs** by total received | `signal_<lig>-<rec>.png` | 3/core |
| `sender_receiver_map` *(ours)* | same, **pathway** aggregates — labelled so they cannot be read as pairs | `pathway_<pw>.png` | 5/core |
| `sender_receiver_map` *(ours)* | the standing requested LRIs | `requested_rank<N>_<LR>.png` | 2/core |
| `communication_direction` → **`ct.pl.plot_cell_communication`** | signalling vector field, **tutorial arguments** (`plot_method='grid'`, `normalize_v=True`, `normalize_v_quantile=0.995`, `grid_density=0.4`, `ndsize=8`). Two backgrounds × sender/receiver × (5 pathways + 2 requested) | `native_vf_{sender,receiver}_*_{sig,ct}.png` | **28/core** |
| **`ct.pl.plot_cluster_communication_network`** | cell-type network + colour legend | `native_cluster_network*.pdf` | 2/core |
| `ours_dotplot` *(ours — substitute)* | cell-type × cell-type dots, size = −log10(p), colour = strength | `ours_dotplot_top.png` | 1/core |
| `ct.pl.plot_cluster_communication_dotplot` | — | ❌ **unrunnable, version gap** (see Gotchas) | 0 |
| `ct.pl.plot_communication_impact` | downstream-impact heatmap | ⚠️ partial — see *Downstream impact* | — |
| `ct.pl.plot_communication_dependent_genes` | signalling-DE gene heatmap | ❌ **blocked** — its input needs R/tradeSeq | 0 |

**`scale` is the one vector-field value that cannot be copied verbatim** — it reaches
`quiver(scale_units='x')`, so arrow length is `|v|/scale` in *data* units, and the tutorial's
`0.00003` is tuned to ~9,000-unit Visium pixel coordinates against our ~1,800–2,400 µm cores.
Rescaled per core to hold arrow length at the tutorial's fraction of the field (1.1e-4 – 1.5e-4).
Same class of unit trap as `dis_thr`.

### Multi-sample / differential mode

**None native.** COMMOT has no cross-sample differential test; `communication_deg_detection`
relates signalling to gene expression *within* a sample. A grade contrast would have to be
hand-rolled from per-core outputs.

### Gotchas

- **Unusable on numpy ≥ 2** — `np.Inf` as a module-level default in `_usot.py`. Requires `numpy<2`.
  (Fixed upstream, but not in the 0.0.3 release we run.)
- **A dense N×N distance matrix is materialised** → the full slide would need **80.3 GB**.
  Per-core is mandatory, not optional.
- **The tutorial's prose contradicts its own units** — `dis_thr=500` is described as "500 µm" but
  Visium coordinates are full-resolution pixels at 0.73 µm/unit, so the real constraint is
  **365 µm**. Copying 500 onto micron data gives a 37% wider neighbourhood.
- **`min_cell_pct=0.05` does not transfer to single cells** — 5% of multi-cell Visium spots vs 5%
  of individual cells retains 0.9–1.8% of pairs instead of the tutorial's 20.9%.
- **Runtime is not predictable from size** — the OT solver runs to convergence and iteration count
  varies per core. At a uniform 671 pairs, core 10 (17,435 cells) took **91.1 min** while core 1
  (26,456 cells) took **47.2 min**. **Do not project runtime from cell count.**
- **`ct.pl.plot_cluster_communication_dotplot` cannot run on a current matplotlib/seaborn.**
  commot 0.0.3 was written against matplotlib <3.9 and seaborn <0.13; installed are **3.10.9** and
  **0.13.2**. Two breakages: `_plotting.py:788` reads `legend.legendHandles` (removed in mpl 3.9),
  and past that seaborn 0.13 hands back `Line2D` rather than `PathCollection`. Recorded, not
  patched. `ours_dotplot_top.png` substitutes.
- **`plot_cluster_communication_network` needs pygraphviz**, not a COMMOT dependency. Installed
  1.14 into `comp-commot`; without it the call raises `ImportError` and the plot is **silently lost
  to a try/except**.
- **The local clone `/Users/jiayifan/tansey_lab/COMMOT/` is NOT the installed package** — it
  carries two post-0.0.3 commits in `_cot.py`. Verify signatures against the installed package.
- **The OT is normalised per run.** `_cot.py:269-271` divides both marginals by
  `max_amount = max(S.sum(), D.sum())` and `:335` multiplies the plan back, so units are restored —
  but `eps_p` and `rho` are fixed constants acting on the *normalised* masses, so the shape of the
  solution depends on which cells are in the run. **Ranks are comparable across cores, magnitudes
  are not.**
- No LR-pair-level significance; magnitudes only. Cell-type-level p-values exist for whatever
  `cluster_communication` was pointed at — and pointing it at the *alphabetically* first pathways
  instead of the strongest ones is an easy and silent mistake.

### Runs on our data

| Dataset | Tier | Path | Status | Key numbers |
|---|---|---|---|---|
| GBM | `cellchatdb2` | `commot/GBM/cellchatdb2/<core>/` | ✅ **296.1 min, 13/13 cores**, zero failures | 2,259 input pairs → **671 used in every core**; 771 obsp keys/core (671 pairs + 99 pathways + total); obsp up to **63.8M** nonzeros; peak RSS **11.5 GB**; **2.4 GB**, **533 figures**, 262 CSVs, 13 `adata_commot.h5ad` |
| GBM | `default` | — | ❌ not run | bundled CellChatDB v1 tier deferred |
| LUAD | both | — | ❌ not run | deferred |
| GBM | impact | `commot/GBM/impact/` | ⚠️ **PARTIAL — do not cite** | see *Downstream impact* |

**All 13 cores are analysable.** Under a per-core expression filter, cores 2 and 6 retained zero
pairs and were reported as a method limitation on sparse tissue. Both are **low-grade**, so
dropping them left 7 high + 4 low against the TMA's true 7 + 6 — a grade-biased analysis set.
Evaluating the filter once on all 100,197 cells restores full coverage **and removes the
pair-count confound**: every core now uses exactly 671 pairs, so ranks are cross-core comparable.
Magnitudes still are not (per-run OT normalisation). **Compare ranks, not magnitudes.**

### Requested LRIs — the strongest single signal of the benchmark

Every core ranks the same 671 pairs, so these ranks share one denominator.

| LRI | tested | median rank | range | |
|---|---|---|---|---|
| **GRN → SORT1** | **13/13 cores** | **2** (of 671) | **1–12** | #1 in 5 cores, top-3 in 9, **top-5 in 12** |
| ANXA1 → FPR1 | **13/13 cores** | 21 | 4–126 | never #1; top-5 in 2 |

COMMOT puts GRN→SORT1 in the top 5 of **12 of 13 cores** out of 671 candidates. The one exception
is core 2 (rank 12), the smallest at 819 cells. That is far higher than CytoSignal (66/895) or
stLearn (21/526) place it, but **the statistic differs**: COMMOT ranks by transported signal
magnitude, which favours abundantly expressed ligands, whereas the others rank by significance
against a null. **The agreement worth noting is directional** — all methods place GRN→SORT1 well
above ANXA1→FPR1 — and COMMOT's margin should not be read as stronger evidence.

Rank-1 LR pair per core: **FGF1–FGFR2** (7 cores), **GRN–SORT1** (5), WNT5A–FZD3 (1). Most frequent
members of a core's top 5: **GRN–SORT1 (12/13)**, FGF1–FGFR2 (7), PDGFA–PDGFRA (6), ANGPTL2–TLR4
(5), PDGFC–PDGFRA (5), CSF1–CSF1R (4). Top-5 pathways by frequency: **FGF 13/13, PDGF 13/13**,
Glutamate 8, GAS 6, BMP 6, COMPLEMENT 5, ANGPTL 3, GRN 3, ncWNT 3, GALECTIN 2, IGF 2, PROS 1.

`cluster_communication` was additionally run with `lr_pair=` for both requested LRIs, so each core
has cell-type × cell-type permutation p-values for them — **the only significance COMMOT can give
a named pair**.

### What COMMOT found in this tumour — plain reading

**1. Growth-factor axes dominate, unambiguously.** FGF and PDGF are top-5 pathways in 13 of 13
cores — concretely `FGF1–FGFR2` (rank 1 in 7 cores), `PDGFA–PDGFRA` and `PDGFC–PDGFRA`. The least
surprising and most reassuring result in the benchmark: PDGFRA is a canonical glioma driver, and a
method with no knowledge of this tumour recovered both axes from 671 candidates. **This is the
sanity check that the run is not mis-specified.**

**2. Glutamatergic signalling is the most interesting hit.** Glutamate is top-5 in 8 of 13 cores.
It comes from v2's `Non-protein Signaling` class, which we added on top of the tutorial's
Secreted-only restriction — so the deviation paid off; it **independently corroborates SpatialDM**,
which ranked glutamatergic pairs as core 13's entire top-6 under a completely different statistic;
and it points at neuron–glioma synaptic signalling.

**3. Myeloid cells are the communication hub; tumour cells are comparatively quiet** — mGAM tops
**both** directions (received 0.804, sent 0.777, mean within-core percentile) while every tumour
state sits near the median (0.455–0.483). ⚠️ **Discount this.** COMMOT's magnitudes grow with the
number of cells inside the 365 µm radius, so "myeloid is the hub" is partly a density artifact and
COMMOT cannot separate the two.

**4. The requested pair is near the top everywhere, but carries no grade signal.** Neither pair's
rank differs by grade (`GRN→SORT1` rank 1 in high vs 2 in low, p = 0.413), whereas ALARMIST's
motif-1 loading does (p = 0.022). **COMMOT sees the interaction but not its clinical association.**

### The vector-field figures — what they are, and whether they mean anything

`communication_direction` collapses the transport plan to one 2-D arrow per cell — the k=5-weighted
average *direction* in which that cell's signal flows — and `plot_cell_communication` interpolates
onto a grid. **One arrow = "in this patch of tissue, this pathway's signal is on average flowing
that way."**

⚠️ **The most common misreading: the arrow is a SPATIAL direction, not a cell-type direction.** It
says "towards the upper right of the image", **not** "from mGAM to MES-like". For sender→receiver
relationships between cell types you need `cluster_communication`'s matrix and its p-values.

**On this dataset they are presentational rather than analytical.** (1) **No test** — neither
direction nor length has a null. (2) **They largely trace tissue geometry** — transported mass
flows from dense to sparse regions. (3) **Edge artifacts** — the longest arrows sit at the tissue
perimeter, where one side of the interpolation kernel is empty. (4) **Wrong setting** — the
tutorial's mouse-brain section has strong anatomical gradients; a 1–2 mm TMA punch of relatively
homogeneous tumour has no comparable macro-gradient.

Still useful as a fast visual check of whether a pathway's signal is diffuse or concentrated, and —
for the `_ct` variant, whose background is cell type — of what the high-signal regions are made of.
**If only one background is kept, keep `_ct`.** **Do not** cite them as evidence that communication
is directional.

### Does COMMOT independently corroborate ALARMIST motif 1? — mostly **no**

`compare_motif1_commot.py` → `commot/GBM/vs_alarmist/`. Joins COMMOT's per-cell sent/received
amounts to `results/GBM/single_cell/cell_loadings.npy` (100,197 × 20, motif 1 = column 1). All 13
cores, all cells; join key is `obs_names` and every COMMOT cell id resolves. **Alignment guard:**
the script refuses to run unless motif 1 peaks on mGAM — it does (5.15e-3), and is the only column
checked that does, corroborating both the row order and the column index. *(Motif 1's margin over
Lymphoid, 5.12e-3, is razor-thin — motif 1 is high in mGAM, Lymphoid and Vascular alike, so "the
mGAM motif" is not a clean per-cell-type label.)*

**The decisive control kills the naive reading.** COMMOT's grand totals over all 671 pairs
correlate with motif-1 loading **better than any motif-1-specific quantity does**:

| COMMOT quantity | median Spearman vs motif 1 (13 cores) | best-matching motif |
|---|---|---|
| **`s-total-total`** (all signalling sent) | **0.167** | motif 1 |
| **`r-total-total`** (all signalling received) | **0.156** | motif 1 |
| `composite-mGAM-side` | 0.119 | motif 1 |
| `s-ANXA1-FPR1` | 0.089 | motif 9 |
| `r-ANXA1-FPR1` | 0.084 | motif 1 |
| `s-GRN-SORT1` | 0.079 | motif 1 |
| `r-GRN-SORT1` | 0.024 | motif 6 |
| `r-`/`s-FGF1-FGFR2` *(unrelated-pair control)* | −0.014 / 0.004 | motif 6 |

So motif-1 loading tracks **"how much signalling this cell does at all"** better than it tracks
COMMOT's estimate of the two LRIs that *define* motif 1. The FGF1–FGFR2 control does its job, but
the total-total control shows the shared variance is a **generic communicativeness / local-density
axis**, not pair-specific agreement. All effect sizes are small (|ρ| ≤ 0.17).

**Directional test, before and after subtracting the matching `total-total` profile:**

| Prediction | raw rank | baseline-corrected rank | residual |
|---|---|---|---|
| `s-GRN-SORT1` → mGAM | 1/9 ✅ | **7/9 ❌** | −0.121 |
| `r-GRN-SORT1` → MES-like | 4/9 ❌ | 2/9 (top: Glial-Neuronal) | +0.012 |
| `s-ANXA1-FPR1` → MES-like | 3/9 ❌ | **1/9 ✅** | +0.075 |
| `r-ANXA1-FPR1` → mGAM | 1/9 ✅ | **8/9 ❌** | −0.187 |

**Correction reverses the raw reading**: the mGAM sides, which looked confirmed, fall below what
mGAM's general communicativeness predicts; the MES-like sides, which looked refuted, rise. Only one
prediction lands first — MES-like as the ANXA1 sender — and even that is suggestive rather than
decisive, since MES-like sits mildly above baseline for *all four* quantities.

**Bottom line — the two methods agree about *which LRIs matter*, not about *which cells are doing
it*.** COMMOT puts GRN→SORT1 in the top 5 of 12/13 cores out of 671, which is strong independent
support for the pair's importance. It provides **no independent per-cell support for motif 1's
cell-type attribution**. Report the concordance at the LRI level; report cell-type attribution and
clinical association as ALARMIST's claims.

### Downstream impact — half of COMMOT's chain is unreachable on this machine

The tutorial's chain is `communication_deg_detection` → `communication_deg_clustering` →
`plot_communication_dependent_genes` → `communication_impact` → `plot_communication_impact`.

**Step 1 needs R.** `communication_deg_detection` imports `rpy2` + `anndata2ri` **inside its own
body** and drives `tradeSeq::fitGAM` + `associationTest` and `clusterExperiment`. It is the
**only** function in this otherwise pure-Python package that touches R — `pip show commot` declares
no R dependency and every other call works without it. On this machine `rpy2`/`anndata2ri` are
absent and **none of the four R installs has tradeSeq**. Consequently `communication_deg_clustering`
and `plot_communication_dependent_genes` — both pure Python — are starved of input.
**`SETUP_tradeSeq.md`** is the plan for enabling it in an isolated `comp-commot-r` env (nothing in
it has been executed); **`run_commot_deg.py`** implements the authors' full chain against it — also
never executed, written from the installed source, syntax-checked only.

**Steps 4–5 are runnable**, because `communication_impact`'s docstring sanctions a substitute
("A list of genes … for example, the highly variable genes"). `run_commot_impact.py` runs it with
every tutorial argument intact and two `ds_genes` variants — `hvg` (method-faithful) and `alarmist`
(top genes from `results/GBM/impact/motif_1_celltype_mGAM_de_results.csv`, so both methods are
asked about the *same* genes).

> ⚠️ **The substitute run was deliberately stopped.** 9 of 13 cores of the `GRN × alarmist` variant
> completed (~50 min) before the decision to install the R dependency and run the authors' actual
> chain instead. The partial output is kept, clearly marked, at `commot/GBM/impact/`
> (`PARTIAL_RUN_DO_NOT_USE_AS_FINAL.md`); it has no `run_manifest.json`/`impact_summary.csv`, which
> is the marker of an incomplete run. The `hvg` control and the proposed motif-14 control never ran.
> **Nothing here is a finished result.**
>
> One uncontrolled observation from the completed cores: on core 6 the receiver-side scores sat far
> above the 0.5 null (`r-GRN-SORT1` median 0.888, 28/30 genes > 0.8) while the sender side sat on it
> (median 0.519). Direction is what one would expect, but **the confound is not excluded**: motif-1's
> mGAM genes are myeloid genes and GRN→SORT1 receipt concentrates in myeloid-dense regions, so any
> mGAM-associated gene set might score the same. The control that would settle it — `ds_genes` from
> motif 14's mGAM impact table, same method and cell type but only 6/30 gene overlap — was designed
> and not run. **Do not cite the agreement until it exists.**

Three implementation facts worth keeping:

- **It requires `adata.raw`.** `run_commot.py` drops `.raw` before writing `adata_commot.h5ad`, so
  `run_commot_impact.py` rebuilds it from the source h5ad's `layers['counts']` — verified genuine
  integer counts (max 128). Still no OT re-run.
- **`pathway_name='GRN'` is a degenerate, convenient case.** The GRN pathway contains exactly
  **one** LR pair, so `df_impact` is 4 × (N+1) and the analysis is literally "what does GRN→SORT1
  signal explain". The figure needs `cluster_knn` below the row count: the default 5 raises
  `IndexError`; the tutorial hits the same wall and passes 2; GRN's 2 rows need **1**. The runner
  falls back 5 → 2 → 1.
- **The score is a percentile against background genes, and its null is 0.5 — not 0.**
  `_similarity.py:85` fits `n_repeat` random forests on `[signal feature | ~500 background genes]`,
  takes the signal's rank among the importances, and returns `mean((n_bg - rank)/n_bg)`. So 0.97
  means "more important than ~97% of background genes" and a useless feature lands near **0.5**.
  This is why the tutorial's own table sits entirely in 0.78–0.999 and looks uniformly high.
  **Always read these against 0.5.**

Like ALARMIST's GLM, `communication_impact` excludes the LR genes themselves
(`exclude_lr_genes_list`) — the same anti-circularity guard. Cost: **3.9 min per core at 3,092
cells / 30 genes**; `RandomForestRegressor` already runs at `n_jobs=-1`, so running variants in
parallel buys nothing.

### What COMMOT's output objects are, beside ALARMIST's

All shapes **measured on disk** (COMMOT: core 13, 9,126 cells, one such object per core ×13;
ALARMIST: `results/GBM/`).

**The one structural fact that explains everything else: ALARMIST has a latent axis; COMMOT does
not.** Every ALARMIST matrix carries **motifs** (K = 20) on one side — a learned axis that does not
exist in the data. Every COMMOT object is indexed by **observed LR pairs** (671) and by **cells**.
There is no matrix anywhere in COMMOT whose axis corresponds to ALARMIST's motif axis, and no way
to construct one.

A second difference: **ALARMIST's "LRI" axis already contains the cell-type direction.** Its 25,271
columns are `celltype1|celltype2|ligand|receptor|signaling_type`. COMMOT's 671 LRIs are *only*
ligand×receptor; cell type enters exclusively at `cluster_communication`. That is why ALARMIST's
factor matrix can express "mGAM→MES-like via GRN→SORT1" as a single weighted feature and COMMOT
cannot.

| ALARMIST | COMMOT | Comment |
|---|---|---|
| `cell_loadings` **100,197 × 20** (cell × motif) | `sum-sender`/`sum-receiver` **9,126 × 771** per core (cell × LRI) | closest pair. But: latent motif vs observed LRI; one number vs a sent/received split; one object over all cells vs 13 objects with **no shared cell axis** |
| `lri_factors` **20 × 25,271** (motif × LRI) | — **nothing** | COMMOT never groups or weights LRIs. This is the missing piece that makes motifs impossible |
| impact **gene × motif × cell type** (180 tables) | — **nothing produced** | `communication_impact`/`deg_detection` exist but are per-pathway within one sample |
| cell-type heatmap per motif (derived from V) | `cluster_communication` **9 × 9** + p-values, per key | ALARMIST's is per *motif* and has no p-value; COMMOT's is per *LR pair or pathway* and does. Complementary |
| — **nothing** | `obsp` transport plans **671 × 9,126 × 9,126** | ALARMIST produces no cell-to-cell edges at all. COMMOT's unique contribution |
| patch loadings **13,113 × 20** | — **nothing** | no spatial-unit abstraction in COMMOT |

**In one line:** ALARMIST factorises into *programs* and is fine-grained along the program axis;
COMMOT resolves individual *cell-to-cell edges* and is fine-grained along the cell axis. **Neither
can be converted into the other.**

### Methods paragraph

> For COMMOT (v0.0.3), we performed spatial communication inference with the
> `spatial_communication` function, following the authors' tutorial. Counts were normalised and
> log-transformed, and ligand–receptor pairs were taken from CellChatDB v2 restricted to the
> diffusible signalling classes, since COMMOT applies a single distance constraint to all pairs;
> heteromeric complexes were handled natively with `heteromeric=True` and the minimum rule. Pairs
> were filtered with `filter_lr_database` requiring at least 100 cells expressing each side, an
> absolute criterion chosen because the tutorial's 5%-of-spots threshold does not transfer to
> single-cell resolution; the filter was evaluated once across the whole specimen and the resulting
> 671 pairs applied to every core, so that pair sets and hence rankings are comparable between
> cores. The spatial distance constraint was set to 365 µm, reproducing the physical constraint of
> the authors' own example after converting their Visium pixel coordinates to microns.
> Cell-type-level communication and its permutation p-values were obtained with
> `cluster_communication` (`n_permutations=100`, fixed random seed) for the five pathways carrying
> the most received signal in each core and, separately, for each ligand–receptor pair of interest.
> Signalling directions were interpolated with `communication_direction` and visualised with
> `plot_cell_communication`. Each of the 13 tissue microarray cores was analysed separately, as
> COMMOT materialises a dense pairwise distance matrix that is intractable at whole-slide scale.

---

## LIANA+ — Python, v1.8.1, env `comp-liana`

Tutorials: `/Users/jiayifan/tansey_lab/liana-py/docs/notebooks/{bivariate,inflow_score,inflow_mofaflex,misty,sma}.ipynb`.
Contract: `liana/NOTES.md`. Deviations: `liana/DEVIATIONS.md` (1,143 lines — the fullest deviation
record of any method here, and the authority for everything marked "see DEVIATIONS.md" below).

**LIANA+ is not one method but a decision tree.** Eight of its branches were run. This section is
organised as: the shared model → a branch map → one subsection per branch → the tier comparison →
the relation to ALARMIST → the two open contract deviations.

### Core algorithm

LIANA+'s spatial mode computes a **spatially-weighted bivariate similarity** between a ligand and
a receptor — by default a weighted cosine, one of six metrics (`li.mt.bivariate.show_functions()`).
One call returns two levels simultaneously:

- **Local** — a score *per cell per LR pair* (`lrdata.X`), with permutation p-values from shuffling
  cell labels (`layers['pvals']`) and categorical labels (`high-high`, `low-low`, …) in
  `layers['cats']`.
- **Global** — per-pair summaries in `lrdata.var`: `mean`, `std`, and **bivariate Moran's R**
  (Lee's statistic) with `morans_pvals`.

**LIANA+ is the only method in this benchmark that natively returns both a per-cell score and a
per-pair global statistic with p-values from a single call.** Note what Moran's R measures: not
"is co-expression high" but "are ligand and receptor *spatially arranged* together" — a pair that
is ubiquitous and uniform scores low by design.

### Spatial model

Gaussian kernel `exp(-d²/(2·bandwidth²))`, weights below `cutoff` zeroed, on a
`max_neighbours`-capped KNN graph. Effective radius = `bandwidth × sqrt(-2·ln cutoff)`.

**Ours: bandwidth 13.1454 µm, cutoff 0.1 → support radius 28.2 µm**, **median 14** neighbours per
cell, **max 52** (verified over all 100,197 cells). Derivation — an **equal-area correspondence to
a 50 µm square patch**:

| step | value |
|---|---|
| `k = sqrt(-2·ln cutoff)`, cutoff = 0.1 | 2.14597 |
| equal-area disk radius of an s×s patch, `s/√π` | **28.2095 µm** ← support radius |
| `bandwidth = R / k = s / 3.804` | **13.1454 µm** |

This is the **smallest neighbourhood of any method here** (SpatialDM 135 µm, CytoSignal 200,
stLearn 250, COMMOT 365), appropriate to single-cell resolution — the relevant range is tens of µm,
whereas on Visium it is ~100 µm because each unit is a 55 µm multi-cell spot. The `max_neighbours`
cap was left at LIANA's default 100 because it is **already non-binding** (max 52); raising it
would have been a gratuitous deviation.

> ⚠️ **Provenance, stated plainly because a reviewer will probe it.** The 50 µm is ALARMIST's patch
> edge length. Only the *geometry* is borrowed — an area-preserving square→disk conversion — and no
> ALARMIST output enters at any point. It nevertheless means **LIANA's spatial scale is set by
> reference to the method it is being compared against**, which is neither LIANA's own default nor
> the tissue's. This is contract deviation **CD-1**, still ❌ OPEN. Anyone reporting these results
> must say so.

**The exploration step cannot settle the bandwidth on this tissue** (`choose_bandwidth.py` →
`GBM/bandwidth_choice/`, 2 figures + 4 CSVs; outcome: the number **does not move**):

- **`li.ut.query_bandwidth` returns `ceil(MEDIAN) − 1`, not the mean**, despite the variable being
  named `avg_nn` (`query_bandwidth.py:71-72`). This reconciles three numbers that otherwise
  disagree: a BallTree **mean** of 14.65 neighbours at R = 28.21, the **13** read off the curve
  there, and the manifests' **median 14**. `ceil(14) − 1 = 13`.
- **⚠️ Its x-axis is a HARD query radius, not a σ.** `query_bandwidth` uses `BallTree.query_radius`,
  whereas `spatial_neighbors(bandwidth=…)` takes a Gaussian **σ** truncated at `cutoff`, reach =
  σ × 2.145966. **A value read off that curve must be divided by 2.146 before being passed as
  `bandwidth=`**; passing it straight through inflates the neighbourhood **area by 4.6×**. This is
  an inconsistency in the tutorial itself. Our figures do not inherit it —
  `run_inflow_downstream.py:163-167` draws its guide at the **support radius**, correctly.
- **The curve has no plateau, elbow or inflection between 5 and 120 µm**, so the exploration step
  selects nothing. Mean neighbours per cell: 1.75 (R=10), 7.49 (20), **14.65 (28.21)**, 28.88 (40),
  59.16 (57.94), 85.21 (70), 240.06 (120).
- **Cell spacing gives a floor, not a choice** — median nearest-neighbour distance **7.86 µm**
  (IQR 6.30–10.47), so a strictly juxtacrine reach is ~8–10 µm (σ ≈ 3.7–4.9).
- **There is no characteristic signalling length scale in this tissue.** Pooling the LRIC radial
  profiles over **all 1,088 resolvable LR pairs × 13 punches = 11,795 pair×punch observations**,
  the median *g(r)* by annulus is **1.459** (0–50 µm) declining to **1.395** (200–225 µm) — a
  **4.4% decline across the entire range**. The tutorial's biological criterion ("reflect the
  typical range of molecular signaling") has nothing to bind to.
- **`max_neighbours=100` supplies the upper bound**: at R = 70 the mean is already 85.2. The
  defensible window is R ∈ [20, 58], σ ∈ [9.3, 27].
- **The alternative not taken:** `inflow_mofaflex.ipynb`'s own value is **σ = 27 µm → R = 57.94 µm**.
  Ours is half that spatial scale, a quarter of the area.
- **The bandwidth and the QC attrition are coupled.** A wider kernel would by itself have prevented
  much of the MOFA-Flex view loss below: mGAM reachability rises **0.319 → 0.712** going from
  R = 28.21 to 57.94, Vascular 0.157 → 0.417.

### LR database

Default is LIANA's **consensus** resource (4,624 unique pairs). `li.mt.bivariate` also accepts
`resource: pd.DataFrame` with `['ligand','receptor']`, and LIANA joins heteromeric subunits with
`_` — exactly our CellChatDB v2 encoding, so the v2 handover is **direct and lossless**, with
native complex support. Pairs are named `<ligand>^<receptor>`; `^` and `_` mean different things.

3,218 unique v2 pairs → 27 self-interactions removed → **131 pairs** survive `nz_prop = 0.02`.

#### ⚠️ A silent join trap against the ALARMIST GBM run

**The LIANA runs consumed the re-exported (post-2026-07-28) `CellChatDBv2.0.human.csv`** (proven by
complex-subunit ordering: 450 × `TGFBR1_TGFBR2`, 0 × `TGFBR2_TGFBR1`; 633/633 distinct LR keys
match the current file, only 496/633 the `.old.csv`). **But `results/GBM/patch_lri_columns.csv` —
the ALARMIST run this is compared against — has 210 × `TGFBR2_TGFBR1` and 0 reversed**, i.e. the
**old** export. Neither run is wrong, but the two label heteromeric complexes in **different
subunit order**, so a raw-string join silently drops keys:

| join of LIANA's 633 LR keys against ALARMIST's 712 | shared |
|---|---|
| raw string | **495** |
| after canonicalisation | **545** |
| of ALARMIST's 205 heteromeric keys: raw / canonical | **94** / **144** |

**50 heteromeric keys (24% of them) vanish from a naive join** — silently, with no error.
**Remedy: canonicalise each side before joining**, applying `'_'.join(sorted(s.split('_')))` to the
ligand and the receptor *separately*:

```python
def canon(lr):                      # 'TGFB1^TGFBR2_TGFBR1' -> 'TGFB1^TGFBR1_TGFBR2'
    l, r = lr.split('^')
    return '_'.join(sorted(l.split('_'))) + '^' + '_'.join(sorted(r.split('_')))
```

**`GRN^SORT1` and `ANXA1^FPR1` are identical in both files**, so the two requested LRs are unaffected.

### Input

`layers['counts'] = X.copy()`; `sc.pp.normalize_total(target_sum=1e4)`; `sc.pp.log1p` — the
tutorial's exact recipe on genuine raw counts. ⚠️ **`run_liana.py`'s bivariate branch skips the QC**
that `run_inflow.py` applies, which is why that branch has 100,197 cells against the others'
100,190 — **cross-branch joins must be on cell ID, not position.**

### Branch map

| branch | entry point | output dir under `liana/GBM/` | wall | script |
|---|---|---|---|---|
| bivariate (spot-based) | `li.mt.bivariate` | `cellchatdb2/`, `default/` | 2.0 / 5.6 min | `run_liana.py` |
| bivariate, local = Moran's R | same, `--local-name morans` | `cellchatdb2_morans/` | 3.6 min | `run_liana.py` |
| **inflow** (single-cell) | `li.mt.inflow` | `cellchatdb2_inflow/`, `default_inflow/` | 1.0 / 1.3 min | `run_inflow.py` + `run_inflow_downstream.py` |
| NMF on either | `li.multi.nmf` | `nmf_{bivariate,inflow}{,_default}/` | ~0.9 min each | `run_nmf.py`, `plot_liana_full.py` |
| **MOFA-Flex** on inflow | `mofaflex` | `mofaflex_inflow/` (+`sensitivity_nzf0.001/`) | 76.0 min | `run_mofaflex.py` |
| MOFA-Flex, reachability-normalised QC | `--nzf-mode reachability` | `mofaflex_inflow_reachnorm/` | 41.5 min | `run_mofaflex.py` |
| **LRIC / cross-PCF** | `li.mt.lric`, `li.mt.cross_pcf` | `lric_percore/` | 2.6 min | `run_lric.py` |
| **LR-MISTy** | `li.mt.lrMistyData` + `misty` | `misty/linear_fullslide/` | 2.9 min | `run_misty.py` |
| factor annotation | `dc.mt.mlm` / `ulm` | `factor_annotation/` | 0.22 min | `annotate_factors.py` |
| vs ALARMIST | — | `vs_alarmist/`, `bandwidth_choice/` | — | `cosine_factors_vs_motifs.py`, `compare_programs_to_alarmist.py`, `why_no_mgam_motif.py`, `choose_bandwidth.py` |

**Whole-slide was verified safe rather than assumed** for the bivariate/inflow branches: a 28.2 µm
support against a measured 222.9 µm minimum inter-core distance gives **zero** cross-core pairs
(still zero at 200 µm), a **7.9× margin**. **That argument is per-branch** — LRIC's density
normalisation and LR-MISTy's 607 µm nominal reach are separate exposures; see their subsections.

### Branch: bivariate

**Workflow.** `li.ut.query_bandwidth` → `li.ut.spatial_neighbors(bandwidth, cutoff, kernel,
set_diag)` → `li.mt.bivariate(resource, local_name='cosine', global_name='morans', n_perms=100,
mask_negatives=False, add_categories=True, nz_prop, use_raw=False)`.

**Data outputs** (`cellchatdb2/data/`): `global_scores.csv` (131 × 10 — ligand/receptor, `*_means`,
`*_props`, **`morans`**, `morans_pvals`, `mean`, `std`), `local_scores.npz`, `local_pvals.npz`,
`local_categories.npz` (all cells × pairs, float32), `cell_meta.csv`.

**Image outputs** — 35 PNGs (27 top level + 8 in `requested/`): top-6 by Moran's R as
`local_*.png` / `local_*_pvals.png` / `local_*_cats.png` / `genes_*.png`, plus `top_morans.png`,
`bandwidth_query.png`, `connectivity_idx57404.png`. Requested pairs are segregated into
`plots/requested/` per `SKILL.md:56-61`. Not produced, recorded rather than run: the
`mask_negatives=True` re-plot (reproducible from the persisted matrices with no re-fit) and the
cell2location / decoupler MuData route (**N/A** — single-cell Xenium has neither modality).

**Top 8 by Moran's R:** COL4A1^CD44 (0.1155), COL4A2^CD44 (0.1051), APP^TNFRSF21 (0.0900),
DLL3^NOTCH1 (0.0895), C3^C3AR1 (0.0870), TNC^SDC4 (0.0752), HLA-DQA1^CD4 (0.0739), IGF2^IGF2R
(0.0694) — **ECM–CD44 adhesion dominates the head of the list.**

| LRI | rank | Moran's R | p | ligand prop | receptor prop |
|---|---|---|---|---|---|
| GRN^SORT1 | **33 / 131** | 0.035702 | 0.0 | 0.109 | 0.188 |
| ANXA1^FPR1 | 60 / 131 | 0.013804 | 0.0 | 0.121 | 0.029 |

Both **spatially significant** but with modest effect sizes. GRN→SORT1 sits at the **25th
percentile** — better than chance but firmly mid-table, and well below COMMOT's rank 1–4.

⚠️ **Merging the cores dilutes core-specific signal — measured.** A smoke run on cores 13+14 alone
ranked **`SLC17A7_GLS^GRIA1` (glutamatergic) first**. On the full 13-core slide that pair drops out
of the top entirely, replaced by ECM–CD44 adhesion. Nothing is wrong — Moran's R is computed over
all cores jointly, so a program strong in two cores is averaged against eleven where it is not.
This is the concrete cost of the merged design, and it is worth weighing against **SpatialDM and
COMMOT, both run per-core, independently identifying glutamatergic signalling as top in core 13**.
If core-level heterogeneity is of interest, LIANA+ should be re-run per core.

### Branch: local metric `cosine` vs `morans`

`cellchatdb2_morans/`, produced by `run_liana.py --local-name morans` (3.6 min, 35 figures, same
131 pairs, same resource, same bandwidth). `local_name` is the **only** thing that changed.

**The global statistic is unchanged; the local one is a different object entirely.** `global_name`
was already `morans`, so global Moran's R, its p-values and the whole pair ranking are
**bit-identical** across the two runs (max |Δ| = 0.0 over all 131 pairs). Only `lrdata.X` differs:

| | `cosine` | `morans` |
|---|---|---|
| range over 100,197 × 131 | **[0, 1.0000002]** | **[−14.084, 73.273]** |
| zeros | **83.6214%** | **0.0000%** |
| negatives | **0** | **4,378,791 (33.4%)** |

Per-pair Pearson correlation between the two local matrices: **median 0.195** (range 0.095–0.337);
`GRN^SORT1` 0.189, `ANXA1^FPR1` 0.253.

**"Is this interaction spatially structured overall" is robust to the local metric; "which cells
carry it" is almost entirely determined by it.** Any figure, niche assignment or downstream
factorisation built on the *local* scores is a statement about the chosen metric as much as about
the tissue.

Two consequences. **Moran's R is not NMF-admissible** (33.4% negative, none zero), so the
communication-program branch could not have been run on it — a *post hoc* justification for
`cosine`, not the original reason, and recorded as such. And **the cross-punch reproducibility
filter is only correct for a non-negative local metric**: it counts a feature as present in a punch
if non-zero, implemented as `X > 0`. For `cosine` (0 negatives) `X > 0` and `X != 0` coincide; for
`morans` they do not, and a third of the matrix would be silently treated as absent.

⚠️ The `requested/` filenames in this tree still carry **rank 33 / 60** — global Moran's R sets the
ranking and is unchanged. **Do not read those ranks as a `morans`-*local* result.**

**Cross-method:** all 131 LIANA pairs appear in SpatialDM's union of 1,662 across the 13 cores.
LIANA's whole-slide Moran's R vs SpatialDM's median per-core z gives **Spearman 0.616, Pearson
0.588**. ⚠️ **The caveat is not optional:** this is agreement between two implementations **at
their own default spatial scales** — SpatialDM per-core at ~709 neighbours/cell, LIANA whole-slide
at 28.2 µm support with median 14 — because `SKILL.md:45-46` forbids harmonising them. The residual
disagreement mixes *implementation* and *spatial scale* and **cannot be attributed to either**. An
r of 0.62 across a ~50× difference in neighbourhood size is the headline; it is not an
implementation-agreement coefficient.

### Branch: inflow — and why it is the right branch

LIANA's README decision tree routes **Spatially-resolved → Single-cell → Interaction scoring →
Inflow Score**, while `li.mt.bivariate` sits under the **Spot-based** branch. Xenium is
single-cell. The package encodes the same judgement in its defaults: `_inflow.py` sets
`nz_prop = 0.001` against `bivariate`'s `0.05` — **the authors addressed single-cell sparsity by
writing a different method, not by asking users to lower a threshold.**

| | `li.mt.bivariate` | `li.mt.inflow` |
|---|---|---|
| branch | spot-based | **single-cell** |
| `nz_prop` | 0.05 (we used 0.02) | **0.001** |
| features | **131** | **4,608** (35.2×) |
| unique LR interactions | 131 | **633** (4.83×) |
| feature identity | LR pair | **sender cell type × LR pair** |
| sparsity | 83.6214% zeros | **99.4538%** zeros |
| runtime | 2.0 min | 1.0 min |

**The inflow feature space is ragged, not a grid.** 4,608 = **9 senders × 633 distinct LR pairs**,
with **250 to 602 LR pairs per sender** (Lymphoid 250, non-mGAM 382, mGAM 508, Vascular 520,
Glial-Neuronal 559, AC-like 590, OPC-like 596, NPC-like 601, MES-like 602) — each sender keeps only
the pairs that clear `nz_prop` for *that* sender.

Inflow carries **sender identity inside the feature** —
`Inflow_{i,s,l,r} = (Σ_j W_ij L_{j,l} C_{j,s} / Σ_j W_ij) · R_{i,r}` with `C_{j,s}` a hard
cell-type indicator — so a feature reads as "cell type *s* sending ligand *l* to receptor *r* here".
`bivariate` has no such axis.

**Both of the tutorial's optional spatial filters are no-ops on this data, and for the same
reason.** The SVG gene pre-filter: 5,097 → 5,097 genes. The SVI interaction filter: 4,608 → 4,608,
reproduced twice, 19.4 s over all 100,190 cells. In both cases **the `I > 0.01` half carries it**
(genes min I = 0.159; interactions min I = 0.0578, median 0.320) while **the FDR half is vacuous** —
all `pval_norm` values are **exactly 0.0**, which at n = 100,190 is floating-point underflow of the
normal approximation, not evidence. **Do not quote the FDR half as if it discriminated.** The SVI
filter was left off deliberately (pre-selecting the feature space for spatial structure would
confound the NMF comparison) and is exposed as `--svi-filter`.

**Downstream (`run_inflow_downstream.py`, 1.8 min, 80 figures).** The tutorial continues well past
`li.mt.inflow`, through `compute_global_specificity(groupby='cell_type')` → `groupby=<region>` →
`cell_type::region` composite → `li.ut.spatial_pair_proximity` → `li.mt.rank_aggregate`. All five
now run. `compute_global_specificity` on 100,190 × 4,608 at `n_perms=1000` takes **8 s**.

**Requested LRIs — inflow corroborates ONE arm of the loop, not both.** `global_interactions`:
41,472 source×target×LR rows, 5,417 at p<0.05; 616 LRs significant somewhere (of 633 — 17 are
significant nowhere). p is permutation-based at `n_perms=1000`, so **p = 0.000999 is the floor**.
**All of these p-values are cell-level (n = 100,190), not per-core — see CD-2.**

| ALARMIST motif-1 arm | `lr_mean` | p | rank within that LR (of 81 pairs) |
|---|---|---|---|
| **ANXA1→FPR1, MES-like → mGAM** | 0.1263 | **0.000999** (floor) | **2 / 81** |
| **GRN→SORT1, mGAM → MES-like** | 0.0202 | **1.0** | 20 / 81 |

**The MES-like→mGAM arm is strongly corroborated; the mGAM→MES-like arm is not corroborated at
all.** As whole LRs both rank mid-table among the 616 — GRN^SORT1 54th, ANXA1^FPR1 72nd — so
neither is a headline interaction for the method; the corroboration is specific to the *directed
cell-type pair*, which is exactly what ALARMIST claims.

⚠️ **Caveat that limits how far this can be pushed.** For both LRs the ranking is dominated by
**self-self pairs** — GRN→SORT1's top four are Vascular→Vascular, mGAM→mGAM,
Glial-Neuronal→Glial-Neuronal, non-mGAM→non-mGAM; ANXA1→FPR1's top is mGAM→mGAM. The inflow score
multiplies a neighbourhood-averaged ligand term by the receiver's own receptor expression, so a
cell type that both expresses the ligand and clusters with itself scores highly without any
cross-type signalling. **Read the off-diagonal cells; treat the diagonal as a co-localisation
baseline.**

⚠️ **`lr_ranking_by_lr_mean.csv` is a *significance-conditioned maximum*, not a ranking by
`lr_mean`.** It computes `gi[gi.pval < 0.05].groupby('lr')['lr_mean'].max()` — a maximum over the
81 source→target pairs. Hence 616 rows against 633 distinct LRs. Its columns are now self-documenting
(`lr, max_lr_mean_over_signif_pairs, n_signif_pairs`) and an unconditional companion
**`lr_ranking_all_pairs.csv`** (633 rows) sits beside it; both put GRN^SORT1 54th and ANXA1^FPR1 72nd.

### Branch: NMF communication programs

**NMF is demonstrated only in `bivariate.ipynb`; `inflow_score.ipynb` has no NMF section.** So
NMF-on-inflow is **our composition**, not an author-demonstrated path — the decision tree's own
answer for unsupervised decomposition at single-cell resolution is MOFA-Flex (next subsection).
Both NMFs are kept: on `bivariate` as the tutorial-sanctioned composition, on `inflow` as the
resolution-appropriate input.

**Choosing `k_range`.** `li.multi.nmf` defaults to `range(1, 11)`. **Kneedle normalises the curve
over the range it is given, so a rank obtained on 1..11 and one on 1..41 are not comparable
numbers.** Three configurations were run; **`range(1,21)` is in force**, because both branches were
fitted on the *same* window: bivariate rank **6**, inflow rank **7**, neither at a boundary. The
runner warns when a rank lands on an endpoint. NMF cost is not the constraint (fits get *faster*
with k: 18.7 s at k=5, 10.8 s at k=40).

⚠️ **Never quote "rank 6" / "rank 7" bare — the elbow is weakly supported.** The elbow metric is
**MAE**, and on a matrix that is 83.6% / 99.45% zeros the zero-predictor is a formidable baseline:

| branch | rank | relative Frobenius error | fraction of SS captured | achieved MAE | zero-predictor MAE |
|---|---|---|---|---|---|
| bivariate | 6 | **0.7607** | **42%** | 0.080583 | 0.072522 |
| inflow | 7 | **0.7771** | **40%** | 0.015287 | 0.011817 |

**Neither fit beats the trivial all-zero predictor on the metric the elbow is computed with.** That
is a property of an L1 elbow on an extremely sparse non-negative matrix rather than a defect in the
run — but **the rank must always be reported with `rel_frob` and the zero baseline alongside it,
not as a count of programs the tissue "has".**

**How much does the program structure change? — a lot.** Spatial correlation of the two `NMF_W`
matrices over the 100,190 shared locations: best match **r = 0.43**, and **5 of 7 inflow factors
have max |r| < 0.3 to any bivariate factor**. Three separate bivariate factors map onto the same
inflow F3, so the mapping is not even one-to-one. **The decompositions are not nested** — bivariate
is not a coarse summary of inflow, it is a different partition. The organising axis differs
qualitatively:

| | top features per factor |
|---|---|
| **bivariate — by LR family** | F1 JAG1^NOTCH1/2 · F2 immune/complement · F3 ECM (COL4A1/A2^CD44) · F4 TNC^SDC4, APP^SORL1 · F5 neurexin · F6 DLL3/DLL1^NOTCH1/2 |
| **inflow — by sender cell type** | F1 Glial-Neuronal · F2 OPC-like · F3 MES-like · F4 NPC-like · F5 AC-like · F6 Glial-Neuronal · F7 NPC-like |

**Every one of inflow's 7 factors is anchored to a single sender cell type**; not one bivariate
factor is. Inflow also splits the NOTCH ligands by sender where bivariate collapses them — a genuine
resolution gain, not a re-labelling.

**The finding to carry forward: the recovered "communication programs" are strongly determined by
the upstream expression filter and feature construction, not only by the tissue.** A reviewer shown
only the bivariate NMF would conclude six programs organised by ligand family; shown only inflow,
seven organised by sender. Both are defensible runs of the same package. **Any claim about
"communication programs" from this class of method must state the upstream configuration.**

**Sensitivity: the rank moves with the bandwidth.** Both rows fitted at `k_range = range(1,41)`, so
bandwidth is the only difference: at 18.75 µm support 40.2 µm, inflow rank **11**; at 13.1454 µm
support 28.2 µm, inflow rank **7** — a 30% reduction in support radius took the recovered program
count down 36%, while bivariate stayed at 3 throughout (pinned by its 131-feature ceiling, itself
an artifact of `nz_prop`). **Neither rank is a property of the tissue without stating the bandwidth.**

**Outputs.** `nmf_{bivariate,inflow}/`: `elbow.png`, `factor_maps.png`, `NMF_W_factor_scores.csv`,
`NMF_H_loadings.csv`, `nmf_WH.npz`, `nmf_errors.csv`, `top10_loadings_per_factor.csv`,
`punch_presence.csv`, plus `plots_full/{factors,global,interactions}/`, and
`nmf_factor_correlation.csv` + `nmf_error_vs_k.png` at the GBM level.

⚠️ **`plots/` and `plots_full/` are DISJOINT, and `plots_full` does not mean "the full set".**

| directory | written by | contents |
|---|---|---|
| `nmf_*/plots/` | `run_nmf.py` | exactly **2** files: `elbow.png` and `factor_maps.png` |
| `nmf_*/plots_full/` | `plot_liana_full.py` | **30** files (bivariate: 25 PNG + 4 CSV + manifest) / **62** (inflow: 55 PNG + 6 CSV + manifest) |

**`plots_full` means "the output of `plot_liana_full.py`" — a second, *additive* pass.** Not one
filename appears in both. **`elbow.png` and `factor_maps.png` exist only under `plots/`; deleting
`plots/` as "redundant" would destroy the entire rank-selection evidence.** Both visualise the same
fit — verified: same `nmf_WH.npz`, factor counts agree across manifests and panel counts, mtimes
place `plots_full` strictly after each fit with no refit between. **There is no stale-figure problem.**
(The name is also overloaded across this document — stLearn uses `plots_full` to mean "went beyond
the tutorial".)

⚠️ **Reading traps in this tree**, none of which changes a result:

- **Two `global_interactions.csv` with different p-values.** The copy under `nmf_inflow/plots_full/global/`
  was written at `--n-perms 100`, **unseeded**; `cellchatdb2_inflow/data/` at `n_perms=1000, seed=1337`.
  Both 41,472 rows with identical `lr_mean`, but **393 rows flip across p < 0.05** (5,290 vs 5,417).
  **`cellchatdb2_inflow/data/global_interactions.csv` is canonical**; the other is now named
  `global_interactions.SUPERSEDED_nperms100_unseeded.csv`.
- **`plots_full/global/interaction_total_inflow.csv` has 4,608 rows** inside a tree whose NMF used
  only the 2,704 post-punch-filter features. Verified that no co-indexed array mixes the two spaces
  and all five plotted top interactions are in the kept set — **a latent hazard, not a current
  error**. The bivariate copy is also mis-named: its `total_inflow` column holds summed **cosine**.
- **`nmf_factor_correlation.csv` has no producing script in the repo.** `grep -rn` matches only this
  document. The **values are correct** (the stored 6×7 matrix reproduces the cell-ID-aligned
  correlation; max |r| = 0.4299 matches the 0.43 above) — **only the code is missing, so the table
  cannot currently be regenerated.** Its GBM-level neighbour `nmf_error_vs_k.png` *does* have one
  (`plot_nmf_errors.py`), which performs no stability selection or rank test LIANA does not itself perform.
- **`pair_proximity.csv`'s `interacting` column is `1` for all 81 pairs** — no information, expected
  at the deliberately coarse 100 µm proximity bandwidth used for that call (distinct from the
  13.1454 µm scoring kernel). **The usable output is the continuous `proximity` column** (0.0435–0.9933).
- **`inflow_lrdata.h5ad` still carries ALARMIST outputs** — `obs['motif']` (15 categories),
  `obs['patch_id']`, `uns['motif_colors']`, plus orphaned `obsp['connectivities']`/`['distances']`.
  **Nothing in the LIANA pipeline ever read them** (`run_inflow_downstream.py` hard-refuses them as
  `--region-col`), so no result is contaminated — but it is a **leakage trap** for any future script
  that iterates `adata.obs` keys generically. **Drop `motif` and `patch_id` at load.**

**Two plotting traps, both hit here.** `li.multi.nmf`'s `_plot_elbow`, `li.pl.dotplot` and
`li.pl.connectivity` return **plotnine** objects (save with `.save()`; `plt.gcf()` captures
*nothing* after them), while `li.pl.feature_by_group` returns a **matplotlib `(fig, ax)` tuple**.
The three entry points disagree on return type, so a single save helper does not work. This
produced four blank all-white PNGs across three scripts before being fixed; all now route plotnine
returns through a dedicated `save_gg` helper and refuse to write a figure with no drawn content.
**Colour scaling on the inflow score maps** uses `vmin=0, vmax='p99.5', sort_order=True`, because
the matrix is 99.45% zeros with a long tail and full-range scaling renders an all-black map;
percentile scaling is the tutorial's own idiom.

### Branch: MOFA-Flex on inflow (the authors' prescribed unsupervised route)

`run_mofaflex.py` → `mofaflex_inflow/`, following `inflow_mofaflex.ipynb` cell by cell. The
decision tree's answer for *Single-cell → Unsupervised* is **Communication Programs = Inflow +
MOFA-Flex**, not NMF.

⚠️ **CD-1 applies to this branch.** MOFA-Flex re-specifies no spatial model of its own: it
factorises `cellchatdb2_inflow/data/inflow_lrdata.h5ad`, so all spatial weighting already happened
inside `li.mt.inflow` at the ALARMIST-derived 13.1454 µm bandwidth. Unlike LRIC (own `cKDTree`) and
LR-MISTy (tutorial `bandwidth=200`), which CD-1 exempts, **this branch inherits the deviation** —
record it wherever a MOFA-Flex factor count is quoted.

#### ⚠️ The single most reportable result: the authors' own QC deletes both arms of ALARMIST motif 1

Two consecutive tutorial cells do it, and **neither errors**:

| tutorial cell | filter | effect on this dataset |
|---|---|---|
| 19 | `nonzero_fraction > 0.01` | **4,608 → 447 features.** Removes **every** `ANXA1^FPR1` feature: its maximum non-zero fraction across all 9 senders is **0.009422** (MES-like); mGAM is 0.004911. All 9 fail |
| 23 | `lrdata_to_mudata(min_features=25)` | drops the **entire mGAM view**, which retained exactly **24** features — **one short**. (`non-mGAM` also drops, at 4.) So `mGAM^GRN^SORT1` is gone even though its non-zero fraction, 0.031770, comfortably cleared cell 19 |

**Run exactly as its authors demonstrate, LIANA+ structurally cannot see either arm of motif 1 on
this dataset.** A finding about the comparator, not about the biology.

**And the cell-19 cut is worse than "a strict threshold" — it is an abundance filter in disguise.**
An inflow feature can only be non-zero where its sender is *reachable*, so a single global threshold
imposes a different effective stringency per sender: surviving-feature count vs sender abundance is
**Spearman ρ = 0.917** (p = 5.1e-4, n = 9). Reachability at R = 28.2096 µm: NPC-like **0.761**,
AC-like 0.750, OPC-like 0.748, MES-like 0.727, Glial-Neuronal 0.379, mGAM **0.319**, non-mGAM 0.205,
Vascular **0.157**, Lymphoid **0.025**. **Lymphoid cannot pass the cut for arithmetic rather than
biological reasons** — its largest `nonzero_fraction` is 0.004332, below 0.01, so **zero** Lymphoid
features can survive whatever the biology is. Expressed as a fraction of the cells that could
*receive* it, the same 0.01 cut is **1.3%** for NPC-like and **39%** for Lymphoid. **Vascular is
penalised for clustering, not rarity** — 3.2% of cells but only 15.7% reachable, because vessels
are spatially aggregated.

#### The three fits

| | tutorial QC (`global`) | sensitivity `nzf>0.001` | **reachability-normalised** (ours) |
|---|---|---|---|
| features after cut | 4,608 → **447** | 4,608 → **1,550** | 4,608 → **779** |
| views ≥ `min_features=25` | **6 of 9** | 8 of 9 | **9 of 9** |
| mGAM / non-mGAM / Lymphoid | 24 → dropped / 4 → dropped / 0 → absent | — | **66 / 34 / 39** |
| Vascular | 26 | — | **138** (5.3× recovery) |
| active factors | **17 / 20** | **20 / 20** | **19 / 20** |
| epochs (cap 1000) | 632 | 199 | 294 |
| fit / wall | 70.5 / 76.0 min | 33.7 / 35.3 min | 40.1 / 41.5 min |
| peak RSS | 4.1 GB | 8.3 GB | 6.78 GB |
| punch-level grade test | **0 / 17**, min raw p 0.013986 → q 0.218 | **0 / 20**, min p 0.008159 → q 0.163 | **0 / 19**, min p 0.022145 → q 0.332 |

`nzf_norm = nonzero_fraction / reach[sender]`, keep `> 0.01` — "is this feature non-zero in an
appreciable share of the cells that *could* receive it", which is scale-free in sender abundance.
**The default stays `--nzf-mode global`, so the tutorial-faithful run remains reproducible; the
reachability fit is an additional fit, not a replacement.** The SVI filter is a no-op on it too
(779 → 779, measured). Primary fit: max single-view R² 0.1681, max |r| between factor scores 0.243
— **zero** pairs exceed the tutorial's 0.6 redundancy flag. A 20-epoch determinism probe on the
primary is bit-identical (`max_abs_weight_diff: 0.0`); no probe was run for the reachability fit.

R² per sender view on the reachability fit — **and the view the tutorial deleted is the second best**:
Glial-Neuronal **0.5058**, **mGAM 0.4453**, NPC-like 0.3931, MES-like 0.3542, OPC-like 0.3254,
Vascular 0.2708, AC-like 0.2629, non-mGAM 0.1397, **Lymphoid 0.0013**.

⚠️ **Caveat that must travel with the reachability run: Lymphoid is admitted but uninformative.**
R² 0.0013, 39 features resting on 2.5% reachability, largest single factor×view cell 0.00021. The
criterion lets it in; the model explains essentially nothing there. **Do not interpret Lymphoid
factors.** No absolute cell-count floor was applied, so nothing else guards against a thin cell base.

#### ⚠️ `n_factors = 20` is a BINDING CEILING, not a selection

**MOFA-Flex does not choose K.** `n_factors` is a hard ceiling: the model fits exactly that many
and inactive ones are pruned *afterwards* by a 2%-R²-in-one-view floor. The 20 comes from
`inflow_mofaflex.ipynb` cell 25. **It coincidentally equals ALARMIST's K = 20; that coincidence is
not the justification and must never be presented as one.**

The ceiling **binds in all three fits** and there is no taper at the bottom of the list — the three
weakest active factors hold 6.7% / 8.3% / 7.7% of total R², and the reachability fit's weakest
factor clears the 0.02 floor by 44% (0.0287). **So 17 / 19 / 20 are artefacts of the argument, not
discovered numbers**, and none may be compared to ALARMIST's K = 20 as though both were fitted.
Same failure class as the `k_range` confound. `run_mofaflex.py` now warns when ≥ 90% of factors are
active and records `n_factors_requested` / `n_factors_active` / `ceiling_binding`. **A K-sweep
(20/30/40/60) has NOT been run**, so where the active count saturates is unknown.

*(For symmetry: `results/GBM/analysis_parameters.csv` records no selection criterion for ALARMIST's
K = 20 either. That is a statement about what is recorded, not a demonstration that the two choices
are equally unprincipled — nothing further was checked.)*

#### Where the two motif-1 arms land

Sensitivity fit (`nzf>0.001`, 1,541 features): `mGAM^GRN^SORT1` peaks on **Factor 19** (+0.651,
rank 73/1541), `MES-like^ANXA1^FPR1` on **Factor 7** (−0.278, rank 136/1541). **Different factors —
the GP-prior factorisation does not put the two arms of the loop on one program.** Factor 7 does
carry three of the four loop features with concordant sign and is myeloid-anchored (largest weight
`mGAM^C3^C3AR1` at −2.892) — but Factor 7 is only the **eighth** strongest factor for
`mGAM^GRN^SORT1` (ordered by |weight|: F19 0.651, F15 0.590, F1 0.558, F16 0.556, F17 0.531,
F9 0.451, F6 0.413, F7 0.398). **The concordant signs are real; the preference is not.**

Reachability fit (779 features), ranks within-factor by |weight|:

| feature | peak factor | weight | rank | next strongest |
|---|---|---|---|---|
| `mGAM^GRN^SORT1` | **Factor 19** | **+0.329** | **67 / 779** (top 8.6%) | F7 −0.327 (80), F9 −0.300 (62), F16 −0.286 (119), F1 +0.279 (99) |
| `MES-like^ANXA1^FPR1` | **Factor 1** | **+0.078** | **276 / 779** (top 35.4%) | F7 −0.071 (305), F15 −0.054 (308) |
| `MES-like^GRN^SORT1` (reverse) | Factor 4 | +0.258 | 121 / 779 | F12 +0.148 (169) |
| `mGAM^ANXA1^FPR1` (reverse) | Factor 1 | +0.165 | 164 / 779 | F7 −0.150 (176) |

1. **`mGAM^GRN^SORT1` is spread flat across factors with mixed signs** — its top **six** factors lie
   within **0.057** of each other and their signs alternate — so **no factor claims it**.
2. **`MES-like^ANXA1^FPR1` is nowhere**: roughly two-thirds of all features load more strongly on
   its own best factor than it does.
3. **The autocrine-ish `mGAM^ANXA1^FPR1` outranks the biologically meaningful `MES-like^ANXA1^FPR1`
   on 17 of the 20 factors.**

**Neither arm is a top-10 feature of any factor**, so neither ever appears in `top_weights.png`. The
best joint factor by worst-of-the-two rank is Factor 1 (99 and 276) — and all four best candidates
carry the **same sign** for both arms (15 of 20 factors do). **The model is not placing the two arms
at opposite poles; it simply places neither anywhere prominent.** **Restoring the mGAM view changed
neither the co-loading nor the grade result.**

#### ⚠️ What the three MOFA-Flex figures actually show — read before quoting one

All three were read out of the installed `mofaflex/pl/_plotting.py`, because two are easy to
describe wrongly.

**`top_weights.png`** — the x axis is **`| Weight |`**, the *absolute* loading. **Sign is carried
only by the glyph** (⊕ / ⊖) and `scale_shape_manual(..., guide=None)` means **no legend for it is
ever drawn**; a reader who has not been told this will read every bar as positive. Within each facet
the top *n* are taken by |weight| then sorted **ascending**, so the **largest sits at the TOP**.
`facet_wrap(scales="free")` — **every panel has its own x scale, so bar lengths are not comparable
across factors.** These are **raw** weights, no prevalence normalisation. Now at top-10.

**`variance_explained.png`** — **it is a `geom_tile` HEATMAP, not a bar chart.** Rows = factors,
columns = views = sender cell types, fill = R², one facet per group (there is a single group here).
Rows ordered by total R² descending with the **largest at the BOTTOM**. `r2_per_view.csv` is exactly
the **column sum over factors** (verified to 4 d.p.). On the reachability fit the darkest cell is
Factor 5 × Glial-Neuronal = 0.1641, and the **Lymphoid column is effectively blank** (largest 0.00021).

**`circle_plot_Factor*.png`** — **the edges are factor-independent.** The factorisation never sees
the receiver. `source` comes from the feature name (model), and *which* 10 edges are drawn comes
from |loading| (model) — but **`target` and the edge weight come from
`inflow_means = lrdata.to_df().groupby(obs[cell_type]).mean()`**, i.e. the *receiving* cell's own
annotation and the **raw** mean inflow, not anything factor-weighted. Verified: **Factor 1 and
Factor 6 share 108 of their 225 (source, LR, target) rows with edge weights identical to 0.000e+00**,
while their loadings on those same rows differ by up to **1.9001**. So two factors that share a
feature display the *identical* sub-network, the edges ignore the factor's sign, and the apparent
"per-factor network" is mostly not the factor. **Read it as "the top-10 interactions this factor
selects, and where those interactions generally go", NOT as "this factor's sender→receiver
structure".** Now emitted for all active factors.

**Version gap — recorded, not guessed.** `inflow_mofaflex.ipynb` cell 5 states verbatim that it
"uses the MOFA-Flex `0.2.0` API, which is not yet released on PyPI" and instructs installing from
git main. The installed build is **`0.1.0.post2.dev179+g9792b435f`** from git main. **Every symbol
the notebook uses exists in the installed build with matching argument names**, so nothing had to be
guessed — but the version strings do not match.

### Branch: LRIC and cross-PCF (spatial co-occurrence)

`run_lric.py` → `lric_percore/`. The decision tree's **Single-cell → Spatial co-occurrence → LRIC**
branch, and **the one place in this benchmark where a method resolves direction (sender → receiver)
as an explicit argument** rather than as a feature name.

**Core algorithm.** `li.mt.cross_pcf` is a cross-type pair-correlation function: for an annulus
[r, r+w) it counts observed *A*-near-*B* pairs and divides by the count expected under complete
spatial randomness, giving g(r) with g = 1 meaning "no more co-located than chance". `li.mt.lric`
is the same construction restricted to the *ligand-expressing* members of *A* and the
*receptor-expressing* members of *B*. **The LRIC / cross-PCF ratio is therefore the whole point:
it asks whether expression adds anything over cell-type co-location.**

**Spatial model.** Annuli, no kernel and no connectivity graph — both functions build their own
`cKDTree`, so **`li.ut.spatial_neighbors` is never called and the repo's 13.1454 µm bandwidth does
not enter this branch at all** (CD-1 exempt). Realised bins: `[0,50) [50,75) … [200,225)` — 8 bins,
the first double-width because `extend_first_annulus` merges the contact band into it, and
`max_radius=200` is the **inner** edge of the last bin so the true reach is **225 µm**. Gates:
`min_cells=50` (a cell type is dropped from a punch) and `min_expressing=20` (an LR pair is NaN'd
for that directed type pair). Both bite.

**Why per punch — measured, not asserted.** Both functions normalise by density computed as
`n_points / bounding-box area`. On a TMA the global bounding box is mostly empty: global bbox
**123,358,214.8 µm²** against the sum of the 13 punch bboxes **52,094,949.7 µm²** — **42.2306%
occupancy**, so whole-slide density is understated by **2.3679×**. A pooled control was run anyway
(`--whole-slide-check`), and the **measured** inflation exceeds that lower bound: LRIC g(r)
**3.4284×**, cross-PCF **2.9445×**. **And it distorts the ratio, which is the quantity under test** —
whole-slide gives `GRN^SORT1` mGAM→MES-like ratios of 1.212–1.284 per bin and `ANXA1^FPR1`
1.737–2.261, against per-punch medians of 0.951–1.073 and 0.926–1.134. **A whole-slide run would
have manufactured a false positive — "expression adds ~2× over co-location" — for exactly the claim
being tested.** Note the contrast with bivariate/inflow, which *were* safe whole-slide: **density
normalisation is a different exposure and is not protected by the 222.9 µm argument.**

**Cost and support.** 13 punches at 816–26,456 cells, 3,418–5,092 genes, 299–1,087 LR pairs, 6–9
cell types kept. Per-punch LIANA call time **65.3 s**; whole script **155.7 s wall**, peak RSS
**17.017 GB** (dominated by the pooled 100k-cell control). mGAM clears `min_cells=50` in **11/13**
punches; `min_expressing=20` removes more — `GRN^SORT1` **8/13** both directions, `ANXA1^FPR1`
**8/13** (MES-like→mGAM) and **6/13** (mGAM→MES-like). **Every gate is attributed per punch** in
`combined/target_availability.csv` and `target_expression_support.csv` — a design feature: the
reason a number is missing is on disk next to the numbers that are present.

**Result 1 — direction asymmetry holds for GRN→SORT1, not for ANXA1→FPR1.** Punches, not cells, are
the replicate unit.

| LR | direction | n | median mean-g | IQR |
|---|---|---|---|---|
| `GRN^SORT1` | **mGAM→MES-like** | 8 | **1.566** | 1.424–1.668 |
| `GRN^SORT1` | MES-like→mGAM | 8 | 1.410 | 1.312–1.475 |
| `ANXA1^FPR1` | mGAM→MES-like | 6 | 1.524 | 1.459–1.626 |
| `ANXA1^FPR1` | MES-like→mGAM | 8 | 1.560 | 1.438–1.685 |

`GRN^SORT1`: forward > reverse in **7 of 8** punches, paired Wilcoxon **p = 0.0156** (the minimum
attainable at n = 8), median paired difference **+0.122**; within each punch's own valid pairs the
forward direction has median rank-pct 0.41 and the reverse 0.72. `ANXA1^FPR1`: median difference
**+0.013 in favour of MES-like→mGAM**, **p = 0.844**, 4/6 punches — **the single-punch asymmetry
does not replicate.**

**Result 2 — the LRIC / cross-PCF ratio is ~1 everywhere. This is the central finding, and it is
negative.**

| LR | direction | per-punch median ratio | IQR | Wilcoxon vs 1 |
|---|---|---|---|---|
| `GRN^SORT1` | mGAM→MES-like | **1.043** | 0.965–1.077 | p = 0.383 |
| `GRN^SORT1` | MES-like→mGAM | **0.963** | 0.900–1.004 | p = 0.109 |
| `ANXA1^FPR1` | mGAM→MES-like | **1.040** | 0.974–1.097 | p = 0.438 |
| `ANXA1^FPR1` | MES-like→mGAM | **1.033** | 0.999–1.058 | p = 0.313 |

Per bin the median ratio stays in **0.93–1.13** everywhere except one dip. **So across punches the
spatial co-occurrence LIANA reports for these pairs is fully explained by cell-type co-location plus
the underlying point pattern; ligand/receptor expression contributes nothing detectable.**

**This also refutes the "the 28.2 µm kernel missed a longer-range interaction" hypothesis.** The
ratio is ~1 in *every* annulus out to **225 µm** — 8× the bivariate/inflow support radius. There is
no longer-range LR-specific structure for a wider kernel to have found.

**Result 3 — the grade contrast is not testable, and the loss is systematic.** Of the 8 punches
supporting `GRN^SORT1`, **7 are high grade and exactly 1 is low grade**. The loss is not random:
the low-grade punches are precisely the ones where mGAM falls below `min_cells` (2, 6) or the
relevant expressing populations fall below `min_expressing` (4, 12, 14). **n = 7 vs 1 supports no
test, and none was run.** The cell-type-**agnostic** fallback (13/13 punches) is null:
`GRN^SORT1` p = **0.234** (median 1.556 high vs 2.053 low), `ANXA1^FPR1` p = **0.788**. **Do not
read the direction** — the low-grade punches pushing the agnostic median up are the small, sparse
ones, so this is more likely residual within-punch density structure than biology.

**Reproducibility.** Punch 3 was re-run inside the 13-punch set; the re-run is **bit-identical** to
the earlier benchmark (max absolute difference **0.0** across 58,744 summary rows, NaN pattern
identical). One off-by-one survives **in the benchmark's derived CSV only**:
`_benchmarks/lric_punch3/TARGETS_lric_vs_crosspcf.csv` lists `rank_by_mean = 42` for `GRN^SORT1`
mGAM→MES-like where the underlying summary ranks it **41 of 60** (verified, no ties). **The g values
are unaffected.**

**Outputs.** `combined/` (13 CSVs incl. `statistics.csv`, `target_availability.csv`,
`whole_slide_vs_perpunch.csv`), `punches/punch_<id>/` × 13, `whole_slide_check/` (labelled a
**control**, not a result, in the code, the manifest and the figure legend), and `figures/` —
**12 figures × png + pdf + svg**, the only LIANA branch that exports all three formats.

### Branch: LR-MISTy (multi-view modelling)

`run_misty.py` → `misty/linear_fullslide/`, following `misty.ipynb`'s own LR-MISTy configuration.

**Core algorithm.** `lrMistyData` builds two views over the same cells: an **intra** view of
receptor expression (the targets) and an **extra** view of *spatially-weighted neighbourhood ligand*
expression (the predictors). Each target is regressed on the extra view under 10-fold CV; the
statistic is the CV R². **It is a *predictive* model, not a test — no null, no p-value.**

**Configuration**, all tutorial or package values: `lrMistyData(bandwidth=200, set_diag=False,
cutoff=0.01, nz_threshold=0.1, kernel='misty_rbf')`, then `misty(bypass_intra=True,
model=LinearModel, k_cv=10, seed=1337)`. **Bandwidth 200 is the tutorial's, not this repo's
13.1454** — stated explicitly, because **this is the one LIANA branch that does not carry CD-1.**

| | value |
|---|---|
| cells | 100,197 → **100,190** after `filter_cells(min_genes=10)` |
| genes | 5,119 → **5,097** (21 control probes dropped, then `filter_genes(min_cells=3)`) |
| **intra view** | 100,190 × **382** receptor targets |
| **extra view** | 100,190 × **37** ligand predictors |
| connectivity | 10,119,190 nnz, degree **median = max = 100**, min 16, **cap binds for 99.7% of cells** |
| fit | **2.69 min** = 0.423 s/target × 382; **2.91 min** total; peak RSS **6.12 GB**; bit-identical across two runs |

⚠️ **Fact 1 — `bypass_intra=True` means `gain_R2` is not a gain.** `_Misty.py` sets
`intra_r2 = ... if not bypass_intra else 0`, so `intra_R2` is **0.0 for all 379 non-NaN rows** and
`gain_R2 == multi_R2` **exactly** — it is the CV R² of the ligand extra view *alone*, with nothing
subtracted. **The question "does the spatial extra view add anything over an intrinsic model" is
never computed by the tutorial configuration.**

⚠️ **Fact 2 — with one view the contribution column is 1.0 by construction.** The meta-model
coefficient is normalised to itself, so `contributions.csv` reads `extra = 1.0` for every target and
`plots/contributions_top.png` is **uninformative by design, not broken**.

**Target metrics — the hypothesised structural null does NOT hold.** Median `gain_R2` **0.00414**;
q75/q90/q95/q99 = 0.01539 / 0.03382 / 0.05109 / 0.11107; max **0.57675** (MET); **133/382 (34.8%)
above 0.01**; 20 above 0.05; 6 above 0.1; **57 exactly 0**; **3 NaN**. Top targets: MET 0.5767,
CNTN2 0.1779, CADM4 0.1212, GRIK3 0.1155, ADCYAP1R1 0.1098, CD44 0.1061. The **57 exact zeros** are
`clip(min=0)` on a *negative* mean CV R² — ultra-sparse receptors (SSTR3 at 0.22% non-zero, SELE
0.09%, …). The **3 NaNs** are heteromeric complexes where a complex is the elementwise **min** over
subunits, no cell co-expresses all three, and the variance is 0.

**But the signal is spatial context, not ligand–receptor mechanism.** Three pieces of evidence, the
third decisive: (1) **MET, the best target at R² = 0.577, is predicted without its own ligand** —
top predictors THY1 (t = −38.1), TNR (+32.0), TNC (+28.8), NTN1 (−27.0), CD99 (−25.0), while **HGF
is absent from the model entirely** (non-zero in only 3.77% of cells, fails `nz_threshold=0.1`).
(2) **Nothing is subtracted** (Fact 1), so `gain_R2` cannot separate "spatial ligand information"
from "anything correlated with position". (3) **The extra view is 37 ligands wide out of 3,218
pairs**, and the 100-NN cap binds for 99.7% of cells — what is being regressed is a heavily smoothed
100-cell neighbourhood average of 37 abundant genes, **a niche-composition proxy**.

**Net: LR-MISTy on this slide is a spatial-context predictor of receptor abundance, not an LRI
detector, and its `gain_R2` ranking must not be read as an LRI ranking.**

**The requested LRIs — a genuine partial positive.** `nz_fraction` gates whether a ligand is in the
model at all: GRN **10.88% ✅**, ANXA1 **12.06% ✅** on the full slide (against 7.27% and 2.52% on
punch 4 — **on that punch neither ligand cleared the threshold, so the earlier punch-4 negative was
uninformative**). On the full slide both clear it, so this is a fair test. Predictor ranks below are
by **|importance|** — MISTy importances are signed, and this asks "how much does this predictor
matter, in either direction"; the report README ranks the same rows by **signed** importance and
therefore gives different numbers (GRN 2nd and ANXA1 15th for FPR1, GRN 6th for SORT1). **Neither is
a correction of the other**; both are recomputed from `data/interactions.csv`, `view == 'extra'`.

| target | `gain_R2` | rank | GRN as predictor | ANXA1 as predictor |
|---|---|---|---|---|
| **SORT1** | 0.06210 | **17 / 382** (top 4.5%) | **t = +5.46, rank 9 / 37** by \|imp\| (6 signed) | t = **−3.50**, rank 11 / 37 (33 signed) |
| **FPR1** | 0.02615 | **54 / 382** (top 14%) | **t = +7.08, rank 3 / 37** by \|imp\| (2 signed) | t = **+1.34**, rank 22 / 37 (15 signed) |

Top-5 predictors of SORT1: SEMA4D +8.8, NTN4 −6.7, SEMA4C +6.5, PRNP +6.2, DHCR24 +6.2. Of FPR1:
C3 +14.1, APP −8.4, **GRN +7.1**, TGFB2 +6.5, DHCR24 +4.1. **GRN→SORT1 and GRN→FPR1 both come out
with strong positive importance and both targets sit in the upper half of the ranking. ANXA1→FPR1
does not replicate (t = +1.34, unranked), and ANXA1 enters SORT1 with the wrong sign.** Read against
the caveat above: **this is evidence that neighbourhood GRN tracks these receptors, not that a
receptor–ligand mechanism was detected.**

**Deviations from `misty.ipynb`:**

| # | Item | Tutorial | Ours | Number that forced it |
|---|---|---|---|---|
| **DEV-1** | HVG pre-step `adata[:, hvg]` | applied | **omitted** | Measured: full slide intra **100,190 × 82** / extra **× 13**, against 382 / 37 without it — **78% of receptor targets and 65% of ligand predictors discarded**. On a single 5,363-cell punch it collapses the extra view to **one** predictor, i.e. vacuous. The tutorial frames HVG as "for the sake of computational speed" on a genome-wide Visium slide; this is a 5,119-gene targeted panel and the full-panel fit takes 2.7 min |
| DEV-3 | bandwidth | 200 | **200 — the tutorial's** | this branch deliberately does not inherit 13.1454 µm, so **CD-1 does not apply** |
| DEV-4 | `max_neighbours` | not exposed by `lrMistyData` | package default 100 | cap binds for **99.7%** of cells, so the 200 µm bandwidth is largely inert. Tutorial behaviour, not a change |
| DEV-5 | preprocessing | tutorial QC | matches `run_inflow.py` | ⚠️ `run_liana.py`'s bivariate branch **skips this QC** — 100,197 vs 100,190 cells |
| DEV-6 | resource | `consensus` | CellChatDB v2 | consistency with the rest of the tree |
| DEV-7 | RandomForest secondary | run | **not run** | a 4-target full-slide probe took **5.85 min = 87.7 s/target** → 382 targets is **9.31 h**, ~9× the budget. The probe agrees with the linear model on ordering (MET 0.5950, CD44 0.1169, SORT1 0.0444, FPR1 0.0190) |

**Multi-sample / differential mode — none was run, deliberately.** LR-MISTy returns a
cross-validated R² **per target**, not a per-cell or per-core score, so there is no quantity to
aggregate to the 13 cores and no native contrast in the package. Fitting per core would give 13
incomparable models, each with its own `nz_threshold`-determined predictor set. Verified on disk:
the branch writes no grade or per-punch file at all. **This branch contributes no grade result by
design** — recorded because its silence could otherwise be mistaken for an oversight.

**Outputs.** `data/`: `target_metrics.csv`, `target_metrics_ranked.csv`, `interactions.csv`
(14,014 rows), `interactions_{SORT1,FPR1}.csv`, `contributions.csv`, `view_features.csv`,
`nz_fractions.csv`. `plots/`: **7 PNGs**. All three previously-unexercised `li.pl.*` MISTy plotting
functions ran clean against liana 1.8.1. ⚠️ **PNG only — no PDF/SVG**, unlike the LRIC branch.

### Branch: annotating the NMF factors (PROGENy + CellChatDB pathways)

`annotate_factors.py` → `factor_annotation/` (0.22 min, **no re-fit** — reads `NMF_H_loadings.csv`
only). `dc.mt.mlm` (the authors' choice in `mofatalk`) and `dc.mt.ulm` for robustness, decoupler
2.2.0, BH across pathways within each factor. Deviation from `mofatalk`: the LR universe is
CellChatDB v2 rather than `select_resource('consensus')`.

⚠️ **PROGENy is a poor fit here and must be reported as one.** It covers only **47/131 (35.9%)** of
bivariate and **160/524 (30.5%)** of inflow LR pairs. And most of what it returns is an artifact of
a single adhesion pair: of the inflow MLM hits at FDR < 0.05, **6 are WNT and 5 are p53, every one
with negative t** — and for F2/F4/F6 these are restatements of *"this factor loads `NCAM1^NCAM1`"*
(PROGENy weight `WNT ← NCAM1^NCAM1` = **−0.877**; those factors' NCAM1 shares are 23.5% / 39.9% /
30.5%). For F7 the driver is `DLL3^NOTCH1` instead, since **F7's NCAM1 share is exactly 0%**.
**Only the Hypoxia hits are defensible:** bivariate F4 (t = 3.55, FDR 3.2e-3), inflow F5 (3.42,
3.9e-3), inflow F3 (2.74, 2.7e-2).

**CellChatDB pathway composition carries the annotation instead — 100% coverage:** bivariate F1
NOTCH/JAG1 lateral induction (NOTCH 28.7%, ULM t = 4.01), F2 complement + galectin, myeloid-immune
(composition only, no FDR < 0.05 hit), F3 basement-membrane COLLAGEN (33.2%, t = 9.65, FDR
**2.3e-15** — the strongest call anywhere), F4 FGF–tenascin (FGF 24.4%, TENASCIN 11.8%, t = 5.03;
⚠️ the third-largest component is **APP at 10.5%, not THBS**), F5 synaptic NRXN (30.7%), F6
DLL3–NOTCH (NOTCH **55.5%**, t = 6.90); inflow F7 NOTCH (50.6%, t = 9.73, FDR 1.1e-17).
⚠️ **The CellChatDB column is named `pathway`, not `pathway_name`.**

⚠️ **Structural finding — the inflow NMF factors are mostly sender identity plus one adhesion pair.**
**Six of seven inflow factors are ≥ 75% a single sender** (F1 Glial-Neuronal 92.9%, F2 OPC-like
92.2%, F3 MES-like 82.5%, F4 NPC-like 86.4%, F5 AC-like 75.6%, F6 Glial-Neuronal 93.5%; only F7 is
mixed at NPC-like 63.2%). **And five of seven share the same top feature, `X^NCAM1^NCAM1`.** So
these are largely **sender identity plus one dominant adhesion pair, not distinct communication
programs** — a stronger statement of the same point the program-structure comparison makes.

**Neither branch co-loads the two arms of motif 1.** Inflow: GRN→SORT1 on F1 and F6 (both ≥92.9%
Glial-Neuronal), ANXA1→FPR1 on F3 (82.5% MES-like). Bivariate: F1 (loading 3.042, rank 7/131) and
F2 (0.851, rank 19/131). **And the strongest `GRN^SORT1` feature anywhere in the inflow
decomposition is `Glial-Neuronal^GRN^SORT1` — rank 9 of 2,704 on F1 — not the mGAM one.**

### The `default` tier — what the LR database does and does not confound

Driver: `run_default_tier.sh` (5 steps, one log each). Every parameter except `--db` is identical to
`cellchatdb2`. **The two resources are a real comparison, not a formality:** LIANA's `consensus` has
**4,624** unique pairs against CellChatDB v2's **3,218**, sharing only **1,663 — 36.0% of consensus,
51.7% of CellChatDB**. Both required LRs are present in both.

| | `cellchatdb2` | `default` (consensus) |
|---|---|---|
| bivariate pairs after `nz_prop=0.02` | 131 | **388** |
| inflow features | 4,608 | **9,448** (1,217 LR pairs × 9 senders, ragged 650–1,173) |
| inflow % zeros | 99.4538% | **99.3555%** |
| inflow features after ≥5/13 punch filter | 2,704 (58.7%) | **6,178 (65.4%)** |
| global-specificity rows | 41,472 → 5,417 at p<0.05 | **85,032 → 12,902** |
| NMF rank (bivariate / inflow), `k_range` 1..20 both | 6 / 7 | **6 / 5** |
| figures (bivariate / inflow) | 35 / 80 | **35 / 80** |
| NMF figures | 27 / 64 | **2 / 2** — see below |

#### ✅ The key result: the LR database is NOT a confounder for any per-interaction statistic

Both tiers joined on their shared entries and compared value-by-value:

| comparison | shared | max &#124;difference&#124; |
|---|---|---|
| bivariate global Moran's R | **79 pairs** | **0.000e+00** (means also identical) |
| bivariate `morans_pvals`, `mean`, `std` | 79 pairs | **0.0** each |
| inflow `lr_mean`, source×target×LR | **23,787 rows** | **0.0** |
| inflow `pval`, source×target×LR | 23,787 rows | **0 rows differ** |

Not luck — it follows from LIANA scoring each pair independently, with no cross-pair normalisation
anywhere. Swapping the resource changes exactly three things: **which pairs are tested**; **the
multiple-testing / ranking denominator** (`GRN^SORT1` moves from rank **33/131** to **77/388** with
an *unchanged* Moran's R of 0.035702; `ANXA1^FPR1` 60/131 → 187/388, R unchanged at 0.013804); and
**the feature space the factorisation sees**.

**Consequence for every rank quoted from this method: a LIANA rank is a statement about the
resource, not about the interaction. The score is not.**

#### The factor count is driven by the window and the bandwidth, not the resource

| upstream choice | measured effect on the inflow NMF rank | status |
|---|---|---|
| bandwidth 18.75 → 13.1454 µm (`k_range` 1..40 both) | **11 → 7** | clean — only bandwidth differs |
| `k_range` 1..10 vs 1..40 vs 1..20 (bivariate) | 3 / 3 / 6 | clean — only the window differs |
| LR database, `cellchatdb2` → `default` (`k_range` 1..20 both) | **7 → 5** (inflow); **6 → 6** (bivariate) | clean — refitted on the matched window |

**The database effect, once measured cleanly, is the weakest of the three**: it does not move the
bivariate rank at all (6 → 6 across a 3× change in feature count) and moves inflow only 7 → 5.
Both default fits still fail the zero-predictor check (rel-Frobenius 0.7582 / 0.8294, MAE 0.085256
and 0.016617 vs zero-predictor 0.077797 and 0.012344), exactly as the `cellchatdb2` fits do.

⚠️ **`nmf_*_default/` hold 6 data files and exactly 2 PNGs each, and no `plots_full/` tree** —
`plot_liana_full.py` was never run on the default tier. Against the `cellchatdb2` NMFs' 27 and 64
PNGs that is a 25 / 62 figure deficit, and it is a **deliberate scope limit, not a loss**: the
default tier exists to answer whether the resource confounds the statistics, and the data files
answer that.

⚠️ **Small provenance blemishes in the default tier — none affects a number.**
`default/run_manifest.json` records `"dataset": "default"` (the dataset is GBM; the string is the
output-dir name) and `resource_fingerprint: null`; all four default manifests record
`"tier": "cellchatdb2"` because the tier string is hardcoded in `run_nmf.py` /
`run_inflow_downstream.py` and was not parameterised — **the `resource` fields are correct**
(`"LIANA consensus"`, `resource_n_pairs: 4624`), so the runs are not misidentifiable;
`default_inflow/downstream_manifest.json` records `"db": ".../consensus"`, a path that does not
exist (the loader correctly branches on `a.db == "consensus"`). **Do not read that key as a file.**

### Relation to ALARMIST — the converging conclusion

**Four independent LIANA analyses agree that LIANA does not reconstruct the ALARMIST motif-1
mGAM ⇄ MES-like loop as a single program.**

| analysis | what it says about the loop |
|---|---|
| inflow global specificity | **ANXA1→FPR1 (MES-like→mGAM) rank 2/81, p = 0.000999 (floor); GRN→SORT1 (mGAM→MES-like) p = 1.00000.** One arm only |
| NMF (bivariate and inflow) | the arms land on **different factors with different dominant senders**; the strongest `GRN^SORT1` feature is Glial-Neuronal, not mGAM |
| MOFA-Flex | **different peak factors** (F19 vs F7) at `nzf > 0.001`; at the tutorial's own QC **both arms are deleted outright**. With the reachability-normalised QC they *still* land differently (F19 vs F1) and **neither is a top-10 feature of any factor** |
| LRIC / cross-PCF | GRN→SORT1 is genuinely directional (7/8 punches, p = 0.0156) but the **LRIC/cross-PCF ratio is ~1 in every bin**, so co-occurrence is fully explained by cell-type co-location |

#### ⚠️ The dominant reason is the UNIT OF ANALYSIS — feature indexing is secondary

The four rows above describe *symptoms* of the factorisations; the cause sits one level upstream, in
what a **row** of LIANA's matrix is (`why_no_mgam_motif.py` → `vs_alarmist/why_no_mgam_motif.json`).

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

**Pearson rises 26-fold purely by aggregating cells into patches. The loop is a property of a
*neighbourhood*, and a neighbourhood is not a row in LIANA's matrix** — so no factorisation applied
to that matrix can recover structure the input does not contain. Feature indexing (sender inside the
feature name, receiver nowhere) is real and compounds this, but it is the **secondary** effect:
correcting the QC to restore every sender view changed neither the co-loading nor the grade result.

Two further points that only fall out of the whole set: **LRIC is the only branch that resolves
direction**, and it is the only one that supports the forward arm as directional at the correct
replicate unit; and **LRIC also refutes the "wider kernel" explanation**, its ratio staying at ~1
out to 225 µm.

The one dissent is **LR-MISTy**, which finds GRN a top-10 predictor of both SORT1 and FPR1 — but
that branch is a niche-composition predictor rather than an LRI detector, so it does not overturn
the above. **CellChat recovers both arms** (see its section), so **the disagreement is between
comparators, not between ALARMIST and the field.**

#### Factor-vs-motif cosine — they agree on vocabulary, not on programs

`cosine_factors_vs_motifs.py` → `vs_alarmist/`, following the matching procedure in
`.claude/skills/alarmist`. Compares the **reachability-normalised** MOFA-Flex loadings against
`results/GBM/bptf/lri_motifs.csv`. Three confounds had to be handled, and each costs something:

1. **The feature spaces are not the same object.** ALARMIST is (sender, receiver, ligand, receptor,
   contact mode) — **25,271** features per motif; MOFA-Flex is (sender, ligand, receptor) — **779**.
   The only common space is the latter, so **ALARMIST must be collapsed by SUMMING over receiver and
   contact mode**: 25,271 rows → **4,756** keys, median **5** rows merged per key. *Worked example —
   `mGAM|GRN|SORT1` on motif 1:* the MES-like receiver arm scores **3.091** of a summed **12.369**,
   i.e. **25%**. The biologically meaningful direction is a quarter of the number that enters the
   cosine, and mGAM→mGAM autocrine is summed in as if equivalent. **The comparison is therefore
   deliberately biased in LIANA's favour — an UPPER BOUND on agreement, not a neutral measurement.**
   ALARMIST mass retained after collapse + join: 90.9% on raw `V`, 74.0% on `V*`.
2. **DB export mismatch** (the subunit-ordering trap above). Canonicalising subunits takes the
   overlap from **713** to **742** keys.
3. **⚠️ Signed vs non-negative — this one silently manufactured a result.** MOFA-Flex weights are
   signed; BPTF factors are not. Three modes are implemented, `--sign-mode {poles,abs,signed}`:

| mode | what it does | verdict |
|---|---|---|
| `abs` | `abs(weight)` | **WRONG, and it is what was originally used.** **57.2%** of weights are negative and the minor pole holds a median **38.5%** of a factor's mass, so `abs` **merges two anti-correlated poles** and manufactures similarity. Kept only for provenance |
| `poles` | split each factor into `max(w,0)` and `max(−w,0)` → 40 non-negative vectors | **default.** Each pole is a genuine non-negative feature set, directly comparable to a motif |
| `signed` | signed weights as-is; cosine may be negative | plotted on a diverging scale centred at 0 |

Results, ALARMIST scored on **`V*` = V/(mean_LR + 1)** (the prevalence-normalised column the skill
says to rank on), 742 shared keys, 200-permutation null preserving sparsity and magnitude:

| sign mode | max cosine | median best-match | motifs with a match > 0.5 | motifs beating the null at p < 0.05 | null median max-cosine |
|---|---|---|---|---|---|
| `abs` | 0.671 | 0.432 | **3 / 20** | 18 / 20 | 0.241 |
| **`poles`** | **0.743** | **0.517** | **13 / 20** | **20 / 20** | 0.221 |
| `signed` | +0.478 (max \|cos\| 0.643) | 0.345 | **0 / 20** | 20 / 20 | 0.088 |

**Reading:** the two methods **agree substantially on LR VOCABULARY** (13/20 above 0.5, 20/20 above
the null under `poles`) but **not on how that vocabulary groups into programs** (the cell-space
Spearman ceiling below is |ρ| ≈ 0.46, with 9 of 20 motifs collapsing onto one hub factor) and **not
on receiver or direction**, which are absent from LIANA by construction.

Three guards on that conclusion. **It is not a multiple-comparisons artefact** — going from 20 to 40
candidate vectors *lowered* the permutation null (0.241 → 0.221), because pole vectors are sparser.
**`aggfunc` is not load-bearing** — summing / max / mean when collapsing ALARMIST moves the third
decimal, not the conclusion (0.743 / 0.783 / 0.760 with 13 / 15 / 13 motifs above 0.5 under `poles`).
**On raw `V` instead of `V*`, 19/20 motifs clear 0.5 under `abs`** — that is a shared-prevalence
artefact, both methods being dominated by the same ubiquitous adhesion pairs, **not agreement**.

#### Per-cell Spearman — one hub factor, and it is an activity axis

`compare_programs_to_alarmist.py` → `vs_alarmist/comparison_summary.json` + `data/rho_*.csv`. Both
methods emit a per-cell × per-program matrix over the **same 100,190 cells**, so Spearman is
directly computable — no collapsing, no key repair.

| LIANA program set | programs | max \|ρ\| | motifs whose best match is the SAME factor | hub vs total ALARMIST loading | motif 1's best match |
|---|---|---|---|---|---|
| MOFA-Flex, tutorial QC | 17 | **0.518** | **13 / 20** → Factor 18 | **+0.506** (+0.476 vs total inflow) | Factor 11, ρ **−0.222** |
| MOFA-Flex, sensitivity | 20 | 0.458 | 6 / 20 → Factor 18 | +0.400 | Factor 7, ρ −0.266 |
| MOFA-Flex, reachability-normalised | 19 | 0.456 | **9 / 20** → Factor 18 | **+0.422** (+0.362 vs total inflow) | Factor 18, ρ **+0.214** |

**The hub is a general-activity axis, not a program.** It correlates +0.42 to +0.51 with a cell's
*total* ALARMIST loading and +0.30 to +0.48 with its *total* inflow — so most of what looks like
motif↔factor agreement in cell space is "this cell has a lot of signalling in it". By contrast
**ALARMIST motif 1 vs total inflow is only +0.146**, i.e. **motif 1 is specific, not an activity
axis, which is precisely why it has no good match.**

**Cell-type placement disagrees at the top.** ALARMIST motif 1 is highest in **mGAM (0.794)** and
lowest in Glial-Neuronal (0.307). Its best-matching factor (reach-norm Factor 18) peaks instead in
**NPC-like (0.654)** and puts mGAM **fifth** (0.504); the two agree only at the bottom.

⚠️ **The NMF row of that script is UNVERIFIED — do not quote it.** For `kind == "csv_nmf"` the script
re-reads `nmf_inflow/data/NMF_W_factor_scores.csv` **without `index_col`**, so the numeric `cell_id`
column survives `select_dtypes("number")` and is scored as if it were a factor: the summary reports
**8 factors for a rank-7 NMF**, and `data/rho_nmf_8f.csv` has 8 columns. **Open issue; the three
MOFA-Flex rows are unaffected.**

#### The factor count must never be set beside ALARMIST's K

ALARMIST used **K = 20** BPTF motifs. Inflow+NMF gives **7**, bivariate+NMF **6**, the `default`
tier **6 / 5**, MOFA-Flex **17 / 19 / 20**. The spread is real (5, 6, 7, 17, 19 or 20 on the same
tissue) but it is driven by the **factorisation, the elbow window and the bandwidth**, not by the
resource — and the MOFA-Flex numbers are not outputs at all, since `n_factors` is a binding ceiling.
**The factor count is not a property of the tissue.** The *structural* comparison stands; the
numerical one does not. And per *Factor annotation*, six of the seven inflow NMF factors are ≥ 75%
one sender, so even structurally they are closer to cell-type identity than to ALARMIST motifs.

### ⚠️ Contract deviations — OPEN, awaiting sign-off

Departures from `.claude/skills/comparator-benchmark/SKILL.md`, **not** from the LIANA tutorials.

| # | Invariant | What SKILL.md says | What was done | Status |
|---|---|---|---|---|
| **CD-1** | kernel scale | `SKILL.md:45-46`, verbatim: *"Keep each method's own neighborhood/kernel definition at its default. Do NOT harmonize spatial scale across methods, and do not match it to ALARMIST's patch size."* `SKILL.md:105` also lists an unpinned kernel scale under **STOP-and-ask** | bandwidth **13.1454 µm** derived from ALARMIST's 50 µm patch edge. Five `run_manifest.json` carry it (`cellchatdb2`, `cellchatdb2_inflow`, `cellchatdb2_morans`, `default`, `default_inflow`); the four `nmf_*` ones carry no `bandwidth` key. **Scope: CD-1 does NOT apply to LRIC** (own `cKDTree`) **or LR-MISTy** (tutorial `bandwidth=200`); it **DOES** apply to MOFA-Flex, which factorises the inflow matrix built at that bandwidth | ❌ **OPEN — violates the invariant.** Not a defensible "default", because LIANA pins none. **Load-bearing:** the sensitivity table measures this one parameter moving the inflow NMF rank **11 → 7**. Kept by user decision after the 2026-08-06 evidence review — **not resolved by it** |
| **CD-2** | native multi-sample / differential mode | `SKILL.md:47-49`: GBM → split by `obs['grade']`, with the **13 `obs['tma_id']` cores as the units** | `li.mt.compute_global_specificity(groupby='grade')` only. Reading the installed source, that is a one-sided **per-group specificity** test permuting labels across **cells** — not a contrast. `region_global_interactions.csv` is 9,216 rows = 4,608 × {high, low} with **no contrast column**, and `grep -n tma_id run_inflow_downstream.py` returns **nothing** — the cores never enter the grade analysis | ❌ **OPEN — requirement NOT satisfied, and inapplicable.** LIANA has no native spatial differential mode, and `SKILL.md:49` says *"If the method has no multi-sample mode, say so — do not hand-roll one."* The 5,417 rows at p<0.05 are **cell-level p-values (n = 100,190)**, pseudoreplicated by ~4 orders of magnitude relative to the 13 cores |

**Wherever a grade-associated p-value from this method is quoted, carry the CD-2 caveat**: it is a
cell-level p-value, not a per-core one, and its magnitude is not interpretable. **Treat the
*direction* of an effect as the result, not the *p*.**

**A punch-level test was run as an *additional* analysis, and it is null everywhere.**
`analyse_existing.py` → `nmf_inflow/punch_level/` aggregates the inflow `NMF_W` by `obs['tma_id']`
and runs a two-sided Mann-Whitney over **7 high vs 6 low** punches, BH-corrected:

| test | signif. at BH q<0.05 | smallest raw p | smallest q |
|---|---|---|---|
| 7 inflow NMF factors vs grade | **0 / 7** | 0.013986 (Factor4, log2FC +1.51) | 0.0816 |
| 20 required-LR features vs grade | **0 / 20** | 0.013986 (`MES-like^GRN^SORT1`, log2FC +2.62) | 0.1399 |
| MOFA-Flex primary (17 active factors) | **0 / 17** | 0.013986 (Factor 18) | 0.2179 |
| MOFA-Flex sensitivity (20 factors) | **0 / 20** | 0.008159 (Factor 11) | 0.1632 |
| MOFA-Flex reachability-normalised (19 active) | **0 / 19** | 0.022145 (Factor 18) | 0.3322 |
| LRIC / cross-PCF, cell-type-resolved | ⚠️ **not testable** — 8 of 13 punches informative and **7 of those 8 are high grade**, the loss systematic. Agnostic fallback null (p = 0.234 / 0.788) | | |

**Five punch-level tests, five nulls (one untestable)** — consistent, but **not** evidence of no
grade effect: every one is bounded by the 7-vs-6 floor of p = 0.0011655. This test is
**hand-rolled**, which `SKILL.md:49` discourages; it is recorded as an additional analysis and does
**not** convert CD-2 to satisfied.

### Gotchas

- **`query_bandwidth` returns a column spelled `bandwith`** (missing the second `d`). Indexing
  `'bandwidth'` raises `KeyError`.
- **`max_neighbours=100`** is a KNN ceiling — the same class of cap that binds destructively in
  SpatialDM. Non-binding for bivariate/inflow (max 52), **binding for 99.7% of cells in LR-MISTy**.
- **`nz_prop=0.2` in the tutorial is 20% of Visium *spots*.** On single-cell data the median gene is
  detected in **4.3%** of cells, so 0.2 is above the 90th percentile of genes and keeps only **13**
  pairs. See `DEVIATIONS.md` for the binomial spot→cell conversion.
- **Gene names containing `_` collide with the complex-subunit separator.** LIANA only *warns*
  (stLearn hard-errors); the 21 `Intergenic_Region_*` control probes are dropped anyway.
- **With `n_perms=100` the p-value grid is coarse** — the floor is 0, and **106/131** pairs report
  `morans_pvals = 0`. p-values here separate "spatially structured" from "not"; **they do not rank.**
- Moran's R rewards spatial *structure*, not abundance. A ubiquitous, uniformly expressed pair scores
  near zero **by design**.
- **A replot must not clobber fit provenance.** A 2026-08-04 regeneration pass overwrote both
  MOFA-Flex manifests, replacing `fit_seconds` / `determinism_probe` with the *replot's* cost. The
  keys were restored from `logs/mofaflex_{primary,sensitivity}.log`, the only surviving source, and
  each manifest carries a `provenance_note` saying so. `run_mofaflex.py` now carries fit provenance
  forward when reusing a cached model and records the replot separately under `last_replot`.

### Methods paragraph

> For LIANA+ (v1.8.1), we identified spatial neighbours using the `spatial_neighbors` function and
> computed bivariate scores using the `bivariate` function, following the authors' tutorial. The
> Gaussian kernel bandwidth was set to 13.1454 µm with a cutoff of 0.1, corresponding to a support
> radius of 28.2 µm; this value was fixed by an area-preserving correspondence between the kernel's
> support disk and a 50 µm square patch, because the inflow branch of LIANA+ that applies to
> single-cell data specifies no numeric bandwidth rule. The resulting graph gave a median of 14
> neighbours per cell (maximum 52, below the `max_neighbours` ceiling of 100). Local scores were
> computed with the weighted cosine metric and global scores with bivariate Moran's R, using 100
> permutations and a minimum non-zero expression proportion (`nz_prop`) of 0.02; the latter replaces
> the tutorial's 0.2, which refers to multi-cell Visium spots. Ligand–receptor pairs were taken from
> CellChatDB v2, whose heteromeric complexes LIANA+ represents natively; 131 of 3,218 pairs passed
> the expression filter. All thirteen tissue microarray cores were analysed together, after
> verifying that the minimum distance between cells of different cores (222.9 µm) exceeds the
> kernel's support radius by 7.9-fold. Because LIANA+ provides no native spatial differential mode,
> no differential test between tumour grades was performed on the bivariate scores. We additionally
> ran the single-cell branch of the authors' decision tree (`inflow`), decomposed it both with
> non-negative matrix factorisation (`li.multi.nmf`) and with the authors' prescribed MOFA-Flex
> route, and ran the spatial co-occurrence (`lric`, `cross_pcf`) and multi-view (`lrMistyData`,
> `misty`) branches; the latter two were run at their own default spatial supports.

---

## NICHES — R, v1.2.4, env `comp-niches`

Local source `/Users/jiayifan/tansey_lab/NICHES` @ `d698e37b` (2026-01-29). R 4.3.3, Seurat 5.3.0.
Contract: `niches/NOTES.md`. Deviations: `niches/DEVIATIONS.md`. Code: `prepare_gbm_input.py`,
`run_niches.R` (per core), `run_all_cores.sh` (the driver over 13 cores × 2 imputation sub-runs),
`niches_io.R` (persistence), `analyze_niches.R` (merge + embed + differential),
`summarize_niches.py`, `audit_edges.py` (the neighbourhood audit), plus the vendored
`vendor/seurat-wrappers/` ALRA source.

### Core algorithm

**NICHES is a *transformation*, not a test.** For a directed edge (sending cell *i* → receiving
cell *j*) and mechanism *L→R* it computes the **product of normalized ligand expression in *i* and
normalized receptor expression in *j***, with multi-subunit complexes multiplied across subunits.
`NeighborhoodToCell` averages that product over every edge landing on *j* (`blend = "mean"`),
giving one **niche vector per cell**; `CellToCellSpatial` keeps one vector **per edge**. The result
is handed back as an ordinary Seurat assay (mechanisms as "genes", cells or edges as "cells").

**There is no null model, no permutation and no p-value anywhere in the scoring step.** Every
inferential statement comes from whatever Seurat test you then run — here the vignettes'
`FindAllMarkers(test.use = "roc")`. So a "significant interaction" in NICHES is a *downstream*
claim about a group contrast, and the unit it attaches to is the receiving cell or the cell–cell
edge, never the interaction itself.

### Spatial model

**Mutual k-nearest-neighbour graph, `k = 4` (package default), no kernel and no distance cutoff.**
`ComputeEdgelist.R:46-56` ranks every cell's neighbours by euclidean distance, keeps the `k+1`
nearest, and symmetrises with `adj & t(adj)` — an edge survives only if each cell is in the other's
top 4. `rad.set` offers a hard radius instead but is ignored whenever `k` is non-NULL, and
`RunNICHES` always passes `k`. Three consequences, all measured on our data
(`summary_edge_audit.csv`, `audit_edges.py`):

- **`k` is unitless, so NICHES has no radius to quote** — the neighbourhood must be measured after
  the fact. Pooled over all 13 cores the real (non-self) edges have **median 10.1 µm, p95 26.6 µm,
  p99 37.8 µm, max 243.1 µm**. That is **the tightest neighbourhood of any method here** (LIANA+
  28.2, SpatialDM 135, CytoSignal 200, stLearn 250, COMMOT 365) — NICHES is effectively scoring
  *directly abutting* cells.
- **The radius is density-dependent, and density tracks grade.** Because `k` is fixed, a denser core
  gets a physically *smaller* neighbourhood: per-core median edge length runs from **8.1 µm (core 1,
  26,456 cells) to 22.2 µm (core 6, 3,092 cells)**, a 2.7× spread. **The high- and low-grade arms
  are not being measured at the same physical scale.**
- **Every cell is its own neighbour.** `order(dis_vec)[1:(k+1)]` includes the cell itself at distance
  0 and the mutual-NN step always keeps it, so the edge list contains one self-edge per cell:
  **100,197 of 419,299 edges (23.9%)** are `cell—itself`. `NeighborhoodToCell` therefore mixes each
  cell's own **autocrine** ligand×receptor product into its "neighbourhood" average. **Package
  behaviour, not a configuration choice, and not mentioned in the vignettes.**

### LR database

Default is **FANTOM5** (`ncomms8866_human`, bundled `.rda`); OmniPath is the other built-in. For
`cellchatdb2` we pass `LR.database = "custom"` with `custom_LR_database` = our CellChatDB v2 CSV.
**The conversion is trivial** — `LoadCustom` wants exactly a 2-column data.frame whose first column
is the ligand subunits and second the receptor subunits, each `_`-separated, which is already our
export's format, so the "conversion" is a straight `db[["ligand","receptor"]]` selection with no
complex re-encoding. Complexes are handled by **multiplying subunit expression**, and
`FilterGroundTruth` keeps a mechanism only if **every** subunit is present in the object. On the
5,119-gene panel that takes CellChatDB v2 from **3,218 unique pairs → 1,088 mechanisms**, identical
in every core. `species` is silently ignored in custom mode (`LoadCustom.R:12`). The FANTOM5
`default` tier was **not run**.

### Input

A Seurat object with (a) **raw counts**, normalized by `NormalizeData` inside the runner (the
vignettes normalize themselves, so counts must be handed in un-normalized — we read
`layers['counts']`, since this h5ad's `X` is already log-normalized); (b) **`x` and `y` as ordinary
`meta.data` columns**, not an `@images` object — units are irrelevant because `k` is unitless, but
ours are µm; (c) a **cell-type column** (`cell_type`, 9 types), which `RunNICHES.Seurat` copies into
`Idents`. No QC filter is applied: `min.cells.per.ident` and `min.cells.per.gene` are left at their
`NULL` defaults, so all 100,197 cells and all 5,119 genes enter.

### Workflow

| Step | Call with argument values | Produces |
|---|---|---|
| 1 | `CreateSeuratObject(counts, meta.data)` — per TMA core | 5,119 × n_cells object |
| 2 | `NormalizeData(obj)` — `LogNormalize`, `scale.factor = 1e4` | `RNA` `data` layer |
| 3 | *(alra sub-run only)* `RunALRA(obj)` — `k = NULL` (auto), `q = 10`, `use.mkl = FALSE` | `alra` assay |
| 4 | `RunNICHES(assay, LR.database = "custom", custom_LR_database, species = "human", cell_types = "cell_type", position.x = "x", position.y = "y", k = 4, rad.set = NULL, blend = "mean", min.cells.per.ident = NULL, min.cells.per.gene = NULL, CellToCellSpatial = TRUE, NeighborhoodToCell = TRUE, CellToCell = FALSE, output_format = "seurat")` | 2 Seurat objects per core |
| 5 | tag `Condition <- grade`, `Core <- tma_id`; `merge()`; `JoinLayers()` | one 13-core object |
| 6 | *(CellToCellSpatial only)* `subset(nFeature_CellToCellSpatial > 5)` | vignette 04's low-information filter |
| 7 | `ScaleData` → `FindVariableFeatures(selection.method = "disp")` → `RunPCA(npcs = 100)` → `RunUMAP(dims = 1:50)` | embedding |
| 8 | `FindAllMarkers(min.pct = 0.25, only.pos = TRUE, test.use = "roc")` with `Idents` = ReceivingType, then grade | marker tables |
| 9 | per ReceivingType: `subset` → `ScaleData` → `FindVariableFeatures` → `RunPCA(npcs = 50)` → `RunUMAP(dims = 1:40)` → `FindAllMarkers` → `DoHeatmap(top_n(20, myAUC))` | per-population grade contrast |

### Data outputs

Per core, per imputation sub-run:

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
`summary_requested_lr.csv`, `summary_edge_audit.csv`, `run_timings.csv`. Analysis-level
(`_analysis/<org>/`): `merged.rds`, `mechanism_detection.csv`, `markers/markers_*.csv`,
`differential_summary.csv`, `composition_*.csv`, `requested_*_by_ReceivingType_grade.csv`,
`analysis_manifest.json`.

### Image outputs

All as `.png` + `.pdf` (+ `.svg` under 50k points; skipped above that and logged, since an SVG of a
100k-point UMAP is hundreds of MB).

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

Not produced: `SpatialFeaturePlot` (needs a Seurat `@images` object the Xenium h5ad has none of —
replaced by an equivalent ggplot scatter) and `ALRAChooseKPlot` (diagnostic for ALRA's rank choice;
the chosen `k` is recorded per core instead).

### Multi-sample / differential mode

**Native, and it is the split-run-merge pattern**, not a joint model. Vignettes 04 and 07 both say
to *split the object by condition first*, run NICHES separately on each, tag the outputs, `merge()`,
then do ordinary Seurat differential testing on the merged matrix. We apply it at the TMA-core level
(13 runs) and contrast `grade` on the merged object. **There is no shared latent space, no batch
term and no random effect for core — every cell is an independent observation in the marker test,
which is pseudoreplication with respect to the 13 cores. That is the method's own design, not our
choice.**

### Gotchas

1. **`ComputeEdgelist` is dense O(N²) and there is no fast path.** `apply(df, 1, ...)` materialises
   the full N × N double distance matrix, then `adj_mat`, `t(adj_mat)` and `1*(...)` — roughly
   4 × N² × 8 bytes. At 100,197 cells that is **~80 GB for the distance matrix alone**. The
   `nn.method = 'aoz'` spNNGP fast path exists in the signature but its body is **commented out**
   (`ComputeEdgelist.R:76-131`) and `RunNICHES` never passes the argument. **This is why the slide
   is split per core.**
2. **Every cell is its own spatial neighbour** — 23.9% of edges are self-edges, so
   `NeighborhoodToCell` blends autocrine signal into the niche. Undocumented in the vignettes.
3. **`k` is unitless, so the physical neighbourhood shrinks as density rises** — 8.1 µm in the
   densest core vs 22.2 µm in the sparsest.
4. **NICHES' `data` slot is not log-space.** It holds raw LR products in both `counts` and `data`
   (`RunNeighborhoodToCell.R:77-81`). Any downstream Seurat call that assumes log1p — **including
   the vignette's own `FindVariableFeatures(selection.method = "disp")`** — is being fed the wrong
   space. On the imputed tier the products reach 4,933.5 and `exp()` overflows, crashing
   `CalcDispersion` outright. **Never `NormalizeData` a NICHES assay.**
5. **Em-dash, not hyphen.** Mechanisms are `LIGAND—RECEPTOR` with U+2014. Every `FetchData` /
   `FeaturePlot` lookup must use it.
6. **Seurat rewrites underscores in feature names** — `CreateSeuratObject` warns and replaces `_`
   with `-`, so a multi-subunit mechanism is stored as `TGFB1—TGFBR1-TGFBR2`. Both requested LRs are
   single-subunit and unaffected.
7. **`species` is silently ignored** when `LR.database = "custom"`.
8. **`BisRNA` is declared in `Imports` but never called** anywhere in `R/`, and it is archived on
   CRAN — install from the archive purely to satisfy `R CMD INSTALL`.
9. **`SeuratWrappers` 0.4.0 will not install** because of its `Banksy` dependency chain; `RunALRA`
   itself needs none of it.
10. **ALRA picks its rank per object.** Across our 13 cores `RunALRA` chose **k = 20–56**, and higher
    rank means less smoothing: core 6 (k = 56) ends up with **611 / 1,088** mechanisms detected while
    every comparable core reaches ~1,000. **Since cores are fit independently, imputation strength is
    not comparable across cores** — a confound layered on top of the density confound.
11. **No memory cliff, but a wall-clock one.** Peak RSS was well below the naive estimate (14.5 GB
    noimpute / 19.8 GB alra on the 26,456-cell core, vs ~26 GB projected) because R reclaims the
    distance matrix before the scoring stage.

### Deviations from the vignettes

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| Unit of a run | one section | one **TMA core**, 13 runs, merged | dense O(N²) edgelist = ~80 GB whole-slide; cores are physically disjoint anyway |
| Imputation | vignettes 01/07 use ALRA; 04 does not | **both**, as parallel sub-runs | requested; and it is decisive — mean mechanisms detected per core **526 vs 966** of 1,088, and the whole grade differential is **null without imputation** (6 vs 229 mechanisms reach `min.pct = 0.25`) |
| Source of `RunALRA` | `library(SeuratWrappers)` | authors' **unmodified** `alra.R` + `internal.R` @ `8df8343`, vendored | SeuratWrappers 0.4.0 won't install (Banksy's Bioc chain fails, 16 pkgs) |
| LR database | `"fantom5"` | `"custom"` = CellChatDB v2.0 human | the `cellchatdb2` tier; format already matches `LoadCustom` |
| `meta.data.to.map` | vignette 07 omits (all cols) | names 6 columns | guarantees x/y/grade carry through to the later stages |
| ReceivingType loop | one hand-picked population | **all 9** | no a-priori population; hand-picking would be a silent choice |
| `SpatialFeaturePlot` | Visium `@images` | equivalent `ggplot` scatter | Xenium h5ad has bare coordinates, no image object to dispatch on |
| `nFeature > 5` filter | vignette 04, on CellToCell | CellToCellSpatial only | vignette 07 applies none to NeighborhoodToCell |
| `selection.method = "disp"` | vignettes 01/04/07 | `"disp"` on noimpute; **`"vst"` on alra** | the vignette's own call **crashes** on imputed data (gotcha 4) |
| `JoinLayers` after `merge` | n/a (pre-v5 vignettes) | added | Seurat v5 keeps one layer per merged object |

### Runs on our data

| Tier | Sub-run | Status | Key numbers |
|---|---|---|---|
| `cellchatdb2` | `noimpute` | ✅ 13/13 cores, 3.9 min | 100,197 cells → **419,299 edges** (319,102 real + 100,197 self); **1,088** mechanisms; **526 detected/core on average** (112–702); mean niche density **0.73%**; peak RSS 14.5 GB |
| `cellchatdb2` | `alra` | ✅ 13/13 cores, 9.2 min | same 419,299 edges; **966 detected/core on average** (529–1,088); mean niche density **9.81%**; ALRA rank k = 20–56 per core; peak RSS 19.8 GB |
| `cellchatdb2` | `noimpute` differential | ✅ **grade contrast null** | only **6 / 1,088** mechanisms reach `min.pct = 0.25`; **zero grade markers** globally and in all 9 niches. The cell-type contrast is not quite empty: **1** marker survives (`CNTN2—CNTN2` for Glial-Neuronal, `myAUC` 0.795) |
| `cellchatdb2` | `alra` differential | ✅ | **229 / 1,088** mechanisms reach `min.pct = 0.25`; **180** cell-type markers, **28** global grade markers, 3–98 per-niche grade markers (Glial-Neuronal: 0) |
| `cellchatdb2` | `noimpute` CellToCellSpatial | ✅ | the `nFeature > 5` filter keeps **89,648 / 419,299** edges (21.4%); **13 / 1,088** mechanisms reach `min.pct = 0.25`; **53** VectorType markers, **0** grade markers |
| `cellchatdb2` | `alra` CellToCellSpatial | ✅ | filter keeps **381,690 / 419,299** edges (91.0%); **116 / 1,088** mechanisms; **215** VectorType markers, **1** global grade marker; `disp` overflowed again (max product **9,693.1**) → `vst` |
| `default` (FANTOM5) | — | ❌ not run | this pass is `cellchatdb2`-only |
| LUAD | — | ❌ not run | not attempted |

⚠️ **Two confounds specific to running NICHES on this TMA**, both quantified in
`summary_confound.csv` / `summary_edge_audit.csv`, and **both push in the same direction as grade**,
so a naive reading of the high-vs-low contrast will overstate it: **detection density tracks core
cellularity** (r = 0.881 alra / 0.410 noimpute; for essentially *every* mechanism `frac_high` is
about twice `frac_low` — NCAM1—NCAM1 55.3% vs 39.0%, DLL3—NOTCH1 32.5% vs 14.7%), and **the
neighbourhood is physically smaller in the denser high-grade cores**. The **cell-type** contrasts
are far less exposed, since they compare populations *within* the same cores.

### Requested LRIs

Both are in CellChatDB v2.0 human and all four genes are on the Xenium 5K panel, so both survive
`FilterGroundTruth` and are scored **in all 13 cores in both sub-runs** — nothing is missing.

| LR | Sub-run | Detection | Rank by detection | Verdict |
|---|---|---|---|---|
| **GRN—SORT1** | `alra` | **69.2%** of cells (mean 58.1% per core) | mean rank **22** of 1,088; top-10 in 6/13 cores | **strongly recovered** |
| **GRN—SORT1** | `noimpute` | 8.5% of cells | mean rank 28 | recovered but sparse |
| **ANXA1—FPR1** | `alra` | 12.1% of cells | mean rank 307 | present, low-ranked overall — **but see below** |
| **ANXA1—FPR1** | `noimpute` | 1.2% of cells | mean rank 127 | barely detectable |

Two results worth stating plainly, both from the `alra` sub-run:

- **`ANXA1—FPR1` is the single best marker of the mGAM niche** in the vignette-01 cell-type contrast
  — `myAUC = 0.866`, `avg_log2FC = +7.60` for `ReceivingType == mGAM`. **Its globally low rank is a
  *prevalence* statement** (FPR1 is a rare myeloid receptor); **its specificity to mGAM is the
  strongest of any mechanism NICHES scored.** An independent corroboration of the ANXA1→FPR1 arm.
- **`GRN—SORT1` is a high-grade marker in 5 of the 9 niches** — AC-like (0.787), MES-like (0.806),
  NPC-like (0.763), OPC-like (0.747), Vascular (0.708), all up in high grade. It is also a cell-type
  marker of the Vascular niche (0.739). Read against the density confound, **the direction is
  consistent but the effect size is inflated by core cellularity.**

One cross-method concordance: **`JAM3—F11R` is the top low-grade mGAM marker** here
(`myAUC = 0.291`, i.e. down in high grade), matching CytoSignal's `JAM3–F11R` (log-FC −1.10), one of
only four FDR-significant interactions in its grade test.

### Methods paragraph

> Cell–cell signalling was inferred with NICHES v1.2.4 (Raredon et al.) in R 4.3.3 / Seurat 5.3.0.
> Each of the 13 TMA cores was processed independently: raw counts were normalized with
> `Seurat::NormalizeData` (`LogNormalize`, scale factor 1e4) and, in a parallel sub-run, imputed
> with `RunALRA` at its default automatic rank. `NICHES::RunNICHES` was then applied per core with
> `LR.database = "custom"` supplying CellChatDB v2.0 human (3,218 ligand–receptor pairs, of which
> 1,088 had every subunit represented on the 5,119-gene Xenium panel), `cell_types = "cell_type"`,
> spatial coordinates in microns, and the package-default mutual k-nearest-neighbour graph (`k = 4`,
> `blend = "mean"`), requesting the `CellToCellSpatial` and `NeighborhoodToCell` organizations.
> Per-core outputs were tagged with tumour grade, merged, and embedded with `ScaleData`,
> `FindVariableFeatures`, `RunPCA` (100 PCs) and `RunUMAP` (50 PCs). Differential signalling was
> tested with `Seurat::FindAllMarkers` (`test.use = "roc"`, `min.pct = 0.25`, `only.pos = TRUE`)
> both across the nine receiving cell types and, within each, between high- and low-grade cores.

---

## CellChat — R, v2.2.0.9001, env `comp-cellchat`

Local source `/Users/jiayifan/tansey_lab/CellChat` @ `75253cd0` (2026-03-04). R 4.3.3, Seurat 5.1.0,
presto 1.0.0, NMF 0.27. Contract: `cellchat/NOTES.md` (which also carries this method's deviations —
it has no separate `DEVIATIONS.md`). Code: `prepare_gbm_input.py`, `audit_db_equivalence.R`,
`run_cellchat.R` (inference), `cellchat_io.R` (persistence + the `save_all_formats()` R saver),
`plot_cellchat.R` (the full plot suite), `rebuild_chords.R` / `rebuild_river_plots.R` (targeted
replots), `control_nonspatial.R`, `finalize_manifest.R`, `install_env.R`, `run_all_gbm.sh` (the
driver), `build_report.py`.

Our data is spatial, single-cell resolution, multi-section, two conditions, and **no single vignette
covers that**, so the contract is assembled from five in a stated precedence order:
`..._multiple_spatial_transcriptomics_datasets.Rmd` (**V2, primary** — the one that matches our
shape: many sections in one object via `meta$samples`), `..._spatial_transcriptomics_data.Rmd`
(**V1**), `FAQ_on_applying_CellChat_to_spatial...Rmd` (**VF** — the Xenium `spatial.factors`),
`Comparison_analysis_of_multiple_datasets.Rmd` (**VC** — high vs low grade) and
`CellChat-vignette.Rmd` (**VB** — the downstream plot inventory both spatial vignettes defer to).

### Core algorithm

For each L-R pair and each ordered pair of **cell groups** (i → j), CellChat computes a Hill
function of the product of group-average ligand expression in i and group-average receptor
expression in j: `P = L·R / (Kh^n + L·R)` with `Kh = 0.5`, `n = 1`. L and R are computed over
`data/max(data)` with a **10% truncated mean** per group; multi-subunit complexes enter as the
**geometric mean** across subunits; the receptor term is further multiplied by co-activation and
divided by co-inhibition receptor terms. Significance is a **permutation test**: cell-group labels
are shuffled `nboot = 100` times and `pval = #{P_boot >= P_obs}/nboot`.

**The unit of a CellChat result is a (sender cell type, receiver cell type, L-R pair) triple** — not
a cell, not a spot, not an edge. `netP` sums the probabilities of all L-R pairs in a pathway;
`aggregateNet` counts significant links and sums probabilities per cell-type pair.

### Spatial model

Space enters at the **cell-group** level, not the cell level. For each sample k and ordered group
pair (i,j), `computeRegionDistance` takes every cell of group i, finds its **1-nearest neighbour**
in group j (`BiocNeighbors::queryKNN`, Annoy), converts to µm by `× ratio[k]`, and takes the **10%
trimmed mean** → `d.spatial[i,j,k]`. Then `adj.spatial[i,j,k] = 1` iff ≥ `k.min = 10` distinct
group-j cells lie within `interaction.range + tol` (250 + 5 µm); `adj.contact[i,j,k] = 1` iff ≥ 10
lie within `contact.range + tol` (10 + 5 µm); across samples `d.spatial` is **averaged** and the
adjacency matrices are 1 if **any** sample says 1, then symmetrised. **Group pairs with
`adj.spatial == 0` become `NaN` and are excluded entirely.** With `distance.use = FALSE` (V2's
value, ours) distance is a **hard filter only**; with `TRUE` (V1) it additionally down-weights by
`1/(d × scale.distance)`.

**Consequence for this benchmark: CellChat produces no per-cell and no per-spot output at all.** Its
finest spatial statement is "these two cell types have ≥10 mutual neighbours within 250 µm somewhere
in this section" — **coarser than every other method here**, and the reason its counts are not
comparable to CytoSignal's per-cell or stLearn's per-spot tallies.

### LR database — and why the two tiers cannot differ by resource

**The bundled `CellChatDB.human` *is* CellChatDB v2** (3,233 rows, 338 complexes, 32 cofactors).
`audit_db_equivalence.R` re-derives the repo CSV's flattening from it using the mapping CLAUDE.md
documents and diffs the key sets:

| | bundled `CellChatDB.human` | `data/LRdatabase/CellChatDBv2.0.human.csv` |
|---|---|---|
| rows | 3,233 | 3,233 |
| unique ligand\|receptor keys | 3,218 | 3,218 |
| Secreted / Non-protein / Cell-Cell Contact / ECM-Receptor | 1280 / 994 / 535 / 424 | 1280 / 994 / 535 / 424 |

**Shared keys 3,218 — only-bundled 0, only-repo 0, Jaccard 1.0000.** The 3,233 vs 3,218 gap is 15
keys listed twice under two annotations (the three `POMC|OPR*` pairs appear as both Secreted and
Non-protein Signaling, and the merge cross-products them), which is the entire content of the "6
signaling_type disagreements" the audit reports. Both sides carry both rows.

So the tiers differ by **annotation scope**, the only knob the vignettes expose:

- **`default`** — `subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")`, the
  literal call in V1:127 and V2:131 → **1,280 interactions, 158 pathways**.
- **`cellchatdb2`** — `subsetDB(CellChatDB)`, the alternative commented in on the next line of the
  same vignette → **2,239 interactions, 252 pathways**, i.e. Secreted + ECM-Receptor + Cell-Cell
  Contact. **This is the scope closest to the resource ALARMIST uses.**

Re-importing the flat repo CSV through `updateCellChatDB` was **rejected, not overlooked**: it would
discard 22 columns including `agonist`, `antagonist`, `co_A_receptor` and `co_I_receptor` —
precisely the terms `computeCommunProb` multiplies into `dataRavg` — while adding zero interactions.
**Degrading the method to re-supply a database it already ships is not a second tier.**

**Widening the DB does not change any individual pair's probability.** The requested LRIs come out
bit-identical between tiers, because the Hill function for a pair depends only on that pair's
ligand/receptor expression. **The tier changes *coverage* and *pathway-level aggregation*, nothing
else.**

### Input

| requirement | detail |
|---|---|
| `data.input` | **normalized** data, genes × cells (V1:49). `adata.X` is already `log1p(CP10K)` — `prepare_gbm_input.py` asserts it and records the implied library size (median exactly 10,000.0). `layers['counts']` is **not** used |
| `meta` | `labels` ← `obs['cell_type']` (9 types), and a `samples` column ← `obs['tma_id']`. Both factors; **unused levels abort `computeCommunProb`** |
| `coordinates` | `obsm['spatial']`, already µm |
| `spatial.factors` | `ratio = 1`, `tol = 5`, **one row per sample in `levels(samples)` order** |

**Object layout.** The vignettes are explicit that several sections of *one* condition go into
**one** object via `samples`, but comparison across conditions needs a **separate object per
condition**. GBM is therefore two objects: `high` (7 cores — tma_id 1,3,5,8,10,11,13; 79,998 cells)
and `low` (6 cores — 2,4,6,9,12,14; 20,199 cells). All 9 cell types occur in all 13 cores, so
`mergeCellChat` and the **functional**-similarity manifold analysis are both applicable.

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
| `quant/requested_lr_status.csv`, `requested_lr_in_DEA.csv`, `<c>_requested_lr_plot_status.csv` | — | where the requested LRIs landed, **including when they landed nowhere** |
| `run_manifest.json` | — | DB, every parameter, per-condition counts, wall time, versions, git SHA |

### Image outputs

Every plot goes through one `save_all_formats()` in `cellchat_io.R` → png + pdf + svg (CLAUDE.md
rule). Per condition: `plots/<cond>/aggregate/` (circle count/weight, heatmap count/weight, one
circle per sender), `plots/<cond>/pathways/<pathway>/` (circle, chord, hierarchy, heatmap, LR
contribution, signalingRole network, chord_cell, gene-expression violin, plus circle+chord for every
enriched L-R pair), `plots/<cond>/spatial/core<id>/<pathway>/` (spatial network and the
incoming-weighted variant — **every pathway × every core**), `plots/<cond>/systems/` (role scatter,
outgoing/incoming/all role heatmaps, selectK curve, river + dot pattern plots, functional and
structural embeddings). Cross-condition: `plots/comparison/` (compareInteractions, diffInteraction
count/weight, differential heatmaps, per-condition circles, role scatter, `signalingChanges` scatter
for all 9 cell types, pairwise functional+structural embeddings, rankSimilarity, rankNet
stacked/unstacked, role heatmaps side by side, bubble comparison/increased/decreased, bubble+chord
for up- and down-regulated pairs, enrichment wordclouds, split gene-expression violins).

Per-LR plots are kept in **two separate trees**, as the benchmark requires:
`plots/top_lr/<cond>/<LR>/` for the pairs CellChat itself ranks highest by summed communication
probability, and `plots/requested_lr/<cond>/<LR>/` for **GRN_SORT1 and ANXA1_FPR1 regardless of
rank** (circle, chord, bubble, pathway contribution, and per core: gene spatial maps, continuous LR
spatial map, binary LR spatial map).

**Totals: `default` 1,632 plots and `cellchatdb2` 3,666 plots.**

> ⚠️ **The oversized SVGs were deleted 2026-08-12.** Written in all three formats these 5,298 plots
> occupied **18.5 GB**, of which **8.7 GB was 2,518 SVGs larger than 1 MB** — almost all of them the
> per-pathway × per-core `spatial/` maps, which are 20k–80k-point scatters whose SVG is neither
> openable nor usefully editable. Those 2,518 files were deleted; the **2,780 SVGs under 1 MB are
> intact** (circle, chord, heatmap, bubble, role and comparison figures — the ones that are genuinely
> editable vector art for a manuscript), as is **every PNG and PDF**. Regenerate with
> `bash scripts/comparators/cellchat/run_all_gbm.sh` if a large SVG is ever needed.

**Not produced, with reasons** — written to `plots_not_produced.txt` on every run, never dropped
silently. Exactly **one** entry remains in each tier: `netVisual_chord_gene` at *all* sources × *all*
targets, which is not renderable — the high-grade network has **161 distinct ligand/receptor
sectors** and circlize fails to allocate them at every `small.gap` tried (1, 0.5, 0.2, 0.1). The
tutorial never draws it at that scope either (both VB and VC use a single sender), so the covering
set is **one chord per sender**, which renders for all 9 cell types in both conditions. Also
deliberately not run: `runCellChatApp` (interactive Shiny, not a file artifact);
`netVisual_embeddingZoomIn` (attempted, skipped when a group is too small); `projectData`/PPI
smoothing (commented out in both spatial vignettes); `mergeInteractions` coarse-cell-type regrouping
(our 9 labels are already coarse and no defensible 3-way grouping exists without asking); and any
network view of a requested LR with **no significant link**, which is itself recorded as a result.

### Multi-sample / differential mode

**Native, and used in both forms.** *Within* a condition, all cores live in one object as
`meta$samples` levels and CellChat computes cell-group distances per sample then averages them — so
**physically separate punches are never treated as neighbours**. *Across* conditions,
`mergeCellChat` + the whole of VC: `compareInteractions`, `netVisual_diffInteraction`,
`rankNet(do.stat=TRUE)`, `netAnalysis_signalingChanges_scatter`, pairwise manifold learning, and the
presto-backed cross-condition DEA.

### Gotchas

1. **`spatial.factors` is indexed positionally by `levels(meta$samples)`** — one row per sample in
   level order. **Wrong order silently applies the wrong µm conversion. No vignette says so.**
2. **`contact.range` is mandatory when `contact.dependent = TRUE`.**
3. **`scale.distance` is validated, not defaulted** — with `distance.use = TRUE`, if
   `min(d × scale.distance) < 1` CellChat aborts and prints the value to use. **Never copy `0.01`
   from the Visium tutorial.**
4. **Unused factor levels abort the run**; `droplevels` is required after subsetting to one grade.
5. **presto changes the numbers, not just the speed.** `do.fast = TRUE` is the default and falls
   back to `stats::wilcox.test` **silently** if presto is missing, giving larger logFC while VC's
   `thresh.fc = 0.05` was tuned *for* presto. presto 1.0.0 is installed, asserted at startup and
   recorded in the manifest.
6. **`computeCellDistance` is dense O(N²)** — it needs ~51 GB at 79,998 cells and returns NA.
   Replaced by `BiocNeighbors::findKNN(k = 1)`, which gives the identical quantity in O(N log N).
7. **`sample.use` is mandatory** for `netVisual_aggregate(layout="spatial")` and
   `spatialFeaturePlot` in multi-sample mode.
8. **`netEmbedding` prefers python `umap-learn` via reticulate**; we pass `umap.method = "uwot"` so
   the run does not depend on a reticulate python.
9. **`selectK` cannot be made to work inside a long plotting session.** It calls
   `NMF::nmfEstimateRank` with NMF's default parallel foreach backend and exposes no override; in a
   session that has already run Seurat/ComplexHeatmap work, every run dies with "All the runs
   produced an error", while the *identical call succeeds in a clean session*. We compute the
   measures with `.pbackend = "seq"` and redraw the curve from them.
10. **`identifyCommunicationPatterns` requires `k`** and the vignette picks it by eye. We apply the
    vignette's stated rule — "the one at which Cophenetic and Silhouette begin to drop suddenly" —
    programmatically to the same measures, and persist the measures so the choice is auditable.
    NMF's rank must also stay below both matrix dimensions, so `k.range` is capped at
    `min(10, n_celltypes-1, n_pathways-1)`.
11. **The `future.rng.onMisuse` warnings are benign.** CellChat draws its permutation matrix once
    under `set.seed(seed.use)` *before* the parallel loop. Two runs with identical arguments produced
    **bit-identical** `net_full`, `net_count` and `net_weight` — verified, not assumed.
12. **`netAnalysis_river` needs `ggalluvial` ATTACHED, not just installed.** It builds a ggalluvial
    plot with `stat = "stratum"`, and ggplot2 resolves stats by name off the search path, so
    `ggalluvial::` is not enough — it dies with `Can't find stat called "stratum"`. **Cost us all 8
    river plots on the first pass.**
13. **`netVisual_chord_gene` has a hard scale ceiling** — at 161 sectors circlize cannot allocate the
    layout at any `small.gap`. Restricting to one sender always works.

### Deviations from the tutorial

| item | tutorial | ours | why |
|---|---|---|---|
| `contact.range` | `100` (V1:182, V2:158) | **`10`** | 100 is the 10X Visium spot pitch. VF:39 / V1:60 / `man/computeCommunProb.Rd` pin `10` for single-cell-resolution platforms. **Measured on our data: median nearest-neighbour distance 11.5 µm (low) and 7.45 µm (high)** — so 10 µm is the right scale and 100 µm would have made "contact" meaningless |
| `spatial.factors` | Visium `ratio = 65/spot_diameter_fullres`, `tol = 32.5` | `ratio = 1`, `tol = 5` | VF:78-83, the FAQ's own **Xenium** row: coordinates already µm, `spot.size = 10` = typical human cell |
| normalization | `GetAssayData(slot="data", assay="SCT")` | `adata.X` as-is | already `log1p(CP10K)`, asserted at export; re-normalizing would log twice |
| `distance.use` | V1 `TRUE` + `scale.distance=0.01`; V2 `FALSE` + `NULL` | **`FALSE`** | V2 governs our data shape; V1's `0.01` is Visium-specific and would abort (gotcha 3) |
| `variable.both` | V1 passes `F`, V2 passes nothing (`TRUE`) | `TRUE` | V2 is the governing vignette and `TRUE` is the package default |
| `umap.method` | not passed (→ `umap-learn`) | `"uwot"` | avoids a reticulate python dependency; the package's own documented alternative |
| `sources.use`/`targets.use` | hardcoded indices (`4`, `5:11`) | all 9 cell types | the indices are specific to VC's 12-cluster skin dataset |
| `pathways.show` | one hand-picked pathway | **every** pathway in `netP$pathways` | benchmark invariant: produce every plot the standard workflow can produce |
| `selectK` | called directly | measures recomputed with `.pbackend="seq"`, curve redrawn | gotcha 9 |
| `k` for patterns | read off the curve by eye | same rule applied programmatically, measures persisted | reproducibility |
| `computeCellDistance` | used as-is | `BiocNeighbors::findKNN(k=1)` | gotcha 6 — dense O(N²) OOMs at 79,998 cells |
| `mergeInteractions` (VC option D) | 12 clusters → 3 coarse types | skipped | our 9 labels are already coarse; no defensible 3-way grouping without asking |
| `netVisual_chord_gene` scope | one sender | all-sources attempted, **plus one chord per sender** | gotcha 13; per-sender is the tutorial's own scope and covers all 9 cell types |
| env: CRAN source | CRAN HEAD | dated snapshot `2024-06-01` | the env is R 4.3.3 and current CRAN sources need ≥ 4.4 |
| env: igraph | stock build | `--disable-graphml` + `xml2-config` shadowed | base anaconda's `xml2-config` leaks in and igraph links `libxml2.2.dylib`, absent from the env |
| smoke test | on the tutorial's demo data | on a 819-cell GBM core | **CellChat ships no demo expression data** — only `CellChatDB.*.rda` and `PPI.*.rda`; the vignettes load from the author's local OneDrive paths |

### Runs on our data

| dataset | tier | status | key numbers |
|---|---|---|---|
| GBM | `default` | ✅ | DB 1,280 interactions / 158 pathways (Secreted only); 525 signaling genes on the panel. **low**: 20,199 cells / 6 cores, 355 LR tested, **270 significant links**, 23 pathways, 511 s. **high**: 79,998 cells / 7 cores, 382 LR tested, **398 significant links**, 25 pathways, 778 s |
| GBM | `cellchatdb2` | ✅ | DB 2,239 interactions / 252 pathways (Secreted + ECM + Contact); 722 signaling genes. **low**: 568 LR tested, **927 significant links**, 61 pathways, 1,257 s. **high**: 601 LR tested, **1,191 significant links**, 67 pathways, 1,854 s |
| GBM | db audit | ✅ | Jaccard 1.0000 vs the repo CSV | `cellchat/db_audit/` |
| LUAD | — | ❌ not run | |

Aggregate comparison (`default`): low 270 links / total strength 3.074 / 23 pathways versus high 398
links / 3.227 / 25 pathways — i.e. **high grade has 47% more significant links but only 5% more total
interaction strength, so the difference is breadth, not intensity.** The presto DEA found **170 up**
and **65 down** L-R pairs in high grade; the up list is led by BMP (32 pairs), COMPLEMENT (26), PDGF
(17), TGFb (15) and IGF (15), the down list by COMPLEMENT (12), GAS (12), GRN (10), FGF (9) and PDGF
(7). **COMPLEMENT, PDGF and GRN appearing on both sides is the effect VC warns about**: the DEA is
run per cell group, so one pair can be up in one group and down in another.

### Requested LRIs — CellChat recovers both arms, and only in high grade

Both are `Secreted Signaling`, so both are in **both** tiers, and both survive
`identifyOverExpressedInteractions` in both conditions. **Tier makes no difference to the numbers.**

| LR | condition | significant cell-type pairs | max prob | ALARMIST direction | that direction's prob |
|---|---|---|---|---|---|
| `GRN_SORT1` (pathway GRN) | low | 15 | 0.04775 | mGAM → MES-like | **absent (n.s.)** |
| `GRN_SORT1` | high | **25** | 0.04434 | mGAM → MES-like | **0.01128, p < 0.001** |
| `ANXA1_FPR1` (pathway ANNEXIN) | low | 8 | 0.02521 | MES-like → mGAM | 0.00140 |
| `ANXA1_FPR1` | high | **10** | 0.02493 | MES-like → mGAM | **0.01624, p < 0.001** (3rd-ranked pair of that LR) |

**The complete bidirectional mGAM ⇄ MES-like loop is significant only in high grade**: in low grade
the GRN→SORT1 arm is not called at all and the ANXA1→FPR1 arm is ~12× weaker. **This is an
independent corroboration of ALARMIST motif 1 being grade-associated, from a method with a
completely different inference target (cell-type pairs, permutation null) and no knowledge of the
motif decomposition.**

Two caveats. **CellChat's own ranking does not put the ALARMIST direction on top**: the strongest
GRN_SORT1 pair in both conditions is mGAM → Glial-Neuronal, and the strongest ANXA1_FPR1 pairs are
Vascular → mGAM and Lymphoid → mGAM. The ALARMIST direction is present and significant but not the
maximum. And in the cross-condition DEA `GRN_SORT1` appears in **both** `net.up` and `net.down` (the
per-cell-group artefact above), while `ANXA1_FPR1` appears only in `net.up`.

Also of note: ANXA1_FPR1's receivers in high grade are **exactly `mGAM` and `non-mGAM`** and nothing
else, consistent with FPR1 being mGAM-restricted on this panel.

### Reproducing the GBM run

```bash
source scripts/comparators/cellchat/activate_env.sh
bash scripts/comparators/cellchat/run_all_gbm.sh all
```

`prepare_gbm_input.py` (env `bptf`) regenerates the input tree; `run_all_gbm.sh` runs the DB audit,
then both tiers (inference + plots). Env: `env.lock.yml` (conda) + `r_packages.lock.csv` (249 R
packages with versions); `install_env.R` rebuilds the R library.

### Methods paragraph

> Cell–cell communication was inferred with CellChat v2.2.0.9001 (Jin et al.) in R 4.3.3. For each
> tumour grade a separate CellChat object was created with `createCellChat(datatype = "spatial")`
> from log-normalized Xenium expression, the nine annotated cell types as `labels`, and the
> constituent TMA cores as `samples`, so that cell-group distances are computed within each core and
> averaged across cores; spatial factors were set to the authors' recommended values for
> single-cell-resolution platforms (`ratio = 1`, `tol = 5` µm). The ligand–receptor database was
> CellChatDB v2 (human), used both as the tutorial's `Secreted Signaling` subset (1,280 interactions)
> and as the full protein-coding database (`subsetDB`, 2,239 interactions). Over-expressed genes and
> interactions were identified with `identifyOverExpressedGenes` and
> `identifyOverExpressedInteractions` at package defaults (Wilcoxon test via presto), and
> communication probabilities were computed with `computeCommunProb` using a 10% truncated mean
> (`type = "truncatedMean"`, `trim = 0.1`), a 250 µm diffusion range as a hard spatial filter
> (`distance.use = FALSE`), a 10 µm contact range, and a 100-permutation label-shuffling test
> (`seed.use = 1`). Communications were filtered with `filterCommunication(min.cells = 10)`,
> aggregated to pathway level with `computeCommunProbPathway` and to cell-type level with
> `aggregateNet`, and network centralities computed with `netAnalysis_computeCentrality`. High- and
> low-grade objects were combined with `mergeCellChat` and contrasted with `compareInteractions`,
> `netVisual_diffInteraction`, `rankNet` (paired Wilcoxon), `netAnalysis_signalingChanges_scatter`
> and a cross-condition differential expression analysis (`identifyOverExpressedGenes` with
> `group.dataset = "datasets"`, `thresh.fc = 0.05`) mapped back onto the inferred communications
> with `netMappingDEG`.

---

## Figure 6 / Figure 7 — the synthesis layer

**This is what the seven method sections above are for**, and it was documented nowhere until
2026-08-12. Both figures **read only files that already exist on disk; they run no inference and
refit nothing.** Env: `bptf`. Outputs: `results/comparators/_benchmark/`.

| script | reads | writes |
|---|---|---|
| `figure6_gather.py` | every method's persisted ranking / significance / manifest | `figure6/panel_{a,b,c,d,e,f}*.csv` + `panel_f_meta.json` |
| `figure6_plot.py` | those CSVs only | `figure6/figure6_comparator_benchmark.{png,pdf,svg}` |
| `figure6_supp_density.py` | `panel_f_*.csv` | `figure6/figure6_supp_density_collinearity.{png,pdf,svg}` |
| `figure7_gather_spatial.py` | each method's per-cell / per-spot map for the two arms, all 13 punches | `figure7/spatial_maps_all_cores.npz`, `provenance_all_cores.json`, `within_method_rho_all_cores.json` |
| `figure7_plot_spatial.py` | that npz only | `figure7/punches/figure7_core<N>.{png,pdf,svg}` ×13 + `figure7_rho_summary.*` |

`figure6_plot.py` also **exports `set_style()` and `save_all_formats()`, which the other two figure
scripts import from it** — it is the style authority for this figure family, and since 2026-08-12
it delegates to `scripts/comparators/_common/plotting.py`.

### Figure 6 — what each panel asserts

| panel | claim | source |
|---|---|---|
| **a** | **every comparator detects both arms of motif 1** — 16 rows, one per (method × arm), each carrying its rank, its denominator, its percentile and its `rank_provenance` (`native` where the method ships a ranking, `sorted` where we sorted its own statistic) | `panel_a_recovery.csv` |
| **b** | **…but in tested sets whose sizes differ ~30-fold**, and each row states *why* — stLearn 526 because it cannot encode heteromeric complexes; SpatialDM 1,546 per core because the expression filter is re-run each core; CytoSignal 895 of 1,088 scored | `panel_b_denominators.csv` |
| **c** | **no comparator puts the two arms in one data-derived object.** Each method is typed `flat` / `curated` / `learned`, and both arms' object paths are given so the claim is checkable | `panel_c_objects.csv` |
| **d** | **why: at the cell the two arms barely co-occur, at the 50 µm patch they do** — the unit demonstration, 2 rows (cell: both 0.095%, Pearson 0.0177, enrichment 3.17×; patch: both 1.235%, Pearson 0.4562, enrichment 14.23×) | `panel_d_unit.csv` |
| **e** | **and only two comparators ship a between-condition test at the punch.** All rows are at the **TMA-punch** replicate unit, n = 13 | `panel_e_grade.csv` |
| **supp** | **grade and cellularity are collinear in this TMA** — ρ = 0.784, p = 0.0015; high-grade punches hold **11,428** cells on average against **3,366**; the 7-vs-6 Mann-Whitney floor is **p = 0.0011655** | `panel_f_*.csv`, `panel_f_meta.json` |

The supplement's own docstring states its limit correctly and it should be quoted with the figure:
**it does not show the grade result is an artifact.** At ρ = 0.78 and n = 13 grade and cellularity
cannot be separated, and residualising on density removes the grade signal for **every** motif —
which is what collinearity does regardless of the truth.

### Figure 7 — one spatial panel per punch

Columns are methods, rows are the two arms, and **ALARMIST occupies a single cell spanning both
rows, because for it the two arms are not two objects.** Three properties are built in rather than
annotated:

- **Colour is the percentile within each map.** Values are **not comparable across methods** —
  LRscore / co-expression / Moran z / transport mass / mechanism score / cosine are different
  quantities on different supports.
- **ρ under each column is the Spearman correlation between that method's OWN two maps, in that
  method's OWN units, on that punch** — a *within*-method statistic, so it carries none of the
  cross-method scale problems the maps have.
- **Whole-slide methods are cropped to a punch, never recomputed on it** (`stlearn_core_assignment`,
  `whole_slide_methods_cropped_not_refit` in `provenance_all_cores.json`).

Two absences are recorded as results, not gaps: **CellChat has no spatial output at all** (its
finest unit is a cell-type pair within a section), and **SpatialDM contributes a map only where the
pair was globally selected** — GRN in 7 cores, ANXA1 in 6, both in only 1
(`spatialdm_pairs_not_globally_selected`).

---

## Code inventory

`scripts/comparators/` is a flat tree of per-method directories plus:

| path | role |
|---|---|
| `METHODS.md` | this file |
| `_common/plotting.py` | **the single figure saver and style setter for this tree** (added 2026-08-12). `apply_publication_style(**overrides)` + `save_all_formats(fig, stem, dpi=…, close=…, verbose=…)`. Seven scripts had each re-declared their own, with three saver signatures and three default dpi. **Anything new here imports from it — do not write an eighth copy.** Deliberately *not* in `src/alarmist/plotting/`, per CLAUDE.md |
| `figure6_*.py`, `figure7_*.py` | the synthesis layer (above) |
| `_archive/METHODS.2026-08-12.md` | the verbatim 4,370-line predecessor of this file |
| `_archive/cytosignal/` | 7 one-off CytoSignal comparison scripts, superseded by `figure6_*`/`figure7_*` but retained because their outputs are still under `cytosignal/GBM/cellchatdb2/run_full/plots/` |

Per method: `NOTES.md` (call contract), `DEVIATIONS.md`, `env.lock.yml`, `run_<method>.{py,R}`, a
plotting script, and an `activate_env.sh` for the R methods (because degraded conda makes
`conda activate` fall through to system R).

| | NOTES.md | DEVIATIONS.md | env.lock.yml | activate_env.sh |
|---|---|---|---|---|
| cytosignal | ❌ | ❌ | ❌ | ✅ |
| stlearn | ✅ | ✅ | ✅ | n/a |
| spatialdm | ✅ | ✅ | ✅ | n/a |
| commot | ✅ | ✅ | ✅ | n/a |
| liana | ✅ | ✅ | ✅ | n/a |
| niches | ✅ | ✅ | ✅ | ✅ |
| cellchat | ✅ | (in NOTES.md) | ✅ + `r_packages.lock.csv` | ✅ |

---

## Open issues

Consolidated from every section. Ordered by how much they would change a claim.

### Scope not yet covered

- ❌ **LUAD is essentially untouched.** `SKILL.md:47-48` names GBM **and** LUAD. Only CytoSignal ever
  ran on it, on a 2.5 mm crop of P21_LUAD; the AIS↔LUAD contrast across `P{17,21}` has not been
  attempted by any method. Full-section CytoSignal is blocked at ~57 GB and has a portable bundle
  waiting for a ≥64 GB node.
- ❌ **The `default` tier exists for 3 of 7 methods** (SpatialDM, LIANA+, CellChat). `SKILL.md:51-54`
  asks for both tiers per method. Where it *was* run, the answer is consistent and reassuring — the
  LR resource changes which pairs are tested and the ranking denominator, **not any per-pair
  statistic** (SpatialDM Pearson r = 1.000 on shared pairs; LIANA max |Δ| = 0.0 on 23,787 rows;
  CellChat bit-identical between tiers) — so the missing four are a completeness gap rather than a
  suspected confound.

### Known defects, not yet fixed

- ❌ **`compare_programs_to_alarmist.py` NMF index bug** — reads 8 numeric columns from a rank-7
  `NMF_W_factor_scores.csv` written with `index=False`, so `cell_id` is scored as a factor. **That
  row of `vs_alarmist/comparison_summary.json` is unverified and must not be quoted.** The three
  MOFA-Flex rows are unaffected.
- ❌ **`nmf_factor_correlation.csv` has no producing script in the repo.** Values verified correct;
  **the table cannot currently be regenerated.**
- ⚠️ **SpatialDM's `default` tier `plots/` still holds 143 blank 7,544-byte PNGs** from the first
  plotting pass. Cleaned from the `cellchatdb2` tree, not from this one. Read `plots_full/`.
- ⚠️ **COMMOT `GBM/impact/` is a PARTIAL run** (9 of 13 cores, one variant, no manifest, no controls),
  marked on disk with `PARTIAL_RUN_DO_NOT_USE_AS_FINAL.md`. **Nothing in it is a finished result.**
- ⚠️ **`reports/liana_plus_GBM_cellchatdb2/liana_plus_GBM_methods.html` is stale** relative to this
  document. It is built from `_liana_report_sections.json` by `build_liana_report.py`.

### Contract deviations awaiting sign-off

- ❌ **CD-1 — LIANA's bandwidth (13.1454 µm) is derived from ALARMIST's 50 µm patch**, violating
  `SKILL.md:45-46`. Kept by user decision after the evidence review, **not resolved by it**. It is
  load-bearing: it alone moves the inflow NMF rank 11 → 7. Exempt: LRIC, LR-MISTy. Inherited:
  MOFA-Flex.
- ❌ **CD-2 — no native multi-sample / differential mode for LIANA.** The package has none for the
  spatial branches; the punch-level Mann-Whitney is a documented hand-rolled substitute, null in
  five places, and does **not** satisfy the requirement.
- ❌ **No K-sweep for MOFA-Flex.** `n_factors = 20` is a binding ceiling in all three fits and nobody
  has looked for where the active count saturates.
- ❌ **No prevalence normalisation of the MOFA-Flex loadings.** ALARMIST has `V*`; LIANA has no
  equivalent, so the cosine comparison normalises one side only.

### Repository hygiene

- 🔴 **`scripts/comparators/` is not tracked in git.** `git ls-files` returns nothing and
  `git check-ignore` says it is not ignored either — it is simply never added. Every
  `run_manifest.json` records `git_sha: 95208de`, which pins the **package**, not the comparator
  scripts that produced 44 GB of results. An edit can be confirmed to have the intended content, but
  there is no way to prove no other line moved. **This is the highest-value single fix here.**
- ⚠️ **CytoSignal has no `NOTES.md`, `DEVIATIONS.md` or `env.lock.yml`**, although `SKILL.md:71`
  names it the reference implementation the other methods were matched against.
