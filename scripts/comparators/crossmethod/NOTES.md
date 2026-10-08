# crossmethod — the 1-D classifier benchmark

**Question (Wesley, 2026-09-10 meeting).** Use the pathologist's tissue annotation as the
label. Fit the simplest possible classifier — a one-dimensional logistic regression — on one
feature at a time. ALARMIST contributes motif 24 and motif 10; every comparator contributes
each of its individual LR pairs. That gives, per method, the whole *distribution* of "how well
can a single LR pair predict this label", and the question is where ALARMIST's motif falls in it.

**Nothing here re-runs any comparator.** Every score matrix is read exactly as its own run
wrote it. The only new computation on a comparator's output is CytoSignal's `.rds → .f32`
re-serialisation (`02_export_cytosignal.R`), which changes container, not values.

## Files

| file | what it does |
|---|---|
| `_config.py` | host paths, the two label definitions, per-method input locations |
| `auc.py` | the scoring primitive; two implementations + their cross-check |
| `01_build_labels.py` | `obs` → per-section label table + ALARMIST features + reference features |
| `02_export_cytosignal.R` (+`.sbatch`) | banked `score_*.rds` → flat float32 for the Python side |
| `03_run_auc.py` (+`.sbatch`) | every method × section × label → one tidy row per feature |
| `04_summarise.py` | where the motif falls in each method's distribution |
| `05_plot.py` | the two figures |

Outputs: `/data1/tanseyw/projects/fanj2/results/_crossmethod/label_auc/`.

## The design decisions that matter

**The label is model-independent.** `obs['tissue_layer']` was drawn by a professional
pathologist and is not derived from ALARMIST, which is the only reason this comparison is
worth anything. Five levels; two binarisations, both run:

| label | positive | negative | excluded |
|---|---|---|---|
| `clean` (main figure) | `Tumor_region` | `Normal_peri`, `Normal_distal` | both interface layers |
| `full` (sensitivity) | `Tumor_region`, `Interface_tumor` | `Normal_peri`, `Normal_distal`, `Interface_normal` | none |

**AUC *is* the 1-D logistic regression's AUC.** With one predictor the fitted model is monotone
in that predictor, so the ROC of the fitted probabilities equals the rank ROC of the raw
feature. We score by rank AUC — and *prove* the identity numerically against sklearn rather
than asserting it. Recorded in `checks.json` as `logistic_identity`.
Measured 2026-09-10: **max |difference| = 1.71e-08**.

**Direction-agnostic.** Everything is ranked on `auc_abs = max(AUC, 1 − AUC)`, so an LR pair
that is strongly *down* in tumour gets full credit. This makes the test harder for ALARMIST,
not easier.

**Every comparator is paired to its own denominator.** Methods drop cells for their own
reasons — CytoSignal's `removeLowQuality` takes out ~13%, SpatialDM leaves a handful with
non-finite local z, stLearn scores a grid. So the motif AUC is *recomputed on exactly the units
each method retained* and written as `alarmist@<method>`; `04_summarise.py` refuses to compare a
method against a paired row with a different `n_units`. Comparing everything against one
all-cells motif number would have put different denominators on the two sides.

**Scale is not harmonised.** `comparator-benchmark` SKILL.md:45-46. Each method keeps its own
kernel/neighbourhood, each is compared only against its own distribution, and LR-pair *counts*
are never compared across methods.

## Per-method contract

| method | unit | features | source file | note |
|---|---|---|---|---|
| ALARMIST | cell | 2 | `obs['motif_{10,24}_loading']` | the claim |
| SpatialDM | cell | 579 | `data/local_z.npz` | **globally selected pairs only** (579 of 1,692) |
| CytoSignal | cell | 1,083 | `quant/score_{diffusion,contact}_*.rds` | post-QC subset; `Raw` slot exported but not scored |
| LIANA | cell | 9,313 | `cellchatdb2_inflow/data/inflow_scores.npz` | features are `celltype^ligand^receptor`, finer than an LR pair |
| stLearn | **spot** | 532 | `data/spot_lr_scores.csv.gz` | ~51 µm grid; spots labelled by majority vote |
| CellChat | cell-type-constant | 302 | `quant/{AIS,LUAD}_net_full.csv.gz` | **no per-cell score exists** — see below |
| MOFA-Flex | cell | 19 | `mofaflex_inflow_joint/data/factor_scores.csv.gz` | the equal-complexity baseline |
| reference | cell | 21 | built here | local density, total counts, cell-type indicators |

**CellChat is a unit-of-analysis result, not a defect.** Its output is
(LR pair, sender type, receiver type) within a *group* (AIS or LUAD), so within a section it
carries no spatial variation at all. The most generous per-cell feature derivable from it is the
communication probability attached to a cell's own type, built both as incoming and outgoing so
it gets twice the features and direction-agnostic credit on each. Such a feature takes at most
19 distinct values, so its AUC is bounded by what cell-type identity alone achieves — which the
`best_cell_type_indicator` reference measures directly. Say this; do not present it as CellChat
underperforming.

**The reference features are an addition, not part of the takeaway.** `local_density_80um`
(cells within ALARMIST's 80 µm patch radius) is there because the first thing a reviewer will
say is "a feature smoothed over 80 µm will track a large contiguous region regardless of what it
measures". `best_cell_type_indicator` bounds what cell identity alone can do. They are labelled
`kind == 'reference'` and can be dropped without touching anything else.

## Traps hit here

- **`local_z.npz` has no names.** It was saved from a bare ndarray, so `pairs`/`cells`
  degenerated to `'0','1','2',…` (`run_spatialdm.py:258-260`). Real names come from
  `local_n_spots.csv` (index = `global_res[selected]`, verified equal) and `cell_meta.csv`
  (row order). The loader asserts the shapes match rather than assuming.
- **Non-finite local z is a property of *cells*, not pairs.** Every affected pair is NaN on
  exactly the same handful of cells (10 of 182,378 in P17_AIS, 3 in P21_AIS — isolated cells
  with no usable neighbourhood). Dropping pairs instead would have thrown away 440 of 579 pairs
  over 10 cells. Drop the cells.
- **`obs_names` are `'0','1','2',…`** in `adata_with_region.h5ad` and carry no identity. The
  join key is `sample_id + '_' + cell_id`, which reproduces every comparator's cell ids exactly
  (182,378/182,378 for P17_AIS, zero unmatched either way).
- **CytoSignal scores are 70-80 % dense** despite being `dgCMatrix`, so MatrixMarket would be
  ~100 M lines for one section. Exported as flat float32, written column by column, i.e.
  **Fortran order** — the Python side reads with `order='F'`.
