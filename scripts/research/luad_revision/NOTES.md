# NOTES — `inventory_anndata.py`

**What this is.** A structural inventory of the four prepped LUAD/AIS Xenium AnnData objects. It is
descriptive only: no normalization, no module scores, no differential testing, no plots. Output goes
to `inventory.json`; the written-up numbers are in `../../../RESULTS.md` § 1.

**Contract it follows.** There is no upstream tutorial for this — it is a repo-local audit. The
things it is designed to answer are the ones later spatial analyses silently get wrong:

1. *Is `X` counts or normalized?* Decided on values (min, max, integer-valued, per-cell sum median),
   never on the layer's name. Both `X` and `layers['counts']` are profiled.
2. *Are the spatial coordinates in µm?* Decided on the median nearest-neighbour distance between
   cells, which must land in the single-digit-to-low-teens µm range for a cell-resolution platform.
   Extent alone cannot distinguish µm from pixels.
3. *Do the sections share a coordinate frame?* Decided by comparing per-section bounding boxes.
4. *What exactly are the cell-type strings?* Read out of the categorical, verbatim, with counts, so
   nothing downstream has to guess at `pDC` vs `PDC` vs `Myeloid_C11_pDC`.

**Method notes.**

- Opened with `h5py`, never `anndata.read_h5ad`. No matrix is ever materialized: min/max and
  integer-ness stream the CSR `data` array in 8M-element chunks, and per-cell sums use
  `np.add.reduceat` over 40k-row blocks of `indptr`.
- Per-cell-type gene statistics (pDC, cDC) make one blocked pass over `layers['counts']`, touching
  only rows in the mask.
- Nearest-neighbour distances use `scipy.spatial.cKDTree`. The tree is built on **all** cells; only
  the query points are subsampled, to 200,000, seeded `default_rng(0)`, for sections above that size.
  P17_AIS (182,378) uses all cells.
- "Genes detected ≥5%" means ≥1 raw count in ≥5% of that cell type's cells in that section. It is a
  detection rate, not an expression level, and is confounded with sequencing depth.

**Scope choice.** The `prepped` files were inventoried rather than the original `xenium_mm.h5ad`,
by explicit instruction on 2026-08-24. The two differ: prepped drops unannotated cells (~20% of each
section) and replaces `X` with log1p. See RESULTS.md § 1.0 for both file sets and their sizes.

**Cost.** ~1 minute wall clock for all four sections on `iscf031`, I/O bound, well under 4 GB RSS.

## pdc_neighborhood_composition.py (2026-08-26)

Cell-type composition of the 30 µm neighbourhood of every pDC, against (a) the whole-tissue
composition and (b) the mean neighbourhood of a random reference cell. Contexts: each of the four
sections, AIS pooled, LUAD pooled, all four.

- Input: `/data1/tanseyw/projects/fanj2/alarmist_luad/adata_with_region.h5ad`, `obs/cell_type`,
  `obs/sample_id`, `obsm/spatial` only. No ALARMIST output is read; this analysis is model-free.
- Neighbour counts are **pooled over pDC** (a cell in range of two pDC counts twice). The
  alternative -- the mean of per-pDC fractions -- is in the JSON as `nbhd_frac_mean_per_pdc` and
  agrees to within ~0.02 everywhere.
- The centre cell is removed from its own count, so `pDC` enrichment measures pDC--pDC clumping and
  nothing else.
- Two nulls. The cell-level one (2,000 draws) is anticonservative because pDC are clumped; the
  clump-matched one (200 draws) reproduces the observed DBSCAN(eps=30 µm) clump-size distribution by
  expanding random anchor cells to the matching sizes. Report the clump-matched p; its floor is
  1/200 = 0.005.
- Runtime ~90 s on a login node, single core. Reference pool 30,000 cells/section.

## vasendo_neighborhood_reciprocal.py + far_stratum_composition.py (2026-08-26)

Answers "does §8 contradict §5.3". Three parts, all model-free:
- A: §8 with the centres swapped to Vas_Endo. Cell-level permutation only (2,000 draws) -- the
  clump-matched null of §8 is not meaningful for 137,915 endothelial cells forming contiguous
  vessels, so the p-values here are anticonservative and are not quoted in RESULTS.md.
- B: absolute mean neighbour counts per cell type (capped at 8,000 cells/type/section) against a
  30,000-cell random reference. This is the metric that is comparable to §5.3; §8's fractions are
  compositional and cannot be.
- C: the identity E1/E2 = mean-nbrs-Vas_Endo / mean-nbrs-pDC, verified to 3 dp with identical pair
  counts from both sides. Consequence: the reciprocal direction is NOT independent evidence.
- far_stratum_composition.py: shows the ">30 um from a vessel" stratum is 26-64% of non-endothelial
  cells and is tumour-core/airway-dominated, which is why §5.3's ratio is inflated.

Runtime ~4 min (the per-type sampling dominates).

Label trap (found 2026-08-26, fixed): the grey reference points are cell types **as the centre
cell**, and an early version labelled "whichever grey point is highest". In panel a that is Vas_Endo
and reads naturally; in panel b it is Myeloid_other, which made the panel look as if it were about
Myeloid_other rather than pDC. Both panels now label the target type's own self-adjacency as the
ceiling, and any cell type with n < 50 in that section is drawn hollow and never labelled
(Myeloid_other is n = 45 in P21_AIS). The numbers were never affected -- the magenta point has
always been the pDC-centred value, re-verified against an independent recompute:
mean T within 30 um of a pDC = 3.385 / 7.071 / 2.841 / 4.188 across the four sections.

Redesigned 2026-08-26 on user feedback: the earlier strip-plot version (sections on x, 19 grey dots
per column) was unreadable and its "label the top grey point" annotation misread as if the panel
were about Myeloid_other. Now one row per cell type, one marker per section, tab20 colours shared
with the composition figures, pDC enlarged and its label bold pink. Rows are ordered by the
four-section mean, so the row order is itself the ranking. Titles state the question directly.

## pdc_isg_by_motif_group.py (2026-08-26)

RESULTS.md §5.2 re-cut by motif 10 / both / motif 24 state instead of by AIS/LUAD, on user request.
Everything except the grouping is copied verbatim from `perivascular_pdc.py`: the per-section
KD-tree distances, the ISG_SAFE and BLEED_CONTROL sets, `score_genes(ctrl_size=50, n_bins=25,
random_state=0)` over all 2,590 pDC at once, the six distance bands, MIN_N=50 / DRAW_MIN=10.

Section A prints the group × section table before any statistic, because `motif 10 only` pDC are
95.8% LUAD and `motif 24 only` are 68.4% AIS -- the motif split is largely a re-labelled stage
split. Section E is the decisive one: every group comparison repeated within each section. The
pooled ISG difference (p = 4e-4) does not survive it and reverses sign in P21_LUAD.

## plot_pdc_vessel_vs_tcell_by_motif.py (2026-08-26)

§9.5's figure with the pDC row split by motif 10 / both / motif 24 / neither. Computes only the
four pDC subgroup neighbourhood counts (2,590 cells, seconds); the 18 reference rows and both
random-cell denominators are read from `vasendo_neighborhood_reciprocal.json`, so the new rows are
on exactly the same scale as the old ones. Re-run that script first if it is ever recomputed.
Colours are the user's: motif 10 #0033cc blue, both #7a0ecc purple, motif 24 #cc0000 red,
neither grey, all-pDC magenta as in §9.5. Deliberately saturated so they do not read as tab20
Tumor_epi blue / Plasma red.

Section F / the per-patient figure added 2026-08-26: same script, P17 = both P17 sections, P21 =
both P21 sections. This is the cleaner cut than the per-section one -- the pooled ISG difference is
driven by a 3.5x difference in baseline interferon tone between the two patients (+0.06 vs +0.21)
crossed with the opposite patient composition of the two motif groups. Rows of the figure share a
y-axis across patients so that offset is visible rather than normalised away.

Section G / the per-section figure added 2026-08-26 (user asked for tissue-level, not patient-level).
Sample-size reality: 0 of 24 section x group x band cells reach n>=50 in either P17 section, and
`motif 10 only` never reaches 50 in any band of any section (largest bin n = 28). Only P21_AIS and
P21_LUAD have drawable curves, and only for `both` and `motif 24 only`.

## plot_pdc_vessel_vs_tcell_motif_split.py (2026-08-27)

Three separate figures instead of the single stacked one above -- `pdc_vessel_vs_tcell_ranked_`
`{motif10only,both,motif24only}` -- because the four subgroup rows in one panel could not be
compared against anything. Each figure keeps the 18 unrestricted reference rows, so those 18 and
their order are identical across the three and only the pDC point moves. Statistic, radius and
random-cell denominators are §9.5's, read from `vasendo_neighborhood_reciprocal.json`.

Two rows were added on top of what the user asked for, and they are the point of the figure:

- `pDC` (magenta) -- all pDC of the section, so the subgroup can be read against the whole.
- `any cell · <state>` (grey) -- all cell types inside the same motif region, 20,000 sampled per
  section, seed 0. **Motif states are regions.** Without this row a coloured point that moves
  between the three figures is unattributable: motif-`both` regions are 2.38x vessel-dense in
  P21_LUAD and motif-10-only P21_AIS is 0.74x, so most of the between-figure spread in panel a is
  the region, not the pDC. With it, panel a reads "pDC at or below their own region" in 6 of 7
  drawable cells and panel b reads "pDC 1.3-4.7x above their own region" in 7 of 7.

x limits are shared across the three figures so they can be flipped through. The banner is placed
from the drawn tight bbox of the axes rather than a `tight_layout` rect -- the rect and the real
title height disagree by ~8% of the figure height and leave a visible band.

Trap repeated from the by_motif script: y-tick labels are recoloured by **zipping with `order`**,
never by matching `lbl.get_text()`. The text has had `_` replaced by a space and would silently
miss or mis-hit rows.

## pdc_neighborhood_composition_by_motif.py (2026-08-27)

S8's composition analysis re-run per pDC motif state (motif 10 only / both / motif 24 only /
neither), written to answer "look at the Vas_Endo column". Two figures:
`pdc_neighborhood_composition_by_motif` (the full 19-type stacked composition, 4 sections x
5 bars, **Vas_Endo outlined in black in every bar** so the column can be found) and
`pdc_neighborhood_vasendo_by_motif` (that column alone, three panels: raw % with a dumbbell to
its region-matched reference, log2 vs whole tissue, log2 vs region-matched with stars).

Reference design is the point of the script. S8 compared the pDC neighbourhood to the whole
section; that cannot work here because **motif states are regions with different vascular
content** -- P21_AIS motif-24 regions are 16.1% Vas_Endo, P21_LUAD motif-neither regions 3.0%.
So three references are carried: `tissue_frac`, `region_frac` (all cells inside the motif region)
and `region_matched_frac` (mean 30 um neighbourhood of a random cell inside the motif region,
<=20,000 sampled). The last is the one to read. Without it the P21_LUAD `both` group looks
+0.814 enriched (vs tissue) when it is -0.108 against its own region.

Null: the S8 clump-matched permutation, but anchored in the motif region rather than the section.
200 draws, so p floors at 0.005. Anchors are region cells; expansion to the observed DBSCAN
clump size takes nearest neighbours in the whole section, so a pseudo-clump can spill out of the
region exactly as a real one can.

Runs in ~15 s -- neighbourhood counting is over <=2,590 pDC plus <=20,000 reference cells per
section x group, not over the 1.68M cells.

Result: Vas_Endo ranks 8th-17th of 19 by log2-vs-region-matched in all nine drawable groups;
pDC itself is first in all nine (+2.65 to +4.57), then cDC / T / Myeloid_other. Three of the four
significant Vas_Endo cells are depletions. Same conclusion as S11's absolute version.

## plot_pdc_vessel_vs_tcell_by_motif_pooled.py (2026-08-27)

`pdc_vessel_vs_tcell_by_motif` with the four sections collapsed to one marker per row. Plot-only --
reads `pdc_vessel_vs_tcell_by_motif.json` plus S9.5's per-section denominators, recomputes nothing,
so re-run `plot_pdc_vessel_vs_tcell_by_motif.py` first if anything upstream changes.

**Estimator, second version (the one plotted).** "One big tissue": numerator is
`sum(n_s*mean_s)/sum(n_s)` over all four sections, denominator is `sum(N_s*ref_s)/sum(N_s)`, i.e.
a random cell drawn uniformly from all 1,676,162 cells (1.292 Vas_Endo and 2.447 T within 30 um).
Neighbourhoods are still built inside a section; only the averaging is pooled. No section is
dropped for small n. The first version averaged the four per-section ratios geometrically over
sections with n >= 50 -- that is still computed and printed next to the pooled value, along with
the four per-section ratios, because the two disagree in an instructive way.

**Where they disagree and why.** The four sections differ 5x in local T density (random-cell T
within 30 um: 3.088 / 4.328 / 0.965 / 0.862) and 2x in vascular density, and the pDC motif
subgroups are heavily section-skewed (`motif 24 only` is 94% P21, `motif 10 only` 96% LUAD). So
`pDC · motif 24 only` on T cells is 1st of 23 by section-mean (4.150x) and 9th (1.457x) pooled --
same 3.57 T cells within 30 um either way, different denominator. Whenever a row's section
composition differs from the tissue's (the table is in RESULTS.md 13.3), the pooled ratio is
partly reporting which sections that row lives in.

Both numbers are in the json and the log. Say which one a figure is showing.

## spatial_maps_pdc_irf7.py (2026-08-27)

`spatial_maps_immune.py` with the pDC layer coloured by IRF7 instead of flat magenta, on request.
Four maps, one per section, `figures/luad_spatial_maps/spatial_<section>_pdc_IRF7.{png,pdf,svg}`.

Choices, all printed on the figure so a reader does not have to open the script:
- value is `X`, which in this h5ad is **log-normalised** (raw counts are `layers/counts`).
- colour scale is **shared across all four sections**, 0 to the 99th percentile of IRF7 over all
  2,590 pDC (5.11), so the maps are comparable; cells above p99 clip to the top colour.
- pDC are drawn dimmest-first so the brightest sit on top, at s=22 with a black edge, so an
  IRF7-zero pDC is still visible as a dark point rather than disappearing.
- context layers changed colour from the original map: Vas_Endo red, T light blue, **cDC brown**.
  cDC was green there and viridis runs purple -> green -> yellow, so green had to go.
- the colour bar and the per-section IRF7 histogram sit in a reserved band **below** the map.
  Putting them in an inset overlapped the scale bar and the data in two of the four sections.

Reading the gene for 2,590 rows is done with per-row CSR slice reads, not a full-matrix scan --
seconds instead of reading 427M nonzeros.

Numbers: IRF7 is detected in 54.7% of pDC overall; per section 35.9% (P17_AIS, mean 1.37) /
28.3% (P17_LUAD, 1.06) / 68.3% (P21_AIS, 2.61) / 47.7% (P21_LUAD, 1.83). The P17-vs-P21 gap is
the same patient-level interferon-tone difference recorded in RESULTS.md 10.6.

## vasendo_pdc_adjacent.py (2026-08-27)

The Vas_Endo-centred analysis. Grouping variable chosen by the user after seeing the counts:
**>= 2 pDC within 50 um** (52 / 60 / 1,000 / 847 per section). Readouts stay at 30 um.

Two control sets throughout, per the user's choice: `rest` (all other Vas_Endo of the section) and
`matched` (non-adjacent Vas_Endo matched on cells-within-50um, 20 quantile bins, <=5 per case,
seed 0). **The two disagree enough that reporting only one would have been wrong**: T/cDC
enrichment is 1.4-4.8x vs rest and 1.0-2.9x vs matched; against `rest` nearly every one of the 25
motifs looks "up", against `matched` the pattern resolves into motif 24 up / motif 10 down.

Runtime ~10 min, of which the 8 x 5,096 Mann-Whitney loops are most of it. The expression matrix is
pulled by streaming 40,000-row CSR blocks and keeping only Vas_Endo rows (29.5M nnz of 427M).

Two things to know when re-reading the outputs:
- the DE csvs have an empty-string `marker` column for non-marker genes, which **pandas reads back
  as NaN**. `fillna('')` before filtering, or every gene looks like a marker.
- `p` from Fisher on the two P17 sections is uninformative (52 and 60 cases). Those columns exist
  for direction only.

Headline: pDC-adjacent endothelium is venous (2.6-4.0x vs matched, p<1e-18, both P21 sections),
motif-24-positive / motif-10-negative (4/4 sections agree in direction), at the tumour interface
rather than the core, IL-6/NF-kB/AP-1 activated with CCL19 and CCL14 -- and **not** interferon
activated (IRF7, STAT1, TAP1, IRF1, MX1 all flat). 32 genes are consistent across the two P21
sections and none of them is a pDC/T/cDC marker, so it is not segmentation bleed.

## vasendo_pdc_adjacent_threshold_sweep.py (2026-08-27)

Robustness of s.14's "pDC-adjacent = >= 2 pDC within 50 um" over ten (radius, min count)
definitions spanning 0.01%-10.1% of Vas_Endo. Every s.14 readout is recomputed against
density-matched controls, with the matching radius following the definition's radius.

What makes it affordable: a **vectorised tie-corrected Mann-Whitney** (`mw_matrix`) that ranks the
whole 5,096-gene matrix in one pass instead of 5,096 scipy calls. Validated against
`scipy.stats.mannwhitneyu(method='asymptotic')` to 3e-6 relative on tie-heavy sparse data. The
tie term uses the identity sum_runs(t^3 - t) == sum_elements(t_i^2 - 1), which is what lets the run
lengths be computed with two accumulate passes and no per-column loop. Whole sweep: 3.6 min.

Findings that changed the s.14 write-up rather than just confirming it:
- **motif 24 is significant in all four sections** once the definition is loose enough to have
  cases (P17_LUAD +1.04*** at 30 um >= 1, n = 294). s.14 called P17 "direction only".
- **P17 is a real null for venous fraction and for expression**, not a power failure: at 50 um >= 1
  it has 489 and 838 cases and still returns 0-2 significant genes and no venous enrichment.
- **motif 10 flips sign with radius in P17_LUAD** (+0.26** at 30 um >= 1, -0.52*** at 80 um >= 2).
  That conclusion must always be quoted with its definition.
- **The gene list is less stable than the program.** A different random matched-control draw at the
  same definition changes 6 of 32 genes. Quote the 14-gene core (>=15/20 replication and present in
  both draws), not "32 genes".

`50 um >= 1` is the better default if one definition has to be carried forward: loosest definition
that still shows every effect, and the only one with >= 400 cases in all four sections.

## DE switched to the package function (2026-08-27)

`vasendo_pdc_adjacent.py` and `vasendo_pdc_adjacent_threshold_sweep.py` now call
**`alarmist.core.glm.differential_expression`** for every gene-level test. Before this they used
hand-rolled Mann-Whitney implementations (a per-gene `scipy.stats.mannwhitneyu` loop in the first,
a vectorised whole-matrix version in the second) written without checking whether the repo already
had one. It does -- from commit 95208de, "sparse tie-corrected Mann-Whitney marker DE".
**Anything new under `scripts/research/` that needs marker DE should call it, not re-implement it.**

`check_de_implementations.py` + `.log` record the three-way comparison on the s.14 reference
comparison (P21_AIS, 1,000 vs 5,000):

| | max |Δp| | Spearman | genes p<0.05 |
| --- | --- | --- | --- |
| package vs scipy per-gene loop | 3.2e-04 | 0.999995 | 870 / 870 / 870, identical set |
| package vs vectorised, continuity ON | 3.2e-04 | 0.999996 | identical set |
| package vs vectorised, continuity OFF | **9.2e-15** | 1.000000 | identical set |

So the only real difference is the **continuity correction**: the package omits it, scipy applies
it by default. Runtime: package 7.1 s, scipy loop 20.5 s, vectorised 5.2 s.

Three things that DID change with the switch, all documented in the scripts:

1. **The expressed-fraction filter sets the BH denominator.** At the package default 1e-4 a gene
   needs at least one expressing cell in *each* group. In P21_AIS vs matched that keeps 4,680 of
   5,096; in P17_AIS vs matched, where there are 52 cases, only 2,539. A smaller denominator makes
   q slightly smaller, which is why the significant-gene counts rose a little (P21_AIS matched
   214 -> 220, P21_LUAD 284 -> 293) and the cross-P21 consistent set went 32 -> 33 (a strict
   superset; ALPL is the addition, all 14 core genes unchanged).
2. **Two effect sizes now coexist.** `diff` = mean(adj) - mean(ctrl) on log-normalised X (the mean
   log fold change) is what RESULTS.md quotes and what DIFF_MIN thresholds. The package's
   `logfoldchanges` = log2(mean) - log2(mean) is carried in the csvs. On log-normalised data these
   are **not** the same quantity -- Spearman 0.72 between them on real data. Do not mix them.
3. **The hand-rolled BH is still used** for the neighbourhood-count and motif families (it is
   NaN-safe, `scipy.stats.false_discovery_control` is not). It agrees with scipy to 2.2e-16.

## A latent bug in the package DE, found and FIXED 2026-08-27

`alarmist.core.glm.mann_whitney_u_sparse_nonneg`, `src/alarmist/core/glm.py:1084`:

```python
r1 = n_zeros1 * ranks[0.0] + np.sum([ranks[v] for v in d1])
```

If a gene is detected in **every** cell of both groups, `total_zeros == 0`, so `value_counts[0.0]`
is never created and `ranks` has no `0.0` key. `ranks[0.0]` is evaluated regardless of whether
`n_zeros1` is zero, so the call dies with `KeyError: 0.0`.

**Fixed on the user's instruction**: `r1` is now built from the stored values and only adds
`n_zeros1 * ranks[0.0]` when `n_zeros1` is nonzero. `tests/test_lri.py` gained
`test_mann_whitney_u_sparse_nonneg_all_cells_detected`, which constructs a fully detected gene,
asserts no exception, and checks the p-value against `scipy.stats.mannwhitneyu(...,
use_continuity=False)` -- the package deliberately omits the continuity correction, so that is the
right reference. There is no pytest in any env on this machine; the tests were run as plain
functions (4 passed). A randomised check over 200 columns (40 fully detected, 10 all zero, 150
sparse) matches scipy to 2.2e-16, exactly 0 on the fully detected columns.

The temporary `de_guarded()` workaround has been **removed** from both analysis scripts, which now
call `differential_expression` directly. Both were rerun against the fixed package and **every
number in RESULTS.md s.14 and s.15 is unchanged** -- the workaround had dropped exactly one gene
(CCL5) in 4 of 37 DE cells, all P17_LUAD, and it was never significant.

### …and the verification that it does not matter (check_pkg_bug_impact.py, 2026-08-27)

Two questions, answered separately, log in `check_pkg_bug_impact.log`:

**Can the bug produce a wrong number rather than a crash?** No. The only silent-failure route in
`mann_whitney_u_sparse_nonneg` is `n_zeros = n - len(stored data)`, which undercounts zeros if the
matrix stores *explicit* zeros. `X/data` has **0 stored zeros in 427,196,462 nonzeros**, so it
cannot fire. Every other path either returns the tie-corrected statistic or raises `KeyError`, and
a `KeyError` aborts the run -- the sweep died and wrote no output at all rather than a wrong one.
Nothing already on disk can have been silently corrupted by it.

**Does the guard change any result?** No. All four affected comparisons were rerun against a
locally corrected copy of the function (`src/alarmist` untouched):

| comparison | dropped gene | its q if tested | significant genes guarded / fixed | max abs dq on genes tested both ways |
| --- | --- | ---: | ---: | ---: |
| P17_LUAD 15 um >=1 (77 vs 385) | CCL5 | 0.992 | 0 / 0 | 3.3e-04 |
| P17_LUAD 30 um >=2 (10 vs 50) | CCL5 | 0.974 | 3 / 3 | 8.3e-04 |
| P17_LUAD 50 um >=2 (60 vs 300) | CCL5 | 0.985 | 0 / 0 | 3.5e-04 |
| P17_LUAD 100 um >=5 (11 vs 55) | CCL5 | 0.909 | 0 / 0 | 8.1e-04 |

Symmetric difference of the significant sets is 0 genes everywhere. The residual dq comes only from
the BH denominator moving by one gene and never crosses 0.05.

## plot_vasendo_pdc_adjacent_neighbourhood_all.py (2026-08-27)

Panel a of `vasendo_pdc_adjacent_neighbourhood`, on request: **one panel** (the `rest` comparison,
not the density-matched one), **all 19 cell types**, **radius 50 um** instead of 30. New file
`vasendo_pdc_adjacent_neighbourhood_all_50um.{png,pdf,svg}` + `..._all.csv`; the original two-panel
three-type figure is untouched.

Two annotations that the extra types and the wider radius force onto the figure:

- **`rest` is not density matched.** pDC-adjacent endothelium has 51-67 neighbours within 50 um
  against 34-48 for the rest, so part of every point above 1 is local density. The density-matched
  version exists only for pDC / T / cDC (RESULTS.md s.14.1) and roughly halves those folds.
- **pDC is definitional at this radius** -- the group IS ">= 2 pDC within 50 um" -- so its 34-102x
  is a sanity check. It is drawn in a grey off-scale band at the top with its values printed,
  because on the shared log axis it compressed the other 18 types into the bottom third.

What the extra 16 types add over the three-type version:

1. The whole top of the ranking is **immune**: Myeloid_other, T, cDC, B, Langhans_cell, Macro.
2. **B cells 1.30x / 2.10x / 4.97x / 2.31x, significant in all four sections** -- not visible in
   the old figure. B + T + CCL19 + CCL14 at a venule is the textbook composition of a tertiary
   lymphoid structure, so this materially strengthens the s.14 HEV/TLS reading.
3. **SMC is at or below 1 in three of four sections (0.56 / 0.60 / 0.55)**, and Pericyte hovers at
   1. pDC-adjacent vessels are mural-cell-poor, i.e. not arterioles or large muscular vessels --
   independent support for "post-capillary venule" from a direction the venous `cell_state`
   annotation does not cover.
4. Structural types (Fibro, Pericyte, Vas_Endo) sit near 1, so this is not "more of everything".

### v2, 2026-08-27: broken y axis instead of the off-scale band

The pDC column (34-102x, definitional) is now handled with a **broken y axis** -- two stacked axes
sharing x, each on its own true log scale, diagonal break marks on the y spine -- instead of a grey
band with hollow markers. **pDC markers are drawn identically to every other point**: same section
colours, same fill, same size. Nothing about the pDC column is encoded in colour any more; the
break carries it.

### Which plotting skill governs this, and the palette check (2026-08-27)

`nature_publication_figures` is listed in CLAUDE.md but is **not present in this checkout** --
`.claude/skills/` here holds only `alarmist`, `comparator-benchmark`, `iris-compute`. Per
`VENDORED.md` it is a gitignored vendored copy sourced from `~/tansey_lab/es_xenium/.claude/skills/`
on the Mac, so it cannot be invoked from iris. What every figure in this directory does follow is
the **Plotting section of CLAUDE.md directly**: Arial, `pdf.fonttype = 42`, `svg.fonttype = 'none'`,
png + pdf + svg through the single shared saver `scripts/comparators/_common/plotting.py`.

The harness-provided `dataviz` skill was applied on 2026-08-27. Its runnable palette validator
could not execute (node here is v10, the script is an ES module), so its six checks were ported to
Python and cross-checked against the skill's own documented default palette, which passes.

**Result on the four section colours (`#4c72b0 #dd8452 #55a868 #c44e52`, seaborn deep), `--pairs
all`: FAIL.**

| pair | normal | protan | deutan | verdict |
| --- | ---: | ---: | ---: | --- |
| P17_LUAD orange / P21_AIS green | 19.7 | **4.5** | 6.5 | FAIL (below the 6.0 floor) |
| P21_AIS green / P21_LUAD red | 26.2 | 18.9 | **7.3** | WARN (6-8 floor band) |
| P17_LUAD orange / P21_LUAD red | **13.7** | 15.0 | 12.3 | normal-vision floor FAIL (needs >= 15) |

Mitigating: every figure here uses a **distinct marker shape per section** (o / s / ^ / D), which is
the secondary encoding the skill requires for the WARN band -- identity is never colour-alone. Not
mitigated: the protan 4.5 and the 13.7 normal-vision floor are hard gates in the skill's terms.

A passing four-colour alternative was found and validated: **Okabe-Ito
`#0072b2 #e69f00 #009e73 #cc79a7`** (RESULT: PASS, one WARN-band pair covered by the marker shapes).
**Not applied** -- the user asked for no colour change, and repainting would desynchronise a dozen
figures from s.9 to s.15 that all share `SECCOL`. If the palette is ever changed, change it in one
place across all of them.

### Palette swapped to Okabe-Ito across every figure in this directory (2026-08-27)

`SECCOL` is now `P17_AIS #0072b2, P17_LUAD #e69f00, P21_AIS #009e73, P21_LUAD #cc79a7` in all four
scripts that define it (`vasendo_pdc_adjacent.py`, `vasendo_pdc_adjacent_threshold_sweep.py`,
`plot_vasendo_pdc_adjacent_neighbourhood_all.py`, `pdc_neighborhood_composition_by_motif.py`).
Validator result on the new palette, `--pairs all`: **PASS**, worst normal-vision pair 18.7 (was
13.7), worst CVD pair 7.6 in the WARN band (was 4.5 FAIL) and covered by the marker shapes. Hue
identity is preserved for three of four sections; P21_LUAD moves red -> pink.

The Python port of the validator now lives here as `validate_palette.py` rather than in a session
scratchpad, and is referenced from CLAUDE.md. Cross-check it against the skill's own default
palette (must be PASS) before trusting it after any edit.

**Two palettes in this directory were deliberately NOT changed**: the 19-cell-type tab20 ramp (no
categorical palette passes at 19 slots -- the skill's own rule is to fold past 8 into "Other" or
facet, which is a design change, not a recolour) and the motif colours blue/purple/red, which the
user chose explicitly.

### The Vas_Endo self-adjacency column, checked (2026-08-27)

Asked why P17 shows no significant Vas_Endo clustering around pDC-adjacent endothelium. Premise
needed correcting and the check changed the reading of that column.

**P17_AIS is the strongest of the four** (10.94 vs 5.54 other Vas_Endo within 50 um = 1.98x,
q = 6.4e-11). The one below 1 is **P17_LUAD** (3.98 vs 5.07 = 0.79x, p = 0.042, q = 0.081 after BH
over the 19 types, i.e. not significant).

Why P17_LUAD is low: its pDC-adjacent endothelium has the **highest total neighbour count in the
table (66.8)** but the **lowest endothelial count (3.98)**, so Vas_Endo is 6.0% of that
neighbourhood against 10.7% for the rest. An isolated vessel segment inside very dense tissue --
consistent with s.14.4, where 93.3% of P17_LUAD pDC-adjacent endothelium is in `Tumor_region`.
**Not a region effect**: restricted to `Tumor_region` cells on both sides it is still 0.79x
(3.77 vs 4.75, p = 0.045).

**The important result is the threshold check.** Self-adjacency at a fixed 50 um readout, across
six grouping definitions:

| definition | P17_AIS | P17_LUAD | P21_AIS | P21_LUAD |
| --- | --- | --- | --- | --- |
| 15 um >=1 | 0.99 | 0.97 | 0.99 | 1.03 |
| 30 um >=1 | 1.05 | 0.96 | 0.99 | 1.16*** |
| 50 um >=1 | 1.11** | 0.98 | 1.01** | 1.21*** |
| **50 um >=2** | **1.98*** | **0.79*** | **0.94*** | **1.32*** |
| 80 um >=2 | 1.66*** | 1.04 | 0.97** | 1.28*** |
| 100 um >=3 | 1.34*** | 1.12* | 0.97** | 1.36*** |

At tight definitions **all four sections are ~1.0**. The 1.98 / 0.79 spread appears only at
50 um >= 2, which is exactly the row where the **grouping radius equals the readout radius** so the
two are coupled, and where P17 has only 52 and 60 cases; both extremes collapse at the neighbouring
definitions. Only P21_LUAD is consistently above 1 (1.16-1.36 over five definitions, all ***).

**Reading: pDC do not sit at vascular hubs.** The right statement is not "P17 fails to cluster" but
"none of the four clusters, except mildly P21_LUAD". That supports the s.14 interpretation rather
than weakening it -- a single post-capillary venule wrapped in immune cells raises the total
neighbour density (34-48 -> 51-67) without raising the endothelial share.

Note this column is the most demanding comparison in the figure: the reference group is also
endothelium, every one of which already lines a vessel (baseline 5.07-7.13 endothelial neighbours),
so it is really asking "is the pDC on a bigger vessel", not "is there a vessel".

### vasendo_pdc_adjacent_motifs restyled to match the neighbourhood figure (2026-08-27)

Rebuilt in the grammar of `vasendo_pdc_adjacent_neighbourhood_all_50um`: motifs on **x** with
rotated labels instead of on y, log2 on **y**, one **filled** marker per section, **BH stars above
the points** instead of filled/hollow, light vertical separators, dashed reference line, footers
under the axes.

Two deliberate departures from that figure:

- **Two rows, not one panel.** `rest` (row a) and `matched` (row b) share x and the ordering. Row a
  is kept because it is what makes the density confound visible -- against `rest` almost every one
  of the 25 motifs reads "up" -- and row b is the one to read. Dropping row a would remove the
  evidence for why row b is needed.
- **Ordered by the density-matched effect, not by motif index.** Motif number is an identifier, not
  a scale, so ranking is free and makes the pattern legible; both rows use the same order so a
  motif can be tracked between them. Motif 24 lands 3rd from the left, motif 10 five from the
  right -- the two vasculature motifs sit at opposite ends, which is the result.

The two vasculature motifs keep a shaded column (green / amber) with the label written inside row
a's axes; the x tick labels are uniform `motif N`, bold for 10 and 24. Nothing about the statistics
changed -- the figure is a restyle of the same `mo` table, and the rerun reproduced every number.

## vasendo_pdc_adjacent_motif_loadings.py (2026-08-27)

The motif figure on the continuous **loadings** instead of the binary ON/OFF state, same cells,
same two controls, same 25 motifs, same figure grammar. Grouping and matching are reproduced
bit-for-bit from `vasendo_pdc_adjacent.py` (same constants, one RNG consumed across sections in the
same order) and **asserted** against that script's published counts (52/260, 60/300, 1000/5000,
847/4235) before anything is computed -- the assert is in the script, not a comment.

**`state` is a hard threshold on `loading` for all 25 of 25 motifs** (verified: every
state-positive cell has a strictly larger loading than every state-negative one). The thresholds
are ~1e-10, i.e. the ON call is essentially "the loading is numerically nonzero". So the state
figure asks *what fraction is over the line* and this one asks *by how much* -- related, not
redundant: loading-log2 vs state-log2 across the 100 section x motif cells is Spearman **0.666**.

Two effect measures, because loadings are non-negative, heavily right-skewed and differ ~15x in
scale between motifs:

- **log2( mean loading adjacent / mean loading control )** -- plotted. Comparable in spirit to the
  state figure's log2 fraction ratio, but tail-sensitive.
- **rank-biserial r** = `2U/(n1*n2) - 1` from the same Mann-Whitney -- scale-free and
  tail-insensitive, in the csv.

**They disagree in sign in 22 of 100 matched cells (Spearman 0.669), so do not quote a single
motif from the plotted value alone** -- check `rank_biserial` in
`vasendo_pdc_adjacent_motif_loadings.csv` first. The disagreements are motifs where a handful of
high-loading cells move the mean against the bulk (e.g. P17_LUAD motif 22: log2 -3.35 but
rank-biserial +0.094). **The two vasculature motifs are not among them** -- both measures agree in
all eight of their cells.

Result on the two vasculature motifs, vs density-matched, with the state figure's value alongside:

| motif | section | log2 loading | rank-biserial | q | state log2 (s.14) |
| --- | --- | ---: | ---: | ---: | ---: |
| 24 healthy | P17_AIS | +0.836 | +0.312 | 2.4e-03 | +0.171 |
| 24 healthy | P17_LUAD | +0.522 | +0.404 | 2.0e-05 | +0.781 |
| 24 healthy | P21_AIS | -0.027 | -0.013 | 0.60 | +0.118 |
| 24 healthy | P21_LUAD | +0.890 | +0.398 | 3.5e-74 | +0.821 |
| 10 tumor | P17_AIS | -0.933 | -0.200 | 6.4e-02 | -0.334 |
| 10 tumor | P17_LUAD | -0.963 | -0.076 | 0.52 | -0.320 |
| 10 tumor | P21_AIS | -0.360 | -0.080 | 1.0e-04 | -0.404 |
| 10 tumor | P21_LUAD | -1.603 | -0.324 | 1.3e-49 | -0.435 |

**motif 10 is negative in 4/4 sections and the effect is larger on loadings than on the state**
(-1.60 vs -0.44 in P21_LUAD). **motif 24 is positive in 3/4** and, importantly, **now reaches
significance in BOTH P17 sections** (q = 2.4e-3 and 2.0e-5), which the state version could not do
at 52 and 60 cases -- a continuous outcome recovers power a binary one throws away. The one cell
that flips is **P21_AIS motif 24: state +0.118*** but loading -0.027, n.s.** There, more cells are
over the line but each carries no more loading; worth stating rather than smoothing over.

## vasendo_pdc_adjacent_motif_groups.py (2026-08-27)

The motif comparison on the **four mutually exclusive** vasculature-motif states instead of the 25
motifs tested one at a time -- the same partition used for pDC in RESULTS.md s.10-s.13:
`motif 10 only` / `both` / `motif 24 only` / `neither`. Statistic unchanged: fraction of Vas_Endo in
the state, pDC-adjacent vs control, log2 ratio, Fisher exact, BH over the four states. Grouping and
matching reproduced from `vasendo_pdc_adjacent.py` with an `assert` on its published counts.

`neither` is drawn although the question is about the other three: the four partition every
Vas_Endo cell, so hiding one makes the remaining percentages unreadable.

**Result, vs density-matched, adjacent% / control%:**

| state | P17_AIS (52) | P17_LUAD (60) | P21_AIS (1,000) | P21_LUAD (847) |
| --- | --- | --- | --- | --- |
| motif 10 only | 0.0/6.5 (no log2) | 25.0/41.3 -0.73 | 0.2/1.2 **-2.55**** | 7.6/34.5 **-2.19***** |
| both | 63.5/73.5 -0.21 | 16.7/10.7 +0.64 | 31.1/40.2 **-0.37***** | 43.1/34.0 **+0.34***** |
| motif 24 only | 36.5/15.4 **+1.25**** | 20.0/10.7 +0.91 | 63.0/46.5 **+0.44***** | 38.3/12.1 **+1.66***** |
| neither | 0.0/4.6 (no log2) | 38.3/37.3 +0.04 | 5.7/12.1 **-1.08***** | 11.1/19.5 **-0.81***** |

**This is a much cleaner statement than the 25-motif version.** `motif 24 only` is up in 4/4
sections (significant in 3), `motif 10 only` is down in 4/4 (significant in 2, and one of the other
two is 0 of 52 cells, i.e. down as far as it can go). `both` is inconsistent (-0.21 / +0.64 /
-0.37 / +0.34), so the effect is a **switch between the two exclusive states**, not a general rise
in vasculature-motif activity. The raw percentages carry it: in P21_LUAD, pDC-adjacent endothelium
is 38.3% motif-24-only against 12.1% in matched controls, and 7.6% motif-10-only against 34.5%.

**Zero-count cells have no log2.** P17_AIS has 0 of 52 in `motif 10 only` and 0 of 52 in `neither`;
with the 1e-4 pseudocount those became about -9 and flattened the whole axis in the first draft.
They are now drawn as an open down-triangle on the bottom spine labelled with the raw `0/52`, and
excluded from the y limits. Any future figure using a pseudocount log ratio should do the same.

The loadings variant (`vasendo_pdc_adjacent_motif_loadings.py`) stays on disk but the motif
comparison in the write-up uses positive fractions, per the user's instruction.

## Figures stripped to publication minimum (2026-08-27)

Every figure in this directory had accumulated methods prose: a sentence-long suptitle, panel
titles that explained the comparison, one or two footer lines, and in-plot annotations. **That
belongs in the caption, not in the figure.** Removed from all five:

`vasendo_pdc_adjacent_motifs`, `vasendo_pdc_adjacent_motif_groups`,
`vasendo_pdc_adjacent_neighbourhood_all_50um`, `vasendo_pdc_adjacent_threshold_sweep`,
`vasendo_pdc_adjacent_expression`.

What a figure here keeps: data, axes with short labels, a legend, significance stars, and a bare
panel letter via `ax.set_title(lab, loc='left')`. Where a panel needs a one- or two-word identity
(`all endothelium` / `density-matched`) it goes as small grey text inside the top-right corner, not
as a title. Nothing else.

What was deleted and now has to live in the caption: the group definition (>=2 pDC within 50 um),
the per-section n, what the two controls are, the test and the BH family, the meaning of the open
down-triangle (zero cells, no log2), and the note that panel a is density-confounded. Draft
captions for the two motif figures are in the conversation; write them into the manuscript rather
than back into the figures.

**Where the line is -- I got this wrong once and had to walk it back.** Delete the *methods
narration* (what the control is, the per-section n, the test and its BH family, "panel a is
density-confounded", "pDC is definitional"). Keep the *panel identity* -- a one- to three-word
noun phrase saying what the panel plots. Stripping `vasendo_pdc_adjacent_threshold_sweep` to bare
`a`-`h` made panels c-g unreadable, because their y labels are all `adjacent / matched` or
`log2 vs matched`; nothing on the figure said which was T cells, which was motif 24, which was
venous. They are back as `a cases`, `b share of endothelium`, `c T cells`, `d cDC`, `e motif 24`,
`f motif 10`, `g venous`, `h agreement with 50 um >=2`.

Two things now exist ONLY in the caption and will mislead if the caption is not written: the four
points above the axis break in `..._neighbourhood_all_50um` (pDC, 34-102x, definitional), and the
coloured points in `..._expression` (pDC/T/cDC markers, i.e. transcript bleeding).

Note for the next figure written here: `set_title(lab, loc='left')` for the panel letter, not
`ax.text(0, 1, ...)` in axes coordinates -- the latter collides with the topmost y tick label.

## vasendo_pdc_adjacent_venous.py (2026-08-27)

The Q4 venous result of `vasendo_pdc_adjacent.py` as a figure. That script computed and printed it
but drew nothing, so it was the only one of the four s.14 conclusions with no figure -- it lived
only in the log, the json `context` field, and RESULTS.md s.14.4. Same grammar as
`vasendo_pdc_adjacent_motif_groups`: two rows for the two controls, log2 on y, section markers,
BH stars, no prose. x is the four sections (a single binary feature has nothing else to put there),
so the tick labels carry section identity and there is no legend.

**One correction to s.14 made here: BH over the four sections, within each control.** The original
stored raw Fisher p only, inconsistent with every other family in that script. It changes nothing --
the two significant cells go from p = 1.3e-18 / 4.1e-25 to q = 2.6e-18 / 1.6e-24 -- but the table
should not have carried bare p.

| section | adjacent venous % | matched % | rest % | log2 vs matched | q |
| --- | ---: | ---: | ---: | ---: | ---: |
| P17_AIS | 7.69 (4/52) | 10.00 | 4.44 | -0.38 | 0.80 |
| P17_LUAD | 15.00 (9/60) | 21.33 | 12.89 | -0.51 | 0.40 |
| P21_AIS | 8.00 (80/1,000) | 2.04 | 0.96 | **+1.97** | **2.6e-18** |
| P21_LUAD | 20.78 (176/847) | 7.96 | 7.97 | **+1.38** | **1.6e-24** |

Note the shape the figure makes visible and the RESULTS.md table did not: **the two P17 sections
are slightly NEGATIVE against matched controls** (-0.38, -0.51), not merely non-significant. With
4 and 9 venous cells they carry no weight, but "no effect in P17" is the honest phrasing, not
"weaker effect".

Provenance guards in the script, because this readout depends entirely on someone else's
annotation: it asserts that Vas_Endo splits into exactly `Stromal_C0_Vas_endo_capillary` and
`Stromal_C7_Vas_endo_venous`, and that the substring `venous` matches exactly one of the 56
`cell_state` categories. `cell_state` came from `Jiayi_sample.8.21.2026.rds` via
`add_cell_state_to_adata.py` and is untouched by ALARMIST, which is what makes this
model-independent.

## §16 — pooled expression DE, and marker exclusion done the GBM.ipynb way (2026-09-01)

Two requests, one script plus one edit.

**What §14 Q3 actually was** (checked before changing anything): per section, four
`differential_expression` calls per control set, `alarmist.core.glm.differential_expression` with
its default `min_in/out_group_fraction=1e-4`. Not scanpy, not pooled. The volcano plotted mean
log-norm difference vs −log₁₀ p and *coloured* a 31-gene hand list of pDC/T/cDC markers.

**Marker exclusion.** Now `tutorials/GBM.ipynb`'s route, not a hand list:
`al.load_exclusion_mask('/data1/tanseyw/projects/fanj2/alarmist_luad/markers/exclusion_matrix.csv')`
→ `glm._filter_genes_for_volcano`'s two rules (expressed fraction ≥ 0.02 within the focal cell
type; drop `exclusion_mask[other].any(0)`). Focal type = Vas_Endo. 3,425/5,096 genes are another
type's marker; 1,021 survive both filters.

- **The mask was NOT recomputed.** `compute_exclusion_mask` subsamples with an **unseeded**
  `np.random.choice`, so recomputing gives a different mask than the one this dataset's own run
  saved in Feb 2026 and that the published impact figures used. Notebook cell 86 is exactly
  `load_exclusion_mask` on that file — reuse is the documented path, not a shortcut.
- There are three files in that directory: `exclusion_matrix.csv` (Feb 2026, used),
  `exclusion_matrix_question.csv` (Jan 29) and `exclusion_matrix_wrongde.csv` (Jan 19). Only the
  first is canonical; the third differs (e.g. A2ML1 is flagged a Tumor_epi marker there and not in
  the canonical one). **Do not grab one by tab-completion.**
- The mask catches all 31 hand-picked genes and 3,394 more. Two genes that §14 was reading as
  biology are now gone: **CCL19** (T marker) and **APP** (pDC marker), both in §14's old top-25.
- Filtering happens **after** testing, as `glm_volcano` does, so BH stays over all tested genes.
  Both counts are printed everywhere ("N all genes -> M after the filter").

**Volcano style** is now the package's own `alarmist.plotting.glm_plots.volcano_plot`, called with
the notebook's `fdr_threshold=0.05` / `lfc_threshold=0.2` and
`_create_volcano_figure_for_motif`'s `n_top=30`, `fontsize=5`, 4×4-inch panels. Consequences worth
knowing: **x becomes log₂FC and y becomes −log₁₀(q)** (§14 plotted mean log-norm difference and
−log₁₀ p, so panel positions are not comparable to the old figure); `volcano_plot` calls
`sns.set(...)`, which **clobbers rcParams globally**, so `apply_publication_style()` has to be
re-called before `save_all_formats` or the PDF/SVG lose `fonttype 42` / `svg.fonttype none`; and
its axis labels are hard-coded at `fontsize + 14` = 19 pt, which is proportionate only at the
notebook's 4-inch panel size. The `>= 300` jitter it inherits is stochastic — `np.random.seed(0)`
before plotting.

**Pooling.** Grouping stays per section (coordinates don't cross sections) and the script
**asserts** §14's published 52/60/1,000/847 and 260/300/5,000/4,235 before doing anything. Only
labels are pooled. Cases are 94% P21; `matched` is 94% P21 too (5× per section, by construction),
`rest` is 57% — so the rest panel is section-confounded and the matched panel is not.

**The answer:** 43 pooled hits vs matched after filtering (36 up), and they are an NF-κB/AP-1
activation program (SOCS3, JUNB, JUN, NFKB1, NFKB2, RIPK1, NR4A1, IL6, **SELE**, CDKN1A, HIF1A,
PIM1, PER1, MYC). SELE is the post-capillary-venule adhesion molecule, i.e. the expression-side
echo of §14 Q4's venous result. But Spearman(pooled, per-section) is +0.71/+0.70 for the two P21
sections and +0.07/+0.11 for the two P17 ones, and P17 gives **0** significant genes on its own:
the pooled number is a P21 number. 37/43 also pass a Stouffer meta that never pools cells, so
pooling is not manufacturing the signal — it is just letting P21 speak for everyone.

**Trap avoided:** the q values in the first draft of RESULTS §16.4 were transcribed by eye from a
log that printed `+0.000`, and five of them were wrong by orders of magnitude. The table is now
generated from `vasendo_pdc_adjacent_expression_pooled.json` by a script. Never hand-copy a number
a float-formatted log has already rounded away.

**Both control sets are now saved as separate figures (2026-09-01).**
`vasendo_pdc_adjacent.py`'s fig 3 loops over `('matched', 'rest')` and writes
`vasendo_pdc_adjacent_expression_matched.*` and `vasendo_pdc_adjacent_expression_rest.*`. The
old single file `vasendo_pdc_adjacent_expression.*` was the *matched* one and nothing in its name
said so — the rename exists so nobody reads an uncontrolled volcano as a controlled one. The old
files are left on disk (identical content to `_matched`) rather than deleted; delete them when
convenient. The two figures are meant to be read as a pair: their difference is the density
effect, which in P17 is the whole result (11 and 5 hits vs rest, 0 and 0 vs matched).

**`adjust_text` is not reproducible run to run.** Re-running `vasendo_pdc_adjacent.py` with the
same seed reproduced the P17_AIS and P17_LUAD volcano panels **byte for byte** but moved label
positions in P21_AIS and P21_LUAD (3.6% and 6.5% of pixels; all differing pixels are label ink,
none are data points, and every printed statistic is identical). `np.random.seed(0)` does not fix
it — the label solver's placement drifts on its own. So: a volcano PNG is not a reproducibility
check, and if a byte-stable figure is ever needed, lower `N_TOP` until labels stop colliding.

**Volcano relabelling, and why `_volcano.py` exists (2026-09-01).** `glm_plots.volcano_plot` picks
its own top-n labels and runs `adjust_text` on them before returning, so any gene that must always
be labelled can only be added afterwards, in a second pass that knows nothing about the first
pass's boxes — they collide. `_volcano.py::volcano_panel` copies the package's grammar
(colours, guides, label boxes, `fontsize+14` axis labels, limit rules) and does **one**
`adjust_text` over top-n *plus* forced labels. Both `vasendo_pdc_adjacent.py` and
`vasendo_pdc_adjacent_expression_pooled.py` import it — one implementation, not two copies.

Two things learned drawing it:

1. **A single cell can set a panel's axis.** P17_AIS CXCL11 is expressed in **1** of 52
   pDC-adjacent cells, giving log₂FC +3.69 at q = 0.025; bold-labelling it widened that panel to
   ±4 and flattened everything real. Hence `MIN_DET_CELLS = 5` on forced labels, with every skip
   printed. Generalise the habit: anything force-drawn should carry a detection floor.
2. **The chemokines are not in the plotted set.** They are all other cell types' markers, so §16's
   filter removes them; drawing them means overlaying genes the figure otherwise excludes. Open
   black diamonds + a one-line legend, never as ordinary points — otherwise a reader assumes they
   passed the same filter as everything else.

**Spatial maps of the §14 case set (2026-09-01).** `spatial_maps_pdc_adjacent_vasendo.py`, colours
inherited verbatim from `spatial_maps_immune.py`. Three things that had to be got right:

1. **Draw order is the whole figure.** pDC must be drawn *after* the orange pDC-adjacent discs.
   In the first version pDC were at zorder 5 and the discs at 10, so a pDC sitting under a case
   cell vanished — inset i of P17_AIS looked like one pDC when the log says two, i.e. the figure
   silently contradicted the ≥2 rule it was meant to illustrate.
2. **Insets are mandatory, not decoration.** 52 cells in a 9 × 10 mm section are invisible; without
   the 250 µm zooms the figure shows nothing. Centres are picked deterministically (densest
   clusters ≥800 µm apart), never by eye.
3. `'ivx'[j]` gives **i, v, x** — roman 1, 5, 10, not 1, 2, 3. Shipped that in the first draft.
   Use an explicit `['i', 'ii', 'iii']`.

Also: I again waited on `pgrep -f <script name>`, and again the waiting shell's own command line
matched, so the loop ran the full 10-minute timeout while the job had finished in ~2. CLAUDE.md
documents this exact trap. Capture the PID (`cmd & PID=$!`) and wait on `kill -0 "$PID"` — which
is what the later reruns in this session did, correctly.

## spatial_explorer.py — interactive HTML explorer for all four sections (2026-09-02)

Builds ONE offline HTML holding **all 1,676,162 cells** of P17/P21 AIS+LUAD at full resolution,
recolourable by cell type (tab20, 19), cell state (56, checkbox-subset), motif loading
(continuous, one motif per panel), motif state, `tissue_layer`, and the curated H&E blood-vessel
lumens. Modelled on `/data1/tanseyw/projects/fanj2/alarmist_luad/08_spatial_explorer.py` (a
different dataset — ES Xenium TMA punches, 20 motifs), but the rendering had to be rewritten.

```bash
export PYTHONNOUSERSITE=1
/data1/tanseyw/fanj2/envs/comp-liana/bin/python \
  scripts/research/luad_revision/spatial_explorer.py \
  --out /data1/tanseyw/projects/fanj2/alarmist_luad/spatial_explorer_AIS_LUAD.html
# ~30 s, <1 GB RSS, 21.0 MB out.  --motifs 10,24 gives a ~6 MB page instead.
```

**Output goes to `/data1`, deliberately.** `*.html` is **not** in this repo's `.gitignore`
(`git check-ignore foo.html` → rc 1), so a 21 MB file written anywhere in the tree except
`figures/` would be committable. Note also that on this iris checkout `scripts/research/` is
*not* locally excluded, contrary to CLAUDE.md:441 — `.git/info/exclude` is empty here.

### The four decisions that were measured rather than guessed

1. **`fillRect` per cell does not scale.** The reference draws one `fillRect` with a per-point
   `fillStyle` string; its worst panel is 27k cells, ours is 640,739 — 23.7×. Measured at
   **81 ms per panel with the canvas call stubbed out to a no-op**. Replaced by writing into a
   `Uint32Array` aliased onto an `ImageData` buffer (**5.9–15 ms** for the same 640k at
   1800×1800) plus a `Uint8` per-pixel depth buffer. The depth buffer is both faster than a
   16-bucket counting sort (3.2 vs 5.7 ms) *and* exactly correct, where the reference's
   two-pass `v<8` / `v>=8` is only approximately ordered. Endianness is detected, not assumed.
2. **Linear percentile scaling of the loadings is unusable.** Positive loadings span ~4 orders
   of magnitude, so linear p1–p99 puts a **mean of 55 % and a worst case of 81 % (motif 2) of
   the positive cells into colour level 1** — a flat wash with a few bright dots. **Log,
   anchored at p5/p99 of the positive cells**, gives a modal bin of 11.9 % (worst 15.4 %), no
   bin under 0.8 %, and 3.72 of 3.91 available bits. p1 is the wrong anchor: for motifs
   3,4,7,8,9,14,19,20,21 the 1st percentile of positives sits on the ~1e-10 numerical floor.
   Cost: the motif planes grow from 11.5 MB to 16.2 MB of base64. Worth it.
3. **Plain `magma` is the wrong direction on a white page.** Its high end `#fcfdbf` has
   **1.05:1 contrast against white** — the strongest-loading cells become invisible — while its
   near-black low end makes the 35–65 % weakest cells the most salient thing in the panel. The
   ramp is `magma_r` truncated to [0.20, 0.90], with level 0 (motif-negative) held **out** of
   the ramp as `#e8e8e8`, so "no signal" and "weakest signal" are never the same pixel.
4. **Motif state costs zero bytes.** Every motif's `positive` state is exactly a threshold on
   its loading, so the 4-bit level is built as 0 for negatives and 1..15 for positives and
   `level > 0 ⇔ positive` is recovered exactly. The builder asserts this per motif against the
   stored `motif_k_state`, which matters because the underlying threshold is razor-thin
   (motif 0: max negative 1.899527e-10 vs min positive 1.899533e-10 — a gap in the 7th
   significant digit). Saves 3.7 MB.

### Traps that bit, or nearly did

- **A real deadlock, caught only by running the JS.** `pump()` stops scheduling animation
  frames once every remaining job is waiting on a motif plane — but the job only *asked* for
  its plane when it was dequeued, so the last motif panels were never dequeued, never
  requested, and stayed blank forever. `needMotif()` is now called in `enqueue()`. Syntax
  checks and visual PNG renders both missed this; a node harness that stubs
  DOM/Canvas/Blob/Response/DecompressionStream and drives the real code found it immediately.
  That harness is the reason to write one: `node harness.js <the html>`.
- **`<SEC>_lumens.parquet` and `<SEC>_vessels.parquet` do not share a row order.** P17_AIS
  disagrees (lumens row 1 is `P17_AIS__30`); the other three happen to line up, so a row-order
  join passes 3 of 4 spot checks and silently corrupts P17_AIS's bbox fallbacks. Join on the
  integer suffix of `vessel_id` and assert against `cx_um`.
- **`poly_world_json` can be the literal string `'null'`.** `json.loads` returns `None` without
  raising, and `pd.isna()` does not catch it. Two rows (both P21_AIS) need the rotated-bbox
  fallback.
- **`plt.get_cmap('tab20', 19)` is not the first 19 tab20 colours** — it resamples and drops
  `#c5b0d5`, shifting 10 of 19 assignments. Both were run through `validate_palette.py
  --pairs all`: neither can pass (19 levels vs a normal-vision floor of 15), but the first-19
  assignment has worst pair 5.5 and **no FAIL**, while the resampled one FAILs on
  `#c7c7c7 / #9edae5`. The first-19 version also separates the two types this dataset most
  needs separated, Tumor_epi `#c7c7c7` from Vas_Endo `#bcbd22`, where the resampled one
  collides them at `#bcbd22 / #dbdb8d`. Kept the first 19.
- **"Rare categories on top" is off by default for cell type.** It is a genuine aid for finding
  a 2,590-cell population, but at 550 px roughly one cell owns each device pixel and drawing
  the rarest one on top turns P17_LUAD into confetti that hides the tumour mass. It stays on
  permanently for `tissue_layer` only, where the levels are contiguous regions so the ordering
  decides nothing but the boundary pixel — which is exactly where the thin interface bands are.
- **56 cell states cannot be told apart by colour.** The measured ΔE ceiling for *any* 56
  colours on a light page is 8.77 against a required floor of 15; the shipped lineage palette
  reaches 5.51. The checkbox subsetting is the mitigation and is load-bearing, not a
  convenience — which is why the layer was specified that way in the first place.

### Verification actually performed

- Numeric round-trip: the emitted HTML is re-decoded and the kernel re-implemented in numpy
  (`verify_render.py`), asserting per section that the range is contiguous and `qy`-sorted with
  the declared extent, and per motif that `level>0` count equals the stored `npos` and that no
  colour level is empty. Panels rendered to PNG and inspected.
- Vessel frame: the curated lumens were **never documented as sharing `obsm['spatial']`'s
  frame**, so it was tested rather than assumed — cell density inside the lumens is 0.03–0.22×
  the section mean and the vascular-cell fraction within 40 µm of a lumen centroid is 0.32–0.73
  against 0.11–0.26 at random points, in all four sections. A zoom crop shows the endothelium
  lying exactly on the traced rim. **This belongs in RESULTS.md and is not yet there.**
- End-to-end: `node harness.js` boots the real JS, decodes all 1.68M cells, builds the
  controls, renders 24 panels across 6 views, drives the IntersectionObserver and the lazy
  per-motif decode, and opens the zoom modal.

## motif_vessel_alignment.py (2026-09-02)

"Does any motif's spatial structure align with the curated blood-vessel annotation?" Signed
distance from every cell to the nearest curated H&E lumen **boundary** (negative inside),
then per section × motif: enrichment of motif-positive cells at 0–25 µm, AUROC of the
continuous loading for perivascular (≤25 µm) vs far (>100 µm), the radial profile in seven
distance bands, and the fraction of the 309 lumens whose 0–25 µm ring is enriched.

**Answer: motif 1**, in all four sections. It is the only motif enriched in **9/9**
section × tissue-layer strata (median 2.33×); motif 23 manages 5/9 and everything else less.

Three things this analysis had to get right, and one thing it revealed:

1. **Endothelium is trivially at vessels**, so every statistic is also computed on
   non-vascular cells only (dropping Vas_Endo, Lym_Endo, SMC, Pericyte — 16.5% of cells).
   Motif 1 keeps AUROC 0.695–0.763 there, so the alignment is not just "SMC sit on vessels".
2. **Tissue layer is the real confound** — vessels are not uniformly distributed across
   tumour/normal/interface. Stratifying by `tissue_layer` is what makes the result safe.
3. **The annotation is 71–83 large lumens per section**, so only 3.2–5.2 k of 180–640 k cells
   are within 25 µm of one. Motif 1's *specificity* is therefore tiny (0.4–3.2% of its
   positives are perivascular). The claim is "enriched at curated large vessels", never
   "only at vessels", and the two must not be confused.
4. **The two motifs the paper calls vasculature do not align with the curated lumens, and
   that is not a contradiction.** Motif 10 (tumour vasculature) and motif 24 (healthy
   vasculature) have radial profiles that *rise with distance* (peak/far 0.12–0.55× and
   0.44–0.64×) and are enriched in 1/9 strata each. Their top LRIs say why: motif 10 is
   `VEGFA→KDR Tumor_epi→Vas_Endo` plus pericyte→tumour laminins — the tumour–neovessel
   interface, i.e. capillaries inside tumour, which the curator did not outline. Motif 24 is
   `COL4A3→SDC4` / `VEGFA→KDR Alveolar_epi→Vas_Endo` — the alveolar–capillary unit. Neither
   is a muscular vessel wall. Motif 1 is: `JAG1→NOTCH3 SMC→SMC`, `GJA5→GJA5 Vas_Endo→Vas_Endo`
   (connexin-40, arterial endothelium), `ANGPT1→TEK SMC→Vas_Endo`, `EDN1→EDNRA/B
   Vas_Endo→SMC`, `PDGFA→PDGFRB` — textbook endothelium–smooth-muscle crosstalk, and SMC is
   its top enriched cell type (+1.51 log2). So the geometry and the biology agree
   independently, which is the reason to believe it.

Figures: `mv_<section>.png` in the session scratchpad (700 µm crops of the three largest
lumens per section, motif 1 / 23 / 10 / 24 loading + cell type, curated outline in cyan).
Numbers: `/data1/tanseyw/projects/fanj2/alarmist_luad/motif_vessel_alignment.json`.
**Not yet written into RESULTS.md.**

## Adversarial review of `spatial_explorer.py` (2026-09-02)

Six independent reviewers over the builder and its embedded JS, each finding adversarially
verified. Only one finding completed verification before the session limit killed 26 of the
52 agents, so the rest were triaged by hand against the code and the data. Everything below
was checked directly before it was acted on; nothing was applied on an agent's word alone.

**Two real bugs the earlier verification had missed, and why it missed them.**

1. **All 11 per-family "all / none" links in the cell-state sidebar were dead.**
   `buildControls()` uses one function-scoped `var h` as a moving cursor over the sidebar
   hosts and reassigns it to `#f_opt` at the end; the per-family click handlers close over
   it, so by the time one fires `h.querySelector('input[data-v=…]')` returns null and
   `boxes.every(b => b.checked)` throws. The node harness had exercised only the three
   *global* all/none links, which read their host from `data-all`. Fixed by capturing
   `hState` before the reassignment. The harness now clicks all 11 (2 → 56 states selected,
   0 throws) — a stub bug hid this too: `el().className` was a plain property, so
   `classList.contains('fam')` never matched and the test silently found 0 headers.
2. **The zoom modal's "N cells in view" was the y-band width, not the canvas count.**
   `paint()` returned `i1 - i0` from the `lowerBound` search, which restricts `qy` only and
   therefore spans the full width of the section; `renderSlice` culls x per pixel. At the
   opening zoom the two agree (640,739 vs 640,734), which is why it looked right; at 16×
   it overstated by 25–70×. `renderSlice` now counts what it actually stamps. Verified:
   full view 640,734, 16× zoom 2,094.

**Also fixed, all verified against the data first:**

- A rejected motif-plane decode left the rejected promise in `PENDING`, so every later
  `needMotif` returned it, `pump` was never rescheduled and the page stuck on
  "rendering N / M" forever. Now caught, the motif is marked failed, its jobs are dropped
  and the failure is surfaced.
- The `k >= 3` stamp branch rejected any cell whose disc crossed the panel edge instead of
  clipping it, blanking a `k`-pixel band. Invisible at `k=3` in the grid, but `k` reaches 14
  in the zoom modal. Interior cells keep the fast path; border cells take a clipped one.
- The vessel view hard-coded the context grey, so "Faint grey context cells" did nothing
  there; the motif-loading legend kept its grey "motif negative" swatch even when negatives
  were painted white; the motif-state legend showed a meaningless two-stop gradient when
  "coloured by cell type" was on, where the panel actually carries all 19 tab20 colours.
- The status line reported all 1,676,162 cells regardless of which sections were ticked.
- Vessel vertex counts were accumulated **after** Douglas-Peucker, so the log claimed
  simplification removed nothing (`995 -> 995`). It removes 15.8% (4,482 → 3,773).
- `assert (vv > 0) == pos` is true by construction and proved nothing. Replaced by the
  property actually claimed: `min(positive loading) > max(negative loading)`, per motif.

**Five statements on the page were wrong, and are corrected.** Each was re-measured:

| said | measured |
|---|---|
| loadings span "about four orders of magnitude" | **8.1–9.3 decades** over the positive cells; p5–p99 is the ~4-decade *window shown*, and it clips |
| state/loading "gaps of ~1e-16" | **4.4e-16 to 9.3e-15** |
| "309 hand-curated" outlines | 309 curated, of which **246 hand-edited, 63 seed-accepted unmodified**, 5 flagged, 2 bbox fallbacks |
| "one micron per pixel everywhere" | one *common* µm-per-pixel, **~21 µm per CSS pixel** at the 550 px default |
| tab20 first-19 has "no outright FAIL" | **both** assignments return RESULT: FAIL with 33 failing rows; the tie is broken on Tumor_epi/Vas_Endo separation, not on the score |

**And one that had gone stale within the session:** the page said motifs 10 and 24 were the
only characterised ones and the other 23 unknown. Motif 1 is now characterised (above), so
it is labelled, and a caveat states what the vessel analysis found — that the two motifs
*named* "vasculature" are named for their LRI content and are not the ones sitting on the
curated lumens.

**Deliberately not changed**, with the reason: the zoom modal ignores "same micron scale"
(fitting each section is the right behaviour for a detail view); the `k=2` stamp is anchored
top-left rather than centred (a half-device-pixel offset against the vector overlays); a
panel that scrolls out of view before its queued job runs is repainted rather than
cancelled; and the same panel can be enqueued twice, inflating `TOTAL`. All are cosmetic or
performance-only and none changes a number a reader sees.

## plot_vessel_strata.py (2026-09-02)

The geometry figure for the vessel-motif analysis: what "near a blood vessel" actually
means. Same crop style as the `mv_<section>.png` panels, but the strata are drawn rather
than left implicit. Row a colours every cell by the stratum the test assigns it — inside
the lumen, RING 0–25 µm, the unused 25–100 µm gap, SURROUND 100–250 µm, beyond 250 µm, and
cells thrown out of the surround for being within 25 µm of a *different* curated lumen —
with the 25 / 100 / 250 µm contours drawn from an exact distance field (`contourf`, so
concave lumens are handled). Row b is motif 1 loading on the identical crops.

Three choices that make the figure say something rather than decorate:

1. **Fill the bands, do not outline them.** A 25 µm ring against a ~1 mm field is 2% of
   the width; the first draft drew it as a line and it read as a line, not as the area it
   stands for. Translucent `contourf` fills plus a thin contour is what makes it legible.
2. **Show a large, a median and a small lumen, not the three largest.** The whole point of
   measuring distance to the lumen *boundary* rather than the centroid is that one 25 µm
   ring means the same thing for a 315,219 µm² artery and a 147 µm² venule. Three big
   vessels cannot show that; one of each does.
3. **Exclude the two rows with no traced lumen from the size ranking.** Their
   `lumen_area_um2` is None, which sorted to 0.0 and made a rotated-bbox fallback the
   "smallest vessel" example in P21_AIS — a misleading illustration of a real vessel.

Reading it: the ring is thin and the surround is large, which is why the per-vessel counts
are lopsided (P21_LUAD__0: 40 ring cells vs 1,140 surround) and why no single vessel
reaches significance. P21_LUAD panel 1 is the clearest example of the exclusion rule doing
work — a cluster of orange cells sits in the surround band but next to a neighbouring lumen.

`--section {P17_AIS,P17_LUAD,P21_AIS,P21_LUAD} --motif K --n-vessels N`. Rendered without
Arial (not installed on the compute nodes) — re-render before manuscript use.

## Ring radius widened to 80 µm (2026-09-02)

25 µm is 2–3 cell layers and was too tight, so the ring/surround radii are now CLI
parameters on all three scripts (`--ring / --sur-lo / --sur-hi`, `--near / --far` for
`motif_vessel_alignment.py`) and every output carries the radius in its filename. The
25 µm results are kept, renamed `*_ring25*`; the current run is `*_ring80*`.

**The surround has to move with the ring.** It must start beyond where the perivascular
signal has decayed, or the contrast is diluted and the effect is *under*-stated. The
measured motif-1 radial profile is back at section baseline by 100–200 µm, so a 25 µm ring
pairs with 100–250 µm and an 80 µm ring needs 200–400 µm.

**What widening buys** (`vessel_motif_ring_sweep.csv`, surround pinned at 200–400 µm so
only the ring changes):

| ring µm | vessels scorable | median ring cells | m1 | m23 | m10 |
|---:|---:|---:|---:|---:|---:|
| 10 | 115 | 30 | +0.58 | +0.05 | −1.96 |
| 25 | 283 | 45 | +0.44 | +0.27 | −1.54 |
| 50 | 308 | 102 | +0.52 | +0.26 | −1.23 |
| **80** | **308** | **189** | **+0.48** | **+0.14** | **−0.89** |
| 120 | 308 | 339 | +0.37 | +0.08 | −0.50 |
| 160 | 308 | 529 | +0.32 | +0.06 | −0.34 |

**Motif 1 is rank 1 at every radius from 10 to 160 µm.** The conclusion does not depend on
where the ring is drawn, which is the point of running the sweep rather than asserting a
choice. 80 µm buys full vessel coverage (308 of 309, against 283 at 25 µm) and 4× the cells
per ring, at the cost of ~8% of motif 1's effect size.

**Widening also separated two signals that 25 µm had conflated.** Motif 1 decays slowly
(+0.58 → +0.32 from 10 to 160 µm) — a broad perivascular program. Motif 23 peaks at 25–50 µm
(+0.27) and collapses by 120 µm (+0.06), dropping from rank 2 to rank 5: a *tight*
endothelial-adjacent signal. That difference in spatial scale is real and was invisible at a
single radius.

At 80 µm, five motifs are enriched at q < 0.05 (Wilcoxon signed-rank of the per-vessel
scores, BH over 25): m1 +0.48, m23 +0.14, m16 +0.13, m17 +0.12, m18 +0.08 — motif 1 is 3.3×
the next. Per-vessel significance improves but stays low: **13 of 308** vessels, up from
8 of 283. Motif 10 remains the most depleted (−0.89, q 5e-28).

Two things the figure caught that trusting the code would not have: the panel-a title still
said "only motifs 1 and 23 are enriched", which is false at 80 µm, so it is now derived from
the data; and `plot_vessel_strata.py` had four legend/title strings hard-coded to 25/100/250
while the geometry was correctly 80/200/400. **Look at the output.**

## Why the surround, and what the background reference gives (2026-09-02)

Asked directly: why compare the ring to an outer annulus rather than to everything else on
average? Two separate things had been conflated, and the honest answer is that only one of
them was a reason.

**Not a reason.** The v1 section-wide comparison collapsed under proper null calibration
(null SD 4.1–4.7, zero surviving calls). That was a failure of the *inference*, not of the
*estimand* — cells in a ring are spatially autocorrelated whatever you compare them with.
Calibrating against spatially matched random rings fixes it for either reference. Preferring
the surround on those grounds was wrong.

**A reason, but a trade.** The surround removes regional confounding: vessels sit
non-randomly (81 of 84 tumour-region vessels in one cluster), so a section-wide reference
partly measures where vessels are rather than what is around them. The price is that it is
conservative — any perivascular signal reaching into the annulus is differenced away — and
that it changes the question from "is this enriched at vessels?" to "is this enriched at
vessels beyond its neighbourhood?". The composition adjustment on
(section × cell type × tissue layer) already absorbs much of the regional effect, so the
background reference is better controlled than it first looks.

`--reference {surround,background}` now selects it; `background` = every cell in the section
more than `--ring` µm from **any** curated lumen (96.3–99.4% of cells), written to
`*_bg.*`.

**At 30 µm the two references agree, Spearman ρ = +0.93 over the 25 motif medians:**

| | background | own surround |
|---|---:|---:|
| m1 SMC / NOTCH3 vessel wall | **+0.481** | **+0.493** |
| m23 endothelial identity | **+0.327** | **+0.296** |
| m17 plasma / BAFF-APRIL | +0.132 | +0.046 |
| m3 / m12 / m19 (most depleted) | −3.10 / −3.05 / −3.02 | −2.00 / −1.59 / −1.74 |
| enriched at q<0.05 | m1, m23, m17, m18 | m1, m23 |
| vessels significant on their own | 0 of 300 | 8 of 299 |

Two systematic differences, both expected:

1. **Depletions are ~1.5× larger against the background**, because the background includes
   tumour cores and airways where those motifs are high, whereas a vessel's own annulus is
   more similar to its ring. The enrichments barely move, which is the reassuring part.
2. **The surround is slightly more powerful per vessel** (8 vs 0), because differencing the
   local region removes variance from the observed score and the null in a matched way.
   Neither is anywhere near enough for per-vessel calls.

So the headline is reference-independent: **motifs 1 and 23 are the only ones enriched
around annotated vessels, and no individual vessel can be called.** Both are shipped
(`vessel_motif_local_ring30.*` and `..._ring30_bg.*`) and panel e of each figure is the
agreement scatter.

## vessel_motif_cooccurrence.py — do vasculature motifs co-occur? (2026-09-02)

New question, not answerable from the earlier work: `vessel_motif_local.py` asked motif by
motif whether each is *enriched* near vessels; this asks whether, among perivascular cells,
the motifs are on *together*. Cells within 30 µm of a curated lumen (20,263 of 1,676,162,
0.60–3.55% per section), six motifs the user named (0, 1, 6, 10, 23, 24), full 25×25 also
computed.

**The user's motif list checks out on a loosened criterion.** By share of a motif's total
LRI score involving Vas_Endo as sender or receiver, the six rank 1, 2, 3, 4, 5 and 9 of 25:
m23 85.5%, m10 56.5%, m24 28.8%, m1 26.2%, m0 14.4%, m6 9.8%. My first reading — that m0
(T/T) and m6 (Macro/Macro) "are not vasculature motifs" — was based only on their top-6
LRIs, i.e. the dominant signal, not the whole motif. On the %Vas_Endo criterion m17 (12.1),
m16 (11.3) and m22 (10.0) would also qualify, above m6.

**Method.** Loadings are zero-inflated and span 8–9 decades, so each is converted once to a
within-section normal score (rank → probit) and every correlation is a Pearson on those —
computed globally, not re-ranked inside each subset, so perivascular and background numbers
stay comparable. Everything is reported raw and as a partial correlation with the
(cell_type × tissue_layer) stratum mean removed, because perivascular cells are
SMC/endothelium-rich and two motifs could correlate merely through shared composition.

**The trap, and it nearly produced the wrong headline.** Against matched pseudo-rings the
raw reading is that m0/m6/m10/m23/m24 all correlate *more* near vessels (z = +2 to +5).
That is almost entirely a **global shift**: over all 300 motif pairs the median z is
**+3.6** and **60% of pairs exceed +3**. Perivascular cells are globally more coherent than
matched tissue patches. A pair is therefore only notable relative to the other pairs, and on
that basis **not one of the positive focus pairs reaches even the top decile** — the best,
m6–m10, sits at percentile 65.

**The real result is anti-co-occurrence, and it is strong.** Of the 15 focus pairs, six are
in the bottom decile of all 300:

| pair | r perivascular | r background | z | percentile of 300 |
|---|---:|---:|---:|---:|
| m1–m24 | −0.22 | +0.14 | −8.0 | **0.0** (most anti-associated pair of all 300) |
| m10–m23 | +0.26 | +0.39 | −7.3 | **0.3** |
| m1–m23 | +0.02 | +0.26 | −5.1 | 2.0 |
| m0–m1 | −0.14 | +0.16 | −4.7 | 2.7 |
| m1–m10 | −0.11 | +0.07 | −4.2 | 3.3 |
| m1–m6 | −0.08 | +0.04 | −1.8 | 6.3 |

So near vessels the vascular programs **segregate into different cells** rather than
stacking. Motif 1 (SMC/NOTCH3 wall) de-couples from all five others; and m10 (tumour–
neovessel, VEGFA→KDR) vs m23 (endothelial identity, CDH5/CD34/JAM2) is the second most
anti-associated pair in the entire matrix despite correlating at +0.39 in background tissue.
m23–m24 and m10–m24 are also below the 12th percentile. The wall, the endothelial identity,
the tumour-facing neovessel program and the alveolar–capillary program occupy distinct
perivascular cells.

Panel h of the figure exists solely to stop panel c being misread; it is the control, not
decoration. **Do not quote a bare z from this analysis** — quote the percentile or the
excess over the +3.6 median, both of which are columns in `*_pairs.csv`.

## CORRECTIONS to the two sections above, after independent verification (2026-09-02)

Four checkers re-derived every quoted number from the saved files. The core result survives;
several statements around it do not. Corrections, in order of how badly they mislead:

**1. "The vascular programs occupy different perivascular cells" is WRONG.** Occupying
different cells requires mutual exclusivity — O/E well below 1 and negative correlation.
The files show the opposite for most flagged pairs. Perivascular O/E: m1-m10 0.753,
m0-m1 0.871, m1-m24 0.906, m1-m6 0.916, m1-m23 0.986, **m10-m23 1.074**. Perivascular
partial correlations: only m1-m24 (−0.133) and m1-m10 (−0.010) are negative; m0-m1 +0.060,
m1-m6 +0.016, m1-m23 +0.112, **m10-m23 +0.233** are positive.

**"m10 and m23 segregate at vessels" is WRONG.** They co-occur at vessels 7% *above*
independence, and **95.6% of m10-positive perivascular cells are also m23-positive**. The
z = −7.3 says only that they are *less* correlated at vessels than in 4,000 matched
pseudo-rings (+0.26 vs a null mean near +0.37). The defensible sentence is
**"less coupled at vessels than expected"**, never "segregate" or "mutually exclusive".

A structural limit that makes the strong reading impossible anyway: the ON-state marginals
in the ring are 89.1% (m23), 78.2% (m1), 43.9% (m0), 40.0% (m24), 28.9% (m6), 20.2% (m10).
At those marginals, negative association is ceiling-compressed — a near-ubiquitous motif
*cannot* show strong anti-association whatever the biology.

Also: the bottom decile holds **7** focus pairs, not 6 (m23-m24 at percentile 9.67 was
omitted). And "survives removing cell type × tissue layer" is overstated — m1-m6 falls below
|z| = 2, m0-m1 drops 38%, m1-m23 drops 20%; only m10-m23 and m1-m24 strengthen.

**RESOLVED (2026-09-03).** The checker could not reproduce `R_peri`/`R_bg` because the
transform is not one of the six they tried. It is a **Pearson on within-section normal
scores** — rank → quantile → probit (van der Waerden) — with the ranks taken over *every*
cell of the section and deliberately **not** recomputed inside the subset, so perivascular
and background stay on one scale. Reproduced exactly (max |difference| **0.0**).

Two things fell out of doing this:

1. **The transform is doing real work.** Across the 300 pairs the shipped values agree with
   Spearman at rho = 0.965, but with a plain Pearson on the raw loadings at only **0.647**.
   On raw loadings m0–m6 is **−0.004** where the shipped value is **+0.462**: the zero
   inflation and the 8–9 decade range dominate a raw Pearson. Never quote these as "the
   correlation" without naming the transform.
2. **The builder casts loadings to float32**, which collapses many near-floor values into
   ties — motif 1 in P17_AIS goes from 1,213 duplicate values in float64 to **101,670**
   (55% of the section). Unintended. Recomputing in float64 moves the 15 focus pairs by
   **≤0.001 with no sign changes** (max 0.030 anywhere among the 300 pairs; rank agreement
   rho = 0.9996), so no conclusion depends on it — but a future run should use float64.
   Saved as `vessel_motif_cooccurrence_ring30_float64check.npz`.

**2. Motif 23 is badly under-sold by calling it "second".** It is **rank 1 of 25 in
Vas_Endo LRI content (85.5%)** where motif 1 is only rank 4 (26.2%); it is positive at more
vessels than motif 1 (85.2% vs 77.4% at ring 25); its per-section signed-rank test is
significant in **4/4** sections where motif 1 manages **3/4**; and it beats motif 1 outright
in P21_LUAD (AUROC 0.865 vs 0.695). The honest framing is a **division of labour** — motif 1
is stronger where present, motif 23 is present at more vessels and in more sections — not a
ranking.

**3. "Motif 1 in all four sections" is true for direction, false for primacy.** Positive in
all four, but **rank 1 in only three**: in P21_LUAD it is rank 4 behind m23 (0.865),
m17 (0.738) and m16 (0.733). The per-vessel signed-rank fails in P21_LUAD at ring 25
(p = 0.414) and ring 30 (p = 0.193), and in P17_AIS against the background reference
(p = 0.621).

**4. "Motifs 10 and 24 do not track the lumens" — right for 10, overstated for 24.**
Motif 24's per-vessel median (−0.249) is rank 8 of 25, above the grand median; its
non-vascular AUROC is 0.597 in P21_LUAD (above 0.5); and **inside the lumens it is enriched
in both LUAD sections** (22.8% vs 12.5% baseline, 1.83×; 34.3% vs 22.3%, 1.54×). Motif 24 is
inconsistent, not absent. Only motif 10 is depleted everywhere (AUROC 0.429/0.223/0.344/
0.288, all < 0.5).

**5. "Motif 10 is the most depleted" is wrong.** Its median of −0.89 at ring 80 is correct
but it ranks **8th** of 25; m5 (−1.69), m12 (−1.64), m19 (−1.51), m7, m8, m21 and m4 are all
more depleted. Seven motifs avoid vessels more strongly than the one named "tumour
vasculature".

**6. The per-vessel count "0-13 depending on radius" is wrong: the sweep gives 5-36**
(7/6/5/13/26/36 at 10/25/50/80/120/160 µm). And the rise is a **power artefact** — median
ring cells go 30 → 529 over that span while motif 1's own median score *falls* past 50 µm.

**The omission that matters most in Q2: not one vessel is individually significant for motif
1 or motif 23.** `m1_active` and `m23_active` sum to 0 in every file. The 8-13 significant
calls are m4, m7, m13, m15, m19, m20, m2, m3, m8, m9, m11 — motifs whose *overall* vessel
medians are strongly negative. Those calls are not evidence for the vascular story.

**7. The alignment statistics are cell-level, i.e. the same pseudoreplication that retired
`vessel_motif_profiles_ring25.csv`.** AUROC, enrichment and the band profiles pool
130,064-528,946 cells from only 58-83 lumens. **Lead with the radius sweep instead** — motif
1 rank 1 of 25 at every radius from 10 to 160 µm, over 115-308 vessels as the replicate unit.
That is the single best-supported statement in the whole analysis. Related: the 0-10 µm
non-vascular band holds only 180-412 cells; distance is confounded with tissue layer (the
0-10 µm band is 31% Tumor_region against 58% at >200 µm in P17_AIS); and "peak/far ratio"
was never defined — only one of ~10 plausible definitions reproduces the quoted 3.35-5.12×
(on all cells it is 4.7-8.7×).

## vessel_motif_clustering.py — can vessels be classified by motif composition? (2026-09-02)

Unit is the VESSEL. Feature = fraction of cells in the 30 µm ring positive for each of the
25 motifs, over 300 curated lumens. Plain methods only, as asked: variance ratios, z-scores,
PCA, k-means, bootstrap ARI.

**Step 1 — do vessels differ at all, beyond sampling noise?** A ring holds 20–190 cells, so
identical vessels would still scatter. Observed between-vessel variance ÷ binomial variance
expected at those ring sizes: **15–35× for all 25 motifs, median 23×**. Unambiguous yes.
This had to be checked first — with rings this small it could easily have been noise.

**Step 2 — k-means.** k = 2 is the only supported solution: silhouette 0.440 and bootstrap
ARI 0.963. Every k ≥ 3 collapses (silhouette 0.14–0.19, ARI 0.45–0.66). Two classes, not a
taxonomy.

| | class 0 (n=56) | class 1 (n=244) |
|---|---|---|
| motifs | m15 78%, m22 77%, m14 70%, m10 70%, m7 60% | m1 82% |
| motif 1 | 45% | 82% |
| lumen area | 835 µm² | 2,438 µm² |
| ring SMC / Vas_Endo / Tumor | .08 / .15 / .08 | .30 / .25 / .00 |

**Step 3 — and here is the finding that matters.** The split is nearly one-dimensional:
PC1 is 42% of variance, **24 of 25 PC1 loadings share a sign**, and 22 of 25 motifs are
higher in class 0. So the dominant axis is not "vessel type A vs B" — it is **overall motif
activity**, with **motif 1 the single motif that runs the other way**.

And that axis is largely technical-adjacent: it correlates with **local cell density at
ρ = −0.59** (distance to the 10th neighbour; class 0 rings 20.2 µm vs class 1 27.3 µm,
Kruskal p = 2.4e-21). Denser tissue puts more cells in a BPTF patch, so more LRIs are
detected and more motifs come out positive. **Sequencing depth is NOT the driver**
(ρ = +0.01, p = 0.81; Kruskal p = 0.18) — that was the other candidate and it is ruled out.
Also correlated: ring %Tumor_epi +0.57, %SMC −0.55, %Vas_Endo −0.45, lumen area −0.44.

**Does anything survive?** Yes, partly. The same split reappears when each section is
clustered independently — |cos| with the global axis 0.96 (P17_AIS), 0.90 (P17_LUAD), 0.82
(P21_LUAD), but only 0.45 in P21_AIS. And after regressing out size, ring composition and
section, 54% of the variance remains and k=2 still gives silhouette 0.248 (ARI 0.39 with the
original) — right at the "continuum, not separated groups" line.

**Answer to the question as asked.** Vessels genuinely differ. They can be classified, into
two groups. But the hoped-for pattern — some vessels carrying {A, B}, others {A, C} without
B — is **not** what the data shows. There is one gradient (how much motif activity, driven
substantially by cell density) plus **one genuine opposition: motif 1 against everything
else**. Motif 1 is high at large, SMC-rich, sparser mural vessels and low at small, dense,
tumour-embedded ones. That is the only vessel-intrinsic motif contrast the 300 vessels
support.

## Clustering redone on the SIX vasculature motifs only (2026-09-02)

The section above clustered on all 25 motifs, which was wrong for the question asked: the
request was the six vasculature-related motifs (0, 1, 6, 10, 23, 24), and letting 19
unrelated motifs drive the split buries the vascular structure. `--motifs 0,1,6,10,23,24`
now restricts the feature set; the overdispersion report still covers all 25 for context.

**The six-motif split is more biologically legible than the 25-motif one.** k = 2 again
(silhouette 0.263, bootstrap ARI 0.870), but the classes are now **condition**-associated
(Cramér's V 0.43 for AIS/LUAD, up from 0.25 on 25 motifs):

| | class 0 (n=120) | class 1 (n=180) |
|---|---|---|
| m24 / m0 / m6 / m10 | 76 / 69 / 50 / 41 % | 21 / 32 / 19 / 7 % |
| m1 | 60% | **85%** |
| sections | P17_AIS 59, P21_AIS 37 (AIS-leaning) | P17_LUAD 63, P21_LUAD 54 (LUAD-leaning) |
| lumen area | 1,304 µm² | 3,460 µm² |
| ring SMC | 0.12 | 0.36 |

**m23 is useless for classifying vessels.** It is positive in 92% of ring cells and its
median across vessels is 100%, so a median split gives 0% — it is on in essentially every
vessel and carries no between-vessel information. The script now detects saturation, drops
such motifs from the pattern enumeration and says so, rather than letting them contribute a
constant bit.

**The question as actually asked — "some vessels have A and B, others A and C but not B" —
is answered by enumerating combinations, not by clustering.** Each of the 5 informative
motifs was split at its across-vessel median (so every marginal is 50% by construction and
no motif dominates through prevalence), giving 32 patterns over 300 vessels:

| pattern | observed | expected | O/E | sections |
|---|---:|---:|---:|---|
| **m1 alone** | 48 | 9.4 | **5.1** | P17_LUAD 29, P21_AIS 12 |
| **m0+m6+m10+m24 (no m1)** | 36 | 9.4 | **3.8** | P17_AIS 24, P21_AIS 6 |
| (none) | 18 | 9.4 | 1.9 | P21_LUAD 12, P17_LUAD 6 |
| m1+m6+m10+m24 | 15 | 9.4 | 1.6 | |

χ² = 312, df = 26, **p = 8e-51** — the six motifs are emphatically not independent across
vessels.

**So the answer is yes, and it is one clean opposition rather than a rich combinatorics:**
the two most over-represented patterns are **motif 1 on its own** (5.1× expected) and
**everything-except-motif-1** (3.8× expected), and they are mutually exclusive by
construction. There is also a real third group of 18 vessels with **none** of the six on.
The AIS/LUAD leaning of the two classes makes this the most promising result of the vessel
work so far.

**Caveats that stand.** Silhouette 0.263 is at the "continuum, not separated groups" line —
read it as a gradient with two ends, not two species. And the class axis still tracks local
cell density (Kruskal p = 4.4e-17), ring %SMC (6.1e-15) and lumen area (5.8e-10); sequencing
depth remains ruled out (p = 0.096). The density confound has NOT gone away by restricting
the motif set.

## Transform sensitivity of the co-occurrence matrix (2026-09-03)

`--transform {normal,log1p,log10,raw}` on `vessel_motif_cooccurrence.py` and the panel
script. Three full runs (each with its own 4,000 pseudo-ring null) are on disk.

**log1p is a near-no-op on these loadings and should not be used.** `log(1+x) ≈ x` for
`x << 1`, and the loadings have median 2.3e-10, 99.9th percentile 0.042, max 0.31. Pearson
on log1p vs Pearson on raw: max difference 0.0019 across the 300 pairs, **Spearman
0.999993**. It was run because it was asked for; it is the no-transform baseline under
another name.

**The conclusion is robust to normal-scores vs log10, and not to log1p/raw.** Spearman of
the partial-correlation z across all 300 pairs:

| | log1p | log10 |
|---|---:|---:|
| normal | **+0.10** | **+0.72** |
| log1p | — | +0.28 |

Percentile of each focus pair within its own 300 (lower = more anti-associated):

| pair | normal | log10 | log1p |
|---|---:|---:|---:|
| m1–m24 | 0.0 | 3.0 | 7.3 |
| m10–m23 | 0.3 | 0.3 | **27.7** |
| m1–m23 | 2.0 | 1.3 | **37.7** |
| m0–m1 | 2.7 | 6.0 | **51.3** |
| m1–m10 | 3.3 | 7.3 | **30.0** |
| m6–m10 | 64.7 | 48.3 | **2.7** |

Normal scores and log10 agree on the whole story: motif 1 anti-associates with the others,
and m10–m23 is the most anti-associated pair (percentile 0.3 under both). Under log1p/raw
that structure disappears and m6–m10 — a mid-ranking *positive* pair under the other two —
becomes the most anti-associated. **A raw-scale Pearson on a variable spanning 8–9 decades
is dominated by a handful of high-loading cells, so pairs reorder.** That is the reason the
default is a rank-based transform, and the reason log10 corroborates it while log1p does not.

The global shift also changes: median z over the 300 pairs is +3.60 (normal), +1.07 (log10),
+0.87 (log1p), with 60% / 22% / 23% of pairs above +3. So the "perivascular cells are
globally more coherent" effect is itself partly a property of the rank transform, and is
weaker but still present on a log scale.

**Recommendation: quote the normal-score version, cite log10 as corroboration, and do not
use log1p or raw for this quantity.**

## vessel_motif_kmeans.py — rerun on max-normalised cell-level loadings (2026-09-03)

Replaces `vessel_motif_clustering.py` at the user's instruction; the old figures are deleted,
not archived. Unit is still the vessel, 30 µm ring, 300 vessels, motifs 0/1/6/10/23/24.
What changed: the feature is now the mean over ring cells of the **cell-level loading divided
by that motif's maximum over all cells**, replacing the fraction-of-cells-positive; and K is
set by **elbow** on the within-cluster sum of squares (kneedle: furthest point from the chord),
not by silhouette. Every figure is written separately.

**Vessels still differ far beyond sampling noise** — between-vessel variance of the per-vessel
mean ÷ the sampling variance of a mean at these ring sizes: m1 167, m0 155, m10 125, m23 107,
m24 91, m6 21.

**Elbow gives k = 4, but it is shallow and the choice is close.** Chord distances: k=2 0.248,
k=3 0.298, **k=4 0.302**, k=5 0.279. k=3 and k=4 differ by 0.004 — this is a smooth curve, so
"k = 4" should be quoted as "3 or 4", not as a determination.

**The clusters are extremely unbalanced: 225 / 6 / 10 / 59**, and this is a direct consequence
of the requested normalisation. (Cluster membership is deterministic within a process but
moved by 2 of 300 vessels between two runs of the script; treat marginal membership as
approximate.)

| cluster | n | m0 | m1 | Σ over the six | dominant motif |
|---|---:|---:|---:|---:|---|
| 0 | 225 | 0.0032 | 0.0134 | 0.033 | mixed: m1 40%, m24 29%, m23 16% |
| 1 | **6** | 0.0069 | **0.2941** | **0.320** | **m1, 92%** |
| 2 | **10** | **0.1161** | 0.0020 | 0.134 | **m0, 87%** |
| 3 | 59 | 0.0053 | 0.0880 | 0.116 | **m1, 76%** |

**The composition view is what makes the problem legible.** Clusters 1 and 3 have the *same*
composition — both m1-dominated (92% and 76%) — and differ almost entirely in magnitude
(Σ 0.320 vs 0.116). So k=4 is really two compositions (m1-high, m0-high) each split by level,
plus the mixed low-signal bulk. PC1 carries 73% of the variance and is essentially an
m1-magnitude axis; the PCA is an L, one arm m1 and one arm m0, with the bulk at the origin.
That is what k-means does to unstandardised, right-skewed features: it partitions magnitude,
not profile.

Dividing by the maximum puts each motif in [0, 1] but leaves the heavy right skew intact: the
mean normalised loading is 0.004–0.018 for a typical vessel, so a handful of vessels carrying
a few very high-loading cells sit 20–35× above the rest and k-means spends its clusters on
them. Cluster 1 (n=6) is defined by m1 at 0.294 against 0.0135 in the main cluster; cluster 2
(n=11) by m0 at 0.113 against 0.0032. **At this feature scaling k-means is doing outlier
detection, not finding vessel classes.** Standardising each feature (z-score) instead of, or
in addition to, the max division would change this — it was not requested and was not done.

Confounds are unchanged in kind but weaker than before: ring cell type Cramér's V 0.44,
cell density Kruskal p = 1.4e-14, ring %SMC p = 3.3e-16; section drops to V = 0.27 (was 0.48).

**The combination enumeration is the part that reproduces cleanly.** Independence is rejected
at chi2 = 343, df = 57, **p = 1.6e-42**, and the single most over-represented pattern is again
**m1 alone** (29 vessels, O/E = 6.2). So the m1-versus-the-rest opposition survives the change
of feature, even though the clustering itself does not.

Figures (separate, no composites): `elbow`, `silhouette_ari`, `pca`, `cluster_profiles`,
`confounds`, `combinations` — in `vessel_kmeans_ring30_m0-1-6-10-23-24/`.

## Re-run on the updated vessel annotation (2026-09-03)

The curated lumens were edited on 2026-09-03 (files rewritten 16:29–16:59). **What actually
changed is very small, and worth recording because the file timestamps overstate it.**

- The lumen count is **unchanged**: 309 total, 83 / 74 / 81 / 71 per section, same
  `vessel_id`s, same flagged count (5). Status moved by one lumen, seed-accepted → edited.
- The edit log shows **116 lumens re-saved** on 09-03 (24 / 14 / 57 / 23). But reconstructing
  the pre- and post-edit polygon from the append-only jsonl, **114 of those 116 are
  byte-identical in geometry** — they were re-confirmed, not reshaped.
- **Only 2 polygons genuinely changed**, both in P21_AIS: `P21_AIS__68`
  315,219 → 333,400 µm² (+5.8%, the largest lumen in that section, and the "large" example in
  the strata figure) and `P21_AIS__58` 48,320 → 47,932 µm² (−0.8%).

`.bak` is **not** usable as the "before" state — it is written after the edit and matched the
new file exactly. The append-only `<SEC>_lumens.jsonl`, keyed on the `ts` field, is the only
way to recover the prior geometry.

**Nothing in the conclusions moves.** Per-vessel motif-1 median: ring25 +0.466 (was +0.461),
ring30 +0.499 (+0.493), ring30_bg +0.481 (unchanged), ring80 +0.482 (unchanged); motif 1 is
still rank 1 of 25 in all four, and rank 1 at every sweep radius from 10 to 160 µm with the
sweep medians unchanged to four decimals. Co-occurrence: perivascular cells 20,403 (was
20,263), and every focus-pair percentile moves by at most 2.7 points — m1–m24 still 0.0,
m10–m23 still 0.3. k-means still gives the same four-cluster structure (m1-92% micro-cluster,
m0-dominated small cluster, m1-dominated mid cluster, mixed bulk) with membership shifting by
~7 vessels, which is within the marginal instability already recorded for this clustering.

**Two process fixes made at the same time.**

1. `plot_motif_vessel_crops.py` is new: it replaces the untracked scratchpad script that made
   the earlier `mv_<section>.png` and was lost when that session's scratchpad was deleted —
   exactly the risk flagged when those figures were first produced. It reads the h5ad and the
   lumen parquets directly, so it no longer depends on the built explorer HTML.
2. `rerun_vessel_analyses.sh` drives the whole re-run from tracked scripts (analyses, then
   sweep table, then figures) so the next annotation update is one command.

The `log1p` / `log10` co-occurrence variants were **removed** from the bundle rather than
re-run: they were built on the old annotation, and the conclusion was already that the
normal-score version is the one to quote.

## Cluster example vessels — spatial crops (2026-09-04)

`plot_cluster_vessel_examples.py` (tracked). Answers "give a few example spatial plots of
vessels from each k-means cluster, plotted like the `mv_<section>` crops".

Design: for each cluster, take the `--n-examples` (default 3) vessels whose six-motif
feature vector is closest in Euclidean distance to that cluster's centroid — i.e. the most
*typical* members, not a random or extreme sample. One figure per cluster; rows = motifs
0, 1, 6, 10, 23, 24 plus a cell-type reference row; columns = the example vessels. Subject
lumen outlined solid cyan (`#00e0ff` lw 1.8 over `#003c46` lw 0.6), every other lumen in
the field thin grey (`#9e9e9e` lw 0.7) so neighbouring vessels are visible but never
confused with the subject. Column headers give vessel_id / ring cell count / lumen area;
100 µm scale bar in the bottom-left panel. Reads the h5ad + parquets directly — it does
**not** depend on the explorer HTML (the earlier `mv_*` generator was an untracked
scratchpad script and was lost with its session; this replaces that capability properly).

Two field-of-view sets, because the clusters differ hugely in vessel size:

| dir suffix | FOV | note |
| --- | --- | --- |
| `_examples` | 700 µm | matches the original `mv_<section>` crops |
| `_examples_fov350` | 350 µm | **more legible**; cluster 1's vessels are 436–814 µm² and are nearly invisible at 700 µm |

Selected examples (identical in both sets — FOV does not affect selection):

| cluster | n | dominant | example vessels |
| --- | --- | --- | --- |
| 0 | 230 | m1 42%, Σ 0.035 | P17_AIS__75, P17_AIS__61, P17_AIS__18 |
| 1 | 6 | m1 92%, Σ 0.320 | P21_LUAD__0, P21_AIS__23, P21_AIS__71 |
| 2 | 52 | m1 79%, Σ 0.117 | P21_AIS__73, P21_AIS__61, P21_AIS__59 |
| 3 | 12 | m0 81%, Σ 0.134 | P21_AIS__41, P17_AIS__71, P17_AIS__56 |

Read visually: cluster 1's motif-1 row is conspicuously darker than its other motif rows;
cluster 3 is the mirror image — motif 0 dense and dark around the lumen, motif 1 nearly
empty. That is the intended check, and it passes: the clusters look like what the
composition barplot says they are.

**Caveat, unchanged from the clustering itself.** These are *illustrations of the
clusters*, not evidence for them. The k-means structure is still dominated by overall
loading magnitude (PC1 73%, cluster sizes 230/6/52/12, cluster axis ρ = −0.59 with local
cell density) — so a large part of what looks different between these crops is how *much*
signal the neighbourhood carries, not only which motif carries it. Do not present a crop
as a discovered vessel subtype.

Rendered on iris → **DejaVu Sans, not Arial** (Arial is not installed on the compute
nodes). Re-render before manuscript use, as with every figure in this directory.

## Cluster composition stacked bars — section and cell type (2026-09-04)

`plot_cluster_composition.py` (tracked). Two SEPARATE single-panel stacked bars, cluster on
x, no panel letters, short title each. Rings are rebuilt with the same geometry as
`vessel_motif_kmeans.py` (same H5, polygons, R = 30 µm, `MIN_RING = 20`) and joined to that
script's `vessels.csv` on `vessel_id` — nothing is re-clustered. 300 vessels, 20,412 ring
cells. Output `vessel_kmeans_ring30_m0-1-6-10-23-24_composition/`, numbers in
`cluster_composition.csv`.

**Cell-type bar uses the mean of per-vessel fractions, not pooled counts** — the vessel is
the replicate, so pooling would let the largest rings speak for the cluster. Both are in
the CSV (`ct_meanfrac_*` vs `ct_pooledfrac_*`).

Colours: sections use the fixed Okabe–Ito `SECCOL`; cell types use tab20 indexed by position
in the h5ad category list, the same assignment the interactive explorer uses, so a colour
means the same cell type in both. With 19 categories no palette is colourblind-safe — the
stack order is fixed and identical across clusters, and the legend is ordered to match.
`other` is drawn as a **hatched white** wedge, not a third grey: at first render it was
`#c8c8c8` and collided with tab20's `Tumor_epi` (`#c7c7c7`).

### Cell-type composition (mean per-vessel fraction of ring cells)

| cluster | n | SMC | Vas_Endo | T | Fibro |
| --- | --- | --- | --- | --- | --- |
| 0 (m1 42%, Σ 0.035) | 230 | 0.216 | 0.233 | 0.091 | 0.165 |
| 1 (m1 92%, Σ 0.320) | 6 | **0.617** | 0.221 | 0.030 | 0.042 |
| 2 (m1 79%, Σ 0.117) | 52 | **0.497** | 0.270 | 0.040 | 0.064 |
| 3 (m0 81%, Σ 0.134) | 12 | 0.021 | 0.152 | **0.507** | 0.059 |

This is the clearest read yet of what the clusters are. The motif-1 clusters (1 and 2) are
**muscularised vessels** — SMC 50–62% of ring cells against 22% in the background cluster —
which is exactly what motif 1 is (SMC/NOTCH3 arterial-mural: `JAG1→NOTCH3`, `GJA5→GJA5`,
`ANGPT1→TEK`, `EDN1→EDNRA/B`, SMC +1.51 log2). Cluster 3 is the opposite: **SMC 2%,
T cells 51%** — lymphocyte-cuffed vessels. So motif 0, which I had only ever justified as
"vasculature-adjacent" by its Vas_Endo LRI weight, is here carried by a T-cell-rich
perivascular compartment. Worth checking motif 0's LRI content against that directly; not
done yet.

### Section composition (fraction of vessels)

| cluster | P17_AIS | P17_LUAD | P21_AIS | P21_LUAD |
| --- | --- | --- | --- | --- |
| 0 | 0.278 | 0.287 | 0.196 | 0.239 |
| 1 | 0 | 0 | 0.333 | 0.667 |
| 2 | 0.154 | 0.077 | 0.538 | 0.231 |
| 3 | **0.917** | 0 | 0.083 | 0 |

**The clusters are section-confounded, and this needs to be said wherever they are shown.**
Crosstab 4×4: χ² = 63.5, df = 9, asymptotic p = 2.8e-10, permutation p < 5e-5 (20,000 label
shuffles), Cramér's V = 0.266. Cluster 3 is **11 of 12 vessels from P17_AIS alone**, and
cluster 1 is P21-only. Only the 230-vessel background cluster 0 is spread evenly.

So the honest statement is: the small clusters are **section-specific** groups of vessels
that differ in perivascular cell composition. With four sections from two patients there is
no way to separate "this vessel type occurs in AIS" from "this is what P17's AIS section
looks like" — n = 2 patients, and the replicate unit for any condition claim is the patient,
not the vessel. Do **not** write cluster 3 up as an AIS-specific lymphocyte-cuffed vessel
state; it is one section's worth of vessels.

Rendered on iris → DejaVu Sans, not Arial. Re-render before manuscript use.

## K sweep for the vessel clustering, K = 4–8 (2026-09-04)

The elbow was doing the choosing; this sweeps K instead. `vessel_motif_kmeans.py` gained
`--k`, which forces K and writes to a `_k<K>`-suffixed directory so the elbow run is never
overwritten; `plot_cluster_composition.py` and `plot_cluster_vessel_examples.py` gained the
same flag so each K gets the *same full set* of downstream figures. Driver:
`sweep_vessel_kmeans_k.sh` (K=5–8, two at a time; K=4 is the elbow run already on disk and
was not recomputed). All four completed, 0 tracebacks.

**Two bugs the sweep exposed, both fixed:**

- `CCOL` was a hard-coded **five**-colour list, so the PCA panel raised
  `IndexError: list index out of range` for every K ≥ 6 — K=6 died and produced no
  `vessels.csv`, which then cascaded into three `FileNotFoundError`s downstream. Extended to
  ten and indexed modulo its length.
- Past ~5 categories **no palette passes `validate_palette.py`**: Okabe–Ito(8) clears the
  normal-vision floor (worst pair 15.6 ≥ 15.0) but still fails the colourblind pairs; the
  naive extension of the old five failed both (worst normal pair 9.3). Shape is already
  spent on section in that figure, so cluster identity no longer rests on colour alone —
  **each centroid is now labelled with its cluster number**. Do not "fix" this by hunting
  for a better 8-colour palette; there isn't one.
- `plot_cluster_composition.py` now scales figure width with cluster count (89 mm holds four
  bars; eight collide). K=4 is unaffected — its width still evaluates to 89 mm.

The elbow/silhouette table is identical for every run (same feature matrix): elbow at
**k = 4** (chord distance 0.301 at k=4 vs 0.300 at k=3, 0.280 at k=5 — the elbow is
essentially flat over 3–5, which is itself a reason not to trust it). Silhouette peaks at
k=2 (0.608) and never recovers; bootstrap ARI is 0.75–0.88 throughout and *rises* to 0.88
at k=8, so stability does not pick a K either.

### What each K gives (Σ = summed mean max-normalised loading; cell-type = mean per-vessel ring fraction)

| K | clusters (n, dominant motif, Σ, SMC, T, top section) |
| --- | --- |
| 4 | 230 m1 .035 SMC .22 T .09 · **6 m1 .320 SMC .62** · 52 m1 .117 SMC .50 · **12 m0 .134 SMC .02 T .51 (P17_AIS 92%)** |
| 5 | 196 m24 .030 · 76 m1 .085 SMC .50 · 14 m1 .196 SMC .56 · 2 m1 .461 SMC .60 · **12 m0 .134 T .51** |
| 6 | 181 m1 .025 · 67 m1 .086 SMC .51 · **24 m24 .086 SMC .06** · 14 m1 .196 · 2 m1 .461 · **12 m0 .134 T .51** |
| 7 | 157 m1 .023 · 65 m1 .056 · 34 m1 .124 · **23 m24 .087 SMC .07** · 12 m0 .134 T .51 · 7 m1 .234 · 2 m1 .461 |
| 8 | 154 m1 .022 · 63 m1 .055 · 34 m1 .124 · **23 m24 .087** · 15 m0 .104 T .41 · 7 m1 .234 · 2 m1 .461 · 2 m0 .234 T .70 |

**Read across K, there are exactly two things that are not the magnitude gradient:**

1. **The m0 / T-cell cluster is the one genuinely stable structure.** n = 12, SMC 0.02,
   T 0.51, 92% P17_AIS — *bit-identical* at K = 4, 5, 6 and 7. At K=8 it splits into 15 + 2,
   the pair being an extreme (T 0.70). If anything here is a real vessel type, it is this.
2. **A motif-24 cluster appears only at K ≥ 6** and is then stable: n = 24 → 23 → 23,
   SMC 0.06–0.07 (vs 0.50 for the m1 clusters), Vas_Endo 0.15, P21_AIS 62–65%. This is the
   sweep's one new finding — the elbow's K=4 hides it inside the 230-vessel background
   cluster. Motif 24 was previously characterised as *healthy* vasculature, and a
   low-SMC/low-Vas_Endo perivascular compartment is consistent with that; worth a look.

Everything else is the **same magnitude gradient being sliced finer**: the m1/SMC clusters
form a clean Σ ladder (.022 → .055 → .124 → .234 → .461) with SMC essentially constant at
0.50–0.60 across all of them. Raising K peels off ever-smaller extreme groups (7 vessels,
then 2) rather than finding new kinds of vessel. That is the same conclusion the K=4 run
reached (PC1 73%, cluster axis ρ = −0.59 with cell density) — the sweep confirms it rather
than overturning it.

**The section confound gets worse, not better, as K rises.** The n=2 cluster is 100%
P21_LUAD, the n=2 m0 cluster is 100% P17_AIS, the m0 cluster is 92% P17_AIS, the m24 cluster
is ~65% P21_AIS. With 4 sections from 2 patients none of these can be called a condition- or
patient-level finding; the replicate unit is the patient (n = 2).

Outputs: `vessel_kmeans_ring30_m0-1-6-10-23-24_k{5,6,7,8}{,_examples,_examples_fov350,_composition}/`,
bundled as `figures_vessel_kmeans_{4,5,6,7,8}/`. The bundle grew 244 MB → 604 MB, almost all
of it the vector (`pdf`/`svg`) example crops — those are scatter plots of tens of thousands
of points, ~6–9 MB each. PNG-only would be about a fifth of that.

Rendered on iris → DejaVu Sans, not Arial. Re-render before manuscript use.

### The motif-24 cluster is alveolar, and that is the sweep's real result

Checked after the fact, because `Alveolar_epi` is not in the K=4 figure's top-10 cell types
and only becomes visible once the cluster exists:

| K | cluster | n | top ring cell types |
| --- | --- | --- | --- |
| 6 | 2 | 24 | **Alveolar_epi 0.17**, Vas_Endo 0.15, Fibro 0.15, T 0.10, Macro 0.08 |
| 7 | 5 | 23 | **Alveolar_epi 0.18**, Vas_Endo 0.16, Fibro 0.15, T 0.10, Macro 0.08 |
| 8 | 4 | 23 | **Alveolar_epi 0.18**, Vas_Endo 0.16, Fibro 0.15, T 0.10, Macro 0.08 |

Background (K=4 cluster 0, n=230): `Alveolar_epi` **0.03**. So this cluster is ~6× enriched
for alveolar epithelium, and alveolar epithelium is its *most abundant* ring cell type —
these are vessels sitting in alveolar parenchyma rather than in tumour or stroma. That is
independently consistent with motif 24 having been called **healthy** vasculature, and it is
reached here from cell-type composition alone, without using the motif's LRI content.

So the sweep is worth having: K=4 was hiding a small but coherent alveolar-vessel group
inside its 230-vessel background cluster. The caveat stands — n = 23 vessels, ~65% from
P21_AIS, so it is not a condition-level claim.

## New feature definition: GMM positive FRACTION, not continuous loading (2026-09-04)

`vessel_motif_kmeans.py --features fraction`. Cells are binarised with the GMM ON/OFF call
**already stored in obs** — `motif_<k>_state`, a categorical with categories
`['negative', 'positive']` (checked explicitly; the script refuses to run if that is not the
encoding). Each vessel's 30 µm ring is then summarised as the **fraction of its cells
positive** for each of motifs 0, 1, 6, 10, 23, 24, and those fraction vectors are clustered.
No normalisation — the features are already fractions on a common 0–1 scale, so dividing by
the max (which is 1) would be a no-op that only obscures what was done.

Cell-level ON fractions overall: m0 0.487, m1 0.350, m6 0.461, m10 0.420, m23 0.633,
m24 0.370. Output `..._frac{,_k<K>}/`; sweep driver `sweep_vessel_kmeans_frac.sh`
(K = 2–10, K=3 unforced/elbow). All nine completed, 0 tracebacks. The example crops needed
no change: they already grey out GMM-negative cells and shade only the positives, so they
read correctly for this feature definition.

### This is a different partition, and a better-behaved one

**ARI(loading K=4, fraction K=3) = 0.006.** The two methods are not variations on a theme —
they agree no better than chance. Crucially the fraction clustering is **not** the magnitude
artefact the loading one was:

| confound | loading K=4 | fraction K=3 |
| --- | --- | --- |
| `mean_counts` (depth) | ρ = +0.01, p = 0.18 | Kruskal **p = 0.31** |
| `n_ring` (ring size) | — | Kruskal **p = 0.50** |
| cluster sizes | 230 / 6 / 52 / 12 | 100 / 149 / 51 |

and it stays balanced all the way to K=10 (smallest cluster 8, most 20–70) instead of peeling
off n=2 outlier groups. Combination enumeration is also far more interpretable on binary
features: "m1 alone" O/E 5.1 (48 vessels) and "m0+m6+m10+m24" — everything *except* m1 —
O/E 3.9 (36 vessels), χ² = 313, df = 26, p = 5.3e-51. That is the complementary-pattern
structure the original question asked for ("some vessels have A and B, others A and C but
no B"), and the loading features never produced it.

K=3 (elbow; chord distance 0.286 at k=3 vs 0.282 at k=4 — flat again, hence the sweep):

| cluster | n | m0 | m1 | m6 | m10 | m23 | m24 | SMC | lumen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 100 | .54 | .92 | .46 | .27 | .95 | **.81** | .21 | 1,532 µm² |
| 1 | 149 | .33 | .85 | .18 | .06 | .90 | **.11** | .37 | 3,683 µm² |
| 2 | 51 | .76 | **.13** | .39 | .53 | .92 | .64 | **.00** | 890 µm² |

m23 is saturated (~0.9 in every cluster) and carries no information. The discriminating axes
are **m1** (which tracks SMC content: m1-present clusters have SMC 0.30–0.53, m1-absent
0.00–0.03) and **m24**.

### The result: motif 24 in the perivascular ring separates AIS from LUAD

Across **all 54 clusters** produced at K = 2…10:

| | clusters | %LUAD range | median |
| --- | --- | --- | --- |
| m24 positive fraction > 0.5 | 29 | **0 – 35%** | 12% |
| m24 positive fraction ≤ 0.5 | 25 | **65 – 100%** | 72% |

**Zero overlap** — the highest %LUAD of any m24-present cluster (35%) is below the lowest of
any m24-absent cluster (65%), at every K from 2 to 10. Per vessel, m24 ring fraction
separates the conditions **within each patient separately**:

| patient | AIS median | LUAD median | AUROC |
| --- | --- | --- | --- |
| P17 | 0.758 (n=83) | 0.000 (n=70) | **0.960** |
| P21 | 0.749 (n=76) | 0.094 (n=71) | **0.795** |

Motif 1 does **not** replicate — it rises in P17 (0.645 → 0.843) and falls in P21
(0.920 → 0.606); m10 falls in P17 (0.466 → 0.089) and is flat in P21 (0.122 → 0.115).
**m24 is the only one of the six that moves the same way, and strongly, in both patients.**
Motif 24 was independently characterised as *healthy* vasculature, and the K≥6 loading-based
sweep separately found its cluster to be alveolar (Alveolar_epi 0.18 vs 0.03 background). All
three lines agree: the healthy/alveolar vascular program is switched off around vessels in
invasive tumour.

**Statistical caveat, and it is the usual one.** Every p-value above is computed over
*vessels*, and vessels within a section are spatially autocorrelated — this is the
pseudoreplication trap that has already produced two wrong answers in this project. The
Mann-Whitney over clusters (p = 3.3e-10) is worse, not better: clusters at different K are
not independent observations. **The defensible claim is the effect size and its direction
replicating independently in two patients**, plus the complete separation across nine
clusterings. The replicate unit for a condition claim is the patient, and n = 2.

### Stability across K

Consecutive-K ARI: 0.659, 0.620, 0.586, 0.694, **0.895, 0.904, 0.937**, 0.857 for
K=2→3 … 9→10. Past K=6 the partitions are refinements of each other rather than reshuffles.
Recurring archetypes, stable in size wherever they appear:

| archetype | signature | n across K≥6 |
| --- | --- | --- |
| pure arterial/mural | m1+m23 only, SMC 0.50–0.53 | 71, 69, 65, 66, 64 |
| non-muscular angiogenic | m0+m10+m23+m24, m1 off, SMC 0.00 | 29, 28, 26, 24, 24 |
| all-on | m0+m1+m6+m23+m24 | 32, 32, 31, 23 |
| m23-only (silent) | m23 alone, 100% LUAD | 25, 25, 25 |

Bundled as `figures_vessel_kmeans_frac_{2..10}/`. Example crops are **fov350 for every K**;
fov700 only for K=3, because the vector crops are 6–9 MB each and shipping both for nine K
would add ~1 GB. Any missing set is one command:
`plot_cluster_vessel_examples.py --ring 30 --features fraction --k <K>`.

Rendered on iris → DejaVu Sans, not Arial. Re-render before manuscript use.

### Controls on the m24 result

Per-vessel Spearman of m24 ring fraction against the covariates:

| covariate | ρ | p |
| --- | --- | --- |
| `n_ring` (ring size) | −0.011 | 0.86 |
| `mean_counts` (depth) | +0.091 | 0.12 |
| `frac_SMC` | −0.282 | 6.7e-07 |
| `frac_VasEndo` | −0.250 | 1.2e-05 |
| `frac_Tumor` | **+0.453** | 1.4e-16 |
| `nn10_um` (local density) | −0.372 | 2.7e-11 |
| `lumen_area` | −0.319 | 1.8e-08 |

**Not** a ring-size or sequencing-depth effect — those are the two that would have made it an
artefact, and both are flat. The `frac_Tumor` correlation is *positive*, which looks wrong
until you remember AIS is adenocarcinoma **in situ** and therefore carries `Tumor_epi` of its
own; m24-high rings are AIS rings, so they have more annotated tumour epithelium, not less.

Restricting to the 239 vessels whose rings contain **no** `Tumor_epi` at all, m24 still
separates the conditions essentially unchanged — P17 AUROC 0.934 (was 0.960), P21 0.789 (was
0.795). So the signal is not tumour cells displacing normal cells in the ring.

The ring cell-type composition of the three K=3 clusters is also strikingly *undramatic* —
cluster 0 (m24-high, AIS) SMC 0.25 / Vas_Endo 0.22 / Alveolar_epi 0.06, cluster 1 (m24-low,
LUAD) SMC 0.36 / Vas_Endo 0.27 / Tumor_epi 0.02. Both are SMC/endothelial rings. **These
vessels differ in which signalling program is on, not in what cells surround them**, which is
the more interesting version of the result and the opposite of what the loading-based
clusters showed (there, composition was the whole story). The remaining real confound is
local density (ρ = −0.372) and lumen size (−0.319): LUAD vessels here are larger and sit in
sparser tissue.

## Niche clustering of the perivascular cells, by cell type only (2026-09-04)

`vessel_niche_cluster.py`, adapted from `/home/fanj2/cluster_analysis.py` with the smallest
set of changes that lets it run on this question. **No motifs are used.** The method is that
script's: count cell types in each cell's r-µm neighbourhood, z-score the columns, k-means
with K chosen by `KneeLocator` over a 10-restart inertia curve.

Seven changes, all listed in the script's docstring; the ones that matter:

- **Scope** — only cells within 30 µm of a curated lumen boundary and outside the lumen are
  clustered, the same perivascular set every other analysis here uses. Neighbourhood counts
  are still searched against **all** cells, so a ring cell's composition is its true local
  composition, not a composition of ring cells only. (I build them only for the ring rows;
  each such row is identical to what the original, which built one per cell, would give.)
- **Per section** — one KDTree per section. The four sections share a coordinate frame, so a
  single global tree would merge cells millimetres apart in different tissue. This forced
  `neighbor_cell_types` to take an explicit `cell_type_map`, so all sections yield the same
  columns in the same order.
- **Seed** — `autotuned_kmeans` now fixes `random_state` (the original left KMeans unseeded).
- `estimate_cluster_associations` read `n_clusters` / `cluster_labels` from module globals
  and **could not run as written**; they are arguments now, and its permutation loop shuffles
  once per trial and scores all K×K pairs from that shuffle instead of reshuffling per pair.
  Each pair's marginal null is unchanged; it turns an hour into seconds.
- **Dropped `plot_meta_cluster_volcanos`** — it calls `differential_expression` and
  `volcano_plot`, neither of which is defined or imported in `cluster_analysis.py`, so that
  function cannot run as given. It also needs the expression matrix, which this question does
  not.

**Result.** 20,279 perivascular cells over 300 vessels, 19 cell types, neighbourhood 30 µm.
KneeLocator picked **K = 11**.

| niche | n | top cell types |
| --- | --- | --- |
| 4 | 6,272 | **SMC 0.67**, Vas_Endo 0.15 |
| 5 | 4,296 | Fibro 0.25, Plasma 0.16, Vas_Endo 0.16 |
| 8 | 4,274 | **Vas_Endo 0.39, SMC 0.36** |
| 1 | 1,268 | **Alveolar_epi 0.24**, Vas_Endo 0.16 |
| 6 | 1,209 | **T 0.39**, Neutro 0.10 |
| 0 | 920 | Plasma 0.21, T 0.18 |
| 3 | 663 | **Tumor_epi 0.38** |
| 2 | 428 | T 0.32 |
| 7 | 409 | T 0.20, Fibro 0.16 |
| 9 | 398 | T 0.24, Tumor_epi 0.12 |
| 10 | 142 | **Airway_epi 0.54** |

These are the expected perivascular compartments and they are cleanly separated: a muscular
wall (4), an endothelial/mural wall (8), fibrous-plasma stroma (5), alveolar (1), lymphoid
(6, 2), tumour (3), airway (10). Niches are spatially coherent — mean self-association
log2FC **+3.53**, and 26 of 121 adjacency pairs reach p < 0.05.

Rolled up per vessel, **189 of 300 vessels (63%) have a ring that is >50% a single niche**,
dominated by niche 5 (95 vessels), niche 4 (85) and niche 8 (46) — so most vessels do have a
dominant perivascular niche rather than a mixture.

**Same confound as everything else here:** niche vs section Cramér's V = 0.35, and niche vs
condition 0.35 (both p ≈ 0). Niche 10 is 99% P17_AIS, niche 3 is 83%, niche 6 is 84%. With
4 sections from 2 patients the replicate unit for a condition claim is the patient.

Outputs in `vessel_niche_ring30_nb30/`: `niche_labels.csv` (one row per perivascular cell,
with its `vessel_id`), `niche_cell_type_proportions.csv`, `neighborhood_counts.npy`,
`elbow_scores.csv`, and the reference script's four figures — spatial assignments per
section, per-niche cell-type composition, per-section niche proportions, and the niche
adjacency heatmap. `--k` forces K if the elbow is not wanted.

### Per vessel, not per cell (the unit that was actually wanted)

`--unit vessel` (now the default; `--unit cell` reproduces the run above). One row per
vessel: **the vessel's neighbourhood IS its 30 µm ring**, so its feature is simply the
cell-type composition of the ring cells, and `--neighborhood` is unused. Each ring cell is
assigned to its nearest lumen, so no cell is counted for two vessels. 300 vessels x 19 types.

**One deliberate deviation from the reference, flagged because it changes the answer.** Ring
size spans **20-262 cells, a 13x range**, so raw counts would mostly encode how big the ring
is -- the exact artefact that sank the loading-based motif clustering. `--composition`
therefore defaults to `proportion` for `--unit vessel`; `--composition count` restores the
reference's raw-count behaviour. (For `--unit cell` the default stays `count`, as in the
original.)

`estimate_cluster_associations` also had to find neighbours **within a section** for this
unit -- with 300 vessels across four sections sharing a coordinate frame, a global k=20
would have made vessels in different tissue "adjacent".

**KneeLocator picked K = 9, and should not be trusted here:** the inertia curve is close to
linear (5700, 4994, 4503, 4232, 3997, 3765, ...), and K=9 yields clusters of **n = 1 and
n = 2**. Swept K=3-8:

| K | usable? | clusters |
| --- | --- | --- |
| **3** | **yes, the only K with no degenerate cluster** | 141 / 110 / 49 |
| 4 | no | 143 / 109 / 47 / **1** |
| 5 | no | 133 / 109 / 45 / 12 / **1** |
| 6 | no | 134 / 109 / 42 / 12 / **2** / **1** |
| 7 | no | 106 / 101 / 41 / 26 / 13 / 12 / **1** |
| 8 | no | 104 / 74 / 39 / 39 / 31 / 10 / **2** / **1** |

The n=1 cluster is the same single airway-adjacent vessel (Airway_epi 0.45) at every K >= 4;
it is a real outlier, not a bug.

**K = 3 — three perivascular niches:**

| niche | n | composition | sections |
| --- | --- | --- | --- |
| 2 | 110 | **SMC 0.55, Vas_Endo 0.28** — muscular wall | spread |
| 0 | 141 | Vas_Endo 0.24, Fibro 0.22, Plasma 0.11 — fibrous/plasma stroma | spread |
| 1 | 49 | **T 0.28**, Neutro 0.07 — immune-cuffed | **76% P17_AIS** |

K=8 adds interpretable types worth knowing about even though it carries two tiny clusters:
an **alveolar** niche (n=31, Alveolar_epi 0.18, 52% P21_AIS), a **plasma-cell** niche (n=39,
Plasma 0.26, **72% P21_LUAD**) and a **tumour** niche (n=2). The alveolar and plasma niches
are the per-vessel echo of what the motif work found — alveolar rings in AIS, and the LUAD
rings differing.

Same confound throughout: the immune-cuffed niche is 76-82% P17_AIS at every K, the plasma
niche 72% P21_LUAD. n = 2 patients.

Outputs: `vessel_niche_ring30_pervessel{,_k3..k8}/` and `vessel_niche_ring30_nb30_percell/`.

Figures here follow the **reference script's** style (its `plt.savefig` calls, PDF only), not
this directory's png+pdf+svg saver — that was the "don't change too much" instruction. Say
the word and they can go through `save_all_formats` instead.

## Composition figures unified with the niche clustering's format (2026-09-04)

`plot_cluster_composition.py` now **imports** `plot_cluster_cell_types` and
`plot_meta_clusters_per_punch` from `vessel_niche_cluster.py` instead of drawing its own, so
`cluster_celltype_composition` and `figures_niche/*/neighborhood_cluster_cell_types` are the
same figure produced by the same code in two analyses, not two similar figures. Both are
rendered to png + pdf + svg by calling the function once per extension — the functions
themselves are left untouched (they take an `outfile`, as in the reference).

Consequences of adopting that format, all deliberate: **all 19 cell types are shown** (no
top-10 + hatched "other"), the x axis is the plain cluster index rather than
`cluster c (n = …)`, and the legend sits outside right titled "Cell type". The statistic is
unchanged and already matched the niche script's — the mean of the **per-vessel** fractions,
so the vessel stays the replicate and big rings do not speak for the cluster.

**New figure, `neighborhood_cluster_proportions`** — the reference's per-section view:
x = section, stacks = cluster, heights = proportion of that section's vessels. It is the
transpose of `cluster_section_composition` (x = cluster, stacks = section), which is kept.
On the fraction clustering at K=3 it states the AIS→LUAD result in one panel: cluster 1
(m1+m23, m24 off) is 87% of P17_LUAD and 70% of P21_LUAD, while cluster 0 (m24-high) is 46%
of P17_AIS and 66% of P21_AIS — the same direction in both patients.

Re-rendered for all 18 runs (`regen_composition_figs.sh`), 0 tracebacks.

## CORRECTION (2026-09-08) — audit of the positive-fraction vessel clustering

Re-checked on request. **The method is exactly as described and reproduces bit-for-bit; the
biological headline drawn from it was wrong.**

### What the method actually does (verified against the code and re-derived from the h5ad)

1. `L[i] = (motif_<k>_state code == 1)` — the ON/OFF call already in `obs`, categories
   `['negative','positive']`, asserted at load. **This state is a clean GLOBAL threshold on
   the loading**: for all six motifs `max(loading | negative) < min(loading | positive)`,
   the gap sitting in the ~7th significant figure (m0 1.899527e-10 vs 1.899533e-10). So it
   is one cut per motif over all 1,676,162 cells — *not* per section, per patient or per
   cell type. Cell-level ON: m0 .487, m1 .350, m6 .461, m10 .420, m23 .633, m24 .370.
2. **No normalisation** (`Ln, mx = L, ones`) — deliberate, they are already 0/1.
3. Ring geometry per curated lumen: KD-tree over the densified polygon boundary; a cell is
   in the ring if `0 < d_boundary <= 30 µm` and it is not inside the polygon; the vessel is
   skipped if fewer than `MIN_RING = 20` cells qualify. 300 vessels survive.
4. Feature = `Ln[:, ring].mean(1)` → **the fraction of that ring's cells that are positive**.
5. K by elbow (`elbow_k`, max distance to the chord) over `KMeans(k, n_init=25,
   random_state=0)`, k = 1..10; final fit `KMeans(kbest, n_init=50, random_state=0)`.

Reproduction: rebuilt (300, 6) from the h5ad, `max |rebuilt − features.npy| = 0.0`, vessel
order identical, `n_ring` identical, and `KMeans(3, n_init=50, random_state=0)` reproduces
the stored labels exactly. Ring overlap is negligible — 133 of 20,279 cells (0.7%) fall in
more than one vessel's ring, so the 20,412 ring-cell slots inflate the cell count 1.007×.

### The bug in the interpretation: the features are absolute, never contrasted

A vessel's feature is its ring's raw level, so it inherits its section's overall level.
Section-wide ON fraction over **all** cells:

| section | m0 | m1 | m6 | m10 | m23 | **m24** |
| --- | --- | --- | --- | --- | --- | --- |
| P17_AIS | .527 | .293 | .493 | .610 | .731 | **.819** |
| P17_LUAD | .692 | .403 | .428 | .406 | .632 | **.120** |
| P21_AIS | .416 | .419 | .681 | .324 | .780 | **.877** |
| P21_LUAD | .278 | .271 | .374 | .424 | .527 | **.243** |

The m24 AIS→LUAD collapse is a property of the **whole tissue**. The perivascular rings
(.702/.089, .648/.218) merely track it, and slightly *below* it. Across vessels,
`corr(ring m24, its section's background m24) = 0.683`, the highest of the six.

Subtracting each vessel's own section background:

| motif | raw AUROC P17 / P21 | section-centred P17 / P21 |
| --- | --- | --- |
| m24 | **0.960 / 0.795** | **0.487 / 0.400** |

The separation is gone — at or below chance in both patients.

### The local control, which is the design this directory already established

Each vessel's 0–30 µm ring against its **own** 100–250 µm annulus (other lumens' 30 µm
neighbourhoods excluded), 300 vessels, paired Wilcoxon:

| motif | ring | surround | ring − surround | p | AIS/LUAD AUROC on the contrast, P17 / P21 |
| --- | --- | --- | --- | --- | --- |
| m1 | .752 | .359 | **+0.392** | 1.8e-38 | 0.527 / 0.830 |
| m23 | .920 | .693 | **+0.227** | 2.7e-40 | 0.410 / 0.400 |
| m0 | .470 | .503 | −0.033 | 9.4e-02 | 0.813 / 0.473 |
| m24 | .431 | .528 | **−0.097** | 1.1e-08 | 0.439 / 0.267 |
| m10 | .208 | .354 | −0.146 | 2.1e-22 | 0.625 / 0.296 |
| m6 | .310 | .474 | −0.164 | 2.4e-20 | 0.573 / 0.219 |

**Motif 24 is DEPLETED around vessels**, not enriched. And on the vessel-vs-surround
contrast **no motif gives a replicated AIS/LUAD difference** — every one lands on opposite
sides of 0.5 in the two patients, m24 included (0.439 / 0.267, both the *reverse* of the
raw claim).

### What stands and what does not

- **STANDS:** motifs **1** (SMC/mural) and **23** (endothelial) are genuinely perivascular,
  +0.392 and +0.227 over their own local surround. This is the same result the ring-radius
  sweep gave (motif 1 rank 1 of 25 at every radius 10–160 µm) and it does not depend on the
  feature definition.
- **STANDS:** the clustering is a correct, reproducible implementation of what was asked.
- **WITHDRAWN:** "motif 24 in the perivascular ring separates AIS from LUAD (AUROC 0.960 /
  0.795), the only motif replicating in both patients." The rings do separate, but only
  because their sections do; vessels are not different from their surroundings, and the
  effect reverses under the local control. The earlier reading — "the alveolar program is
  switched off around vessels in invasive tumour" — should be **"the alveolar program is
  switched off across the invasive tumour section, vessels included, and slightly more so
  away from vessels than at them."**
- Also withdrawn by implication: the "perfect no-overlap m24 / %LUAD dichotomy across all
  clusters at K=2..30". It is real as arithmetic and meaningless as biology — it is
  restating the section identity of each cluster.

Root cause, and it is the same one as the two earlier errors in this file: **a between-group
comparison was made on an absolute per-unit quantity with no matched local reference.** The
fix is the design already in `vessel_motif_local.py` — ring vs its own surround — which
should be applied to the fraction features too if a condition claim is ever wanted.

## `cluster_motif_composition_absolute` — the stack in true values (2026-09-08)

`cluster_motif_composition` normalises each cluster's bar to 100%, which throws away how
much signal a cluster carries. The new `cluster_motif_composition_absolute` stacks the raw
centroid values instead, so **the bar's length is the cluster's summed motif level**. Same
colours, same motif annotations, same legend, still no per-segment numbers and no total
annotation — the length is the total.

For `--features fraction` each segment is a positive fraction in [0, 1], so the bar lies in
**[0, n_motifs] = [0, 6]**, and the axis is fixed to that full range so bars stay comparable
across clusters and across K. This is the figure that matches the "let each motif have its
own [0,1] value, so the sum of each vector is [0,6]" framing directly. At K=3 the totals are
cluster 0 = 3.93, cluster 2 = 3.37, cluster 1 = 2.42 — cluster 1 is not just differently
composed, it has substantially *less* motif activity overall, which the percentage version
hides completely.

For `--features loading` the values are means of loading/max and occupy only a sliver near
zero, so a fixed 0–6 axis would render everything invisible; the axis is left automatic
there. That version is worth keeping precisely because it makes the earlier diagnosis visual:
the four loading clusters form a clean magnitude ladder (0.035 → 0.117 → 0.134 → 0.320) with
m1 filling almost the whole bar — one gradient, not a taxonomy.

Re-ran `vessel_motif_kmeans.py` for all 18 clusterings (`regen_kmeans_figs.sh`, does not
touch the example crops), 0 tracebacks; clustering is deterministic so `vessels.csv` and
`features.npy` came back unchanged (fraction K=3 still 149/100/51).

## Figure layout fixes for high K (2026-09-08)

Skills: `nature_publication_figures` is **not vendored in this checkout**; `dataviz` is, and
was followed (form, fixed-order categorical hues, validated palette, 2 px surface gap between
stacked fills, recessive text, render-and-look-at-it).

- **Y labels collided at K > 10.** Figure height was fixed at `SINGLE * .78` regardless of
  cluster count. Now `H_STK = max(45, 12 + 5.0·K) mm` for the stacked figures and
  `H_LEV = max(45, 12 + 1.9·n_motifs·K) mm` for the grouped one (which needs a row per
  cluster × motif bar). At K=30 that is 162 mm and 354 mm.
- **Cluster tick labels** went from two lines (`cluster 0` / `(n=100)`) to one (`0 (n=100)`)
  with `ylabel = "cluster"`, halving the vertical demand per row.
- **Legend offset must be physical, not fractional.** An axes-anchored legend at a fixed
  axes fraction cleared the x label at K=30 and sat on top of it at K=3. It is now a
  figure-level legend anchored 5 mm below the figure, computed from `fig.get_figheight()`.
- **Text cut to labels, not sentences:** "Motif level per cluster" → "Motif level";
  "Motif composition per cluster" → "Motif share"; "Motif level per cluster, stacked" →
  "Motif level, stacked"; "fraction of ring cells positive (GMM)" → "positive fraction";
  "share of the cluster's total positive fraction (%)" → "share (%)"; "sum of the 6 motif
  positive fractions (0-6)" → "Σ positive fraction"; "Cluster quality (does not set k)" →
  "Cluster quality"; "Vessels: shape = section, colour and label = cluster" → "Vessels";
  "What else does the cluster axis track?" → "Cluster axis vs covariates"; "observed ÷
  expected if motifs were independent" → "observed / expected"; "elbow, k = N" → "k = N".
- **2 px surface gap** (`edgecolor="white", linewidth=.5`) added between stacked segments.

`MCOL` (`#0072b2 #e69f00 #009e73 #cc79a7 #56b4e9 #d55e00`) re-validated with the repo's
Python port of the dataviz validator, `--pairs all`: **PASS**, normal-vision floor 15.6 vs a
15.0 requirement. Unchanged.

All 18 clusterings re-rendered (`regen_kmeans_figs.sh`), 0 tracebacks, `vessels.csv` /
`features.npy` unchanged (fraction K=3 still 149/100/51).

## Motif clusters (fraction, K=10) vs cell-type niche clusters (2026-09-08)

`compare_motif_vs_niche_clusters.py`. Both partition the same 300 vessels from disjoint
information — motif positive fractions vs ring cell-type composition. Compared with
partition measures, not Pearson: cluster ids are nominal.

| measure | value |
| --- | --- |
| ARI | **+0.158** |
| AMI | +0.231 |
| Cramér's V | 0.387 (χ² = 359, df = 72, p = 1.1e-39) |
| ARI **within section**, weighted | **+0.134** (P17_AIS +0.151, P17_LUAD +0.157, P21_AIS +0.092, P21_LUAD +0.138) |
| permutation p, labels shuffled **within section** | **< 5e-5** (20,000 shuffles; null ARI 0.033 ± 0.009) |

**The association survives the section confound** — within-section ARI 0.134 against a
within-section null of 0.033 — which is the opposite of what happened to the motif-24
AIS/LUAD claim. Magnitude is still modest: ARI 0.158 means the two partitions agree far
above chance but are nowhere near the same partition, i.e. motif usage carries information
that ring cell-type composition does not, and vice versa.

Strongest pairs (n >= 5), motif signature against niche composition:

| n | log2 O/E | motif cluster | niche |
| --- | --- | --- | --- |
| 6 | +2.91 | 5 (8) m6+m23+**m24** | 5 (30) Vas_Endo .21, **Alveolar_epi .19** |
| 11 | +2.10 | 4 (24) **m0**+m10+m23+m24 | 0 (32) **T .28**, Vas_Endo .15 |
| 13 | +2.08 | 6 (25) m23 only | 6 (37) **Plasma .27**, Vas_Endo .25 |
| 7 | +1.58 | 7 (22) **m0**+m1+m10+m23+m24 | 0 (32) **T .28** |
| 6 | +1.38 | 3 (23) m0+m1+m6+m23+**m24** | 5 (30) **Alveolar_epi .19** |
| 12 | +1.35 | 9 (47) m1+m23+**m24** | 5 (30) **Alveolar_epi .19** |
| 47 | +1.08 | 1 (64) **m1**+m23 | 3 (104) **SMC .57** |

**This independently validates the motif annotations.** m1 = SMC, m0 = T cell, m24 =
alveolar were read off each motif's LRI factor weights in `lri_motifs.csv` (2026-09-04); here
the motif clusters land in exactly the matching cell-type niche, from data the annotation
never touched — m24-carrying clusters in the alveolar niche (3 of them), m0-carrying clusters
in the T-cell niche, the m1 cluster in the SMC niche (the single largest cell of the table,
47 vessels).

Figure caveat: the colour limit is the 95th percentile of |log2 O/E| over cells with >= 5
vessels and is clipped there — the niche partition has an n=1 and an n=2 cluster whose O/E
reaches ~2^5 and would otherwise flatten everything. Marginal counts are on the tick labels.

Outputs `motif_k10_vs_niche/`: heatmap, `contingency_k10.csv`, `association_k10.csv`,
`vessel_labels_k10.csv`. `--k` compares a different motif clustering.

## Explorer: colour the vessel outlines by cluster (2026-09-08)

Skill: `interactive-spatial-plot` (tracked here). Its rule "never fork
`scripts/build_explorer.py`" is already broken in this project and predates this change —
the shipped explorer is built by the bespoke `spatial_explorer.py` (1,605 lines), not the
skill's `build_explorer.py` (1,385). Consequence, hit here: **the skill's
`verify_explorer.js` cannot boot this page** (it expects `META.samples/layers/overlays` and
`#f_sam`; this page has `META.sections/celltypes/vessels` and `#f_sec`). Left as is rather
than widened, because the two pages have genuinely different APIs, not just different ids.

**New:** a "Vessel colour" select in the sidebar — *outline only* (default, unchanged
behaviour), *Motif k-means (K=10)*, *Cell-type niche*. It applies to the blood-vessel view
and to the vessel overlay on every other view.

- Labels join on `vessel_id` against the lumens parquet, never row order. **300 of 309
  lumens carry a label** under both schemes; the 9 without are the lumens whose 30 µm ring
  holds fewer than `MIN_RING = 20` cells, so neither clustering scored them. They draw in
  `#c9c9c9` and the legend counts them.
- Colours are tab10. At K=9–10 no categorical palette is identity-safe; the mitigations are
  that a coloured lumen is **filled** (0.45 alpha for real outlines, 0.90 for the sub-6 px
  ring markers) rather than only outlined, so the hue actually reads.
- Zero size cost that matters: 21.02 MB, unchanged to 2 dp.

**`verify_explorer_luad.js` is new and tracked** — the harness for this builder. The
original was untracked and was lost with its session, exactly the risk flagged then. It
stubs DOM/Canvas/streams and drives the real page: boot and blob decode, the three
categorical planes, every view × every section through `paint()`, every legend, the 11
per-family all/none links, the render queue, lazy motif decode, the modal count, and the new
control under every scheme on both the vessel view and an overlaid view. **ALL CHECKS
PASSED** (71 putImageData, 4,944 strokes, 1,236 fills, 1,768 ring markers).

Two harness bugs found while writing it, both mine not the page's: `String.prototype.matchAll`
does not exist on iris's node 10.24, and `CT_LV` is a per-level rarity lookup rather than the
per-cell plane (`CT`).

Not verified: anything visual. The page needs opening in a real browser once.

### Row-normalised heatmaps and a Sankey (2026-09-08)

Added to `compare_motif_vs_niche_clusters.py`, same 300 vessels, same contingency table.

- `rownorm_motif_rows_k10` — rows = motif cluster, each row sums to 100 %: where that
  motif cluster's vessels land among the niches.
- `rownorm_niche_rows_k10` — rows = niche, each row sums to 100 %: the transpose question.
- Both use a **sequential single hue** (Blues, 0–100 %), not a diverging map: a row share
  is a magnitude with no meaningful midpoint. Marginal counts sit on both tick labels so a
  percentage can always be turned back into a count.
- `sankey_k10` — bipartite flow, ribbon width = vessels, ribbon colour = motif cluster.

**Read the two 100 % cells in the niche-row version as noise, not signal.** Niche 4 is
n = 1 and niche 7 is n = 2 (the airway singleton and the tumour pair from the KneeLocator
K=9), so their rows are 100 % by construction. The marginal counts on the tick labels are
there to make that immediately checkable.

Sankey layout notes: node order comes from **6 alternating barycentre sweeps** (Sugiyama) —
at 10 × 9 the numeric order is a hairball. Order is layout only; identity is the node label
and the ribbon colour, both tied to the cluster id and assigned in fixed order. Ribbons are
drawn widest-first per source so thin ones stay visible, and node labels are pushed apart to
a 0.026 minimum spacing because the n=1 and n=2 niches are slivers whose labels otherwise
overlap their neighbour's.

The Sankey makes the headline visible directly: motif cluster 1 (m1+m23, n=64) is one thick
band into niche 3 (SMC .57, n=104), 47 of 64 vessels. Everything else genuinely mixes —
which is what ARI 0.158 means.

### Each vessel classification against the pathologist's region (2026-09-08)

Two more row-normalised heatmaps in `compare_motif_vs_niche_clusters.py`: rows = vessel
class, cols = pathology region, each row summing to 100 %. A vessel has no intrinsic
region — `layer` is the **modal `tissue_layer` of its 30 µm ring cells**, as
`vessel_motif_kmeans.py` records it. Columns run normal → tumour (`PATH_ORDER`), not
alphabetically, because the levels are ordinal.

| classification | Cramér's V | χ² | df | p |
| --- | --- | --- | --- | --- |
| motif cluster (K=10) | **0.370** | 164 | 36 | 3.2e-18 |
| niche cluster | **0.322** | 125 | 32 | 7.4e-13 |

Both associate with region; the motif partition slightly more than the niche partition.
Region marginals are 159 Tumor_region / 61 Normal_distal / 59 Normal_peri / 13
Interface_normal / 8 Interface_tumor, so **53 % is the uninformative baseline for the
Tumor_region column** — a row near 53 % there says nothing.

What stands out (motif): clusters 4 and 6 are **100 % Tumor_region**, cluster 0 is 91 %;
at the other end cluster 9 is 51/40 % Normal_distal/Normal_peri and cluster 3 is 48 %
Normal_distal. (Niche): niche 6, the plasma-cell niche, is **95 % Tumor_region**; niche 3,
the SMC niche, is 67 % normal; niche 5, the alveolar niche, is 77 % normal.

**Caveat that has to travel with these two figures.** Pathology region is largely a
restatement of condition: `layer` vs `condition` Cramér's V = **0.576** (Normal_distal is
60 AIS / 1 LUAD; Tumor_region 44 AIS / 115 LUAD), and `layer` vs section 0.500. So
"cluster × region" is substantially the same fact as "cluster × section", which the earlier
audit showed the raw fraction features carry by construction. Unlike the motif-vs-niche
association, **these two were not tested against a within-section null** — do not read them
as evidence that a cluster is a tumour-region vessel type.

Rows with n = 1 or 2 (niche 4, niche 7) read 100 % by construction; marginal counts are on
the tick labels. Long region names are rotated 40°, which the shared `rownorm_heatmap`
now does whenever a column label exceeds 6 characters.

## Vessel size against both classifications (2026-09-08)

`plot_vessel_size_by_cluster.py`. Lumen area spans **147–333,400 µm², 3.4 orders of
magnitude**, so every axis is log10. 298 of 300 vessels carry an area.

**Not a section artefact.** Lumen area is only weakly section-dependent (Kruskal p = 0.019
across the four sections; medians 1,450 / 2,322 / 1,759 / 3,422 µm²), and after centring
log10 area on each section's own median the cluster effect survives:

| classification | Kruskal on log10 area | section-centred |
| --- | --- | --- |
| motif K=10 | 1.7e-11 | **7.5e-09** |
| niche K=9 | 2.8e-08 | **2.4e-07** |

Motif clusters by median lumen area — a **13× span** across clusters:

| cluster | n | median µm² | IQR | identity |
| --- | --- | --- | --- | --- |
| 1 | 64 | **8,086** | 2,056–23,205 | m1/SMC only — large arteries |
| 2 | 36 | 3,253 | 1,380–6,385 | m1+m0 |
| 6 | 25 | 2,976 | 1,296–6,871 | m23 only, plasma-cuffed |
| 5 | 8 | 2,021 | 1,115–3,478 | m24 alveolar |
| 8 | 29 | 1,851 | 790–4,247 | m0+m6+m1 |
| 9 | 47 | 1,810 | 1,200–4,900 | m1+m24, normal lung |
| 0 | 22 | 1,388 | 743–2,477 | m0 high, m1 low — T-cuffed |
| 7 | 22 | 1,229 | 687–1,929 | m10+m1 — arterialised angiogenic |
| 3 | 21 | 1,041 | 849–1,709 | m6+m0+m24 |
| 4 | 24 | **604** | 442–963 | m10 high, m1 low — neovessels |

Niche clusters span 909–4,053 µm² over the non-degenerate groups (niche 7, n=2, is 185 µm²)
— a narrower spread than the motif partition, on the same vessels.

**Mechanism, from `lumen_area_vs_motif_k10`:** m1 (SMC) is the **only** motif positively
correlated with calibre (ρ = +0.222 raw, +0.262 section-centred); every other motif is
negative, m23 most strongly (−0.372 / −0.407). So the cluster ordering is not an extra fact
— it is the SMC program tracking vessel wall calibre, which is what a mural programme should
do, and the k-means ordering falls out of it.

That also sharpens the earlier cluster 4 vs 7 contrast: both carry m10, but 4 is m1-low and
604 µm² while 7 is m1-positive and 1,229 µm² — de-novo neovessels versus angiogenesis on an
existing muscular vessel, separated by twofold calibre.

Figures order rows by median (layout only; colour is the fixed tab10 cluster identity shared
with the Sankey and contingency figures). The correlation figure's legend keys are neutral
grey — taking them from the bars would borrow m0's hue and read as a row label.

Caveats: descriptive, no multiple-comparison correction across the 10 pairwise size
contrasts; lumen area comes from the curated H&E polygons, so it inherits whatever the
curation does with obliquely-sectioned vessels; and 2 of 300 vessels have no recorded area.

## Motif loading vs distance to the curated vessels (2026-09-08)

`plot_motif_distance_profile.py`. Reuses `signed_distance` / `load_polys` from
`motif_vessel_alignment.py` rather than reimplementing the geometry. That script already
stores a profile, but only in 7 coarse bands with an open-ended ">200 µm" bucket that has no
plottable x, so the bins are recomputed at **10 µm** over −40…300 µm. Not split by condition.

**y = van der Waerden normal score** of the loading (rank → quantile → probit), per motif
**within section over all of that section's cells**, never re-ranked inside a subset — the
convention already used by `vessel_motif_cooccurrence.py`. Raw loadings cannot share an axis
(orders of magnitude apart, extremely right-skewed). 0 = section average. Sections are
averaged with **equal weight**, not pooled, so the two large LUAD sections do not set the
curve; the shaded band is ±1 SD across the four. Bins with <50 cells in a section contribute
nothing.

**The result: motif 1 is the only motif that rises at vessels.** Ranking all 25 by the
0–10 µm score minus the 150–200 µm plateau:

| | all cells | non-vascular cells only |
| --- | --- | --- |
| **m1 SMC** | **+1.135** (rank 1, the only positive in the top 8) | **+0.723** (still the only positive in the top 8) |
| next 7 | m14 −0.911, m3 −0.885, m19 −0.883, m21 −0.865, m7 −0.833, m15 −0.832, m8 −0.825 | m7 −0.779, m8 −0.753, m20 −0.744, m19 −0.742, m4 −0.729, m14 −0.727, m21 −0.723 |

m1 peaks at +1.14 in the first 10 µm and decays to the section average by ~150 µm — a
~100 µm length scale, consistent with a mural programme on the vessel wall.

**The `--nonvasc` control matters and is shipped alongside.** Endothelium, SMC, pericytes
and lymphatic endothelium occupy the wall by construction, so on all cells any motif carried
by another cell type looks depleted at the lumen for purely compositional reasons. Averaging
over non-vascular cells only (1,399,750 of 1,676,162):

- m1 survives, weaker (+1.135 → +0.723) — it is not only the SMCs.
- m23 (endothelial) **does not**: its peak drops from +0.53 to a spike confined to the first
  ~10 µm, and it crosses below zero by ~50 µm. On non-vascular cells it is close to a null
  motif. Its ±1 SD band is the widest of the six in both versions — the least consistent
  across sections.
- m6, m10, m24 stay depleted at the lumen but much less so; m10 remains the most depleted.
- m0 (T cell) is non-monotone in both: a dip right at the boundary, then a plateau slightly
  **above** the section average from ~50 µm out.

### Caveats

1. **The far plateau is 0 by construction.** The normal score is centred within section and
   most cells are far from a lumen, so every curve must return to ~0. Only the near-vessel
   deviation carries information; the plateau is not a measurement.
2. **The annotation is sparse and biased to large vessels** (71–83 lumens per section), so
   "far from a curated lumen" is not "avascular" — a real capillary bed sits in the far bins.
   This is the same caveat `motif_vessel_alignment.py`'s docstring carries.
3. Per-cell means over 1.4–1.7 M spatially autocorrelated cells: the curves are descriptive,
   and no p-value is attached to them here. The vessel-level version of this question, with
   each ring against its own surround, is `vessel_motif_local.py`.

Outputs `motif_distance_profile{,_nonvasc}/`: the two figures, `profile_mean.csv`,
`profile_sd.csv`, `distance_dependence.csv`, `bin_cell_counts.csv`.

## figures_interpretation/ — supporting plots for the biological read (2026-09-09)

**Corrected the same day.** The first build interpreted the **motifs** (m1 = SMC, m24 =
alveolar …). The ask was to support the interpretation of the **vessel clusters** — "what is
cluster 0" — which is what the text answer of 2026-09-08 had given. Rebuilt around the
clusters; the motif figures are kept as `motif_background/` because they are the key needed
to read the cluster figures, not the answer to the question.

### `clusters/` — the answer

`motif_cluster_profiles` and `niche_cluster_profiles` are the main figures: one row per
cluster (ordered by median lumen area, n printed), columns in three blocks — **signalling**
(the six motifs, section-centred), **ring cell type**, **size / context** (log10 lumen area,
ring cells, nn10 density, % Tumor_region, % LUAD). Every cell prints its raw value; colour is
that column's z across clusters, purely so one scale carries fractions, calibre and
percentages together. Numbers in the matching CSVs.

Reading a row gives the vessel type directly, e.g.

| cluster | n | reading |
| --- | --- | --- |
| 0 | 22 | m0 +0.36, m1 −0.52; SMC 0.01, T 0.25; 1,477 µm²; 91 % Tumor_region → T-cell-cuffed non-muscular tumour vessels |
| 1 | 64 | m1 +0.18 (raw 0.99), rest down; SMC 0.53; **7,687 µm²** → large arteries, signalling-pure |
| 4 | 24 | m10 +0.42, m1 −0.54; SMC 0.00; **699 µm²**, densest tissue; 100 % Tumor_region → de-novo neovessels |
| 7 | 22 | m10 +0.28 **and** m1 +0.21; SMC 0.27 → angiogenesis on a muscularised vessel |
| 6 | 25 | every motif ≤ 0; Plasma 0.20; sparsest tissue; 100 % LUAD → plasma-cuffed, no programme on |

`niche_cluster_profiles` makes the earlier text observation visible: **niche 2 (n=74) has an
all-white signalling block** — a quarter of the vessels form a cell-type group with no
signalling identity.

Also in `clusters/`: `cluster_identity_matrix` (the signalling block alone, larger),
`neovessel_vs_arterialised` (the two clusters niche 0 merges), and
`examples_motif_clusters/` (three example vessels per cluster, copied from the K=10 sweep).

`CAPTIONS.md` carries the per-figure claim and caveat, and closes with what the folder does
**not** show — no condition claim, and no claim that the motif classification is better.

---

### First build (motif-level, retained as `motif_background/`)

`build_interpretation_figures.py`. One script for all six, so style, palettes and the saver
are shared; `CAPTIONS.md` in the output directory carries per-figure claim + caveat.

| figure | claim it supports |
| --- | --- |
| `motif_identity_celltype` | the annotations are read off the LRI factor weights, not assumed |
| `motif_identity_top_lri` | the same at gene level (`JAG1→NOTCH3`, `VEGFA→KDR`, `PECAM1`, `SFTPD`) |
| `perivascular_enrichment` | **m1 +0.491 (p=2e-38) and m23 +0.245 (p=3e-40) are the only two of 25 enriched at vessels**; each vessel is its own control |
| `m1_m23_distance_control` | drop the wall's own cells and m1 survives (+1.135 → +0.723), m23 collapses to a ~10 µm spike |
| `cluster_identity_matrix` | each vessel cluster is a distinct section-centred signalling identity, ordered by calibre |
| `neovessel_vs_arterialised` | the worked case: cluster 4 vs 7 both sit in niche 0, motif separates them |

New number computed for the caption of figure 4, and it is the justification for the whole
`--nonvasc` control: **83.2 % of cells in the 0–10 µm bin are vessel-wall cells** (Vas_Endo,
Lym_Endo, SMC, Pericyte), against 17.2 % at 200–300 µm and a 16.5 % tissue baseline. Distance
to a vessel is very nearly the same variable as cell type, which is why an uncontrolled
"loading vs distance" curve cannot be read as a spatial gradient.

The ring-vs-own-surround numbers here (m1 +0.491, m23 +0.245) are on **positive fractions
over 300 vessels with a clean 100–250 µm surround**, slightly different from the earlier
audit's +0.392 / +0.227 because that run used a different surround-exclusion pass; both give
the same two motifs and the same ordering.

Three figures needed layout fixes found by rendering and looking: the rotated per-motif tags
in figure 2 landed on top of the pair labels (replaced with a legend + block hairlines), and
figure 5's x labels collided (rotated 35°).

`CAPTIONS.md` closes with what these figures do **not** show — no condition claim (the m24
AIS/LUAD result stays withdrawn), and no claim that the motif classification is better than
the niche one, which still needs the held-out tests proposed above.

### Clustered versions of the two row-normalised heatmaps (2026-09-09)

`clustermap_motif_rows_k10` and `clustermap_niche_rows_k10`, added to
`compare_motif_vs_niche_clusters.py`. The un-clustered `rownorm_*` versions are **kept** —
they preserve numeric cluster order, which is what you want when checking a specific id.

Method: each margin is clustered on **its own** profile normalisation — rows on the row-%
matrix, columns on the column-% matrix — so neither dendrogram is driven by how big the
other margin's groups happen to be. Average linkage, Euclidean. The heatmap displayed is
always the row-% one, identical to the un-clustered figure. Row labels sit on the right
because a left dendrogram occupies exactly where they would otherwise be.

**Structure the clustering exposes (motif rows):**

| block | motif clusters | what they share |
| --- | --- | --- |
| muscular | 8, 2, 9, 1 | concentrate in niche 3 (SMC .57) and niche 2 (Fibro) — cluster 1 is 73 % niche 3 |
| plasma / fibrous | 6, 0 | niche 6 (Plasma .27) and niche 2 — cluster 6 is 52 % niche 6 |
| immune | 7, 3, 4 | spread over niches 6, 0, 8 — the T-cell-cuffed side |
| alveolar | 5 | alone: 75 % in niche 5 (Alveolar_epi .19) |

Niche columns cluster as {2, 3} (the two large composition buckets), {6, 0, 8, 1} (the
immune-ish niches), {5} (alveolar) and {4, 7} on their own.

**Caveat:** niches 4 (n=1) and 7 (n=2) branch off first, which is their size, not their
biology — a row that is 100 % one column is maximally far from everything else by
construction. Read the dendrogram over the non-degenerate groups only.

The pathology heatmaps were deliberately **not** clusterised: their columns are ordinal
(normal → interface → tumour) and reordering them would destroy that.

## vessel_explorer_AIS_LUAD.html — the cluster-driven vessel browser (2026-09-09)

`vessel_explorer.py` + `verify_vessel_explorer.js`, both tracked. **7.69 MB**, offline,
double-click.

Skill: `interactive-spatial-plot` (tracked here). This is a **new tool, not a fork** — the
skill's `build_explorer.py` builds a section-panel explorer from a config, and a
cluster-driven crop browser is not expressible in that schema. It is also not a copy of
`spatial_explorer.py`: it imports that script's helpers (`gzb64`, `load_vessels`,
`load_vessel_clusters`, `TAB20`, `TISSUE_COLORS`, `read_cat`, `QMAX`) and ships its own page.
Neither existing builder was modified.

**What the page does.** Pick a classification (motif k-means K=10 / cell-type niche K=9) and
one cluster; you get a per-section panorama with that cluster's lumens picked out in the
cluster colour (all other lumens thin grey), and a grid of zoom-in crops, one per vessel in
the cluster. Crop size 700 (default) / 500 / 350 / 250 µm. Cells recolour by cell type,
pathology layer, motif loading (25 motifs, any of 6 ramps) or motif state. First 24 crops,
with a "show all" link.

**Size lever:** only cells within **520 µm** of a curated lumen are shipped — 650,940 of
1,676,162 (39 %). 520 µm is what a 700 µm crop centred on any lumen can reach
(700/2·√2 = 495). The consequence is stated on the page: the panorama shows perivascular
tissue, not the whole section. `--ctx` changes it.

**Two bugs the harness caught**, neither visible without it:

1. The panorama cell layer is cached so that clicking through clusters does not re-stamp
   650 k points. The first cache was a flat map cleared on every call, so **each section
   evicted the previous one and nothing was ever reused** — total fillRect over the test
   sweep was 44 M. Keyed per section it is 6.2 M, and a warm panorama re-stamps 8 instead of
   650,952.
2. My own first version of the "canvas actually drew" check counted `fill()` (polygon fills)
   and not `fillRect()` (cell stamps), so a page that drew **no cells at all** would have
   passed. The stub now counts `fillRect` and asserts >100 k, plus a cold/warm/invalidate
   triple on the cache.

Harness coverage: boot and blob decode, the three categorical planes, section ranges
contiguous and y-sorted (the crop lookup binary-searches on that), a label per lumen in
range for both classifications matching META's counts, **every classification × every
cluster × every colour layer rendering without throwing** (19 clusters × 4 layers), the
y-band lookup against a brute scan, every crop size, and the cache. **ALL CHECKS PASSED.**

Not verified: anything visual. Open it in a browser once.

Not shipped, to keep the file at 7.7 MB: `cell_state` (56 levels). Say the word and it goes
in at roughly +1 MB.

### Linear-axis versions of the vessel-size box plots (2026-09-09)

`plot_vessel_size_by_cluster.py` now emits three versions of each size figure; the log one
is unchanged and still the default read.

| suffix | axis | note |
| --- | --- | --- |
| *(none)* | log10 | as before |
| `_linear` | raw µm², full range 0–333,400 | honest, and mostly unreadable: cluster 1's tail pushes every other cluster into the left edge |
| `_linear_x25000` | raw µm², clipped at `--xmax` | **the usable linear one**; "21 vessel(s) beyond the axis" is printed on the figure |

**The box statistics are not the same between the two scales.** On the log figure the
quartiles and whiskers are computed on log10 area; on the linear ones they are computed on
the raw values. They are two different summaries of the same data, not a relabelling of one
axis — a 1.5×IQR whisker in log space is a multiplicative span, in linear space an additive
one. Row order is identical in both because the median is order-preserving under log.

Clipping never drops a vessel silently: `--xmax` counts what falls beyond and prints it.
Without `--xmax` only the two unclipped versions are written.

Fixed while doing this: the linear branch used `WIDE`, which this script had never defined
(only `SINGLE`) — `NameError` on the second figure.

### Globally-centred versions of the cluster figures (2026-09-09)

Asked for after questioning why the motif columns were centred on each vessel's own section
rather than on all vessels. `build_interpretation_figures.py` now takes `centre` and emits
both; the section version keeps its filenames, the other gains `_globalcentre`:
`cluster_identity_matrix_globalcentre`, `motif_cluster_profiles_globalcentre`,
`niche_cluster_profiles_globalcentre` (+ CSVs).

Numbers that justify keeping section-centring as the default — between-section share of the
per-vessel positive-fraction variance: **m24 48 %, m10 30 %**, m1 14 %, m6 12 %, m0 7 %,
m23 7 %. And under global centring a cluster's m24 value is essentially a restatement of its
section composition (**R = 0.976** predicting it from section composition alone).

Mean |difference| between the two centrings, per motif, tracks that exactly:

| motif | mean \|change\| | between-section variance share |
| --- | --- | --- |
| m24 | **0.180** | 48 % |
| m10 | **0.085** | 30 % |
| m6 | 0.053 | 12 % |
| m1 | 0.051 | 14 % |
| m0 | 0.033 | 7 % |
| m23 | 0.014 | 7 % |

So the global version inflates precisely the two motifs whose section offsets the audit
already showed to be whole-tissue effects: cluster 4's m10 +0.42 → +0.62, its m24
+0.23 → +0.42; cluster 5's m24 +0.27 → +0.52. Cluster ordering barely moves (m24 column,
Spearman 0.96 between the two).

The honest bottom line, now written into `CAPTIONS.md`: with 4 sections and clusters that are
70–100 % one section, section and cluster are partly collinear — centring moves the confound,
it does not remove it.

### The z composite retired; three raw-value blocks instead (2026-09-09)

Asked for after questioning what an across-cluster z does to interpretation.
`motif_cluster_profiles` / `niche_cluster_profiles` (one wide z-scored composite) are moved
to `clusters/superseded/`, along with `cluster_identity_matrix{,_globalcentre}` which the new
signalling figure duplicates. Nothing deleted.

Replaced by `cluster_blocks()`, three figures per classification, **coloured by the raw
value**:

| figure | colour |
| --- | --- |
| `<cls>_signalling` (+`_globalcentre`) | the centred positive fraction itself, one symmetric scale across all six motifs |
| `<cls>_celltype` | the mean per-vessel ring fraction itself, one 0→max scale across all seven types |
| `<cls>_context` | mixed units, so colour is each column's own min→max position among clusters — a shading, not a magnitude; the colourbar says so |

**The three ways the z misled**, all measured:

1. Its denominator is the SD of **cluster means**, which is 0.25–0.98 of the **vessel** SD
   column by column (`ring cells` 0.25, `log10 area` 0.44, `frac_VasEndo` 0.58, m24 0.64,
   m1 0.98). z = +2 on `ring cells` was half a vessel SD, not two.
2. Its zero is the **unweighted** centre of the 10 cluster means, so n=8 and n=64 count
   equally. For `frac_SMC` that centre sits **0.42 cluster-SD below** the true vessel mean —
   a cell shown as neutral white was in fact below-average.
3. Per-column standardisation made colour incomparable **between** columns: `frac_VasEndo`
   spans 0.21 absolute but got a **wider** z range (3.64) than m1 SMC, which spans 0.75
   (2.38). The figure implied VasEndo varied more than the SMC motif; the opposite is true.

Fixed while doing this: the context block's header said `log10 area` while the cells print
raw µm².

Colour change (2026-09-09): the two sequential blocks, `<cls>_celltype` and `<cls>_context`,
use **Reds** rather than Blues on request. Both are single-hue light→dark, which is the
sequential rule; the diverging `<cls>_signalling` keeps RdBu_r with its neutral midpoint.

Row order flipped (2026-09-09): the size-ordered **heatmaps** now put the **largest lumen at
the top**. `imshow` draws row 0 at the top, so the sort is now descending by median area;
the y-axis label says so. The box plots in `plot_vessel_size_by_cluster.py` were already
largest-at-top (matplotlib's y axis increases upward, so an ascending sort puts the largest
last, i.e. highest) and were not touched. `cluster_identity()` is no longer called from
`main()` — `<cls>_signalling` supersedes it — so the superseded figures are not regenerated.

Colour change (2026-09-09, cont.): the row-% panels in `compare_motif_vs_niche_clusters.py`
— `rownorm_motif_rows`, `rownorm_niche_rows`, `rownorm_motif_vs_pathology`,
`rownorm_niche_vs_pathology` and both `clustermap_*` — switched from Blues to **Reds**, to
match the `<cls>_celltype` / `<cls>_context` blocks. All of them show the same kind of
quantity (a row share, 0–100 %), so they share one single hue. `motif_k10_vs_niche_heatmap`
keeps RdBu_r: log2 observed/expected is signed and needs a neutral midpoint.

## Per-motif: positive fraction vs vessel size, and vs distance to the vessel (2026-09-10)

`plot_motif_size_distance.py`. Back off the clusters, one figure per motif, 25 motifs, two
families. `figures_motif_profiles/`.

**`size_vs_positive/m<k>`** — one point per vessel (298 with a recorded lumen area): x = log
lumen area, y = fraction of that vessel's 30 µm ring cells positive. Colour **and** marker
encode section. LOWESS trend. Spearman reported overall and per section, because section is
a known confound.

**`distance_vs_positive/m<k>`** — one point per **cell**: x = distance to the nearest lumen
boundary (0–300 µm), y = the GMM state. Smoothed with a **binomial GLM on a natural cubic
spline basis (df = 6)**, fitted per section. Cells are binned at 2 µm and the GLM is fitted
to the per-bin (positive, total) counts — binomial sufficiency makes that the same likelihood
as one row per cell, and it turns 100 fits on ~85 k rows into 100 fits on ~150 rows. Four
thin section curves + a thick equal-weight mean + a **dashed non-vascular-only mean**.

### Both figures agree, and they agree with the ring-vs-surround test

Size (Spearman with log lumen area): **m1 is the only motif that goes UP with calibre**
(+0.22; per section +0.55, +0.33, −0.19, +0.20). All 24 others are negative or flat, most
negative m16 −0.47, m9 −0.38, m0 −0.37, m23 −0.37.

Distance, on **non-vascular cells only**, P(positive) at 10 µm minus at 150 µm:

| motif | Δ non-vasc | P at 10 µm | size ρ |
| --- | --- | --- | --- |
| **m1 SMC** | **+0.406** | 0.70 | +0.22 |
| **m23 endothelial** | **+0.285** | 0.94 | −0.37 |
| m17 | +0.086 | 0.58 | −0.19 |
| m16 | +0.075 | 0.77 | −0.47 |
| m18 | +0.073 | 0.75 | −0.22 |
| … | | | |
| m3 | −0.244 | 0.14 | −0.28 |

Only m1 and m23 rise appreciably toward the vessel once the wall's own cells are dropped;
everything else is flat or falls. That is the same answer the paired ring-vs-own-surround
test gave (m1 +0.491, m23 +0.245, everything else ≤ 0.035), reached by a different route —
per-cell smoothing over a continuous distance instead of a per-vessel paired contrast.

**m23 behaves differently between the two smoothings, and the reason matters.** On the
*normal score of the loading* (`figures_distance_profile/`) m23's non-vascular peak collapsed
to a ~10 µm spike; on *P(GMM positive)* here its non-vascular Δ is +0.285, larger than the
all-cells Δ (+0.226). These are different quantities — a rank-based mean of a continuous
loading versus the probability of crossing the GMM threshold — and m23 is saturated
(P = 0.94 at 10 µm, 0.68 at 150 µm), so the loading can fall a long way while the state stays
on. Read m23 as "almost every cell near a vessel is m23-positive", not as a strong gradient.

Caveats: the per-section curves are per-cell fits over spatially autocorrelated cells, so no
p-value is attached; the annotation is sparse and biased to large vessels, so the far bins
are not "avascular"; and the 0–10 µm bin is 83 % vessel-wall cells, which is exactly why the
dashed line is on every figure.

Dashed line removed (2026-09-10): `distance_vs_positive/m<k>` now shows only the four
section curves and their equal-weight mean. The non-vascular-only fit is **still computed
and still in `distance_vs_positive.csv`** (`p_at_10um_nonvasc`, `delta_nonvasc`) — it is only
no longer drawn. That column is what separates "a real gradient in the surrounding tissue"
from "the curve is tracking cell-type composition", since 83 % of cells 0–10 µm from a lumen
are vessel-wall cells against a 16.5 % tissue baseline. Reading these figures without
checking it will over-read every motif that rises toward the lumen.

Confidence areas added (2026-09-10). `distance_vs_positive/m<k>` now carries **two** bands,
and they are not the same quantity — the figure says so under the title:

- **coloured band** = that section's spline fit 95 % pointwise CI, from
  `GLMResults.get_prediction().summary_frame()` (statsmodels applies the inverse link, so it
  is on the probability scale). These come out as **hairlines**: each section contributes
  70–110 k cells, and the CI does not account for spatial autocorrelation between
  neighbouring cells, so it understates the real uncertainty. A narrow band here is not
  evidence of a reliable curve.
- **grey band on the mean** = ±1 SD **between the four sections**. This is far wider and is
  the honest uncertainty for any statement about "the" profile — on m10, for instance, the
  four sections span 0.06 to 0.40 at the boundary, so the mean is nearly meaningless there
  while each section's own curve is tight.

That contrast is the point of showing both: the fit is precise, the biology is not
reproducible across sections for most motifs. m1 is the exception where all four sections
move the same way.

### Pooled six-motif distance figures (2026-09-11)

`plot_motif_size_distance.py --pooled` → `distance_vs_positive_pooled/m{0,1,6,10,23,24}`,
six figures, alongside the 25 per-section ones (which are unchanged).

Differences from `distance_vs_positive/`:

- **only the six Vas_Endo-related motifs**, not all 25;
- **all four sections pooled into one curve** instead of four + a mean;
- the **binned observations are drawn** — open circles, 10 µm bins, area ∝ √n. The FIT still
  uses the 2 µm grid; the dots are a coarser read of the same data, there so the smoother
  can be checked against what it is smoothing;
- **every plotted observation gets a graduation on the x axis** (minor ticks at the dot
  positions, major every 50 µm);
- the 95 % CI band is still drawn but is a hairline at n = 341,289.

| motif | P at 10 µm | P at 150 µm | Δ |
| --- | --- | --- | --- |
| **m1 SMC** | 0.819 | 0.362 | **+0.458** |
| **m23 endothelial** | 0.912 | 0.670 | **+0.243** |
| m24 alveolar | 0.384 | 0.440 | −0.056 |
| m0 T cell | 0.413 | 0.547 | −0.134 |
| m10 tumour-endo | 0.188 | 0.366 | −0.178 |
| m6 macrophage | 0.269 | 0.479 | −0.210 |

**Pooling weights sections by cell count** — P17_LUAD contributes 109,749 of the 341,289
cells, so the pooled curve leans on it. The per-section figures are what show the
disagreement between sections, and for most motifs that disagreement is larger than the
gradient; m1 is the one where all four sections move the same way. Nothing here is corrected
for cell-type composition either: 83 % of cells 0–10 µm from a lumen are vessel-wall cells,
and `distance_vs_positive.csv`'s `delta_nonvasc` column is still the check for that.

**Conventions note.** These were written under the old layout
(`/data1/tanseyw/projects/fanj2/alarmist_luad/` + the `download_20260902/` bundle) and with
`/data1/tanseyw/fanj2/envs/comp-liana/bin/python`, to stay byte-consistent with their 50
siblings. The updated CLAUDE.md asks for `scripts/_common/paths.py` (`results()` /
`deliverables()`), says not to keep a second downloadable copy, and pins ALARMIST work to
`/home/fanj2/.conda/envs/spatial/bin/python`. Migrating this whole output tree is a separate
decision — flagged, not done.

## Per-vessel ring positivity, per motif — the distributions before any cutoff (2026-09-11)

`plot_vessel_positivity_hist.py` → `figures_vessel_positivity/m{k}_vessel_positive_hist`,
six histograms over the 300 curated vessels, x = fraction of that vessel's 30 µm ring
positive. Bars stacked by section; dashed lines mark each section's **whole-section
cell-level** positive rate, i.e. what a ring would score if it were a random patch of that
tissue. Numbers in `vessel_positivity.csv`.

| motif | median | q75 | q90 | =0 | =1 | ≥0.5 | ≥0.75 | ≥0.9 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| m0 T cell | 0.433 | 0.762 | 0.989 | 8.7 % | 10.0 % | 138 | 79 | 51 |
| **m1 SMC** | **0.978** | 1.000 | 1.000 | 7.0 % | **48.7 %** | 230 | 200 | 175 |
| m6 macrophage | 0.250 | 0.502 | 0.715 | 11.7 % | 1.7 % | 76 | 21 | **9** |
| **m10 tumour-endo** | **0.051** | 0.370 | 0.700 | **42.7 %** | 3.7 % | 52 | 26 | 18 |
| **m23 endothelial** | **1.000** | 1.000 | 1.000 | 0.3 % | **55.3 %** | 290 | 269 | 228 |
| m24 alveolar | 0.358 | 0.807 | 1.000 | 19.3 % | 13.0 % | 129 | 89 | 60 |

Whole-section cell-level baselines (the random-patch score):

| motif | P17_AIS | P17_LUAD | P21_AIS | P21_LUAD |
| --- | --- | --- | --- | --- |
| m0 | 0.527 | 0.692 | 0.416 | 0.278 |
| m1 | 0.293 | 0.403 | 0.419 | 0.271 |
| m6 | 0.493 | 0.428 | 0.681 | 0.374 |
| m10 | 0.610 | 0.406 | 0.324 | 0.424 |
| m23 | 0.731 | 0.632 | 0.780 | 0.527 |
| m24 | **0.819** | **0.120** | **0.877** | **0.243** |

### What this means for defining a "positive vessel"

**m1 and m23 cannot carry a specificity definition.** Half the vessels sit at *exactly* 1.0
(m1 48.7 %, m23 55.3 %); 230 and 290 of 300 are ≥ 0.5. A cutoff anywhere selects most of the
cohort, and the ties at 1.0 mean it cannot even rank them. Their signal is real — both are
the two motifs enriched over a vessel's own surround — but it is not *discriminating between
vessels*.

**m10 is the opposite and the most usable.** 42.7 % of rings are exactly 0 and the median is
0.051, so a vessel at ≥ 0.75 (26 of 300) is genuinely unusual.

**m6, m0, m24 sit in between**, and m24 is bimodal with the two modes split by section, which
is the m24 section effect already on record — its baselines differ 7× between sections
(0.877 vs 0.120), so a global m24 cutoff selects AIS vessels almost by construction.

**A global cutoff is not comparable across sections for any of these.** The baselines move
by 2.5× (m0) to 7× (m24). If the goal is a specific *vessel*, the cutoff should be relative
to that vessel's own section — or better, to its own 100–250 µm surround, which is the
design `vessel_motif_local.py` already implements and which removes regional variation too.
Say which of the three definitions you want and I will apply it:

1. absolute (`fraction ≥ c`) — simple, but not comparable between sections;
2. section-relative (`fraction − section background ≥ c`);
3. local (`ring − its own 100–250 µm surround ≥ c`) — strictest, and the one whose null we
   have already characterised.

## The cleanest single-motif vessels (2026-09-11)

`find_motif_specific_vessels.py` → `figures_motif_specific_vessels/`: a table, an example
table, and one crop figure per motif (rows = the six motifs + cell type, columns = the top 3
vessels, target row labelled in red, 350 µm fields).

### The naive definition returns almost nothing, and that is the finding

"f_k ≥ 0.75 and every other motif ≤ 0.25" gives **0 vessels for m0, m6, m10 and m24**, and
only 2 each for m1 and m23. The reason is in `figures_vessel_positivity/`: m1 and m23 are
positive nearly everywhere (median ring fraction 0.978 and 1.000; 48.7 % and 55.3 % of
vessels exactly 1.0). Taking `f_k − max(other)` on raw fractions, **m6's best margin over
all 300 vessels is +0.00 and m10's is +0.07** — those motifs are never a vessel's top motif
by any margin on the raw scale.

### The definition actually used

Two corrections, both forced by the data:

1. **m23 is excluded from the competitor set** for the other five motifs. It is an identity
   marker, not a discriminator — "clean" can only mean "clean apart from m23". m23's own
   score is still against all five others.
2. **Each motif is converted to a within-section percentile rank first**, because the
   marginals differ wildly (m23 median 1.00 vs m10 median 0.05) and so do the section
   baselines (m24: 0.877 in P21_AIS, 0.120 in P17_LUAD). Then

   `margin_k(v) = rank_k(v) − max_{j ≠ k, j ≠ 23} rank_j(v)`

| motif | margin ≥ 0.20 | is the top motif | strict absolute | best vessel | its fractions |
| --- | --- | --- | --- | --- | --- |
| m0 T cell | **11** | 50 | 0 | P17_LUAD__10 | m0 1.00, others(−m23) 0.00, m23 0.86 |
| m1 SMC | **50** | 83 | 2 | P17_AIS__25 | m1 1.00, others 0.08, m23 0.62 |
| m6 macrophage | **8** | 42 | 0 | P21_LUAD__61 | m6 0.70, others 0.38, m23 1.00 |
| m10 tumour-endo | **18** | 57 | 0 | P21_LUAD__54 | m10 1.00, others 0.04, m23 1.00 |
| m23 endothelial | **3** | 35 | 2 | P21_LUAD__52 | m23 1.00, others 0.19 |
| m24 alveolar | **10** | 53 | 0 | P21_AIS__6 | m24 1.00, others 0.86, m23 0.95 |

The crops confirm the good ones visually: for P21_LUAD__54 the m10 row is dark around the
lumen while m0/m1/m6/m24 are almost entirely grey; same for P17_LUAD__10 on m0.

### What is weak, stated plainly

- **m6 and m24 have no clean example.** m6's best still has m1 = 0.38 and m23 = 1.00;
  m24's best has m1 = 0.86. Their "specificity" is a rank effect, not a clean single-motif
  vessel, and the crops show it.
- **m23 has only 3 vessels at margin ≥ 0.20** — because it is near-saturated, it can rarely
  out-rank the others by much even when it is 1.00.
- The third example for a motif with few qualifiers is genuinely weak (m10's third has
  m10 = 0.23). That is the tail of the list, not a curation error.
- Counts are sensitive to the margin cutoff: at 0.30 they fall to 2 / 30 / 2 / 4 / 2 / 1.
  `--margin` changes it; `motif_specific_vessels.csv` has the counts.

### Crop figures reworked: state not loading, field scaled to the vessel (2026-09-14)

Three changes to `find_motif_specific_vessels.py`. **The vessel selection, the margins and
all counts are unchanged** — only the drawing.

1. **The motif rows now show the GMM ON/OFF state**, not the loading: positive cells solid
   `#3b0f70`, negative `#e2e2e2`, positives drawn on top. The 15-level magma percentile ramp
   is gone. This matches what the selection is actually made on (the per-vessel score is a
   positive *fraction*), so the picture and the number are now the same quantity.
2. **The field of view is per vessel**: `clip(lumen extent + 2 × 150 µm, 350, 1200)`. The
   examples span a 19–470 µm lumen extent, so a fixed 350 µm box cropped the two largest
   (m1's `P17_AIS__25`, 470 µm) while wasting the frame on the smallest. Fields now range
   350–770 µm. Because panels differ in scale *within* a figure, **every column carries its
   own scale bar** (length picked as the largest round number under a third of the field)
   and prints its field size in the title. Dot size scales as `√(350/fov)` so density stays
   comparable.
3. **The cell-type row already existed**; it now has a **legend** — without one, 19 tab20
   colours were unreadable. Listed types are those ≥ 2 % of the cells drawn on that figure,
   with their share; the rest are summarised as a count with no swatch.

The "no swatch" is deliberate. The first attempt recoloured sub-2 % types to `(.80,.80,.80)`,
which is indistinguishable from tab20's own `Tumor_epi` grey `#c7c7c7` **in the same row** —
the same collision already fixed once in `cluster_composition`. Every cell now keeps its
tab20 colour and the *legend* is trimmed instead of the data.

Reading the results: `P21_LUAD__54` (m10) has a dense positive ring with m1/m6/m24 almost
entirely grey, and its cell-type row is 51 % `Tumor_epi` — consistent with the
tumour-endothelial annotation. `P17_AIS__25` (m1) shows positives hugging the lumen wall
across a 770 µm field. The weak tail is still visible and still honest: m10's third example
has m10 = 0.23, m6's best still has m23 = 1.00.

### Pooled figures: raw per-cell points, not binned means (2026-09-14)

I had misread "points" as the binned observations. Corrected: `distance_vs_positive_pooled/`
now shows **one point per cell** — y = its GMM state (0 or 1), y-jittered by ±0.045 — plus a
**rug of the same points' x-projections** on the axis. The 10 µm binned means are gone; they
sat on the fitted curve by construction and showed nothing the curve did not.

Two rendering decisions, neither of which is on the figure:

- **The cloud is a seeded random subsample of 20,000** of the 341,289 cells in range
  (`--n-points` via `NPOINTS`). All 341 k render as two solid bars and a ~50 MB vector file.
  **The fit always uses every cell.** `n_drawn` is a column in
  `distance_vs_positive_pooled.csv` so the subsample is on record.
- The cloud and rug are `rasterized=True`, so the SVG embeds a bitmap for the points while
  axis text stays vector (CLAUDE.md § Plotting requires editable vector *text*).
- The rug is at alpha 0.025. At 20 k marks across 300 µm any opaque rug is a solid black
  bar; this low it reads as a density instead — visibly darker to the right, i.e. most cells
  are far from a curated lumen.

**Figure style.** `nature_publication_figures` is **not in this checkout** (CLAUDE.md
§ Skills), so CLAUDE.md § Plotting is what was followed. The descriptive annotations were
removed from all three families — the pooled one ("dots: observed fraction in 10 µm bins …"),
the per-section one ("coloured band: 95 % CI … grey band: spread between sections") and the
size one (the Spearman line). That text belongs in a caption, not on the axes; it is in this
file and in the CSVs. The per-section and size figures were re-rendered too, so the whole
directory is consistent — say if you wanted the Spearman value kept on the size figures.

Reading them: m1's y = 1 band is dense at the boundary and thins with distance; m10's does
the opposite; the curve is the same fit as before, unchanged.

### Per-section figure text trimmed (2026-09-14)

Two removals from `distance_vs_positive/`, and the same title change applied to
`size_vs_positive/` for consistency:

- title `m1 SMC · distance to vessel` → `m1 SMC`. The x-axis label already states the
  variable; the title now carries only the motif identity, matching the pooled figures.
- legend entry `mean of sections (±1 SD)` → `mean ± 1 SD`.

What remains on the distance figure is the minimum needed to read it: the motif identity,
both axis labels, the four section names, and the meaning of the black line and its band.
The band definitions (per-section colour band = that section's 95 % fit CI; grey band =
between-section SD) are in this file and were removed from the axes on 2026-09-11.

All 56 figures in the directory were re-rendered. No data, fit, or reported value changed.

### size_kde/ — lumen area of positive vs negative vessels (2026-09-14)

A third family in `figures_motif_profiles/`: `size_kde/m<k>`, 25 figures. Vessels are split
at a ring positive fraction of **0.5**, and the density of log lumen area is drawn for each
group on the same x axis as `size_vs_positive/`. Rugs show every vessel. A group below
`KDE_MIN = 15` gets the rug only — a Gaussian KDE over ten points draws a shape the data
does not support, and m23's negative group is exactly that (n = 10 of 298). Colours
`#b2182b` / `#2166ac`, validated PASS with `--pairs all`. Numbers in `size_kde.csv`.

AUROC is P(a positive vessel is larger than a negative one); 0.5 is no difference.

| motif | n pos / neg | median area pos | median area neg | AUROC | p |
| --- | --- | --- | --- | --- | --- |
| **m1 SMC** | 228 / 70 | **2,375** | 1,188 | **0.672** | 2e-05 |
| m17 | 144 / 154 | 1,769 | 2,338 | 0.403 | 0.004 |
| m18 | 212 / 86 | 1,812 | 3,672 | 0.394 | 0.004 |
| m6 macrophage | 74 / 224 | 1,423 | 2,365 | 0.354 | 3e-05 |
| m0 T cell | 136 / 162 | 1,458 | 3,395 | 0.332 | 7e-07 |
| m23 endothelial | 288 / 10 | 1,973 | 18,123 | 0.314 | 0.046 |
| m24 alveolar | 127 / 171 | 1,373 | 3,368 | 0.306 | 1e-08 |
| m10 tumour-endo | 52 / 246 | 891 | 2,469 | 0.251 | 1e-08 |
| m13 | 10 / 288 | 510 | 2,105 | 0.126 | 4e-05 |

**m1 is the only motif whose positive vessels are larger**; for 22 of 25 the positive group
is smaller, several strongly so. This is the same result the Spearman correlations gave
(m1 the only positive ρ with calibre, +0.22) reached by a coarser route, so the thresholded
view adds clarity rather than a new claim.

Two limits on how far this can be read:

- **The 0.5 threshold is not comparable between motifs.** From
  `figures_vessel_positivity/`, the per-vessel median ranges from 0.05 (m10) to 1.00 (m23),
  so 0.5 selects the top ~17 % of vessels for m10 and 97 % for m23. "Positive" therefore
  means something different in each panel, and the group sizes above show it.
- p-values are over vessels, which are spatially autocorrelated, and the whole-section
  baselines differ up to sevenfold between sections. Read the effect size and the group
  sizes, not the exponent.

A defect worth recording: the first attempt anchored the insertion on `def size_vs_motif`,
a function that lives in `plot_vessel_size_by_cluster.py` and not in this file. `str.replace`
matched nothing and returned the input unchanged, `ast.parse` still passed, and the failure
only surfaced as a `NameError` at run time. Every replacement in this file is now asserted
to have changed the source.

### Mann-Whitney U on the size_kde figures (2026-09-14)

The test was already being computed when `size_kde/` was added — `mannwhitneyu(x[pos],
x[~pos])` on `x = log10(lumen area)`, with `U / (n1 n2)` in `size_kde.csv`. What was missing
is that it was not on the figure. Each panel now carries `AUC …  P = …` right-aligned on the
title line, and `size_kde.csv` reports `mannwhitney_u`, `auc` and `p` explicitly. The column
`auroc_pos_larger` is renamed `auc`.

Two properties worth recording, both verified rather than assumed:

- **U and P are identical on raw area.** Checked across all 25 motifs: `U(log10 area)` equals
  `U(area)` exactly, because U is rank-based and log is monotone. The log sets the axis and
  the KDE bandwidth, nothing else.
- **AUC is the two-sample counterpart of the Spearman ρ on the scatter**, as expected. Signs
  and ordering agree across the six vasculature motifs:

| motif | U | AUC | 2·AUC − 1 | Spearman ρ(positive, log area) |
| --- | --- | --- | --- | --- |
| m1 SMC | 10,729 | **0.672** | +0.344 | **+0.253** |
| m0 T cell | 7,311 | 0.332 | −0.336 | −0.290 |
| m6 macrophage | 5,870 | 0.354 | −0.292 | −0.218 |
| m23 endothelial | 905 | 0.314 | −0.372 | −0.116 |
| m24 alveolar | 6,637 | 0.306 | −0.389 | −0.333 |
| m10 tumour-endo | 3,212 | 0.251 | −0.498 | −0.327 |

They are not numerically equal — ρ against a binary indicator is scaled by the group
proportions, so 2·AUC − 1 and ρ differ in magnitude while agreeing in sign and rank.

The caveats on the threshold stand: 0.5 selects the top ~17 % of vessels for m10 and 97 %
for m23, so AUC is comparable *within* a panel, not across panels; and vessels are spatially
autocorrelated, so the effect size and group sizes carry the argument, not P.

### Motif-specific vessels — the selection standard, visualised (2026-09-17)

`plot_motif_specificity_standard.py` → `figures_motif_specific_vessels/` (recomputes the picks
and asserts they equal `example_vessels.csv`; no change to `find_motif_specific_vessels.py`).
`m<k>_standard_rank` = the decision rule itself (x = within-section rank of k, y = max rank
of the competitors, shaded = margin ≥ 0.20); `m<k>_standard_raw` = the same vessels on raw
positive fractions with the absolute box (f_k ≥ 0.75, max other excl. m23 ≤ 0.25);
`standard_margin_counts` = n vessels with margin ≥ t per motif; `standard_picks_{raw,rank}` =
the 18 picks × six motifs, target cell outlined; `standard_counts.csv`.

Absolute-rule counts with m23 excluded from "others": m0 1, m1 45, m6 0, m10 1, m23 2, m24 0
(the earlier 0/2/0/0/2/0 included m23 as a competitor). **A rank margin does not imply raw
dominance.** Of the 18 picks, 4 are not the raw top motif (m0 #2 0.83 vs 0.96; m6 #3 0.75 vs
0.79; m10 #3 0.23 vs 0.95; m24 #2 0.83 tied with m1), and m24 #3 is top only at 0.14 — picked
because the target is high *for its section*. 5 picks pass the absolute rule: m0 #1, m1 #1,
m1 #3, m10 #1, m23 #1. m6 and m24 have none.

**Rule replaced (2026-09-17, same day).** The rank margin is out; selection is now on raw
fractions, per the user: target motif's ring positive fraction >= 0.5, then the other motifs
as low as possible. `find_motif_specific_vessels.select()`: eligible = f_k >= --min-frac
(0.5, the size_kde positive cut), ordered by max-other ascending (m23 excluded unless it is
the target), ties -> higher f_k, more ring cells, vessel id. `plot_motif_specificity_standard.py`
now draws only raw space (`m<k>_standard`, `standard_max_other_counts`, `standard_picks`);
the rank figures and the old picks moved to `superseded_rank_margin/`.

Picks (target / max other): m0 P17_LUAD__10 1.00/0.00; m1 P21_AIS__0, P17_LUAD__0,
P21_AIS__66 all 1.00/0.00; m6 P21_LUAD__61 0.70/0.38; m10 P21_LUAD__54 1.00/0.04;
m23 P17_LUAD__42 1.00/0.14; m24 P17_AIS__64 0.97/0.43. Eligible n = 138/230/76/52/290/129;
n with max-other <= 0.25 = 1/47/0/1/2/0. **m6 and m24 have no clean vessel at any threshold**
(lowest max-other 0.38 and 0.43) and all three m24 picks are P17_AIS, where m24 is high
section-wide — quote them as least-contaminated examples, not as specific vessels.

### size_kde: positive-only, four motifs, one axes (2026-09-17)

`plot_motif_size_distance.py --kde-combined` (`--kde-thresh` 0.9, `--kde-motifs` 1,10,23,24)
→ `figures_motif_profiles/size_kde/combined_thresh0.9.*` + `size_kde_combined.csv`. It reuses
the per-vessel ring fractions and returns before the distance field, so it touches none of
the other figures or CSVs. Curve colours are the VASC_MOTIFS colours of the pooled distance
figure; the 4-colour subset passes `validate_palette.py --pairs all` (worst normal pair 15.6).

At fraction >= 0.9: m1 175 vessels (median lumen 2,456 um2), m10 18 (610), m23 227 (1,788),
m24 58 (1,062) of 298. m1 is the only motif whose positive vessels reach the 10^4-10^5 um2
tail; m10's curve rests on 18 vessels and its shape should not be read closely.

`combined_thresh0.9_violin` is the same four groups as violins (log area on y, every vessel
drawn jittered, median + IQR in black). It is the honest version for m10: 18 points, so the
violin's shape is kernel smoothing over almost nothing, which the dots make visible.

### Cell-type composition of motif-positive cells (2026-09-17)

`plot_motif_celltype_composition.py` → `figures_motif_celltype/`: `m<k>_celltype_pie` (m1, m10,
m23, m24), `motif_celltype_stacked` (all four in the stacked format of
`vessel_niche_cluster.plot_cluster_cell_types`, tab20 by cell-type category order — the same
mapping as the vessel crops), `motif_celltype_composition.csv` (proportions + counts, with an
"all cells" background row). `--cells ring` restricts to the 30 um rings; default is all cells.

Positive cells are a large share of the tissue (m1 586 k, m10 704 k, m23 1.06 M, m24 620 k of
1.68 M), so every composition sits close to the background (Tumor_epi 23.5 %, Fibro 17.4 %,
T 10.9 %, Macro 8.9 %, Vas_Endo 8.2 %, Alveolar_epi 6.0 %). Read the ratio, not the share:
m1 SMC 14 % = 2.84x background; m24 Alveolar_epi 16 % = 2.69x; m23 Vas_Endo 12 % = 1.48x;
m10 has no type above 1.5x (Pericyte 1.44x, Tumor_epi 28 % = 1.19x). The composition figures
alone therefore under-state how motif-specific m1/m24 are and over-state m10.

`combined_thresh0.9_violin_sig` adds pairwise brackets to that violin (the plain version is
kept). Two-sided Mann-Whitney on log area for all 6 pairs, Holm-corrected with the
`multipletests` used in `core/glm.py`; only pairs with Holm p < 0.05 get a bracket
(* <0.05, ** <0.01, *** <0.001). 5 of 6 are significant: m1 vs m24 Holm p 4e-5, m1 vs m10
7e-4, m23 vs m24 5e-4, m1 vs m23 6.5e-3, m10 vs m23 6.5e-3; m10 vs m24 is not (0.44).
AUC (a > b): m1 vs m10 0.77, m1 vs m24 0.73, m23 vs m24 0.67, m1 vs m23 0.59, m10 vs m23 0.29.
Numbers in `size_kde_combined_pairwise.csv`. **Caveat to carry:** vessels are not independent
replicates -- 298 vessels from 4 sections / 2 patients -- so these P values are
pseudoreplicated; the AUCs and the medians carry the argument.

**2026-09-18.** `combined_thresh0.9_violin_sig` now sorts the violins by median lumen area,
largest left (m1 2,456 > m23 1,788 > m24 1,062 > m10 610 um2), and labels them with the
user's motif names (`NAME` in `plot_motif_size_distance.py`): m1 Arterial mural motif, m10
Tumor vasculature motif, m23 Vascular homeostasis motif, m24 Healthy alveolar motif. Same
tests, same Holm P; `size_kde_combined_pairwise.csv` now lists pairs in the sorted order with
name columns, so AUC is P(left group > right group) -- the two m10 pairs flip (m23 vs m10
0.71, m24 vs m10 0.56). The plain `_violin` keeps its m1/m10/m23/m24 order and short labels.

### Blood vessel or lymphatic? (2026-09-18)

`classify_vessel_lymphatic.py` → `figures_vessel_lymphatic/`. Per vessel (same 300 rings as
vessels.csv, asserted): share = Lym_Endo/(Lym+Vas_Endo) in the 30 um ring; M = PROX1|RELN vs
VWF|PLVAP detection over the ring's endothelial cells, cutoff 0.295 = midpoint of all-Vas_Endo
(0.005) and all-Lym_Endo (0.584). Lymphatic = share > 0.5 AND M >= cutoff.
**Result: 0 of 300 lymphatic** -- blood 292, discordant 1 (P21_LUAD__70: 2 Lym / 2 Vas, M 0.33),
indeterminate 7 (< 3 endothelial cells). No vessel has share > 0.5; max M 0.33. Holds for
every zone (30 / 15 / 10 um), marker set (PROX1, +RELN, +FLT4, +PDPN) and AND-rule tested by
the 2026-09-18 verification workflow (23 agents, all 300 calls reproduced independently).
Rule calibration on endothelial neighbourhoods (rule_calibration.csv): specificity 0.95,
**sensitivity only 0.58** -- PROX1|RELN is detected in 26 % of Lym_Endo, so the marker line
misses many true lymphatic neighbourhoods; the 0-lymphatic conclusion rests on share, which
never exceeds 0.5.
Vas_Endo PROX1+: 257 / 137,915 (0.19 %, mostly 1 count) -- "almost none", not none.
**This says the curated set is blood vessels, not that the tissue lacks lymphatics**: only
0.93 % of Lym_Endo lie within 30 um of a curated lumen (Vas_Endo 3.70 %). The Lym_Endo band at
the right of P21_LUAD is peribronchial (airway wall; the channel is an Airway_epi-lined airway,
FOXJ1/PIFO+), NOT a lymphatic lumen -- the maps now draw Airway_epi. Uncurated Lym_Endo-lined
lumen candidates from the verification (obsm/spatial um): P21_LUAD (10669, 6873) 15 Lym/1 Vas;
P21_LUAD (8369, 7599) 14/1; P21_AIS (6882, 998) 24/12; P21_LUAD (8676, 6820) 9/0.
P21_AIS__16 / __60 have no traced lumen (4-corner vessel-box fallback, src == 1); marked
`box_fallback` and "(vessel box)" in titles; calls blood under every geometry tested.
Old M over all ring cells kept as `M_allcells` (it capped a pure Lym_Endo ring near 0.55).
`plot_vessel_zoom.py --vessel P21_LUAD__70` → `figures_vessel_lymphatic/zoom_P21_LUAD__70`: 200 um
field, three panels (cell type / PROX1|RELN / VWF|PLVAP), 30 um ring outlined, the 4 ring
endothelial cells circled (ring asserted equal to the calls table). The two Vas_Endo sit on
the lower-right end of the lumen and are both VWF|PLVAP+; the two Lym_Endo sit together at
the upper-left end, one PROX1+, neither VWF|PLVAP+. Read: a blood vessel with a small
lymphatic profile abutting one end, not a lymphatic lumen.
`combined_thresh0.9_violin_sig` now brackets every tested pair; the one non-significant pair
(Healthy alveolar m24 vs Tumor vasculature m10, Holm p 0.44, AUC 0.56) reads "ns".
`distance_vs_positive/m<k>` (per motif, not pooled) now display 0-100 um (`XMAX_DIST`); the
splines are still fitted on 0-300 um (`DMAX`), so the curves and distance_vs_positive.csv are
unchanged (CSV byte-identical to the previous run). Pooled figure keeps 0-300 um.

### Tip-cell module in Vas_Endo (2026-09-18)

`plot_vas_endo_tip_module.py`: scanpy `score_genes` (COL4A1, KDR, ESM1; ctrl 50, 25 bins,
random_state 0) on the 137,915 Vas_Endo only, using X (log1p of counts per 1e4, verified).
Maps `figures_vas_endo_tip/tip_module_map_<section>` (one global colour scale, 1st-99th pct),
`tip_module_by_layer.csv`; per-cell scores in `vas_endo_tip_module/vas_endo_tip_scores.csv.gz`.
Detection in Vas_Endo: COL4A1 42 %, KDR 28 %, ESM1 9.7 %; score Spearman with COL4A1 0.76,
KDR 0.62, ESM1 0.46 (COL4A1 drives it most). Within every section the Tumor_region Vas_Endo
score higher than the Normal layers (mean: P17_AIS 0.39 vs -0.01..0.08; P17_LUAD 0.12 vs
-0.26..-0.14; P21_AIS 0.40 vs 0.05..0.21; P21_LUAD 1.42 vs 0.34..0.55). P21_LUAD is higher
in every layer, so the between-section level is confounded with section; compare within a
section. Large-vessel linings (the curated lumens) read low; high scores sit in the
capillary bed, in P21_LUAD along the tumour nests.

### Motif-specific example crops in morphology (2026-09-24)

`find_motif_specific_vessels.py --render morphology` (now the default; `--render points` keeps
the old dots) draws Xenium cell and nucleus polygons instead of centroids, following the
DSRCT `04_morphology_plots.py` pattern (cell body then nucleus, same colour per cell;
alpha 0.45 / 0.95, white 0.12 pt cell outline).

Geometry: `{XEN}/{section}_Xenium/{cell,nucleus}_boundaries.parquet` under
/data1/tanseyw/projects/spatial_data/linghua/j.ccell.2025.10.004. **Same micron frame as
obsm/spatial** -- verified by joining on `cell_id` per section (memory: never obs_names):
all 1,676,162 cells match with centroid difference exactly 0.0, so no transform is applied.
Read per crop with a pyarrow bbox filter (+30 um pad, 0.1-0.2 s per file); nuclei are grouped
by (cell_id, label_id) because 16/469 cells in a test crop are multinucleate.
Coverage over the 18 drawn crops: 11,970 cells, 0 without a cell boundary; a few have no
nucleus polygon (0-91 per crop, most in P21_LUAD__54), which are drawn as cell body only.
Runtime 57 s for the six figures. Selection rule and every number are unchanged.

### Example vessels on H&E (2026-09-24)

`plot_vessel_he_morphology.py` → `figures_motif_specific_vessels/m<k>_specific_examples_he`
(separate figures; the `m<k>_specific_examples` grids are untouched). One row per motif,
same picks and same fields as those grids: aligned H&E underneath, the Xenium cell polygon of
every cell filled translucent (alpha 0.68, no outline) in its tab20 cell-type colour; the
curated lumen is cyan at alpha 0.6 (neighbouring lumens grey, same alpha). **Nuclei are not drawn** and **Tumor_epi is recoloured** to #54278f (CT_COL) --
tab20's grey vanishes on H&E; the deep purple stays clear of tab20's #9467bd (Myeloid_other)
and #c5b0d5 (NK). The other figures keep grey for Tumor_epi. (Earlier version: cell boundary
solid + nucleus boundary dashed as lines.) `load_shapes` now takes `kinds`, so this figure
reads only the cell parquet.

H&E: `{XEN}/{section}_Xenium/HE/aligned_HE_rgb.ome.tif` (present for all four sections,
1.4-3.0 GB, pyramidal 9 levels, OME PhysicalSize 0.2125 um/px). It is already registered to
the Xenium frame, so micron = pixel x 0.2125 with origin (0,0) and no further transform;
the overlay lands exactly on the H&E nuclei. `he_crop()` picks the pyramid level giving
~2000 px per crop and reads only the crop's window (~1 s each; 60 s for all six figures).
Axes use `set_ylim(hi, lo)` because image rows run downward in y.

### Per-motif columns (2026-09-24)

*(Superseded within the same day: these two figures were first built per tissue section; the
user retracted that and asked for one column per motif. The per-section PNG/PDF/SVGs and
`size_kde_density_by_section.csv` were moved to `superseded_per_section/` in each folder;
`celltype_positive_fraction.csv` still carries the per-section numbers.)*

**Density instead of violins.** `plot_motif_size_distance.py --kde-combined` now also writes
`size_kde/density_thresh0.9` + `size_kde_density.csv`: one column per motif (m1 / m23 / m24 /
m10 by median), the violin cut in half and laid on its side, vessels as ticks and the median
as a black bar, all sections pooled. `sharex`/`sharey`, so every column has the same scale.
The pooled `_violin` and `_violin_sig` figures are unchanged. (Per section the groups would be
too thin to draw: m10 n = 13 / 1 / 1 / 3 and m24 n = 27 / 0 / 28 / 3 -- see
`superseded_per_section/`.)

**Cell-type bars, one column per motif.** `plot_motif_celltype_composition.py` now also writes
`celltype_positive_fraction` + `celltype_positive_fraction.csv`: one column per motif, four
horizontal rows each (Tumor_epi, Vas_Endo, SMC, Alveolar_epi), all sections pooled. No axes:
each bar is a thin black 0-100 % frame filled to **% of that cell type's cells that are motif
positive**, and the fill uses that cell type's tab20 colour (the mapping shared with the pies,
the stacked bar and the vessel crops). Only the first column carries the row labels --
with `sharey` the ticks are shared, so setting blank labels on a later column wipes them all. The other reading (share
of positive cells that are that type) is the `share_of_positive_cells` column of the CSV,
which also keeps the per-section breakdown. Pooled: SMC 99.5 % for m1, Vas_Endo 93.7 % for
m23, Alveolar_epi 99.4 % for m24; m10 is the flat one (Tumor_epi 50 %, Vas_Endo 43 %,
Alveolar_epi 43 %, SMC 14 %). Per section, SMC-m1 is 98.8-99.7 % everywhere while
Alveolar_epi-m24 is 99 % in the AIS sections but 82 % in P17_LUAD.

### Whole-section morphology maps, cells coloured by motif state (2026-09-24)

`plot_section_motif_morphology.py` → `figures_section_motif_morphology/<section>_m<k>`
(4 sections x m1/m10/m23/m24 = 16 figures; `--motifs`, `--sections`, `--dpi`) +
`section_motif_cell_counts.csv`. Every cell drawn as its Xenium cell polygon, filled dark
(#3b0f70) when `motif_<k>_state` is positive and pale grey otherwise; curated lumens outlined
cyan on top. No H&E.
Geometry: whole `cell_boundaries.parquet` per section, split into per-cell arrays via
factorize + `np.split` on the contiguous cell_id runs (~0.5 s for 230k cells); the cell layer
is `rasterized=True` so the PDF/SVG stay ~0.5 MB instead of carrying 600k vector polygons.
3.5 min for all 16 figures; total 24 MB.
Cells segmented but absent from the adata are not drawn (P17_AIS 47,652 of 230,030; the adata
is the high-quality subset) -- the counts in the CSV are over drawn cells only.
Positive fractions are whole-section, so they are much higher than the ring fractions:
m23 0.53-0.78, m24 0.12 (P17_LUAD) to 0.88 (P21_AIS), m1 0.27-0.42, m10 0.32-0.61.
m1 visibly rings the cyan lumens; in P21_LUAD m10 fills the tumour front and avoids the
vessel-rich centre.

### UpSet of the four vessel motifs (2026-09-25)

`plot_motif_upset.py` → `figures_motif_upset/`: `motif_upset_cells`, `motif_upset_vessels`,
`motif_upset_counts.csv` (all 16 combinations, both units, plus a 0.9 vessel column).
Hand-rolled (no upsetplot/pyupset in the env, no prior code in the repo): intersection bars
on top, set sizes left, membership matrix below, combinations **exclusive** and sorted by
size. Labels are the user's motif names only, no m-numbers on the figure.
**m1's display name changed** to "SMC vascular stabilization motif" (was "Arterial mural
motif"); the other three are unchanged. Only this script uses the new name so far.

Cells (1,676,162): sets m23 1.06 M > m10 704 k > m24 620 k > m1 586 k; 299,524 cells (17.9 %)
carry none of the four; the largest combination is m23 + m10 (200,099), then
m23 + m24 + m10 (167,884), then m23 + m1 (146,998); all four together 91,924 (5.5 %).
Vessels (300, ring fraction >= 0.5): **every vessel carries at least one** (none = 0);
m23 + m1 122, m23 + m1 + m24 77, m23 alone 32, m23 + m24 + m10 25, all four 18.
**m10 never occurs alone in a vessel** and only ever with m23; m23 is in 290 of 300.

### 16-state OKLab palette, shared across sections and figures (2026-09-25)

`motif_state_palette.py` is the single source: MOTIFS/NAME/ANCHOR, OKLab conversion
(Ottosson matrices; round trip 4e-15), `state_code` (bit-packed, bit j = MOTIFS[j]),
`state_palette()` → {0..15: hex}, `state_label`, `draw_state_legend`, `palette_separation`.
Consumers: `plot_section_motif_states.py` (the maps) and `plot_motif_upset.py` (intersection
bars now carry the state colour).

Palette (anchors re-chosen 2026-09-25): m1 #d94527 red-orange, m23 #2b86d8 blue,
m24 #74a54b green, m10 #e98aea orchid; a multi-positive state is the OKLab mean of its
anchors, minus LSTEP = 0.13 lightness per extra motif; all-negative #ededed.
**Why they changed:** the first version reused the established motif colours, but m1 #e69f00
and m24 #d55e00 are both orange, so every mix involving them landed in one brown band. The
new four were found by random + local search maximising the closest of the 16 state colours,
subject to: all 16 in gamut, anchor lightness inside validate_palette.py's light band,
chroma >= 0.10, and every anchor pair clearing the normal floor (>= 15) and the colour-vision
floor (min(protan, deutan) >= 6). Closest pair of the 16 went 0.054 -> 0.083 (0.080 as drawn
after the hex rounding), and the four anchors now come back RESULT: PASS from the validator
(worst normal pair 24.8; one WARN, red-orange vs green under deuteranopia 6.2).
**Cost:** these two figure families no longer share single-motif colours with the density,
cell-type-bar, violin and per-motif section figures, which still use the old four.
**Why the lightness step:** the plain mean cannot separate 16 states -- for symmetric anchors
the mean of two opposite motifs equals the mean of all four (measured 0.000). LSTEP 0.13
breaks that tie (`palette_separation`) and makes "more motifs" read as "darker".
Chroma-restore (rescaling mixed a/b to the mean anchor chroma) measured worse, 0.031; kept as
`chroma="restore"`, and `lstep=0` gives the plain recipe. Even at 0.054, triples and the
quadruple are brown and hard to tell apart in the map -- the legend's dot matrix carries
identity, not the colour alone.

`plot_section_motif_states.py` → `figures_section_motif_morphology/<section>_motif_states`
(+ `motif_state_palette.csv`, `motif_state_counts_by_section.csv`). Cells drawn in order of
how many motifs they carry so multi-positive cells sit on top; lumens outlined black (cyan
would clash with the m23 anchor). 70 s for the four sections.
Regional structure is obvious in P21_LUAD: the tumour centre is sky blue (Vascular
homeostasis alone) while the periphery is dark (three or four motifs); 177,106 of 560,183
cells there carry none of the four, against 10,212 of 182,378 in P17_AIS.

**Venn version (2026-09-25).** `plot_motif_venn.py` → `figures_motif_upset/motif_venn_cells`,
`motif_venn_vessels`. Four ellipses (four circles cannot make all 15 regions); each region is
filled with its state colour from `motif_state_palette`, so a region matches the same
combination in the UpSet and the section maps, and the background is the all-negative colour.
Regions are rasterised on a 1400x1400 grid and each count sits at that region's pole of
inaccessibility (distance transform), so no label escapes its region. The ellipse geometry was
searched for a layout where all 15 regions exist and the smallest is as large as possible
(6,137 px on a 700-grid; the usual pyvenn layout gives 3,399 and my first hand-set one missed
two regions entirely). `plot_motif_upset.py` gained `load_cells()` / `load_vessels()` so both
scripts share one definition of the sets.
**Areas are not proportional** -- impossible for four sets; the Venn shows which combinations
exist, the UpSet is what sizes should be read from.
Vessels: the zeros are informative -- Tumor vasculature alone 0, with Healthy alveolar alone 0,
SMC + Healthy alveolar 0; every vessel carrying Tumor vasculature also carries Vascular
homeostasis.

**Per-section UpSet (2026-09-25).** `plot_motif_upset.py` also writes
`motif_upset_cells_by_section` + `motif_upset_by_section.csv`: same 15 columns and the same
membership matrix, but one bar row per tissue section and the bar is that combination's share
of **that section's** cells (%), on a shared y scale; each row carries its own "none of the
four" share and n. `load_cells(with_sections=True)` returns the section label per cell.
The sections disagree strongly: Vascular homeostasis + Healthy alveolar + Tumor vasculature is
the top state in P21_AIS (24 %) and P17_AIS (30 %) but near zero in P17_LUAD, whose top state
is Vascular homeostasis + Tumor vasculature (18 %, and 15 % in P21_LUAD). "None of the four"
runs 2 % (P21_AIS), 6 % (P17_AIS), 16 % (P17_LUAD), 32 % (P21_LUAD).

**Y axis flipped (2026-09-28).** Both whole-section families (`<section>_m<k>` and
`<section>_motif_states`) now plot y increasing upward (`set_ylim(lo[1], hi[1])`) instead of
image order, matching the lymphatic maps and the vessel crops; the 1 mm bar moved to the
bottom-left inside the axes in both. Note the sections are therefore mirrored top-to-bottom
against the H&E crops (`m<k>_specific_examples_he`), which stay in image order because the
raster is tied to it. Nothing but orientation changed -- all counts and colours are identical.

### Cell-type colours corrected to the run's own palette (2026-09-28)

The run's figures colour cell types with `alarmist.plotting.utils.get_cell_type_colors`, which
is `plt.get_cmap("tab20", len(cell_types))` -- tab20 **resampled** to 19 -- indexed over the
h5ad category order. Every luad_revision figure until now used plain `tab20[i % 20]`, which
shifts the whole assignment: SMC came out #f7b6d2 pink instead of #7f7f7f grey, Tumor_epi
#c7c7c7 grey instead of #bcbd22 olive, Vas_Endo #bcbd22 instead of #dbdb8d.
New `celltype_palette.py` calls the package function (adds <repo>/src to sys.path; the package
is not pip-installed in comp-liana) and is now used by `find_motif_specific_vessels.py`,
`plot_vessel_he_morphology.py`, `plot_motif_celltype_composition.py` and `plot_vessel_zoom.py`;
those figures are regenerated. The H&E overlays keep the requested Tumor_epi -> #54278f purple
override on top of the run palette -- which now differs from the olive Tumor_epi everywhere
else, so drop `CT_COL` if that inconsistency is not wanted.
**Still on plain tab20** (older folders, not regenerated): `cluster_celltype_composition` in
every `vessel_kmeans_*_composition/`, `neighborhood_cluster_cell_types` in the
`vessel_niche_*` folders, and the cluster example crops.

**2026-09-28 (colour alignment).** `celltype_positive_fraction` (and the pies / stacked bar
from the same script, which share its palette) now draw Tumor_epi in #54278f, the purple the
H&E overlays use; every other cell type keeps the run palette. The script also imports `NAME`
from `motif_state_palette`, so m1 reads "SMC vascular stabilization motif" there too.
`size_kde/density_thresh0.9` now takes its four colours from `motif_state_palette.ANCHOR`
(the UpSet / Venn / state-map anchors: m1 #d94527, m23 #2b86d8, m24 #74a54b, m10 #e98aea),
is long and flat (62 x 30 mm per panel), titles carry only the m-number and n (no motif
name), every panel repeats the x label, and Q1 / median / Q3 are vertical ticks (median taller);
`size_kde_density.csv` gained `q1_area` and `q3_area`. The violin figures keep MCOL.
**Still on the old m1 name**: `plot_section_motif_morphology.py` (the 16 per-motif maps).

**Presentation overrides (2026-09-28).** `celltype_palette.PRESENTATION` now holds the two
deliberate departures from the run palette, shared by `plot_vessel_he_morphology.py` and
`plot_motif_celltype_composition.py`: Tumor_epi #54278f (ColorBrewer Purples-8; the run's olive
#bcbd22 is hard to read on H&E) and **SMC <-> Plasma swapped** (SMC #f7b6d2 pink, Plasma
#7f7f7f grey) on request. Every other figure keeps `get_cell_type_colors` untouched, so the
point/morphology crops (`m<k>_specific_examples`) and `zoom_P21_LUAD__70` still show grey SMC.

### H&E overlay palette: distance from the measured stain (2026-09-29)

`make_he_palette.py` (rewritten), dict in `he_overlay_palette.HE_CELLTYPE_COLORS`, consumed
by `plot_vessel_he_morphology.py` and, from the same day, by
`plot_motif_celltype_composition.py` (pies, stacked bar, celltype_positive_fraction), so
the two families share one palette; `celltype_palette.PRESENTATION` is no longer used. `coloraide` 8.13 pip-installed into comp-liana.
Supersedes the fixed-hue-band version from earlier the same day.

Method: 50,000 tissue pixels sampled from the 18 crops the figure draws, balanced across the
four section images, OKLab L kept in 0.25-0.92 (drops glass and near-black nuclei). k-means
k=8 on those pixels: **all eight centres are pink/purple/magenta, hue 298-341** (largest
#dddae7 pale 29 %, then #a33e84, #882d74, #b6569b, #691e62, #c676b3, #47114c, #d4a2cc) --
eosin here reads magenta-pink, there is no separate red cluster. Candidate grid L 0.35-0.85 /
C 0.10-0.22 / H 0-358 step 2 = 6,174 in gamut; each scored by OKLab distance to its 50th
nearest sampled pixel (KD-tree). Survivors: 2,604 at 0.10, **2,013 at 0.12 covering 124 of 180
hue steps**, 1,478 at 0.14, 1,064 at 0.16. **Layout rule (revised 2026-09-29, second pass):** steps 1-4 now only decide which hues are
allowed, and the colours follow a fixed rule -- chroma 0.14 for all of them, the 19 types laid
evenly along the allowed hues in group order (within-group 9.5 deg, between-group 19.0 deg),
lightness cycling 0.55 / 0.75 / 0.65 continuously (restarting per group made the boundary
pairs the closest in the palette, 0.034). A colour failing the stain check is moved in
LIGHTNESS ONLY, to the passing lightness that stays furthest from the rest of the palette:
only Tumor_epi needed it (0.109 -> 0.121 at L 0.69).
**The allowed-hue mask must be taken at the fixed chroma**: over all chromas it admits hues
that only clear the threshold at C 0.20+, which C 0.14 cannot reach -- that bug left five
colours below the stain threshold. At C 0.14 the usable range is 46-250 deg (101 steps).
Chroma was checked: 0.16 and 0.18 give FEWER usable hues (75 and 46 steps) because high
chroma falls out of sRGB gamut, so 0.14 is the better choice.
Result: every colour is at least 0.120 from the stain; the closest pair in the palette is 0.057
(Vas_Endo/SMC) and nine pairs sit below 0.08. Those are same-lightness pairs three positions
apart: at chroma 0.14 a 28 deg hue gap is only ~0.07 in OKLab, so with three lightness levels
they cannot clear 0.08. Measured alternative, same rule otherwise: a four-level cycle
(0.50 / 0.75 / 0.62 / 0.85) gives closest pair 0.069 and only three pairs below 0.08.
Diagnostics: `he_stain_clusters`, `he_palette_candidates` (hue vs lightness of survivors,
chosen circled), `he_palette_swatches`, `he_palette.csv`.
**Known compromise:** the greedy piles picks at the far end of each segment (vascular/stromal
lands at 124-146 of its 84-148 band), so five colours sit in the green band. Spreading picks
evenly in hue, or maximising mutual distance subject to the stain threshold, would fix it.

### H&E overlay palette rebuilt in OKLCH (2026-09-29, superseded the same day)

`make_he_palette.py` generates it, `he_overlay_palette.HE_CELLTYPE_COLORS` holds it (one dict,
hand-editable), `plot_vessel_he_morphology.py` consumes it. **Only that figure changed** --
the composition figures and the crops keep the run palette (+ `celltype_palette.PRESENTATION`).
`coloraide` 8.13 was pip-installed into the comp-liana env for this.

Stain band measured, not guessed: 3.9 M pixels from the 18 crops the figure actually draws,
converted to OKLCH, chroma-weighted hue mass (chroma <= 0.02 dropped as it is unstained glass).
90 % of the mass in 316-353 deg, peak 339, 99 % arc 288-5 deg; padded by 20 deg the excluded
band is 268-25 deg, so hue runs 25-268 deg -- close to the 280-40 guess, a little warmer.
Mean stained pixel #9b508b, used as the swatch background.

Layout: hue = group (epithelial / vascular-stromal / lymphoid / myeloid, order as given),
within-group step 11.6 deg, between-group 23.1 deg (GAP_FACTOR 2), chroma fixed 0.15,
lightness cycling 0.55 / 0.78 / 0.66 across the whole ramp so neighbours differ in both;
sRGB via coloraide `.fit()`.
Checks: every pair >= 0.08 by coloraide delta_e "ok". Three pairs needed a +-0.06 lightness
nudge (pDC vs Macro 0.066, Vas_Endo vs SMC 0.070, cDC vs Myeloid_other 0.078); after repair the
closest pair is 0.085 (Langhans_cell/Neutro, T/Plasma, Macro/Mast), median pair 0.24, none
below the floor. Swatches on the sampled background and on white:
`figures_motif_specific_vessels/he_palette_swatches`, values in `he_palette.csv`.
**Weakest against the background: Tumor_epi #b94642 at delta_e 0.12** -- it sits at the warm
edge of the allowed band. Raising PAD in make_he_palette.py pushes the whole ramp away from
the stain at the cost of tighter spacing.

**Third pass (2026-09-29): vivid instead of flat chroma.** Chroma is no longer fixed -- each
colour takes 90 % of the max in-gamut chroma at its own L and hue (`CHROMA_FRAC`, binary search
in `max_chroma`), so chroma now runs 0.085 (teals, where sRGB simply has little) to 0.283
(magenta). Over the yellow/olive band 70-125 deg the lightness cycle is shifted up to
0.72 / 0.90 / 0.81 (`YELLOW_LIGHTNESS`), because dark yellow reads brown; the repair keeps
L >= 0.70 there. The hue scan is redone under the chroma rule: **144 of 180 hues now clear 0.12**
(against 101 at fixed chroma 0.14), so the ramp spans 44-330 deg with 13.6 deg within a group
and 27.2 deg between -- wider spacing, as asked.
Result: closest to the stain 0.120, closest pair 0.065 (Macro/Myeloid_other), four pairs below
0.08 (was nine). Seven colours needed the lightness-only repair, two of them large (cDC
0.75 -> 0.42, Langhans_cell 0.65 -> 0.46), which breaks the neat alternation in the myeloid run.
**Watch this:** the widened range puts Myeloid_other / Neutro / Mast at 302-330 deg, i.e. inside
the stain's own hue band; they clear the check only because they are far more saturated than the
tissue. If vivid magenta on pink reads badly, cap the ring (the allowed mask) below ~300 deg.
`he_palette_swatches` now shows the new palette above the previous one on both backgrounds.

### Niche composition in absolute cells (2026-10-01)

`plot_niche_cell_counts.py` → in every `vessel_niche_ring30_pervessel*` folder:
`neighborhood_cluster_cell_counts` (same stacked bars as `neighborhood_cluster_cell_types`,
but bar height = number of ring cells, so niche size is visible), `niche_cell_type_counts.csv`
(counts per niche x cell type, plus n_vessels and n_cells), and `neighborhood_cluster_cell_types`
re-rendered so both use the run's own cell-type palette instead of plain tab20.
`neighborhood_counts.npy` holds proportions in the per-vessel runs, so counts are recovered as
proportion x n_ring from niche_labels.csv -- exact, max deviation from integer 1.4e-14 and the
total matches sum(n_ring) = 20,279 ring cells over the 300 vessels.
The size skew the percentage figure hides: in the auto-K run niche 3 holds 9,280 of the 20,279
ring cells while niches 4 and 7 hold 97 and 121; at K=3 the split is 9,511 / 7,562 / 3,206.
The per-cell run (`_nb30_percell`) was left alone -- there a "count" would sum overlapping
neighbourhoods, which is not the same quantity.

**Extended to the motif k-means runs (2026-10-01).** `plot_niche_cell_counts.py` also writes
`cluster_celltype_counts.{png,pdf,svg}` + `cluster_celltype_counts.csv` into every
`vessel_kmeans_ring30_*_composition` folder (18 runs: the two defaults, frac k2-k10/k15/k20/
k25/k30, loading k5-k8), and re-renders `cluster_celltype_composition` on the run palette.
Those runs cluster the same 300 vessels on the same rings, so the per-vessel counts come from
the niche table joined on vessel_id; the vessel sets are asserted identical in every run.
Size skew at k10 (frac): cluster 1 holds 5,939 of the 20,279 ring cells, cluster 5 only 407;
at k25 the smallest cluster is 32 cells, which the percentage figure drew full height.

**Colour swap + vessel counts (2026-10-01).** The count figures label each bar with the number
of vessels in that cluster, and `plot_niche_cell_counts.SWAP` is now
`celltype_palette.PRESENTATION` itself (reused, not copied): Tumor_epi #54278f deep purple plus
the 3-cycle SMC #f7b6d2 pink (was grey), Plasma #2ca02c green (was pink), Langhans_cell #7f7f7f
grey (was green). No other cell type changes, and the three swapped colours are permuted among
themselves so none is duplicated. These figures now match the H&E overlays' overrides. Applied to every figure this script writes --
the 7 niche runs and the 18 motif k-means composition folders, counts and percentage versions
alike, so the two families stay consistent.

**The singleton niche, looked at (2026-10-01).** `plot_vessel_he_morphology.py --vessels
P17_AIS__6 --fov 700` (new `--vessels` / `--fov` flags draw any vessel, not just the motif
picks) → `he_vessels_P17_AIS__6`. It settles the question: the curated lumen is a small vessel
embedded in a bronchiole wall, with Airway_epi epithelium wrapping it on every side --
52 % of the cells in the 700 um field are Airway_epi, against 0.36 % of all ring cells. So the
singleton cluster at every K >= 4 is one genuinely peribronchial vessel, not a mis-segmentation,
but with n = 1 it is an outlier rather than a niche.

### Defending the choice of K for the vessel niches (2026-10-01)

`niche_k_robustness.py` → `figures_niche_robustness/`: `niche_k_stability`,
`niche_k_silhouette`, `niche_consensus_k3`, `niche_robustness.csv`,
`niche_perturbations_k3.csv`. Features are the per-vessel ring proportions z-scored and
**clipped at +-5**, which is what stops the 16-sigma peribronchial vessel owning a cluster.

| K | silhouette | stability (ARI, 200 resamples) | smallest cluster |
|---|---|---|---|
| 2 | 0.152 | 0.86 +- 0.21 | 128 |
| **3** | 0.156 | **0.91 +- 0.04** | 49 |
| 4 | 0.170 | 0.74 | 43 |
| 9 | 0.164 | **0.57** | 9 |
| 10-12 | 0.18-0.20 | 0.61-0.65 | 7-10 |

So K=3 is the defensible choice: highest stability by a wide margin and the tightest spread,
while the KneeLocator's K=9 is the *least* stable setting tested. Consensus at K=3:
within-cluster co-clustering 0.93, between-cluster 0.03.
Perturbations at K=3 (ARI vs the reference labels): seeds 0.95-1.00, raw z features 0.96,
rare types pooled 0.90, 80 % subsample 0.94 -- but **leave-one-section-out 0.56-0.90**, worst
when P17_AIS is dropped. That is the honest limit: the partition is reproducible against noise
and feature choices, not against removing a section, with 4 sections from 2 patients.
**What can be claimed:** silhouette is 0.15-0.20 at every K, so the vessels are a continuum;
K=3 is a reproducible coarse partition of that gradient, not evidence of three discrete niches.

## 2026-10-02 — K selection replaced by per-cluster bootstrap Jaccard (`niche_k_jaccard.py`)

The global-index battery (bootstrap ARI, silhouette, consensus matrix, perturbation table) was
replaced at the user's request with the standard per-cluster protocol. Old outputs moved, not
deleted: `figures_niche_robustness/superseded_ari_consensus/`.

Protocol: scan K = 2..9 (9 = the elbow); at each K, 100 subsamples of 80 % of the 300 vessels
*without replacement*, recluster, score each reference cluster by its best Jaccard against the
subsample's clusters; choose the **largest** K at which **every** cluster's mean Jaccard ≥ 0.75
(not the most stable K — that is always K = 2). Features: ring proportions, z-scored, clipped
at ±5.

Result (`niche_k_choice.csv`): min per-cluster Jaccard K2 0.951, **K3 0.903**, K4 0.591,
K5 0.546, K6 0.465, K7 0.495, K8 0.382, K9 0.372. Only K = 2 and K = 3 clear 0.75, so **K = 3**.
The drop is a cliff between 3 and 4, not a slope.

Per-cluster at K = 3 (`niche_jaccard_per_cluster.csv`): c0 n=109 J=0.99 (SMC 0.56 / Vas_Endo
0.28), c1 n=142 J=0.96 (Vas_Endo 0.24 / Fibro 0.22), c2 n=49 J=0.90 (T 0.27 / Vas_Endo 0.13).
All three span 4 sections, but **c2 is 79.6 % P17_AIS** — the same section whose removal gave the
worst leave-one-section-out ARI (0.56) in the superseded analysis. Report it as section-skewed.

Small clusters, per the same rule: the only cluster ≤ 10 vessels in the scan is K = 9 cluster 5
(n = 9, J = 0.536 ± 0.432) — low Jaccard, i.e. noise, which is evidence for a smaller K, not a
rare population. Without the z-clip the earlier singleton returns: raw-z K = 4 cluster 3 is
`P17_AIS__6` alone, J = 0.637 ± 0.473, **1 section** — high-variance, single-section, n = 1, the
textbook artefact case (it is the vessel embedded in a bronchiole wall, 45 % Airway_epi in its
ring). Raw-z min Jaccard per K: K3 0.901, K4 0.558, K5 0.437, K6 0.410, K7 0.296, K8 0.334,
K9 0.501 — the choice of K = 3 does not depend on the clip.

Clustree (`niche_clustree`, nodes = clusters sized by vessel count, coloured by Jaccard, edges =
share of a K cluster going to K+1, columns ordered by parent barycentre): from K = 3 on, the
large SMC-high cluster keeps ~105 vessels and stays yellow at every K while new clusters are
shaved off the other two — new nodes are outlier-shaving, not new structure.

Figures/tables: `figures_niche_robustness/{niche_jaccard_by_k,niche_clustree}.{png,pdf,svg}`,
`niche_jaccard_per_cluster.csv`, `niche_k_choice.csv`; synced to
`download_20260902/figures_niche_robustness/`.

**Consensus added for every K (same run, same subsamples).** `niche_consensus_all_k` (2 × 4
panel grid) plus `niche_consensus_k2 … k9` individually; within/between columns added to
`niche_k_choice.csv`: within-cluster co-clustering K2 0.959, K3 0.953, K4 0.837, K5 0.758,
K6 0.729, K7 0.681, K8 0.606, K9 0.623; between-cluster 0.02–0.08 throughout. The within column
falls with the Jaccard cliff (0.95 → 0.84 between K=3 and K=4), so consensus and per-cluster
Jaccard agree on K = 3; the K ≥ 4 panels show pale, speckled diagonal blocks rather than the
solid blocks at K = 2–3.

**Niche × motif correlation, every K (`plot_niche_motif_correlation.py`).** Per-vessel Pearson r
between each niche's 0/1 membership and each of the four picked motifs' 30 µm ring positive
fraction (features.npy of `vessel_kmeans_ring30_m0-1-6-10-23-24_frac`, joined to the niches by
vessel_id — the two tables are NOT in the same row order). Reds, vmin 0, shared vmax 0.60;
negative r left white with the number printed; each niche's best motif boxed. Motif axes are
labelled by NAME, not index (SMC vascular stabilization / Vascular homeostasis / Healthy
alveolar / Tumor vasculature), rotated 38 deg, coloured with `motif_state_palette.ANCHOR`, in
the usual m1/m23/m24/m10 order; the CSV also carries the names, not m-numbers.

Mean best-match r by K: K2 0.42, K3 0.38, K4 0.44, K5 0.35, K6 0.28, K7 0.30, K8 0.26, K9 0.25.
Splitting past K=4 does not buy motif correspondence — it falls. At K=3: niche 0 (n=109,
SMC-high) ↔ m1 r=0.52; niche 1 (n=142) ↔ m23 r=0.17 only; niche 2 (n=49, T-high) ↔ m10 r=0.46
(m24 0.33). K=4 is marginally higher than K=3 because its extra cluster separates an
m24-leaning group (r=0.54), worth one sentence in the paper rather than a change of K: K=4 fails
the Jaccard rule (0.59).

Honest reading: the niches are correlated with, not equivalent to, the motifs — the strongest
cell-type niche ↔ motif correspondence in the whole scan is r ≈ 0.5, so cell-type composition
and LRI motif state are related but not interchangeable descriptions of a vessel.

Outputs: `figures_niche_robustness/niche_motif_corr_k2 … k9`, `niche_motif_corr_all_k`,
`niche_motif_corr.csv`; synced to the download bundle.

## 2026-10-06 — vessel size: perimeter diameter instead of lumen area (`lumen_geometry.py`)

Raised by the user: a lumen is flattened when the block is cut and mounted, and veins collapse
more than arteries, so **area** is a biased size metric. The invariant under flattening is the
traced **perimeter**, so the metric to prefer is `d_perim = P / π` (the diameter of the circle
with the same circumference). Note `d_area = 2√(A/π)` is a monotone function of A, so on a log
axis it is the area plot with a relabelled axis and *identical* Spearman/Mann-Whitney values —
only `d_perim` can change a ranking.

`lumen_geometry.py` (new) computes per curated lumen: shoelace area (max 5.9e-10 relative
difference from the curated `lumen_area_um2`, so the stored areas are reproduced exactly),
perimeter, d_area, d_perim, Feret max/min, circularity 4πA/P², aspect ratio. Written to
`lumen_geometry.csv` (309 lumens; 2 fall back to the vessel bounding box, `src = 1`, excluded
from the summaries because a rectangle's perimeter is not the vessel's).

How squashed are they: **median circularity 0.70, 10th percentile 0.38, 20 % below 0.5**; median
Feret aspect 1.83, 42 % above 2. So flattening is real and worth taking seriously.

Does the metric change anything: Spearman(area, d_perim) = 0.976, median rank shift 10 places
(max 110 of 307). Re-running the figures with `--size diameter`: per-motif size correlations move
by ≤ 0.04 (m10 −0.303 → −0.263 is the largest), every sign and every call is unchanged; the
positive-vs-negative AUCs move by ≤ 0.03; all six pairwise violin comparisons keep the same
significance call after Holm. **One borderline flip: m23's size_kde p goes 0.046 → 0.056**, i.e.
it was never a real result. Conclusion: the diameter version is the defensible one to publish,
and it does not change any claim.

Does squishing track arteriality? Circularity vs ring SMC fraction is ρ = −0.20 (p = 6e-4) —
the *opposite* of the expectation — but inside area quartiles the correlation vanishes
(−0.17, −0.14, +0.03, −0.02, all n.s.): it is a size confound, big lumens are less circular
(ρ(d_perim, circularity) = −0.57). Note also that circularity conflates flattening with a wavy
traced outline, which large arteries have anyway; the Feret aspect ratio separates the two and is
in the CSV. So: no evidence here that arteries resist squashing, and no need to claim it.

Scripts changed: `plot_motif_size_distance.py` and `plot_vessel_size_by_cluster.py` gained
`--size {area,diameter}` (axis labels, stat column names and file stems follow the metric; the
diameter run writes to `figures_motif_profiles_diam/` and `vessel_size_by_cluster_k10_diam/`, so
the area versions are untouched). Also fixed a stale display name in both
`plot_motif_size_distance.py` and `plot_section_motif_morphology.py`: m1 was still
"Arterial mural motif", now "SMC vascular stabilization motif"; the area-metric KDE figures were
re-rendered to pick it up.

Outputs: `figures_vessel_size/lumen_size_metrics`, `lumen_geometry.csv`,
`figures_motif_profiles_diam/{size_vs_positive,size_kde}/`, `vessel_size_by_cluster_k10_diam/`;
all synced to `download_20260902/`.

## 2026-10-06 (later) — the z-clip dropped, and the labels now match the stored niche runs

Two defects the user caught: the Jaccard / consensus / clustree / motif-correlation figures
were computed on z CLIPPED at +/-5, while everything on disk
(`vessel_niche_ring30_pervessel*`) was cut on unclipped z, so the singleton cluster that
appears from K = 4 in the stored folders was invisible in those figures. Fixed by dropping the
clip, after checking it changes nothing that matters.

Clip sweep, smallest cluster per K (only 2 of the 5,700 z values exceed |z| = 10, both in
P17_AIS__6 and P17_AIS__62):

    clip    K3   K4   K5   K6   K7   K8   K9
    none    53    1    1    1    1    1    1
    12      51   13    3    9    8    2    2
    10      47   13   12   13    3    8    8
    8       51   12   12    9    9    9    8
    5       49   43   37   12   11   13    9

Clipping at 10 removes the 1-vessel cluster but still leaves a 3-vessel one at K = 7; only the
heavy clip at 5 gives an even spread. The K choice, however, is the same everywhere -- min
per-cluster Jaccard at K=3 is 0.90-0.91 and <= 0.6 from K=4 under no clip, clip 10 and clip 5
alike, so all three choose K = 3. Unclipped is now the default in `niche_k_jaccard.py` and
`plot_niche_motif_correlation.py`: one less free parameter, and the singleton is itself an
argument for the small K.

Second fix: `fit()` now reproduces `vessel_niche_cluster.autotuned_kmeans` exactly (10 restarts
at random_state seed+i, keep the lowest inertia) instead of `KMeans(n_init=10)`. With that,
the refit labels are IDENTICAL to the stored folders -- ARI = 1.000 at K = 3, 4 and 9 (it was
0.95 / 0.96 / 0.89 before). Every figure in figures_niche_robustness now describes the
clustering that is actually on disk.

Unclipped results: min Jaccard K2 0.953, K3 0.905, K4 0.652, K5 0.166, K6 0.218, K7 0.376,
K8 0.428, K9 0.529 -> still K = 3 (clusters n=141 J=0.96, n=110 J=0.98, n=49 J=0.91).
Consensus within/between: K3 0.95/0.02, falling to 0.72-0.75 / 0.05-0.08 by K7-K9.
The singleton P17_AIS__6 now shows up in the small-cluster table with Jaccard 0.65 (K=4)
rising to 1.000 (K=9) -- the textbook "high Jaccard but tiny" case, resolved by the next two
columns: n_sections = 1, top_section_frac = 1.00, 45 % Airway_epi. Artefact, not a rare niche.

New figure: `niche_motif_best_match` -- rows K = 2..9, columns the four motifs, each cell the
Pearson r of whichever niche matches that motif best at that K, annotated with which niche and
its size (table: niche_motif_best_match.csv). It reads as: the SMC-stabilization motif always
finds a ~105-110 vessel niche at r ~ 0.51-0.56 regardless of K; Tumor vasculature peaks at
K = 3 (r = 0.41, the 49-vessel T-high niche) and degrades as that niche is split; Healthy
alveolar improves slightly with K (0.26 -> 0.46) by carving out a 26-31 vessel subset;
Vascular homeostasis never exceeds r = 0.21 at any K, because m23 is ON in 290 of 300 vessels
and so has almost no variance to explain.

**Slim variant added**: `niche_motif_best_match_slim` is the same K x motif best-match
heatmap with only the r and the niche id in each cell (no cluster size); the annotated
version `niche_motif_best_match` is kept.
