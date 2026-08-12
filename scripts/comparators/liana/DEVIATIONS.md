# LIANA+ — deviations, and a corrected method choice

## ⚠️ We initially used the wrong branch of the authors' own decision tree

LIANA+'s README ships a decision tree. Its spatial half reads:

```
Spatially-resolved -> Resolution?
  Single-cell -> Interaction scoring -> INFLOW SCORE
              -> Interaction scoring -> Standard LR methods (spatially-constrained)
              -> Spatial co-occurrence -> LRIC
              -> Unsupervised        -> Communication Programs (Inflow + MOFA-Flex)
  Spot-based  -> Bivariate -> Local / Global Bivariate Metrics
              -> Unsupervised -> MISTy
```

**`li.mt.bivariate` sits under the *spot-based* branch. Xenium is single-cell.** The
single-cell route is `li.mt.inflow` (or spatially-constrained standard methods, or `li.mt.lric`).

The package makes the same point in its defaults: `_inflow.py` sets **`nz_prop = 0.001`** against
`bivariate`'s **`0.05`** — a 50× lower threshold. The authors addressed single-cell sparsity by
**writing a different method**, not by asking users to lower a threshold. Our earlier binomial
spot→cell derivation of `nz_prop = 0.02` for `bivariate` was therefore solving a problem that
should not have been posed; reading the decision tree first would have avoided it.

**Both runs are kept.** They answer different questions and neither is discarded:

| Run | Path | Rationale |
|---|---|---|
| `cellchatdb2/` | `li.mt.bivariate`, `nz_prop=0.02` | the spot-based branch, run on single-cell data — retained as the comparison point and because the *bivariate* tutorial is the only place the authors demonstrate NMF |
| `cellchatdb2_inflow/` | `li.mt.inflow`, `nz_prop=0.001` | the branch the authors' tree actually prescribes for this resolution |

**Since 2026-08-04 four more of the tree's leaves have been exercised**, so the "wrong branch"
framing above now covers a minority of what is on disk:

| Run | Decision-tree leaf | Path | Cost |
|---|---|---|---|
| `cellchatdb2_morans/` | Bivariate, second local metric | `li.mt.bivariate(local_name='morans')` | 3.6 min |
| `mofaflex_inflow/` | **Single-cell → Unsupervised → Communication Programs** | `run_mofaflex.py`, `inflow_mofaflex.ipynb` | 76.0 min |
| `lric_percore/` | **Single-cell → Spatial co-occurrence → LRIC** | `run_lric.py`, `li.mt.lric` + `li.mt.cross_pcf`, **13 per-punch runs** | 2.6 min |
| `misty/linear_fullslide/` | **Spot-based → Unsupervised → MISTy** | `run_misty.py`, `misty.ipynb` LR-MISTy config | 2.9 min |
| `default*/`, `nmf_*_default/` | both branches at LIANA's own `consensus` resource | `run_default_tier.sh` (now passes `--k-max 21` explicitly) | 5.6 + 1.3 + 3.4 + 3.3 + 2.2 min — the two NMF terms are the `--k-max 21` **refits**; this row read `+ 1.8 + 1.4` until 2026-08-04, those being the discarded `k_range` 1..10 fits |
| `factor_annotation/` | annotation only, no re-fit | `annotate_factors.py` | 0.22 min |

`li.mt.lric` is now the **only** branch here that resolves *direction* as an argument, and the only
one that supports the mGAM→MES-like arm as directional at the punch replicate unit. Note also that
neither LRIC nor LR-MISTy inherits the 13.1454 µm bandwidth, so **D1 below does not apply to
either** — see the tutorial-deviation table.

## Consequences, measured

| | `bivariate` | `inflow` |
|---|---|---|
| features | **131** | **4,608** (35.2×) |
| unique LR interactions | 131 | **633** (4.83×) |
| feature identity | LR pair | **sender cell type × LR pair** (`<cell_type>^<lig>^<rec>`) |
| senders | — | **9** |
| sender × LR expansion | — | **RAGGED, not a complete grid** — 250–602 LR pairs per sender |
| non-negative | — | **yes**, 0 negatives of **461,675,520** (= 100,190 cells × 4,608 features) |
| sparsity | — | **99.4538%** zeros |

> **Corrected 2026-08-04.** This table previously read "535 unique LR interactions (4.1×)" and
> "0 negatives of 482,414,850". Both were wrong. Recomputed by splitting the 4,608 feature names
> in `cellchatdb2_inflow/data/inflow_scores.npz` on the first `^`: **9 senders, 633 distinct LR
> pairs**, so the expansion over bivariate's 131 is **4.83×**. The 535 has no basis in this run —
> it appears to have been carried over from the CellChatDB Cell-Cell-Contact row count in
> `../METHODS.md`. The cell count is 100,190 (not 100,197: `run_inflow.py` applies
> `sc.pp.filter_cells(min_genes=10)`, which `run_liana.py` does not), so the matrix is
> 100,190 × 4,608 = 461,675,520 entries. 99.4538% zeros is `run_manifest.json`'s own
> `pct_zeros`, which is authoritative.

**The expansion is ragged, and that is structural, not noise.** Per-sender LR-pair counts are
MES-like 602, NPC-like 601, OPC-like 596, AC-like 590, Glial-Neuronal 559, Vascular 520, mGAM 508,
non-mGAM 382, Lymphoid 250. A complete 633 × 9 grid would be 5,697 features; only 4,608 exist,
because `nz_prop` is applied within each sender's own expression profile. So **feature count per
sender is partly an abundance/expression-breadth statistic**, not purely a signalling one — the
rarest cell type (Lymphoid) contributes the fewest features. Never describe the inflow feature
space as "512 interactions × 9 senders" or as any other product; it is not a product.

`inflow` carries **sender cell-type identity in the feature itself** (`C_{j,s}` is a hard
indicator in the score), which `bivariate` does not. That is a qualitative difference, not just
a count: an inflow factor can be read as "cell type X signalling into these locations".

## 🔴 Deviations from the BENCHMARK CONTRACT — OPEN, awaiting sign-off

**These are a different category from the tutorial deviations below, and must not be merged into
that table.** The table below records where we depart from *LIANA's authors*. This section records
where we depart from **our own benchmark contract**, `.claude/skills/comparator-benchmark/SKILL.md`
— i.e. from the rules that make the seven comparators comparable *to each other and to ALARMIST*.
Both were made without being flagged, and both are still **OPEN**. Recorded 2026-08-04; no
parameter was changed and nothing was re-run in recording them.

### D1 — the kernel bandwidth is derived from ALARMIST's patch size

`SKILL.md:45-46` reads verbatim:

> **Keep each method's own neighborhood/kernel definition at its default.** Do NOT harmonize
> spatial scale across methods, and do not match it to ALARMIST's patch size.

`SKILL.md:105` additionally lists "any patch/kernel/neighborhood scale parameter the tutorial does
not pin down" under STOP-and-ask. The *inflow* tutorial — our branch — pins nothing.

What we did is the opposite of the rule: `run_liana.py:30-31` documents 13.1454 µm as `s/3.804`
with `s = 50 µm` = **ALARMIST's patch edge**, via equal-area correspondence. Every
`run_manifest.json` in the bivariate / inflow trees — **both tiers**, **five** of them
(`cellchatdb2`, `cellchatdb2_inflow`, `cellchatdb2_morans`, `default`, `default_inflow`) — carries
`"bandwidth": 13.1454`. *(Corrected 2026-08-04: this read "**eight** of them" and included the NMF
trees. Verified by parsing all 16 `run_manifest.json` under `results/comparators/liana/GBM/`: the
four `nmf_*` manifests carry **no** `bandwidth` key at all, and the eight `run_manifest.json` that
contain the literal string `13.1454` include `lric_percore` and both `misty` manifests, whose
`params_provenance` says the opposite — "NOT the 13.1454 um bandwidth used by this repo's
bivariate/inflow LIANA runs" / "the 13.1454 um bandwidth is irrelevant here". The old count was a
grep that counted the negative statements.)*

**Scope, added 2026-08-04: D1 applies to the bivariate, inflow and NMF branches only.**
`li.mt.lric` / `li.mt.cross_pcf` consume no connectivity graph at all (they build their own
`cKDTree`), and LR-MISTy was run at the **tutorial's** `bandwidth=200`. Neither inherits the
ALARMIST-derived number, so neither is in violation.

Why this is not cosmetic: `../METHODS.md`'s *Sensitivity: the number of communication programs
moves with the bandwidth* table measures that this parameter alone moves the inflow NMF rank from
**11 to 7**. The headline "LIANA finds 7 programs" is therefore a function of
a number imported from the method it is being compared against. The choice is defensible on its
own terms — equal-area correspondence is a principled mapping, and it is transparently derived —
but *defensible* is not *signed off*, and it is exactly the harmonization the contract forbids.

**Status: OPEN.** Either (a) get explicit sign-off to keep it, recording that the comparison is
scale-harmonized by design, or (b) add a bandwidth arm at a LIANA-native default. **Do not silently
change the parameter** — every result on disk was produced at 13.1454.

### D1, update 2026-08-06 — the bandwidth STAYS at 13.1454 µm by user decision, and D1 STAYS OPEN

The user decided to keep σ = 13.1454 µm / support R = 28.2096 µm. **That decision does not close
D1**: the value is still derived from ALARMIST's patch size, which `SKILL.md:45-46` forbids, and it
is still not signed off as a deliberate scale-harmonization. **Do not mark this resolved.** What
changed is that the tutorial's own exploration step has now been run properly
(`choose_bandwidth.py` → `results/comparators/liana/GBM/bandwidth_choice/`, no re-fit), so the
choice is at least documented rather than merely asserted. Four measured findings:

**1. `li.ut.query_bandwidth` returns `ceil(MEDIAN) - 1`, not the mean.**
`liana/utils/query_bandwidth.py:71-72` computes `avg_nn = np.ceil(np.median(num_neighbors))` and
then writes `avg_nn - 1` — despite the variable name. This reconciles three numbers that otherwise disagree: our own BallTree **mean** of 14.6 at
R = 28.21, the value **13** read off the plotted curve, and the **median** of 14. `ceil(14) - 1 = 13`.
Verified against `bandwidth_choice/data/query_bandwidth_tutorial_5_35.csv`, whose row at
x = 28.077 reads 13.

**2. The tutorial mixes two different quantities on one axis.**
`query_bandwidth`'s x-axis is a **hard query radius** (`BallTree.query_radius`), while
`li.ut.spatial_neighbors(bandwidth=)` is a **gaussian σ** truncated at `cutoff`, with reach
`R = σ·sqrt(−2 ln cutoff) = 2.145966·σ`. A value read off the curve is an **R** and must be
**divided by 2.146** before being passed as `bandwidth=`; passing it straight through inflates the
neighbourhood **area by 4.6×**. This is an inconsistency in the tutorial itself, not in our port.
> ⚠️ **Retracted 2026-08-07.** The 2026-08-06 pass wrote here that *"the pre-existing
> `cellchatdb2_inflow/plots/global/bandwidth_query.png` draws its vertical guide at the σ (13.1454)
> on a radius axis and is therefore wrong"*, and carried it into the open-issues table below. **That
> is false, and it was never checked against the code or the image.**
> `run_inflow_downstream.py:163-167` reads
> `R = a.bandwidth * np.sqrt(-2 * np.log(a.cutoff))` → `geom_vline(xintercept=float(R))`, i.e. the
> guide is at the **support radius 28.2096 µm** — the correct scale — labelled
> `support radius = 28.2 um (gaussian sigma = 13.1454)`. The PNG on disk agrees: the dashed line
> sits at ≈28 µm and crosses the curve at ≈14 neighbours. **No figure inherits the tutorial's unit
> confusion; our port got this right the first time.** The only defect is cosmetic — the rotated
> annotation is anchored at `y = neighbours.max()` with `va="bottom"` and is clipped above the
> panel, leaving `su` visible. Kept on the record because the wrong claim was plausible enough to
> survive three documents. The two `bandwidth_choice/figures/query_bandwidth_{tutorial_5_35,extended_5_120}.png`
> draw the same guide, unclipped and with a second guide at the tutorial's σ = 27.

**3. The exploration step does not select a bandwidth on this tissue.** Mean neighbours vs hard
radius is smooth and monotonic over 5–120 µm with **no plateau, elbow or inflection**
(`bandwidth_choice/data/neighbours_vs_radius.csv`, interpolated at the quoted radii):

| R (µm) | 10 | 20 | 28.21 | 40 | 57.94 | 70 | 120 |
|---|---|---|---|---|---|---|---|
| mean neighbours | 1.7 | 7.5 | **14.6** | 28.9 | **59.2** | 85.2 | 240.1 |

So the technical half of the tutorial's criterion ("too wide blurs, too narrow misses") has nothing
to grip on here.

**4. The biological half cannot be satisfied either — there is no characteristic length scale.**
LRIC `g(r)` over **all 1,088 resolvable LR pairs × 13 punches = 11,795 pair-punch observations**
(`lric_percore/punches/*/lric_agnostic_matrix.csv.gz`), median `g` per annulus:

| bin (µm) | 0–50 | 50–75 | 75–100 | 100–125 | 125–150 | 150–175 | 175–200 | 200–225 |
|---|---|---|---|---|---|---|---|---|
| median g | 1.459 | 1.419 | 1.403 | 1.394 | 1.401 | 1.393 | 1.381 | 1.395 |

A **4.4 % decline across the entire 225 µm range** — co-occurrence is essentially flat, so the
tutorial's biological criterion ("reflect the typical range of molecular signaling") **does not
constrain the choice on this tissue**. *(An earlier version of this claim rested on only the **2**
required LR pairs, i.e. `lric_percore/combined/aggregate_per_bin.csv`, 32 rows. It now rests on
1,088 pairs and must be stated with that denominator.)*

**What the evidence does bound.** Median nearest-neighbour distance is **7.86 µm** (IQR 6.30–10.47,
p95 19.97), so a strictly juxtacrine reach is ~8–10 µm (σ ≈ 3.7–4.9) — a **floor**, not a choice.
LIANA's `max_neighbours=100` default gives an upper bound: at R = 70 the mean is already 85.2. That
leaves a defensible window of **R ∈ [20, 58] µm, σ ∈ [9.3, 27]**; 13.1454 sits inside it, which
makes the value *defensible* but still not *derived from LIANA*.

**The alternative not taken:** the inflow tutorial's own value is **σ = 27 µm → R = 57.94 µm**
(59.2 neighbours). Ours is **half that spatial scale and a quarter of the area**. Choosing it would
have been the LIANA-native option and would have closed D1.

⚠️ **The bandwidth and the QC attrition are coupled**, so this is not an isolated parameter: at
R = 57.94 µm mGAM reachability rises **0.319 → 0.712** (measured, `bandwidth_choice.json`). Since
the global `nonzero_fraction` cut is bounded above by reachability, a wider kernel would have
reduced the view attrition recorded in the reachability section above — plausibly enough to carry
the mGAM view over `min_features=25` without deviating from the authors' QC at all. **That last
step is an inference, not a measurement: no fit at σ = 27 has been run.**

### D2 — no native multi-sample / differential mode was used

`SKILL.md:47-49` requires: *"Use the method's native multi-sample / differential mode when it has
one: GBM → split by `obs['grade']` (high vs low; 13 TMA cores in `obs['tma_id']` are the units)."*

`../METHODS.md`'s *Multi-sample / differential mode* text used to mark that requirement
**satisfied**. It is not, and both files now say so.

What was actually run is `li.mt.compute_global_specificity(groupby='grade')`. Reading the installed
implementation (`liana/method/sp/_compute_global_specificity.py`), that is a one-sided **per-group
specificity** test which permutes labels across **cells** — it is not a contrast. The output bears
this out: `cellchatdb2_inflow/data/region_global_interactions.csv` is **9,216 rows = 4,608 × 2**
(one block for `high`, one for `low`), with columns `source, ligand_complex, receptor_complex,
target, lr_mean, pval` and **no contrast column at all** — no effect size, no direction, no
high-vs-low statistic.

The cores never enter it: `grep -n tma_id scripts/comparators/liana/run_inflow_downstream.py`
returns **nothing**. `obs['tma_id']` is read only by `run_nmf.py`'s ≥5/13 cross-punch presence
filter (see the last section of this file), which is a reproducibility filter, not a test.

**Consequence, and it must be stated wherever these p-values are quoted:** the 5,417 rows at
p<0.05 in `global_interactions.csv` (and everything derived from them, including
`lr_ranking_by_lr_mean.csv` and the plotted top-6) are **cell-level permutation p-values at
n = 100,190**, pseudoreplicated relative to the **13 TMA cores** that are the actual experimental
units. They are a ranking statistic, not evidence of a grade difference. This is the standard
pseudoreplication trap that `spatial-workflow` owns; LIANA is not exempt from it.

`grep -n by_sample scripts/comparators/liana/*.py` returns no matches — none of `by_sample`,
`dotplot_by_sample`, `lrs_to_views` or `to_tensor_c2c` was ever called.

**Status: still OPEN as a contract item, but no longer unaddressed.** LIANA has no native spatial
differential mode, so `SKILL.md:47-49`'s requirement cannot be met on its own terms — and
`SKILL.md:49` explicitly says *"If the method has no multi-sample mode, say so — do not hand-roll
one."* So the honest reading is: **the requirement is inapplicable, and must be marked ❌/N-A, not
✅.**

**A punch-level test has nonetheless been run** (`analyse_existing.py` → `results/comparators/liana/GBM/nmf_inflow/punch_level/`,
2026-08-04, 0.12 wall-min, no re-fit). It aggregates `nmf_WH.npz`'s `W` by `tma_id` and runs a
two-sided Mann-Whitney over the **7 high vs 6 low punches**, BH-corrected. It is recorded as an
*additional* analysis, not as the method's own mode. Result, verified from
`punch_level/data/punch_factor_tests.csv` and `punch_requiredLR_tests.csv`:

| test | n signif. at BH q<0.05 | smallest raw p | smallest q |
|---|---|---|---|
| 7 inflow NMF factors vs grade | **0 / 7** | 0.013986 (Factor4, log2FC +1.51) | 0.0816 |
| 20 required-LR features vs grade (9 senders × 2 LRs + 2 aggregates) | **0 / 20** | 0.013986 (`MES-like^GRN^SORT1`, log2FC +2.62; `AC-like^GRN^SORT1`) | 0.1399 |

**Nothing survives correction at the correct replicate unit.** Note the hard floor: a two-sided
rank test on 7 vs 6 units cannot go below **p = 0.001166**, so this design has little power and
the null result is weak evidence of absence, not evidence of no effect. Contrast this with the
5,417 "significant" cell-level rows — that gap *is* the pseudoreplication, made concrete.

## Deviations from the tutorials

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| bandwidth | inflow: 27 µm (MERFISH mouse brain); bivariate: 200 (Visium px) | **13.1454 µm** | Set by equal-area correspondence to an s x s patch: support = s/sqrt(pi) = 28.2 µm, bandwidth = s/3.804 = 13.1454 µm. The earlier 18.75 (bivariate tutorial's ~6-neighbour "first ring" rule) is VOID -- first ring is a hexagonal-lattice concept with no referent on irregularly packed single cells. |
| resource | `consensus` (+ mouse ortholog translation) | **CellChatDB v2**, no translation | Our data is human; `resource` accepts a `['ligand','receptor']` frame and LIANA joins complex subunits with `_`, so the handover is direct and lossless. |
| SVG gene filter | Moran's I, FDR<0.05 & I>0.01 | applied — **no-op** | 5,097 → 5,097 genes. **The no-op is carried entirely by the `I > 0.01` half** (min I on the panel is **0.159**, 16× the threshold): the FDR half is an *underflow artifact* — all 5,097 `pval_norm` and `pval_norm_fdr_bh` values in `cellchatdb2_inflow/data/gene_moranI.csv` are **exactly 0.0** at n = 100,190, so that criterion is vacuous and cannot exclude anything at this sample size. Net: every gene on the Xenium panel is spatially autocorrelated well above the effect-size threshold — expected for a targeted panel on structured tissue — but do not cite the FDR as evidence. |
| SVI interaction filter | marked *optional*, Moran's I, FDR≤0.05 & I>0.01 | **not applied** — and measured to be a **no-op** | Left off so the feature space is not pre-selected for spatial structure before factorisation, which would confound the NMF comparison. Exposed as `--svi-filter`. **Measured post-hoc (2026-08-04) so the deviation is numbered, not argued: 4,608 → 4,608 interactions** — applying it would have changed nothing. Same structure as the SVG row above: **the `I > 0.01` half carries it** (min I is **0.0578**, `NPC-like^WNT4^FZD6_LRP6`; max **0.9962**, `Glial-Neuronal^CXCL2^CXCR1`; median 0.320), and the FDR half is again an *underflow artifact* — all 4,608 `pval_norm` and `pval_norm_fdr_bh` are **exactly 0.0** at n = 100,190, so it is vacuous. `sq.gr.spatial_autocorr(mode='moran', use_raw=False)` on `cellchatdb2_inflow/data/inflow_lrdata.h5ad`, its stored `spatial_connectivities` (no rebuild), 19.4 s / 1.15 GB peak RSS; two replicates agree to 1.1e-16. Table: `cellchatdb2_inflow/data/interaction_moranI.csv`. |
| control probes | n/a | 21 `Intergenic_Region_*` dropped | `_` collides with LIANA's complex-subunit separator; LIANA warns rather than errors. |
| core split | single section | **whole slide** (bivariate / inflow only) | Verified: effective radius 28.2 µm vs a measured 222.9 µm minimum inter-core cell–cell distance → **zero** cross-core pairs. ⚠️ **This clearance is branch-specific and does not transfer.** LR-MISTy at the tutorial's `bandwidth=200, cutoff=0.01` has a 607.0 µm nominal support and **does** produce cross-core edges (2,520 of 10,119,190, 0.0249%, touching 232 cells); LRIC/cross-PCF is run **per punch** for a different reason again (density normalisation). See the two rows below. |
| local metric | `cosine` (bivariate tutorial) | **also run with `morans`** | `cellchatdb2_morans/`, 2026-08-04, `--local-name morans` the only change. Not a deviation so much as a second arm: both are in `li.mt.bivariate.show_functions()`. Global Moran's R is bit-identical across the two; the *local* scores agree at median per-pair r = **0.195**. Moran's R is **not NMF-admissible** (33.4% negatives, 0 zeros), which is why `cosine` remains the NMF input. |
| LR resource tier | `consensus` | **both tiers now run** | `run_default_tier.sh` (2026-08-04) adds `default/`, `default_inflow/`, `nmf_bivariate_default/`, `nmf_inflow_default/` at LIANA's own `consensus` (4,624 pairs, 36.0% overlap with CellChatDB v2). The `resource` row above therefore describes the `cellchatdb2` tier only. Measured: the resource changes **which** pairs are tested and the ranking denominator, but **not one score** — see the new section below. |
| LRIC / cross-PCF spatial support | n/a — no kernel | **13 per-punch runs**, annuli to 225 µm | `li.mt.lric` / `li.mt.cross_pcf` normalise by `n_points / bbox area`, and the TMA's global bbox is only **42.2306%** occupied → whole-slide density understated **2.3679×**. **Measured** with a pooled control: LRIC g(r) inflated **3.4284×**, cross-PCF **2.9445×**, and the LRIC/cross-PCF ratio distorted from ~1.0 to 1.2–2.3, i.e. whole-slide would have produced a false positive. `li.ut.spatial_neighbors` is deliberately **not** called (neither function consumes a graph), so the 13.1454 µm bandwidth does not enter — **CD-1 does not apply to this branch**. |
| LR-MISTy bandwidth | 200 (Visium tutorial) | **200 — unchanged** | This branch keeps the tutorial's value rather than inheriting 13.1454 µm, so **CD-1 does not apply to it either**. Consequence: nominal support 200·√(−2·ln 0.01) = **607.0 µm**, above the 222.9 µm inter-core floor. What limits the damage is `max_neighbours=100`, which `lrMistyData` does not expose — degree median = max = 100, cap binding for **99.7%** of cells, so the extra view is in practice a 100-NN neighbourhood, not a 200 µm kernel. |
| LR-MISTy HVG pre-step | `adata[:, hvg]` | **omitted** | Measured (`--hvg 1 --n-top-genes 1521 --construct-only 1`): intra **100,190 × 82** / extra **100,190 × 13** with it, against **382 / 37** without — 78% of receptor targets and 65% of ligand predictors discarded; on one 5,363-cell punch the extra view collapses to **1** predictor. The tutorial frames HVG as *"for the sake of computational speed"* on a genome-wide Visium slide; this is a 5,119-gene targeted panel and the full-panel fit takes **2.7 min**. |
| LR-MISTy secondary model | `RandomForestModel` | **not run** | Measured on the full slide, not extrapolated: **87.7 s/target** (4-target probe, 5.85 min) → **9.31 h** for 382 targets, ~9× the budget and 2.8× worse than the punch-4 scaling law predicted. Kept as `misty/rf_rate_probe/`. |
| MOFA-Flex feature QC | cell 19 `nonzero_fraction > 0.01`; cell 23 `min_features=25` | **applied as written**, plus a `nzf > 0.001` sensitivity arm | Applying it as written **deletes both arms of ALARMIST motif 1** — 4,608 → 447 features removes every `ANXA1^FPR1` (max nzf **0.009422**), and `min_features=25` drops the whole mGAM view at exactly **24** features. The sensitivity arm exists because a filter that removes the quantity under test cannot be the only run. |
| MOFA-Flex feature QC — **reachability normalisation** *(new 2026-08-06)* | cell 19 applies **one global** `nonzero_fraction > 0.01` to every feature, whatever its sender | **third arm, opt-in**: `--nzf-mode reachability` divides each feature's `nonzero_fraction` by its **sender's reachability** before the cut. **The default stays `global`**, so the tutorial-faithful run remains reproducible. | A `<sender>^<lig>^<rec>` feature can only be non-zero in a cell that has that sender inside the kernel support, so its `nonzero_fraction` is **bounded above by the sender's reachability**. One global threshold therefore ranks senders by abundance, not by signalling: surviving-feature count vs sender abundance is **Spearman ρ = 0.917** (p = 5.1e-4, n = 9). Effect: 4,608 → **779** features, **9 / 9** views kept (was 6 / 9). Full evidence in the section below. |

## The tutorial's global `nonzero_fraction` cut is an abundance filter — measured, and deviated from 2026-08-06

**This is a tutorial deviation, not a contract one** — it belongs to the table above, and the
default is left at the tutorial's behaviour (`--nzf-mode global`) so nothing already on disk is
orphaned. `results/comparators/liana/GBM/mofaflex_inflow_reachnorm/` is the third arm.

**The arithmetic.** An inflow feature `<sender>^<lig>^<rec>` scores the *receiving* cell, and it can
only be non-zero for a cell that has at least one cell of that **sender** type inside the kernel
support. So

```
nonzero_fraction(feature)  <=  reach(sender) = P(>=1 neighbour of type s within R)
```

and `nonzero_fraction > 0.01` is not one criterion applied nine times — it is nine different
criteria. Expressed as a fraction of the cells that *could* carry the feature, the same 0.01 cut is
**1.3 %** for NPC-like and **39 %** for Lymphoid.

**Lymphoid cannot pass, for arithmetic and not biological reasons.** At R = 28.2096 µm its
reachability is **0.0254**, and its highest observed `nonzero_fraction` over all 250 of its features
is **0.004332** — below 0.01 by construction. **Zero** Lymphoid features can survive the tutorial
cut at any expression level whatsoever. **Vascular is penalised for clustering, not for rarity:**
it is **3.22 %** of cells but only **15.7 %** reachable, because vessels are spatially aggregated.

| sender | share of cells | reachability at R = 28.2096 µm | features kept, tutorial `global` | features kept, `reachability` |
|---|---|---|---|---|
| NPC-like | 22.78% | 0.761 | 98 | **114** |
| OPC-like | 18.91% | 0.748 | 78 | **92** |
| AC-like | 18.12% | 0.750 | 92 | **111** |
| MES-like | 17.49% | 0.727 | 97 | **122** |
| Glial-Neuronal | 12.27% | 0.379 | 28 | **63** |
| mGAM | 4.02% | 0.319 | **24** — one short of `min_features=25` | **66** |
| Vascular | 3.22% | 0.157 | 26 | **138** (5.3×) |
| non-mGAM | 2.89% | 0.205 | 4 | **34** |
| Lymphoid | 0.30% | 0.025 | **0** — arithmetically impossible | **39** |
| **total** | | | **447**, 6 views | **779**, 9 views |

Views the tutorial recipe drops: **mGAM** (24 features, one short of `min_features=25`),
**non-mGAM** (4), **Lymphoid** (0). Reachabilities read off
`results/comparators/liana/GBM/bandwidth_choice/bandwidth_choice.json`; kept-feature counts off the
two `run_manifest.json`; the ρ and the Lymphoid maximum recomputed from
`cellchatdb2_inflow/data/inflow_lrdata.h5ad`.

**The fix** (`run_mofaflex.py:316-357`) is `nzf_norm = nonzero_fraction / reach[sender]`, keep
`> 0.01`. It is a **deviation from the authors**, recorded as such, and it is opt-in.

**The reachability-normalised fit.** `--nzf-mode reachability --nzf-threshold 0.01`, everything else
identical to the primary arm: CPU, K = 20, batch 2048, lr 0.005, patience 50, seed 0. The **SVI
filter is again a measured no-op (779 → 779)**. Converged after **294** of 1000 epochs; fit
**40.1 min**, wall **41.5 min**, peak RSS **6.78 GB**; **19 / 20** factors active. *(As with the
other two fits, the epoch count is **not** in `run_manifest.json` — it lives only in
`logs/comparators/liana-mofaflex-reachnorm.log`. Same provenance gap as recorded for the primary
and sensitivity arms.)*

| view (SENDER) | R² | note |
|---|---|---|
| Glial-Neuronal | 0.5058 | |
| **mGAM** | **0.4453** | **second highest — the view the tutorial recipe deletes** |
| NPC-like | 0.3931 | |
| MES-like | 0.3542 | |
| OPC-like | 0.3254 | |
| Vascular | 0.2708 | |
| AC-like | 0.2629 | |
| non-mGAM | 0.1397 | |
| Lymphoid | **0.0013** | ⚠️ admitted, but the model explains essentially nothing here |

⚠️ **The Lymphoid caveat must travel with this result.** Lymphoid clears the normalised criterion
with 39 features and then reaches **R² = 0.0013** — it is *admitted but uninformative*, and its
features rest on a **2.5 %** reachable cell base. **Do not interpret Lymphoid factors.** The
deviation implemented is the normalisation *alone*; **no absolute cell-count or reachability floor
was added**, so a view can now enter the model on a base too small to fit. That is a known,
deliberate limitation of this arm, not an oversight.

**Restoring mGAM did not restore motif 1.** Punch-level grade test on the 19 active factors:
**0 / 19** at BH q < 0.05, smallest raw p **0.0221** (Factor 18) → q **0.332**, against the same
7-vs-6 floor of p = 0.0011655. And the two motif-1 arms **still peak on different factors** — see
*Where the two motif-1 arms actually rank* below. The tutorial's QC is a real defect in the
comparison, but it is **not** what prevents LIANA from seeing the loop — the **unit of analysis**
is, and that is a property of `li.mt.inflow` no feature filter can undo. See *Why LIANA cannot
recover the mGAM loop*.

## ⚠️ LR-database provenance — these runs and `results/GBM/` used DIFFERENT CellChatDB exports

The `resource` row above says "CellChatDB v2". That is not specific enough, and the ambiguity is a
join trap. **CLAUDE.md records that the 2026-07-28 re-export changed 1,120 of 3,218 human LR keys**
by reordering complex subunits (`TGFBR2_TGFBR1` → `TGFBR1_TGFBR2`, `RAMP2_CALCR` → `CALCR_RAMP2`)
and by turning some ligands into homo-complexes.

**Which export each side used — verified 2026-08-04 by subunit order and by key matching:**

| | export used | evidence |
|---|---|---|
| the LIANA runs (all four dirs) | **current / re-exported** `data/LRdatabase/CellChatDBv2.0.human.csv` | `cellchatdb2_inflow/data/global_interactions.csv` has **450 ×** `TGFBR1_TGFBR2` and **0 ×** `TGFBR2_TGFBR1` |
| `results/GBM/` (the ALARMIST run being compared against) | **old** `CellChatDBv2.0.human.old.csv` | `results/GBM/patch_lri_columns.csv` has **210 ×** `TGFBR2_TGFBR1` and **0** reversed |

Key-matching confirms it independently. Taking the run's LR keys as `(ligand, receptor)` tuples:

| key set | vs current export | vs `.old.csv` |
|---|---|---|
| inflow, 633 distinct pairs | **633 / 633** | 496 / 633 |
| bivariate, 131 pairs | **131 / 131** | 115 / 131 |

**Neither run is wrong** — CLAUDE.md directs fresh runs at the current file, and `results/GBM/`
predates the re-export. The trap is downstream: **a raw-string LR join between the LIANA outputs
and `results/GBM/` silently drops heteromeric keys**, with no error and no warning — they simply
fail to match and vanish from the intersection.

**Remedy — canonicalise both sides before joining**, applying it to ligand and receptor
*separately* (never to the whole `lig^rec` string):

```python
canon = lambda s: "_".join(sorted(s.split("_")))
key = (canon(ligand), canon(receptor))
```

**The two required LRs are unaffected**: `GRN^SORT1` and `ANXA1^FPR1` are single-subunit on both
sides and identical in both files. So motif-1 conclusions do not depend on this — but any
*panel-wide* overlap statistic between LIANA and ALARMIST does, and would be biased downward
without the fix.

## ⚠️ NMF on inflow is OUR composition, not an author-demonstrated path

`inflow_score.ipynb` contains **no NMF section**. NMF (`li.multi.nmf`) appears only in
`bivariate.ipynb`. Applying it to inflow output is therefore our own construction and must not be
presented as following the authors' workflow.

Moreover, the decision tree shows the authors' *own* answer for unsupervised decomposition at
single-cell resolution: **"Communication Programs — Inflow + MOFA-Flex"**, demonstrated in
`inflow_mofaflex.ipynb`. So for this exact goal there is an author-sanctioned route, and it is a
**different factorisation** (MOFA-Flex, not NMF).

We run NMF on both anyway, deliberately:

- **NMF on `bivariate`** — the tutorial-sanctioned composition. Evidence that we exercised the
  authors' demonstrated path.
- **NMF on `inflow`** — the resolution-appropriate input, factorised with the only decomposition
  the package demonstrates for CCC output.

Neither is the authors' prescribed unsupervised route. ~~**MOFA-Flex on inflow remains the
untested, author-sanctioned alternative**~~

**✅ CLOSED 2026-08-04 — MOFA-Flex on inflow has been run.** `run_mofaflex.py` →
`results/comparators/liana/GBM/mofaflex_inflow/`, following `inflow_mofaflex.ipynb` cell by cell:
CPU, K = 20, batch 2048, lr 0.005, patience 50, seed 0, early-stopped at **632 / 1000** epochs;
fit **70.5 min**, wall **76.0 min**, peak RSS **4.11 GB**; **17 / 20** factors active. **The NMF runs
remain our own composition** — that part of this section stands unchanged.

**What running it revealed is a finding about the comparator, and it is the most reportable single
result of this pass: the authors' own QC deletes both arms of ALARMIST motif 1.**

| tutorial cell | filter | effect here |
|---|---|---|
| 19 | `nonzero_fraction > 0.01` | 4,608 → **447** features; removes **every** `ANXA1^FPR1` feature (max nzf across 9 senders = **0.009422**, MES-like; mGAM 0.004911) |
| 23 | `lrdata_to_mudata(min_features=25)` | drops the **whole mGAM view**, which retained exactly **24** features — one short — taking `mGAM^GRN^SORT1` with it despite its nzf of 0.031770 |

So **LIANA+, run exactly as its authors demonstrate, structurally cannot see either arm of motif 1
on this dataset.** A sensitivity fit at `nzf > 0.001` (1,550 → **1,541** features, 8 views,
early-stopped **199 / 1000**) was therefore added, and there the two arms peak on **different**
factors — `mGAM^GRN^SORT1` on Factor 19 (**+0.651**, rank 73/1541), `MES-like^ANXA1^FPR1` on
Factor 7 (**−0.278**, rank 136/1541). Factor 7 carries three of the four loop features with
concordant sign and is myeloid-anchored (largest weight `mGAM^C3^C3AR1` = **−2.892**), but
⚠️ **Factor 7 is only the *eighth* strongest factor for `mGAM^GRN^SORT1`** (|w| order: F19 0.651,
F15 0.590, F1 0.558, F16 0.556, F17 0.531, F9 0.451, F6 0.413, F7 0.398) — the concordant signs are
real, the preference is not.

Punch-level grade test: **0 of 20** factors significant in either fit (primary smallest raw
p = 0.013986 → q = 0.217949, Factor 18; sensitivity p = 0.008159 → q = 0.163170, Factor 11), against
the same 7-vs-6 floor of **p = 0.0011655**.

**Version gap, recorded rather than guessed:** `inflow_mofaflex.ipynb` cell 5 states it targets the
**MOFA-Flex 0.2.0 API, "not yet released on PyPI"**; the installed build is
`0.1.0.post2.dev179+g9792b435f` from git main. Every symbol the notebook uses exists in the
installed build with matching argument names.

**Provenance gap — and it widened on 2026-08-04.** It was originally recorded here as "`run_manifest.json`
records `fit_seconds` but **not the stopping epoch** — 632 / 199 live only in
`logs/mofaflex_{primary,sensitivity}.log`". It no longer records `fit_seconds` either: the
blank-dotplot regeneration pass (see **Defect 4** below, `refit: false`) **overwrote both
manifests**, which now carry `fit_seconds: null` and `determinism_probe: null`, and whose
`wall_seconds` / `peak_rss_gb` describe the *replot* — **73.6 s / 3.24 GB** (primary,
`timestamp 2026-08-04T22:36:20`) and **81.5 s / 4.24 GB** (sensitivity, `22:37:54`). So the fit's
70.5 min, the 76.0 min wall, the 4.1 GB peak, the stopping epochs 632 / 199 **and** the determinism
probe now all live only in `logs/mofaflex_{primary,sensitivity}.log`. The manifest should carry the
stopping epoch, and a replot must not clobber fit provenance.

**RESOLVED 2026-08-04.** The clobbered keys were restored into both manifests from `logs/mofaflex_{primary,sensitivity}.log`, the only surviving source: `fit_seconds` 4,230 s / 2,022 s (70.5 / 33.7 min), `wall_seconds` 76.0 / 35.3 min, `peak_rss_gb` 4.1 / 8.3, `n_epochs` **632** / **199** against a 1000 cap (the log says *Training converged after N epochs*, so this is convergence, not a cap hit), and the primary's `determinism_probe` `{'epochs': 20, 'max_abs_weight_diff': 0.0, 'bitwise_identical': True}`. `run_mofaflex.py` now carries `fit_seconds` / `n_epochs` / `determinism_probe` / `peak_rss_gb` / `wall_seconds` forward when it reuses a cached model, and records the replot separately under `last_replot`, so this cannot recur. Each manifest also carries a `provenance_note` saying the values were restored rather than remeasured.

## Why the NMF comparison matters

NMF decomposes whatever is in `lrdata.X`. Run on `bivariate` output, the resulting programs are
defined over the 131 features that survived `nz_prop` — i.e. **partly an artifact of an
upstream filter calibrated for a different data modality**. Run on `inflow`, they are defined
over 4,608 features chosen by a threshold intended for this modality. Comparing the two measures
**how sensitive the recovered program structure is to the upstream filter**, which is itself a
result worth reporting rather than a nuisance.

This is also the **only comparator output that shares ALARMIST's shape** — locations × factors
and features × factors — so it is the only place a structural comparison to ALARMIST's BPTF
motifs is even possible.

## Why LIANA cannot recover the mGAM loop — it is the UNIT OF ANALYSIS, not the feature indexing

**Correction of emphasis, 2026-08-06.** This file's *Consequences, measured* table frames the
LIANA/ALARMIST difference as one of **feature identity** — LR pair vs `sender × LR pair`. That
difference is real, but it is **secondary**. The dominant reason LIANA cannot see ALARMIST motif 1
is that **inflow scores the receiving CELL**, so the two arms of a bidirectional loop land on
**disjoint populations of rows**. Measured by `why_no_mgam_motif.py` →
`vs_alarmist/why_no_mgam_motif.json`:

| | non-zero cells | share | dominant receiver |
|---|---|---|---|
| `mGAM^GRN^SORT1` | 3,183 | 3.18% | MES-like 19.2%, AC-like 17.0%, NPC-like 16.7% |
| `MES-like^ANXA1^FPR1` | 944 | 0.94% | **mGAM 45.2%** |
| **both** | **95** | **0.095%** | — |

Pearson r between the two arms across cells is **+0.0177**, Spearman **+0.0385**. The two arms are,
at cell resolution, essentially independent.

**The same two interactions, the same tissue, the same DB — only the unit changes:**

| | cells (LIANA inflow) | 50 µm patches (ALARMIST) |
|---|---|---|
| rows | 100,190 | 13,113 |
| rows carrying BOTH arms | 0.095% | **1.235%** |
| Pearson r | +0.0177 | **+0.4562** |
| Spearman ρ | +0.0385 | **+0.4044** |
| P(arm 2 \| arm 1) ÷ marginal | 3.2× | **14.2×** |

**Pearson rises 26× purely by aggregating cells into patches.** The loop is a **neighbourhood**
property, and a neighbourhood is not a row in LIANA's matrix. **No factorisation can recover
structure the input does not contain** — which is why neither the sensitivity arm nor the
reachability-normalised arm (which *does* restore the mGAM view) puts the two arms on one factor.

### Where the two motif-1 arms actually rank — reachability-normalised fit, 779 features

Recomputed from `mofaflex_inflow_reachnorm/data/mofaflex_loadings.csv`:

| feature | peak factor | loading | rank | |
|---|---|---|---|---|
| `mGAM^GRN^SORT1` | Factor 19 | **+0.329** | **67 / 779** | top 8.6% |
| `MES-like^ANXA1^FPR1` | Factor 1 | **+0.078** | **276 / 779** | top 35.4% — **two-thirds of all features load more strongly on Factor 1 than it does** |

`mGAM^GRN^SORT1` is **spread flat with mixed signs** and no factor claims it: F19 +0.329 (67),
F7 −0.327 (80), F9 −0.300 (62), F16 −0.286 (119), F1 +0.279 (99). **Neither arm is a top-10 feature
of any factor**, so **neither ever appears in `top_weights.png`**.

The reverse directions outrank the biologically meaningful one: `MES-like^GRN^SORT1` peaks on
Factor 4 (+0.258, rank 121) and `mGAM^ANXA1^FPR1` on Factor 1 (+0.165, rank 164). The
autocrine-ish myeloid direction `mGAM^ANXA1^FPR1` carries a larger |weight| than
`MES-like^ANXA1^FPR1` on **17 of the 20** factors — including **every one of the seven** on which
either of them ranks in its own top five. The three exceptions (Factors 12, 13, 18) are factors
where neither is prominent.

The best joint factor by worst-of-the-two rank is **Factor 1** (ranks 99 and 276), same sign. The
two arms in fact share a sign on **15 of the 20** factors, so **the model is not placing them at
opposite poles — it simply places neither anywhere prominent.** *(A shorter list — "Factors 1, 6, 7,
15" — circulated during this pass; the full recomputation gives 15 of 20. The conclusion is
unchanged.)*

## Factor-vs-motif cosine — defects in the COMPARISON itself, recorded 2026-08-06

`cosine_factors_vs_motifs.py` follows the matching procedure in `.claude/skills/alarmist`. Three
confounds have to be handled or the heatmap is an artefact. Key canonicalisation was handled
correctly from the start; the other two were not. **The sign handling was outright WRONG (C2), and
the cost of the feature-space collapse was never stated (C1)** — and correcting the first changes
the headline reading.

### C1 — the (sender, L, R) collapse is biased in LIANA's favour, and is an UPPER BOUND

The two feature spaces are not the same object: ALARMIST is
`(sender, receiver, ligand, receptor, contact mode)` — **25,271** features per motif — and MOFA-Flex
is `(sender, ligand, receptor)` — **779**. The only common space is the latter, so **ALARMIST must
be collapsed by SUMMING over receiver and contact mode**: 25,271 rows → **4,756 keys**, a median of
**5** rows merged per key (mean 5.31, max 10; median 5 receivers).

**The worked example is the number to quote.** For `mGAM|GRN|SORT1` in motif 1, the collapsed
`score` is **12.369**, of which the biologically meaningful MES-like arm contributes **3.091 —
25 %**. mGAM→mGAM autocrine (1.144 + 0.758) is summed in as if equivalent. So the coordinates
ALARMIST has and LIANA does not are **discarded before the cosine is taken**, and they are
discarded in the direction that helps LIANA. **This comparison is an upper bound on the agreement,
not a neutral measurement, and must never be quoted as though it were symmetric.**

ALARMIST mass retained after collapse + join: **90.9 %** on raw `V`, **74.0 %** on `V*`.

**`aggfunc` is not load-bearing** (recomputed, not on disk): switching sum → max → mean moves
max|cos| from 0.643 → 0.678 → 0.660 in `signed` mode and 0.743 → 0.783 → 0.760 in `poles` mode, with
**1 / 20** motifs above 0.5 by |cos| in all three signed variants and 13 / 15 / 13 in poles. The
qualitative reading survives; the poles count moves by 2, so quote the `sum` figures.

*(Key canonicalisation — the second confound — is the DB-export trap already recorded above: sorting
subunits before joining lifts the raw key overlap from **713** to **742**.)*

### C2 — 🔴 `abs()` was used to reconcile signed weights with a non-negative factorisation. That was WRONG.

A BPTF motif is non-negative and means *"these interactions are high together"*. A MOFA-Flex factor
is an **axis**: its two poles are **anti-correlated**, and its global sign is arbitrary. Taking
`|weight|` **merges the two poles**, scoring a factor that contrasts A against B as if it contained
A and B together.

**On this run that is a large distortion, not a rounding:** **57.2 %** of the 779 × 20 weights are
negative, and the **minor pole holds a median 38.5 % of each factor's mass** (range 13.9–49.3 %).

Three modes are now implemented, `--sign-mode {poles,abs,signed}`:

- **`poles` (DEFAULT)** — split each factor into `Factor_k (+) = max(w,0)` and
  `Factor_k (−) = max(−w,0)`, giving **40** genuinely non-negative vectors that preserve the
  contrast.
- **`abs`** — `|weight|`. **Wrong, and what was originally used.** Retained only for provenance.
- **`signed`** — signed weights as-is; cosine may be negative, plotted on a diverging `RdBu_r` scale
  centred at 0.

**Results** (ALARMIST scored on `V* = V/(mean_LR+1)`, 742 shared keys, 200-permutation null that
permutes each MOFA column over the shared keys, preserving its sparsity and magnitude):

| sign mode | max cos | median best-per-motif | motifs > 0.5 | motifs beating the null at p<0.05 | null median max-cos |
|---|---|---|---|---|---|
| `abs` | 0.671 | 0.432 | **3 / 20** | 18 / 20 | 0.241 |
| **`poles`** | **0.743** | **0.517** | **13 / 20** | **20 / 20** | 0.221 |
| `signed` | 0.478 (\|cos\| max 0.643) | 0.345 | **0 / 20** | 20 / 20 | 0.088 |

⚠️ **The `abs` null figures are a RECOMPUTATION, not a file on disk.**
`vs_alarmist/cosine_mofaflex_reachnorm_vs_alarmist.json` predates the permutation block and carries
**no** `null_*` keys; the 0.241 / 18-of-20 above were reproduced here with the script's own
procedure and seed (the same recomputation returns the stored `poles` values 0.2211 / 20-of-20
exactly, which is the check that it is faithful). The signed cosine range across both value columns
is **−0.679 … +0.633**; a quoted range of "−0.55 … +0.63" does **not** reproduce.

On **raw `V`** instead of `V*`, `abs` gives 19 / 20 motifs above 0.5 — that is a **shared-prevalence
artefact** (both sides are dominated by the same high-prevalence adhesion pairs), not agreement.
Report `V*`.

**The poles/abs difference is not multiple comparisons.** Going from 20 to 40 candidate vectors
**LOWERED** the null (median max-cosine 0.241 → 0.221), because pole vectors are sparser. The
improvement is therefore not bought by having more chances.

**✏️ Correction to the record.** An earlier statement from this pass — *"3 / 20 motifs above 0.5, so
the methods are not recovering the same programs"* — was an **artefact of the `abs()` handling**.
The corrected reading is more specific and less flattering to both sides:

- **they agree substantially on LR VOCABULARY** — 13 / 20 motifs match a pole above 0.5, and 20 / 20
  beat the permutation null;
- **they do not agree on how that vocabulary groups into programs** — in cell space the ceiling is
  |ρ| ≈ 0.46 and **9 of 20 motifs collapse onto a single hub factor** (Factor 18, which correlates
  +0.422 with total ALARMIST loading and +0.362 with total inflow per cell — it is a general-activity
  axis, not a program). ALARMIST motif 1 correlates only **+0.146** with total inflow, i.e. it is
  specific, and its best MOFA-Flex match is Factor 18 at **ρ = +0.214**;
- **they do not agree on receiver or direction at all**, because LIANA does not represent either.

*(Cell-space numbers from `compare_programs_to_alarmist.py` → `vs_alarmist/comparison_summary.json`.
⚠️ **The NMF row of that file is UNVERIFIED and must not be quoted:** the script reads 8 numeric
columns from `NMF_W_factor_scores.csv`, which was written with `index=False`, while the `nmf_inflow`
rank is **7** — an off-by-one index bug. The three MOFA-Flex rows are unaffected.)*

Figures: `vs_alarmist/figures/cosine_*` (sorted by best match) and `clustermap_*` (clustered on both
axes), `Reds` for the two non-negative modes and `RdBu_r` for `signed`, PNG + PDF + SVG, each
suffixed by sign mode.

## `nz_prop` — the four different values in play

| source | value | note |
|---|---|---|
| `_bivariate.py` signature default | **0.05** | the package default for the spot-based branch |
| `_inflow.py` signature default | **0.001** | 50x lower; the authors' answer to single-cell sparsity |
| `bivariate.ipynb` (tutorial call) | **0.2** | overrides the default; 20% of Visium SPOTS |
| LIANA+ paper text | **10%** | a third value again, stated in the manuscript |

Four different numbers for the same parameter across code, tutorial and paper. We use **0.02**
for bivariate (binomial spot->cell conversion of the tutorial's 0.2: 20% of ~15-cell spots is
~1.5% of cells) and **0.001** for inflow (its own default, untouched). Any statement about how
many interactions a LIANA run "finds" is meaningless without saying which of these was used.

## ⚠️ The inflow branch was truncated — corrected 2026-08-04

`run_inflow.py` ended after `li.mt.inflow`. It created `plots/` and wrote nothing into it (no
plotting code existed), and it never ran the five downstream steps of `inflow_score.ipynb`:
`compute_global_specificity` (cell-type, region, and `cell_type::region` variants),
`spatial_pair_proximity`, and the spatially-constrained `rank_aggregate`. **Nothing errored —
the calls were simply absent**, which is why the gap survived until the empty `plots/` directory
was noticed. Without `compute_global_specificity` the branch had no source→target significance
at all. `run_inflow_downstream.py` now performs all five and produces 80 figures (0 blank).

Deviations introduced by that script:

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| region variable | `major_brain_region` | **`grade`** (high/low) | The GBM TMA has no anatomical-region annotation. `grade` is the closest tissue-class analogue. `obs` also carries `motif` and `patch_id`, which are **ALARMIST outputs** — the script hard-refuses them as `--region-col` and never reads them. |
| proximity bandwidth | 100 µm | 100 µm (unchanged) | The tutorial deliberately uses a coarser scale for proximity than for the inflow kernel; kept as-is rather than harmonised to 13.1454. **Consequence:** at 100 µm the `interacting` column of `pair_proximity.csv` is **1 for all 81 source→target pairs** — it is saturated and carries no information. Use the continuous **`proximity`** column (range 0.0435–0.9933) instead; that is the usable output. |
| interactions plotted | one hand-picked (`Ntf3^Ntrk3`) | top-6 by `lr_mean` **+ GRN^SORT1, ANXA1^FPR1** | Standing project rule: always plot the two ALARMIST motif-1 LRs alongside the method's own top hits. |
| score-map colour scale | default | `vmin=0, vmax='p99.5'` | 99.45% zeros + long tail ⇒ default scaling gives an all-black map. Percentile scaling is the tutorial's own idiom (`percentile_scaling=(1,97)`). |

## `plot_liana_full.py` and the `plots_full/` trees — three defects, all FIXED 2026-08-04

*(A fourth, **Defect 4**, is filed at the end of this section. It is not a `plot_liana_full.py`
defect — it is `run_mofaflex.py` — but it is the same trap as Defect 1 and it is what showed
Defect 1's blank-figure audit had not been exhaustive.)*

Neither this file nor `NOTES.md` mentioned `plot_liana_full.py` or `plots_full/` before
2026-08-04, which is why the three defects below survived. Recorded now.

### What `plots_full/` actually is — and what it is NOT

`plots_full/` means **"the output of `plot_liana_full.py`"**, a second additive plotting pass. It
does **not** mean "the full set of plots", and it is **not** a superset of `plots/`.

**The two directories are strictly DISJOINT — not one filename appears in both:**

| directory | written by | contents |
|---|---|---|
| `nmf_*/plots/` | `run_nmf.py` | exactly 2 files: `elbow.png` (rank selection) and `factor_maps.png` |
| `nmf_*/plots_full/` | `plot_liana_full.py` | **30** real files (bivariate: 25 PNG + 4 CSV + `plot_manifest.json`) / **62** (inflow: 55 PNG + 6 CSV + manifest) under `global/`, `factors/`, `interactions/` (incl. `interactions/requested/`). *(22 / 31 before the 2026-08-04 regeneration — the increase is Defect 3's `requested/` dirs and the `feature_by_group` loop.)* ✅ **The manifests now match disk exactly: `n_files` 30 / 62 plus a separate `n_png` 25 / 55.** *(They previously said 31 / 63, one too many in each tree, because the closing glob counted `.DS_Store`. Fixed at `plot_liana_full.py:479-480`, which skips dotfiles and counts PNGs separately; re-verified with `find`.)* |

**Do not delete `plots/` as redundant** — it holds the only rank-selection evidence
(`elbow.png`) and the main factor overview, neither of which has any counterpart under
`plots_full/`. The name is also overloaded across comparators: `../stlearn/DEVIATIONS.md:26` uses
`plots_full` to mean "went beyond the tutorial". A less confusing name (`plots_extra/`) is worth
adopting at the next rebuild, but renaming now would break every path already written down.

### Defect 1 — two blank `connectivity.png` (FIXED)

`nmf_bivariate/plots_full/global/connectivity.png` and `nmf_inflow/plots_full/global/connectivity.png`
were both written as **pure white** — 6,365 bytes, 900 × 1188 px, exactly **1 unique colour**,
mean pixel value 255.0. A pixel test over the ~140 PNGs that existed in the LIANA tree **at the
time** flagged only these two — but that audit was **not** exhaustive over the final tree. A later
sweep over all **420** PNGs found **two more**, from the identical trap in a *third* script; see
**Defect 4** below. Four blanks total, across three scripts; the tree is now **420 PNGs, 0 blank**.

**Cause:** `li.pl.connectivity` draws with **plotnine**, not matplotlib. `plot_liana_full.py`
called it without `return_fig=True` and then saved `plt.gcf()` — an empty matplotlib canvas. No
exception was raised, so the script's `guard()` logged "ok" and the blank was written.

The correct idiom already existed 200 lines away in `run_inflow_downstream.py:167-173`
(`return_fig=True` + `fig.save(...)`, with an inline comment naming the trap); its outputs
`cellchatdb2_inflow/plots/global/connectivity_idx{25047,50095}.png` are real figures.

**Fixed:** `plot_liana_full.py` now takes `return_fig=True` and saves the ggplot through a
plotnine-aware saver, and it carries a `_is_blank()` detector — the same one
`run_inflow_downstream.py::_save` uses — so a figure with no drawn content is *refused and logged*
rather than written. Both `plots_full/` trees **have been regenerated (2026-08-04 18:04)** and both
`connectivity.png` are now real figures — 1080 × 1080, ~2,485 unique colours, mean pixel value
217.1 — re-verified with a PIL pixel test. Impact was cosmetic: connectivity is illustrative and
nothing downstream reads it.

### Defect 2 — row-order misalignment corrupted the SENDER and LR-pair axes (FIXED)

**Cause:** `li.ut.get_variable_loadings` re-sorts its rows by `|Factor1|` descending
(`liana/utils/_getters.py:149`), so `data/NMF_H_loadings.csv` is **not** in
`nmf_WH.npz['features']` order. `plot_liana_full.py` built its *annotation* vectors
(`sender_of`, `lr_of`) from the npz but its *value* vectors from the CSV, then masked one by the
other — silently pairing each feature's **name** with a **different feature's loading**. The two
orderings agreed at **0 of 2,704** positions on inflow and **1 of 131** on bivariate, so the
mislabelling was near-total, not marginal.

It was a pure permutation, not a data error — the underlying fit was never wrong.

**Outputs it corrupted** (all under `plots_full/factors/`):

| branch | corrupted |
|---|---|
| `nmf_inflow` | `factor_by_SENDER_celltype.csv`, `factor_by_LRpair.csv`, `heatmap_factor_by_SENDER.png`, `heatmap_factor_by_LRpair_top30.png`, `celltype_communication_by_factor.png`, `identity_Factor1..7.png` |
| `nmf_bivariate` | `factor_by_LRpair.csv`, `heatmap_factor_by_LRpair_top30.png` |

The corruption was visible *inside* a single figure: `nmf_inflow/.../identity_Factor1.png` showed
all ten top-10 bars as `Glial-Neuronal^…` under a suptitle reading "top sender: AC-like".

**Not affected** (verified): everything `run_nmf.py` itself wrote — `NMF_H_loadings.csv`,
`top10_loadings_per_factor.csv`, `nmf_WH.npz`, `plots/elbow.png`, `plots/factor_maps.png`; and
within `plot_liana_full.py` the RECEIVER matrices (built from `W`, npz-ordered), the top-10 bar
panels, and `top_lri_dot_by_factor.png` (label-indexed).

⚠️ **`../METHODS.md`'s prose sender list (in *How much does the program structure change?*) is RIGHT** — it was derived from
`top10_loadings_per_factor.csv`, not from the broken figures. Do **not** "correct" it to match the
old figures; the figures were what was wrong.

**Fixed:** `plot_liana_full.py` now reindexes once, immediately after load
(`loadings = loadings.loc[list(feats)]`), asserts the alignment, and raises if any npz feature is
absent from the CSV — so the two files can never again be assumed to be co-ordered. Both trees
**have been regenerated (2026-08-04 18:03–18:04)**, so the corrupted CSVs and figures listed above
no longer exist on disk; what is described here is the pre-fix state, kept for provenance.

### Defect 3 — the required-LR rule was unmet in both NMF branches (FIXED)

`SKILL.md:56-61` requires that **wherever a plot is per-LR**, two sets are produced: the method's
own top hits, **and always GRN→SORT1 and ANXA1→FPR1**, "in separate output dirs so the requested
pair is never mistaken for the method's own ranking".

`plot_liana_full.py` had **no required-LR mechanism at all** — selection was purely
`top_int = strength.index[:n_top_interactions]`. So:

| directory | status before the fix |
|---|---|
| `cellchatdb2/plots/requested/` | ✅ compliant — `rank33_GRN-SORT1{,_pvals}.png`, `rank60_ANXA1-FPR1{,_pvals}.png`, in their own dir |
| `cellchatdb2_inflow/plots/interactions/` | ⚠️ present but **not segregated** — the 16 required-LR files sit in the same 64-file dir as the top-6 |
| `nmf_bivariate/plots_full/interactions/` | ❌ **absent entirely** — 6 files, all top hits |
| `nmf_inflow/plots_full/interactions/` | ❌ **absent entirely** — 11 files, all top hits |

This mattered most where it was missing: the NMF branch is the **only** comparator output that
shares ALARMIST's shape, and it had **no factor-level per-LR view of either arm of motif 1**. The
data was on disk all along — both plain pairs survive into `nmf_bivariate/data/nmf_WH.npz`, and
**18** sender-resolved variants survive the ≥5/13 punch filter into `nmf_inflow`.

**Fixed:** `plot_liana_full.py` now takes `--required-lrs` (default `GRN^SORT1,ANXA1^FPR1`),
resolves each to every matching feature, plots them into a **separate** `requested/` directory,
and — per the contract — logs an explicit "ABSENT from the NMF feature set … that is itself a
result" when a required LR does not survive upstream filtering.

**Verified on disk after regeneration (2026-08-04):** `nmf_bivariate/plots_full/interactions/requested/`
holds 5 files (`violin_` + `feature_by_group_` for each of the two LRs, plus
`requested_lr_ranks.csv`) and `nmf_inflow/.../requested/` holds 25 (6 plotted sender-resolved
features × `violin_`/`feature_by_group_`/`panel_`/`sender_receiver_`, plus the CSV). Both LRs are
**present, not absent**: all **18** sender-resolved variants (9 senders × 2 LRs) survive the ≥5/13
punch filter into `nmf_inflow`'s 2,704-feature set, and both plain pairs are in
`nmf_bivariate`'s 131. Ranks recorded by `requested_lr_ranks.csv` — bivariate `GRN^SORT1` **25 /
131**, `ANXA1^FPR1` **97 / 131**; inflow best `MES-like^GRN^SORT1` **153 / 4,608**,
`MES-like^ANXA1^FPR1` **352 / 4,608**.

### Defect 4 — two blank `dotplot_focus_factors.png` in the MOFA-Flex tree (FIXED 2026-08-04)

**Not a `plot_liana_full.py` defect** — recorded here because it is the *same* trap as Defect 1 and
because it shows Defect 1's audit was not exhaustive.

`mofaflex_inflow/plots/dotplot_focus_factors.png` and the same file under
`sensitivity_nzf0.001/plots/` were **4,377-byte all-white PNGs**. Cause identical to Defect 1:
`li.pl.dotplot` returns a **plotnine** ggplot and `run_mofaflex.py` saved `plt.gcf()`. That makes it
the **third** appearance of the plotnine return-type trap in this codebase — after
`plot_liana_full.py`'s `connectivity.png`, and after `run_inflow_downstream.py`, which had handled
it correctly from the start.

**Fixed:** `run_mofaflex.py` now routes plotnine returns through `save_gg()` (`:106`, `:139`), and its
`save_current()` (`:118-131`) **refuses to write a blank canvas** rather than logging "ok", naming
the trap in the failure message. Both figures were regenerated **from the cached model — no refit**
(`logs/comparators/liana-mofaflex-replot{,-sens}.log`; the manifests record `refit: false`), and are
now 154,414 and 174,660 bytes. ⚠️ That regeneration is also what clobbered the fit provenance in both
MOFA-Flex manifests — see the *provenance gap* note in the MOFA-Flex section above.

**Final sweep over the whole LIANA tree: 420 PNGs, 0 blank** (`cellchatdb2` 35, `cellchatdb2_inflow`
80, `cellchatdb2_morans` 35, `nmf_bivariate` 27, `nmf_inflow` 64, `default` 35, `default_inflow` 80,
`nmf_bivariate_default` 2, `nmf_inflow_default` 2, `mofaflex_inflow` 30, `lric_percore` 12, `misty`
7, `factor_annotation` 10).

### Provenance gap — the original invocations were never recorded (since largely closed)

**The `plot_liana_full.py` command lines that produced the on-disk `plots_full/` trees are
recorded nowhere.** The literal string `plots_full` does not appear in the script; it came from
`--out-dir` on a command line that was not saved, and `plot_manifest.json` does not store argv.
The invocations had to be *reconstructed* from the manifests in order to regenerate.

Worse, **exactly one log exists for five scripts**: `logs/liana_inflow_downstream.log`. There is
no log for `run_liana.py`, `run_inflow.py`, `run_nmf.py` or `plot_liana_full.py`. Since
`plot_liana_full.py`'s `guard()` catches and logs exceptions rather than raising, **any figure it
silently skipped is unrecoverable** — we cannot now tell whether anything failed and was swallowed.
The **four** blank PNGs — two `connectivity.png` (Defect 1) and two `dotplot_focus_factors.png`
(Defect 4) — are the cases we caught, and only because a blank file is visible on disk; a
`guard()`ed failure leaves no file and no trace at all. *(This read "the two blank connectivity PNGs
are the one case we caught" until 2026-08-04, when the second pair was found.)*

Going forward: tee every comparator invocation to `logs/comparators/<method>-<tier>-<stage>.log`,
and have each script record its own `argv` into its manifest.

**Status after the 2026-08-04 regeneration — mostly closed, verified on disk:**

- `plot_manifest.json` now records enough to reconstruct the invocation exactly: `script`,
  `input`, `dataset`, `tier`, `inflow_dir`, `nmf_dir`, `h5ad`, `out_dir`, `required_lrs`,
  `global_specificity`, `bandwidth`, `cutoff`, `seed`, **`git_sha`**, `n_files`, `wall_min`, plus
  the two fix-confirmation flags `blank_figures_suppressed` and
  `loadings_reindexed_to_npz_order`. It still does **not** store a literal `argv` key.
- `logs/comparators/` now holds **18 LIANA logs** (plus `liana-envexport.err`) — it had five when
  this bullet was first written. In full: `liana-bivariate-plotsonly.log`,
  `liana-plotfull-{bivariate,inflow}.log` **×2 each** (one plain, one with the argv suffixed into
  the filename), `liana-inflow-downstream-v2.log`, `liana-punch-level.log`,
  `liana-GBM-cellchatdb2-morans.log`, the **seven** `liana-GBM-default-*.log`
  (`bivariate`, `inflow`, `downstream`, `nmf-biv`, `nmf-biv-k21`, `nmf-inf`, `nmf-inf-k21` — the two
  `-k21` logs being the `--k-max 21` refits), `liana-mofaflex-replot{,-sens}.log` (the Defect 4
  blank-dotplot regeneration) and `comp-liana-install-mofaflex.log` — alongside the original
  `logs/liana_inflow_downstream.log`. The `-k21` and `-replot` logs are the evidence for the two
  post-documentation changes recorded in this file. There is
  still **no log for the original `run_liana.py` / `run_inflow.py` / `run_nmf.py` fits**, so any
  `guard()`-swallowed failure from *those* remains unrecoverable; the paragraph above stands for
  the pre-2026-08-04 runs.

## Three MOFA-Flex figures were described wrongly — corrected 2026-08-06

**No figure changed and no number moved; what was wrong is the reading.** Each semantic below was
verified against the installed `mofaflex` source and against the CSVs the same call writes.

### `top_weights.png` — the x-axis is |weight|, and the sign is in the glyph with no legend

`mfl.pl.top_weights`, `mofaflex/pl/_plotting.py:1118-1171`:

- **x is `| Weight |`** (`:1165`), the **absolute** loading. The sign is carried **only by the
  marker shape** — `$\oplus$` for w ≥ 0, `$\ominus$` for w < 0 — and `scale_shape_manual(...,
  guide=None)` at `:1162` means **no legend is drawn for it**. A reader who does not know this reads
  a signed factorisation as if it were non-negative.
- Within each facet the top *n* are taken **by |weight|** (`:1141`, `weightabs.argsort()`) and then
  sorted ascending, so the **largest bar is at the TOP** of each panel.
- `facet_wrap("factor", scales="free")` (`:1169`): **every panel has its own x-scale.** Bar lengths
  are **not comparable across factors**.
- These are **raw weights** — no prevalence normalisation of any kind (LIANA has no analogue of
  ALARMIST's `V* = V/(mean_LR+1)`).
- Now regenerated at **top-10** per factor (`--top-weights-n`, `run_mofaflex.py:107`, default 10);
  the tutorial's cell 43 uses 5.

### `variance_explained.png` — it is a HEATMAP, not bars ⚠️

⚠️ **A description of this figure as "bars labelled by sender" circulated during this pass and is
wrong.** `mfl.pl.variance_explained`, `_plotting.py:477-527`, builds
`p9.ggplot(df_r2, p9.aes(x=x, y="factor", fill="R2")) + p9.geom_tile()` (`:513-515`):

- **rows = FACTORS, columns = views = SENDER cell types, fill = R²**, one facet (`group_1`, the
  single group). Rows are ordered by total R² descending, which puts the **largest at the BOTTOM**.
- `data/r2_per_view.csv` is the **column sum over factors** — verified to equal
  `r2_per_factor_view.csv` grouped by view exactly (max |difference| 0.0).
- Single darkest cell in the reach-norm fit: **Factor 5 × Glial-Neuronal = 0.1641**. The **Lymphoid
  column is blank** (largest cell 0.00021), which is the visual form of the R² = 0.0013 caveat above.
- `variance_explained_by_view.png` is the *same data* with `group_by="view"` (x = group, faceted by
  view). With one group that is **nine one-column facets** — the less useful of the two layouts.

### `circle_plot_Factor*.png` — the edges are not the factor ⚠️

The factorisation **never sees the receiver**. In `run_mofaflex.py:624-676`, only two things come
from the model: `source` (parsed out of the feature name) and *which* 10 edges are drawn
(`|loading|`). **`target` and the edge weight come from `inflow_means`** —
`lrdata.to_df().groupby(obs[cell_type]).mean()` — i.e. the *receiving* cell's own annotation and the
**raw mean inflow**, not anything factor-weighted, not even factor-signed.

**Verified:** `factor_interactions_Factor1.csv` and `factor_interactions_Factor6.csv` share **108**
`(source, ligand, receptor, target)` rows whose edge weights are identical to **0.000e+00**, while
their loadings on those same rows differ by up to **1.90**. Consequences:

- two factors that share a feature draw the **identical sub-network**;
- the edges ignore the factor's **sign**, so a `−` pole and a `+` pole look the same;
- the apparent "per-factor network" is mostly **not the factor**.

Read it as *"the top-10 interactions this factor selects, and where inflow of those interactions
generally goes"* — **not** as *"this factor's sender→receiver structure"*. Now emitted for **all 19
active factors** (`--circle-all-factors` / `--circle-top-n`, `run_mofaflex.py:109-113`); the tutorial
draws only its two hand-picked focus factors.

## Cross-punch reproducibility filter (applied before NMF)

The LIANA+ paper retains interactions present in **>=10 of 28 sections**. Scaled proportionally
to this dataset's **13 TMA punches** that is **>=5/13**. A feature counts as present in a punch if
it is non-zero in at least one cell of that punch. Applied to **both** branches immediately
before `li.multi.nmf`, so the factorisation runs on interactions that reproduce across punches
rather than ones driven by a single core. Per-feature punch counts are written to
`punch_presence.csv`.

## Smaller items on record — recorded 2026-08-04

None of these changes a result. All were recomputed from the files on disk; nothing was re-run.

- **`lr_ranking_by_lr_mean.csv` is misnamed — it is NOT a ranking by `lr_mean`.**
  `run_inflow_downstream.py:187-188` computes
  `gi[gi.pval < 0.05].groupby('lr')['lr_mean'].max()`: a **significance-conditioned maximum** over
  the 81 source→target pairs, not a mean and not an unconditional ranking. It has **616 rows
  against 633 unique LRs** — the 17 missing LRs are absent *by construction* (no source→target
  pair of theirs reached p<0.05), not because they scored low. The 6 "top" LRs that were plotted
  were selected by that conditional max. Read the file as "best significant sender→receiver score
  per LR"; and note the p<0.05 gate inherits the pseudoreplication caveat in **D2** above.

  **Since 2026-08-04 this is self-documenting**: the columns are now
  `lr, max_lr_mean_over_signif_pairs, n_signif_pairs` instead of a bare `lr_mean`, and an
  unconditional companion **`lr_ranking_all_pairs.csv`** (633 rows, `max_lr_mean_all_pairs`) is
  written alongside. Both rank `GRN^SORT1` 54th and `ANXA1^FPR1` 72nd — verified — so the
  conditioning never moved the two required LRs; the filename remains the only misleading part.

- **NMF rank selection is weakly supported — report rank 6 / rank 7 with these numbers, not bare.**
  `li.multi.nmf`'s elbow metric is **MAE** (`liana/multi/_nmf.py`), and an L1 elbow on a
  99.45%-zero matrix is a fragile criterion. Recomputed from the stored `W`, `H` and the matching
  input columns:

  | branch | rank | rel. Frobenius error | SS captured | achieved MAE | zero-predictor MAE |
  |---|---|---|---|---|---|
  | `nmf_bivariate` | 6 | 0.7607 | **42%** | 0.080583 | 0.072522 |
  | `nmf_inflow` | 7 | 0.7771 | **40%** | 0.015287 | 0.011817 |

  Both fits have a **worse MAE than predicting zero everywhere** — which is a property of L1 on an
  extremely sparse non-negative matrix, not a defect in the run (the Frobenius fit is respectable
  and the factors are spatially coherent). But "LIANA finds 6 / 7 programs" should never be quoted
  without the caveat that the selecting statistic does not beat a trivial baseline. Combined with
  **D1** — where bandwidth alone moves the inflow rank 11 → 7 — the rank is the *least* robust
  number in this comparator.

- **`cellchatdb2_inflow/data/inflow_lrdata.h5ad` still carries ALARMIST outputs.** Verified
  present: `obs['motif']`, `obs['patch_id']`, `uns['motif_colors']`, plus orphaned
  `obsp['connectivities']` / `obsp['distances']` (leftovers from a scanpy neighbours graph, not
  LIANA's — LIANA's is `obsp['spatial_connectivities']`). **Nothing in the shipped scripts reads
  them**: `run_inflow_downstream.py` hard-refuses them as `--region-col` and
  `downstream_manifest.json` records `region_col: grade`. So this is not a leak that has happened —
  it is a **trap for any future script that iterates `obs` keys generically** (a loop over
  `adata.obs.columns`, an automatic covariate sweep, a generic `groupby` scan) and would silently
  test a comparator against ALARMIST's own labels, i.e. circularly. Drop `motif` and `patch_id` at
  load in anything new.

## ✅ D3 (was: the `default` tier was never run) — CLOSED 2026-08-04, with one confound created and resolved the same day

`SKILL.md:51-54` requires the **`default`** tier (the method's own LR resource) **first** and
`cellchatdb2` second. Only `cellchatdb2` existed. `run_default_tier.sh` closes it: five steps,
`--db consensus`, everything else identical, one log each under
`logs/comparators/liana-GBM-default-*.log`.

**The two resources are not interchangeable:** LIANA's `consensus` has **4,624** unique
`(ligand, receptor)` pairs against CellChatDB v2's **3,218**, sharing only **1,663 — 36.0% of
consensus.** Both required LRs are in both.

**Result — and it settles a question this document had left open: the LR database is not a
confounder for any per-interaction statistic.** Joined on shared entries:

| comparison | shared | max &#124;difference&#124; |
|---|---|---|
| bivariate Moran's R (and `morans_pvals`, `mean`, `std`) | **79 pairs** | **0.000e+00** |
| inflow `lr_mean`, source×target×LR | **23,787 rows** | **0.0** |
| inflow `pval` | 23,787 rows | **0 rows differ** |

That follows from LIANA scoring each pair independently with no cross-pair normalisation. What the
resource *does* change: which pairs are tested (388 vs 131 bivariate; 9,448 vs 4,608 inflow
features), and the **ranking denominator** — `GRN^SORT1` moves from **33/131** to **77/388** with an
unchanged Moran's R of 0.035702. **So a LIANA rank is a statement about the resource; the score is
not.** Report ranks with the resource named.

### ✅ `k_range` confound — introduced by this run, then RESOLVED the same day

**What happened.** `run_default_tier.sh` called `run_nmf.py` **without `--k-max`**, and
`run_nmf.py:31` defaults it to **11**. So the first `nmf_bivariate_default` / `nmf_inflow_default`
carried `"k_range": [1, 11, 1]` while `nmf_bivariate` / `nmf_inflow` carry `[1, 21, 1]`.

This walked straight into a trap **this repo had already documented**: `../METHODS.md`'s
*"Choosing `k_range`"* section states that Kneedle normalises the elbow curve over the window it is
given, so **ranks from different windows are not comparable.** The default tier's first **4 / 4**
was fitted on 1..10 against `cellchatdb2`'s **6 / 7** on 1..20, so the reported "6/7 → 4/4 is a
database effect" was a **database × window** effect and could not be attributed to either.

**Resolution.** Both default-tier NMFs were refitted with `--k-max 21`, and
`run_default_tier.sh` now passes it explicitly with a comment saying why it must never be omitted.
On the matched window:

| branch | `cellchatdb2` | `default` | features |
|---|---|---|---|
| bivariate | rank **6** | rank **6** | 131 → 388 |
| inflow | rank **7** | rank **5** | 2,704 → 6,178 |

So the LR database's effect on the factor count is **much weaker than first reported**: **zero for
bivariate**, despite tripling the feature count, and **7 → 5 for inflow**. The headline "the
database changes the factor count" does not survive; what survives is that the *feature space* the
factorisation sees changes, and that per-interaction statistics are unaffected either way.

**Lesson worth keeping**: the confound was created by omitting one argument, and it produced a
plausible, quotable, wrong result that survived one round of documentation before being caught.
Any cross-tier rank comparison must assert the windows are equal before comparing.

Both default fits also fail the zero-predictor check, as the `cellchatdb2` fits do:
`nmf_bivariate_default` rel-Frobenius **0.758200** / **42.5%** SS (MAE **0.08525634** vs
zero-predictor 0.07779735), `nmf_inflow_default` **0.829425** / **31.2%** SS (**0.01661722** vs
0.01234411). *(Corrected 2026-08-04: this paragraph carried 0.770519 / 40.6% / 0.08690641 and
0.841216 / 29.2% / 0.01677688 — the **pre-refit** `k_range` 1..10 numbers — for one documentation
pass after the table above had already been corrected to rank 6 / rank 5, so the section
contradicted itself. Values above read off `nmf_*_default/run_manifest.json`, i.e. the `--k-max 21`
fits; only the two zero-predictor baselines are unchanged.)*

### Provenance blemishes in the default tier — none affects a number

- All four default-tier manifests, and `default_inflow/downstream_manifest.json`, record
  `"tier": "cellchatdb2"`: the string is hardcoded in `run_nmf.py` / `run_inflow_downstream.py`.
  The `resource` fields are correct (`"LIANA consensus"`, `resource_n_pairs: 4624`).
- `default/run_manifest.json` records `"dataset": "default"` (it is GBM) and
  `resource_fingerprint: null` (the fingerprint block only runs for file-backed resources).
- `default_inflow/downstream_manifest.json` records `"db": "<repo>/consensus"`, a path that does not
  exist — the literal `consensus` was joined to the repo root by the manifest writer. The loader is
  correct (`run_inflow_downstream.py:387-388` branches on `a.db == "consensus"`). Do not read that
  key as a file.

## 🔴 `n_factors = 20` is a CEILING, not a selection — new defect, recorded 2026-08-06

**Same failure class as the `k_range` confound above, and it must be read alongside it.** There the
rank came out of the *window* the elbow was searched in; here the factor count comes out of the
*argument*, and in both cases a plausible, quotable number was produced that is a property of the
call and not of the data.

**MOFA-Flex does not choose K.** `n_factors` is a hard ceiling: the model fits exactly that many,
and inactive factors are pruned *afterwards* by our own floor (`--r2-floor`, 2 % R² in at least one
view). The default of 20 was taken because `inflow_mofaflex.ipynb` cell 25 uses `n_factors=20`. It
**coincidentally equals ALARMIST's K = 20, and that is not the justification** — nothing about this
dataset was consulted in choosing it.

**The ceiling binds.** Verified from the three `run_manifest.json`:

| fit | active / requested | 3 weakest factors' share of total R² | weakest **kept** factor (max view R²) |
|---|---|---|---|
| tutorial QC (`mofaflex_inflow`) | **17 / 20** | 6.7% | — |
| reachability-normalised | **19 / 20** | 7.7% | Factor 15 at **0.0265**, vs the 0.02 floor |
| sensitivity `nzf>0.001` | **20 / 20** | 8.3% | — |

**There is no taper.** The three weakest factors still carry 6.7–8.3 % of total R² in every run, and
in the reach-norm fit the weakest survivor clears the floor by only 1.3× while the single pruned
factor sits at 0.0098 — i.e. the model is cut off mid-distribution, not converging on a natural
number of programs. **So 17 / 19 / 20 are artefacts of the argument and must not be compared to
ALARMIST's K = 20 as though both were fitted quantities.**

`run_mofaflex.py:505-519` now warns when ≥ 90 % of factors are active and records
`n_factors_requested`, `n_factors_active` and `ceiling_binding` in the manifest. ⚠️ **None of the
three manifests on disk carries those keys** — the guard post-dates every fit, so the table above
was recomputed from `r2.n_active_factors` and `active_factors.csv`, not read from a
`ceiling_binding` flag.

**Status: OPEN — a K-sweep (20 / 30 / 40 / 60) has NOT been run**, so where the active count
saturates is unknown.

⚖️ **State this symmetrically:** `results/GBM/analysis_parameters.csv` records **no selection
criterion for ALARMIST's K = 20 either**, so this may well be a shared weakness rather than a
LIANA-specific one. That is as far as the check went — nothing was measured about how ALARMIST's K
was chosen, and nothing more should be asserted.

## The `X > 0` cross-punch filter is only correct for a non-negative local metric

The *Cross-punch reproducibility filter* section above defines presence as "non-zero in at least one
cell", and `run_nmf.py` implements it as `X > 0`. The `cellchatdb2_morans` arm shows those are **not
the same predicate in general**: the `cosine` local matrix has **0 negatives and 83.6214% zeros**,
so `X > 0` ≡ `X != 0`; the `morans` local matrix has **4,378,791 negatives (33.4%) and 0.0000%
zeros**, so `X > 0` would silently treat a third of it as absent. It happened to be applied only to
`cosine`. Moran's R is also **not NMF-admissible** at all (negatives), which retroactively justifies
`cosine` as the NMF input — *post hoc*, not the original reason, and recorded as such.

## 🔴 OPEN — `scripts/comparators/` is not tracked in git

`git ls-files scripts/comparators` returns **nothing**, and `git check-ignore` confirms it is not
ignored either — it has simply never been added. Every run manifest records `git_sha: 95208de`,
which pins the **package**, not the comparator scripts that produced the results.

**This has cost verifiability three separate times during the 2026-08-04 pass**: an edit can be
confirmed to have the intended content, but there is no way to prove that no other line moved, and
no way to tie a manifest to the exact script revision that wrote it. Combined with the
*Provenance gap* section above (no `argv` in most manifests, logs missing for the original fits),
this is now the largest reproducibility hole in this comparator. **Fix: `git add scripts/comparators/`.**

## 🔴 Still open after the 2026-08-06 pass — consolidated

Nothing below is done. Each line points at the section that owns it; the section is authoritative.

| # | Open item | Owner section |
|---|---|---|
| 1 | **bandwidth 13.1454 µm remains a `SKILL.md:45-46` contract violation.** Kept by user decision; that is a decision, not a resolution. | *D1* + *D1, update 2026-08-06* |
| 2 | **No native multi-sample / differential mode.** LIANA has none; the punch-level Mann-Whitney is a documented substitute, marked ❌/N-A rather than ✅. | *D2* |
| 3 | **No K-sweep for MOFA-Flex.** `n_factors=20` is a binding ceiling and where the active count saturates is unknown. | *`n_factors = 20` is a CEILING* |
| 4 | **No prevalence normalisation of the MOFA-Flex loadings** — the analogue of ALARMIST's `V*` was not attempted, so the cosine on `V*` compares a normalised side against an unnormalised one. | *C1* |
| 5 | **Index bug in the NMF arm of `compare_programs_to_alarmist.py`** — reads 8 numeric columns where `nmf_inflow`'s rank is 7 (`NMF_W_factor_scores.csv` written with `index=False`). That row is unverified. | *C2* |
| 6 | ~~`cellchatdb2_inflow/plots/global/bandwidth_query.png` draws its vline on the wrong scale~~ — **NOT AN ISSUE, retracted 2026-08-07.** The guide is at the support radius 28.2096 µm, which is correct; see the retraction in *D1, update 2026-08-06*. The only defect is a clipped annotation. | *D1, update 2026-08-06* |
| 7 | **`scripts/comparators/` is still UNTRACKED IN GIT** — re-verified 2026-08-06, `git ls-files scripts/comparators` returns nothing. | *🔴 OPEN — not tracked in git* |
| 8 | **The HTML report is now STALE.** `reports/liana_plus_GBM_cellchatdb2/liana_plus_GBM_methods.html` is built from `_liana_report_sections.json` by `build_liana_report.py` and was **not** updated by this pass, so it still carries the pre-correction readings — including the `abs()` cosine result and the old figure semantics. Rebuild before circulating it. | this table |
| 9 | **No absolute cell-count floor accompanies `--nzf-mode reachability`.** Lymphoid enters the model at 2.5 % reachability and R² = 0.0013. | *The tutorial's global `nonzero_fraction` cut* |

## Housekeeping — 2026-08-04

- `env.lock.yml` re-frozen after installing torch / decoupler / mofaflex: **167 lines**, pinning
  `decoupler==2.2.0`, `torch==2.13.0`, `pyro-ppl==1.9.1`, `gpytorch==1.15.2`,
  `mofaflex==0.1.0.post2.dev179+g9792b435f`.
- `run_inflow_downstream.py --db` now accepts the literal **`consensus`** (previously a file path
  only) — that change is what let the default tier reach `rank_aggregate`.
- `plot_liana_full.py` no longer counts dotfiles in `n_files` and reports `n_png` separately
  (`:479-480`). Both `plots_full` manifests now match disk exactly.
- ⚠️ A stray **`.DS_Store` was written into `cellchatdb2_inflow/`** during the SVI measurement
  (19:40), and one sits at `results/comparators/liana/GBM/`. They are what inflated the old
  `n_files` counts — delete them rather than count around them.
- `decoupler` **is** installed (2.2.0) and is used by `annotate_factors.py`; `../METHODS.md`'s
  bivariate "plots we did not produce" table said otherwise and has been corrected.
