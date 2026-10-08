# crossmethod — deviations

One row per departure from a documented contract, with the number that forced it.
`XM-n` ids are referenced from `METHODS.md` and `RESULTS.md`.

| id | item | contract says | we did | why |
|---|---|---|---|---|
| **XM-1** | Arial | CLAUDE.md: every publication figure in Arial | figures render in **DejaVu Sans** on iris | No Arial or Helvetica font file exists anywhere on iris (`/usr/share/fonts` and all six conda envs searched, 2026-09-10). `_common/plotting.py` sets the family to `["Arial","Helvetica","DejaVu Sans"]` and `svg.fonttype='none'`, so the **svg carries the Arial-first family string and renders Arial when opened on a machine that has it**; only the png/pdf rasterise DejaVu here. Re-export the pdf on a Mac before submission. |
| **XM-2** | SpatialDM feature set | all LR pairs the method tested | **579 of 1,692** pairs | SpatialDM computes local statistics only for pairs its global stage *selected* (FDR<0.1). The unselected 1,113 have no per-cell score to give a classifier. This is generous to SpatialDM — only its own best pairs compete. |
| **XM-3** | SpatialDM cell set | all labelled cells | 10 cells dropped (P17_AIS), 3 (P21_AIS) | Non-finite local z on a fixed handful of isolated cells, identical across every affected pair. The paired `alarmist@spatialdm` row is recomputed on the same reduced set, so both sides share the denominator. |
| **XM-4** | CytoSignal cell set | all labelled cells | **~13 % dropped** (20,074 of 153,730 for P17_AIS/clean) | `removeLowQuality(counts>=100, genes>=20)` inside `run_cytosignal.R`. Not re-run and not overridden — it is the method's own workflow. Paired ALARMIST row recomputed on the post-QC cells. |
| **XM-5** | CytoSignal slots scored | — | `diffusion` + `contact` (1,083 CCIs); `Raw` exported but **not** scored | `Raw_Raw_smooth` covers the same 171 interactions as `contact` without imputation; scoring both would double-count those interactions in the distribution. The export keeps it so the choice can be revisited. |
| **XM-6** | LIANA feature granularity | "a single LR pair" | `celltype^ligand^receptor` (9,313 features) | That is the unit LIANA's inflow branch actually writes; there is no per-LR inflow output to aggregate to without recomputing, which would breach "do not re-run". These features are *finer* than an LR pair and therefore a **stronger** baseline, not a weaker one. |
| **XM-7** | stLearn unit | cell | **~51 µm grid spot** | stLearn was run on its own grid (`n_row=209, n_col=181, spot_size≈51.4×51.2 µm`). Harmonising it away is forbidden (SKILL.md:45-46). Instead spots are labelled by majority vote of the pathologist labels of the cells in them, and the ALARMIST motif is averaged onto the *same* spots — the aggregation is applied to both sides or neither. Cells further than one spot diagonal from every centre are off-grid and dropped. Majority-vote ties (19 spots in P17_AIS/full) are counted negative. |
| **XM-8** | CellChat unit | per-cell score | piecewise-constant on cell type | CellChat has no per-cell score at all; its unit is (LR, sender type, receiver type) within a group. The derived feature takes ≤19 distinct values per section, so its AUC is bounded by cell-type identity. Built as incoming **and** outgoing (302 features) to be maximally generous. Report as a unit-of-analysis fact. |
| **XM-9** | `empirical_p` | — | reported, but **not** a p-value against a sampling null | Features are not exchangeable draws and the cells are spatially autocorrelated. The defensible statement is the rank: "k of N single LR pairs match or beat the motif". The column is kept because it is the compact form of that rank. |
| **XM-10** | reference features | not in the meeting takeaway | added (`kind == 'reference'`) | `local_density_80um` and the cell-type indicators exist to answer the obvious objection that a spatially smoothed feature tracks a large contiguous region for free. Clearly labelled; droppable. |

## Not deviations, but worth stating

- `scripts/comparators/` is still **not tracked in git** (METHODS.md open issue, unchanged by
  this work). These scripts are therefore unversioned like the rest of the tree.
- METHODS.md's *Open issues* still says "**LUAD is essentially untouched** … 0 of 7 methods".
  That is stale as of the iris port: five methods have complete four-section LUAD runs under
  `/data1/tanseyw/projects/fanj2/results/`. Not edited here — flagged for whoever owns that file.
- COMMOT and NICHES have **no LUAD run**, so they are absent from this benchmark. That is a
  coverage gap in the LUAD tier, not something this analysis chose.
