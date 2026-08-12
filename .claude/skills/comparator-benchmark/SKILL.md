---
name: comparator-benchmark
description: Run a cell-cell communication comparator (CytoSignal, stLearn, SpatialDM, COMMOT, LIANA+, NICHES, CellChat) on the GBM/LGG and LUAD datasets following the authors' default workflow, then document what the method does, inputs, and every data/image output in scripts/comparators/METHODS.md. Use when onboarding, running, re-running, or writing up any comparator method.
---
# Comparator benchmark

Goal for this phase: **understand what each method says about our data on its own terms.**
Not comparing embeddings to ALARMIST yet, not forcing a common "comparable unit".
Announce `comparator-benchmark` when you start using it (CLAUDE.md rule).

## Sources (already local — do NOT WebFetch)

- Upstream repos + tutorials:
  `/Users/jiayifan/tansey_lab/{stLearn,SpatialDM,COMMOT,liana-py,NICHES,CellChat}/`
- CytoSignal has **no local repo** — its tutorials are vendored HTML in
  `results/comparators/cytosignal/reference_notebook/` (4 files: main workflow, differential
  multi-dataset, custom LR DB, container conversion).
- CellChat = v2.2.0.9001, tutorials in `CellChat/tutorial/*.Rmd` (+ rendered `.html`). Use the
  **spatial** ones: `CellChat_analysis_of_spatial_transcriptomics_data`,
  `CellChat_analysis_of_multiple_spatial_transcriptomics_datasets` (multi-sample),
  `FAQ_on_applying_CellChat_to_spatial_transcriptomics_data`, `Update-CellChatDB`.
  Note its **default DB already is CellChatDB v2**, so its two tiers may collapse to one —
  confirm bundled-vs-`data/LRdatabase/` equivalence and say so rather than faking a second run.
- Dataset paths, cell-type columns, coordinate units: **CLAUDE.md**. Never infer them.

## Paths (repo root = /Users/jiayifan/tansey_lab/alarmist)

- Code:        `scripts/comparators/<method>/`        (tracked in git)
- Outputs:     `results/comparators/<method>/<dataset>/<tier>/`   (gitignored)
- Living doc:  `scripts/comparators/METHODS.md`
- Reference implementation: `scripts/comparators/cytosignal/` — match its structure.
- `<dataset>` is **`GBM`** or **`LUAD`** (existing convention — not `LGG`, not `AIS_LUAD`).
- `<tier>` is **`default`** or **`cellchatdb2`**.
- Plural `comparators`, singular `<method>` lowercase. Do not invent new spellings.

## Invariants

- **One conda env per method**, never base, never shared:
  `CONDA_SOLVER=classic /Users/jiayifan/anaconda3/bin/conda create -n comp-<method> -c conda-forge --override-channels ...`
  (libmamba is broken on this box — see CLAUDE.md). R methods additionally need a
  `scripts/comparators/<method>/activate_env.sh`, because `conda activate`/`conda run` fall
  through to the *system* R. Freeze with `conda env export --no-builds > env.lock.yml`.
- **Follow the authors' default workflow.** Same functions, same order, same argument values
  as the tutorial. No "improvements", no extra preprocessing, no reordering.
- **Keep each method's own neighborhood/kernel definition at its default.** Do NOT harmonize
  spatial scale across methods, and do not match it to ALARMIST's patch size.
- **Use the method's native multi-sample / differential mode** when it has one:
  GBM → split by `obs['grade']` (high vs low; 13 TMA cores in `obs['tma_id']` are the units).
  LUAD → AIS vs LUAD (the four `data/linghua/P{17,21}_{AIS,LUAD}_Xenium.h5ad` sections).
  If the method has no multi-sample mode, say so — do not hand-roll one.
- **Two tiers per method, in this order:** `default` (the tutorial's own DB) first, then
  `cellchatdb2` (`data/LRdatabase/CellChatDBv2.0.human.csv`, the LR resource ALARMIST uses).
  `default` proves we used the method as its authors recommend; `cellchatdb2` removes the
  DB as a confounder. Never mix the two trees.
- **Generate every plot the standard workflow can produce**, both tiers, and save all of them.
- **Whenever a plot is per-LR, produce two sets:** (a) the method's *own* top / most significant
  LRs by whatever statistic it ranks with, and (b) always **GRN→SORT1** and **ANXA1→FPR1**
  (the ALARMIST motif-1 mGAM loop), whatever their rank. Keep the two in separate output dirs
  so the requested pair is never mistaken for the method's own ranking. If a requested LR is
  absent from the method's DB or called non-significant, **say so** — that is a result.
  Other datasets' LRs of interest will be named by the user; do not invent them.
- **Persist quantitative outputs, not just plot images and interaction names.** Score
  matrices, p-values, and per-cell/per-spot assignments must land on disk in a re-readable
  form (see `quant_io.R` in the cytosignal reference).
- **Verify every call signature against the INSTALLED package**, not the tutorial, not memory.
- **Code and outputs are separate trees.** Never write generated artifacts into `scripts/`;
  never write `.py`/`.R` into `results/`. Run scripts take an output-dir argument.
- Write a `run_manifest.json` per run: dataset, tier, DB + row count, n_cells in/out of QC,
  every non-default parameter, seed, package version, wall-time, peak RSS, git SHA.

## Procedure — one method at a time, stop for sign-off before the next

1. **Read the tutorial from the local source.** List files first; never assume a path.
2. Write `scripts/comparators/<method>/NOTES.md`: one row per tutorial call —
   step | function + argument values | tutorial file/cell. This is the contract for step 4.
3. Build the env, install, **smoke test on the tutorial's own demo data** if it ships one.
   Freeze `env.lock.yml`.
4. Write `run_<method>.{py,R}` one-to-one against the NOTES.md table. Parameterize
   input dir / output dir / DB; hardcode nothing.
5. Smoke-run on a **crop** of the dataset first, then the full section.
   GBM (100k cells) before LUAD (up to 670k cells).
6. Run tier `default` → `results/comparators/<method>/<dataset>/default/`.
7. Run tier `cellchatdb2` → `.../cellchatdb2/`. Record every DB-conversion step as a deviation.
8. Update `scripts/comparators/METHODS.md` with this method's section (template below).
9. Report in chat: what ran, what the outputs are, every deviation + justification,
   and anything you had to guess. **Then wait for sign-off.**

Only after every method is signed off: build one self-contained HTML report **per dataset**
(`/interactive-report`) covering what was run and what each method found.

## Deviations

Any departure from the tutorial goes in NOTES.md **and** METHODS.md as a row
(what / tutorial did / we did / why) **and** is reported in chat. Never deviate silently.
If a tutorial call does not match the installed signature, that is a version gap: record it,
do not silently "fix" it. Parameter changes forced by our data (panel size, coordinate units,
QC thresholds) are deviations too — justify each with the number that forced it.

## STOP and ask — do not decide these yourself

- Tutorial assumes Visium/Slide-seq-specific structure and needs non-trivial rewriting.
- The method wants a cell-type annotation granularity ours doesn't match.
- Any LR-database choice, or more than one reasonable way to map CellChatDB v2 onto the
  method's DB format (multi-subunit complexes → simple LR pairs especially).
- Any patch/kernel/neighborhood scale parameter the tutorial does not pin down.
- Smoke test passes but the full run OOMs (LUAD at 670k cells is **known** to blow up —
  CytoSignal already hit this at ~57 GB).
- The method's output cannot be mapped onto the output contract in METHODS.md.
- Anything requiring a rerun/overwrite of an existing `results/comparators/` directory.

## METHODS.md section template

One `##` section per method, appended to `scripts/comparators/METHODS.md`:

```markdown
## <Method> — <language>, v<version>, env `comp-<method>`
**Core algorithm** — what it actually computes, in 3-6 sentences. Statistical model, what the
  null is, what a "significant interaction" means, what unit it is assigned to (cell / spot /
  cell-type pair / edge).
**Spatial model** — neighborhood definition, kernel, distance parameter + its default value
  and units. State explicitly if the method is non-spatial.
**LR database** — default resource + size; how CellChatDB v2 was converted for the
  `cellchatdb2` tier; how complexes are handled.
**Input** — exact objects/matrices/columns required, counts vs normalized, coordinate units.
**Workflow** — table: step | call with argument values | what it produces.
**Data outputs** — table: file | shape/schema | meaning. Every persisted file.
**Image outputs** — table: plot function | what it shows | file written. **Every plot the
  standard workflow can produce**, including ones we chose not to run (say why).
**Multi-sample / differential mode** — the native mechanism, or "none".
**Gotchas** — traps, silent failures, version gaps, memory behaviour, scale limits.
**Deviations from the tutorial** — table: item | tutorial | ours | why.
**Runs on our data** — table: dataset | tier | status | key numbers | output path.
**Methods paragraph** — 3-6 sentences in journal methods-section voice, naming the exact
  functions used, for pasting into the manuscript.
```
