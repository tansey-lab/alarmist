# CytoSignal — deviations and findings

Created **2026-08-18** during the iris LUAD port. `METHODS.md` previously recorded that CytoSignal
had no `DEVIATIONS.md`; this is it. Tutorial-call deviations for the GBM run live in
`METHODS.md` § CytoSignal — what is here is what the LUAD port established.

---

## CS-1. `rankIntrSpatialVar` OOM-killed P17_LUAD at 96 GB — measured, not inferred

`P17_LUAD` (640,739 source cells, **536,388** post-QC) ran for 7 h 35 min under
`srun --cpus-per-task 8 --mem=96G` and was **SIGKILLed** inside `rankIntrSpatialVar` →
`SPARK::sparkx`, at Gaussian kernel 2 of 11 in the first of three `sparkx` calls.

**The measurement that settles it** — `sstat -j 8869477.0`:

| | |
|---|---|
| `MaxRSS` | 100,662,144 KiB = **95.9989 GiB** |
| `--mem=96G` cgroup cap | 100,663,296 KiB = 96.0000 GiB |
| headroom remaining | **1,152 KiB = 1.12 MiB** |

`sacct` alone cannot show this: only the R child died, so the job step is still `RUNNING` and has
no recorded exit state. Site config confirms the mechanism — `TaskPlugin=task/cgroup,task/affinity`,
`ProctrackType=proctrack/cgroup`, `cgroup.conf: ConstrainRAMSpace=yes`, and **no** `OverMemoryKill`,
so it is the kernel's cgroup OOM killer, which kills the largest task and leaves the step alive.

**Two independent corroborations, both from this repo rather than from logs:**

1. `bash` prints `Killed` only for SIGKILL. The stage-05 LIANA failure four days earlier printed
   `Segmentation fault (core dumped)` at `MaxRSS 47.85 GB` — different signal, different wording,
   different cause (an ARPACK int32 overflow). The wording is diagnostic.
2. `run_cytosignal.R` wraps the SPARK call in `tryCatch(..., error = ...)`. An R
   `cannot allocate vector of size N Gb` is a *condition*, so it would have been **caught**, and
   execution would have reached `save_lrscore_quant` and written `quant/`. The output directory is
   **empty**. Therefore it was not an R allocation error — it was a signal. `tryCatch` cannot catch
   SIGKILL.

**96 GiB is a FLOOR, not the peak.** The process was capped and killed there; the true requirement
is unknown and unbounded above by any measurement. It also died only *two kernels into the first
of three* `sparkx` calls, so most of the step's work had not yet been attempted. Any re-run must be
sized generously — the partition allows ~975 GB/node and another failure costs another ~7.5 h.

**Not the lever:** `numCores`. `run_cytosignal.R:57` hardcodes `numCores = 4` (which is why the log
says "Running with 4 cores" inside an 8-CPU allocation), but the `mclapply` children in `sparkx.sk`
close over length-912 vectors; fork is copy-on-write and what they dirty is R's small-object node
heap, not the multi-GB matrices. 4 → 1 buys ~1 GiB of 96.

**Also not the lever:** SPARK-X densifying. It does not. On the sparse, no-covariate path CytoSignal
uses, the `as.matrix()` branch is guarded by `if (ncol(counts) < 30000)` and n = 536,387 ≫ 30,000;
the 11 kernels are held one at a time. SPARK-X's own marginal cost is single-digit-to-mid-teens GiB.
**It was the last straw on an object that was already ~86 GiB resident**, not the cause.

---

## CS-2. `result.spx` is a RE-ORDERING of `result.hq`, not a filter — `METHODS.md` said otherwise

`METHODS.md:270` described `result.spx` as "+ spatially variable by SPARK-X (the headline tier)",
and the methods paragraph at `:438-439` said interactions "were reported at the `result.spx` level,
that is, significant, quality-controlled **and spatially variable by SPARK-X**". That is a sentence
written for the manuscript, and it is wrong.

**Measured on the completed `P17_AIS` outputs** (`quant/reslist_*.rds`, all three slots):

| slot | `result.hq` | `result.spx` | same interactions | same order | per-interaction cell sets identical |
|---|---:|---:|---|---|---|
| `contact_Raw_smooth` | 146 | 146 | yes | **no** | **yes** |
| `diffusion_Raw_smooth` | 642 | 642 | yes | **no** | **yes** |
| `Raw_Raw_smooth` | 68 | 68 | yes | **no** | **yes** |

`rankIntrSpatialVar` sorts `result.hq` by SPARK-X adjusted p-value and **discards the p-values** —
only a console message reports how many were significant. Nothing is removed.

Three consequences:

1. **The methods sentence must be reworded** to "reported at the `result.hq` level (significant and
   quality-controlled), ordered by SPARK-X spatial variability". Corrected in place at
   `METHODS.md:441`.
2. **`result = hq = spx` at `METHODS.md:401-402` is not evidence of anything.** `hq = spx` holds for
   *every* interaction in every slot, so it does not show that GRN→SORT1 or ANXA1→FPR1 cleared an
   additional bar.
3. **Skipping the step costs far less than it appears to.** Lost: the SPARK-X ordering, the `n_spx`
   column (it becomes 0), the `result.spx` list element, and `signif_plots/` — though the plot block
   is gated at `run_cytosignal.R` to ≤ 200,000 cells and all three remaining sections are far above
   it, so those figures would not have been produced either way. **Not** lost: every score matrix,
   every significance call and its cell sets, the `n_hq` ordering the summary CSVs are actually
   sorted by, and the entire NEBULA differential (`run_nebula_stage.R` rebuilds from `input/` and
   never reads these outputs). It is still a departure from the authors' documented step order and
   must be declared if used — the opt-out is `CS_SKIP_SPX=1`, default off.

---

## CS-3. The runner banked its results *after* the step that OOMs — fixed

`run_cytosignal.R` ran `inferSignif` (6 h 42 min on P17_LUAD), then `rankIntrSpatialVar`, and only
then `save_lrscore_quant`. **The finished 6 h 42 min result was four statements from disk when the
process was killed.**

The checkpoint flag is not the fix: `save_rds == "save"` writes its checkpoint *before*
`inferSignif`, so even in `save` mode it would have protected only the preceding ~53 min and left
all 6 h 42 min exposed — for ~45-50 GB of gzip.

**Fixed by reordering, additively:** `intr_names` is hoisted above the SPARK block and
`save_lrscore_quant` is now called **twice** — once before `rankIntrSpatialVar` and once after. The
pre-SPARK write lands everything except `n_spx`; `save_lrscore_quant` tolerates an absent
`result.spx` (`quant_io.R` → `cnt()` returns `integer(0)` → the summary column fills with `0L`) and
the summary CSV is sorted by `-n_hq`, not `n_spx`, so the early file is already correctly ordered.
The post-SPARK call overwrites it with `n_spx` filled in. Cost: one extra write of `quant/` per
section. Nothing was deleted.

---

## CS-4. Runner ergonomics for resuming — `03_cytosignal.sh`

Added 2026-08-18, all additive; default behaviour with no arguments is byte-identical (verified:
4 significance sections + 1 differential, same command lines).

| flag | why |
|---|---|
| positional `SECTION …` | The script looped `$LUAD_SECTIONS` unconditionally, so resuming after a mid-stage failure would redo the completed `P17_AIS` **and overwrite its good output**. Now takes section arguments the way `01_stlearn.sh` already did. |
| `--skip-differential` | The NEBULA differential is one call over **all four** sections (`run_nebula_stage.R` hardcodes the list), so it cannot run while a subset is still missing — which is exactly the situation when resuming. |
| `--skip-export` | The four `input/` dirs already exist and check out against their own `provenance.json`. Re-exporting costs ~95 s and rewrites 3.4 GB of `counts.mtx` for nothing. |

**`--skip-export` was a no-op until 2026-08-21 — fixed.** As written on 2026-08-18 it gated only the
`step` banner:

```bash
[ "$SKIP_EXPORT" = "1" ] && step "skipping the CytoSignal input export (--skip-export)"
[ "$SKIP_EXPORT" = "1" ] || step "export the selected sections ..."
for S in $SECS; do run "$PY_PREP" .../export_cs_input.py ...; done   # <-- never gated
```

so it printed *"skipping"* and exported anyway. The sibling flags `--skip-significance` and
`--skip-differential` were proper `if/else` blocks; only this one used the `[ ] && / [ ] ||` idiom.

No number is compromised — `export_cs_input.py` is deterministic and the re-written bytes are
identical — but two things follow. It rewrote ~5.7 GB on every resume, and, more importantly,
**every file under `cytosignal/LUAD/input/` now carries an mtime LATER than the section runs that
consumed it** (MEASURED: `input/*` at 2026-08-20 12:51, after the runs finished). Anyone auditing
provenance by timestamp will read that as "the inputs changed after the analysis". They did not.
Now an `if/else` that also refuses when `$IN/$S` is absent, matching `04_cellchat.sh --skip-export`.

---

## CS-5. What is measured and what is modelled — read before sizing an allocation

Measured, and safe to quote: post-QC cell counts (`P17_AIS` 158,263 · `P21_AIS` 264,769 ·
`P21_LUAD` 498,422 · `P17_LUAD` 536,388); `P17_AIS` whole-pipeline wall clock ≈ 50 min;
`P17_LUAD` `inferSignif` = **6 h 42 min**; the 95.9989 GiB cap-hit; the interaction counts per slot
(912 / 171 / 171, panel-fixed).

**Modelled, and contested between analyses — do not quote:** the per-section peak-memory
projections. Two independent passes disagreed on SPARK-X's marginal cost (≈10 GiB vs ≈15-20 GiB),
on the score-matrix density (a claimed 1.0 was measured at 0.698-0.797), and on the neighbour term
(no fixed coefficient exists — it varies 483-816 across sections). The adversarial pass also noted
that anchoring any model on "96 GiB" is circular, because 96 GiB is the **cgroup ceiling**, not a
measured peak. Treat every projected GB figure as an upper-bounded guess and size generously.

`inferSignif` cost is **super-linear in cell count** (`graphSpatialFDRNew` materialises a dense
length-n vector per column, once per slot; a data-only bound gives an exponent ≥ 1.78). It is *not*
driven by the interaction count. Halving cells roughly quarters this stage.

---

## CS-6. Three of four LUAD sections have no figures — and `inferSignif` need not run again

`run_cytosignal.R:105` wraps **every** plotting call in one gate:

```r
## -------- plots (only when small enough; heavy at >200k cells) --------
if (ncol(cs@counts) <= 200000) { ... plotCluster ... plotSignif ... }
```

`quant/` is written before it and unconditionally, so a section over the threshold ships complete
quantitative output and zero images. Against the measured post-QC counts in CS-5:

| section | in | post-QC | ≤ 200,000 | figures on disk |
|---|---:|---:|---|---|
| `P17_AIS` | 182,378 | **158,263** | yes | `cluster_map.png` + `signif_plots/` (6) |
| `P21_AIS` | 292,862 | 264,769 | no | none |
| `P21_LUAD` | 560,183 | 498,422 | no | none |
| `P17_LUAD` | 640,739 | 536,388 | no | none |

Stage 03 also runs `nosave`, so there is no `cs_checkpoint.rds` / `cs_result.rds` to replot from.
This is declared behaviour of the tracked runner (`run_luad/03_cytosignal.sh:31-34`), not a failure.

**`inferSignif` does not have to run again to get the figures.** On the completed `P21_LUAD` run
(`cyto_P21_LUAD2.log`): findNN 5 min · imputeLR 1 min · `inferIntrScore` 43 min · **`inferSignif`
6 h 03 min** · SPARK 5.5 min · `quant/` writes 2 × 12 min — 7 h 23 min total, 82 % of it in one
stage. And that stage's entire output is `@res.list`, which `quant_io.R` already banks verbatim as
`quant/reslist_<tag>.rds` (`result` / `result.hq` / `result.spx`). Verified populated for all four
sections — `n_spx` is non-zero in every `signif_summary_*.csv` — and because the post-SPARK write
is the one that survives, the names carry the SPARK-X ordering.

`scripts/comparators/cytosignal/replot_from_quant.R` therefore rebuilds only as far as
`inferIntrScore` and **assigns the banked `res.list` back into the object**: ~50 min instead of
~7.5 h, and it also skips the stage where the 192 GB OOM occurred.

`inferIntrScore` is **not** skippable, and this is the reason: `cytosignal:::getIntrValue` — which
`plotSignif` → `plotIntrValue` go through — opens with

```r
sample.index <- sample(ncol(score.obj@lig.null), nCells)
```

unconditionally, and its default `type` includes `ligand_null` / `receptor_null` / `score_null`.
`@lig.null`, `@recep.null` and `@score.null` are produced by `inferIntrScore` and are not in
`quant/`. Only `plotIntrValue`/`plotSignif` need them; `plotEdge`, `plotSigCluster`,
`plotCircosNIntr` and `plotCluster` do not.

**Deviation, declared:** the permutation NULL shown in the `*_null` display panels is a fresh draw,
not the original run's. `inferIntrScore` permutes under `numCores = 4`, so bit-identical
reproduction was never on the table — `plot_signif_rerun.R`'s own header concedes as much. It
touches no statistic, ranking or significance call: every one of those is the authors' own
`inferSignif` + `rankIntrSpatialVar` output, read back unmodified. `CS_MODE=full` re-runs both
stages and removes even that, at the cost of the 6 h. Neither mode writes `quant/`.

Injection is verified before anything is plotted, and every check is **fatal**, because a silent
mismatch would yield confident wrong figures: (a) every banked barcode must exist among this
rebuild's post-QC cells; (b) `n_hq` must agree with `signif_summary_*.csv` on every row; (c) the
recomputed `@score` is compared against the banked one (`@score` is a deterministic function of the
imputation — only the null is permuted). Results land in `plots/injection_check.csv`.

Run it with `scripts/comparators/run_luad/sbatch/cytosignal_replot_<SECTION>.sbatch`
(`submit_all.sh` submits the three that have no figures; `--all` adds `P17_AIS`, which has only the
inline block's two plot types and never got the wider surface). `CS_PREFLIGHT_ONLY=1` validates
paths, env and `quant/` without spending the allocation.

---

## CS-7. `circlize` was never installed — `plotCircosNIntr` has been failing silently all along

`circlize` is in cytosignal's DESCRIPTION under **Suggests**, and the env recipe
(`_common/install_envs_iris.sh:62-67`, `build_cytosignal`) does not list it. `plotCircosNIntr` —
the chord diagram of interaction counts between cell types, and the only cross-cell-type summary in
CytoSignal's plot surface — therefore stopped with the package's own
`stop("Package 'circlize' is required for this function.")`.

Both `plot_signif_rerun.R` and `replot_from_quant.R` wrap that call in `tryCatch`, so it failed
**silently**. Confirmed by absence rather than by log: no circos output exists anywhere under
`results/`, including from the completed GBM run, whose completeness pass calls it.

Fixed additively by `scripts/comparators/cytosignal/install_plot_deps.sh` (idempotent;
`install.packages` into the env's own library against `cloud.r-project.org`, matching
`install_envs_iris.sh:69`). Installed **circlize 0.4.18** on 2026-08-20 and confirmed
`plotCircosNIntr` resolves.

Deliberately **not** folded into `install_envs_iris.sh`: that script records what actually built the
existing envs, and editing its recipe would make a fresh build diverge from the env every current
result came from.

The other Suggests the plotting path reaches — `cowplot` (`plotSignif`) and `scattermore`
(`plotIntrValue`, `raster = TRUE`) — are both present. `plot3D` / `png` / `pheatmap` are **not**
needed despite appearing in a naive grep: they occur only as string literals in
`plot.fmt = c("png","pdf","svg")`-style arguments, never as namespace calls.

---

## CS-8. `plot_signif_rerun.R` used a macOS `df` flag — fixed

The disk-space guard before the optional `cs_result.rds` write was:

```r
free_gb <- tryCatch(as.numeric(system(sprintf("df -g %s | ...", pdir), intern = TRUE)), error = ...)
if (!is.na(free_gb) && free_gb >= 20) { ... }
```

`df -g` is a BSD/macOS flag. On Linux GNU coreutils rejects it (`df: invalid option -- 'g'`), and
`system(intern = TRUE)` responds to a non-zero exit with a **warning, not an error** — so the
`tryCatch` never fires and `free_gb` becomes `numeric(0)`. The `if` then raises *"missing value
where TRUE/FALSE needed"*, killing the script **after every figure is written but before the
manifest**. Latent until now because the GBM run was on macOS. Fixed 2026-08-20 to `df -BG` with an
explicit length-0/NA guard; `replot_from_quant.R` avoids the construct entirely (`CS_SAVE_RDS`).

---

## CS-9. `print(plotIntrValue(...))` kept 1 panel of 8 — fixed in both scripts

`plotIntrValue` does not return one figure. Its last statement is

```r
names(plotList) <- type
outputList[[intrName]] <- plotList
```

so it returns a **nested** list: `outputList[[intrName]]` is itself eight ggplots, one per
`type` — `ligand`, `ligand_ori`, `ligand_null`, `receptor`, `receptor_ori`, `receptor_null`,
`score`, `score_null`.

`plot_signif_rerun.R`'s completeness pass did `print(plotIntrValue(...))` inside a single
`png()` device. Printing a list of ggplots draws each one in turn, so all eight went to the same
single-page device and **only the last (`score_null`) survived in the file** — the seven panels
that carry the actual ligand, receptor and score signal were silently discarded. No error, no
warning, a plausible-looking PNG.

Never caught because `plot_signif_rerun.R` has never been run on iris: there is no
`plot_signif_manifest.json` and no `intrValue_*` file anywhere under `results/`.

Both scripts now recover the inner list and combine it with `cowplot::plot_grid(plotlist = ...,
nrow = 2, labels = names(pl))`, at 2600 × 1400 rather than 1600 × 1200 so eight panels are legible.
`cowplot` was already a dependency of `plotSignif`, so nothing new is required.

Note this is the *opposite* failure mode to CS-7: circos failed loudly into a `tryCatch` and left
no file, whereas this one produced a file that looked fine. When auditing a plot suite, count the
panels, do not just count the files.

---

## CS-10. Sizing `--mem` for the replot — the bound is CS-1, not the full-pipeline peak

The first version of the replot sbatch files asked for **448 / 384 / 224 / 160 GB**. That was wrong
by roughly 4×, and wrong for an avoidable reason: it was anchored on the *full pipeline's* peak
(the 192 GB OOM, and the 256 GB allocation that completed `P21_LUAD`) even though the entire point
of `CS_MODE=inject` is that it does **not run** the stages that produce that peak.

The correct bound is CS-1. `P17_LUAD` — the **largest** section, 536,388 post-QC cells — completed
`findNN` + `imputeLR` + `inferIntrScore` + `inferSignif` inside a **96 GiB** cgroup, and was killed
only afterwards, inside `rankIntrSpatialVar` → `SPARK::sparkx`. Corroborated by the 192 GB run
(`results/_logs/cyto_P21_LUAD.log`, job 9051862): it finished `inferIntrScore` at 18:13:11, 50 min
in, and was OOM-killed at 18:34:15 — **21 minutes into `inferSignif`**, not before it.

`CS_MODE=inject` stops after `inferIntrScore`. It therefore never reaches **either** of the two
stages that have ever OOM-killed this pipeline, and its peak on the largest section is bounded
comfortably under 96 GiB. Current requests: **128 / 128 / 96 / 64 GB** for
`P17_LUAD` / `P21_LUAD` / `P21_AIS` / `P17_AIS`.

Do not treat those as final either — `replot_from_quant.R` logs its own peak RSS (`/proc/self/status`
`VmHWM`) at each stage and records it in `plots/replot_manifest.json`. Read it after the first run
and cut further. CS-5's "size generously" applies to the *full* pipeline, whose memory model is
genuinely contested; it does not license ignoring a measured bound for a shorter one.

---

## CS-11. `sbatch` runs a COPY of the script — `${BASH_SOURCE[0]}` does not point at the repo

The first three replot jobs (9177138-40) died in **3 seconds** with

```
bash: /var/spool/slurmd/job9177138/_replot_body.sh: No such file or directory
```

Slurm copies the batch file to `/var/spool/slurmd/job<id>/slurm_script` and executes it *there*, so
`$0` and `${BASH_SOURCE[0]}` resolve into the spool directory, not the checkout. The
`HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"` idiom that every other script in this repo
uses is correct for `bash script.sh` and **wrong for `sbatch script.sbatch`**. The `.sbatch` files
now use an absolute `REPO` (overridable with `sbatch --export=ALL,REPO=/other/checkout`).

**`sbatch --test-only` does not catch this.** It validates partition, account, QOS and resource
limits against the scheduler and reports when the job would start — it never executes a line of the
script, so it passed on all four files while every one of them was broken. To actually verify a
batch script end to end without spending an allocation, submit it for real in preflight mode:

```bash
sbatch --export=ALL,CS_PREFLIGHT_ONLY=1 --mem=8G --time=00:05:00 \
       scripts/comparators/run_luad/sbatch/cytosignal_replot_P21_AIS.sbatch
```

`_replot_body.sh` then resolves the env, the interpreter, the input dir, the DB and `quant/`,
prints them, and exits 0 in seconds. Verified as job 9178866 (`COMPLETED`, `0:0`, 9 s).

This is the first `sbatch` in the repo; every earlier comparator run was an interactive `srun`,
which is why the trap had not been hit before.

---

## CS-12. `CS_EXTRA` defaulted to an LGG hypothesis and leaked it into all four LUAD sections — fixed

`replot_from_quant.R` took an env knob `CS_EXTRA`, a comma-separated list of CCI ids to plot
into `plots/requested_<tag>/` regardless of rank. Its default was **`CCI-01109` = GRN → SORT1**,
and **nothing anywhere in the repo ever set `CS_EXTRA`** (verified by grep over `scripts/`). So
every LUAD/AIS replot silently produced:

```
cytosignal/LUAD/cellchatdb2/P17_AIS/plots/requested_diffusion_Raw_smooth/Rank_259_GRN-SORT1.png
                            P17_LUAD/...                                  Rank_149_GRN-SORT1.png
                            P21_AIS/...                                   Rank_144_GRN-SORT1.png
                            P21_LUAD/...                                  Rank_318_GRN-SORT1.png
```

**Why this is wrong.** `GRN_SORT1` and `ANXA1_FPR1` are the two arms of ALARMIST motif 1 in the
**LGG/GBM** TMA. They are a hypothesis about that dataset. Naming one of them on a LUAD/AIS run
imports a hypothesis this data was never designed to test and makes the output read as a targeted
test nobody specified.

The old comment justified the default like this: *"GRN → SORT1, the arm of ALARMIST motif 1 that
is testable here … ANXA1 is absent from the LUAD Xenium 5K panel (FPR1, GRN and SORT1 are all
present), so it is not merely non-significant, it was never testable."* **That reasoning is wrong
and was rejected 2026-08-24.** Panel membership means the assay *could* measure the pair; it says
nothing about whether the hypothesis belongs to the dataset. The comment has been removed so it
cannot mislead the next reader — it already misled one, who went on to propose repointing
`run_cellchat_luad.R`'s `REQUESTED_LR` at `GRN_SORT1` for exactly the same bad reason.

**Fixed in three layers, 2026-08-24:**

1. `replot_from_quant.R:97` — `CS_EXTRA` now defaults to **empty**. A dataset-agnostic script must
   not carry any dataset's hypothesis as a default. A GBM run that wants it passes
   `CS_EXTRA=CCI-01109` explicitly. Safe to change: the script is untracked, LUAD-era, and no GBM
   result was ever produced by it (no `replot_manifest.json` outside the LUAD tree).
2. `_common/luad_config.sh` — a documented **INVARIANT** block plus `CS_EXTRA=""; export CS_EXTRA`.
   Every LUAD stage sources this file, so an inherited `CS_EXTRA=CCI-01109` is forced back to empty
   (verified end to end). Any future knob of this kind defaults to empty and gets added here.
3. `03_cytosignal.sh` — `export CS_EXTRA=""` at the call site as well.

**The other four runners were already correct** and are recorded here so nobody "fixes" them:
`01_stlearn.sh` and `02_spatialdm.sh` pass `--requested-lrs ""`; `05_liana.sh:137` passes a bare
`--lr-of-interest` immediately followed by `--tag`, which argparse (`nargs="*"`) resolves to `[]`
(verified); `run_cellchat_luad.R` sets `REQUESTED_LR <- character(0)` via `make_luad_variants.sh`.

**Therefore, on a LUAD run these are CORRECT OUTPUT, not defects:**
`quant/requested_lr_status.csv` 3 bytes · `quant/requested_lr_in_DEA.csv` 3 bytes ·
`mofaflex_inflow_joint/data/lr_of_interest_loadings.csv` 1 byte · `plots/requested*/` empty or
absent. An audit that flags these as silent failures is wrong; do not repopulate them.

**Still on disk:** the four `Rank_N_GRN-SORT1.png` files listed above, from replots run before the
fix. They are stale artefacts of the old default. Deleting them is a results change and needs the
owner's say-so.
