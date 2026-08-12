# CLAUDE.md — project instructions for Claude Code

This repo wears two hats: it is the **ALARMIST method** (the `src/alarmist` Python
package + the nf-core Nextflow pipeline) *and* the **analysis repo** where that method
is applied to real datasets (`scripts/research/`, `results/`, `reports/`). Package code
gets software standards; analysis code gets reproducible-science standards. Say which
one you are in when it matters.

ALARMIST = **A**ssessment of **L**igand **A**nd **R**eceptor Interaction **M**otifs
**I**n **S**patial **T**ranscriptomics — LRI motif discovery by Bayesian Poisson tensor
factorization on a patch × LRI matrix, projected to single cells.

## Announce every skill you use

**Whenever you invoke or follow a skill, tell the user which one, explicitly, at the
point you start using it** — e.g. "using `/alarmist`" or "applying
`/nature_publication_figures` for the figure". If several apply, name all of them. If a
skill looks relevant but you decide it does *not* apply, say that too, and why. Never
apply a skill's rules silently — the user needs to know which playbook produced an
answer in order to trust or correct it.

## Working rules

- Treat this as reproducible science: correctness, provenance, and readability for
  future scientific review matter more than cleverness.
- **Do not silently guess** when project structure, paths, dataset parameters, or data
  provenance are ambiguous — ask. If a dataset is not pinned in the table below, its
  patch size / `K` / cell-type column is *not* known; ask instead of inferring.
- **Ask before risky changes** (deleting/overwriting data or results, rewriting shared
  logic in `src/alarmist/`, changing pipeline semantics).
- **Do not run long analysis jobs** unless explicitly asked. BPTF runs are minutes to
  hours; run them in the background and summarise.
- **Do not copy large `.h5ad` files or old results into the repo** unless asked.

## Environment

Every Python command for this repo runs in the `bptf` conda env, invoked **by absolute
path to the interpreter**:

```bash
/Users/jiayifan/anaconda3/envs/bptf/bin/python run_pipeline_more_plots.py ...
```

Never base Python (`Numba needs NumPy 1.24 or less`, `ModuleNotFoundError: alarmist`).

> **Conda on this machine is degraded.** `source "$(conda info --base)/etc/profile.d/conda.sh"`
> fails (permission denied on `__conda_exe`), `libarchive.13.dylib` is missing so
> libmamba/mamba fail, and `conda activate` / `conda run` silently fall through to the
> *system* interpreter rather than the env's. So: always call the env's binary directly.
> For `conda create` / `install`, prefix `CONDA_SOLVER=classic` and use
> `/Users/jiayifan/anaconda3/bin/conda ... -c conda-forge --override-channels`.

**Comparator methods** (the CCC benchmark) each get their **own** isolated env named
`comp-<method>` — never co-install two methods, their scanpy/numpy/Seurat pins conflict.
ALARMIST itself stays in `bptf`. R comparators are driven via
`source scripts/comparators/<method>/activate_env.sh` (plain `conda activate` gives you
system R). Runners live in `scripts/comparators/<method>/`, outputs in
`results/comparators/<method>/`.

## Running ALARMIST

`run_pipeline_more_plots.py` (repo root) is the **canonical, dataset-agnostic driver**.

> ⚠️ **Two drivers sit in the repo root and only one is current.** `run_pipeline.py` is the
> **stale** original (9 arguments, hardcoded paths) and is still tracked; use
> `run_pipeline_more_plots.py` (39 arguments) for everything. The `/alarmist` skill's
> `SKILL.md:74` and `references/running.md:18` still name `run_pipeline.py` — that skill is a
> synced copy of `~/tansey_lab/spatial_analysis_skills/skills/alarmist` and is deliberately
> **not** edited here, so fix it at the source if you fix it at all. *(The canonical driver was
> called `run_pipeline_new.py` until 2026-08-12; every reference in this file was updated with
> the rename, but a shell history or a `results/*/run_manifest.txt` may still carry the old name.)*

Stages: patchify → bptf → plots-bptf → project →
plots-project → markers → glm → plots-glm, with presets (`full`/`bptf`/`plots`),
`--from`/`--to`, sentinel-based resume, `--dry-run`, `--force`, and patchify reuse.
Nothing is hardcoded: `--data-file` is required, `--cellchatdb` optional, interpreter
defaults to `sys.executable`. Two figure families are produced: `plots_original/` (the
`tutorials/GBM.ipynb` suite) and `plots/` (stock `alarmist.cli.visualize`). The
slide-tags dataset has its own driver, `scripts/research/slide_tags/run_pipeline_adapt.sh`
(+ `prepare_input.py`, which collapses the 16 neuron subtypes).

Requirements and traps:

- **Raw counts must be in both `X` and `layers['counts']`** — patchify binarizes `X`, the
  GLM reads `layers['counts']`, marker DE reads `X`. If `X` is log-normalized, pass
  `--count-layer layers:counts`.
- **Fix seeds:** `np.random.seed` before BPTF (BPTF has no `random_state` argument —
  reproducibility rides on this), `random_state` on the GMM.
- **Non-human data hangs unless you pass the right LRI DB.** `alarmist-bptf`'s
  `process_bptf_results()` annotates pathways against a hardcoded *human* CellChatDB, and
  `annotate_pathways` has an O(rows × pairs) fuzzy fallback — on mouse data nothing
  matches, every row takes the slow path (12+ min, all "Unknown"). `run_pipeline_more_plots.py`
  works around it by running BPTF in-process and forwarding `cellchatdb_path`; if you
  drive the CLI yourself, don't.
- **Legends duplicate multi-word cell types.** The plot helpers loop every key of the
  `ct_colors` dict, and `with_sanitized_aliases` adds an underscored alias per name
  (needed only for GLM volcano/forest lookups), so "White Blood Cells" appears twice.
  Pass the clean palette to bptf/project figures; add aliases only around GLM figures.
  Cosmetic only — colours and data are unaffected.
- Motif indices are **arbitrary across runs**. Any cross-condition comparison must come
  from a single joint run.

## LRI databases

`CellChatDBv2.0.{human,mouse}.csv` exist in **two** places —
`src/alarmist/config/lri_databases/` is what the package loads via
`_get_bundled_database_path`; `data/LRdatabase/` is the working copy and the archive.

**Sync state (verified 2026-07-30):** every file present in both is byte-identical
(`CellChatDBv2.0.human.csv`, `CellChatDBv2.0.human.old.csv`, `CellChatDBv2.0.mouse.csv`,
`CellPhoneDBv5.0.human.csv`). The bundled directory is a strict **subset** — these live
only in `data/LRdatabase/`: `CellChatDBv2.0.mouse.old.csv`, `CellPhoneDBv4.0.human.csv`,
`CellPhoneDBv5.0.human_old.csv`, `cellchatv2.csv`, `signaling_type_difference.csv`.
After re-exporting a DB, copy it into **both**.

Two packaging notes: `_get_bundled_database_path` is only ever called with the two
**human** names, so the mouse DB is never auto-resolved — you must pass `--cellchatdb` /
`cellchatdb_path` explicitly (this is what the non-human hang above is about). And
`pyproject.toml` ships `config/lri_databases/*.csv` by glob, so the pre-fix
`CellChatDBv2.0.human.old.csv` currently goes into the wheel — the archive belongs in
`data/LRdatabase/` only.

They are exported from the CellChat R package v2.2.0. **Derive subunits from the
`complex` table, never from `ligand.symbol` / `receptor.symbol`** — the `*.symbol`
columns are display-only and contain errors (`CXCL8_CXCR1` has
`ligand.symbol == "CXCL1"`, which yields a duplicate `CXCL1|CXCR1` row and loses
CXCL8→CXCR1 entirely; `CD99_PILRB` carries mouse casing in the human DB). Mapping:
`ligand`/`receptor` ← `complex` subunits joined with `_`; `pathway` ← `pathway_name`;
`signaling_type` ← `annotation`; `version` ← `version`. Current state (re-exported
2026-07-28): human 3233 rows / 3218 unique pairs / 1383 genes, mouse 3379 rows / 3359
unique pairs / 1439 genes, all genes validating against `geneInfo$Symbol`; pre-fix
originals kept as `*.old.csv`.

**The re-export broke join-compatibility with every existing result — read this before
rerunning or joining anything.** `CellChatDBv2.0.human.old.csv` is byte-identical to
`git show HEAD:data/LRdatabase/CellChatDBv2.0.human.csv` (md5 `ca044e67…`), i.e. it is the
DB that produced every result currently on disk. The 2026-07-28 re-export reorders complex
subunits (`RAMP2_CALCR` → `CALCR_RAMP2`, `TGFBR2_TGFBR1` → `TGFBR1_TGFBR2`) and turns some
ligands into homo-complexes (`ALDH1A1` → `ALDH1A1_ALDH1A1`): **1,120 of 3,218 LR keys
changed; only 2,097 are shared.** Stored LRI column names use the old ordering
(`results/spatch/coad_xenium/patch_lri_columns.csv`: 546 × `TGFBR2_TGFBR1`, 0 ×
`TGFBR1_TGFBR2`). So: to reproduce or to join against saved matrices, point
`cellchatdb_path` at `data/LRdatabase/CellChatDBv2.0.human.old.csv`; use the current file
only for fresh runs, and never mix results across the two.

Three mouse rows (`Rarres2→Cmklr2`, `H2-D→Klra`, `H2-K1→Klra`) resolve to symbols real
data never uses; they are kept for fidelity and simply never match. `_split_gene_complex`
in `src/alarmist/core/lri.py` now strips whitespace (fixed 2026-07-28) — before that a
`", "`-formatted complex silently dropped the whole LRI (776 → 1087 LRIs on the human
Xenium 5K panel).

## Datasets — the facts the skills defer here

The `/alarmist`, `/spatial-*` and `/cohort-explore` skills are deliberately
dataset-agnostic and send you here for patch size `s`, motif count `K`, cell-type column,
LRI resource, and paths.

**The manuscript lives in a sibling repo: `/Users/jiayifan/tansey_lab/alarmist_tex`**
(`results.tex`, `methods.tex`, `discussion.tex`, `supplement/`, figure assets in `fig/`).
The datasets below are ordered **as the paper introduces them**, then the analyses that
are not in the paper. Parameters marked *(paper)* come from the manuscript; everything
else was read off disk. Where the two disagree, the disk wins and the conflict is called
out. **Anything not listed here is not pinned — ask, don't infer.**

### Fig 2 — semi-synthetic benchmark (`fig:benchmark`)

Not a repo dataset and **not reproducible from this repo**. Ten ground-truth motifs are
simulated on top of a real COAD scRNA-seq reference from the SPATCH study (8,288 cells,
16 types, `major_annotation`) — the reference h5ad lives in the *sibling* project at
`/Users/jiayifan/tansey_lab/ebb_flow/data/coad_sc.h5ad`, not under `data/`. K=20,
**CellPhoneDB** (the one analysis in the paper that does not use CellChatDB), patch
`50 × 50` in synthetic units (not µm), patch-count sweep {100, 400, 1000, 5000, 10000} at
~25 cells/patch. `methods.tex:168` states the simulation scripts are provided separately —
they, the ground-truth motif definitions, and the comparator runs (Tensor-cell2cell,
COMPOTES, NICHES, Leiden, Novae) are **absent here**. Do not try to reconstruct them from
`scripts/research/`; get them from whoever ran them.

### Fig 3 — SPATCH cross-platform concordance (`fig:fig3`)

Matched Xenium Prime 5K + CosMx 6K consecutive sections from one COAD, one OV and one HCC
specimen (6 sections), from the SPATCH portal (Ren et al. 2025).

- **Data:** `data/spatch/{coad,ov,hcc}_{xenium,cosmx}.h5ad`. Xenium 5,001 genes
  (405,927 / 410,193 / 271,687 cells); CosMx 6,175 genes (292,371 / 289,272 / 237,030).
  `X` is raw counts in all six; no `layers`, no `.raw`. Coordinates are already µm.
- **Cell types:** `obs['annotation']` — the run scripts **copy it to `obs['cell_type']`**
  after subsetting to `obs['high_quality'] == 1`. COAD 16 types, HCC 17, OV 13; identical
  vocabulary across platforms within a cancer, different across cancers.
- **Params:** patch **50 µm**, **K=20**, human CellChatDB, seed 42, `sample_column` **none**
  — all six are *independent* fits.
- **Results:** `results/spatch/{coad,ov,hcc}_{xenium,cosmx}/` (per-section runs), plus the
  comparison layer `results/spatch/{correlation,leesl,morans_i,motif_alignment_new,motif_alignment_{coad,hcc,ov}}/`.
  Aligned CosMx copies `data/spatch/*_cosmx_aligned.h5ad` carry
  `obsm['spatial_aligned_to_xenium']` + `uns['alignment_info']`; **the ALARMIST fits used
  the unaligned originals** — alignment matters only for the downstream Lee's L / niche grids.
- **Scripts:** `scripts/research/xenium_cosmx_{align,run1,run2,run3,correlation,motif_alignment,overlap_boxplot,lees_l,morans_i,niche_concordance}.py`.
  In `run1.py` the load/patchify/BPTF block is now **commented out** (only the plotting
  tail is live) — read parameters from `analysis_parameters.csv`, not from the script.
- **Trap:** motif index *k* on Xenium has no relation to motif *k* on CosMx. Every
  cross-platform comparison goes through `results/spatch/correlation/{disease}_cosine.csv`
  with the cosine > 0.4 matching rule *(paper)*.

### Fig 4 — AIS → LUAD progression (`fig:fig4`)

- **Data:** four Xenium 5K sections, `data/linghua/P{17,21}_{AIS,LUAD}_Xenium.h5ad` —
  matched precursor AIS + invasive LUAD from two patients. *(paper: >1.6M cells,
  5,096-gene panel, 19 cell types, 772 detectable LRIs.)* The other `data/linghua/` files
  (P4, P11, P24 AAH/LUAD) are **not** part of this run.
- **Cell types:** `obs['annotation_coarse']`, renamed to `cell_type` at load.
- **Params (on disk, `results/AIS_LUAD/analysis_parameters.csv`):** patch **80 µm** — not
  the 50 µm default, because lung tissue is sparse *(results.tex:53)* — `n_samples=4`,
  58,203 patches, 181,424 LRI columns, 5,101 shared genes, CellChatDB.
  **K=25** (`bptf*/factorization_parameters.csv`).
- **Results:** `results/AIS_LUAD/`. Three BPTF dirs, all K=25 on the same 58,203 patches:
  `bptf_25` (Dec 2025) → `bptf` (Jan 2026) → **`bptf_new` (Feb 2026) = canonical**, with
  `bptf_plots/` the matching figures. **Motif 10 = tumor vasculature, motif 24 = healthy
  vasculature** (`de_results_10.csv`, `de_results_24.csv`, `volcano_motif{10,24}_celltype.svg`).
- **Trap:** `scripts/research/luad_1.py` is a copy-paste from the SAHA scripts and still
  carries stale `n_components = 20` and `patch_size_px = 50/0.12028` lines that do **not**
  match the run. Read `analysis_parameters.csv`, never that script.
- **Not the paper run:** `results/AIS_LUAD_new/` is a July 2026 single-sample
  (P21_LUAD only) torch-backend rerun at seed 53 — see its `run_manifest.txt`.
  `results/AIS_LUAD_plot/` is a replot-only dir.

### Fig 5 — IDH-mutant LGG / GBM TMA (`fig:fig5`)

The paper calls this **LGG** (IDH-mutant low-grade glioma with active high-grade
transformation); the repo calls it **GBM**. Same dataset.

- **Data:** `data/xenium_mm_final_cell_id.h5ad` — 100,197 cells × 5,119 genes, human
  Xenium TMA. `X` is log-normalized, raw counts in `layers['counts']`.
  *(This file was recorded as `data/z.h5ad` until 2026-08-12; no such file exists or ever
  did on this machine. Verified: the path above has exactly the shape described here, and
  it is what all 14 comparator runners and `results/GBM/` were produced from.)*
- **Cell types:** `obs['cell_type']` (9 types). Also present: `old_cell_type` (11),
  `subtype` (9), `grade` (high/low), `tma_id`, `motif` (see the trap below).
- **Params:** patch **50 µm**, **K=20**, human CellChatDB
  (`results/GBM/analysis_parameters.csv`: 13,113 patches × 25,271 LRI columns).
  *(paper: 748 candidate LRIs; CellTypist trained on the GBmap core atlas at annotation
  level 2; mGAM called by a FOSL2 module score > 0.25; tumour states from Neftel 2019.)*
- **Results:** `results/GBM/` (`bptf/`, `single_cell/`, `impact/`, `markers/`,
  `bptf_plots_20/`). `results/GBM_*` siblings are reruns/legacy — read the dir name.
- **Core count — the paper is wrong.** `obs['tma_id']` has **13** cores, **7 high-grade +
  6 low-grade** (verified by crosstab; ids run 1–14 with **7 absent**). Fig 5a's caption
  ("n=6 low-grade and n=7 high-grade") is right; `methods.tex:2` ("14 punch cores, 7
  low-grade and 7 high-grade") is wrong on both numbers.
- **Fig 5i** is external TCGA LGG survival on a 20-gene signature (top 20 by log2FC,
  FDR 0.001) — no methods subsection covers it and no code for it is in this repo.
- **Motif loadings — the one that bites.** Per-cell loadings `U` are
  `results/GBM/single_cell/cell_loadings.npy`, shape **(100197, 20)**, 0-based columns, in
  the same cell order as the h5ad. **Motif 1 (mGAM) = column 1.** Do **not** use
  `obs['motif']` — that is a 15-category argmax label, not the 20-motif loadings. Motif 1
  is a bidirectional mGAM ⇄ MES-like loop: GRN→SORT1 (mGAM→MES-like) and ANXA1→FPR1
  (MES-like→mGAM), where GRN and FPR1 are mGAM-specific by expression while SORT1/ANXA1
  are broad. It is grade-associated (`results/GBM/motif1_grade.svg`).

### SAHA — IBD + ileum (**not in the paper**)

A separate CosMx story: does motif usage separate inflamed IBD from ileum, do motifs
replicate across the two tissues, and does a treatment responder differ from a
non-responder? Written up in `reports/saha_ibd_ile_report*.html`.

- **Data lives outside the repo:** `/Users/jiayifan/tansey_lab/saha_ibd/SAHA_{IBD,ILE}_RNA.h5ad`.
- **Sample key:** `obs['section_ID']` — patchify is run in **multi-sample dict mode**
  (`{section: adata}`), 8 IBD sections and 12 ILE sections (`*/sample_info.csv`).
- **Cell-type column differs per run:** IBD → `cell_type`; ILE → `Insitutype_Broad`
  (copied to `cell_type`); `SAHA_IBD_more` → `cell_type_general`, which is why it has
  365,742 LRI columns against 120,971 for the same 8 sections.
- **Patch size is in CosMx pixels, not µm:** `patch_size = 50/0.12028 ≈ 415.697` px = 50 µm
  (CosMx pitch 0.12028 µm/px). Same class of unit trap as MOSTA — do not pass 50.
- **K:** IBD swept 15/20/25 (`bptf_15|20|25`), **20 is the working choice**; ILE and
  `IBD_more` are K=20 only.
- **Results:** `results/SAHA_IBD/`, `results/SAHA_ILE/`, `results/SAHA_IBD_more/`;
  cross-tissue figures in `figures/`; responder summary
  `results/saha_ibd_responder_summary.json` (**responder `IBD_A2_1_EC04`** 35,561 cells vs
  **non-responder `IBD_A2_2_EC04`** 18,124, compared by per-motif positive fraction against
  the 374,881-cell IBD background); motif digest `results/saha_motif_summary.json`.
- **Niche × motif work** is its own directory: `saha_niche_corr/` (`compute_niche_motif_corr.py`,
  `explore_annotation_consistency.py`, `build_niche_report.py`). Note it also audits the
  several competing cell-type vocabularies (`cell_type`, `cell_type_general`,
  `Insitutype_Broad`, `Insitutype_Labelled`, `ct_major`) — check which one a result used.
- **Scripts:** `scripts/research/saha_{ibd,ile,ibd_more}_{1,2,3}.py` (1 = patchify+BPTF,
  2 = projection, 3 = downstream) and `build_saha_report.py`.

### MOSTA Stereo-seq (**not in the paper**)

`data/mosta/*.h5ad` (e.g. `E16.5_E1S3_cell_bin.h5ad`): 281,377 cells × 28,103 genes, mouse
embryo. Cell-type column **`annotation`** (25 types), *not* `cell_type`. Coordinates are
**0.5 µm per unit** (Stereo-seq DNB pitch), so `--patch-size 100` = 50 µm — native units,
not µm. K=20, **mouse** DB, `X` log-norm with raw in `layers['counts']`.
`obs/annotation` uses legacy anndata categorical encoding (int8 codes with labels under
`obs/__categories/annotation`) — a naive h5py read yields `'0','1',…`; scanpy reconstructs
it correctly. At 50 µm the E16.5 embryo gives ~29k patches × ~191k LRI columns. CNGB
downloads can arrive truncated — verify the file opens first.

```bash
/Users/jiayifan/anaconda3/envs/bptf/bin/python run_pipeline_more_plots.py \
  --data-file data/mosta/E16.5_E1S3_cell_bin.h5ad \
  --cellchatdb data/LRdatabase/CellChatDBv2.0.mouse.csv \
  --output-dir results/mosta/E16.5_E1S3_50um_k20 \
  --cell-type-column annotation --n-components 20 --patch-size 100 \
  --seed 0 --preset bptf --network-threshold 50 100
```

### Everything else

`data/` and `results/` also hold ES, AAH, slide-tags, slide-seq and other runs. **None of
their parameters are pinned here — ask.** Their run directories encode parameters in the
name (e.g. `results/slide_tags/run_seed0_15_100um/`); read the directory and the run's
`analysis_parameters.csv`, don't assume.

## Comparator benchmark — Figures 6 and 7

**The question:** ALARMIST motif 1 is a bidirectional mGAM ⇄ MES-like loop
(GRN→SORT1 + ANXA1→FPR1). **Can any competing cell–cell-communication method find it on its
own?** Seven were run on the GBM/LGG TMA, each following *its own authors' default workflow*
— CytoSignal, stLearn, SpatialDM, COMMOT, LIANA+ (8 branches), NICHES, CellChat.

**Start here: `scripts/comparators/METHODS.md`.** It is the living document — one section
per method (algorithm, spatial model, LR database, inputs, workflow, every output, gotchas,
deviations, the numbers) plus the cross-method reading rules and a consolidated *Open issues*
list. `/comparator-benchmark` is the skill that governs how a method is added.

The four results, in one paragraph: **every method detects both arms** and all of them rank
GRN→SORT1 above ANXA1→FPR1; **none puts the two arms in one data-derived object**; the reason
is the **unit of analysis** — the two arms co-occur at Pearson 0.018 per cell but 0.456 per
50 µm patch, a 26× rise from aggregation alone; and **CellChat independently recovers both
arms, significant only in high grade**, corroborating motif 1's grade association. Nearly
every grade test is null because **grade and cellularity are collinear** in this 13-core TMA
(ρ = 0.78) and a 7-vs-6 rank test floors at p = 0.0012.

Working rules specific to this benchmark:

- **One conda env per method**, `comp-<method>`, never shared — their scanpy/numpy/Seurat
  pins conflict. R methods need `source scripts/comparators/<method>/activate_env.sh`,
  because degraded conda makes `conda activate` fall through to *system* R.
- **Never harmonize spatial scale across methods** (`SKILL.md:45-46`). Each method keeps its
  own kernel default, so **LR-pair counts are not comparable across methods** — compare
  fractions or ranks. LIANA+'s bandwidth violates this and is logged as open deviation CD-1.
- Code in `scripts/comparators/`, outputs in `results/comparators/<method>/<dataset>/<tier>/`
  (gitignored, ~44 GB). `<dataset>` is `GBM` or `LUAD`; `<tier>` is `default` or `cellchatdb2`.
- All figures go through `scripts/comparators/_common/plotting.py`. **Do not re-declare
  rcParams or a per-format saver in a new script** — seven scripts once did.

**What is not done** (full list in METHODS.md § Open issues): **LUAD is essentially untouched**
— 0 of 7 methods completed it, though the skill's invariants name it; the `default` tier
exists for only 3 of 7; LIANA's CD-1/CD-2 await sign-off; the COMMOT downstream-impact run is
partial. Where the `default` tier *was* run it showed the LR resource changes which pairs are
tested and the ranking denominator but **not any per-pair statistic**, so the missing four are
a completeness gap, not a suspected confound.

## Plotting

Every publication figure: **Arial**, editable vector text (`pdf.fonttype = 42`,
`svg.fonttype = 'none'`), and saved as **`.png` + `.pdf` + `.svg` through a single
saver** — never hand-roll per-format `fig.savefig`, and never save only one format
unless the user explicitly asks.

There is currently **no canonical saver** in `src/alarmist/plotting/` (the repo has ~174
bare `savefig` calls and each research script re-declares its own rcParams). Until one
exists: define a small `save_all_formats(fig, path_without_ext)` + style setup at the top
of the script and use it consistently. If you are writing the third copy of it, **ask**
before promoting it into `alarmist.plotting.utils` — don't silently grow the package API.

**Exception — the comparator tree already has one.** `scripts/comparators/_common/plotting.py`
(added 2026-08-12) holds `apply_publication_style(**overrides)` and
`save_all_formats(fig, stem, dpi=…, close=…, verbose=…)`; seven scripts under
`scripts/comparators/` had each re-declared their own, with three different saver
signatures and three different default dpi. **Anything new under `scripts/comparators/`
imports from there — do not write an eighth copy.** It is deliberately *not* in the
package, per the rule above. Bootstrap from a method subdirectory with:

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/comparators
from _common.plotting import apply_publication_style, save_all_formats
```

Full layout/legend/statistics/export conventions → the `/nature_publication_figures`
skill. HTML deliverables → `/interactive-report`.

## Repo organization

- `src/alarmist/` — the package (`core/`, `plotting/`, `cli/`, `config/`, `data/`).
- `nextflow/` — the nf-core pipeline.
- `run_pipeline_more_plots.py` — the general driver (see above).
- `scripts/research/` — one-off analysis scripts, **locally excluded from git**
  (`.git/info/exclude`). Treat them as **patterns to adapt, not an API**.
- `scripts/comparators/` — the CCC benchmark (see its own section above). `METHODS.md` is the
  living document; `<method>/` holds one env's runners plus its `NOTES.md` (tutorial call
  contract) and `DEVIATIONS.md`; `_common/plotting.py` is the shared figure saver;
  `figure6_*.py` / `figure7_*.py` build the paper panels from what is already on disk;
  `_archive/` holds superseded scripts and the pre-2026-08-12 `METHODS.md`.
- `tutorials/GBM.ipynb` — the reference notebook the "original plots" suite came from.
- **`/Users/jiayifan/tansey_lab/alarmist_tex`** — the manuscript (sibling repo, not part of
  this one). `/Users/jiayifan/tansey_lab/ebb_flow` holds the benchmark's scRNA-seq reference.
- `data/`, `results/`, `figures/`, `reports/` — gitignored. `*.png`/`*.pdf`/`*.svg` are
  gitignored repo-wide, so a figure you "saved" is not committed; the **script** is the
  tracked artifact.

## Skills in this repo

Registered under `.claude/skills/`. Only **`alarmist`** and **`comparator-benchmark`** are
tracked in git — this repo owns both; the rest are vendored copies and are gitignored (see
`.claude/skills/VENDORED.md` for provenance and licences).

| Skill                          | Owns                                                                                                                                                                                                       |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `alarmist`                   | Interpreting a fitted run: characterize a motif, screen motifs, test condition enrichment, motif→gene impact, signatures, motif niches, reproducibility.`references/running.md` covers producing a run. |
| `comparator-benchmark`       | Running a competing CCC method on our data**on the authors' own terms**, and writing it up. Owns the one-env-per-method rule, the two-tier LR-database design, the do-not-harmonize-kernel-scale invariant, and the`METHODS.md` section template.                     |
| `spatial-workflow`           | Pipeline order, the replicate-unit rule, the silent-failure traps.                                                                                                                                         |
| `spatial-niche`              | Niche clustering / annotation / association / between-condition comparison.                                                                                                                                |
| `spatial-stats`              | Moran's I, Lee's L, landmark alignment.                                                                                                                                                                    |
| `cohort-explore`             | Turning per-sample metadata into scoped, testable comparisons.                                                                                                                                             |
| `scientific-goals`           | Framing an investigation, the hypothesis ledger, dispatching per-hypothesis work.                                                                                                                          |
| `nature_publication_figures` | Figure layout, legends, statistics annotation, vector export.                                                                                                                                              |
| `interactive-report`         | Self-contained HTML reports and interactive explorers.                                                                                                                                                     |

Deliberately **not** vendored: `es-figures` (depends on `es_utils` helpers this repo does
not have — its rules are summarised in *Plotting* above) and `es-de-units` (Tazemetostat /
TMA-specific; its generic core, "the unit is the replicate, cell-level p-values are
pseudoreplication", is owned by `spatial-workflow`).

## Before finishing a task

Always summarize:

1. what files were changed;
2. what was **not** changed;
3. what commands were run;
4. whether any commands failed;
5. any assumptions or unresolved questions.
