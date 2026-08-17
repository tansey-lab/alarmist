**run CytoSignal, stLearn, spatialDM, COMMOT, LIANA+, NICHES, CellChat on the LUAD and LGG datasets to show what they find and if it's possible to see the same signal or if they're deficient somehow. ultimate goal: have a key result plot convincing the reviewers that other methods dont work or cannot find what alarmist can do.**

---

# Comparators plan (competing CCC methods vs ALARMIST)

## 0. Pin the reference signal first (do before any competitor run)
For each dataset, freeze the ALARMIST result the competitors are measured against:
- **GBM** (`data/gbm_tma1.h5ad`, `results/GBM/`): the grade-associated motif(s)
  (`motif1_grade.svg`, `all_motifs_grade.svg`) + the mGAM→MES-like directional motif.
- **LUAD** (`data/linghua/P*_LUAD_Xenium.h5ad`, `results/LUAD/`): pick the headline motif
  (top-activity motif from `bptf_25`) — TODO: confirm which sample(s) the ALARMIST run used.
- Write down, per reference motif: its defining LRIs, sender→receiver cell-type pairs,
  spatial ON region, downstream impact genes, and condition (grade/stage) association.
  This tuple is the target every competitor is scored against.

## 1. Comparison axes (the scorecard — columns of the key figure)
| axis | what ALARMIST does that we test |
|---|---|
| A spatial | uses coordinates, not cell-type-aggregate only |
| B resolution | single-cell / fine output |
| C **co-occurrence** | discovers a *program* of many co-active LRIs (not one LR at a time) |
| D **direction** | sender→receiver cell-type pairs *within* a program |
| E **motif→gene** | downstream transcriptional impact (GLM) |
| F **states** | per-cell ON/OFF state you can map + test |
| G condition | program associates with grade/stage |
| H scale | runs on Xenium (1e5–1e6 cells) |

## 2. Methods — output type, env, expected gap
| method | lang | output | closest to ALARMIST / gap |
|---|---|---|---|
| CellChat v2 | R | cell-type pathway comms (+spatial) | has D; no C/E/F, no per-patch program |
| LIANA+ | Py | consensus LR + spatial bivariate + NMF/MOFA factors | strongest on C; no D/E, factors over LR scores not directional counts |
| NICHES | R | interaction matrices → embed/cluster | multivariate but clusters of pairwise vecs; no D/E/F |
| COMMOT | Py | OT spatial signaling per pathway | strong A/B; pathway-by-pathway, no C |
| CytoSignal | R | single-cell spatial LR scores | strong A/B; per-LR, no C/D/E |
| stLearn CCI | Py | LR hotspots + celltype-pair perm test | per-LR/per-pair; no C/E |
| spatialDM | Py | bivariate Moran's R per LR | per-LR significance only; no C/D/E/F |

Expected story: each recovers ALARMIST's motif LRIs *individually/fragmented* but none
assembles the program + direction + spatial ON region + downstream genes together.

## 3. Fair-comparison rules (shared inputs)
- Same cells, same cell-type labels, same LR resource (CellChatDB v2) across all methods.
- One frozen input per dataset in `results/comparators/_inputs/` (raw counts in X + `layers['counts']`,
  `obsm['spatial']`, matched `obs` cell types). Subset/downsample identically if a method OOMs on
  full Xenium — record the subset; run ALARMIST reference on the same subset.

## 4. Directory layout
```
results/comparators/
  README.md                # scorecard + per-method notes (this table, filled in)
  _inputs/{LUAD,GBM}.h5ad   # frozen shared inputs (+ provenance)
  alarmist/{LUAD,GBM}/      # pointer/copy of the pinned reference motif artifacts
  cellchat/{LUAD,GBM}/  liana/  niches/  commot/  cytosignal/  stlearn/  spatialdm/
  _figures/                # cross-method key plots
scripts/comparators/
  common/prepare_inputs.py # build the frozen inputs
  <method>/run.(py|R)      # one runner per method, reads _inputs/, writes results/comparators/<method>/
  compare/build_scorecard.py + head_to_head.py
```
Envs (conda; R methods can't share the Python env): `comp-r` (CellChat, NICHES, CytoSignal,
Seurat) + isolated Python envs per method (`comp-liana`, `comp-commot`, `comp-spatialdm`,
`comp-stlearn`) since stLearn/COMMOT pin conflicting deps. ALARMIST stays in `bptf`.

## 5. Deliverable — the key plots
1. **Capability scorecard**: methods × axes (A–H) check-matrix.
2. **Head-to-head on one motif** (per dataset): the reference motif's LRIs highlighted vs where
   each competitor ranks/places those same LRIs (fragmented), + spatial panel: ALARMIST motif-ON
   region vs the best competitor spatial output failing to isolate the same program.

## 6. Phasing
0 pin reference → 1 build frozen inputs → 2 Python methods (COMMOT, spatialDM, stLearn, LIANA+)
→ 3 R methods (CellChat, NICHES, CytoSignal) → 4 collate into scorecard + head-to-head figures.

---

# Runs log

## CytoSignal — GBM full run (2026-07-28) — the ALARMIST GBM reference dataset
- **Data:** `data/xenium_mm_final_cell_id.h5ad` — human Xenium GBM TMA (13 cores, global µm coords),
  100,197 cells × 5,119 genes. Already carries ALARMIST outputs (`motif` ×15, MES/AC/OPC/NPC-like,
  mGAM/non-mGAM, `grade`, `patch_id`) → **direct head-to-head substrate**. Raw counts in
  `layers['counts']` (X is log-norm); clusters = `cell_type` (9 types incl. mGAM/non-mGAM).
- **Full run (no subset):** 89,035 cells after `counts>=100` QC; fits 36GB. Human CellChat DB
  (858/1383 genes on panel), nosave, plots on (<200k gate), watchdog. Input `GBM/input_full/`,
  outputs `GBM/run_full/` (+ `quant/`). Exporter used `layers['counts']` + `cell_type`.
  TMA caveat: cores are spatially separated so diffusion (ε-ball 200µm) is clean; a few negligible
  inter-core Delaunay (contact) edges possible.
- **Run:** completed in ~16 min, 89,035 cells, peak swap ~3GB (comfortable on 36GB). 895 significant
  diffusion + 169 contact interactions. Outputs `GBM/run_full/{cluster_map.png, quant/}` (score
  matrices + reslists + summaries). `plotSignif` hit the custom-DB naming bug (skipped) → use the
  quant score matrices for per-LRI spatial maps.
- **★ Head-to-head result:** CytoSignal **recovers ALARMIST's mGAM↔MES-like LRIs as isolated, buried
  entries** — `GRN→SORT1` at **rank 66/895** (27,630 cells), `ANXA1→FPR1` at **rank 255/895** (15,292
  cells). Its top-ranked pairs are generic developmental signaling (7× WNT3-FZD, BMP2/6/8B, GABA),
  with ~50% of cells "significant" for the top pairs (non-discriminative). It gives **no grouping of
  GRN-SORT1 + ANXA1-FPR1 into one motif, no sender→receiver direction (mGAM→MES-like vs reverse), no
  grade association** — exactly axes C (co-occurrence) + D (direction) that ALARMIST provides. CytoSignal
  finds the pieces; ALARMIST finds the mGAM program. This is the LGG/GBM key-figure material.
- **Plots** (`plotSignif` hit the custom-DB naming bug → built directly from the quant score matrices):
  `GBM/run_full/plots/` = `spatial_scores.png` (per-LRI score maps over the 13 cores: GRN-SORT1,
  ANXA1-FPR1, WNT3, BMP2, DLL3-NOTCH1 + an mGAM/MES cell-type panel), `ranking.png` (top-25 all
  WNT/BMP; mGAM LRIs absent), `lr_panels_{GRN_SORT1,ANXA1_FPR1}.png` (ligand | receptor | score). Top-ranked grids
  (`plot_top_lris.py`): `top_diffusion_grid.png` (top 16) + `top_contact_grid.png` (top 12 — OCLN/
  cadherins/CADM3/JAM2/CD99/L1CAM/MPZ adhesion + DLL3→NOTCH1). All generic tissue signaling, none is the mGAM program.
- **★ Head-to-head figure** (`compare_mgam_alarmist.py` → `plots/comparison_mGAM_vs_cytosignal.png`),
  ALARMIST motif 1 (mGAM) vs CytoSignal on the same 89k cells. **Motif loadings come from
  `results/GBM/single_cell/cell_loadings.npy` (100197×20, `al.project_cell_loadings`), column 1 =
  mGAM motif** — verified: col 1 is mGAM's top-loading column, matching the notebook's `motif_idx=1`
  (NOT `obs['motif']`, which is a 15-category argmax label — my first draft used that and was wrong).
  **direction confirmed from expression** —
  GRN (ligand) & FPR1 (receptor) are mGAM-specific → mGAM **sends** GRN→SORT1 (to tumor) and **receives**
  ANXA1→FPR1 (from tumor): a bidirectional mGAM⇄MES-like loop = ALARMIST's one motif. CytoSignal gives
  the two LRIs as **separate undirected score maps ranked #66 & #255 of 895** (top 15 all generic WNT/BMP);
  both scores peak on mGAM (the hub) but are never grouped, directioned, or surfaced. Panels: motif-loading
  map | GRN→SORT1 map | ANXA1→FPR1 map | expression-direction heatmap | ranking | per-celltype score.
- **Reconstruction test** (`reconstruct_motif1.py` → `plots/reconstruct_motif1_from_cytosignal.png`):
  Σ (ALARMIST factor × CytoSignal per-cell LRI score) over motif 1's top-100 LRIs, weight = `factor`
  and `factor_lrnorm`, raw + per-LRI z-scored. **99/100 LRIs exist in CytoSignal** (the pairs are all
  there), yet even with ALARMIST's exact motif weights the reconstruction only reaches Pearson
  **r = 0.17 (factor raw) → 0.27 (factor_lrnorm z-scored)** vs the true motif-1 loading. Two reasons:
  CytoSignal's scores are neighborhood-smoothed/dense (diffuse vs the motif's punctate niche), and they
  are **cell-type-agnostic** (one score per L-R, collapsing ALARMIST's cell-type-pair resolution). So
  the LRIs are present but the motif is recoverable only WITH ALARMIST's weights — and even then only
  partially — reinforcing axes C/D.
- **Top-25 LRI grid** (`plot_motif1_top_lris_cytosignal.py` → `plots/motif1_top25_lris_cytosignal.png`):
  CytoSignal spatial score map for each of motif 1's top-25 LR pairs (by factor_lrnorm, unique L-R),
  each panel labeled with BOTH its ALARMIST motif-1 rank (+sender→receiver+mode+weight) AND its
  CytoSignal rank/mode/cells. Every motif-1-defining LRI is buried in CytoSignal (HLA-DQA1-CD4 ALARMIST
  #1 → CS #40/166; C4A-C3AR1 #2 → #387/895; ANXA1-FPR1 #5 → #255/895; GAS6-MERTK #6 → #202/895); some
  have 0 significant CytoSignal cells (THY1-ADGRE2, PTGES3-PTGER4). Only TGFB2-TGFBR2_TGFBR1 not scored.

### CytoSignal two-condition (high vs low grade) differential — the axis-G comparison (2026-07-29)
- **CytoSignal HAS a native workflow:** `mergeCytoSignal` (per-sample objects + dataset-level metadata)
  → `runNEBULA` (NB mixed model per interaction, core=subject, grade=fixed) → `plotNebulaVolcano`.
  Built the merged 13-core object (`nebula_grade/merged.rds`; 919 diff + 169 cont shared interactions;
  grade constant within core → 7 high, 6 low). `nebula` wouldn't compile in the conda R (4.3.3 / Eigen
  3.4); I first used interim pseudobulk tests (avg-score / fraction per core + MWU) — **now superseded
  and their scripts/CSVs removed** in favour of the real runNEBULA below.
- **ALARMIST test = the notebook's `plot_motif_fraction_by_grade`** (GMM-binarize loadings →
  motif-positive fraction per `tma_id` core → Mann-Whitney high vs low). Reproduced faithfully (scratchpad
  `alarmist_fraction_grade.py`: `al.gmm_binarize_all_motifs`, multi_sample by patient) →
  `alarmist_motif_fraction_grade.csv`. (An earlier draft wrongly used mean-loading+t-test — corrected;
  FDR is our own add-on, not in the notebook, so the figure shows raw p.)
- **runNEBULA (native) DOES run** — via **system R 4.4.2** (`nebula` 1.5.6 installs as a CRAN *binary*,
  no compile; the conda R 4.3.3 source-compile hit an Eigen 3.4 wall). cytosignal won't *link* in system
  R (Fortran/BLAS), so runNEBULA's exact inputs are extracted in conda R via `cytosignal:::.setup.model`
  (`nebula_inputs.rds`) and `nebula::nebula` is called in system R — this IS runNEBULA's computation
  (NB mixed model, core = random effect, `~grade`, offset = counts, cell-level). **Both stages are one
  script now:** `run_nebula_grade.R` runs stage 1 (build+merge+extract) in the conda R, then re-invokes
  itself in system R 4.4.2 (env `NEBULA_STAGE=2`) for the nebula step → `nebula_grade_results.csv`.
- **★ Corrected, nuanced result (2-panel `plots/grade_comparison_2panel.png`, `al.glm_volcano` style —
  runNEBULA volcano | ALARMIST motif-fraction volcano):** The native cell-level NB
  model has **more power than the pseudobulk MWU and DOES find grade-differential LRIs — 15 at BH<0.05**
  (I was wrong to imply CytoSignal finds nothing). BUT its top hits are **generic adhesion/junction/
  metabolite**: CDH2-CDH2 (↑), NECTIN3-NECTIN1 (↓), F11R/JAM3 junctions (↓), SLC-glutamate→GRM3 (↓),
  BMP5 (↑) — **not the mGAM program**. Of motif-1's top-100 LRIs, only **1** (JAM3-F11R) is BH-sig in
  runNEBULA (0 in avg-score MWU). ALARMIST flags the **mGAM motif itself** (`plot_motif_fraction_by_grade`,
  p=0.014, 17/20 motifs). So the fair statement: CytoSignal's per-LRI test (even the powerful native one)
  scatters grade signal across generic interactions and misses the mGAM program; ALARMIST's motif captures
  it. Scripts: `run_nebula_grade.R` (2-stage: conda-R build/merge → system-R nebula), `analyze_motif1_grade.R`,
  scratchpad `alarmist_fraction_grade.py`, `build_grade_2panel.py` (final figure).
  DB builder patched with empty `protein_name_a/b` (fixes getLigandNames /
  showIntr return.name / plotSignif / mergeCytoSignal on the custom DB).


## CytoSignal — first pass, P21 LUAD (2026-07-28)
- **Env:** isolated conda `comp-cytosignal` (R 4.3.3 + toolchain), CytoSignal 0.5.1 + SPARK +
  scattermore. Activate via `source scripts/comparators/cytosignal/activate_env.sh` (this box's
  conda is degraded — `conda activate/run` fall through to system R; see memory).
- **Input:** first-pass **2.5mm × 2.5mm crop** of `P21_LUAD_Xenium.h5ad` = 28,596 cells
  (immune/stroma/epithelial mix, only ~5% Tumor_epi), `annotation_coarse` (19 types), microns.
  Export + runner: `scripts/comparators/cytosignal/{export…,run_cytosignal.R,plot_cytosignal.R}`;
  outputs in `results/comparators/cytosignal/LUAD/`.
- **Params:** `scale.factor=1` (Xenium microns), `r.eps.real=200`, `counts.thresh=100`
  (default 300 drops ~half the panel). QC kept 25,206 cells.
- **Panel coverage:** only **649/977** CytoSignal DB LR-genes are on the Xenium panel (66%) →
  ~1/3 of the database is unmeasurable here (a per-method limitation for the scorecard).
- **Result:** 277 diffusion + 44 contact spatially-variable significant LRIs.
  Diffusion = TGFβ/PDGF/IL6/IFN/HGF (CAF+immune paracrine); contact = integrin-ICAM1, Notch
  (JAG/DLL→NOTCH), CD8/CD2/CD47 (juxtacrine immune). Output = **per-LR-pair lists**, no motif
  grouping, no directional cell-type program, no downstream gene impact — i.e. exactly the
  C/D/E/F axes CytoSignal can't cover vs ALARMIST.
- **Caveats / next:** crop is tumor-poor; redo on a tumor-rich window and/or full tissue before
  the head-to-head vs the ALARMIST LUAD motif. ~2.7GB rds per run (delete checkpoints).

### CytoSignal — CellChatDB v2 + quantitative outputs (2026-07-28, update)
- **DB swap for fair comparison:** replaced CytoSignal's bundled CellPhoneDB-v2 with alarmist's
  `data/LRdatabase/CellChatDBv2.0.human.csv` (same resource ALARMIST uses). Built via CytoSignal's
  own `formatLRDB` using gene symbols as the protein-ID space (identity g_to_u — no UniProt needed);
  `Cell-Cell Contact`→contact, else→diffusion. → 3,217 interactions / 1,384 LR genes (865 on panel).
  Builder `build_cellchat_db.R` → `cellchat_db_human.rds`; `run_cytosignal.R` takes it as 4th arg.
- **Quantitative outputs now persisted** (were missing — only names before): `run/quant/` holds
  `score_<slot>.mtx|.rds` (cells×interactions LR-score), `reslist_<slot>.rds` (significant cells per
  interaction), `signif_summary_<slot>.csv`. Exporter `quant_io.R`; `@score` is ~dense so large runs
  save it as compressed `.rds`.
- **CellChat-DB crop result:** diffusion top = TGFβ→TGFBR complexes, LAMA3→ITGA6_ITGB4, COL4A4→CD44
  (ECM now included); contact top = PECAM1, CDH5, JAM2/3, OCLN junctions.
- **Full-tissue P21 = hardware-blocked here (scale finding, axis H).** 560,183 annotated cells
  (`input_full/`), 498,422 after QC. CytoSignal got through findNN + imputeLR + started
  inferIntrScore (912/2,683 diffusion interactions valid) then blew into swap during null-scoring
  (it floors the permutation count at n_cells, so the null scales with cell count); watchdog killed
  it at ~57GB working set / 21GB swap. Full 498k needs a ≥64GB node — CytoSignal keeps single-cell
  resolution throughout, unlike ALARMIST's patch-aggregated counting that runs full P21 easily.
- **Deferred to a big node:** self-contained bundle at `results/comparators/cytosignal/bundle_bignode/`
  (scripts + env installer + corrected DB rds/CSV + portable exporter + README with resources/params).
  Run there: `Rscript run_cytosignal.R input_full run_full nosave db/cellchat_db_human.rds`.
- **DB freshness:** the CellChat CSV was re-exported/corrected mid-session (complex-table subunits,
  CXCL8 + NPR fixes); the DB rds was rebuilt from the current CSV. The stale-DB crop run was removed.
  Valid local result kept: `run_crop_2p5mm/` (bundled CellPhoneDB-v2 DB, 28k crop).
