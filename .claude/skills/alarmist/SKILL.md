---
name: alarmist
description: >-
  The analysis playbook for interpreting ALARMIST results — reach for it whenever
  someone has a fitted ALARMIST run (motif loadings, motif ON/OFF states, or
  motif→gene impacts) and asks what to do next, even if they never say "ALARMIST".
  Covers: characterizing a motif (its defining LRIs, its cell-type senders/receivers
  and signaling direction, where it acts), screening all motifs to pick the
  biologically relevant ones, testing whether a motif is enriched in a condition
  (treatment vs control, high vs low grade, disease stage), extending the downstream
  impact (combinatorial motif-state programs, Hallmark module scores, GSEA), building
  a motif-derived gene signature for survival/clinical validation, subtyping cells
  into motif-usage niches, and checking whether a motif is reproducible across
  datasets / platforms / cohorts. Trigger on phrases like "what is this motif",
  "is this motif enriched in <group>", "which signaling programs change with
  <condition>", "impact / DEG genes of a motif", "motif signature", "motif niches",
  "is this motif real". Running ALARMIST to produce the outputs is a prerequisite
  owned elsewhere (see body); dataset-specific parameters (patch size, K, cohort,
  paths) live in CLAUDE.md.
---

# ALARMIST — analyzing motif results

ALARMIST discovers **microenvironment motifs**: recurring patterns of ligand–receptor
co-activity across tissue, found by BPTF on a patch × LRI matrix and projected to
single cells. **This skill is the analysis toolbox for motifs that have already been
discovered** — how to get from a fitted run to a biological answer. It does not run
the pipeline (that *produces* the substrate — see below) and holds no dataset facts:
patch size `s`, motif count `K`, cell-type column, LRI resource, which
samples/conditions, install path — all in CLAUDE.md. **Read CLAUDE.md first; don't
guess these.**

The shape of the work: **one substrate, many analyses.** A run yields a fixed set of
outputs; from there you *screen* the motifs, *select* the ones that bear on your
question, then *compose* independent tools. Pick-and-combine — don't run a fixed script.

## The substrate (what a run gives you)

Every tool below consumes these. `al.*` names are the reproducibility anchor and
belong here.

- **Motifs — LRI factors `V`**: which cell-type–LRI combinations define each motif.
  `al.extract_factors(model)` returns the raw `(W, V)` arrays. For interpretation you
  want the **prevalence-normalized** `V*` (so ubiquitous LRIs don't dominate) — that,
  plus the per-LRI ranking columns, lives in the **`lri_motifs` table** written by
  `al.save_bptf_results` / `al.process_bptf_results` (which call `add_normalized_scores`,
  `V* = V / (µ_(L,R) + 1)`). Rank a motif's LRIs off that table, not off `V` directly.
- **Patch loadings `W`**: motif activity per patch.
- **Cell loadings `U`**: motif activity per cell (`al.project_cell_loadings`, `V`
  clamped). The analysis workhorse.
- **Motif states**: `al.gmm_binarize_all_motifs` fits a two-component GMM on log `U`
  and writes per-cell `obs['motif_{k}_state']` (values `'positive'`/`'negative'`, i.e.
  ON/OFF) and `obs['motif_{k}_loading']`. The positive-component posterior probability
  is returned as arrays (`return_arrays=True`), **not** written to `obs`. In
  multi-sample mode (dict input, or `multi_sample=True` + `sample_column`) the GMM is
  fit on **pooled** cells so thresholds are consistent across samples.
- **Motif→gene impact**: per motif × cell type × non-LR gene, a Poisson GLM of raw
  counts on the z-scored log motif loading (`al.run_poisson_glm_analysis`). LR genes
  are dropped **inside that call** (circularity). Separately,
  `al.compute_exclusion_mask` builds each cell type's marker set so that, when reading
  one cell type's impact, **other** cell types' markers (transcript-bleed spillover)
  can be filtered out at the summary/plot stage.

**Motif indices are arbitrary across runs.** Two separate runs' "motif 3" are unrelated,
so any cross-condition comparison must come from a **single joint run** over all samples
(rule owned by `spatial-workflow`).

## Producing the substrate (prerequisite — not this skill's job)

Assumes segmentation → cell-type annotation is done and **checked**; motifs on bad
labels are meaningless (`spatial-workflow` gate). Run via the nf-core Nextflow pipeline
(`cd nextflow && nextflow run main.nf --input samplesheet.csv --outdir results -profile docker`),
the CLI (`alarmist-patchify` → `alarmist-bptf` → `alarmist-project` → `alarmist-glm` →
`alarmist-visualize`), or `run_pipeline.py`. Two things silently invalidate a run: raw
counts must be present in **both** `X` and `layers['counts']` (patchify binarizes `X`;
the GLM reads `layers['counts']`; marker DE reads `X`) — the core pipeline does **not**
need a normalized layer (that is only for the downstream module-score / GSEA analyses
below); and multi-sample mode needs coordinates that actually overlap, disambiguated by
a per-sample key (`sample_column`, CLAUDE.md). Fix seeds: `np.random.seed` before BPTF
(BPTF has no `random_state` argument — reproducibility rides on this), `random_state` on
the GMM. Full stage-by-stage detail, the counting model, and `s`/`K` selection judgment →
`references/running.md`.

## Toolbox (compose per question)

Each tool is independent. `al.*` = packaged function. *(pattern)* = a reference
implementation lives in the repo's `scripts/research/`; **adapt it, it is not an API.**
Start every analysis by screening, then selecting.

### Screen all motifs, then select
Before analyzing one motif, get the lay of all `K`: rank motifs by overall activity
(`al.get_top_motifs` on patch loadings `W` — this ranks *motifs*, not LRIs), read the
cell-type-pair × motif factor summary as a clustermap
(`alarmist.plotting.plot_lri_clustermap`), skim each motif's network and top LRIs, and
pick the motifs that bear on your biological question. Every paper case study starts
here (e.g. spotting the two motifs that captured contrasting vascular states).
`al.plot_lri_networks` across motifs; top-LRI dotplots (`al.plot_top_lri_interactions_dot`);
a single-cell-type-level factor summary *(pattern)*.

### Characterize a motif — "what is this motif?"
- **Defining LR pairs** → `al.plot_top_lri_interactions_dot`,
  `alarmist.plotting.plot_top_lri_interactions_by_pathway`. Both read the `lri_motifs`
  table's factor / `V*` / score columns (run `annotate_pathways` first for the pathway
  variant).
- **Who signals to whom** → `alarmist.plotting.plot_celltype_communication_by_motif`.
  This is **directional** — the heatmap has rows = sender, columns = receiver, and is
  not symmetrized, so read sender→receiver and look for asymmetry (a cell type
  predominantly on the ligand vs receptor side). The mGAM story turned on splitting one
  motif into mGAM→MES-like vs MES-like→mGAM and finding direction-specific LRIs
  (GRN→SORT1 one way, ANXA1→FPR1 the other). *(directional split: pattern)*
- **Cell-type composition of ON cells** → `al.weighted_celltypes_by_motif` +
  `alarmist.plotting.plot_motif_celltype_composition`.
- **Network topology** → `al.plot_lri_networks` / `plot_lri_networks_html`.
- **Two motifs sharing a cell-type pair — shared vs distinct LRIs** → *(pattern)*.
  (How the paper separated an mGAM motif from its near-identical non-mGAM twin.)

### Locate / count motif states
- **Spatial map** → `alarmist.plotting.plot_motif_spatial`. **One panel per sample** —
  pass a dict of AnnData or a `sample_column` so it grids one panel per sample; a single
  merged AnnData with `sample_column=None` superimposes overlapping tissues onto one axes.
- **ON/OFF counts** → `al.compute_motif_state_counts` /
  `alarmist.plotting.plot_motif_state_counts` (state labels are `'positive'`/`'negative'`);
  motifs-per-cell → `al.compute_positive_motifs_per_cell`.
- **Where two motifs co-occur** (both-ON vs either-ON regions) is the entry point to the
  combinatorial programs below. *(pattern)*

### Is a motif enriched in a condition? — compositional, not DE
Compute the **motif-positive fraction per sample**, then test across samples
(Mann–Whitney; paired / mixed-effects when samples match by patient). **The sample is
the replicate — never test per cell** (scoping owned by `spatial-workflow`). Works for
any two-group contrast: treatment vs control, high vs low grade, progression stages.
*(pattern)*

### Extend the motif→gene impact (GLM already fit)
- **Summarize / visualize** → `al.analyze_glm_results` builds/loads the marker-exclusion
  mask and drives the volcano + forest plots (`al.glm_volcano` / `al.glm_forest`); it
  does **not** itself tabulate up- vs down-regulated DEG counts — compute that tally from
  the per-`(motif, cell type)` result tables.
- **Combinatorial motif-state programs.** Define mutually-exclusive groups (A-only /
  B-only / both) and, within each cell type, compare Hallmark **module scores**
  (`sc.tl.score_genes` on normalized log expression) across groups by Mann–Whitney —
  the both-ON group often reads as an **intermediate** state. Follow with pairwise DE
  (per-sample normalized expression, BH-FDR) and **GSEA** (`gseapy` prerank on the
  logFC-ranked list vs Hallmark). DE/GSEA mechanics owned by `spatial-workflow`.
  *(pattern)*

### Motif signature → survival / clinical validation
Turn a motif's impact into a prognostic score: take the top-N genes by logFC among ON
cells of the relevant cell type (all passing the GLM FDR), score them in an **external
bulk cohort**, split patients high/low, test by **Kaplan–Meier + log-rank**. Moves a
motif from descriptive to prognostic. **There is no packaged helper and no reference
script in the repo for this** — implement it yourself (score the signature in the
cohort, stratify, then `lifelines` / `scikit-survival` for KM + log-rank). The cohort is
a dataset fact → CLAUDE.md.

### Motif-defined niches — subtype cells by motif usage
**k-means** on the **per-motif max-scaled** cell loadings `U` (never raw; fixed seed),
K by elbow → niches that are *motif-usage* profiles. These are a **different niche type**
from neighbourhood-composition niches (which `spatial-niche` owns) — but they are still
niches, so its rules apply: **annotate them, never report bare cluster IDs.** To compare
niches across conditions, fix K on one side and apply it to both, then match by Pearson
correlation of per-grid-square proportion vectors. K + grid size → CLAUDE.md.
To relate motifs to composition-niches (are motif-ON cells concentrated in certain
niches?), partial-correlate motif loadings/states against niche membership, controlling
for a confounder. *(pattern — `scripts/research/*_niche_clustering.py`,
`*_niche_concordance.py`)*

### Cross-dataset reproducibility — "is this motif real or an artifact?"
Generalizes to any two comparable datasets (platforms, cohorts, batches).
- **Match motifs:** cosine similarity of LRI factors `V` (from the `lri_motifs` table),
  **restricted to shared cell-type–LRI combinations**; call a motif *well-aligned* if
  its best cross-dataset match exceeds a threshold (paper: 0.4). A complementary check
  aligns on the **impact** vectors (per-gene × cell-type logFC) instead. *(pattern)*
- **Diagnose a mismatch** at three levels — LR genes, LR pairs, cell-type–LRI — the gap
  usually appears only at cell-type–LRI (i.e. *which* cell type carries the LRI).
  *(pattern)*
- **Spatial + niche + impact concordance:** Lee's L on motif-positive fraction per bin
  (`spatial-stats`); niche matching (above); bin genes by cross-dataset expression
  correlation and correlate matched-motif logFCs within bins. *(pattern)*
- **Confounder first:** gene-**panel overlap** bounds achievable agreement — quantify it
  before reading disagreement as biology. If the sections aren't co-registered (serial
  sections don't share a frame), align first via the landmark similarity transform
  (`spatial-stats`). *(pattern)*

### QC on outputs
- Chunked `project_cell_loadings` can inject scale jumps at chunk edges (each chunk
  re-inits its own BPTF projection) — check loading continuity across chunks; per-sample
  projection is the workaround. *(pattern — `scripts/research/es_diag_chunks.py`,
  `cell_loadings_correlation_check.py`)*
- Motifs depend heavily on cell-type labels — verify annotation before trusting them
  (`spatial-workflow` gate).

## Worked chains
- *Which programs change with treatment / progression?* → enrichment test
  (motif-positive fraction) → read the GLM impact for the motifs that moved → GSEA on
  the ranked genes.
- *What is this motif and where does it act?* → top-LRI + directional cell-type
  communication → spatial map per sample.
- *Real signal or platform artifact?* → factor match → panel-overlap control → impact
  concordance.
- *Does this motif carry prognostic weight?* → characterize → enrichment in the
  aggressive group → signature → survival.

## Delegations & pipeline position
ALARMIST is `spatial-workflow` step 5. Owned elsewhere, referenced above:
- `spatial-workflow` — pipeline order, the replicate-unit rule, the joint-run rule,
  DE/GSEA mechanics, the annotation gate.
- `spatial-niche` — neighbourhood-composition niches + niche annotation rules.
- `spatial-stats` — Moran's I, Lee's L, landmark similarity transform.

**Interpretation guardrail:** ALARMIST is hypothesis-generating — motif / impact / LRI
findings assume mRNA tracks protein, so state them as hypotheses and validate key LRIs
(a specific ligand→receptor claim) at the protein level before calling them mechanism.
