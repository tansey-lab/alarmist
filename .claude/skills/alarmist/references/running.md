# Running ALARMIST (producing the substrate)

Read this when you need to *produce* an ALARMIST run, debug a stage, or decide `s`/`K`.
The analysis toolbox in `SKILL.md` assumes the run already exists. Dataset-specific
values (the actual `s`, `K`, cell-type column, LRI resource, sample keys, install path /
conda env) live in **CLAUDE.md** — read it first.

## Run it as a pipeline, not by hand

The workflow is standardized end to end — nothing here is optional or reorderable.
Assemble stage-by-stage in a notebook only to debug one stage. Entry points:

- **Nextflow (preferred for production):** `cd nextflow && nextflow run main.nf --input
  samplesheet.csv --outdir results -profile docker` (the `main.nf` lives under
  `nextflow/`, not the repo root; swap `-profile singularity` on HPC).
- **CLI:** `alarmist-patchify` → `alarmist-bptf` → `alarmist-project` → `alarmist-glm` →
  `alarmist-visualize`.
- **Python wrapper:** `run_pipeline.py` (a thin wrapper that chains the numbered
  `scripts/research/0*.py` stages with conda-env switching; defaults `--patch-size 50.0`,
  `--n-components 15`).

## Stages

1. **Patchify + LRI counting** — `al.PatchLRIAnalyzer.run_patchify` → patch × LRI count
   matrix. Expression is **binarized** at the cell level (`adata.X > 0`); each LRI count
   is a non-negative integer per patch, per signaling mode (autocrine / paracrine /
   juxtacrine).
2. **BPTF factorization** — `al.run_bptf` fits and returns a **fitted BPTF model object**
   (not the factors). `al.extract_factors(model)` → `(patch_loadings W, LRI factors V)`;
   `al.save_bptf_results` / `al.process_bptf_results` writes the full processed outputs to
   disk — raw and rescaled `W`/`V` plus the normalized **`lri_motifs` table** (`.npy` /
   `.csv`). Motifs = LRI factors `V` (what defines each motif) + patch loadings `W`.
3. **Single-cell projection** — `al.NeighborhoodLRIAnalyzer.run_neighborhood` builds a
   per-cell neighborhood-LRI matrix. Match the patch geometry **yourself**: construct it
   with `neighborhood_size=s` (same edge length as the patch `patch_size=s`) and pass
   `required_columns=<patch LRI columns>` so the column space aligns — the code enforces
   neither. Then `al.project_cell_loadings` (`V` **clamped** via the BPTF mode-1 clamp,
   cells chunked, `chunk_size=50000` default) → per-cell loadings `U`.
4. **GMM state calling** — `al.gmm_binarize_all_motifs` fits a two-component GMM on log
   `U` and writes `obs['motif_{k}_state']` (`'positive'`/`'negative'` = ON/OFF) and
   `obs['motif_{k}_loading']`; the positive-component posterior probability is returned as
   arrays (`return_arrays=True`), not written to `obs`. In multi-sample mode (dict input,
   or `multi_sample=True` + `sample_column`) the GMM is fit on **pooled** loadings so
   thresholds are consistent across samples (the single-sample default fits on one sample).
5. **Motif→gene impact (Poisson GLM)** — `al.run_poisson_glm_analysis`. Per motif × cell
   type × non-LR gene, a Poisson GLM of **raw counts** on the **z-scored log motif
   loading** (only strictly-positive loadings kept; z-scoring done within each cell type).
   Significance = Wald test + BH-FDR. **LRI genes are dropped inside this call** (via
   `extract_lri_genes`, circularity). A Spearman prefilter runs first to accelerate it —
   `al.spearman_prefilter_genes`, **on by default** (`prefilter_spearman=True`, p < 0.001).
   Separately, `al.compute_exclusion_mask` builds each cell type's marker set (rank-based
   Mann–Whitney U on `adata.X`, log2FC ≥ 1, p_adj ≤ 1e-5) so that **other** cell types'
   markers (transcript-bleed spillover) can be dropped when reading a cell type's impact.
   Visualize with `al.glm_volcano` / `al.glm_forest`.

## Input contract (get these wrong and the run is silently invalid)

- `obsm['spatial']` in **microns** (assumed, not validated — must be in the same units as
  `--patch-size`), a **cell-type column**, and **raw counts**. Patchify binarizes
  `adata.X` (presence/absence); the GLM reads counts via `--count-layer` (CLI default
  `layers['counts']`, core default `X`); marker detection reads `adata.X`. So put raw
  counts in **both** `X` and `layers['counts']`. **The core pipeline does not consume a
  normalized layer** — a normalized/log layer is needed only for the downstream `SKILL.md`
  analyses (Hallmark module scores via `sc.tl.score_genes`, pairwise DE, GSEA).
- **Coordinates overlap across samples.** Use multi-sample mode **only** when they
  actually overlap, and disambiguate by a per-sample key — the `sample_column` argument
  (CLI `--sample-column`; required when `--multi-sample`). Overlap is *allowed*, not
  validated: patches are built per sample and tagged by sample id, so a shared coordinate
  frame is fine. Multi-sample uses the **gene intersection** and a **unified LRI column
  space** (union of cell types) across samples.
- **Fix seeds:** BPTF has no `random_state` argument — reproducibility comes from
  `np.random.seed(random_state)` called immediately before `model.fit` (`run_bptf` does
  this). The GMM is seeded via its `random_state` (default 42).

## Counting model (why the counts look the way they do)

Each LRI feature is `f = (sender A, receiver B, ligand L, receptor R, mode m)`. Within a
patch, any ligand-expressing cell can signal to any receptor-expressing cell (juxtacrine
is **not** further restricted to adjacent cells — transcript bleeding blurs sub-patch
position). With `nL` type-A ligand+ cells and `nR` type-B receptor+ cells (multi-subunit
complexes use the **min count across subunits**; if any subunit is absent the count is 0):

- **Autocrine** (non-juxtacrine, A = B only): the number of type-A cells co-expressing all
  ligand + receptor genes (strict AND at the cell level) — a **co-expression cell count,
  not a product**.
- **Paracrine / juxtacrine:** the sender×receiver product `nL × nR`, minus the autocrine
  count when A = B (to avoid double-counting cells that express both).

Features with zero counts across all patches are dropped.

## Factor normalization (for interpretation)

Raw `V` mixes motif-specificity with global LR-pair prevalence. Divide by global LR-pair
prevalence (`V* = V / (µ_(L,R) + ε)`, ε = 1) so ubiquitous LRIs don't dominate motif
characterization (`add_normalized_scores`, called by `process_bptf_results`). The model is
scale-invariant, so loadings are rescaled `W̃ ∈ [0,1]` per motif (`c_k = max_i W_ik`,
`W̃ = W/c_k`, `Ṽ = V·c_k`); the combined score `Ṽ / (µ_(L,R) + ε)` reflects both
motif-specific enrichment and absolute magnitude. Always report a prevalence-normalized
version when ranking a motif's LRIs (these columns are already in the `lri_motifs` table).

## Choosing `s` and `K` (judgment — the value for a dataset goes in CLAUDE.md)

- **Patch size `s`.** Default **50 µm** (captures juxtacrine + paracrine). **Increase for
  sparse tissue** (e.g. lung) so each patch has enough cells to estimate reliably. Too
  large mixes distinct microenvironments; too small starves estimation.
- **Motif count `K`.** Analogous to the number of factors in NMF (`--n-components`,
  default 15). Generally robust; tune to balance program complexity against covering the
  biology present in the tissue. Provide a few extra components as slack for
  noise/redundancy when benchmarking.

## Known artifacts

- **Projection chunk-boundary artifact:** chunked `project_cell_loadings` can inject
  scale jumps at chunk edges (each chunk re-initializes its own BPTF projection
  independently). Check loading continuity across chunks; per-sample projection is the
  known workaround. (Diagnosed in `scripts/research/es_diag_chunks.py` and
  `cell_loadings_correlation_check.py`.)
- Motifs depend heavily on cell-type labels — a bad annotation produces meaningless
  motifs. Verify annotation upstream before trusting anything downstream.
