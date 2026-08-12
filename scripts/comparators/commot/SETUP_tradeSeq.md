# Enabling `ct.tl.communication_deg_detection` — the R/tradeSeq half of COMMOT

**Nothing in this file has been executed.** It is a plan written from the installed source, for
you to run. Every claim about what the code does was read out of
`/Users/jiayifan/anaconda3/envs/comp-commot/lib/python3.10/site-packages/commot/tools/_downstream_analysis.py`
(lines 32–183) — *not* from the docs, which are wrong in at least one place (see Trap 1).

---

## 0. Why this is needed at all

COMMOT is a pure-Python package. `pip show commot` lists **no R dependency**:

```
anndata, importlib-metadata, karateclub, leidenalg, networkx, numpy,
pandas, plotly, pot, pysal, python-igraph, python-louvain, scanpy, scikit-learn
```

Exactly **one** of its functions — `communication_deg_detection` — imports `rpy2` + `anndata2ri`
*inside its own body* and drives R's **tradeSeq** and **clusterExperiment**. The imports are
function-local by design, which is why `import commot` and everything else works without R and
why this dependency is invisible until you call that one function:

```python
>>> ct.tl.communication_deg_detection(None)
ModuleNotFoundError: No module named 'rpy2'
```

It is the **discovery** step of the downstream chain — it nominates the genes whose expression
tracks the amount of signal a cell receives. Without it, `communication_deg_clustering`,
`plot_communication_dependent_genes` and the tutorial's `ds_genes` for `communication_impact` all
have no input.

---

## 1. What to install

Build a **separate env**. Do not touch `comp-cellchat`, `comp-niches` or `comp-cytosignal` — a
shared R would put their pinned Seurat/CellChat stacks at risk.

```bash
CONDA_SOLVER=classic /Users/jiayifan/anaconda3/bin/conda create -n comp-commot-r \
  -c conda-forge --override-channels -y \
  python=3.10 'numpy<2' r-base=4.3 rpy2 anndata2ri
```

Then, inside that env's R, the two Bioconductor packages:

```bash
/Users/jiayifan/anaconda3/envs/comp-commot-r/bin/R --no-save -e \
  'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("tradeSeq","clusterExperiment"), ask=FALSE, update=FALSE)'
```

Finally COMMOT itself and the Python side:

```bash
/Users/jiayifan/anaconda3/envs/comp-commot-r/bin/pip install commot==0.0.3 scanpy anndata
```

Freeze it when it works:

```bash
CONDA_SOLVER=classic /Users/jiayifan/anaconda3/bin/conda env export -n comp-commot-r --no-builds \
  > scripts/comparators/commot/env.r.lock.yml
```

### On the version numbers

The docstring says *"tradeSeq version 1.0.1 with R version 3.6.3 has been tested to work … rpy2==3.4.2
and anndata2ri==1.0.6 have been tested to work"*, and a source comment adds
`# !!! anndata2ri works only with 3.6.3 on the tested machine`.

Those are the versions **the authors happened to test in 2021**, not declared minimums. R 3.6.3 is
long unsupported and current tradeSeq targets R ≥ 4.3, so the recipe above deliberately uses a
modern stack. **This combination has not been verified by anyone here.** Section 3 lists the exact
places it can break and what each failure looks like, so you can tell a version gap from a data
problem. If modern versions do fail, the fallback is the authors' pinned set in an isolated env —
but try modern first; do not start by building a 2020 stack.

---

## 2. Verify the install before touching real data

```bash
/Users/jiayifan/anaconda3/envs/comp-commot-r/bin/python - <<'PY'
import rpy2, anndata2ri, rpy2.robjects as ro
print("rpy2", rpy2.__version__, "| anndata2ri", anndata2ri.__version__)
ro.r('library(tradeSeq); library(clusterExperiment)')
print("R:", ro.r('R.version.string')[0])
print("tradeSeq:", ro.r('as.character(packageVersion("tradeSeq"))')[0])
print("clusterExperiment:", ro.r('as.character(packageVersion("clusterExperiment"))')[0])
# the three attribute accesses COMMOT makes -- these are the likely rpy2-version breakages
ro.numpy2ri.activate(); ro.pandas2ri.activate(); anndata2ri.activate()
print("activate() calls OK")
anndata2ri.deactivate(); ro.numpy2ri.deactivate(); ro.pandas2ri.deactivate()
PY
```

If every line prints, the plumbing is sound and the rest is data.

---

## 3. The five places this can break, and how each looks

Read against the installed source; these are not hypothetical.

| # | Where | Symptom | What it means |
|---|---|---|---|
| 1 | `ro.numpy2ri.activate()` (`:76`) | `AttributeError: module 'rpy2.robjects' has no attribute 'numpy2ri'` | rpy2 ≥ 3.5 made `numpy2ri`/`pandas2ri` submodules that must be imported explicitly. COMMOT relies on something else having imported them (`anndata2ri` usually does). Fix on **our** side, not by editing the package: `from rpy2.robjects import numpy2ri, pandas2ri` before calling COMMOT. `run_commot_deg.py` already does this. |
| 2 | `anndata2ri.py2rpy(adata_deg)` (`:143`) | converter / `SingleCellExperiment` errors | anndata2ri's converter API moved between 1.0 and 1.1+. If it fails, pin `anndata2ri==1.0.6`. |
| 3 | `fitGAM(counts=X, pseudotime=…, cellWeights=…, nknots=6, verbose=TRUE)` (`:153`) | unused-argument or signature error in R | tradeSeq's `fitGAM` gained a matrix/SCE dispatch. The positional-name form used here has been stable, but check `?fitGAM` if it errors. |
| 4 | `associationTest(sce, global=FALSE, lineage=TRUE)` (`:155`) | later `KeyError: 'waldStat_1'` in Python | COMMOT hard-codes the column names `waldStat_1`, `df_1`, `pvalue_1`, which only exist when `lineage=TRUE`. A tradeSeq version that renames them breaks the Python side, not the R side — the error surfaces one step later. |
| 5 | `clusterExpressionPatterns(sce, nPoints=50, genes=…, k0s=4:5, alphas=c(0.1))` (`:166`) | RSEC / clusterExperiment errors | `k0s` and `alphas` are forwarded to clusterExperiment's RSEC. This is the least stable call in the chain. |

**None of these should be fixed by editing the installed package.** Record the gap, work around it
in our runner, and write the deviation into `DEVIATIONS.md` — same rule that produced the
`legendHandles` entry there.

---

## 4. Traps in the data contract (all verified from source)

**Trap 1 — the docstring is wrong about the layer name.** It says the counts "should be available
through `adata.layers['count']`" (singular). The code at `:139` reads **`adata.layers['counts']`**
(plural). Follow the code.

**Trap 2 — our saved AnnData has neither `layers['counts']` nor `.raw`.** `run_commot.py` clears
`layers` before the run and drops `.raw` before writing `adata_commot.h5ad`, to halve the file. So
both must be rebuilt from `data/xenium_mm_final_cell_id.h5ad`'s `layers['counts']` — verified to be
genuine integer counts (max 128). `run_commot_deg.py` does this. **It is still not an OT re-run** —
the transport plans in `obsp` are read straight off disk.

**Trap 3 — gene selection is wide open by default.** With `n_var_genes=None` and `var_genes=None`,
the function calls `highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5)` — Seurat's
defaults, tuned for whole-transcriptome data. On a 5,119-gene targeted Xenium panel the number of
genes that pass is unpredictable, and it directly sets how long `fitGAM` runs. **Pass
`--n-var-genes` explicitly** and record it.

**Trap 4 — R gets a DENSE matrix.** `:146` runs `X <- as.matrix( assay(adata,'X') )`. For core 1
(26,456 cells) at 2,000 genes that is ~420 MB in R on top of everything else; the whole slide is not
an option (and COMMOT is per-core here anyway, for the same dense-matrix reason as the main run).

**Trap 5 — `pathway_name` overrides `lr_pair`.** If `pathway_name` is given, `lr_pair` is ignored
(`:145-149`). For `GRN` the two are equivalent anyway: the GRN pathway contains exactly one LR pair,
`GRN→SORT1`, so `pathway_name='GRN'` and `lr_pair=('GRN','SORT1')` select the same column.

---

## 5. Runtime — plan for this, it is the real cost

The tutorial's `fitGAM` took **4 min 13 s on 3,355 Visium spots**. Our cores run 819 → 26,456
cells, 100,197 in total ≈ **30× the tutorial's cell count**, and `clusterExpressionPatterns` runs
on top of that. A whole-cohort estimate in the multi-hour to overnight range is realistic, and it
scales with the gene count you choose in Trap 3.

**Do not start with all 13 cores.** The order that finds problems cheapest:

1. `--cores 2` (819 cells) with `--n-var-genes 500` — proves the R plumbing end to end in minutes.
2. `--cores 13` (9,126 cells, GRN is top-5 there) at the gene count you actually want — gives a
   real per-core timing to extrapolate from.
3. Then the rest, in background, smallest-first (the runner already sorts that way).

`--skip-deg` reuses a previous `deg_<pathway>_core<NN>.pkl`, so the expensive step is paid once
even if a later step fails. Pickle-first is the tutorial's own pattern.

---

## 6. What to run

```bash
PY=/Users/jiayifan/anaconda3/envs/comp-commot-r/bin/python

# 1) smoke: smallest core, few genes
$PY scripts/comparators/commot/run_commot_deg.py --cores 2 --pathway GRN --n-var-genes 500

# 2) timing: a real core
$PY scripts/comparators/commot/run_commot_deg.py --cores 13 --pathway GRN --n-var-genes 2000

# 3) the cohort, background
nohup $PY scripts/comparators/commot/run_commot_deg.py --pathway GRN --n-var-genes 2000 \
  > logs/commot_deg_GRN.log 2>&1 &
```

`--pathway FGF` is worth a second pass: FGF is top-5 in **13/13** cores and has 14 LR pairs, so
`plot_communication_impact` clusters properly there, whereas GRN's single pair needs
`cluster_knn=1` (the runner falls back automatically).

---

## 7. What you get, and the one comparison this unlocks

Per core, `results/comparators/commot/GBM/deg/<pathway>/`:

| File | Meaning |
|---|---|
| `deg_core<NN>.pkl` | `{df_deg, df_yhat}`, the tutorial's own checkpoint |
| `df_deg_core<NN>.csv` | one row per tested gene: `waldStat`, `df`, `pvalue`, sorted by Wald |
| `df_yhat_core<NN>.csv` | smoothed expression, genes × 50 points along the received-signal axis |
| `df_deg_clustered_core<NN>.csv`, `df_yhat_clustered_core<NN>.csv` | after `communication_deg_clustering` |
| `top_de_genes_core<NN>.txt` | what `plot_communication_dependent_genes(return_genes=True)` returns — **this is the list the whole chain exists to produce** |
| `heatmap_deg_core<NN>.png` | expression trends vs increasing received signal |
| `impact_core<NN>.csv` + `_{sender,receiver}.png` | `communication_impact` on those genes, tutorial arguments |
| `examples_core<NN>.png` | the tutorial's 3-panel figure: received signal, one negative DE gene, one positive |
| `run_manifest.json` | versions, parameters, per-step timings, failures |

**The point of doing this.** The benchmark is currently missing its most natural cross-method
comparison: *what genes does COMMOT itself nominate downstream of GRN→SORT1, and how much do they
overlap ALARMIST's motif-1 impact genes?* Only `communication_deg_detection` can answer the first
half. The earlier substitute run (`results/comparators/commot/GBM/impact/`, marked partial) fed
ALARMIST's gene list **into** COMMOT, which tests something narrower and cannot produce an overlap
statistic.

Once `top_de_genes_core<NN>.txt` exists, the comparison is:

```
COMMOT's nominated genes   ∩   results/GBM/impact/motif_1_celltype_mGAM_de_results.csv (qval<0.05)
```

against a background of the 5,119-gene panel, with a hypergeometric test. **That is the number the
benchmark actually wants.** Two controls to run alongside it, both cheap once the lists exist:
the same overlap for a **different ALARMIST motif's** mGAM genes (motif 14 — same method, same cell
type, only 6/30 gene overlap with motif 1), and for a **different pathway's** COMMOT gene list
(FGF), so "these overlap" can be separated from "any two myeloid-flavoured gene lists overlap".

## 8. Interpreting the numbers

- `df_deg.pvalue` is tradeSeq's association-test p-value — **not FDR-corrected by COMMOT**. Apply
  BH yourself before calling anything significant.
- `communication_impact` scores are **not correlations**. `_utils/_similarity.py:85` returns the
  signal feature's importance percentile among the ~500 background genes, averaged over 100 refits,
  so a useless feature scores ≈ **0.5**, not 0. Read every impact number against 0.5. This is why
  the tutorial's own `df_impact_PSAP` table sits entirely in 0.78–0.999.
- The core is the replicate unit (7 high-grade / 6 low-grade). Per-cell statistics inside one core
  are descriptive; anything about grade is a core-level test.
