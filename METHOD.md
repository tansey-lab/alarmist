# ALARMIST — Patch-Level LRI Quantification: Method and Proposed Extension

This document has two parts:

1. **How the current method works** — how ALARMIST quantifies ligand–receptor
   interactions (LRIs) inside a 50 µm spatial patch, why it does **not** physically
   separate paracrine from juxtacrine signaling, and the matrix trick that makes the
   computation fast.
2. **A proposed extension** — how to quantify paracrine and juxtacrine LRIs
   *separately* while still folding both into the **same** patch feature vector,
   *without* sub-dividing the patch into smaller patches, and while keeping the
   computation efficient.

The reference implementation is `PatchLRIAnalyzer` in
[`src/alarmist/core/lri.py`](src/alarmist/core/lri.py) (`create_spatial_patches`,
`create_column_structure`, `build_patch_lri_matrix`).

---

## Part 1 — Current method

### 1.1 Overview

ALARMIST turns a spatial transcriptomics sample into a **patch × LRI-channel count
matrix**. Each row is a square tissue patch; each column is one
`(sender_cell_type, receiver_cell_type, ligand, receptor, mode)` channel. The entry
is an integer count of how many sender–receiver cell pairs in that patch could
support that interaction. This matrix is the input to the downstream Bayesian tensor
/ matrix factorization that discovers signaling *motifs*.

### 1.2 Step 1 — Spatial patches (`create_spatial_patches`)

A regular axis-aligned grid is laid over the tissue at `patch_size = 50 µm`
(default). Cell coordinates come from `adata.obsm["spatial"]`; each cell is binned
into a patch by `np.digitize` on the x and y grid edges, and given a single integer
`patch_id = x_index * n_y_bins + y_index`. Patches are **fixed, non-overlapping,
and independent of cell type or expression** — purely a geometric tiling. Cells that
fall in the same 50 µm × 50 µm tile share a patch.

### 1.3 Step 2 — LR database (`prepare_lri_database`)

Ligand–receptor pairs are read from CellChatDB (default) or CellPhoneDB. A pair is
kept only if **all** of its ligand genes and **all** of its receptor genes are
present in the panel (complexes are `_`/`,`-separated and must be fully present).
Each pair carries a `signaling_type` annotation from the database — importantly, the
value `"Cell-Cell Contact"` is what later marks a pair as "juxtacrine".

### 1.4 Step 3 — Column structure (`create_column_structure`)

For every retained LR pair and every ordered cell-type pair `(A, B)` a set of
columns is emitted. The **mode** of each column is decided purely from the
database annotation:

| Database `signaling_type` | Sender/receiver | Columns emitted |
| --- | --- | --- |
| `Cell-Cell Contact`       | any `(A, B)`   | one `juxtacrine` column |
| anything else (secreted)  | `A == B`       | `autocrine` **and** `paracrine` columns |
| anything else (secreted)  | `A != B`       | one `paracrine` column |

So a column name looks like
`A|B|ligand|receptor|juxtacrine` (or `…|paracrine`, `…|autocrine`).

**Key point:** whether a channel is called *juxtacrine* or *paracrine* is a
**database label**, not a measurement of physical distance. It is fixed before any
coordinates are looked at.

### 1.5 Step 4 — The count matrix (`build_patch_lri_matrix`)

This is where the "matrix computation to simplify/accelerate" happens.

**Binarization.** Expression is thresholded to presence/absence:
`expr_bool = (X > 0)`. Interaction potential is therefore counted in **numbers of
cells**, not expression magnitude.

**Notation.** For a patch `p`, an LR pair with ligand-complex gene set `𝔏` and
receptor-complex gene set `ℜ`, and a directed cell-type pair `A → B`, define

- `s_p` = number of **ligand-positive senders** in `p`
  = #{ cells `i` in `p` : `type(i) = A` and `i` expresses all genes in `𝔏` }
- `r_p` = number of **receptor-positive receivers** in `p`
  = #{ cells `j` in `p` : `type(j) = B` and `j` expresses all genes in `ℜ` }
- `a_p` = number of **self / autocrine** cells in `p` (only when `A == B`)
  = #{ cells in `p` of type `A` that co-express **both** ligand and receptor }

**The all-to-all identity.** The number of *ordered* sender→receiver cell pairs that
co-occur in patch `p` is simply the product of the two marginal counts:

```
number of (sender, receiver) pairs in patch p  =  s_p · r_p
```

This is the crux of the acceleration: you never enumerate the O(n²) cell pairs inside
a patch. Because "every ligand-positive sender can be paired with every
receptor-positive receiver", the pair count collapses to a **product of two
per-patch marginal counts**.

**How the marginals are computed (vectorized).** For each cell type and each
individual ligand/receptor gene, the code builds a per-patch count vector by a sparse
**cell → patch** aggregation (`patch_by_lig`, `patch_by_rec`):

- restrict the binarized expression to cells of type `A` (or `B`),
- for gene `g`, scatter-add the expressing cells into their patches
  (a `coo_matrix` with `sum_duplicates`), giving a length-`n_patches` vector.

A `patch_by_cell` matrix (`n_patches × n_cells`) aggregates any per-cell quantity to
patches with a single sparse mat-vec.

**Values written per mode.** With `s_p`, `r_p`, `a_p` in hand:

```
juxtacrine (A→B):   s_p · r_p − a_p        # a_p ≠ 0 only when A == B
autocrine  (A==B):  a_p
paracrine  (A→B):   s_p · r_p − a_p   if A == B
                    s_p · r_p        if A != B
```

The autocrine term `a_p` (computed from true per-cell co-expression via
`patch_by_cell.dot(coexpr)`) is subtracted so a single cell that happens to be both
ligand- and receptor-positive does not count as signaling to itself.

**Complexes (implementation detail).** For the marginal counts `s_p`, `r_p` a
multi-gene complex is handled as an element-wise `min` of the individual per-gene
patch counts. This is an *upper-bound approximation* of true co-expression at the
patch level (it can exceed the number of cells expressing *all* subunits), and is
exact for single-gene ligands/receptors (the large majority of pairs). The autocrine
term `a_p` uses exact per-cell co-expression.

### 1.6 Why paracrine and juxtacrine are *not* physically distinguished

Notice that the juxtacrine and paracrine formulas are **identical** (`s_p·r_p`, minus
self-pairs when `A == B`). The only thing that differs is the *label*, taken from the
database `signaling_type`. Concretely:

- A "juxtacrine" channel is **not** restricted to physically touching cells. It counts
  every ligand-positive sender against every receptor-positive receiver anywhere in
  the same 50 µm patch — two cells ~70 µm apart across the patch diagonal are counted
  exactly like two adjacent cells.
- The reason the method *cannot* recover contact information is mathematical: the pair
  count `s_p · r_p` is a **product of marginals**. It records *how many* senders and
  *how many* receivers are present, but discards *which sender is next to which
  receiver*. Adjacency lives in the joint distribution, which the marginal product
  throws away.

So within a patch, ALARMIST currently measures **co-presence within 50 µm**, and uses
the database annotation as a proxy for "this pair is contact-dependent". It does not
measure contact itself.

### 1.7 Efficiency of the current method

Let `P` = number of patches, `C` = number of columns
(`≈ n_LR_pairs · n_celltypes² · {1 or 2}`), `nnz` = nonzeros in the expression matrix.

- Building marginals: one pass over `nnz`, plus sparse scatter to patches.
- Main loop: for each of `C` columns, a handful of length-`P` element-wise ops → `O(C·P)`.
- **No** `O(n_cells²)` pairwise enumeration anywhere.

This is what makes 50 µm patches over whole slides tractable. The trade-off is the
loss of intra-patch geometry described in §1.6.

---

## Part 2 — Proposed extension: separate juxtacrine and paracrine *inside* one patch

**Goal.** Keep the 50 µm patch as the unit of the feature vector, but split each
patch's signal into a **contact (juxtacrine)** part and a **diffusion (paracrine)**
part — as *separate columns on the same patch row* — without cutting smaller patches
and without an O(n²) blow-up.

### 2.1 The missing ingredient: a cell-level contact graph

The only thing the current pipeline lacks is *adjacency*. Add it once, globally, as a
sparse symmetric neighbor graph over cells:

```
C ∈ {0,1}^(n_cells × n_cells),   C[i,j] = 1  iff cells i and j physically touch,
                                  C[i,i] = 0  (no self-loops)
```

Build options (all `O(n log n)`; `NeighborhoodLRIAnalyzer` already builds a KD-tree
graph in this codebase, so the machinery exists):

- **Delaunay triangulation**, edges pruned above a max length (~30 µm) — parameter-free,
  adapts to local cell density, good default for segmented single-cell data.
- **Radius graph** at a physical contact radius `r_c` (e.g. 10–30 µm) via
  `sklearn`/`scipy` KD-tree `query_radius`.
- **kNN graph** (e.g. k = 6) — simplest, fixed degree.

Average degree is small (~6 for Delaunay), so `nnz(C) ≈ 6·n_cells` — cheap to store
and multiply.

To keep the within-patch bookkeeping exact (see §2.3), restrict the graph to
**intra-patch edges**:

```
C̃[i,j] = C[i,j] · 1[ patch(i) == patch(j) ]
```

This mask is a single vectorized operation on the edge list.

### 2.2 Contact (juxtacrine) count per patch — vectorized

Reuse the exact same sparse machinery as today, inserting `C̃`. For a directed pair
`A → B` and an LR pair, define per-cell indicator vectors:

```
l ∈ {0,1}^n :  l_i = 1  iff type(i)=A and i expresses all ligand genes
r ∈ {0,1}^n :  r_j = 1  iff type(j)=B and j expresses all receptor genes
P ∈ {0,1}^(n_patch × n_cell) :  cell → patch assignment (already used today)
```

Then the number of **adjacent** sender→receiver pairs, aggregated to the sender's
patch, is three cheap operations:

```
t = C̃ · r                 # per sender cell i: # of touching receptor-positive receivers
J = P · ( l ⊙ t )          # length-n_patch juxtacrine count  (⊙ = element-wise)
```

`C̃` has a zero diagonal, so self-pairs are excluded automatically — no autocrine
leakage. Each ordered contact pair is counted once, at the (shared) patch of its two
endpoints.

### 2.3 Paracrine (diffusion) count = the remainder — exact, non-negative

The within-patch ordered pairs partition **exactly** into three disjoint groups:

```
s_p · r_p   =   a_p   +   J_p   +   Para_p
(all pairs)  (self)  (contact)  (non-contact, same patch)
```

So paracrine needs **no new pass** — it is what is left over after removing self and
contact pairs:

```
Para_p = s_p · r_p − a_p − J_p     (A == B)
Para_p = s_p · r_p − J_p           (A != B)
```

Because `J` uses the *intra-patch* graph `C̃`, contact pairs are a strict subset of
the `s_p·r_p` within-patch pairs, so `Para_p ≥ 0` is guaranteed and the three
channels sum back to today's all-to-all value. `s_p·r_p` and `a_p` are computed
**exactly as today** — this extension is purely additive.

**Result:** for each LR channel you now emit up to three columns —
`…|autocrine`, `…|juxtacrine`, `…|paracrine` — all indexed by the **same patch `p`**.
They merge into one patch feature vector automatically. No smaller patches are ever
created.

### 2.4 Efficiency

Incremental cost over the current pipeline:

- **Build `C` once:** `O(n_cells log n_cells)`; mask to `C̃`: `O(nnz(C))`.
- **Batch the neighbor aggregation by receptor cell type.** For receiver type `B`,
  precompute `T_B = C̃ · M_B` where `M_B ∈ {0,1}^(n × |genes_B|)` marks
  (cell is type B) ∧ (expresses gene g) for every single receptor gene at once. One
  sparse×sparse product per receptor cell type yields, for every single-gene receptor,
  the per-cell adjacent-receiver counts (`t = T_B[:,g]`). This mirrors the existing
  `patch_by_rec` precompute. Only **complex** receptors need a per-column `C̃·r`
  (a small minority).
- **Per column:** one element-wise `l ⊙ t` (`O(n_cells)`) + one `P·(…)` aggregation
  (`O(n_cells)`) — the same order as today's per-column work.

Net: one global graph build + one sparse product per receiver cell type + `O(1)` extra
work per column. **No `O(n²)` enumeration, no sub-patching** — a small constant-factor
overhead on top of the existing method. Extra memory is one `~6·n_cells`-nnz integer
graph.

### 2.5 Design choices to decide

1. **Which pairs get a juxtacrine column?**
   - *(recommended)* **Keep the database gate:** only `Cell-Cell Contact`-annotated
     pairs get a juxtacrine column — but now it is a *true contact* count via `C̃`,
     not an all-to-all count. Secreted pairs stay autocrine/paracrine. Minimal column
     growth and biologically principled (a membrane-bound ligand should not "signal"
     at diffusion range).
   - *(optional)* **Distance-resolve every pair:** emit both juxtacrine and paracrine
     columns for all pairs. Richer, but ~doubles columns and adds noise for
     contact-only pairs. Offer as a flag.
2. **Contact-pair patch assignment.** Assign to the **sender's** patch (used above),
   consistent with the directed L→R semantics. With the intra-patch graph `C̃`, sender
   and receiver share a patch anyway, so this is unambiguous.
3. **Global vs intra-patch graph.** `C̃` (intra-patch) gives an exact, non-negative
   partition but misses two touching cells split by a grid line (a boundary artifact —
   though the patch tiling already has this). Mitigations if that matters: use
   overlapping/halo patches for the contact term only, or use the global `C` and clamp
   `Para_p = max(0, s_p·r_p − a_p − J_p)`. Default to `C̃` for clean bookkeeping.
4. **Graph type / contact scale.** Delaunay (pruned) adapts to density and is
   parameter-free; a fixed radius `r_c` is simpler when you have a physical contact
   scale in mind. Both plug into the same formulas.

### 2.6 Alternatives considered (ranked)

1. **Contact graph inside the existing matrix pipeline (recommended)** — §2.1–2.4.
   Exact partition, one graph built once, negligible overhead, no sub-patching.
2. **Multi-radius / distance-band decomposition.** Build graphs at radii
   `r₁ < r₂ < …`; band `k` count = `graph(r_k) − graph(r_{k−1})`. Generalizes the
   contact/diffusion split into several concentric rings. Same machinery, more
   columns and more graphs — use only if a full distance profile is wanted.
3. **Sub-patching for juxtacrine** — explicitly what we want to avoid; creates a second
   patch grid and misaligns rows.
4. **Explicit per-patch pairwise enumeration** — `O(Σ_p n_p²)`; defeats the purpose of
   the marginal-product trick. Rejected.

### 2.7 Summary

The current method counts `s_p · r_p` co-presence within a 50 µm patch and labels
channels juxtacrine/paracrine from the database, so it cannot physically separate
contact from diffusion. Inserting a single, globally-built, intra-patch-restricted
**cell contact graph `C̃`** lets you compute a true juxtacrine count
`J = P·(l ⊙ (C̃·r))` and recover paracrine as the exact remainder
`s_p·r_p − a_p − J_p`. Both live on the same patch row, so they merge into one patch
feature vector with no additional patching and only a small constant-factor cost.
