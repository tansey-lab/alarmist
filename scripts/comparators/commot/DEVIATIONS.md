# COMMOT — deviations from the tutorial

The tutorial is a **single Visium section** (3,355 multi-cell spots, coordinates in full-resolution
image pixels); ours is **single-cell Xenium on a 13-core TMA** in microns. Every deviation below
follows from that, from CellChatDB v2 being newer than the bundled v1, or from a version gap.

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| numpy | (era-appropriate) | **pinned `numpy<2` (1.26.4)** | **Hard version gap.** `commot/_optimal_transport/_usot.py` uses `np.Inf` as a *module-level default argument* (3 occurrences), removed in NumPy 2.0 — so `import commot` raises `AttributeError` and the package is entirely unusable on numpy ≥ 2. Pinned rather than patching installed code. Anyone reproducing this needs the same pin. |
| `dis_thr` | `500` | **`365`** | **The tutorial's own prose is wrong about its units.** `Basic_usage.ipynb` calls it "a spatial distance constraint of 500 µm", but `sc.datasets.visium_sge` returns `obsm['spatial']` in **full-resolution image pixels**. Measured on the actual tutorial dataset: nearest-neighbour spacing **137.0 units** against Visium's 100 µm hex pitch → **0.7299 µm/unit** → `dis_thr=500` is **365 µm physically**. We use 365 on our micron coordinates to preserve the authors' real constraint. (A cross-check via `spot_diameter_fullres` gives 0.6144 µm/unit, a 19% disagreement; the hex pitch is a hard geometric constant and is the more reliable calibration.) Copying `500` verbatim would have given a 37% wider neighbourhood than the authors used. Deliberately not harmonised with CytoSignal's 200 µm, stLearn's 250 µm or SpatialDM's 135 µm. |
| `signaling_type` | `'Secreted Signaling'` only | **`Secreted Signaling` + `Non-protein Signaling`** | COMMOT applies **one** `dis_thr` to every pair and has no separate contact/adjacency graph, so the authors' restriction to diffusible signalling is methodologically load-bearing — including Cell-Cell Contact or ECM-Receptor would model juxtacrine signalling as acting over 365 µm. We keep that principle but add v2's `Non-protein Signaling`, which **did not exist in the CellChatDB v1 the tutorial used**, is equally diffusible (neurotransmitters/metabolites), and matters in glioma — SpatialDM ranked glutamatergic pairs as core 13's entire top-6. Consistent with the CytoSignal mapping (Secreted/ECM/Non-protein → diffusion). 2,259 pairs before expression filtering; 862 on the Xenium panel. |
| LR filter | `filter_criteria='min_cell_pct'`, `min_cell_pct=0.05` | **`filter_criteria='min_cell'`, `min_cell=100`** | Both are COMMOT function defaults; we use the other criterion. The tutorial's 5% is 5% of **Visium spots**, each aggregating 10–30 cells and therefore with much higher detection rates. At single-cell resolution the same threshold retains only **0.9–1.8%** of pairs (20–40 per core) against the tutorial's **20.9%** (1,199 → 250). `min_cell=100` lands at 69–217 pairs per core, comparable to the tutorial's absolute 250. |
| cores 2 and 6 | n/a | **not analysable — 0 pairs** | Under `min_cell=100` these two cores retain **zero** LR pairs. They are the two least dense (423 and 904 cells/mm²), and they fail under any defensible absolute floor; a percentage floor would only pass them by admitting pairs estimated from <25 expressing cells, i.e. noise. **Reported as a limitation of the method on sparse tissue rather than worked around by lowering the threshold until a result appears.** 11 of 13 cores are analysed. |
| input | `adata.raw = adata`; `normalize_total`; `log1p` | same, applied to `layers['counts']` | Our `.X` is already log-normalised by an unknown recipe. We start from genuine raw counts (verified: `layers['counts']` max = 651, integers) and reproduce the tutorial's exact two steps. The stale `uns['log1p']` key left by upstream processing is popped first — otherwise `sc.pp.log1p` emits a spurious "already log-transformed" warning even though the input is raw. |
| whole-tissue vs per-core | single section | **one run per TMA core** | See the provenance note below — this was **not** a tutorial instruction and **not** a user instruction; it was forced by memory. |
| per-LR plots | top pathway only | top 3 by total received signal **plus** GRN→SORT1 and ANXA1→FPR1 | Standing request; requested pairs are written as `requested_rank<N>_<LR>.png`. |
| `communication_direction` / `plot_cell_communication` | run for one pathway | **not run** | The vector-field plot requires `background='image'` (an H&E image we do not have) or heavy grid interpolation per pathway; the sender/receiver spatial maps carry the same information for our purposes. Noted as an available output not exercised. |
| `communication_deg_detection` | run for PSAP | **not run** | Downstream signalling-DE analysis, out of scope for this pass. |

---

# 2026-08-10 re-run — what changed and why

The table above describes the **2026-08-01** run. That run's artifacts are preserved at
`results/comparators/commot/GBM/cellchatdb2_prefix_backup_20260810/`. Everything below is a
deviation *from that run*, decided with the user. Scope: **GBM / `cellchatdb2` only** — the
`default` (bundled CellChatDB v1) tier and LUAD remain not run.

## Defects in the 2026-08-01 artifacts that forced the re-run

The script was edited ~2 minutes after that run finished (`run_commot.py` mtime 21:07 vs core 1's
outputs at 21:04-21:05) and **never re-executed**, so the on-disk artifacts predate the fix:

| Defect | Evidence | Effect |
|---|---|---|
| `lr_total_received.csv` mixed pathway aggregates in with LR pairs | 261 rows vs 217 in `lr_pairs_used.csv` (core 1); same gap in all 11 cores | `per_split_summary.csv`'s `GRN_SORT1_rank` / `ANXA1_FPR1_rank` were inflated — GRN→SORT1 read 5–16 (median 8) when the LR-only rank is **1–4 (median 2)** |
| `pathway_total_received.csv` never written | absent in all 11 cores | pathway ranking not persisted (recoverable, since the values were mixed into the other file) |
| "top 3" plots were **pathway aggregates**, not LR pairs | 31 of 33 `signal_*.png` are pathway names (BMP, PDGF, ncWNT, FGF, IGF, COMPLEMENT, Glutamate, GAS); only cores 4 and 12 got one real pair | the method's own top LR pairs were never plotted at all |
| `cluster_communication` ran on the **alphabetically** first pathways | core 1 tested ACTIVIN / ANGPTL / ANNEXIN, which rank 29 / 6 / 32 by signal (true top 3: BMP, ncWNT, PDGF); every core shows the same alphabetical pattern | COMMOT's *only* p-values were spent on near-irrelevant pathways |
| `obsp` never persisted | `--save-obsp` was declared at argparse and **never referenced in the body** — a dead flag | every downstream COMMOT tool required a full OT re-run |
| manifest recorded the wrong filter parameter | `"min_cell_pct": 0.05` written, but `filter_criteria=min_cell, min_cell=100` was in force; neither of those was recorded | reading the manifest alone gives the wrong method |

METHODS.md's *narrative* numbers were correct — they had been recomputed by hand from the CSVs —
but nothing was written back to disk, so the artifacts and the write-up disagreed.

## Deviations introduced by the re-run

| Item | 2026-08-01 | 2026-08-10 | Why |
|---|---|---|---|
| **expression-filter scope** | per core | **`--filter-scope global`** — `filter_lr_database` evaluated once on all 100,197 cells, the resulting **671** pairs handed to every core | Per-core filtering gave each core a different pair set (51–217) and dropped cores **2 and 6** entirely. **Both dropped cores are low-grade**, so the analysable set was 7 high + 4 low against the TMA's true 7 + 6. Global scope restores **13/13 cores**, makes ranks comparable across cores, and removes the `corr(n_cells, n_pairs_used) = 0.819` confound. The tutorial is a single section and does not distinguish the two scopes. Cost: ~3× the pairs per core, so runtime/memory/storage all rise. |
| **top-N LR plots** | top 3 of a pathway-contaminated ranking | top 3 of the **LR-pair-only** ranking (`--n-top-lr 3`); pathway aggregates drawn separately as `pathway_<name>.png` | pathway sums are larger than any single pair by construction, so the contaminated ranking could only ever surface pathways |
| **`cluster_communication` targets** | 3 pathways, selected in column (alphabetical) order | **top 5 pathways by total received signal** (`--n-top-pathways 5`, user's choice of N) **plus each requested LR pair** via `lr_pair=` | puts the permutation test on the pathways that carry signal; `lr_pair=` is the only significance test COMMOT can give GRN→SORT1 and ANXA1→FPR1 |
| **`cluster_communication` seed** | not passed (silently used the default `random_seed=1`) | **`random_seed=--seed`** | the permutation p-values were not tied to the run seed |
| **`obsp` persistence** | none | **`--save-adata`** writes `<core>/adata_commot.h5ad` with all transport plans (`.raw` dropped — reconstructible from `--h5ad`) | so `communication_direction`, `communication_impact`, `deg_detection`, `spatial_autocorrelation`, `group_*` never need another OT run |
| **native `ct.pl` figures** | **none produced** — all 54 PNGs came from our own `sender_receiver_map` | `communication_direction` → **`plot_cell_communication`** (vector fields, sender + receiver, for the 5 pathways and both requested pairs) and **`plot_cluster_communication_network`** | see the correction below |
| `plot_cluster_communication_dotplot` | not run | **still not run — version gap**, replaced by `ours_dotplot_top.png` | see below |
| pygraphviz | absent | **added to `comp-commot`** (1.14) | `plot_cluster_communication_network` goes through `networkx.drawing.nx_agraph.to_agraph` |
| `env.lock.yml` | **missing** (the only comparator without one) | frozen | skill invariant |
| manifest | 11 fields, one of them wrong | adds `filter_criteria`, `filter_scope`, `min_cell`, `n_pairs_global`, `cluster_random_seed`, `n_cells_in`/`n_cells_analysed`, `git_sha`, full `versions` block | reproducibility |
| `peak_rss_gb` column | reported per core | renamed **`peak_rss_gb_running_max`** | `getrusage` is process-wide and monotonic, so it was never a per-core peak; the old column read as one |

## Correction: the "`plot_cell_communication` needs an H&E image" justification was wrong

The 2026-08-01 table (row `communication_direction` / `plot_cell_communication`) says the vector
field "requires `background='image'` (an H&E image we do not have)". Checked against the
installed package: `plotting/_plotting.py:35` has **`background: str = "summary"` as the
default** — no image is involved. Confirmed empirically: the re-run produces 14 vector-field PNGs
per core with no background image. The plots were skipped for no valid reason.

## Self-correction: the first vector-field pass used the wrong arguments

`run_commot.py` called `ct.pl.plot_cell_communication` with the **function's** defaults —
`plot_method='cell'`, `scale=1.0`, `ndsize=1`, no `normalize_v`. Those render without error, and
were initially reported as a success on that basis. Inspecting the output showed they are close
to useless: at 3k–26k cells per core the per-cell arrows are invisible and the `coolwarm`
background is swamped by near-zero values.

The **tutorial's** call is a different argument set entirely (`Basic_usage.ipynb`;
`visium-mouse_brain.ipynb` is the same call with `pathway_name='PSAP'`):

```python
ct.pl.plot_cell_communication(adata, database_name=..., lr_pair=..., plot_method='grid',
    background_legend=True, scale=0.00003, ndsize=8, grid_density=0.4, summary='sender',
    background='summary', clustering='leiden', cmap='Reds',
    normalize_v=True, normalize_v_quantile=0.995)
```

`scripts/comparators/commot/plot_commot_vf.py` re-draws every field with those values. It reads
each core's persisted `adata_commot.h5ad`, so **it costs no OT re-run** — the first concrete
payoff of `--save-adata`.

| Item | Tutorial | Ours | Why |
|---|---|---|---|
| `plot_method` | `'grid'` | `'grid'` | (the first pass wrongly used `'cell'`) |
| `normalize_v`, `normalize_v_quantile`, `grid_density`, `ndsize`, `background_legend` | `True, 0.995, 0.4, 8, True` | identical | — |
| `scale` | `0.00003` | **`0.00003 × (9000 / x_extent_µm)`**, i.e. 1.4e-4 for core 13 | `_utils/_plotting.py:318-320` passes `scale` to `quiver(scale_units='x')`, so arrow length is `|v|/scale` in **data units**. The tutorial's value is tuned to Visium full-resolution pixels (~9,000 units across); our cores are microns (~1,000–3,000 across). Copying it verbatim draws arrows orders of magnitude too long. We hold arrow length at the tutorial's *fraction of the field*. **Same class of unit trap as `dis_thr`.** |
| background | one panel `background='image'` (+`clustering`, `cmap='Alphabet'`), one `background='summary'` (`cmap='Reds'`) | `background='cluster'` (+`clustering='cell_type'`, `cmap='Alphabet'`) and `background='summary'` (`cmap='Reds'`) | we have no H&E, and `'cluster'` is the function's own no-image equivalent of the panel the tutorial gets from the image. Note this also shows the earlier "the vector field *requires* an image" claim was doubly wrong: the tutorial itself draws one of its two panels with `background='summary'`. |
| `clustering` | `'leiden'` | `'cell_type'` | we have real annotations; re-clustering would discard them |

## Version gap: `plot_cluster_communication_dotplot` cannot run here

commot 0.0.3 was written against matplotlib <3.9 and seaborn <0.13; installed are **matplotlib
3.10.9** and **seaborn 0.13.2**. Two independent breakages, both observed:

1. `plotting/_plotting.py:788` iterates `g.legend.legendHandles` — matplotlib removed that alias
   in 3.9 (now `legend_handles`).
2. Past that, seaborn 0.13 returns `Line2D` legend handles, not the `PathCollection` the code
   assumes, so `set_edgecolor` does not exist.

(1) is a one-line alias; (2) needs the function's internals rewritten, which would stop being
"the method as its authors wrote it". **Recorded, not patched.** The information is fully
persisted as `cluster_comm_*.csv` / `cluster_pval_*.csv`, and `ours_dotplot_top.png` draws an
equivalent from exactly those matrices. It is prefixed `ours_`, never `native_`, so it cannot be
mistaken for a COMMOT figure.

## Sharpening the cross-core comparability warning

The existing "Confound inherited from the dataset" section is correct but incomplete. Two further
reasons cross-core **magnitudes** are not comparable, both read off the installed solver:

- **The OT is globally normalised per run.** `_optimal_transport/_cot.py:269-271` computes
  `max_amount = max(S.sum(), D.sum())` and divides both marginals by it; `:335` multiplies the
  transport plan back by `max_amount`, so *units* are restored. But `eps_p` (entropy) and `rho`
  (unmatched-mass penalty) are **fixed constants acting on the normalised masses**, so the shape
  of the solution — how much mass goes unmatched, how diffuse the coupling is — depends on the
  cell set. Each core is normalised by its own constant.
- Consequence: **ranks within a core are sound and (with `--filter-scope global`) comparable
  across cores; absolute received-signal magnitudes are not.** Any grade contrast must be built
  on ranks or within-core statistics, not on raw magnitudes.

The "cross-core coupling is only 0.009% of pairs" figure remains true, but it only bounds the
*coupling structure* lost to splitting — it says nothing about the filter scope or the
normalisation, which are the two things that actually make per-core and whole-slide runs
non-equivalent.

## Note on what COMMOT does and does not provide

**There is no per-LR-pair significance test.** COMMOT produces transport plans and per-cell
sent/received amounts; the only p-values come from `cluster_communication`'s permutation test at
the **cell-type-pair** level, per pathway. So COMMOT's LR "ranking" here is by **total received
signal**, which is a magnitude, not a significance — it is not directly comparable to CytoSignal's
significant-cell counts, stLearn's significant-spot counts, or SpatialDM's FDR. State this
explicitly whenever COMMOT ranks are placed beside the others.

## Confound inherited from the dataset

The density–grade correlation documented in `../spatialdm/DEVIATIONS.md` (r = 0.659, p = 0.014;
4× difference in neighbourhood occupancy) applies here too: COMMOT's `dis_thr` is distance-based,
so the number of cells inside the transport radius — and hence the amount of signal that can be
moved — scales with local density. The `min_cell` filter compounds it, since denser cores also
retain more LR pairs (69–217, correlated with core size). **Cross-grade comparisons of COMMOT
signal magnitudes or pair counts are confounded and must not be read as pure biology.**

## Provenance of the per-core decision — recorded honestly

COMMOT was run per TMA core. Two distinct justifications were conflated when this was first
decided, and they deserve to be separated because they carry very different weight:

**FORCED — memory.** COMMOT materialises a **dense N×N distance matrix**
(`scipy.spatial.distance_matrix`, `_spatial_communication.py:390`). At whole-slide scale
(100,197 cells) that is **80.3 GB** for the distances alone, on a 36 GB machine. Per-core, the
largest (core 1, 26,456 cells) needs 5.6 GB and the run peaked at 13.05 GB. **Whole-slide COMMOT
is not possible here regardless of any methodological preference.** This alone settles it.

**NOT forced — consistency with SpatialDM.** SpatialDM was run per-core on an explicit user
instruction (Moran's I across disconnected cores would distort its analytical null). Extending
that to COMMOT was an inference, not an instruction, and COMMOT's statistic is not a Moran's I —
the null-structure argument does **not** transfer. It should not be cited as a reason.

**Neither COMMOT tutorial prescribes per-sample splitting** — both are single contiguous
sections, and COMMOT has no native multi-sample mode. So per-core here is a departure from the
authors' demonstrated usage, adopted purely because the alternative is computationally
impossible.

For completeness, cross-core coupling at `dis_thr = 365 µm` would have been negligible anyway:
8,792 of 94,949,273 pairs = **0.009%**. So the split costs essentially nothing in signal; it is
the *comparability* with the whole-TMA methods (CytoSignal, stLearn) that it costs — see the
warning at the top of `../METHODS.md`.
