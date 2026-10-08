# CellChat — environment / build deviations (iris port)

**Scope note.** CellChat's *tutorial-call* deviations live in `cellchat/NOTES.md`, not here —
`METHODS.md` states explicitly that this method "has no separate `DEVIATIONS.md`". This file is
narrower and newer: it records what was found when porting the **environment build** to the MSK
iris cluster (linux-64) on **2026-08-17**. Nothing about the analysis contract is restated here.

**Status 2026-08-18 (second update): the env is BUILT and CellChat loads.** `comp-cellchat` exists,
CellChat **2.2.0.9001** is installed from the pinned commit and `CellChatDB.human$interaction` has
3,233 rows. Getting there took one more fix than the script had — see **D-8**, which is now folded
back into `build_env_iris.sh`. Stage 04 is runnable.

**Status update 2026-08-18 (first): the build is now scripted and pre-validated, but NOT yet run.**
`cellchat/build_env_iris.sh` was added. It builds the env end to end and resolves D-1, D-3, D-4
and D-5 below; the conda solve was verified by dry run, and the CellChat pin was verified against
the upstream remote. Nothing has been installed and `comp-cellchat` still does not exist.
The original diagnostic text follows; **D-1/D-2/D-3/D-4/D-5 are superseded by that script** and
are kept because they are the reasoning behind it.

**Status (original): stage 04 (CellChat) is PARKED.** `comp-cellchat` was deliberately **not** built in this
session. The findings below are diagnostic only. **No spec in `env.lock.yml` was changed, and no
package was installed.**

---

## D-1. `env.lock.yml` does not solve on linux-64 — three entries block it

`conda env create -f cellchat/env.lock.yml` fails on iris. Determined by iteratively
dry-run-solving a scratchpad copy of the lock's dependency list (`conda create --dry-run`, which
writes nothing); the tracked lock was never modified.

Exactly **three** specs have no linux-64 build on conda-forge:

| spec | what it is |
|---|---|
| `cctools=1030.6.3` | Apple binutils (`ar`, `ranlib`, `lipo`) — the bare macOS metapackage |
| `ld64=956.6` | Apple linker — the bare macOS metapackage |
| `libintl=0.25.1` | macOS gettext runtime; on linux this role is filled inside glibc |

Remove exactly those three from the spec list and **the solve succeeds**. That is the extent of
the finding — the removal was performed only on the scratchpad copy, to identify the blocking set.
Whether to change the lock is an open decision, not a recommendation made here.

## D-2. The osx-arm64 toolchain in the lock is a lock-generation artifact, not a dependency

The lock carries **19** macOS-specific entries:

```
cctools  cctools_impl_osx-arm64  cctools_osx-arm64
clang  clang-19  clang_impl_osx-arm64  clang_osx-arm64
clangxx_impl_osx-arm64  clangxx_osx-arm64
compiler-rt  compiler-rt_osx-arm64
gfortran  gfortran_impl_osx-arm64  gfortran_osx-arm64
ld64  ld64_osx-arm64  libcxx  libcxx-devel  libgfortran-devel_osx-arm64
libsigtool  sdkroot_env_osx-arm64  sigtool-codesign  tapi
```

None of these is a CellChat dependency. They are what `conda env export` captured because the
source environment was built on an **Apple-silicon Mac**, where the `compilers` metapackage
resolves to the native clang/gfortran toolchain. On linux-64 the same metapackage resolves to
`gcc_linux-64` / `gxx_linux-64` / `gfortran_linux-64`. The macOS entries are therefore
**provenance of the machine the lock was made on**, carrying no information about what CellChat
needs to build.

Note that most of them *do* exist for linux-64 — conda-forge ships them as **cross-compilation**
artifacts (targeting osx-arm64 *from* linux). Their presence in a linux-64 solve is therefore not
an error conda will report; see D-3 for why that matters.

## D-3. Activation-order hazard: two toolchains, one set of `CC`/`CXX`/`FC`

With the three D-1 entries dropped, the successful linux-64 solve pulls in **both** toolchains:

- the osx-arm64 cross set — `clang_osx-arm64`, `cctools_osx-arm64`, `ld64_osx-arm64`,
  `sdkroot_env_osx-arm64`, `tapi`, `sigtool-codesign`, `libsigtool`, `libcxx`, `libcxx-devel`
- the linux native set — `gcc_linux-64 14.3.0`, `gxx_linux-64 14.3.0`, `gfortran_linux-64 14.3.0`

Both ship activation scripts under `$CONDA_PREFIX/etc/conda/activate.d/` and **both write
`CC`, `CXX`, `FC`, `AR`, `LD`** (the osx set additionally writes `SDKROOT`). They are sourced in
glob order, so which compiler R sees is determined by filename ordering, not by intent. If the
osx-arm64 cross-compiler wins, every C/C++/Fortran build emits Mach-O objects that the linux
linker cannot use.

This is not hypothetical exposure of one or two packages. `cellchat/r_packages.lock.csv` pins
**250** R packages, and the conda lock supplies almost none of them (see D-4) — so essentially
all of that set is built **from source** through whatever `CC`/`CXX`/`FC` resolve to.
`cellchat/activate_env.sh` sources those same activation scripts.

## D-4. `install_envs_iris.sh`'s `build_cellchat` is not the right build route

Recorded here because it is the path the runbook (§1.3) points at.

- `cellchat/env.lock.yml` contains **none** of CellChat's R dependencies — no Seurat, NMF,
  ComplexHeatmap, igraph, ggalluvial, circlize, patchwork, presto. It has `r-base=4.3.3` plus
  ~60 low-level `r-*` packages only.
- `install_envs_iris.sh:93-118` (`build_cellchat`) does `from_lock` and then
  `remotes::install_local(CELLCHAT_SRC)`. It **never invokes `cellchat/install_env.R`**, which is
  the script that actually installs those packages, from a dated Posit CRAN snapshot
  (`https://packagemanager.posit.co/cran/2024-06-01`, reachable from iris compute nodes).
- `cellchat/r_packages.lock.csv` (250 rows, pinning `CellChat 2.2.0.9001`, `Seurat 5.1.0`,
  `NMF 0.27`, `ComplexHeatmap 2.18.0`, `ggalluvial 0.12.6`, `igraph 2.0.3`, `presto 1.0.0`) is
  referenced by no installer at all.

The intended route is `install_env.R` + `r_packages.lock.csv`.

`install_env.R:23` also defaults `--cellchat-src` to `/Users/jiayifan/tansey_lab/CellChat`,
a path that does not exist on iris; `--cellchat-src` must be passed explicitly.

## D-5. The CellChat source clone is not available on iris

CellChat installs as `RemoteType: local` — no lock file can carry it, and the source machine
(the laptop the lock was generated on) is **in repair and unreachable** as of 2026-08-17. A
filesystem search of `/home/fanj2`, the project directory (then `/data1/tanseyw/fanj2`, moved to
`/data1/tanseyw/projects/fanj2` later the same day) and `/data1/tanseyw/projects` found no
CellChat package tree and no `DESCRIPTION` with `Package: CellChat`.

**The pin does not need to be reconstructed — it is already recorded in this repo**, in two
places, and does not depend on the laptop:

- `cellchat/NOTES.md:3-4` — Version **2.2.0.9001**, git
  **`75253cd0c9e68410e6e721a6d3a0419a1d7e358f`**, 2026-03-04, commit message *"Update analysis.R"*
- `METHODS.md` (CellChat section) — `` `75253cd0` (2026-03-04) ``

So recovery is a clone of the upstream CellChat repository at `75253cd0c9e6…`, verified against
`DESCRIPTION`'s `Version: 2.2.0.9001`, rather than a date-bounded search. The 2026-03-04 commit
date is comfortably inside the pre-2026-07-28 bound implied by the CellChatDB re-export, so that
bound is consistent with, and superseded by, the recorded SHA.

## D-6. Resource note for whoever picks this up

The ~250 source builds of D-3 are the reason stage 04's environment build was deferred rather
than attempted: it does not fit alongside the other five envs in one short interactive
allocation. Build it in its own long allocation, and check `Rscript -e 'Sys.getenv("CC")'`
**before** starting the package installs.


---

## D-7. How the build was resolved (2026-08-18) — `build_env_iris.sh`

Added `cellchat/build_env_iris.sh`, a single idempotent script. What it settles, and how each
was verified:

| D-1/D-2 blocking + artifact entries | **Resolved by dropping all 30 macOS/clang entries, not just the 3 blockers.** The generated spec list is 143 of the lock's 173 dependencies and is written to `_archive/cellchat_linux64_specs.txt` with the dropped list in its header, so the diff from the lock is a tracked artifact rather than a step someone has to remember. |
|---|---|
| **D-3 dual-toolchain hazard** | **Eliminated, not mitigated.** Verified by a real `conda create --dry-run` on 2026-08-18: the resulting solve contains **zero** `clang` / `cctools` / `ld64` / `tapi` / `sigtool` / `osx-arm64` packages and selects `gcc_linux-64 14.4.0`, `gxx_linux-64 14.4.0`, `gfortran_linux-64 14.4.0`, `binutils 2.46.1`, `sysroot_linux-64 2.28`, `r-base 4.3.3` — 171 packages. The script re-asserts this at run time twice: it greps its own dry-run output for macOS packages and refuses, and after activation it requires `CC`/`CXX`/`FC` to match `*linux-gnu*` **and** compiles a trivial C file, checking the output is ELF. D-6's instruction ("check `Sys.getenv("CC")` before starting the package installs") is therefore enforced by the script rather than left to the operator. |
| **D-4 wrong build route** | The script calls `install_env.R --cellchat-src <clone>` explicitly, never `install_envs_iris.sh`'s `build_cellchat`, and passes the `--cellchat-src` that `install_env.R:23` otherwise defaults to a nonexistent Mac path. |
| **D-5 missing source clone** | **Resolved without needing the laptop.** `git ls-remote` on 2026-08-18 shows `refs/heads/main` of **`github.com/jinworks/CellChat`** is exactly `75253cd0c9e68410e6e721a6d3a0419a1d7e358f` — the SHA recorded in `NOTES.md:3-4`. Note **jinworks**, not `sqjin`: `sqjin/CellChat` is the v1 repository and its head is a different commit (`e4f68625…`). The script clones into `/data1/tanseyw/projects/fanj2/src/CellChat`, checks out that SHA, and asserts `DESCRIPTION` reads `Version: 2.2.0.9001` before installing. |
| **`install_env.R` `Ncpus`** | Fixed (additive). `:19` used `parallel::detectCores() - 2`, which reads the **machine**, not the cgroup — 62 on an iris node — so it would launch 62 concurrent package builds inside an 8- or 16-CPU allocation, each spawning gcc. It now prefers `CELLCHAT_NCPUS`, then `SLURM_CPUS_PER_TASK`, and only falls back to the original expression off-cluster. Same class of bug as stLearn's `os.cpu_count()`; performance and memory only, no package version is affected. |

**What this build does NOT reproduce, and must be stated in the write-up:** the lock's
osx-arm64 build strings for low-level libraries. That is not reproducible on linux by
construction. Everything carrying scientific meaning **is** pinned — `r-base 4.3.3`, the ~90
`r-*` conda versions, the CRAN snapshot date `2024-06-01`, the 250 R package versions in
`r_packages.lock.csv`, and the CellChat commit. The script's final step diffs the installed
library against `r_packages.lock.csv` and prints every missing and version-differing package,
so the residual gap is measured rather than assumed.

**Resource:** run it in its own allocation —
`srun --pty -p componc_cpu --cpus-per-task 16 --mem=64G --time=12:00:00 bash`. ~250 source
builds; the 12 h and 64 G are estimates, not measurements, since the build has not been run.


---

## D-8. The first real build produced a HALF-BUILT env that every check called healthy

`build_env_iris.sh` ran, created the conda env, and installed **224** R packages — Seurat, NMF,
ComplexHeatmap, igraph, presto, remotes, all present. **CellChat was not.** Stage 04 then got as
far as writing its input (1,676,162 cells, 2 conditions, 19 cell types — that part worked) and died
on `Error in library(CellChat) : there is no package called 'CellChat'`.

### The chain, root cause last

```
CellChat            ERROR: dependency 'svglite' is not available
  └─ svglite        installation failed
      └─ systemfonts  ft2build.h: No such file or directory
          └─ configure got PKG_CFLAGS=   (empty)
              └─ pkg-config --cflags fontconfig  FAILED
                  └─ "Package 'expat', required by 'fontconfig', not found"
                      └─ expat.pc is not in the env
```

**Root cause: conda-forge splits expat in two.** `expat` carries the headers, the `.so` symlink and
`lib/pkgconfig/expat.pc`; `libexpat` carries only `lib/libexpat.so.1`. The lock records **only
`libexpat=2.8.1`** — verified from `conda-meta/libexpat-2.8.1-hecca717_1.json`, whose entire file
list is `['lib/libexpat.so.1', 'lib/libexpat.so.1.12.1']`. And `fontconfig.pc` declares
`Requires.private: freetype2 >= 21.0.15, expat`. With `expat.pc` absent the whole fontconfig query
fails, so `systemfonts`' configure falls back to an empty `PKG_CFLAGS` and never sees
`$PREFIX/include/freetype2`. On the source Mac the dev half was present, so the lock never recorded
the need.

**Fix: add `expat` to the derived spec list.** Done in `build_env_iris.sh`, as an explicit `added`
entry that is written into the generated spec file as `#   ADDED expat`, so the diff from the lock
stays a tracked artifact. Repair on the already-built env was
`conda install -n comp-cellchat -c conda-forge expat`, after which `systemfonts 1.1.0`,
`textshaping 0.4.0`, `svglite 2.1.3` and `CellChat 2.2.0.9001` all built.

### A false alarm, recorded so nobody re-raises it

Mid-diagnosis `pkg-config --cflags freetype2` appeared to return **`/usr/include/freetype2`** — the
*system* freetype — which would have meant all 224 packages were built against system libraries.
**It is not true.** Re-tested in a clean shell (`env -i` + only `activate_env.sh`) it returns
`$PREFIX/include/freetype2` correctly; the system path came from the diagnosing shell having base
miniforge3 active. `pkg-config --variable pc_path` confirms the compiled-in search path is the env's
own directories. `PKG_CONFIG_PATH` is unset and does not need to be set.

### Two guards added so this cannot recur silently

1. **`build_env_iris.sh` step 4b — a pkg-config assertion before the ~250 builds**, alongside the
   existing toolchain assertion. `fontconfig`, `freetype2`, `harfbuzz` and `fribidi` must each
   resolve to a non-empty `--cflags`, or the script aborts. All four were the ones that mattered:
   they are what `systemfonts` / `textshaping` / `svglite` discover their headers through.
2. **The build script now EXITS NON-ZERO if CellChat is absent** at verification, instead of
   printing `CellChat loadable: FALSE` and saying DONE.

### And the check that let it through

`_common/luad_config.sh`'s `RS_CELLCHAT` test was `[ -x "$RS_CELLCHAT" ]` — it only asked whether an
`Rscript` binary exists. Once the conda env was created that went green, while the R library had no
CellChat. **`luad_check` now probes the package itself** (`requireNamespace("CellChat")`), prints its
version on success, and on failure says the env is half-built and points at the idempotent rebuild.
The probe is skipped silently when the env does not exist, so it costs nothing before the build.

**Generalisable lesson for the other envs:** every remaining interpreter check in `luad_config.sh` is
still "is the binary executable". That is sufficient for the Python envs only because they were
verified by import at build time. If any of them is ever rebuilt, verify by import, not by `-x`.

---

## D-9. `--plan multicore` on LUAD — a backend deviation, proven not to move a number

**The tutorial's call is `future::plan("multisession", workers = 4)`** (`NOTES.md:92`), and that is still
the default in `run_cellchat.R`. The GBM runs used it. **LUAD does not**, because it could not:

| run | plan | outcome |
|---|---|---|
| 2026-08-19 | `multisession`, workers=4, `--mem=192G` | OOM in the **AIS** condition, `MaxRSS` 192.0 GiB |
| 2026-08-20 | `multicore`, workers=4, `--mem=192G` | **AIS finished** (2166.7 s); OOM in **LUAD**, `MaxRSS` 201,326,424K = 191.99989 GiB |

`multisession` starts PSOCK workers, so every `future_sapply` call serialises its globals; future's own
error named the payload — *"The total size of the 11 globals exported is 2.81 GiB. The three largest
globals are 'data.use' (2.63 GiB)"*. `multicore` forks, so nothing is serialised. MEASURED on a
4,200-cell replica: **multisession 172.2 s vs multicore 23.3 s (7.4x)**, sequential 16.8 s.

**Why the backend cannot change a result.** Five independent checks, 2026-08-20:

1. **The only RNG in `computeCommunProb` is in the parent.** Deparse lines 158-159:
   `set.seed(seed.use); permutation <- replicate(nboot, sample.int(nC, size = nC))`. It sits *above*
   both `my.sapply` calls (line 160, once per condition; line 214, once per LR pair). Both worker
   bodies only index it: `group[permutation[, nE]]`.
2. **No RNG is reachable from a worker.** Grepping every function in the CellChat namespace for
   `sample|runif|rnorm|set.seed|replicate|...` finds RNG only in `computeCommunProb` itself and in
   `computeEnrichmentScore` / `netAnalysis_contribution` / `netVisual*` / `runPCA` / `runUMAP`.
   `computeExpr_LR`, `computeExpr_complex`, `computeExpr_coreceptor`, `computeExpr_agonist`,
   `computeExpr_antagonist`, `triMean`, `thresholdedMean`, `computeRegionDistance` — everything the
   two bodies call — contain none.
3. **The fork-vs-PSOCK RNG channel is real, so the test can fail.** MEASURED:
   `future_sapply(1:8, \(i) runif(1))` returns three *different* vectors under
   sequential / multicore / multisession. The plan **would** change numbers if any body drew RNG.
   `formals(future_sapply)$future.seed` is `NULL`, i.e. future actively detects and warns. A
   zero-warning result inside `computeCommunProb` is therefore a real falsification, not a silence.
4. **No BLAS reduction-order channel either.** Every `crossprod` in either body has inner dimension
   **k = 1** (`matrix(x, nrow = 1)`, deparse 169-170 and 227) — a rank-1 outer product with one
   multiply per output element and no accumulation, hence no reduction order to permute (MEASURED
   bitwise identical to `outer()`). `04_cellchat.sh` also exports `OMP/OPENBLAS/MKL_NUM_THREADS=1`
   before R starts, inherited identically by forks and PSOCK workers.
5. **MEASURED end to end on a real CellChat object** — 4,200 cells, 6 groups, 2 samples, spatial
   datatype, 104 over-expressed interactions, **nboot = 100** (the runner's real value), the runner's
   exact `computeCommunProb` arguments, three fresh R subprocesses:

   | comparison | slot | `identical()` | max abs diff |
   |---|---|---|---|
   | sequential vs multisession | prob / pval | TRUE / TRUE | 0 / 0 |
   | sequential vs multicore | prob / pval | TRUE / TRUE | 0 / 0 |
   | multisession vs multicore | prob / pval | TRUE / TRUE | 0 / 0 |

   Non-degenerate: 1860/3744 `prob` > 0, 210 entries `pval < 0.05`, 93 distinct p-values,
   min non-zero p = 0.01 = 1/nboot. The permutation test is live.

**Known coverage gap, stated honestly.** The 104 pairs in the replica are all `Secreted Signaling`
and all simple — **0 agonist, 0 antagonist, 0 co-receptor**, so deparse 181-183 (the contact
`P.spatial * adj.contact` mutation) and 184-203 / 229-248 (the co-factor branches) never executed.
The real tier-B DB is 535 Cell-Cell Contact + 424 ECM-Receptor + 1280 Secreted, with 342 / 459 / 382 /
486 pairs carrying agonist / antagonist / coA / coI. Those branches are covered by arguments 1, 2
and 4 rather than by the end-to-end measurement. The gap is **branch coverage, not scale**.

**How to state it in the manuscript:** a parallel-backend change, not a modelling change. GBM used the
tutorial's `multisession`; LUAD used `multicore` because the LUAD object does not fit otherwise. If
even that is unwanted, `--workers 1` is a stronger claim still — see D-10.

---

## D-10. The `UNRELIABLE VALUE` warnings are `igraph::eigen_centrality`, they pre-date the port, and nothing consumes the affected number

Every CellChat run in this benchmark — GBM and LUAD, `multisession` and `multicore` — emits per
condition:

```
UNRELIABLE VALUE: One of the 'future.apply' iterations ('future_sapply-N') unexpectedly generated
random numbers without declaring so. There is a risk that those random numbers are not
statistically sound and the overall results might be invalid. ... specify 'future.seed=TRUE'
```

It is alarming and it is **not** about the communication probabilities.

**It pre-dates the iris port.** `logs/comparators/cellchat-GBM-{default,cellchatdb2}-run.log` were
produced on macOS (`arm64-apple-darwin20.0.0`) under the tutorial's `multisession`, before `--plan`
existed, and each carries the identical **8** warnings — `future_sapply-1..4` x 2 conditions.

**It is localised.** Per-stage `withCallingHandlers`, MEASURED, identical under both parallel plans:

| stage | RNG warnings |
|---|---|
| `identifyOverExpressedGenes` | 0 |
| `identifyOverExpressedInteractions` | 0 |
| **`computeCommunProb`** | **0** |
| `filterCommunication` | 0 |
| `computeCommunProbPathway` | 0 |
| `aggregateNet` | 0 |
| **`netAnalysis_computeCentrality`** | **4** (= `nbrOfWorkers()`) |

`netAnalysis_computeCentrality` has its own `my.sapply` (its deparse lines 16-18:
`my.sapply <- ifelse(nbrOfWorkers()==1, pbapply::pbsapply, future.apply::future_sapply)`), called
**once per condition** over `1:nrun` pathways and chunked into `nbrOfWorkers()` futures — hence
exactly 4 warnings, and hence the GBM logs' two runs of four *consecutive* future counters
(`…-225,226,227,228` and `…-497,498,499,500`). Four consecutive counters is the signature of one
once-per-condition call; the `nLR`-times call at `computeCommunProb` deparse 214 would have produced
thousands.

**The consumer is one line.** Snapshotting `.Random.seed` around each centrality measure (MEASURED,
igraph 2.0.3 / sna 2.7.2): `hub_score`, `authority_score`, `page_rank`, `betweenness`, `strength`,
`sna::flowbet`, `sna::infocent` all leave it untouched; **`igraph::eigen_centrality` advances it** —
ARPACK with a random start vector. `CellChat:::computeCentralityLocal` deparse line 12 is
`centr$eigen <- igraph::eigen_centrality(G)$vector`. That is the entire story.

**Its numerical consequence, and why it does not matter here.** MEASURED on one post-`aggregateNet`
object: `centr$eigen` differs across plans at **max 3.22e-14** (148-153 of 2442 entries), and is not
even run-to-run reproducible under `multicore` (run1 vs run2, 3.08e-14). Every other measure —
`outdeg`, `indeg`, `hub`, `authority`, `page_rank`, `betweenness`, `flowbet`, `info` — is bitwise
identical across all plans and repeats. And **`eigen` is never reported**:

- `cellchat_io.R:89-92` writes only `outdeg`, `indeg`, `flowbet`, `info` to `<cond>_centrality.csv`.
- `plot_cellchat.R` contains no reference to `centr` or `eigen`.
- In the whole CellChat namespace only `computeCentralityLocal` (which creates it) and
  `netClustering` mention `eigen`, and `netClustering`'s is base R `eigen(L, symmetric = TRUE)` on a
  network-similarity Laplacian — an unrelated quantity.

So `centr$eigen` is computed, stored in `objects/<cond>.rds`, and consumed by nothing.

**Decision: leave the warning in the logs.** Setting `future.rng.onMisuse = "ignore"` would silence
the authors' own diagnostic, it is present in the GBM results too, and it would hide a genuine change
if CellChat ever grows one. Do not let a reader of the log conclude the p-values are unsound — point
them at this section.

---

## D-11. Resume ergonomics — `--reuse-existing`, `--skip-export`, `--skip-plots`, `--conditions`

Added 2026-08-20, all additive; a bare `bash 04_cellchat.sh` is unchanged (verified by dry-run diff).

`run_cellchat.R` had **no skip logic at all**: `for (cond in conditions) results[[cond]] <- run_condition(cond)`
recomputes every named condition and overwrites `objects/<cond>.rds` plus all 12 `quant/<cond>_*`
files. After the 2026-08-20 OOM that meant a naive re-run would have destroyed a completed AIS
(2166.7 s) to get at LUAD. `--conditions LUAD` avoids the recompute but then `run_manifest.json` —
written once, after the loop, from every element of `results` — would describe only one condition.

| flag | file | why |
|---|---|---|
| `--reuse-existing` | `run_cellchat.R` | loads `objects/<cond>.rds` instead of recomputing, and rebuilds `stats` / `req` / `d_obs` / `n_cells` / `n_samples` **from the loaded object**, so the manifest still describes both conditions. Records `reused: true` and `wall_seconds: null`. |
| `--conditions` | `04_cellchat.sh` | was hardcoded `AIS,LUAD` with no passthrough |
| `--skip-export` | `04_cellchat.sh` | `prepare_gbm_input.py` has no exists-check and rewrites 8.7 GB of exchange mtx on every invocation. Refuses if `$IN/<cond>/data.mtx` is absent. |
| `--skip-plots` | `04_cellchat.sh` | `plot_cellchat_luad.R` only does `readRDS(objects/<cond>.rds)` (`plot_cellchat.R:56`), so it runs standalone later. The GBM plot tree was 9,256 files; do not let it take a 3-hour compute down with it. |

MEASURED verification of the reuse path against the completed AIS: **all 12 `quant/AIS_*` files
byte-identical** (md5, gzipped ones included), `19 cell types / 709 LR pairs / 3737 significant links /
70 pathways` reproduced exactly, `observed_min_cell_distance_um` 0.0069 vs the logged 0.007, 41 s.

**One hazard found and fixed during that verification.** `save_cellchat_quant` (`cellchat_io.R:48`)
begins with a bare `saveRDS`, so the first version of the reuse path re-serialised the 329 MB object
it had just read — i.e. a resume feature was rewriting the artefact it exists to protect, and a kill
during that write would have destroyed a completed condition. `save_cellchat_quant` now takes
`save_rds = TRUE` (default; every existing caller unchanged) and the reuse path passes `FALSE`.
Verified: `objects/AIS.rds` mtime/size unchanged across a reuse run, quant still 12/12 identical.

**Provenance note.** During that same first, pre-fix test the canonical
`results/cellchat/LUAD/cellchatdb2/objects/AIS.rds` was re-serialised in place through a symlink
(mtime 2026-08-20 14:00:30, was 00:53). It is a round-trip of the same in-memory object: size
unchanged to the byte (329,468,061), and a subsequent reuse run regenerated all 12 quant files
byte-identically, `sum(net$prob) = 39.6898131431063`, 3737 significant links, 70 pathways,
centrality present. The content is intact; the bytes are not the originals. Recorded here rather than
quietly fixed.
