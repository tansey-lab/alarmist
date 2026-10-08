#!/usr/bin/env bash
# ======================================================================================
# build_env_iris.sh -- build `comp-cellchat` on the MSK iris cluster (linux-64), end to end.
#
#   bash scripts/comparators/cellchat/build_env_iris.sh
#
# RUN IT ON A COMPUTE NODE, IN ITS OWN ALLOCATION. It builds ~250 R packages from source:
#
#   srun --pty -p componc_cpu --cpus-per-task 16 --mem=64G --time=12:00:00 bash
#
# WHY THIS SCRIPT EXISTS AND NOT `install_envs_iris.sh --only cellchat`
# --------------------------------------------------------------------
# install_envs_iris.sh:93-118 (`build_cellchat`) does `from_lock` then
# `remotes::install_local(CELLCHAT_SRC)` and NEVER invokes cellchat/install_env.R -- which is
# the script that actually installs CellChat's ~250 R dependencies from a dated CRAN snapshot.
# env.lock.yml contains none of them (no Seurat, NMF, ComplexHeatmap, igraph, ggalluvial,
# circlize, patchwork, presto): it is r-base plus ~90 low-level r-* packages. So that route
# produces an env that cannot load CellChat. See cellchat/DEVIATIONS.md D-4.
#
# THREE THINGS THIS SCRIPT HANDLES THAT A NAIVE `conda env create -f env.lock.yml` DOES NOT
# -------------------------------------------------------------------------------------------
# 1. env.lock.yml DOES NOT SOLVE ON linux-64. It was exported on an Apple-silicon Mac, so it
#    carries the macOS toolchain. Three specs have no linux-64 build at all -- cctools, ld64,
#    libintl (D-1) -- and ~27 more are osx-arm64 cross-compilation artifacts that conda-forge
#    *does* ship for linux-64, so they would install silently (D-2).
# 2. If those cross-compilers install, BOTH toolchains ship activate.d scripts and BOTH write
#    CC/CXX/FC/AR/LD. They are sourced in glob order, so which compiler R sees is decided by
#    filename ordering. If the osx-arm64 clang wins, all ~250 source builds emit Mach-O objects
#    the linux linker cannot use (D-3). This script drops every osx entry, so `compilers`
#    resolves natively and the hazard does not arise -- VERIFIED: the resulting solve contains
#    zero clang/cctools/ld64/tapi/sigtool packages and selects gcc/gxx/gfortran_linux-64 14.4.0
#    with binutils 2.46.1 and sysroot_linux-64 2.28. It is then ASSERTED again at step 4 below,
#    because a silent Mach-O toolchain is exactly the failure that would waste the allocation.
# 3. CellChat installs as RemoteType: local -- no lock can carry it, and the laptop the lock
#    came from is unavailable (D-5). The pin is recorded in-repo though: version 2.2.0.9001,
#    git 75253cd0c9e68410e6e721a6d3a0419a1d7e358f (cellchat/NOTES.md:3-4, METHODS.md). Verified
#    2026-08-18: that SHA is the current head of `main` in github.com/jinworks/CellChat -- note
#    JINWORKS, not sqjin (sqjin/CellChat is the v1 repo and its head is a different commit).
#
# WHAT IS PINNED AND WHAT IS NOT -- state this in the write-up
# -----------------------------------------------------------
# Pinned exactly: r-base 4.3.3, the ~90 r-* conda packages at the lock's versions, the CRAN
# snapshot date (2024-06-01, in install_env.R:18), the 250 R package versions in
# r_packages.lock.csv, and the CellChat commit. NOT reproduced: the lock's osx-arm64 build
# strings for low-level libraries -- they cannot be, on linux, by construction.
#
# Every step is idempotent: re-running skips what already exists.
# ======================================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../_common/luad_config.sh"

ENV_NAME="comp-cellchat"
ENV_PREFIX="$CONDA_ENVS/$ENV_NAME"
SRC_ROOT="${CELLCHAT_SRC_ROOT:-/data1/tanseyw/projects/fanj2/src}"
CC_SRC="$SRC_ROOT/CellChat"
CC_SHA="75253cd0c9e68410e6e721a6d3a0419a1d7e358f"
CC_REPO="https://github.com/jinworks/CellChat.git"
CC_VERSION="2.2.0.9001"
SPECS="$HERE/../_archive/cellchat_linux64_specs.txt"   # provenance copy, regenerated below

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
die()  { printf '\n\033[1mFATAL: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- 0. preflight
say "0. preflight"
[ -f "$HERE/env.lock.yml" ]        || die "env.lock.yml missing"
[ -f "$HERE/install_env.R" ]       || die "install_env.R missing"
[ -f "$HERE/r_packages.lock.csv" ] || die "r_packages.lock.csv missing"
command -v conda >/dev/null        || die "conda not on PATH"

NCPU=$(nproc)                       # affinity-aware: what the cgroup actually grants
echo "  node        : $(hostname -s)"
echo "  cpus (avail): $NCPU"
echo "  env prefix  : $ENV_PREFIX"
echo "  CellChat src: $CC_SRC @ ${CC_SHA:0:8}"
if [ "$NCPU" -lt 4 ]; then
    echo "  WARNING: only $NCPU CPUs. ~250 source builds will take many hours."
fi

# Both remotes must be reachable FROM THE NODE, not just from the login node.
say "0b. network reachability from this node"
timeout 30 git ls-remote --heads "$CC_REPO" >/dev/null 2>&1 \
    || die "cannot reach $CC_REPO -- check the node has outbound https"
echo "  github.com/jinworks/CellChat  OK"
SNAP=$(grep -oE 'https://packagemanager\.posit\.co/cran/[0-9-]+' "$HERE/install_env.R" | head -1)
timeout 30 curl -fsI "$SNAP/src/contrib/PACKAGES" >/dev/null 2>&1 \
    || die "cannot reach the CRAN snapshot $SNAP"
echo "  $SNAP  OK"

# ---------------------------------------------------------------- 1. spec list
say "1. derive the linux-64 spec list from env.lock.yml"
mkdir -p "$(dirname "$SPECS")"
python3 - "$HERE/env.lock.yml" "$SPECS" <<'PY'
import re, sys, pathlib
lock, out = sys.argv[1], sys.argv[2]
lines = pathlib.Path(lock).read_text().splitlines()
i = next(n for n, l in enumerate(lines) if l.strip() == "dependencies:")
specs = [l.strip()[2:].strip() for l in lines[i + 1:] if l.strip().startswith("- ")]
# macOS-only, or present only because the macOS clang toolchain pulled them in.
DROP = re.compile(
    r"^(cctools|ld64|libintl|clang|clang-19|compiler-rt|libcxx|libcxx-devel|libsigtool|"
    r"sigtool-codesign|tapi|llvm-tools|llvm-tools-19|libclang-cpp[0-9.]*|libllvm[0-9]*|"
    r"libasprintf|libgettextpo|gfortran)=")
keep = [s for s in specs if "osx-arm64" not in s and not DROP.match(s)]
dropped = [s for s in specs if s not in keep]
# ADD: conda-forge splits expat into `expat` (headers + expat.pc) and `libexpat` (runtime .so
# only). The lock carries ONLY libexpat, so fontconfig.pc's `Requires.private: ... expat`
# cannot resolve, `pkg-config --cflags fontconfig` fails outright, systemfonts' configure gets
# an EMPTY PKG_CFLAGS, and its compile dies on `ft2build.h: No such file or directory`.
# svglite depends on systemfonts and CellChat depends on svglite, so the whole install ends
# with "ERROR: dependency 'svglite' is not available for package 'CellChat'" -- which is
# exactly what happened on the first real build, 2026-08-18. See DEVIATIONS.md D-8.
added = ["expat"]
keep += added
pathlib.Path(out).write_text(
    "# Generated by cellchat/build_env_iris.sh from env.lock.yml -- DO NOT EDIT BY HAND.\n"
    "# %d of %d dependencies kept (+%d added); %d macOS/clang-only entries dropped "
    "(DEVIATIONS D-1, D-2, D-8):\n" % (len(keep) - len(added), len(specs), len(added), len(dropped))
    + "".join("#   ADDED %s\n" % a for a in added)
    + "".join("#   %s\n" % d for d in dropped)
    + "\n".join(keep) + "\n")
print(f"  {len(specs)} dependencies -> {len(keep)} kept, {len(dropped)} dropped")
PY
echo "  wrote $SPECS"

# ---------------------------------------------------------------- 2. solve, then create
if [ -d "$ENV_PREFIX" ]; then
    say "2. conda env already exists -- skipping create"
else
    say "2a. dry-run solve (writes nothing; abort here rather than half-build)"
    CONDA_SOLVER=classic conda create --dry-run -n "$ENV_NAME" \
        -c conda-forge --override-channels --file "$SPECS" > /tmp/cc_solve.$$ 2>&1 \
        || { tail -30 /tmp/cc_solve.$$; die "the solve failed -- see above"; }
    BAD=$(grep -icE "osx-arm64|clang|cctools|ld64|tapi|sigtool" /tmp/cc_solve.$$ || true)
    [ "$BAD" -eq 0 ] || {
        grep -iE "osx-arm64|clang|cctools|ld64|tapi|sigtool" /tmp/cc_solve.$$ | head
        die "the solve pulled macOS toolchain packages -- DEVIATIONS D-3 hazard, refusing"
    }
    echo "  solve OK, zero macOS toolchain packages"
    grep -E "gcc_linux-64|gxx_linux-64|gfortran_linux-64|^  r-base" /tmp/cc_solve.$$ | sed 's/^/  /'
    rm -f /tmp/cc_solve.$$

    say "2b. creating $ENV_NAME (this takes a few minutes)"
    CONDA_SOLVER=classic conda create -y -n "$ENV_NAME" \
        -c conda-forge --override-channels --file "$SPECS"
fi
[ -x "$ENV_PREFIX/bin/Rscript" ] || die "$ENV_PREFIX/bin/Rscript missing after create"

# ---------------------------------------------------------------- 3. activate
say "3. activating (plain \`conda activate\` gives you system R -- never use it here)"
# shellcheck source=/dev/null
source "$HERE/activate_env.sh"
echo "  R      : $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
echo "  R_LIBS : $R_LIBS_USER"

# ---------------------------------------------------------------- 4. toolchain assertion
# DEVIATIONS D-6 says to check this BEFORE starting the package installs. A wrong CC does not
# fail loudly -- it fails ~250 times, hours in.
say "4. toolchain assertion"
CC_SEEN=$(Rscript -e 'cat(Sys.getenv("CC"))')
CXX_SEEN=$(Rscript -e 'cat(Sys.getenv("CXX"))')
FC_SEEN=$(Rscript -e 'cat(Sys.getenv("FC"))')
echo "  CC =$CC_SEEN"
echo "  CXX=$CXX_SEEN"
echo "  FC =$FC_SEEN"
for v in "$CC_SEEN" "$CXX_SEEN" "$FC_SEEN"; do
    case "$v" in
        *linux-gnu*) ;;
        *) die "toolchain is not the linux one ('$v'). Refusing to start ~250 source builds." ;;
    esac
done
printf 'int main(void){return 0;}\n' > /tmp/cc_probe.$$.c
$CC_SEEN /tmp/cc_probe.$$.c -o /tmp/cc_probe.$$ 2>/dev/null || die "$CC_SEEN cannot compile a trivial C file"
file /tmp/cc_probe.$$ | grep -q ELF || die "$CC_SEEN did not produce an ELF binary (Mach-O? see D-3)"
rm -f /tmp/cc_probe.$$.c /tmp/cc_probe.$$
echo "  compiles a trivial C file to ELF  OK"

# The second thing that fails silently and expensively. systemfonts / textshaping / svglite all
# discover their headers through pkg-config; if ANY of these four cannot resolve, their configure
# scripts fall back to an empty PKG_CFLAGS and the compile dies on a missing header -- hours in,
# inside a try(), and the only symptom is CellChat quietly absent from the STILL MISSING line.
# Check it here, before the builds, the same way the toolchain is checked.
say "4b. pkg-config assertion (headers for systemfonts / textshaping / svglite)"
for pc in fontconfig freetype2 harfbuzz fribidi; do
    OUT=$(pkg-config --cflags "$pc" 2>&1) || die "pkg-config --cflags $pc FAILED: $OUT"
    [ -n "$OUT" ] || die "pkg-config --cflags $pc returned EMPTY -- systemfonts will not build"
    printf '  %-11s %s\n' "$pc" "$OUT"
done

# ---------------------------------------------------------------- 5. CellChat source
say "5. CellChat source at the recorded pin"
if [ -d "$CC_SRC/.git" ]; then
    echo "  clone exists, checking it out at $CC_SHA"
    git -C "$CC_SRC" fetch --quiet origin
else
    mkdir -p "$SRC_ROOT"
    git clone --quiet "$CC_REPO" "$CC_SRC"
fi
git -C "$CC_SRC" checkout --quiet "$CC_SHA"
GOT_SHA=$(git -C "$CC_SRC" rev-parse HEAD)
[ "$GOT_SHA" = "$CC_SHA" ] || die "checked out $GOT_SHA, expected $CC_SHA"
GOT_VER=$(awk '/^Version:/{print $2}' "$CC_SRC/DESCRIPTION")
[ "$GOT_VER" = "$CC_VERSION" ] || die "DESCRIPTION says Version: $GOT_VER, expected $CC_VERSION"
echo "  $CC_SRC @ ${GOT_SHA:0:8}  Version: $GOT_VER  OK"

# ---------------------------------------------------------------- 6. the R side
say "6. install_env.R -- ~250 R packages from $SNAP (THE LONG PART)"
echo "  parallel installs: $NCPU (CELLCHAT_NCPUS overrides; the script's own default is"
echo "  detectCores()-2 = $(nproc --all) - 2, which ignores the cgroup and would oversubscribe)"
export CELLCHAT_NCPUS="$NCPU"
Rscript "$HERE/install_env.R" --cellchat-src "$CC_SRC"

# ---------------------------------------------------------------- 7. verify
say "7. verification against r_packages.lock.csv"
Rscript - "$HERE/r_packages.lock.csv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
lock <- read.csv(args[1], stringsAsFactors = FALSE)
inst <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
m <- merge(lock, inst, by = "Package", all.x = TRUE, suffixes = c(".lock", ".have"))
missing <- m[is.na(m$Version.have), "Package"]
differ  <- m[!is.na(m$Version.have) & m$Version.lock != m$Version.have, ]
cat(sprintf("  locked %d | installed %d | MISSING %d | version-differs %d\n",
            nrow(lock), sum(!is.na(m$Version.have)), length(missing), nrow(differ)))
if (length(missing)) cat("  MISSING:", paste(head(missing, 40), collapse = " "), "\n")
if (nrow(differ)) {
  cat("  DIFFERS (lock -> installed):\n")
  print(head(differ[, c("Package", "Version.lock", "Version.have")], 40), row.names = FALSE)
}
ok <- requireNamespace("CellChat", quietly = TRUE)
cat("  CellChat loadable:", ok, "\n")
if (!ok) {
  cat("\n  FATAL: CellChat is NOT installed. install_env.R wraps every install in try(), so the\n")
  cat("         real error is further up this log -- grep it for 'ERROR: '. The classic one is\n")
  cat("         a missing header from a failed pkg-config lookup (see DEVIATIONS.md D-8).\n")
  quit(status = 1)
}
if (ok) {
  suppressMessages(library(CellChat))
  cat("  CellChat", as.character(packageVersion("CellChat")), "\n")
  cat("  CellChatDB.human$interaction rows:", nrow(CellChatDB.human$interaction), "\n")
}
RS
[ $? -eq 0 ] || die "verification failed -- the env is NOT usable (see above)"

say "DONE"
cat <<EOF

  env    : $ENV_PREFIX
  source : $CC_SRC @ ${CC_SHA:0:8}  (Version $CC_VERSION)
  specs  : $SPECS

Next, confirm the runner's config self-check now reports RS_CELLCHAT OK:

  bash scripts/comparators/_common/luad_config.sh

Then stage 04 becomes runnable. Record in cellchat/DEVIATIONS.md what this build did NOT
reproduce: the lock's osx-arm64 build strings for low-level libraries. Everything that carries
scientific meaning -- R 4.3.3, the CRAN snapshot date, the 250 R package versions, the CellChat
commit -- is pinned.
EOF
