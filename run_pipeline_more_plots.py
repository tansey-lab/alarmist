#!/usr/bin/env python3
"""
ALARMIST pipeline runner — general, dataset-agnostic, resumable.

A drop-in successor to run_pipeline.py that keeps that script's "one command runs
the whole thing" convenience but adds the workflow features of the es_xenium
driver: an explicit stage list with presets, resume (--from / --to / sentinels),
--dry-run, --force, reuse of an existing patchify, a shared cell-type colour code,
two interleaved figure families, and a provenance manifest.

NOTHING is hard-coded. Every path, the interpreter, and the ligand-receptor
database are taken from the command line; the only defaults are ALARMIST's own
(cell_type column name, patch size, etc.). Point it at any h5ad. Any
dataset-specific preprocessing (e.g. collapsing cell-type subtypes) is done
beforehand — this runner consumes --data-file as-is.

Stages, in order (figures interleaved with compute, so each family of plots is
drawn as soon as its inputs exist):

    patchify        patch x LRI count matrix                  [alarmist-patchify]
    bptf            BPTF factorization -> motifs              [alarmist-bptf]
  > plots-bptf      motif / LRI figures (both families)
    project         per-cell neighbourhood projection         [alarmist-project]
  > plots-project   composition, GMM states, spatial maps
    markers         per-cell-type marker exclusion mask
    glm             motif -> gene impact (Poisson GLM)        [alarmist-glm]
  > plots-glm       volcano / forest (both families)

Two figure families, same run directory, one shared colour code:

    <run>/plots_original/   the GBM.ipynb suite ("original plots")
    <run>/plots/            stock alarmist.cli.visualize ("workflow plots")

Both resolve every cell type to the SAME colour, derived once from the input
h5ad's cell-type categories and cached in <run>/celltype_colors.json.

Compute stages shell out to `<python> -m alarmist.cli.<stage>` (so they run in
whichever interpreter you choose); figure and marker stages run in-process (the
colour-code registration and the stock CLI's in-process call are the whole point).
Run this with your ALARMIST env's python so both halves share one environment.

Examples
--------
    PY=/path/to/env/bin/python

    # see the plan, run nothing (recommended first)
    $PY run_pipeline_more_plots.py --data-file data.h5ad --preset full --dry-run

    # patchify -> bptf + the figures that need nothing else
    $PY run_pipeline_more_plots.py --data-file data.h5ad --n-components 20 --preset bptf

    # a full run with a custom (e.g. mouse) LR database
    $PY run_pipeline_more_plots.py --data-file data.h5ad \
        --cellchatdb data/LRdatabase/CellChatDBv2.0.mouse.csv

    # reuse a deterministic patchify from another run
    $PY run_pipeline_more_plots.py --data-file data.h5ad \
        --patchify-dir results/run_seed0_15/patchify

    # resume from a stage (finished stages skip on their sentinel output)
    $PY run_pipeline_more_plots.py --data-file data.h5ad --from glm

    # redraw an existing run's figures with an LRI-network threshold sweep
    $PY run_pipeline_more_plots.py --data-file data.h5ad --preset plots --force \
        --network-threshold 1 50 100 200
"""

import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ALL_STAGES = [
    "patchify",
    "bptf",
    "plots-bptf",
    "project",
    "plots-project",
    "markers",
    "glm",
    "plots-glm",
]

PRESETS = {
    "full": ALL_STAGES,
    "bptf": ["patchify", "bptf", "plots-bptf"],
    "plots": ["plots-bptf", "plots-project", "plots-glm"],
}

# The stage whose figures a plots-* stage draws.
PLOTS_STAGE = {"plots-bptf": "bptf", "plots-project": "project", "plots-glm": "glm"}

PALETTE = "tab20"  # stock ALARMIST cell-type palette; do not change casually.


def utcnow() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ===========================================================================
# logging: tee everything (this driver's prints, in-process figure output, and
# subprocess output) to <run>/run_<ts>.log
# ===========================================================================
class _Tee:
    def __init__(self, *streams):
        self.streams = [s for s in streams if s is not None]

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()

    def flush(self):
        for st in self.streams:
            st.flush()


# ===========================================================================
# the cell-type colour code (ported from the es_xenium figure engine, generalised)
# ===========================================================================
def sanitize(name) -> str:
    """The GLM's cell-type key transform: 'natural killer cell' -> 'natural_killer_cell'."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_")


def save_all_formats(fig, path_no_ext, dpi=300):
    """png + pdf + svg via matplotlib (no external dependency)."""
    out = Path(path_no_ext)
    out.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "pdf", "svg"):
        p = f"{path_no_ext}.{ext}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        paths.append(p)
    return paths


def read_celltype_categories(h5ad_path, column):
    """Cell-type levels in category order — handles both anndata categorical encodings.

    Modern (>=0.7): obs/<col> is a group holding a 'categories' array + 'codes'.
    Legacy (<0.7): obs/<col> is an int-code dataset and the string labels live in
    obs/__categories/<col> (this is what the MOSTA cell_bin files use — reading the
    codes directly is what produced numeric '0','1',... colour keys).
    """
    import h5py
    import numpy as np

    def _dec(arr):
        return [
            c.decode() if isinstance(c, (bytes, np.bytes_)) else str(c) for c in arr
        ]

    with h5py.File(h5ad_path, "r") as f:
        key = f"obs/{column}"
        if key not in f:
            raise KeyError(f"'{key}' not found in {h5ad_path}")
        node = f[key]
        # modern: group with an explicit categories array
        if isinstance(node, h5py.Group) and "categories" in node:
            return _dec(node["categories"][()])
        # legacy: int codes in obs/<col>, labels in obs/__categories/<col>
        legacy = f.get(f"obs/__categories/{column}")
        if legacy is not None:
            return _dec(legacy[()])
        # plain string column: sorted unique values (what .astype('category') gives)
        return sorted(set(_dec(np.asarray(node[()]).ravel())))


def build_colors(categories):
    """The stock recipe: tab20 resampled to n_types, indexed in category order."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex

    cmap = plt.get_cmap(PALETTE, len(categories))
    return {ct: to_hex(cmap(i)) for i, ct in enumerate(categories)}


def with_sanitized_aliases(colors):
    """Add each cell type's underscored (GLM) key, pointing at the same colour."""
    out = dict(colors)
    for ct, colour in colors.items():
        out.setdefault(sanitize(ct), colour)
    return out


def load_or_build_colors(run_dir, adata_path, column, refresh=False):
    """Return (colors, categories), cached in <run>/celltype_colors.json."""
    cache = Path(run_dir) / "celltype_colors.json"
    if cache.exists() and not refresh:
        meta = json.loads(cache.read_text())
        if meta.get("cell_type_column") != column:
            raise ValueError(
                f"{cache} was built for cell_type_column='{meta.get('cell_type_column')}' "
                f"but --cell-type-column='{column}'. Delete it or pass --refresh-colors."
            )
        print(
            f"colour code <- {cache} ({len(meta['categories'])} cell types)", flush=True
        )
        return meta["colors"], meta["categories"]

    if not adata_path or not Path(adata_path).exists():
        raise FileNotFoundError(
            f"Need --data-file to derive the colour code; {cache} does not exist yet."
        )
    categories = read_celltype_categories(adata_path, column)
    colors = build_colors(categories)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(
        json.dumps(
            {
                "palette": PALETTE,
                "cell_type_column": column,
                "source_h5ad": str(adata_path),
                "recipe": "plt.get_cmap('tab20', n_types)(i) over adata.obs[col].cat.categories",
                "categories": categories,
                "colors": colors,
            },
            indent=2,
        )
    )
    print(f"colour code -> {cache} ({len(categories)} cell types)", flush=True)
    return colors, categories


def plot_color_key(colors, categories, out_no_ext):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(3.4, 0.32 * len(categories) + 0.6))
    ax.legend(
        handles=[
            Patch(facecolor=colors[c], edgecolor="none", label=c) for c in categories
        ],
        loc="center",
        frameon=False,
        handlelength=1.2,
        fontsize=9,
    )
    ax.axis("off")
    ax.set_title(f"cell-type colour code ({PALETTE}, n={len(categories)})", fontsize=10)
    save_all_formats(fig, out_no_ext)
    plt.close(fig)


# ===========================================================================
# a per-figure runner: one figure's failure must not sink the whole stage
# ===========================================================================
class Runner:
    def __init__(self):
        self.ok, self.fail = [], []

    def step(self, name, fn):
        import traceback

        import matplotlib.pyplot as plt

        print(f"\n=== {name} ===", flush=True)
        try:
            fn()
            self.ok.append(name)
            print(f"[OK] {name}", flush=True)
        except Exception as exc:  # noqa: BLE001 - deliberately non-fatal
            self.fail.append(name)
            print(f"[FAIL] {name}: {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
        finally:
            plt.close("all")


# ===========================================================================
# "original plots" family (the GBM.ipynb suite) — in-process
# ===========================================================================
def load_lri_motifs(run_dir, cellchatdb):
    import pandas as pd

    import alarmist as al

    bptf = al.load_bptf_results(f"{run_dir}/bptf")
    lri_motifs = bptf["lri_motifs"]
    fcol = "factor_lrnorm" if "factor_lrnorm" in lri_motifs.columns else "factor"
    print(f"lri_motifs {lri_motifs.shape}, factor_col={fcol}", flush=True)

    lri_pw = lri_motifs
    if cellchatdb and Path(cellchatdb).exists():
        # Re-annotate from the DB the user actually chose: alarmist-bptf annotates
        # against a hard-coded default, which is wrong for a non-default species.
        try:
            from alarmist.plotting.motif_plots import annotate_pathways

            base = lri_motifs.drop(columns=["pathway"], errors="ignore").copy()
            lri_pw = annotate_pathways(base, pd.read_csv(cellchatdb))
            n = int(lri_pw["pathway"].notna().sum()) if "pathway" in lri_pw else 0
            print(f"annotate_pathways OK: {n}/{len(lri_pw)} annotated", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(
                f"annotate_pathways failed ({exc}); pathway plot may be skipped",
                flush=True,
            )
            lri_pw = lri_motifs
    return lri_motifs, lri_pw, fcol


def original_stage_bptf(a, out, colors):
    import matplotlib.pyplot as plt

    from alarmist.plotting import (
        plot_celltype_communication_by_motif,
        plot_lri_networks,
        plot_lri_networks_html,
        plot_top_lri_interactions_by_pathway,
        plot_top_lri_interactions_dot,
    )

    r = Runner()
    lri_motifs, lri_pw, fcol = load_lri_motifs(a.output_dir, a.cellchatdb)

    r.step(
        "celltype_communication_by_motif",
        lambda: plot_celltype_communication_by_motif(
            lri_motifs,
            factor_col=fcol,
            n_cols=5,
            save_path=f"{out}/celltype_communication.pdf",
        ),
    )
    r.step(
        "top_lri_interactions_dot",
        lambda: plot_top_lri_interactions_dot(
            lri_motifs,
            factor_col=fcol,
            top_n=a.lri_top_n,
            ct_colors=colors,
            save_path=f"{out}/top_lri_interactions_dot.pdf",
        ),
    )
    r.step(
        "top_lri_interactions_by_pathway",
        lambda: plot_top_lri_interactions_by_pathway(
            lri_pw,
            factor_col=fcol,
            top_n=a.pathway_top_n,
            ct_colors=colors,
            save_path=f"{out}/top_lri_interactions_by_pathway.pdf",
        ),
    )

    def _networks():
        for i, thr in enumerate(a.network_threshold):
            path = f"{out}/lri_networks_thr{thr:g}.svg"
            plot_lri_networks(
                lri_motifs,
                top_n=a.network_top_n,
                threshold=thr,
                factor_col=fcol,
                ct_colors=colors,
                save_path=path,
            )
            plt.close("all")
            print(f"  threshold {thr:g} -> {path}", flush=True)
            if i == 0:
                shutil.copyfile(path, f"{out}/lri_networks.svg")

    r.step(f"lri_networks (thresholds {a.network_threshold})", _networks)

    r.step(
        "lri_networks_interactive_html",
        lambda: plot_lri_networks_html(
            lri_motifs,
            f"{out}/lri_networks_interactive.html",
            top_n=a.network_top_n,
            factor_col=fcol,
            ct_colors=colors,
        ),
    )
    return r


def original_stage_project(a, out, colors, adata, cell_loadings):
    import matplotlib.pyplot as plt
    import numpy as np

    import alarmist as al
    from alarmist.plotting import plot_motif_spatial

    r = Runner()
    n_motifs = cell_loadings.shape[1]
    motifs = list(range(n_motifs)) if a.motifs is None else list(a.motifs)

    def _composition():
        fig, _ax, df = al.analyze_motif_celltype_composition(
            adata,
            cell_loadings,
            cell_type_column=a.cell_type_column,
            ct_colors=colors,
            output_dir=out,
        )
        df.to_csv(f"{out}/motif_celltype_weighted.csv")
        save_all_formats(fig, f"{out}/motif_celltype_weighted")

    r.step("analyze_motif_celltype_composition", _composition)

    def _gmm():
        np.random.seed(a.seed)
        summary = al.gmm_binarize_all_motifs(cell_loadings, adata, random_state=a.seed)
        try:
            summary.to_csv(f"{a.output_dir}/gmm_summary.csv", index=False)
        except Exception as exc:  # noqa: BLE001
            print(f"  (gmm_summary not written: {exc})", flush=True)
        cols = [f"motif_{k}_state" for k in range(n_motifs)]
        states = adata.obs[cols].astype(str).copy()
        states.index = adata.obs_names.astype(str)
        states.index.name = "cell_id"
        states.to_parquet(f"{a.output_dir}/motif_states.parquet")
        print(f"  motif_states.parquet {states.shape}", flush=True)

    r.step("gmm_binarize_all_motifs -> motif_states.parquet", _gmm)

    def _state_counts():
        fig, _ax, df = al.analyze_motif_state_counts(adata, output_dir=out)
        df.to_csv(f"{out}/motif_state_counts.csv")
        save_all_formats(fig, f"{out}/motif_state_counts")

    r.step("analyze_motif_state_counts", _state_counts)

    def _spatial():
        sp = f"{out}/motif_spatial"
        os.makedirs(sp, exist_ok=True)
        for m in motifs:
            try:
                plot_motif_spatial(
                    adata,
                    motif_idx=int(m),
                    point_size=a.point_size,
                    cell_type_column=a.cell_type_column,
                    ct_colors=colors,
                    output_dir=sp,
                )
                plt.close("all")
                print(f"  motif {m} spatial OK", flush=True)
            except Exception as exc:  # noqa: BLE001
                print(
                    f"  motif {m} spatial FAIL: {type(exc).__name__}: {exc}", flush=True
                )

    r.step(f"plot_motif_spatial ({len(motifs)} motifs)", _spatial)
    return r


def original_stage_glm(a, out, adata):
    import numpy as np

    import alarmist as al

    r = Runner()
    mask_csv = f"{a.output_dir}/markers/exclusion_matrix.csv"
    if not os.path.exists(mask_csv):
        raise FileNotFoundError(f"{mask_csv} missing — run the 'markers' stage first.")
    glm = al.load_glm_results(f"{a.output_dir}/glm")
    cts, genes, excl = al.load_exclusion_mask(mask_csv)
    cts_san = np.array([sanitize(c) for c in cts])

    original = adata.obs[a.cell_type_column].copy()
    adata.obs[a.cell_type_column] = (
        original.astype(str).map(sanitize).astype("category")
    )
    try:
        r.step(
            "glm_volcano (all cell types)",
            lambda: al.glm_volcano(
                adata=adata,
                de_results=glm,
                cell_types=cts_san,
                all_genes=genes,
                exclusion_mask=excl,
                fdr_threshold=a.fdr,
                lfc_threshold=a.lfc,
                min_expression_frac=a.min_expression_frac,
                output_dir=out,
            ),
        )
        r.step(
            "glm_forest (all cell types)",
            lambda: al.glm_forest(
                adata=adata,
                de_results=glm,
                cell_types=cts_san,
                all_genes=genes,
                exclusion_mask=excl,
                min_expression_frac=a.min_expression_frac,
                output_dir=out,
            ),
        )
    finally:
        adata.obs[a.cell_type_column] = original
    return r


# ===========================================================================
# "workflow plots" family — the stock CLI, run in-process so it inherits palette
# ===========================================================================
VISUALIZE_PLOT_TYPES = {
    "bptf": ["heatmap", "motif_summary", "lri_dot", "lri_network"],
    "project": ["spatial"],
    "glm": ["volcano", "forest"],
}


def run_visualize(a, fig_stage):
    import matplotlib.pyplot as plt

    from alarmist.cli import visualize as viz_cli

    plot_types = VISUALIZE_PLOT_TYPES[fig_stage]
    out = f"{a.output_dir}/plots"
    os.makedirs(out, exist_ok=True)
    argv = [
        "alarmist-visualize",
        "--glm-dir",
        f"{a.output_dir}/glm",
        "--bptf-dir",
        f"{a.output_dir}/bptf",
        "--patchify-dir",
        f"{a.output_dir}/patchify",
        "--output-dir",
        out,
        "--plot-types",
        *plot_types,
        "--format",
        a.format,
        "--cell-type-column",
        a.cell_type_column,
        "--network-top-n",
        str(a.network_top_n),
        "--network-threshold",
        str(a.viz_network_threshold),
        "--min-expression-frac",
        str(a.min_expression_frac),
        "--alpha",
        str(a.fdr),
        "-v",
    ]
    if {"spatial", "volcano", "forest"} & set(plot_types):
        argv += ["--project-dir", f"{a.output_dir}/project"]
    mask = f"{a.output_dir}/markers/exclusion_matrix.csv"
    if os.path.exists(mask):
        argv += ["--exclusion-mask", mask]

    print(f"\n=== stock visualize: {' '.join(plot_types)} ===", flush=True)
    print("  " + " ".join(argv), flush=True)
    saved = sys.argv
    try:
        sys.argv = argv
        rc = viz_cli.main() or 0
    finally:
        sys.argv = saved
        plt.close("all")
    print(f"outputs -> {out}", flush=True)
    return rc


def draw_figures(a, fig_stage, colors):
    """Run the requested figure families for one stage (bptf|project|glm).

    `colors` is the CLEAN palette (one key per real cell type). The GLM stores cell
    types with underscores, so volcano/forest need underscored aliases in the global
    registry to find each colour — but those extra keys would also duplicate every
    multi-word cell type in any legend built by looping the colour dict. So we add
    the aliases to the registry (and pass them) ONLY around the glm figures, then
    restore the clean registry.
    """
    import alarmist as al

    families = {
        "both": ("original", "visualize"),
        "original": ("original",),
        "visualize": ("visualize",),
    }[a.figures]

    if fig_stage == "glm":
        al.set_celltype_colors(with_sanitized_aliases(colors))
    try:
        if "original" in families:
            out = f"{a.output_dir}/plots_original"
            os.makedirs(out, exist_ok=True)
            if fig_stage == "bptf":
                r = original_stage_bptf(a, out, colors)
            else:
                import anndata
                import numpy as np

                print("loading projected_adata.h5ad ...", flush=True)
                adata = anndata.read_h5ad(
                    f"{a.output_dir}/project/projected_adata.h5ad"
                )
                if fig_stage == "project":
                    cl = np.load(f"{a.output_dir}/project/cell_loadings.npy")
                    r = original_stage_project(a, out, colors, adata, cl)
                else:  # glm
                    r = original_stage_glm(a, out, adata)
            print(
                f"\noriginal plots [{fig_stage}]: OK={len(r.ok)} FAIL={r.fail} -> {out}",
                flush=True,
            )

        if "visualize" in families:
            run_visualize(a, fig_stage)
    finally:
        if fig_stage == "glm":
            al.set_celltype_colors(colors)


# ===========================================================================
# compute stages: shell out to the alarmist CLI
# ===========================================================================
def run_cmd(cmd, dry):
    print(f"\n>>>> {utcnow()}  {' '.join(cmd)}", flush=True)
    if dry:
        print("     (dry-run: not executed)", flush=True)
        return
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    for line in proc.stdout:
        sys.stdout.write(line)
    proc.wait()
    if proc.returncode != 0:
        raise SystemExit(f"!! command failed (exit {proc.returncode}): {' '.join(cmd)}")


def cli(a, module, *stage_args):
    return [a.python, "-m", f"alarmist.cli.{module}", *stage_args]


def save_cells_per_patch(a):
    """Save the cells-per-patch histogram (al.plot_cells_per_patch) into the run dir.

    patchify runs as a subprocess, so obs['patch_idx'] is not persisted; we
    re-derive the identical grid assignment here (it depends only on
    obsm['spatial'] and patch_size) and hand it to the stock plot. Best-effort:
    a QC plot must never fail the run.
    """
    if a.dry_run:
        print("     (dry-run) cells-per-patch histogram", flush=True)
        return
    try:
        import anndata
        import matplotlib

        import alarmist as al

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        adata = anndata.read_h5ad(a.data_file)
        analyzer = al.PatchLRIAnalyzer(
            patch_size=float(a.patch_size), cell_type_column=a.cell_type_column
        )
        patch_assignments, _ = analyzer.create_spatial_patches(adata)
        adata.obs["patch_idx"] = patch_assignments
        out = f"{a.output_dir}/cells_per_patch.png"
        al.plot_cells_per_patch(adata, save_path=out, show=False)
        plt.close("all")
        print(f">>>> cells-per-patch histogram -> {out}", flush=True)
    except Exception as exc:  # noqa: BLE001 - a QC plot must not sink the run
        print(
            f">>>> cells-per-patch plot skipped: {type(exc).__name__}: {exc}",
            flush=True,
        )


def stage_patchify(a):
    args = [
        "--adata",
        a.data_file,
        "--output-dir",
        f"{a.output_dir}/patchify",
        "--cell-type-column",
        a.cell_type_column,
        "--resource",
        a.resource,
        "--patch-size",
        str(a.patch_size),
        "-v",
    ]
    if a.cellchatdb:
        args += ["--cellchatdb-path", a.cellchatdb]
    run_cmd(cli(a, "patchify", *args), a.dry_run)
    save_cells_per_patch(a)


def reuse_patchify(a):
    """Symlink an existing patchify into this run after verifying its geometry."""
    abs_src = Path(a.patchify_dir).resolve()
    params = abs_src / "analysis_parameters.csv"
    if params.exists():
        got = {}
        for line in params.read_text().splitlines():
            parts = line.split(",")
            if len(parts) >= 2:
                got[parts[0].strip()] = parts[1].strip()
        ps, rs = got.get("patch_size"), got.get("resource_name")
        if ps and abs(float(ps) - float(a.patch_size)) > 1e-4:
            raise SystemExit(
                f"!! reused patchify patch_size={ps} != --patch-size={a.patch_size}"
            )
        if rs and rs != a.resource:
            raise SystemExit(
                f"!! reused patchify resource={rs} != --resource={a.resource}"
            )
        print(
            f">>>> reuse patchify: {abs_src} (patch_size={ps} resource={rs})",
            flush=True,
        )
    else:
        print(
            f">>>> reuse patchify: {abs_src} (no analysis_parameters.csv — NOT verified)",
            flush=True,
        )
    link = Path(a.output_dir) / "patchify"
    if not link.exists():
        if a.dry_run:
            print(f"     (dry-run) ln -s {abs_src} {link}", flush=True)
        else:
            os.symlink(abs_src, link)
    save_cells_per_patch(a)


def stage_bptf(a):
    # Run BPTF IN-PROCESS, not via `alarmist-bptf`. The CLI calls
    # process_bptf_results() without a cellchatdb_path, so it annotates pathways
    # against its hard-coded HUMAN default. annotate_pathways has an
    # O(rows x pairs) fuzzy-match fallback, so on non-human data (where almost no
    # LRI matches the human table) EVERY row hits the slow path and the stage can
    # hang for many minutes producing useless "Unknown" pathways. Calling the
    # library directly lets us pass the DB the run actually uses -> exact matches,
    # fast, correct. (The library function takes cellchatdb_path; the CLI just
    # never forwards it.)
    patchify_dir = f"{a.output_dir}/patchify"
    bptf_dir = f"{a.output_dir}/bptf"
    print(
        f"\n>>>> {utcnow()}  [bptf] in-process  K={a.n_components}  max_iter={a.max_iter}  "
        f"seed={a.seed}  cellchatdb={a.cellchatdb}",
        flush=True,
    )
    if a.dry_run:
        print(
            f"     (dry-run) al.run_bptf(load_patch_lri_results('{patchify_dir}'), "
            f"n_components={a.n_components}, max_iter={a.max_iter}, random_state={a.seed}) "
            f"-> process_bptf_results(output_dir='{bptf_dir}', cellchatdb_path={a.cellchatdb!r})",
            flush=True,
        )
        return
    import alarmist as al

    results = al.load_patch_lri_results(patchify_dir)
    if "patch_lri_matrix" not in results:
        raise SystemExit(
            f"load_patch_lri_results returned no 'patch_lri_matrix' (keys: {list(results)})"
        )
    model = al.run_bptf(
        results["patch_lri_matrix"],
        n_components=a.n_components,
        max_iter=a.max_iter,
        verbose=True,
        random_state=a.seed,
    )
    kwargs = {"cellchatdb_path": a.cellchatdb} if a.cellchatdb else {}
    al.process_bptf_results(model, results, output_dir=bptf_dir, **kwargs)
    print(f">>>> bptf -> {bptf_dir}", flush=True)


def stage_project(a):
    args = [
        "--adata",
        a.data_file,
        "--bptf-dir",
        f"{a.output_dir}/bptf",
        "--patch-lri-dir",
        f"{a.output_dir}/patchify",
        "--output-dir",
        f"{a.output_dir}/project",
        "--resource",
        a.resource,
        "--cell-type-column",
        a.cell_type_column,
        "-v",
    ]
    if a.cellchatdb:
        args += ["--cellchatdb", a.cellchatdb]
    run_cmd(cli(a, "project", *args), a.dry_run)


def stage_markers(a):
    print(
        f"\n>>>> {utcnow()}  [markers] compute per-cell-type exclusion mask", flush=True
    )
    if a.dry_run:
        print("     (dry-run: not executed)", flush=True)
        return
    import anndata
    import numpy as np

    import alarmist as al

    np.random.seed(a.seed)
    adata = anndata.read_h5ad(f"{a.output_dir}/project/projected_adata.h5ad")
    al.compute_exclusion_mask(
        adata,
        marker_lfc=a.marker_lfc,
        marker_pvalue=a.marker_pvalue,
        marker_subsample=a.marker_subsample,
        output_dir=f"{a.output_dir}/markers",
    )
    print(f"exclusion mask -> {a.output_dir}/markers", flush=True)


def stage_glm(a):
    args = [
        "--input-dir",
        f"{a.output_dir}/project",
        "--adata",
        a.data_file,
        "--output-dir",
        f"{a.output_dir}/glm",
        "--patch-lri-dir",
        f"{a.output_dir}/patchify",
        "--cell-type-column",
        a.cell_type_column,
        "--count-layer",
        a.count_layer,
        "--backend",
        a.glm_backend,
        "--device",
        a.glm_device,
        "--gene-tile",
        str(a.gene_tile),
        "--glm-dtype",
        a.glm_dtype,
        "--alpha",
        str(a.fdr),
        "--seed",
        str(a.seed),
        ("--prefilter-spearman" if a.prefilter_spearman else "--no-prefilter-spearman"),
        "-v",
    ]
    run_cmd(cli(a, "glm", *args), a.dry_run)


# ===========================================================================
# stage selection + orchestration
# ===========================================================================
def sentinel(a, stage):
    o = a.output_dir
    return {
        "patchify": f"{o}/patchify/patch_lri_matrix.npz",
        "bptf": f"{o}/bptf/lri_motifs.csv",
        "plots-bptf": f"{o}/plots_original/lri_networks.svg",
        "project": f"{o}/project/cell_loadings.npy",
        "plots-project": f"{o}/motif_states.parquet",
        "markers": f"{o}/markers/exclusion_matrix.csv",
        "glm": f"{o}/glm",
        "plots-glm": f"{o}/plots_original/volcano_plots_filtered.pdf",
    }[stage]


def select_stages(a):
    if a.preset:
        return list(PRESETS[a.preset])
    if a.stages:
        req = [s.strip() for s in re.split(r"[,\s]+", a.stages) if s.strip()]
        for s in req:
            if s not in ALL_STAGES:
                raise SystemExit(f"unknown stage: {s} (valid: {', '.join(ALL_STAGES)})")
        return req
    lo = ALL_STAGES.index(a.__dict__["from"]) if a.__dict__["from"] else 0
    hi = ALL_STAGES.index(a.to) if a.to else len(ALL_STAGES) - 1
    if lo > hi:
        raise SystemExit("--from is after --to")
    return ALL_STAGES[lo : hi + 1]


def should_run(a, stage, selected):
    if stage not in selected:
        return False
    if a.figures == "none" and stage.startswith("plots-"):
        return False
    if not a.force and os.path.exists(sentinel(a, stage)):
        print(
            f">>>> skip  {stage}  (exists: {sentinel(a, stage)}; pass --force to redo)",
            flush=True,
        )
        return False
    return True


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    p.add_argument("--data-file", required=True, help="Input h5ad (used as-is).")
    p.add_argument(
        "--output-dir",
        default=None,
        help="Run directory (default: results/run_seed<seed>_<ncomp>).",
    )
    p.add_argument(
        "--python",
        default=None,
        help="Interpreter for the CLI compute stages (default: this interpreter).",
    )

    # run parameters (defaults are ALARMIST's own, not dataset-specific)
    p.add_argument(
        "--n-components", type=int, default=15, help="Number of BPTF motifs."
    )
    p.add_argument(
        "--patch-size", type=float, default=50.0, help="Patch grid side (um)."
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-iter", type=int, default=1000, help="BPTF max iterations.")
    p.add_argument("--cell-type-column", default="cell_type")
    p.add_argument(
        "--resource", default="cellchatdb", choices=["cellchatdb", "cellphonedb"]
    )
    p.add_argument(
        "--cellchatdb",
        default=None,
        help="Path to a custom LR database CSV (e.g. a mouse CellChatDB). "
        "Passed to patchify/project and used to (re)annotate pathways. "
        "If omitted, ALARMIST's bundled default is used.",
    )
    p.add_argument(
        "--count-layer", default="X", help="GLM count source: X | raw | layers:NAME."
    )
    p.add_argument("--glm-backend", default="torch", choices=["torch", "sklearn"])
    p.add_argument(
        "--glm-device", default="cpu", help="torch device: cpu | auto | cuda | mps."
    )
    p.add_argument("--glm-dtype", default="float64", choices=["float64", "float32"])
    p.add_argument("--gene-tile", type=int, default=2048)
    p.add_argument(
        "--prefilter-spearman",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Spearman pre-filter before the GLMs (default: on).",
    )
    p.add_argument(
        "--patchify-dir",
        default=None,
        help="Reuse this existing patchify directory (skips the patchify stage).",
    )
    # marker stage
    p.add_argument("--marker-lfc", type=float, default=1.0)
    p.add_argument("--marker-pvalue", type=float, default=1e-5)
    p.add_argument("--marker-subsample", type=int, default=50000)

    # stage selection
    p.add_argument(
        "--preset",
        choices=list(PRESETS),
        default=None,
        help="full = every stage; bptf = patchify,bptf,plots-bptf; "
        "plots = plots-bptf,plots-project,plots-glm.",
    )
    p.add_argument("--stages", default=None, help="Explicit comma/space stage list.")
    p.add_argument("--from", dest="from", default=None, help="Start stage (inclusive).")
    p.add_argument("--to", default=None, help="End stage (inclusive).")
    p.add_argument(
        "--force", action="store_true", help="Re-run stages whose outputs exist."
    )
    p.add_argument(
        "--dry-run", action="store_true", help="Print the plan; run nothing."
    )

    # figures
    p.add_argument(
        "--figures",
        default="both",
        choices=["both", "original", "visualize", "none"],
        help="original = plots_original/ (GBM.ipynb suite); "
        "visualize = plots/ (stock CLI); both (default); none.",
    )
    p.add_argument(
        "--format",
        default="pdf",
        choices=["png", "pdf", "svg"],
        help="Format for the stock workflow figures.",
    )
    p.add_argument(
        "--network-threshold",
        nargs="+",
        type=float,
        default=[100.0],
        help="Original LRI-network edge-weight cutoff(s); pass several to sweep.",
    )
    p.add_argument(
        "--viz-network-threshold",
        type=float,
        default=0.0,
        help="Stock network figure cutoff (0 = alarmist's own behaviour).",
    )
    p.add_argument("--network-top-n", type=int, default=200)
    p.add_argument("--lri-top-n", type=int, default=35)
    p.add_argument("--pathway-top-n", type=int, default=40)
    p.add_argument(
        "--point-size", type=float, default=0.2, help="Spatial map marker size."
    )
    p.add_argument(
        "--motifs",
        nargs="+",
        type=int,
        default=None,
        help="Restrict per-motif spatial maps to these motif indices.",
    )
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=0.2)
    p.add_argument("--min-expression-frac", type=float, default=0.02)
    p.add_argument(
        "--refresh-colors",
        action="store_true",
        help="Rebuild celltype_colors.json (recolours the whole run).",
    )

    a = p.parse_args(argv)
    a.python = a.python or sys.executable
    if a.output_dir is None:
        a.output_dir = f"results/run_seed{a.seed}_{a.n_components}"
    return a


def write_manifest(a, selected):
    lines = [
        f"timestamp_utc = {utcnow()}",
        f"stages = {' '.join(selected)}",
        f"seed = {a.seed}   n_components = {a.n_components}   patch_size = {a.patch_size}   max_iter = {a.max_iter}",
        f"data_file = {a.data_file}",
        f"output_dir = {a.output_dir}",
        f"patchify = {a.patchify_dir or (a.output_dir + '/patchify (built here)')}",
        f"cell_type_column = {a.cell_type_column}   resource = {a.resource}   cellchatdb = {a.cellchatdb}",
        f"count_layer = {a.count_layer}   glm_backend = {a.glm_backend}   glm_device = {a.glm_device}   gene_tile = {a.gene_tile}",
        f"prefilter_spearman = {a.prefilter_spearman}   glm_dtype = {a.glm_dtype}",
        f"figures = {a.figures}   network_threshold = {a.network_threshold}   viz_network_threshold = {a.viz_network_threshold}",
        f"python = {a.python}",
    ]
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:  # noqa: BLE001
        commit = "unknown"
    lines.append(f"alarmist_git_commit = {commit}")
    text = "\n".join(lines) + "\n"
    print(text, flush=True)
    if not a.dry_run:
        (Path(a.output_dir) / "run_manifest.txt").write_text(text)


def main(argv=None):
    a = parse_args(argv)
    selected = select_stages(a)
    Path(a.output_dir).mkdir(parents=True, exist_ok=True)

    logf = None
    if not a.dry_run:
        logf = open(
            Path(a.output_dir)
            / f"run_{utcnow().replace(':', '').replace('-', '')}.log",
            "a",
        )
        sys.stdout = _Tee(sys.__stdout__, logf)
        sys.stderr = _Tee(sys.__stderr__, logf)

    print("=" * 78, flush=True)
    print(
        f"ALARMIST  seed={a.seed}  K={a.n_components}  patch={a.patch_size}um  ->  {a.output_dir}",
        flush=True,
    )
    print(f"  stages : {' '.join(selected)}", flush=True)
    print(f"  figures: {a.figures}   data-file: {a.data_file}", flush=True)
    if a.dry_run:
        print("  *** DRY RUN — nothing will be executed ***", flush=True)
    print("=" * 78, flush=True)

    write_manifest(a, selected)

    # Build the colour code up front whenever any figures will be drawn, so it is
    # derived once and shared by every stage (and the swatch key exists early).
    colors = None
    if (
        a.figures != "none"
        and not a.dry_run
        and any(s.startswith("plots-") for s in selected)
    ):
        import alarmist as al

        colors_orig, categories = load_or_build_colors(
            a.output_dir, a.data_file, a.cell_type_column, refresh=a.refresh_colors
        )
        # Register the CLEAN palette (one key per real cell type). The underscored
        # GLM aliases are added only around the glm figures (see draw_figures), so
        # bptf/project legends — which some plot helpers build by iterating the
        # colour dict — don't list each multi-word cell type twice.
        colors = colors_orig
        al.set_celltype_colors(colors)
        key = Path(a.output_dir) / "celltype_colors"
        if a.refresh_colors or not key.with_suffix(".pdf").exists():
            plot_color_key(colors_orig, categories, key)
            print(f"colour key -> {key}.{{png,pdf,svg}}", flush=True)

    for stage in ALL_STAGES:
        # patchify reuse is handled specially: it substitutes for the stage.
        if stage == "patchify" and a.patchify_dir and stage in selected:
            reuse_patchify(a)
            continue
        if not should_run(a, stage, selected):
            continue

        if stage == "patchify":
            stage_patchify(a)
        elif stage == "bptf":
            stage_bptf(a)
        elif stage in PLOTS_STAGE:
            print(f"\n---- {stage} ----", flush=True)
            if not a.dry_run:
                draw_figures(a, PLOTS_STAGE[stage], colors)
            else:
                print("     (dry-run: figures not drawn)", flush=True)
        elif stage == "project":
            stage_project(a)
        elif stage == "markers":
            stage_markers(a)
        elif stage == "glm":
            stage_glm(a)

    print(f"\n==== DONE  {utcnow()}  ->  {a.output_dir} ====", flush=True)
    print(
        "  compute : patchify/ bptf/ project/ markers/ glm/ motif_states.parquet",
        flush=True,
    )
    print("  figures : plots_original/ (original)   plots/ (workflow)", flush=True)
    print(
        "  colours : celltype_colors.json + celltype_colors.{png,pdf,svg}", flush=True
    )
    if logf:
        # Restore the real streams BEFORE closing the log file, so any writes
        # during interpreter shutdown don't hit a closed file through the tee.
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
