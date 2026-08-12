"""Shared figure conventions for the comparator benchmark.

CLAUDE.md requires every publication figure to use Arial, keep vector text
editable (`pdf.fonttype = 42`, `svg.fonttype = 'none'`) and be written as
png + pdf + svg **through a single saver**. Before this module seven scripts
under `scripts/comparators/` each re-declared that themselves, with three
different `save_all_formats` signatures and three different default dpi.

This is deliberately *not* in `src/alarmist/plotting/` — CLAUDE.md says to ask
before growing the package API, and nothing outside the comparator tree needs it.

Usage
-----
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/comparators
    from _common.plotting import apply_publication_style, save_all_formats

    apply_publication_style(**{"font.size": 6, "axes.spines.top": False})
    save_all_formats(fig, out_dir / "panel_a", dpi=450)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

# The CLAUDE.md baseline. Everything here is a house rule, not a taste call:
# Arial for the journal, fonttype 42 so pdf text stays selectable, svg.fonttype
# 'none' so Illustrator sees real text rather than outlined paths.
BASE_STYLE: dict[str, Any] = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}

FORMATS = ("png", "pdf", "svg")


def apply_publication_style(**overrides: Any) -> None:
    """Apply the CLAUDE.md baseline rcParams, then any per-script overrides.

    Overrides are ordinary rcParams keys, so a script that wants 6 pt text or
    no top/right spines passes them here instead of re-declaring the baseline.
    """
    import matplotlib

    matplotlib.rcParams.update(BASE_STYLE)
    if overrides:
        matplotlib.rcParams.update(overrides)


def save_all_formats(
    fig,
    stem,
    *,
    dpi: int | None = 300,
    bbox_inches: str | None = "tight",
    close: bool = False,
    verbose: bool = False,
) -> list[Path]:
    """Write `fig` as png + pdf + svg beside each other. Returns the paths written.

    `stem` is a path **without** an extension; parent directories are created.
    `dpi=None` leaves the figure's own dpi in force (matplotlib's default
    `savefig.dpi='figure'`), which is what a couple of callers relied on.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    kw: dict[str, Any] = {}
    if dpi is not None:
        kw["dpi"] = dpi
    if bbox_inches is not None:
        kw["bbox_inches"] = bbox_inches

    written = []
    for ext in FORMATS:
        path = stem.with_suffix(f".{ext}")
        fig.savefig(path, **kw)
        written.append(path)

    if close:
        import matplotlib.pyplot as plt

        plt.close(fig)
    if verbose:
        print(f"  {stem.name}.{{{','.join(FORMATS)}}}")
    return written
