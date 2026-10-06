#!/usr/bin/env python
"""The cell-type colours of this LUAD ALARMIST run, straight from the package.

`alarmist.plotting.utils.get_cell_type_colors` builds them as
`plt.get_cmap("tab20", len(cell_types))`, i.e. tab20 **resampled** to the number of cell
types (19 here), indexed over the h5ad's category order. That is not the same as plain tab20
by index: resampling to 19 shifts every assignment, so SMC is #7f7f7f grey and Tumor_epi
#bcbd22 olive, where plain tab20 gives SMC pink and Tumor_epi grey. Figures written before
2026-09-28 used the plain version and disagreed with the run's own figures.

Import needs the package on the path; it is not pip-installed in the comp-liana env, so this
module adds <repo>/src itself.
"""
from __future__ import annotations

import sys
from pathlib import Path

from matplotlib.colors import to_hex

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from alarmist.plotting.utils import get_cell_type_colors  # noqa: E402


# Presentation overrides used by the H&E overlays and the cell-type composition figures.
# Tumor_epi: the run's olive #bcbd22 is hard to read on H&E, so ColorBrewer Purples-8.
# A 3-cycle on request, starting from the run palette (SMC grey, Plasma pink, Langhans green):
# SMC takes Plasma's pink, Plasma takes Langhans_cell's green, Langhans_cell takes SMC's grey.
# All deliberate departures -- every other figure keeps get_cell_type_colors unchanged.
PRESENTATION = {"Tumor_epi": "#54278f", "SMC": "#f7b6d2", "Plasma": "#2ca02c",
                "Langhans_cell": "#7f7f7f"}


def celltype_colors(categories, override=None):
    """{cell type: hex} for the given category list, optionally with per-type overrides"""
    col = {c: to_hex(v) for c, v in get_cell_type_colors(list(categories)).items()}
    if override:
        col.update(override)
    return col
