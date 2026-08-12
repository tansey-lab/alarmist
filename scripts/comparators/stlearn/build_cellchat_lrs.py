#!/usr/bin/env python
"""Convert CellChatDB v2 into the flat `L_R` gene-pair list stLearn expects.

stLearn's LR format (see `st.tl.cci.load_lrs`) is a flat array of "LIGAND_RECEPTOR" strings with
**no support for multi-subunit complexes** -- the underscore is the L/R separator, so it cannot
also mean "subunit of". CellChatDB v2 encodes complexes with underscores (TGFBR1_TGFBR2), which
makes 57.5% of its rows unrepresentable.

We therefore DROP every row whose ligand or receptor is a complex. We do not expand complexes
combinatorially: TGFB1 + TGFBR1_TGFBR2 -> {TGFB1_TGFBR1, TGFB1_TGFBR2} would assert two LR pairs
that CellChat does not claim exist independently, and would double-count one interaction. This is
a limitation of stLearn, recorded in DEVIATIONS.md, not a modelling preference.

Usage: python build_cellchat_lrs.py <cellchat_csv> <out_txt> [panel_genes.tsv]
"""
import sys
import pandas as pd

csv = sys.argv[1] if len(sys.argv) > 1 else "data/LRdatabase/CellChatDBv2.0.human.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "results/comparators/stlearn/cellchatdb2_lrs.txt"
panel_f = sys.argv[3] if len(sys.argv) > 3 else None

db = pd.read_csv(csv)
lig, rec = db["ligand"].astype(str), db["receptor"].astype(str)
is_complex = lig.str.contains("_") | rec.str.contains("_")

simple = db.loc[~is_complex].copy()
pairs = (simple["ligand"].astype(str) + "_" + simple["receptor"].astype(str))
pairs = pairs.drop_duplicates().sort_values().tolist()

print(f"CellChatDB v2 rows            : {len(db)}")
print(f"  dropped (complex L or R)    : {int(is_complex.sum())}  ({100 * is_complex.mean():.1f}%)")
print(f"  kept (simple 1:1)           : {len(simple)}")
print(f"  unique L_R pairs            : {len(pairs)}")

if panel_f:
    panel = set(open(panel_f).read().split())
    on = [p for p in pairs if all(g in panel for g in p.split("_", 1))]
    print(f"  on panel ({len(panel)} genes)     : {len(on)}")
    print("  NOTE: stLearn's run() does its own filtering to genes present in the data;")
    print("        this count is reported for the record, the full list is written out.")

with open(out, "w") as fh:
    fh.write("\n".join(pairs) + "\n")
print(f"wrote -> {out}")
