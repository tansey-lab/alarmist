#!/usr/bin/env python
"""Write the gene panel of a prepped h5ad to a one-symbol-per-line file.

Feeds ``stlearn/build_cellchat_lrs.py``'s optional third positional argument, which filters
CellChatDB to the pairs whose genes are actually on the panel.

This matters: the existing GBM list is 1,371 of CellChatDB's 3,218 pairs precisely because it
was filtered to the GBM 5,119-gene panel. LUAD's panel is a DIFFERENT 5,101 genes, so reusing
the GBM list would test the wrong pair set and give a wrong ranking denominator.

Env: bptf.  Usage:
    python dump_panel_genes.py <prepped.h5ad> <out.tsv>
"""

import sys

import anndata as ad

if len(sys.argv) < 3:
    sys.exit(__doc__)

src, out = sys.argv[1], sys.argv[2]
a = ad.read_h5ad(src, backed="r")
genes = [str(g) for g in a.var_names]
with open(out, "w") as fh:
    fh.write("\n".join(genes) + "\n")
print(f"panel genes: {len(genes)} -> {out}")
