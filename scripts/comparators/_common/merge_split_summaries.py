#!/usr/bin/env python
"""Concatenate the per-section copies of SpatialDM's ``per_split_summary.csv``.

``run_spatialdm.py`` writes ``<out-dir>/per_split_summary.csv`` only after ALL of its splits
finish, and it OVERWRITES that file on every invocation. Because the four LUAD sections are
run as four separate invocations into one shared output directory (so a late failure does not
cost the earlier three), each run clobbers the previous summary. ``02_spatialdm.sh`` copies it
aside as ``per_split_summary_<section>.csv`` after each run; this stitches those back into the
single file that ``run_diff_spatialdm.py`` reads for its density sensitivity analysis.

Env: bptf.  Usage:
    python merge_split_summaries.py <out_dir>
"""

import glob
import os
import sys

import pandas as pd

if len(sys.argv) < 2:
    sys.exit(__doc__)

out_dir = sys.argv[1]
parts = sorted(glob.glob(os.path.join(out_dir, "per_split_summary_*.csv")))
if not parts:
    print(f"no per_split_summary_*.csv under {out_dir} -- nothing to merge")
    sys.exit(0)

df = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
dest = os.path.join(out_dir, "per_split_summary.csv")
df.to_csv(dest, index=False)
print(f"merged {len(parts)} per-section summaries ({len(df)} rows) -> {dest}")
for f in parts:
    print(f"    {os.path.basename(f)}")
