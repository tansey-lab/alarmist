#!/usr/bin/env python
"""Stamp the correct dataset label onto a run_manifest.json after the run.

Two of the runners we use hardcode ``"dataset": "GBM"`` in the manifest they write:

    scripts/comparators/spatialdm/run_diff_spatialdm.py:357
    scripts/comparators/liana/run_inflow_downstream.py:435   (also ``"tier": "cellchatdb2"``)

Both are tracked scripts that produced signed-off GBM output, and the instruction for this
work is not to edit them. The label is metadata written AFTER the computation, so correcting
it afterwards is exactly equivalent to having passed a flag -- and it leaves the tracked
scripts byte-identical.

SKILL.md requires the manifest to name the dataset, so shipping LUAD results labelled "GBM"
is not an option; this closes that gap without touching the runner.

The original value is preserved under ``dataset_as_written_by_runner`` so the correction is
never invisible.

Env: any python 3.  Usage:
    python fix_manifest.py --dataset LUAD <run_manifest.json> [...]
    python fix_manifest.py --dataset LUAD --tier cellchatdb2 --find <dir>
"""

from __future__ import annotations

import argparse
import json
import os
import sys


def patch(path, dataset, tier):
    try:
        with open(path) as fh:
            man = json.load(fh)
    except Exception as e:  # noqa: BLE001
        print(f"  SKIP {path}: unreadable ({e})")
        return False
    if not isinstance(man, dict):
        print(f"  SKIP {path}: top level is {type(man).__name__}, not an object")
        return False

    changed = []
    if man.get("dataset") != dataset:
        if "dataset" in man and "dataset_as_written_by_runner" not in man:
            man["dataset_as_written_by_runner"] = man["dataset"]
        changed.append(f"dataset {man.get('dataset')!r} -> {dataset!r}")
        man["dataset"] = dataset
    if tier is not None and man.get("tier") != tier:
        if "tier" in man and "tier_as_written_by_runner" not in man:
            man["tier_as_written_by_runner"] = man["tier"]
        changed.append(f"tier {man.get('tier')!r} -> {tier!r}")
        man["tier"] = tier

    if not changed:
        print(f"  ok   {path} (already correct)")
        return False
    man["dataset_label_corrected_by"] = os.path.abspath(__file__)
    with open(path, "w") as fh:
        json.dump(man, fh, indent=2)
    print(f"  FIX  {path}: " + "; ".join(changed))
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="LUAD")
    p.add_argument("--tier", default=None,
                   help="also stamp the tier (run_inflow_downstream.py hardcodes it too)")
    p.add_argument("--find", default=None,
                   help="walk this directory and patch every run_manifest.json under it")
    p.add_argument("paths", nargs="*", help="explicit run_manifest.json paths")
    a = p.parse_args()

    targets = list(a.paths)
    if a.find:
        for root, _dirs, files in os.walk(a.find):
            for f in files:
                if f == "run_manifest.json":
                    targets.append(os.path.join(root, f))
    if not targets:
        sys.exit("nothing to do: give manifest paths or --find <dir>")

    n = sum(patch(t, a.dataset, a.tier) for t in sorted(set(targets)))
    print(f"{n} manifest(s) corrected out of {len(set(targets))} inspected")


if __name__ == "__main__":
    main()
