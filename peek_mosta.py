#!/usr/bin/env python
"""
Remotely peek at the metadata of a large .h5ad over HTTP range requests.

Reads only the HDF5 metadata plus obs/var/obsm -- never touches the X data
body. Every byte actually pulled over the wire is counted against a hard
budget; the script aborts loudly rather than silently downloading the file.

Run with the bptf env:
    /Users/jiayifan/anaconda3/envs/bptf/bin/python peek_mosta.py
"""

from __future__ import annotations

import sys
import time
import traceback

import aiohttp
import fsspec
import h5py
import numpy as np
import pandas as pd

URL = (
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/"
    "Mouse_embryo_all_stage.h5ad"
)
REPORT_PATH = "peek_report.txt"

BUDGET_BYTES = 500 * 1024 * 1024  # hard stop: ~500 MB over the wire
BLOCK_SIZE = 512 * 1024  # fsspec read-block; smaller => less over-fetch
HTTP_TIMEOUT = 60  # seconds per request
MAX_RETRIES = 3


class BudgetExceeded(RuntimeError):
    """Raised when the download budget would be blown -- never caught silently."""


class Tee:
    """Write to stdout and the report file at once."""

    def __init__(self, path):
        self.fh = open(path, "w")

    def write(self, s):
        sys.__stdout__.write(s)
        sys.__stdout__.flush()
        self.fh.write(s)
        self.fh.flush()

    def flush(self):
        sys.__stdout__.flush()
        self.fh.flush()

    def close(self):
        self.fh.close()


class ByteMeter:
    """Counts bytes genuinely fetched from the network, and enforces the budget."""

    def __init__(self, budget):
        self.budget = budget
        self.total = 0
        self.n_requests = 0

    def charge(self, n):
        self.total += n
        self.n_requests += 1
        if self.total > self.budget:
            raise BudgetExceeded(
                f"Download budget exceeded: {self.total / 1e6:.1f} MB fetched "
                f"(limit {self.budget / 1e6:.0f} MB) over {self.n_requests} range "
                f"requests. Aborting rather than continuing to download."
            )

    def mb(self):
        return self.total / 1e6

    def would_exceed(self, n):
        return (self.total + n) > self.budget

    def remaining(self):
        return max(0, self.budget - self.total)


def instrument(f, meter):
    """Wrap HTTPFile._fetch_range so every real network fetch is metered + retried."""
    original = f._fetch_range

    def fetch(start, end):
        last = None
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                data = original(start, end)
                meter.charge(len(data))
                return data
            except BudgetExceeded:
                raise  # never retry past the budget
            except Exception as exc:  # noqa: BLE001 - report, don't hang
                last = exc
                if attempt < MAX_RETRIES:
                    wait = 2**attempt
                    print(
                        f"    [retry {attempt}/{MAX_RETRIES}] range {start}-{end} "
                        f"failed ({type(exc).__name__}: {exc}); waiting {wait}s"
                    )
                    time.sleep(wait)
        raise RuntimeError(
            f"Range request {start}-{end} failed after {MAX_RETRIES} attempts: "
            f"{type(last).__name__}: {last}"
        ) from last

    f._fetch_range = fetch
    # fsspec's read-cache captured the *original* bound method in
    # AbstractBufferedFile.__init__, so patching the attribute alone is not
    # enough -- the cache must be repointed or the meter silently reads zero.
    if getattr(f, "cache", None) is not None:
        f.cache.fetcher = fetch
    return f


# get_storage_size() walks the whole chunk B-tree. On a multi-GB remote dataset
# that is hundreds of MB of scattered metadata reads and effectively hangs, so
# only ask for it on datasets small enough for the walk to be cheap.
STORAGE_SIZE_LIMIT = 64 * 1024 * 1024  # logical bytes


def stored_bytes(obj, default=0):
    """
    Compressed on-disk size of a dataset (what a full read would actually cost).

    Returns None for datasets too large to interrogate remotely.
    """
    if not isinstance(obj, h5py.Dataset):
        return default
    if obj.nbytes > STORAGE_SIZE_LIMIT:
        return None
    try:
        return obj.id.get_storage_size()
    except Exception:  # noqa: BLE001
        return obj.nbytes


def fmt_stored(d):
    sb = stored_bytes(d)
    if sb is None:
        return "not queried (chunk-index walk too costly remotely)"
    return f"{sb:,}"


def describe_dataset(d, indent="    "):
    lines = [
        f"{indent}shape={d.shape} dtype={d.dtype}",
        f"{indent}chunks={d.chunks} compression={d.compression}"
        f"{'' if d.compression_opts is None else f' opts={d.compression_opts}'}"
        f" shuffle={d.shuffle}",
        f"{indent}stored (compressed) bytes = {fmt_stored(d)}"
        f"  |  logical bytes = {d.nbytes:,}",
    ]
    return "\n".join(lines)


def cost_of(d):
    """Bytes a full read of `d` would pull; falls back to logical size when unknown."""
    sb = stored_bytes(d)
    return d.nbytes if sb is None else sb


def attrs_of(obj):
    out = {}
    for k, v in obj.attrs.items():
        if isinstance(v, bytes):
            v = v.decode()
        elif isinstance(v, np.ndarray):
            v = v.tolist()
        out[k] = v
    return out


def decode(arr):
    """bytes -> str for h5py string arrays."""
    arr = np.asarray(arr)
    if arr.dtype.kind in ("S", "O"):
        return np.array(
            [x.decode() if isinstance(x, bytes) else str(x) for x in arr.ravel()]
        ).reshape(arr.shape)
    return arr


LEGACY_CATS = "__categories"


def categorical_parts(grp, name):
    """
    Return (categories_dataset, codes_dataset) if `name` is categorical, else None.

    Handles both encodings:
      * modern h5ad (>=0.8): a Group holding {categories, codes}
      * legacy h5ad (<0.8):  an integer Dataset whose 'categories' attr is an
        HDF5 object reference into obs/__categories/<name>
    """
    obj = grp[name]

    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        return obj["categories"], obj["codes"]

    if isinstance(obj, h5py.Dataset) and "categories" in obj.attrs:
        ref = obj.attrs["categories"]
        if isinstance(ref, h5py.Reference):
            return grp.file[ref], obj
        if LEGACY_CATS in grp and name in grp[LEGACY_CATS]:
            return grp[LEGACY_CATS][name], obj

    if LEGACY_CATS in grp and name in grp[LEGACY_CATS] and isinstance(obj, h5py.Dataset):
        return grp[LEGACY_CATS][name], obj

    return None


def read_column(grp, name, meter, log):
    """
    Return (pandas Series or None, kind, cost_bytes).

    Reads a categorical column as codes+categories, a plain dataset directly.
    Refuses to read anything that would blow the remaining budget.
    """
    obj = grp[name]

    parts = categorical_parts(grp, name)
    if parts is not None:
        cats_d, codes_d = parts
        cost = cost_of(cats_d) + cost_of(codes_d)
        if meter.would_exceed(cost):
            log(
                f"    !! SKIPPED: would cost ~{cost / 1e6:.1f} MB, only "
                f"{meter.remaining() / 1e6:.1f} MB of budget left"
            )
            return None, "categorical", cost
        cats = decode(cats_d[:])
        codes = codes_d[:]
        ser = pd.Categorical.from_codes(codes, categories=list(cats))
        return pd.Series(ser), "categorical", cost

    dt = obj.dtype
    kind = "string" if dt.kind in ("S", "O", "U") else f"numeric[{dt}]"
    cost = cost_of(obj)
    if meter.would_exceed(cost):
        log(
            f"    !! SKIPPED: would cost ~{cost / 1e6:.1f} MB, only "
            f"{meter.remaining() / 1e6:.1f} MB of budget left"
        )
        return None, kind, cost
    vals = obj[:]
    if dt.kind in ("S", "O"):
        vals = decode(vals)
    return pd.Series(vals), kind, cost


def contiguous_runs(mask):
    """Compact a boolean mask into [(start, end_inclusive), ...] row ranges."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    brk = np.flatnonzero(np.diff(idx) != 1)
    starts = np.concatenate(([idx[0]], idx[brk + 1]))
    ends = np.concatenate((idx[brk], [idx[-1]]))
    return list(zip(starts.tolist(), ends.tolist()))


def main():
    tee = Tee(REPORT_PATH)
    log = lambda s="": tee.write(s + "\n")  # noqa: E731

    meter = ByteMeter(BUDGET_BYTES)
    log("=" * 78)
    log("REMOTE h5ad METADATA PEEK  (range requests only, X body never read)")
    log("=" * 78)
    log(f"URL    : {URL}")
    log(f"budget : {BUDGET_BYTES / 1e6:.0f} MB   block_size: {BLOCK_SIZE // 1024} KiB")
    log("")

    fs = fsspec.filesystem(
        "http",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=HTTP_TIMEOUT)},
        block_size=BLOCK_SIZE,
    )
    size = fs.size(URL)
    log(f"remote size: {size:,} bytes ({size / 1024**3:.2f} GiB)")

    f = instrument(fs.open(URL, "rb", block_size=BLOCK_SIZE), meter)

    with h5py.File(f, "r") as h5:
        log(f"root keys  : {list(h5.keys())}")
        log(f"root attrs : {attrs_of(h5)}")
        log(f"[wire so far: {meter.mb():.1f} MB in {meter.n_requests} requests]")

        # ---------------- shape ----------------
        obs, var = h5["obs"], h5["var"]
        obs_idx_key = obs.attrs.get("_index", b"_index")
        var_idx_key = var.attrs.get("_index", b"_index")
        obs_idx_key = (
            obs_idx_key.decode() if isinstance(obs_idx_key, bytes) else obs_idx_key
        )
        var_idx_key = (
            var_idx_key.decode() if isinstance(var_idx_key, bytes) else var_idx_key
        )
        n_obs = obs[obs_idx_key].shape[0]
        n_var = var[var_idx_key].shape[0]

        log("")
        log("#" * 78)
        log("# 1. SHAPE")
        log("#" * 78)
        log(f"n_obs = {n_obs:,}")
        log(f"n_var = {n_var:,}")

        # ---------------- 5. X storage (before any bulk reads) ----------------
        log("")
        log("#" * 78)
        log("# 5. X STORAGE FORMAT")
        log("#" * 78)
        X = h5["X"]
        log(f"h5['X'] is a {type(X).__name__}")
        log(f"attrs: {attrs_of(X)}")
        if isinstance(X, h5py.Group):
            log(f"members: {list(X.keys())}")
            for k in X.keys():
                log(f"  X/{k}:")
                log(describe_dataset(X[k], indent="      "))
        else:
            log(describe_dataset(X, indent="  "))

        # ---------------- 6. layers / raw ----------------
        log("")
        log("#" * 78)
        log("# 6. layers / raw")
        log("#" * 78)
        for key in ("layers", "raw"):
            if key in h5:
                node = h5[key]
                log(f"{key}: PRESENT ({type(node).__name__})")
                log(f"  attrs: {attrs_of(node)}")
                if isinstance(node, h5py.Group):
                    log(f"  keys: {list(node.keys())}")
                    if key == "layers":
                        for lk in node.keys():
                            lo = node[lk]
                            log(f"  layers['{lk}'] is a {type(lo).__name__}")
                            log(f"    attrs: {attrs_of(lo)}")
                            if isinstance(lo, h5py.Group):
                                for sk in lo.keys():
                                    log(f"    layers/{lk}/{sk}:")
                                    log(describe_dataset(lo[sk], indent="        "))
                            else:
                                log(describe_dataset(lo, indent="    "))
                    if key == "raw" and "X" in node:
                        rx = node["X"]
                        log(f"  raw/X is a {type(rx).__name__}, attrs={attrs_of(rx)}")
                        if isinstance(rx, h5py.Group):
                            for k in rx.keys():
                                log(f"    raw/X/{k}: shape={rx[k].shape} dtype={rx[k].dtype}")
                        else:
                            log(describe_dataset(rx, indent="    "))
                    if key == "raw" and "var" in node:
                        log(f"  raw/var keys: {list(node['var'].keys())}")
            else:
                log(f"{key}: ABSENT")

        # ---------------- 4. obsm ----------------
        log("")
        log("#" * 78)
        log("# 4. obsm")
        log("#" * 78)
        for grp_name in ("obsm", "varm", "obsp", "uns"):
            if grp_name not in h5:
                log(f"{grp_name}: ABSENT")
                continue
            g = h5[grp_name]
            keys = list(g.keys())
            log(f"{grp_name}: {keys}")
            if grp_name in ("obsm", "varm"):
                for k in keys:
                    o = g[k]
                    if isinstance(o, h5py.Dataset):
                        log(
                            f"  {grp_name}['{k}']: shape={o.shape} dtype={o.dtype} "
                            f"chunks={o.chunks} compression={o.compression} "
                            f"stored={fmt_stored(o)} B"
                        )
                    else:
                        log(f"  {grp_name}['{k}']: Group, keys={list(o.keys())}")

        log("")
        log(f"[wire so far: {meter.mb():.1f} MB in {meter.n_requests} requests]")

        # ---------------- 3. var ----------------
        log("")
        log("#" * 78)
        log("# 3. var")
        log("#" * 78)
        log(f"var attrs   : {attrs_of(var)}")
        log(f"var index   : '{var_idx_key}'")
        log(f"var columns : {list(var.keys())}")
        for k in var.keys():
            o = var[k]
            if isinstance(o, h5py.Dataset):
                log(f"  var['{k}']: shape={o.shape} dtype={o.dtype} stored={fmt_stored(o)} B")
            else:
                log(f"  var['{k}']: categorical Group {list(o.keys())}")
        first20 = decode(var[var_idx_key][:20])
        log("")
        log("first 20 var index entries:")
        for i, g in enumerate(first20):
            log(f"  [{i:2d}] {g}")

        # ---------------- 2. obs ----------------
        log("")
        log("#" * 78)
        log("# 2. obs")
        log("#" * 78)
        log(f"obs attrs : {attrs_of(obs)}")
        log(f"obs index : '{obs_idx_key}'")
        col_order = obs.attrs.get("column-order")
        if col_order is not None:
            col_order = [c.decode() if isinstance(c, bytes) else c for c in col_order]
            log(f"column-order ({len(col_order)}): {col_order}")
        all_keys = list(obs.keys())
        log(f"all obs keys ({len(all_keys)}): {all_keys}")
        if LEGACY_CATS in obs:
            log(
                f"NOTE: legacy h5ad encoding -- categories live in "
                f"obs/{LEGACY_CATS}/: {list(obs[LEGACY_CATS].keys())}"
            )
        cols = [c for c in all_keys if c != LEGACY_CATS]

        log("")
        log("--- per-column inventory (cost = compressed bytes a full read pulls) ---")
        plan = []
        for c in cols:
            o = obs[c]
            parts = categorical_parts(obs, c)
            if parts is not None:
                cats_d, codes_d = parts
                cost = cost_of(cats_d) + cost_of(codes_d)
                log(
                    f"  {c:28s} categorical  n_categories={cats_d.shape[0]:<8,} "
                    f"codes_dtype={codes_d.dtype}  cost~{cost / 1e6:8.2f} MB"
                )
                plan.append((c, True, cost))
            else:
                cost = cost_of(o)
                is_str = o.dtype.kind in ("S", "O", "U")
                log(
                    f"  {c:28s} {'string' if is_str else str(o.dtype):12s} "
                    f"shape={str(o.shape):16s} cost~{cost / 1e6:8.2f} MB"
                )
                plan.append((c, is_str, cost))

        total_cost = sum(c for _, want, c in plan if want)
        log("")
        log(
            f"object/category columns to read: "
            f"{sum(1 for _, w, _ in plan if w)}  |  total ~{total_cost / 1e6:.1f} MB"
        )
        log(f"[wire so far: {meter.mb():.1f} MB]")

        log("")
        log("--- value_counts().head(20) ---")
        loaded = {}
        for c, want, cost in plan:
            if not want:
                continue
            log("")
            log(f"### obs['{c}']   (~{cost / 1e6:.2f} MB)")
            try:
                ser, kind, _ = read_column(obs, c, meter, log)
            except BudgetExceeded:
                raise
            except Exception as exc:  # noqa: BLE001
                log(f"    !! FAILED: {type(exc).__name__}: {exc}")
                continue
            if ser is None:
                continue
            loaded[c] = ser
            vc = ser.value_counts()
            log(f"    kind={kind}  n_unique={ser.nunique():,}  n_rows={len(ser):,}")
            for lbl, n in vc.head(20).items():
                log(f"      {str(lbl):40s} {n:>12,}")
            if len(vc) > 20:
                log(f"      ... {len(vc) - 20:,} more distinct values")
            log(f"    [wire so far: {meter.mb():.1f} MB]")

        # ---------------- judgement calls ----------------
        log("")
        log("#" * 78)
        log("# JUDGEMENT CALLS")
        log("#" * 78)
        log(f"n_obs = {n_obs:,}")
        log(
            "MOSTA cell_bin across all 53 sections is ~10^7 cells; bin50 is 1-2 "
            "orders of magnitude smaller (~10^5-10^6)."
        )

        # Hunt for the section / stage columns among what we loaded.
        log("")
        log("--- candidate columns for section id and developmental stage ---")
        for c, ser in loaded.items():
            u = ser.astype(str).unique()
            sample = sorted(u)[:25]
            looks_section = any(
                ("E" in str(x) and "S" in str(x)) for x in u[: min(len(u), 500)]
            )
            looks_stage = any(str(x).startswith("E") and "." in str(x) for x in u[:500])
            tag = []
            if looks_section:
                tag.append("SECTION?")
            if looks_stage:
                tag.append("STAGE?")
            log(f"  {c:28s} n_unique={ser.nunique():<8,} {' '.join(tag)}")
            log(f"      sample: {sample}")

        # E12.5 / E1S1 / E1S2 localisation
        log("")
        log("--- locating E12.5 E1S1 / E1S2 ---")
        found_any = False
        for c, ser in loaded.items():
            s = ser.astype(str)
            for token in ("E1S1", "E1S2"):
                hit = s.str.contains(token, regex=False, na=False)
                n_hit = int(hit.sum())
                if n_hit == 0:
                    continue
                found_any = True
                matched_vals = sorted(s[hit].unique())[:10]
                log("")
                log(f"  obs['{c}'] contains '{token}': {n_hit:,} rows True")
                log(f"    matching values: {matched_vals}")
                runs = contiguous_runs(hit.to_numpy())
                log(f"    {len(runs)} contiguous run(s) of row indices")
                for a, b in runs[:10]:
                    log(f"      rows {a:,} .. {b:,}  ({b - a + 1:,} rows)")
                if len(runs) > 10:
                    log(f"      ... {len(runs) - 10} more runs")
        if not found_any:
            log("  No obs column contained the literal 'E1S1' / 'E1S2'.")

        # The obs index looks like "<x>_<y>-<n>". The "-<n>" suffix is the
        # AnnData.concatenate batch tag, i.e. the only surviving trace of which
        # section each row came from. Cross-tab it against timepoint to see how
        # many sections are actually in this file.
        log("")
        log("--- obs index suffix ('-N' concat batch tag) x timepoint ---")
        if obs_idx_key in loaded and "timepoint" in loaded:
            names = loaded[obs_idx_key].astype(str)
            suffix = names.str.rsplit("-", n=1).str[-1]
            tp = loaded["timepoint"].astype(str)
            log(f"distinct suffixes: {sorted(suffix.unique())}")
            ct = pd.crosstab(suffix, tp)
            log("")
            log(ct.to_string())
            log("")
            log(
                f"=> {suffix.nunique()} distinct batch tags vs "
                f"{tp.nunique()} timepoints"
            )
            log("")
            log("--- contiguous row ranges per timepoint ---")
            for t in sorted(tp.unique()):
                mask = (tp == t).to_numpy()
                runs = contiguous_runs(mask)
                log(f"  {t}: {int(mask.sum()):>8,} rows, {len(runs)} run(s)")
                for a, b in runs[:5]:
                    log(f"      rows {a:,} .. {b:,}  ({b - a + 1:,})")
                if len(runs) > 5:
                    log(f"      ... {len(runs) - 5} more runs")

            log("")
            log("--- E12.5 subset (the stage E1S1/E1S2 belong to) ---")
            m125 = (tp == "E12.5").to_numpy()
            log(f"  boolean mask True count: {int(m125.sum()):,}")
            for a, b in contiguous_runs(m125):
                log(f"  row range: {a:,} .. {b:,}  ({b - a + 1:,} rows)")
            sub = suffix[m125]
            log(f"  batch tags present within E12.5: {sorted(sub.unique())}")
            log(f"  cells per tag within E12.5: {sub.value_counts().to_dict()}")
        else:
            log("  (index or timepoint column not loaded -- cannot cross-tab)")

        log("")
        log("=" * 78)
        log(
            f"DONE. Total fetched over the wire: {meter.mb():.1f} MB "
            f"in {meter.n_requests} range requests "
            f"({100 * meter.total / size:.4f}% of the 21.4 GiB file)."
        )
        log("X data body was never read.")
        log("=" * 78)

    tee.close()


if __name__ == "__main__":
    try:
        main()
    except BudgetExceeded as exc:
        print(f"\n*** ABORTED ON BUDGET ***\n{exc}", file=sys.stderr)
        sys.exit(2)
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        sys.exit(1)
