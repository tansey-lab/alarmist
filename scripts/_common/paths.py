"""Where things go — one resolver, so no script hardcodes a machine path again.

THE RULE: the repo holds **code**. Every byte a script *generates* goes to the data disk.
Three destinations, chosen by how big the output is and who reads it:

    results(proj)       big machine-readable output   — h5ad, .npy, per-cell tables, logs.
                                                        Never downloaded wholesale.
    deliverables(proj)  figures + small summary tables — the stuff you actually look at.
                                                        Kept small on purpose: one
                                                        `rsync` pulls the whole tree.
    scratch(proj)       node-local temp                — deleted whenever; never cite it.

On iris these live on /data1; the repo's `figures/` and `results/` are symlinks into the
first two, so a repo-relative path keeps working and `du` of $HOME stays small. On the Mac
laptop they default to repo-relative directories. Either way the *same* call works:

    from _common.paths import deliverables, results
    out = deliverables("dsrct_slide2")        # -> .../deliverables/dsrct_slide2/
    fig.savefig(out / "motif_06_lollipop.png")

Bootstrap from anywhere under scripts/:

    import sys; from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[N]))   # N up to scripts/
    from _common.paths import deliverables

Override any root with an env var (ALARMIST_RESULTS / ALARMIST_DELIVERABLES /
ALARMIST_SCRATCH) — useful for a one-off run you do not want mixed into the main tree.
"""

from __future__ import annotations

import os
import socket
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# --- per-host defaults -------------------------------------------------------------
# iris: everything generated lives on the group data disk. Note /data1/tanseyw was 91%
# full as of 2026-09-10 — check `df -h /data1/tanseyw` before starting anything large.
_IRIS_ROOT = Path("/data1/tanseyw/fanj2")


def _on_iris() -> bool:
    """True on the MSK iris cluster. Keyed on the data disk existing, not on hostname,
    so it also works from a compute node with a different name."""
    return _IRIS_ROOT.is_dir()


def _root(env_var: str, iris_sub: str, repo_sub: str) -> Path:
    override = os.environ.get(env_var)
    if override:
        return Path(override).expanduser()
    return (_IRIS_ROOT / iris_sub) if _on_iris() else (REPO / repo_sub)


def _make(root: Path, project: str | None, parts: tuple[str, ...]) -> Path:
    p = root if project is None else root / project
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def results(project: str | None = None, *parts: str) -> Path:
    """Big machine-readable output. Created if absent."""
    return _make(_root("ALARMIST_RESULTS", "results", "results"), project, parts)


def deliverables(project: str | None = None, *parts: str) -> Path:
    """Figures + small summary tables — the rsync-home set. Created if absent.

    Keep this tree small (target: the whole thing under a couple of GB). If an output is
    big or nobody reads it by eye, it belongs in results() instead.
    """
    return _make(
        _root("ALARMIST_DELIVERABLES", "deliverables", "figures"), project, parts
    )


def scratch(project: str | None = None, *parts: str) -> Path:
    """Node-local temp. Ephemeral — never cite a path under here in a report."""
    default = Path("/tmp") / f"alarmist-{os.environ.get('USER', 'user')}"
    override = os.environ.get("ALARMIST_SCRATCH")
    return _make(Path(override).expanduser() if override else default, project, parts)


def describe() -> str:
    """One-line summary of where this process will write. Print it at the top of a run
    so the log records the destinations, not just the inputs."""
    return (
        f"host={socket.gethostname()} iris={_on_iris()}\n"
        f"  results      -> {_root('ALARMIST_RESULTS', 'results', 'results')}\n"
        f"  deliverables -> {_root('ALARMIST_DELIVERABLES', 'deliverables', 'figures')}\n"
        f"  scratch      -> {scratch().parent}"
    )


if __name__ == "__main__":
    print(describe())
