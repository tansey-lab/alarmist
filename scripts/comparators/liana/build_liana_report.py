#!/usr/bin/env python
"""Build the self-contained LIANA+ / GBM methods report.

Writes reports/liana_plus_GBM_cellchatdb2/liana_plus_GBM_methods.html --- a tutorial-detail
English walkthrough of every LIANA+ call we made on the GBM Xenium TMA, with the reason for
every departure from the authors' tutorials and from this project's own benchmark contract.

Per `.claude/skills/interactive-report`:
  * single self-contained .html, vanilla JS + inline CSS, no CDN, works offline
  * English only
  * figures embedded as base64 below INLINE_MAX_BYTES, linked file:// above
  * THIS SCRIPT is the tracked artifact; reports/ is gitignored

Figure policy (user instruction, 2026-08-04): only figures that carry the BIOLOGICAL story,
or that stop a biological claim being over-read. QC, calibration, connectivity, elbow,
factor-correlation and plain cell-type maps are excluded.

Run:
    /Users/jiayifan/anaconda3/envs/bptf/bin/python scripts/comparators/liana/build_liana_report.py
"""
from __future__ import annotations

import base64
import html as _html
import json
import subprocess
from pathlib import Path

ROOT = Path("/Users/jiayifan/tansey_lab/alarmist")
REPORT_DIR = ROOT / "reports" / "liana_plus_GBM_cellchatdb2"
OUT = REPORT_DIR / "liana_plus_GBM_methods.html"
SECTIONS_JSON = Path(__file__).with_name("_liana_report_sections.json")

# Inline images up to this size; larger ones link out and lazy-load.
INLINE_MAX_BYTES = 400_000


def img_src(rel: str) -> str | None:
    """data: URI for small images, file:// for large ones, None if absent."""
    p = REPORT_DIR / rel
    if not p.exists():
        return None
    if p.stat().st_size <= INLINE_MAX_BYTES:
        return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()
    return f"file://{p}"


def git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


CSS = """
:root{
  --bg:#ffffff; --fg:#1a1a1a; --muted:#5c6470; --rule:#e3e6ea; --accent:#1f4e79;
  --code-bg:#f6f8fa; --warn-bg:#fff8e6; --warn-bd:#e0b000;
  --dev-bg:#f2f7fb; --dev-bd:#8fb8d8; --bad-bg:#fdf1f1; --bad-bd:#d98a8a;
}
@media (prefers-color-scheme: dark){
  :root{
    --bg:#15181c; --fg:#e6e8ea; --muted:#9aa4b0; --rule:#2b3038; --accent:#7fb3e0;
    --code-bg:#1d2127; --warn-bg:#2b2410; --warn-bd:#8a6d00;
    --dev-bg:#16222c; --dev-bd:#3a6d96; --bad-bg:#2a1a1a; --bad-bd:#8a4a4a;
  }
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:15px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;}
.wrap{display:grid;grid-template-columns:290px minmax(0,1fr);gap:0;max-width:1500px;margin:0 auto}
nav{position:sticky;top:0;align-self:start;max-height:100vh;overflow-y:auto;
  padding:22px 16px 40px;border-right:1px solid var(--rule);font-size:13px}
nav h2{font-size:12px;text-transform:uppercase;letter-spacing:.08em;color:var(--muted);margin:18px 0 8px}
nav a{display:block;padding:3px 8px;color:var(--fg);text-decoration:none;border-radius:4px;
  border-left:2px solid transparent}
nav a:hover{background:var(--code-bg)}
nav a.sub{padding-left:20px;color:var(--muted);font-size:12.5px}
nav a.active{border-left-color:var(--accent);color:var(--accent);font-weight:600}
main{padding:28px 40px 120px;min-width:0}
h1{font-size:27px;line-height:1.25;margin:0 0 6px}
h2{font-size:21px;margin:44px 0 10px;padding-top:10px;border-top:1px solid var(--rule)}
h3{font-size:17px;margin:26px 0 8px;color:var(--accent)}
h4{font-size:15px;margin:18px 0 6px}
p,li{max-width:78ch}
code{background:var(--code-bg);padding:1px 5px;border-radius:4px;font-size:12.8px;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
pre{background:var(--code-bg);padding:13px 15px;border-radius:7px;overflow-x:auto;
  border:1px solid var(--rule);font-size:12.6px;line-height:1.5}
pre code{background:none;padding:0;font-size:12.6px}
table{border-collapse:collapse;margin:14px 0;font-size:13.4px;width:100%}
th,td{border:1px solid var(--rule);padding:7px 10px;text-align:left;vertical-align:top}
th{background:var(--code-bg);font-weight:600}
blockquote{margin:14px 0;padding:10px 16px;background:var(--warn-bg);
  border-left:4px solid var(--warn-bd);border-radius:0 6px 6px 0}
blockquote p{margin:6px 0}
figure{margin:20px 0;padding:0}
figure img{display:block;max-width:100%;height:auto;border:1px solid var(--rule);border-radius:6px}
figcaption{font-size:13px;color:var(--muted);margin-top:7px;max-width:78ch}
.lead{font-size:16px;color:var(--muted);max-width:78ch}
.meta{font-size:12.5px;color:var(--muted);margin:16px 0 0}
.tag{display:inline-block;font-size:11px;letter-spacing:.05em;text-transform:uppercase;
  padding:2px 7px;border-radius:20px;border:1px solid var(--rule);color:var(--muted);margin-right:6px}
.dev{background:var(--dev-bg);border-left:4px solid var(--dev-bd);padding:10px 16px;
  border-radius:0 6px 6px 0;margin:14px 0}
.bad{background:var(--bad-bg);border-left:4px solid var(--bad-bd);padding:10px 16px;
  border-radius:0 6px 6px 0;margin:14px 0}
.missing{font-size:12.5px;color:var(--muted);font-style:italic}
@media (max-width:1000px){
  .wrap{grid-template-columns:1fr}
  nav{position:static;max-height:none;border-right:none;border-bottom:1px solid var(--rule)}
  main{padding:20px}
}
"""

JS = """
// active-section highlighting in the sidebar
const links = [...document.querySelectorAll('nav a')];
const targets = links.map(a => document.getElementById(a.getAttribute('href').slice(1)))
                     .filter(Boolean);
function sync(){
  let best = null, bestTop = -Infinity;
  for (const t of targets){
    const top = t.getBoundingClientRect().top;
    if (top < 120 && top > bestTop){ bestTop = top; best = t; }
  }
  links.forEach(a => a.classList.toggle('active',
    best && a.getAttribute('href') === '#' + best.id));
}
document.addEventListener('scroll', sync, {passive:true});
sync();
"""


import re as _re


def anchor_of(s: dict) -> str:
    """The id the nav should point at.

    Section fragments are authored independently and give their own <h2> an id that need not
    equal section_id. Point at the fragment's first id when section_id is not actually present,
    otherwise the nav link silently goes nowhere.
    """
    sid = s["section_id"]
    if f'id="{sid}"' in s["html"]:
        return sid
    m = _re.search(r'<h2 id="([^"]+)"', s["html"]) or _re.search(r'id="([^"]+)"', s["html"])
    return m.group(1) if m else sid


def build_nav(sections: list[dict]) -> str:
    out = ["<nav><h2>Contents</h2>"]
    for s in sections:
        out.append(f'<a href="#{anchor_of(s)}">{_html.escape(s["title"])}</a>')
        for sub_id, sub_title in s.get("subheads", []):
            # subheads are scraped from the fragments' <h3> text, so they still carry HTML
            # entities (&mdash;, &amp;). Unescape before re-escaping, or they render literally.
            clean = _html.escape(_html.unescape(sub_title))
            out.append(f'<a class="sub" href="#{sub_id}">{clean}</a>')
    out.append("</nav>")
    return "\n".join(out)


def figure_html(rel: str, caption: str) -> str:
    src = img_src(rel)
    cap = _html.escape(caption)
    if src is None:
        return (f'<p class="missing">[figure not found on disk: '
                f'<code>{_html.escape(rel)}</code>]</p>')
    lazy = ' loading="lazy"' if src.startswith("file://") else ""
    return (f'<figure><img src="{src}" alt="{cap}"{lazy}>'
            f'<figcaption><strong>{_html.escape(Path(rel).name)}</strong> — {cap}'
            f'</figcaption></figure>')


def main() -> None:
    if not SECTIONS_JSON.exists():
        raise SystemExit(
            f"missing {SECTIONS_JSON}\n"
            "It holds the report prose (list of {section_id,title,html,figures,subheads}).\n"
            "It is written alongside this builder and is part of the tracked artifact.")
    sections = json.loads(SECTIONS_JSON.read_text())

    body, n_inline, n_linked, n_missing = [], 0, 0, 0
    for s in sections:
        body.append(s["html"])
        for f in s.get("figures", []):
            src = img_src(f["path"])
            if src is None:
                n_missing += 1
            elif src.startswith("data:"):
                n_inline += 1
            else:
                n_linked += 1
            body.append(figure_html(f["path"], f["caption"]))

    n_png = sum(1 for _ in REPORT_DIR.rglob("*.png"))
    head = f"""<h1>LIANA+ on the GBM Xenium TMA — what was run, and why</h1>
<p class="lead">A step-by-step account of every LIANA+ call made on this dataset, at the level of
detail of the authors' own notebooks, together with the reason for every departure from those
notebooks and from this project's benchmark contract.</p>
<p class="meta">
<span class="tag">liana 1.8.1</span><span class="tag">env comp-liana</span>
<span class="tag">CellChatDB v2</span><span class="tag">100,197 cells</span>
<span class="tag">13 TMA punches</span><span class="tag">git {git_sha()}</span></p>
<p class="meta">Companion figure tree: <code>{n_png}</code> PNG in this directory, indexed by
<code>README.md</code> and <code>figure_manifest.json</code>. Rebuild this page with
<code>scripts/comparators/liana/build_liana_report.py</code>.</p>"""

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LIANA+ on the GBM Xenium TMA — methods walkthrough</title>
<style>{CSS}</style></head><body>
<div class="wrap">
{build_nav(sections)}
<main>
{head}
{"".join(body)}
</main></div>
<script>{JS}</script>
</body></html>"""

    OUT.write_text(doc)
    mb = OUT.stat().st_size / 1048576
    print(f"wrote {OUT}  ({mb:.1f} MB)")
    print(f"  sections {len(sections)} | figures inline {n_inline}, linked {n_linked}, "
          f"missing {n_missing}")
    if n_missing:
        print("  WARNING: some figure paths did not resolve — check them, a wrong path is silent")


if __name__ == "__main__":
    main()
