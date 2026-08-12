# Provenance of the skills in this directory

Only **`alarmist/`** and **`comparator-benchmark/`** are tracked in git — this repo owns
both. Everything else here is a vendored copy and is gitignored — update it at its source,
then re-copy, rather than editing in place.

| Skill | Source | Licence |
| --- | --- | --- |
| `alarmist/` | `~/tansey_lab/spatial_analysis_skills/skills/alarmist` (identical) | MIT, Tansey Lab — `LICENSE.spatial_analysis_skills` |
| `comparator-benchmark/` | **written in this repo**, not vendored. It is the contract `scripts/comparators/METHODS.md` cites throughout (`SKILL.md:45-46` on kernel scale, `:47-49` on multi-sample mode, `:51-54` on the two tiers, `:56-61` on requested-LR segregation) — those citations dangle if it is not tracked. Added to git 2026-08-12 | MIT, Tansey Lab |
| `spatial-workflow/` | `~/tansey_lab/spatial_analysis_skills/skills/` | same |
| `spatial-niche/` | same | same |
| `spatial-stats/` | same | same |
| `cohort-explore/` | same | same |
| `scientific-goals/` | same | same |
| `nature_publication_figures/` | `~/tansey_lab/es_xenium/.claude/skills/`, originally the `nature_plot_skills` project by Feiyang Huang | MIT — `nature_publication_figures/LICENSE`, **keep it intact** |
| `interactive-report/` | adapted from `es-interactive-report` in `~/tansey_lab/es_xenium` | — |

`interactive-report/` is the only one that was **modified** rather than copied verbatim
(output path, generalized tooltip axes, large-image `file://` fallback, English-only rule
kept). Its changes are recorded in a comment at the bottom of its `SKILL.md`. Copied
2026-07-30.

Not vendored, on purpose: `es-figures` and `es-de-units` from es_xenium — see the *Skills
in this repo* section of the repo `CLAUDE.md` for why.
