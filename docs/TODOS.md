# ALARMIST: Reproducibility & Portability TODOs

Upgrade plan to make this project reproducible and portable, following the `ebbf` project framework.

**Reference**: `~/Code/ebbf` for patterns and templates.

---

## Phase 1: Project Infrastructure (High Priority)

### 1. Migrate to uv for Dependency Management
Replace conda/pip with uv for fast, reproducible dependency resolution.

- [x] Update `pyproject.toml` with proper dependency groups:
  ```toml
  [project]
  requires-python = ">=3.10"
  dependencies = [...]

  [project.optional-dependencies]
  dev = ["pytest>=7.0", "pytest-cov>=4.0", "ruff>=0.1.0", "pre-commit>=3.0", "nf-core", "multiqc"]
  ```
- [x] Generate lock file: `uv lock`
- [x] Create `uv.lock` for reproducible builds
- [x] Pin BPTF to specific commit: `git+https://github.com/aschein/bptf.git@<sha>`
- [x] Remove `requirements.txt` (superseded by uv.lock)
- [x] Update README with uv install instructions

### 2. Add Makefile
Create a Makefile following ebbf pattern for common tasks.

- [x] `make dev` - Set up full dev environment (uv sync, pre-commit install)
- [x] `make test` - Run pytest
- [x] `make test-nf` - Run Nextflow integration tests
- [x] `make lint` - Run ruff
- [x] `make install-uv`, `make install-nextflow`, etc.

### 3. Add Pre-commit Hooks
- [x] Create `.pre-commit-config.yaml`:
  ```yaml
  repos:
    - repo: local
      hooks:
        - id: pytest
          name: pytest
          entry: uv run pytest
          language: system
          pass_filenames: false
          always_run: true
    - repo: https://github.com/astral-sh/ruff-pre-commit
      rev: v0.4.0
      hooks:
        - id: ruff
          args: [--fix]
        - id: ruff-format
  ```
- [x] Switch from black/flake8 to ruff (faster, unified)
- [x] Add `[tool.ruff]` config to pyproject.toml

### 4. Add Docker Support
- [x] Create `Dockerfile` (use python:3.11-slim base)
- [x] Include BPTF installation in image
- [x] Create `.dockerignore` (exclude .git, notebooks, results, data/*.h5ad)
- [x] Set up automated Docker builds in CI

---

## Phase 2: nf-core Pipeline Setup (High Priority)

### 5. Create nf-core Pipeline Skeleton
Set up Nextflow pipeline structure following nf-core conventions.

- [x] Create `nextflow/` directory structure (via `nf-core create`)
- [x] Use nf-schema plugin for parameter validation
- [x] Define resource profiles (local, docker, singularity, slurm)

### 6. Map Pipeline Steps to Nextflow Modules
The ALARMIST pipeline has 5 steps that need corresponding CLI commands and NF modules:

| Step | CLI Command | NF Module | Status |
|------|-------------|-----------|--------|
| 1. Patchify | `alarmist-patchify` | `ALARMIST_PATCHIFY` | Done |
| 2. BPTF | `alarmist-bptf` | `ALARMIST_BPTF` | Done |
| 3. Project | `alarmist-project` | `ALARMIST_PROJECT` | Done |
| 4. GLM | `alarmist-glm` | `ALARMIST_GLM` | Done |
| 5. Visualize | `alarmist-visualize` | `ALARMIST_VISUALIZE` | Done |

- [x] Complete missing CLI commands (project, glm, visualize)
- [x] Create NF module for each step
- [x] Wire up workflow in `alarmist.nf`

### 7. Add Cirro Integration (Optional)
- [ ] Create `.cirro/` directory with form schemas
- [ ] Define `process-form.json` for UI parameter input
- [ ] Add `process-input.json` and `process-output.json`

---

## Phase 3: CI/CD & Testing

### 8. GitHub Actions Workflows
Create workflows following ebbf pattern.

- [x] `.github/workflows/ci.yml` - Lint and test on push/PR
- [x] `.github/workflows/docker.yml` - Build and push container on tags
- [ ] `.github/workflows/bump_version.yml` - Automated versioning

### 9. Add Test Suite
- [x] Create `tests/` directory
- [x] Unit tests for core modules (lri.py, factorization.py, glm.py)
- [x] Integration test with minimal fixture data
- [x] Add small test AnnData file to `tests/fixtures/`
- [x] Configure pytest in pyproject.toml

---

## Phase 4: Code Cleanup

### 10. Fix Hardcoded Paths
- [x] Use `importlib.resources` for package data (LRI databases)
- [x] Remove hardcoded paths from `src/alarmist/core/lri.py`
- [x] Bundled LRI databases in `src/alarmist/config/lri_databases/`

### 11. Consolidate scripts/ and CLI
- [x] Port remaining logic from `scripts/0*.py` to CLI modules
- [x] Archive `scripts/` directory to `archived_scripts/` (excluded from git)
- [x] All CLI commands now implemented: patchify, bptf, project, glm, visualize

### 12. Documentation
- [x] Update README with nf-core pipeline usage
- [x] Add samplesheet examples in `samplesheets/`
- [x] Document all Nextflow parameters
- [ ] Clean up research notebooks vs tutorial notebooks

---

## Phase 5: Data & Reproducibility

### 13. Data Provenance
- [ ] Document LRI database versions and sources
- [ ] Add checksums for bundled data files
- [ ] Version track container images in `version.config`

### 14. Random Seed Audit
- [ ] Verify all stochastic operations respect `random_state`
- [ ] Document reproducibility guarantees
- [ ] Add reproducibility integration test

---

## Reference: ebbf Project Structure

Target structure after migration:
```
alarmist/
├── src/alarmist/           # Python package
│   ├── cli/                # CLI entry points
│   ├── core/               # Core analysis modules
│   ├── data/               # Data loading
│   └── plotting/           # Visualization
├── nextflow/               # nf-core pipeline
│   ├── main.nf
│   ├── nextflow.config
│   ├── workflows/
│   ├── modules/local/
│   └── conf/
├── tests/                  # pytest tests
├── .github/workflows/      # CI/CD
├── .cirro/                 # Cirro integration (optional)
├── Dockerfile
├── Makefile
├── pyproject.toml
├── uv.lock
└── .pre-commit-config.yaml
```

---

*Generated: 2026-03-09*
