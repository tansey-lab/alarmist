# Notebooks

This directory contains Jupyter notebooks for ALARMIST.

## Structure

```
notebooks/
├── research/           # Research and analysis notebooks (not for distribution)
│   ├── 0112_*.ipynb   # Date-prefixed analysis notebooks
│   ├── 1203_*.ipynb   # Figure generation notebooks
│   └── ...
└── README.md
```

## Tutorials

For user-facing tutorials and examples, see the [`tutorials/`](../tutorials/) directory:

- **GBM.ipynb** - Complete walkthrough of the ALARMIST pipeline on glioblastoma data

## Research Notebooks

The `research/` subdirectory contains internal analysis notebooks used during development.
These notebooks may contain hardcoded paths and are not intended for external use.

To run the pipeline, use either:
1. The CLI commands (`alarmist-patchify`, `alarmist-bptf`, etc.)
2. The Nextflow pipeline (`nextflow run main.nf`)
3. The tutorial notebooks in `tutorials/`
