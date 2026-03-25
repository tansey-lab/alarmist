#!/bin/bash
set -euo pipefail

# Run ALARMIST pipeline on test data
# Usage: ./scripts/run_test.sh [h5ad_path]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default path - override with first argument if provided
H5AD_PATH="${1:-${PROJECT_DIR}/tests/fixtures/ovarian_cancer_sample.h5ad}"
# Convert to absolute path
H5AD_PATH="$(cd "$(dirname "$H5AD_PATH")" && pwd)/$(basename "$H5AD_PATH")"

# Output directories
OUTDIR="${PROJECT_DIR}/results/test"
WORKDIR="${PROJECT_DIR}/work"

# Create samplesheet
SAMPLESHEET="${PROJECT_DIR}/samplesheets/test.csv"
mkdir -p "$(dirname "$SAMPLESHEET")"

cat > "$SAMPLESHEET" << EOF
sample_id,adata_path
ovarian_sample,${H5AD_PATH}
EOF

echo "Created samplesheet: $SAMPLESHEET"
echo "Using h5ad: $H5AD_PATH"
echo "Output directory: $OUTDIR"

# Run Nextflow pipeline
cd "${PROJECT_DIR}/nextflow"

nextflow run main.nf \
    -profile test \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    -work-dir "$WORKDIR" \
    -resume

echo "Pipeline complete. Results in: $OUTDIR"
