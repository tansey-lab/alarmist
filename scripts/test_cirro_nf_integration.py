#!/usr/bin/env python3
"""
Integration test for Cirro config <> Nextflow pipeline.

This script simulates how Cirro invokes the pipeline:
1. Creates test input data that mimics a Cirro dataset
2. Runs the preprocessing script logic
3. Invokes Nextflow with the generated samplesheet
4. Validates the entire flow works correctly

Usage:
    python scripts/test_cirro_nf_integration.py [--dry-run] [--keep-temp]
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / ".cirro"))


class MockLogger:
    """Mock logger that mimics Cirro's PreprocessDataset logger."""

    def info(self, msg: str):
        print(f"[INFO] {msg}")

    def warning(self, msg: str):
        print(f"[WARN] {msg}")

    def error(self, msg: str):
        print(f"[ERROR] {msg}")


@dataclass
class MockPreprocessDataset:
    """
    Mock of Cirro's PreprocessDataset for local testing.

    Simulates the interface that preprocess.py expects from Cirro.
    """

    files: pd.DataFrame
    samplesheet: pd.DataFrame = field(default_factory=pd.DataFrame)
    params: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=lambda: {"inputs": []})
    logger: MockLogger = field(default_factory=MockLogger)
    output_params: dict = field(default_factory=dict)

    def add_param(self, key: str, value: str, overwrite: bool = False):
        """Register a parameter to pass to Nextflow."""
        self.logger.info(f"add_param({key}={value})")
        self.output_params[key] = value

    def remove_param(self, key: str):
        """Remove a parameter."""
        self.logger.info(f"remove_param({key})")
        if key in self.params:
            del self.params[key]


def create_test_data(
    temp_dir: Path, fixture_h5ad: Path | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create test data that mimics what Cirro would provide.

    Creates AnnData (.h5ad) file references that the ALARMIST pipeline expects.

    Args:
        temp_dir: Directory to create test files in
        fixture_h5ad: Optional path to a real h5ad file to use as fixture

    Returns:
        Tuple of (files DataFrame, samplesheet DataFrame) matching Cirro's format
    """
    files_data = []
    sample_names = []

    # Create test directory structure
    inputs_dir = temp_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    for sample_name in ["sample_A", "sample_B"]:
        sample_names.append(sample_name)
        sample_dir = inputs_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Create or link h5ad file
        h5ad_path = sample_dir / f"{sample_name}.h5ad"

        if fixture_h5ad and fixture_h5ad.exists():
            # Symlink to fixture file for realistic testing
            h5ad_path.symlink_to(fixture_h5ad.resolve())
        else:
            # Create placeholder file
            h5ad_path.touch()

        # Each file gets its own row in Cirro's files DataFrame
        files_data.append({"sample": sample_name, "file": str(h5ad_path)})

    files_df = pd.DataFrame(files_data)
    samplesheet_df = pd.DataFrame({"sample": sample_names})

    return files_df, samplesheet_df


def run_preprocess(ds: MockPreprocessDataset, work_dir: Path) -> Path:
    """
    Run the preprocessing logic from .cirro/preprocess.py.

    Args:
        ds: Mock dataset with files to process
        work_dir: Directory to write samplesheet to

    Returns:
        Path to generated samplesheet
    """
    import os

    original_dir = os.getcwd()
    os.chdir(work_dir)

    try:
        # Import the preprocessing functions
        from preprocess import prepare_samplesheet

        ds.logger.info("Creating samplesheet from input files")
        ds.logger.info(f"Input files: {len(ds.files)} rows")
        ds.logger.info(f"Input columns: {list(ds.files.columns)}")

        samplesheet = prepare_samplesheet(ds)
        ds.logger.info(f"Samplesheet created with {len(samplesheet)} samples")

        samplesheet_path = work_dir / "samplesheet.csv"
        return samplesheet_path
    finally:
        os.chdir(original_dir)


def validate_samplesheet_against_schema(samplesheet_path: Path) -> list[str]:
    """
    Validate the generated samplesheet against the Nextflow schema.

    Returns:
        List of validation errors (empty if valid)
    """
    schema_path = PROJECT_ROOT / "nextflow" / "assets" / "schema_input.json"

    with open(schema_path) as f:
        schema = json.load(f)

    samplesheet = pd.read_csv(samplesheet_path)
    errors = []

    # Check required columns from schema
    required_props = schema["items"].get("required", [])
    for col in required_props:
        if col not in samplesheet.columns:
            errors.append(f"Missing required column: {col}")

    # Check for unexpected columns
    expected_cols = set(schema["items"]["properties"].keys())
    actual_cols = set(samplesheet.columns)
    extra_cols = actual_cols - expected_cols
    if extra_cols:
        errors.append(f"Unexpected columns in samplesheet: {extra_cols}")

    missing_cols = expected_cols - actual_cols
    if missing_cols:
        errors.append(f"Missing expected columns: {missing_cols}")

    return errors


def run_nextflow_pipeline(
    samplesheet_path: Path,
    output_dir: Path,
    dry_run: bool = False,
    profile: str = "dev",
) -> bool:
    """
    Run the Nextflow pipeline with the generated samplesheet.

    Args:
        samplesheet_path: Path to the samplesheet CSV
        output_dir: Directory for pipeline outputs
        dry_run: If True, only print the command without executing
        profile: Nextflow profile to use

    Returns:
        True if pipeline succeeded, False otherwise
    """
    nextflow_dir = PROJECT_ROOT / "nextflow"

    cmd = [
        "nextflow",
        "run",
        "main.nf",
        "-profile",
        profile,
        "--input",
        str(samplesheet_path),
        "--outdir",
        str(output_dir),
    ]

    print(f"\n{'=' * 60}")
    print("Running Nextflow pipeline:")
    print(f"  Command: {' '.join(cmd)}")
    print(f"  Working dir: {nextflow_dir}")
    print(f"{'=' * 60}\n")

    if dry_run:
        print("[DRY RUN] Would execute the above command")
        return True

    result = subprocess.run(cmd, cwd=nextflow_dir)
    return result.returncode == 0


def print_summary(
    samplesheet_path: Path,
    schema_errors: list[str],
    pipeline_success: bool | None,
):
    """Print a summary of the integration test results."""
    print(f"\n{'=' * 60}")
    print("INTEGRATION TEST SUMMARY")
    print(f"{'=' * 60}")
    print(f"Samplesheet: {samplesheet_path}")

    # Show samplesheet contents
    print("\nSamplesheet contents:")
    print(pd.read_csv(samplesheet_path).to_string(index=False))

    # Schema validation
    print("\nSchema validation: ", end="")
    if schema_errors:
        print("FAILED")
        for error in schema_errors:
            print(f"  - {error}")
    else:
        print("PASSED")

    # Pipeline execution
    if pipeline_success is not None:
        print(f"\nPipeline execution: {'PASSED' if pipeline_success else 'FAILED'}")

    print(f"{'=' * 60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Integration test for Cirro config <> Nextflow pipeline"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print Nextflow command without executing",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary directories after test",
    )
    parser.add_argument(
        "--skip-pipeline",
        action="store_true",
        help="Skip running the Nextflow pipeline (only test preprocessing)",
    )
    parser.add_argument(
        "--fixture-h5ad",
        type=Path,
        help="Path to a real h5ad file to use as test fixture",
    )
    args = parser.parse_args()

    # Check if preprocess.py exists
    preprocess_path = PROJECT_ROOT / ".cirro" / "preprocess.py"
    if not preprocess_path.exists():
        print(f"ERROR: {preprocess_path} does not exist")
        print("Create a preprocess.py in .cirro/ first")
        return 1

    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="cirro-integration-test-"))
    print(f"Working directory: {temp_dir}")

    try:
        # Step 1: Create test data
        print("\n--- Creating test data (h5ad files) ---")
        files_df, samplesheet_df = create_test_data(temp_dir, args.fixture_h5ad)
        print(f"Created {len(files_df)} test files for {len(samplesheet_df)} samples")

        # Step 2: Run preprocessing
        print("\n--- Running Cirro preprocessing ---")
        ds = MockPreprocessDataset(files=files_df, samplesheet=samplesheet_df)

        try:
            samplesheet_path = run_preprocess(ds, temp_dir)
            print(f"Generated samplesheet: {samplesheet_path}")
        except Exception as e:
            print(f"Preprocessing failed: {e}")
            import traceback

            traceback.print_exc()
            return 1

        # Step 3: Validate against schema
        print("\n--- Validating samplesheet against Nextflow schema ---")
        schema_errors = validate_samplesheet_against_schema(samplesheet_path)

        # Step 4: Run Nextflow pipeline (optional)
        pipeline_success = None
        if not args.skip_pipeline and not schema_errors:
            print("\n--- Running Nextflow pipeline ---")
            output_dir = temp_dir / "output"
            output_dir.mkdir()
            pipeline_success = run_nextflow_pipeline(
                samplesheet_path, output_dir, dry_run=args.dry_run
            )
        elif schema_errors:
            print("\n--- Skipping Nextflow pipeline due to schema errors ---")

        # Print summary
        print_summary(samplesheet_path, schema_errors, pipeline_success)

        # Return appropriate exit code
        if schema_errors:
            return 1
        if pipeline_success is False:
            return 1
        return 0

    finally:
        if args.keep_temp:
            print(f"\nTemp directory preserved: {temp_dir}")
        else:
            shutil.rmtree(temp_dir)
            print(f"\nCleaned up temp directory: {temp_dir}")


if __name__ == "__main__":
    sys.exit(main())
