#!/usr/bin/env python3
"""
Cirro preprocess script for ALARMIST pipeline.

This script prepares the samplesheet.csv input file required by the
nf-core/alarmist Nextflow pipeline.

The ALARMIST pipeline expects a samplesheet with:
- sample_id: Unique identifier for each sample
- adata_path: Path to the AnnData (.h5ad) file
"""

from __future__ import annotations

from urllib.parse import urlparse

import numpy as np
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

SAMPLESHEET_REQUIRED_COLUMNS = (
    "sample_id",
    "adata_path",
)


def is_url(string: str) -> bool:
    """Check if a string is a URL."""
    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


def samplesheet_from_files(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Create a samplesheet from Cirro's files DataFrame.

    Uses ds.files to get h5ad files and creates adata_path from each file.
    Merges with ds.samplesheet to get sample metadata.
    """
    files = ds.files

    ds.logger.info(f"Found files in ds.files: {files.to_dict()}")

    # Filter for .h5ad files only
    files = files.copy()
    files = files[files["file"].str.endswith(".h5ad")]

    if files.empty:
        ds.logger.warning("No .h5ad files found in dataset files")
        return pd.DataFrame()

    # Use the file path directly as adata_path
    # For s3 links, ensure proper formatting
    files["adata_path"] = files["file"].apply(
        lambda x: (
            x.replace("s3:/", "s3://")
            if x.startswith("s3:/") and not x.startswith("s3://")
            else x
        )
    )

    # Rename sample column to sample_id
    files = files.rename(columns={"sample": "sample_id"})
    files = files[["sample_id", "adata_path"]]

    # Deduplicate - keep one row per sample
    files = files.drop_duplicates(subset=["sample_id"])

    # Merge with samplesheet if it has additional metadata
    if not ds.samplesheet.empty and "sample" in ds.samplesheet.columns:
        ds.samplesheet = ds.samplesheet.rename(columns={"sample": "sample_id"})
        samplesheet = pd.merge(files, ds.samplesheet, on="sample_id", how="left")
    else:
        samplesheet = files

    return samplesheet


def samplesheet_from_params(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Create a samplesheet from dataset metadata inputs.

    Falls back to using ds.metadata['inputs'] when no files are found.
    """
    if not ds.metadata.get("inputs"):
        return pd.DataFrame()

    data_params = pd.DataFrame(
        {
            "sample_id": [x["name"] for x in ds.metadata["inputs"]],
            "adata_path": [x["dataPath"] for x in ds.metadata["inputs"]],
        }
    )

    return data_params


def prepare_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Prepare the samplesheet for the pipeline.

    Tries to create from files first, falls back to params if no files found.
    Ensures all required columns are present and cleans up params.
    """
    ds.logger.info(f"Params: {ds.params}")

    samplesheet = samplesheet_from_files(ds)

    # Check if pipeline uses Cirro samplesheet, and if not prepare it from params
    if samplesheet.empty:
        ds.logger.warning(
            "No files found in dataset. Preparing samplesheet from params."
        )
        samplesheet = samplesheet_from_params(ds)
        if samplesheet.empty:
            raise ValueError(
                "No files found in dataset and unable to prepare "
                "samplesheet from params."
            )
        ds.logger.info("Prepared samplesheet from params.")

    # Ensure all required columns are present (populate missing)
    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in samplesheet.columns:
            ds.logger.warning(
                f"Samplesheet is missing required column '{colname}'."
                " Populating with NaN."
            )
            samplesheet[colname] = np.nan

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=False)

    # Clear params that we wrote to the samplesheet
    # cleared params will not overload the nextflow.params
    to_remove = []
    for k in ds.params:
        if k in SAMPLESHEET_REQUIRED_COLUMNS:
            ds.logger.info(
                f"Removing param '{k}' from dataset params"
                " as it is now in the samplesheet."
            )
            to_remove.append(k)

    for k in to_remove:
        ds.remove_param(k)

    ds.add_param("input", "samplesheet.csv", overwrite=True)

    # Log the samplesheet
    ds.logger.info(f"Samplesheet: {samplesheet.to_dict()}")

    return samplesheet


def main():
    ds = PreprocessDataset.from_running()

    ds.logger.info("Creating samplesheet from input files")
    ds.logger.info(f"Input files: {len(ds.files)} rows")
    ds.logger.info(f"Input columns: {list(ds.files.columns)}")
    ds.logger.info(f"Files DataFrame:\n{ds.files.to_string()}")

    samplesheet = prepare_samplesheet(ds)
    ds.logger.info(f"Samplesheet created with {len(samplesheet)} samples")

    # Log final params
    ds.logger.info(f"Final params: {ds.params}")


if __name__ == "__main__":
    main()
