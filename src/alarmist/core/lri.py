"""
Patch-based Ligand-Receptor Interaction Analysis

This module implements spatial patch-based LRI analysis for matrix factorization approaches.
It divides tissue into regular grid patches and counts all-to-all interactions within each patch.
"""

import logging
import os
from abc import ABC
from importlib import resources
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

# from numba import njit, prange
import scipy.sparse as sp
import scipy.sparse as sparse
from liana.resource import select_resource
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse import vstack as sparse_vstack
from sklearn.neighbors import KDTree
from tqdm import tqdm

logger = logging.getLogger(__name__)


def _get_bundled_database_path(database_name: str) -> Path:
    """
    Get the path to a bundled LRI database file using importlib.resources.

    Parameters
    ----------
    database_name : str
        Name of the database file (e.g., 'CellChatDBv2.0.human.csv')

    Returns
    -------
    Path
        Path to the database file
    """
    try:
        # Python 3.9+ approach
        ref = resources.files("alarmist.config.lri_databases").joinpath(database_name)
        # For resources that need to be extracted to filesystem
        with resources.as_file(ref) as path:
            return Path(path)
    except (AttributeError, TypeError):
        # Fallback for older Python versions
        import pkg_resources

        return Path(
            pkg_resources.resource_filename(
                "alarmist", f"config/lri_databases/{database_name}"
            )
        )


def _split_gene_complex(gene_str: str) -> list[str]:
    """
    Split a gene complex string into individual gene names.
    Supports both '_' and ',' as separators.

    Parameters
    ----------
    gene_str : str
        Gene complex string (e.g., "A_B", "A,B", "A_B,C")

    Returns
    -------
    List[str]
        List of individual gene names

    Examples
    --------
    >>> _split_gene_complex("IL12A_IL12B")
    ['IL12A', 'IL12B']
    >>> _split_gene_complex("IL12A,IL12B")
    ['IL12A', 'IL12B']
    >>> _split_gene_complex("A_B,C")
    ['A', 'B', 'C']
    """
    # Replace all commas with underscores, then split
    normalized = gene_str.replace(",", "_")
    return normalized.split("_")


def load_database_genes(
    resource_name: str,
    cellchatdb_path: str | None = None,
    cellphonedb_path: str | None = None,
) -> set[str]:
    """
    Return the full set of ligand/receptor gene names referenced by an LRI
    database, without filtering against any input dataset. Used for assessing
    overall database coverage.
    """
    if resource_name.lower() == "cellchatdb":
        if cellchatdb_path is None:
            cellchatdb_path = str(
                _get_bundled_database_path("CellChatDBv2.0.human.csv")
            )
        resource = pd.read_csv(cellchatdb_path)
    elif resource_name.lower() == "cellphonedb":
        if cellphonedb_path is None:
            cellphonedb_path = str(
                _get_bundled_database_path("CellPhoneDBv5.0.human.csv")
            )
        resource = pd.read_csv(cellphonedb_path)
    else:
        resource = select_resource(resource_name)

    genes: set[str] = set()
    for col in ("ligand", "receptor"):
        if col not in resource.columns:
            continue
        for val in resource[col].dropna():
            genes.update(_split_gene_complex(str(val)))
    genes.discard("")
    return genes


class BaseLRIAnalyzer(ABC):
    """
    Base class for Ligand-Receptor Interaction Analyzers

    Provides common functionality for LRI analysis including database loading,
    column structure creation, and support for ligand/receptor complexes.
    """

    def __init__(
        self,
        resource_name: str = "cellchatdb",
        spliter: str = "|",
        cellchatdb_path: str | None = None,
        cellphonedb_path: str | None = None,
        cell_type_column: str = "cell_type",
    ):
        """
        Initialize the base LRI analyzer.

        Parameters
        ----------
        resource_name : str, default 'cellchatdb'
            LRI database to use from liana
        spliter : str, default '|'
            Separator for column names
        cellchatdb_path : str, optional
            Path to local CellChatDB CSV file. If None, uses bundled database.
        cellphonedb_path : str, optional
            Path to local CellPhoneDB CSV file. If None, uses bundled database.
        cell_type_column : str
            Column name for cell types in adata.obs
        """
        self.resource_name = resource_name
        self.spliter = spliter
        self.cell_type_column = cell_type_column

        # Use bundled databases if paths not provided
        if cellchatdb_path is None:
            self.cellchatdb_path = str(
                _get_bundled_database_path("CellChatDBv2.0.human.csv")
            )
        else:
            self.cellchatdb_path = cellchatdb_path

        if cellphonedb_path is None:
            self.cellphonedb_path = str(
                _get_bundled_database_path("CellPhoneDBv5.0.human.csv")
            )
        else:
            self.cellphonedb_path = cellphonedb_path

        # Analysis results (to be populated)
        self.lr_pairs = None
        self.ligand_genes_list = None
        self.receptor_genes_list = None
        self.signaling_types = None
        self.cell_types = None
        self.column_names = None

    def prepare_lri_database(
        self, adata: anndata.AnnData | None = None, gene_names: list[str] | None = None
    ) -> tuple[list[tuple], list[list[str]], list[list[str]], list[str]]:
        """
        Prepare ligand-receptor pairs from database.

        Parameters
        ----------
        adata : anndata.AnnData, optional
            Spatial transcriptomics data. Used to get gene names if gene_names is not provided.
        gene_names : List[str], optional
            List of gene names to filter LR pairs. If provided, overrides adata.var_names.
            This is useful for multi-sample analysis where we want to use the intersection
            of genes across all samples.

        Returns
        -------
        lr_pairs : List[Tuple]
            List of (ligand_str, receptor_str) tuples
        ligand_genes_list : List[List[str]]
            List of ligand gene lists for each LR pair
        receptor_genes_list : List[List[str]]
            List of receptor gene lists for each LR pair
        signaling_types : List[str]
            List of signaling types for each LR pair
        """
        # Determine which gene set to use for filtering
        if gene_names is not None:
            available_genes = set(gene_names)
        elif adata is not None:
            available_genes = set(adata.var_names)
        else:
            raise ValueError("Either adata or gene_names must be provided")

        logger.debug(f"Loading {self.resource_name} database...")
        logger.debug(f"Filtering with {len(available_genes)} available genes")

        # Load from local CSV if cellchatdb or cellphonedb
        if self.resource_name.lower() == "cellchatdb":
            resource = pd.read_csv(self.cellchatdb_path)
            # Check required columns
            required_cols = ["ligand", "receptor", "signaling_type"]
            if not all(col in resource.columns for col in required_cols):
                raise ValueError(
                    f"CellChatDB CSV must contain {required_cols}. "
                    f"Found columns: {resource.columns.tolist()}"
                )
        elif self.resource_name.lower() == "cellphonedb":
            resource = pd.read_csv(self.cellphonedb_path)
            # Check required columns
            required_cols = ["ligand", "receptor", "signaling_type"]
            if not all(col in resource.columns for col in required_cols):
                raise ValueError(
                    f"CellPhoneDB CSV must contain {required_cols}. "
                    f"Found columns: {resource.columns.tolist()}"
                )
        else:
            # Use liana's select_resource for other databases
            resource = select_resource(self.resource_name)
            # LIANA doesn't have signaling_type, add it as 'Unknown'
            if "signaling_type" not in resource.columns:
                resource["signaling_type"] = "Unknown"

        # Filter pairs where ALL ligand genes and ALL receptor genes exist
        lr_pairs = []
        ligand_genes_list = []
        receptor_genes_list = []
        signaling_types = []

        for idx in range(len(resource)):
            ligand = resource.iloc[idx]["ligand"]
            receptor_str = resource.iloc[idx]["receptor"]
            signaling_type = resource.iloc[idx]["signaling_type"]

            # Skip if ligand or receptor is NaN
            if pd.isna(ligand) or pd.isna(receptor_str):
                continue

            # Convert to string
            ligand_str = str(ligand)
            receptor_str = str(receptor_str)
            signaling_type = (
                str(signaling_type) if not pd.isna(signaling_type) else "Unknown"
            )

            # Parse ligand and receptor genes (supports both '_' and ',' separators)
            ligand_genes = _split_gene_complex(ligand_str)
            receptor_genes = _split_gene_complex(receptor_str)

            # Check if ALL ligand genes and ALL receptor genes exist in available genes
            if all(lig_gene in available_genes for lig_gene in ligand_genes) and all(
                rec_gene in available_genes for rec_gene in receptor_genes
            ):
                lr_pairs.append((ligand_str, receptor_str))
                ligand_genes_list.append(ligand_genes)
                receptor_genes_list.append(receptor_genes)
                signaling_types.append(signaling_type)

        logger.debug(f"Initial L-R pairs in data: {len(lr_pairs)}")
        logger.debug(
            f"  Single ligand: {sum(1 for lg in ligand_genes_list if len(lg) == 1)}"
        )
        logger.debug(
            f"  Multi ligand: {sum(1 for lg in ligand_genes_list if len(lg) > 1)}"
        )
        logger.debug(
            f"  Single receptor: {sum(1 for rg in receptor_genes_list if len(rg) == 1)}"
        )
        logger.debug(
            f"  Multi receptor: {sum(1 for rg in receptor_genes_list if len(rg) > 1)}"
        )

        # Print signaling type distribution
        from collections import Counter

        sig_type_counts = Counter(signaling_types)
        logger.debug("Signaling type distribution:")
        for sig_type, count in sig_type_counts.items():
            logger.debug(f"  {sig_type}: {count}")

        self.lr_pairs = lr_pairs
        self.ligand_genes_list = ligand_genes_list
        self.receptor_genes_list = receptor_genes_list
        self.signaling_types = signaling_types

        return lr_pairs, ligand_genes_list, receptor_genes_list, signaling_types

    def create_column_structure(
        self, adata: anndata.AnnData, signaling_types: list[str]
    ) -> list[str]:
        """
        Create column names for the LRI matrix.

        For Cell-Cell Contact: one 'juxtacrine' column per cell type pair
        For other types: 'autocrine' and 'paracrine' columns (autocrine only for same cell type)

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
        signaling_types : List[str]
            Signaling type for each LR pair (same order as lr_pairs)

        Returns
        -------
        column_names : List[str]
            List of column names in format: ligand_ct|receptor_ct|ligand|receptor|mode
        """
        # Handle both categorical and non-categorical columns
        cell_type_col = adata.obs[self.cell_type_column]
        if hasattr(cell_type_col, "cat"):
            self.cell_types = cell_type_col.cat.categories.tolist()
        else:
            self.cell_types = sorted(cell_type_col.unique().tolist())
        column_names = []

        for idx, (lig, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]

            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == "Cell-Cell Contact":
                        # Cell-Cell Contact: only juxtacrine (no autocrine/paracrine split)
                        column_names.append(
                            f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}juxtacrine"
                        )
                    else:
                        # Non-contact signaling: autocrine + paracrine (or just paracrine)
                        if lig_ct == rec_ct:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}autocrine"
                            )
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )
                        else:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )

        self.column_names = column_names
        return column_names

    def remove_zero_columns(
        self, matrix: csr_matrix, column_names: list[str]
    ) -> tuple[csr_matrix, list[str]]:
        """
        Remove columns that are all zeros from the matrix and corresponding column names.

        Parameters
        ----------
        matrix : csr_matrix
            Sparse matrix to filter
        column_names : List[str]
            Column names corresponding to matrix columns

        Returns
        -------
        filtered_matrix : csr_matrix
            Matrix with zero columns removed
        filtered_column_names : List[str]
            Column names with zero columns removed
        """
        # Calculate sum for each column
        col_sums = np.array(matrix.sum(axis=0)).ravel()

        # Find non-zero columns
        nonzero_cols = np.where(col_sums > 0)[0]

        n_total = len(column_names)
        n_nonzero = len(nonzero_cols)
        n_removed = n_total - n_nonzero

        logger.debug("Filtering zero columns:")
        logger.debug(f"  Total columns: {n_total}")
        logger.debug(f"  Non-zero columns: {n_nonzero}")
        logger.debug(
            f"  Zero columns removed: {n_removed} ({n_removed / n_total * 100:.1f}%)"
        )

        # Filter matrix and column names
        filtered_matrix = matrix[:, nonzero_cols]
        filtered_column_names = [column_names[i] for i in nonzero_cols]

        # Update instance variable
        self.column_names = filtered_column_names

        return filtered_matrix, filtered_column_names

    def _split_adata_by_sample(
        self, adata: anndata.AnnData, sample_column: str
    ) -> dict[str, anndata.AnnData]:
        """
        Split a merged AnnData object into a dictionary of AnnData objects by sample.

        Parameters
        ----------
        adata : anndata.AnnData
            Merged AnnData object containing multiple samples
        sample_column : str
            Column name in adata.obs that identifies different samples

        Returns
        -------
        adata_dict : Dict[str, anndata.AnnData]
            Dictionary mapping sample_id -> AnnData subset
        """
        logger.debug(f"Splitting merged AnnData by '{sample_column}'...")

        # Get unique sample IDs
        sample_ids = adata.obs[sample_column].unique()

        # Handle categorical vs non-categorical
        if hasattr(sample_ids, "categories"):
            sample_ids = sample_ids.categories.tolist()
        else:
            sample_ids = sorted([s for s in sample_ids if pd.notna(s)])

        logger.debug(f"Found {len(sample_ids)} samples: {sample_ids}")

        adata_dict = {}
        for sample_id in sample_ids:
            # Subset adata for this sample
            mask = adata.obs[sample_column] == sample_id
            adata_subset = adata[mask].copy()

            # Ensure cell_type column is categorical with only present categories
            if self.cell_type_column in adata_subset.obs.columns:
                ct_col = adata_subset.obs[self.cell_type_column]
                if hasattr(ct_col, "cat"):
                    # Remove unused categories
                    adata_subset.obs[self.cell_type_column] = (
                        ct_col.cat.remove_unused_categories()
                    )

            adata_dict[str(sample_id)] = adata_subset
            logger.debug(
                f"  {sample_id}: {adata_subset.n_obs} cells, "
                f"{adata_subset.obs[self.cell_type_column].nunique()} cell types"
            )

        return adata_dict


class PatchLRIAnalyzer(BaseLRIAnalyzer):
    """
    Patch-based Ligand-Receptor Interaction Analyzer

    Divides spatial transcriptomics data into regular grid patches and counts
    ligand-receptor interactions within each patch using all-to-all strategy.
    """

    def __init__(
        self,
        patch_size: float = 50.0,
        resource_name: str = "cellchatdb",
        spliter: str = "|",
        cellchatdb_path: str | None = None,
        cellphonedb_path: str | None = None,
        cell_type_column: str = "cell_type",
    ):
        """
        Initialize the patch-based LRI analyzer.

        Parameters
        ----------
        patch_size : float, default 50.0
            Size of spatial patches in micrometers
        resource_name : str, default 'cellchatdb'
            LRI database to use from liana
        spliter : str, default '|'
            Separator for column names
        cellchatdb_path : str, optional
            Path to local CellChatDB CSV file. If None, uses bundled database.
        cellphonedb_path : str, optional
            Path to local CellPhoneDB CSV file. If None, uses bundled database.
        cell_type_column : str
            Column name for cell types in adata.obs
        """
        # Initialize base class
        super().__init__(
            resource_name=resource_name,
            spliter=spliter,
            cellchatdb_path=cellchatdb_path,
            cellphonedb_path=cellphonedb_path,
            cell_type_column=cell_type_column,
        )

        # Patch-specific attributes
        self.patch_size = patch_size
        self.patch_assignments = None
        self.patch_coords = None
        self.patch_lri_matrix = None

    def create_spatial_patches(self, adata: anndata.AnnData) -> tuple[np.ndarray, dict]:
        """
        Divide tissue into regular grid patches.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data

        Returns
        -------
        patch_assignments : np.ndarray
            Array mapping each cell to patch_id
        patch_info : dict
            Dictionary containing patch grid information
        """
        coords = adata.obsm["spatial"][:, :2]

        # Define grid boundaries with padding
        x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
        y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

        # Create grid bins
        x_bins = np.arange(x_min, x_max + self.patch_size, self.patch_size)
        y_bins = np.arange(y_min, y_max + self.patch_size, self.patch_size)

        # Assign cells to patches
        x_indices = np.digitize(coords[:, 0], x_bins) - 1
        y_indices = np.digitize(coords[:, 1], y_bins) - 1

        # Create unique patch IDs
        patch_assignments = x_indices * len(y_bins) + y_indices

        # Create patch coordinate mapping
        patch_coords = {}
        for patch_id in np.unique(patch_assignments):
            x_idx = patch_id // len(y_bins)
            y_idx = patch_id % len(y_bins)

            # Calculate patch center coordinates
            if x_idx < len(x_bins) - 1 and y_idx < len(y_bins) - 1:
                center_x = (x_bins[x_idx] + x_bins[x_idx + 1]) / 2
                center_y = (y_bins[y_idx] + y_bins[y_idx + 1]) / 2
                patch_coords[patch_id] = (center_x, center_y)

        patch_info = {
            "x_bins": x_bins,
            "y_bins": y_bins,
            "patch_coords": patch_coords,
            "n_patches": len(patch_coords),
        }

        self.patch_assignments = patch_assignments
        self.patch_coords = patch_coords

        return patch_assignments, patch_info

    def build_patch_lri_matrix(
        self, adata: anndata.AnnData, signaling_types: list[str]
    ) -> csr_matrix:
        """
        Build the patch-LRI interaction matrix, distinguishing autocrine from paracrine.
        Supports both ligand and receptor complexes (all genes must be co-expressed).
        """
        logger.debug(
            "Building patch-LRI matrix with autocrine/paracrine distinction..."
        )

        # ─── 1) Prepare basics ────────────────────────────────────────────────────────
        unique_patches = np.array(
            [p for p in np.unique(self.patch_assignments) if p in self.patch_coords]
        )
        patch_idx_map = {pid: i for i, pid in enumerate(unique_patches)}
        n_patches = len(unique_patches)
        n_columns = len(self.column_names)
        logger.debug(f"Processing {n_patches} patches × {n_columns} LRI combinations")

        # ─── 2) Index mappings ───────────────────────────────────────────────────────
        ct_to_idx = {ct: i for i, ct in enumerate(self.cell_types)}

        # Check for nan/missing cell types and map them
        cell_type_values = adata.obs[self.cell_type_column].values
        cell_types_idx = []
        for ct in cell_type_values:
            if pd.isna(ct):
                # Assign -1 for nan cell types (will be excluded from analysis)
                cell_types_idx.append(-1)
            elif ct not in ct_to_idx:
                # Cell type not in shared cell types (shouldn't happen in multi-sample mode)
                cell_types_idx.append(-1)
            else:
                cell_types_idx.append(ct_to_idx[ct])
        cell_types_idx = np.array(cell_types_idx, dtype=int)

        # Report cells with missing/invalid cell types
        n_invalid = np.sum(cell_types_idx == -1)
        if n_invalid > 0:
            logger.debug(
                f"  Warning: {n_invalid} cells ({n_invalid / len(cell_types_idx) * 100:.1f}%) have missing/invalid cell types and will be excluded"
            )

        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # Parse column metadata - NOW HANDLES LIGAND COMPLEXES
        col_meta = []
        for idx, (lig_str, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]

            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == "Cell-Cell Contact":
                        # One column: juxtacrine
                        col_meta.append(
                            (
                                len(col_meta),  # Column index
                                ct_to_idx[lig_ct],
                                ct_to_idx[rec_ct],
                                lig_str,  # Store ligand string (may be complex)
                                rec_str,
                                "juxtacrine",
                                sig_type,
                            )
                        )
                    else:
                        # Non-contact
                        if lig_ct == rec_ct:
                            # Two columns: autocrine + paracrine
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "autocrine",
                                    sig_type,
                                )
                            )
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "paracrine",
                                    sig_type,
                                )
                            )
                        else:
                            # One column: paracrine
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "paracrine",
                                    sig_type,
                                )
                            )

        # ─── 3) Binarize expression ───────────────────────────────────────────────────
        X = adata.X
        if sp.issparse(X):
            expr_bool = (X > 0).astype(int).tocsc()
        else:
            expr_bool = csr_matrix((X > 0).astype(int)).tocsc()

        expr_coo = expr_bool.tocoo()

        # ─── 4) Build patch_by_lig for INDIVIDUAL ligand genes ────────────────────────
        patch_by_lig = {}
        logger.debug("Building patch-by-ligand matrices (individual genes)...")
        for ct_idx in range(len(self.cell_types)):
            mask_cells = cell_types_idx == ct_idx
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (
                    expr_coo.data[entry_mask],
                    (expr_coo.row[entry_mask], expr_coo.col[entry_mask]),
                ),
                shape=expr_bool.shape,
            )

            # Collect all individual ligand genes for this cell type (expand complexes)
            all_lig_genes_ct = set()
            for _, lct, _, lig_str, _, _, _ in col_meta:
                if lct == ct_idx:
                    all_lig_genes_ct.update(_split_gene_complex(lig_str))
            all_lig_genes_ct = sorted([g for g in all_lig_genes_ct if g in gene_to_idx])

            # Build matrix for each individual ligand gene
            lig_gene_matrices = {}
            for lig_gene in all_lig_genes_ct:
                lig_gene_idx = gene_to_idx[lig_gene]
                lig_mask = sub.col == lig_gene_idx
                rows_cells = sub.row[lig_mask]
                data_vals_ct = sub.data[lig_mask]
                patch_rows = np.array(
                    [patch_idx_map[self.patch_assignments[c]] for c in rows_cells],
                    dtype=int,
                )

                lig_coo = coo_matrix(
                    (data_vals_ct, (patch_rows, np.zeros(len(patch_rows), dtype=int))),
                    shape=(n_patches, 1),
                )
                lig_coo.sum_duplicates()
                lig_gene_matrices[lig_gene] = lig_coo.tocsr()

            patch_by_lig[ct_idx] = lig_gene_matrices

        # ─── 5) Build patch_by_rec for INDIVIDUAL receptor genes ──────────────────────
        patch_by_rec = {}
        logger.debug("Building patch-by-receptor matrices (individual genes)...")
        for ct_idx in range(len(self.cell_types)):
            mask_cells = cell_types_idx == ct_idx
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (
                    expr_coo.data[entry_mask],
                    (expr_coo.row[entry_mask], expr_coo.col[entry_mask]),
                ),
                shape=expr_bool.shape,
            )

            all_rec_genes_ct = set()
            for _, _, rct, _, rec_str, _, _ in col_meta:
                if rct == ct_idx:
                    all_rec_genes_ct.update(_split_gene_complex(rec_str))

            all_rec_genes_ct = sorted([g for g in all_rec_genes_ct if g in gene_to_idx])

            rec_gene_matrices = {}
            for rec_gene in all_rec_genes_ct:
                rec_gene_idx = gene_to_idx[rec_gene]
                rec_mask = sub.col == rec_gene_idx
                rows_cells = sub.row[rec_mask]
                data_vals_ct = sub.data[rec_mask]
                patch_rows = np.array(
                    [patch_idx_map[self.patch_assignments[c]] for c in rows_cells],
                    dtype=int,
                )

                rec_coo = coo_matrix(
                    (data_vals_ct, (patch_rows, np.zeros(len(patch_rows), dtype=int))),
                    shape=(n_patches, 1),
                )
                rec_coo.sum_duplicates()
                rec_gene_matrices[rec_gene] = rec_coo.tocsr()

            patch_by_rec[ct_idx] = rec_gene_matrices

        # ─── 6) Build patch_by_cell matrix ────────────────────────────────────────────
        n_cells = adata.n_obs
        cell_patches = np.array(self.patch_assignments, dtype=int)
        pb_rows = [patch_idx_map[p] for p in cell_patches]
        pb_cols = list(range(n_cells))
        patch_by_cell = coo_matrix(
            (np.ones(n_cells, int), (pb_rows, pb_cols)), shape=(n_patches, n_cells)
        ).tocsr()

        # ─── 7) Compute LRI interactions ──────────────────────────────────────────────
        logger.debug("Computing LRI interactions...")
        row_inds, col_inds, data_vals = [], [], []

        for j, lig_ct_idx, rec_ct_idx, lig_str, rec_str, mode, sig_type in tqdm(
            col_meta, desc="Processing LRI interactions", unit="interaction"
        ):
            # Get ligand counts (支持配体复合体，AND逻辑)
            lig_genes = _split_gene_complex(lig_str)
            if lig_genes[0] in patch_by_lig[lig_ct_idx]:
                count_lig = np.array(
                    patch_by_lig[lig_ct_idx][lig_genes[0]].toarray()
                ).ravel()
                for lig_gene in lig_genes[1:]:  # 如果是复合体
                    if lig_gene in patch_by_lig[lig_ct_idx]:
                        count_lig = np.minimum(
                            count_lig,
                            np.array(
                                patch_by_lig[lig_ct_idx][lig_gene].toarray()
                            ).ravel(),
                        )
                    else:
                        count_lig = np.zeros(n_patches, dtype=int)
                        break
            else:
                count_lig = np.zeros(n_patches, dtype=int)

            # Get receptor counts (支持受体复合体，AND逻辑)
            rec_genes = _split_gene_complex(rec_str)
            if rec_genes[0] in patch_by_rec[rec_ct_idx]:
                count_rec = np.array(
                    patch_by_rec[rec_ct_idx][rec_genes[0]].toarray()
                ).ravel()
                for rec_gene in rec_genes[1:]:
                    if rec_gene in patch_by_rec[rec_ct_idx]:
                        count_rec = np.minimum(
                            count_rec,
                            np.array(
                                patch_by_rec[rec_ct_idx][rec_gene].toarray()
                            ).ravel(),
                        )
                    else:
                        count_rec = np.zeros(n_patches, dtype=int)
                        break
            else:
                count_rec = np.zeros(n_patches, dtype=int)

            # Compute autocrine for same cell type (支持配体和受体复合体)
            if lig_ct_idx == rec_ct_idx:
                coexpr = np.ones(n_cells, dtype=int)

                # All ligand genes must be co-expressed
                for lig_gene in lig_genes:
                    lig_gene_idx = gene_to_idx[lig_gene]
                    coexpr = coexpr * expr_bool[
                        :, lig_gene_idx
                    ].toarray().ravel().astype(int)

                # All receptor genes must be co-expressed
                for rec_gene in rec_genes:
                    rec_gene_idx = gene_to_idx[rec_gene]
                    coexpr = coexpr * expr_bool[
                        :, rec_gene_idx
                    ].toarray().ravel().astype(int)

                # Only cells of the correct type
                coexpr = coexpr * (cell_types_idx == lig_ct_idx).astype(int)
                auto = np.array(patch_by_cell.dot(coexpr)).ravel()
            else:
                auto = np.zeros(n_patches, int)

            # Fill matrix based on mode
            if mode == "juxtacrine":
                # Cell-Cell Contact: exclude autocrine (same-cell interactions)
                juxta = count_lig * count_rec - auto
                rows = np.nonzero(juxta)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(juxta[rows].tolist())

            elif mode == "autocrine":
                # Non-contact autocrine: same-cell co-expression only
                rows = np.nonzero(auto)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(auto[rows].tolist())

            else:  # paracrine
                # Non-contact paracrine: different-cell interactions
                if lig_ct_idx == rec_ct_idx:
                    para = count_lig * count_rec - auto
                else:
                    para = count_lig * count_rec

                rows = np.nonzero(para)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(para[rows].tolist())

        # ─── 9) Assemble final matrix ─────────────────────────────────────────────────
        patch_lri_matrix = csr_matrix(
            (data_vals, (row_inds, col_inds)), shape=(n_patches, n_columns), dtype=int
        )

        self.patch_lri_matrix = patch_lri_matrix
        logger.debug(
            f"Matrix density: {patch_lri_matrix.nnz / (n_patches * n_columns) * 100:.2f}%"
        )
        return patch_lri_matrix

    def run_patchify(
        self,
        adata: anndata.AnnData | dict[str, anndata.AnnData],
        output_dir: str | None = None,
        multi_sample: bool = False,
        sample_column: str | None = None,
    ) -> dict:
        """
        Run the complete patch-based LRI analysis.

        Supports three input modes:
        1. Single AnnData with multi_sample=False: standard single-sample analysis
        2. Dict of AnnData objects: automatically treated as multi-sample
        3. Single merged AnnData with multi_sample=True: split by sample_column

        For multiple samples, patch IDs are made globally unique and all samples
        share the same column structure based on the intersection of genes.

        Parameters
        ----------
        adata : anndata.AnnData or Dict[str, anndata.AnnData]
            Single AnnData object, or a dictionary mapping sample_id -> AnnData.
            For multi-sample mode, all samples will be processed with a shared
            gene set (intersection) and unified column structure.
        output_dir : str, optional
            Directory to save results. If None, results are not saved to disk.
        multi_sample : bool, default False
            Whether to treat a single AnnData as containing multiple samples.
            Automatically set to True if adata is a dict.
            If True and adata is a single AnnData, sample_column must be provided.
        sample_column : str, optional
            Column name in adata.obs that identifies different samples.
            Required when multi_sample=True and adata is a single AnnData.
            Ignored when adata is a dict.

        Returns
        -------
        results : dict
            - patch_lri_matrix: sparse matrix (n_patches x n_columns)
            - column_names: list of column names
            - sample_info: (multi-sample only) dict mapping sample_id -> {n_patches, n_cells, ...}

        Notes
        -----
        The function adds 'patch_idx' column to adata.obs, mapping each cell to its
        corresponding row in patch_lri_matrix. For multi-sample mode, patch_idx is
        the global index across all samples.

        Examples
        --------
        # Mode 1: Single sample
        >>> results = analyzer.run_patchify(adata)

        # Mode 2: Dict of samples
        >>> results = analyzer.run_patchify({'sample_A': adata_a, 'sample_B': adata_b})

        # Mode 3: Merged AnnData with sample column
        >>> results = analyzer.run_patchify(
        ...     adata_merged,
        ...     multi_sample=True,
        ...     sample_column='patient_id'
        ... )
        """
        # Route to appropriate method based on input type
        if isinstance(adata, dict):
            # Dict input: automatically multi-sample mode
            if len(adata) == 0:
                raise ValueError("adata_dict cannot be empty")
            elif len(adata) == 1:
                # Single sample in dict format - extract and run single mode
                sample_id, single_adata = next(iter(adata.items()))
                logger.debug(
                    f"Single sample detected ('{sample_id}'), running in single-sample mode"
                )
                return self._run_patchify_single(single_adata, output_dir)
            else:
                # Multiple samples
                return self._run_patchify_multi(adata, output_dir)
        else:
            # Single AnnData object
            if multi_sample:
                # Split by sample_column and run multi-sample mode
                if sample_column is None:
                    raise ValueError(
                        "sample_column must be provided when multi_sample=True "
                        "and adata is a single AnnData object"
                    )
                if sample_column not in adata.obs.columns:
                    raise ValueError(
                        f"sample_column '{sample_column}' not found in adata.obs. "
                        f"Available columns: {list(adata.obs.columns)}"
                    )

                # Split merged AnnData into dict
                adata_dict = self._split_adata_by_sample(adata, sample_column)
                return self._run_patchify_multi(adata_dict, output_dir)
            else:
                # Standard single-sample mode
                return self._run_patchify_single(adata, output_dir)

    def _run_patchify_single(
        self, adata: anndata.AnnData, output_dir: str | None = None
    ) -> dict:
        """
        Run patch-based LRI analysis for a single AnnData object.
        This is the original run_patchify logic.
        """
        logger.debug("Starting patch-based LRI analysis (single sample)...")
        logger.debug(f"Patch size: {self.patch_size} μm")
        logger.debug(f"Data shape: {adata.shape}")

        # Step 1: Create spatial patches
        patch_assignments, patch_info = self.create_spatial_patches(adata)

        # Step 2: Prepare LRI database
        lr_pairs, ligand_genes_list, receptor_genes_list, signaling_types = (
            self.prepare_lri_database(adata=adata)
        )

        # Step 3: Create column structure
        column_names = self.create_column_structure(adata, signaling_types)

        # Step 4: Build patch-LRI matrix
        patch_lri_matrix = self.build_patch_lri_matrix(adata, signaling_types)

        # Step 5: Remove zero columns
        patch_lri_matrix, column_names = self.remove_zero_columns(
            patch_lri_matrix, column_names
        )

        # Step 6: Add patch_idx to adata.obs
        # Map each cell to its patch matrix row index
        unique_patches = np.array(
            [p for p in np.unique(patch_assignments) if p in self.patch_coords]
        )
        patch_id_to_idx = {pid: i for i, pid in enumerate(unique_patches)}

        patch_idx_array = np.array(
            [patch_id_to_idx.get(patch_assignments[i], -1) for i in range(adata.n_obs)]
        )
        adata.obs["patch_idx"] = patch_idx_array
        logger.debug(
            f"Added 'patch_idx' to adata.obs ({(patch_idx_array >= 0).sum()} cells assigned to patches)"
        )

        # Step 7: Save results (optional)
        if output_dir is not None:
            logger.debug("Saving results...")
            os.makedirs(output_dir, exist_ok=True)

            # Save sparse matrix
            matrix_file = os.path.join(output_dir, "patch_lri_matrix.npz")
            sparse.save_npz(matrix_file, patch_lri_matrix)

            # Save column names
            columns_file = os.path.join(output_dir, "patch_lri_columns.csv")
            pd.DataFrame({"column_name": self.column_names}).to_csv(
                columns_file, index=False
            )

            # Save analysis parameters
            params_file = os.path.join(output_dir, "analysis_parameters.csv")
            params_df = pd.DataFrame(
                {
                    "parameter": [
                        "patch_size",
                        "resource_name",
                        "n_patches",
                        "n_lri_combinations",
                        "matrix_sparsity",
                    ],
                    "value": [
                        self.patch_size,
                        self.resource_name,
                        patch_info["n_patches"],
                        len(column_names),
                        f"{(1 - patch_lri_matrix.nnz / np.prod(patch_lri_matrix.shape)) * 100:.2f}%",
                    ],
                }
            )
            params_df.to_csv(params_file, index=False)

            # Save patch metadata
            patch_metadata_list = []
            for i, patch_id in enumerate(unique_patches):
                if patch_id in self.patch_coords:
                    x_center, y_center = self.patch_coords[patch_id]
                    patch_metadata_list.append(
                        {
                            "patch_idx": i,
                            "patch_id": patch_id,
                            "x_center": x_center,
                            "y_center": y_center,
                            "x_min": x_center - self.patch_size / 2,
                            "x_max": x_center + self.patch_size / 2,
                            "y_min": y_center - self.patch_size / 2,
                            "y_max": y_center + self.patch_size / 2,
                        }
                    )
            patch_metadata_df = pd.DataFrame(patch_metadata_list)
            patch_metadata_file = os.path.join(output_dir, "patch_metadata.parquet")
            patch_metadata_df.to_parquet(patch_metadata_file)

            logger.debug(f"Results saved to: {output_dir}")
            logger.debug(f"- Patch-LRI matrix: {matrix_file}")
            logger.debug(f"- Column names: {columns_file}")
            logger.debug(f"- Analysis parameters: {params_file}")
            logger.debug(f"- Patch metadata: {patch_metadata_file}")

        return {"patch_lri_matrix": patch_lri_matrix, "column_names": column_names}

    def _get_shared_genes(self, adata_dict: dict[str, anndata.AnnData]) -> list[str]:
        """
        Get the intersection of genes across all AnnData objects.

        Parameters
        ----------
        adata_dict : Dict[str, anndata.AnnData]
            Dictionary mapping sample_id -> AnnData

        Returns
        -------
        shared_genes : List[str]
            Sorted list of genes present in all samples
        """
        gene_sets = [set(adata.var_names) for adata in adata_dict.values()]
        shared_genes = set.intersection(*gene_sets)
        logger.debug(
            f"Gene intersection across {len(adata_dict)} samples: {len(shared_genes)} genes"
        )
        for sample_id, adata in adata_dict.items():
            logger.debug(f"  {sample_id}: {len(adata.var_names)} genes")
        return sorted(shared_genes)

    def _get_shared_cell_types(
        self, adata_dict: dict[str, anndata.AnnData]
    ) -> list[str]:
        """
        Get the union of cell types across all AnnData objects.

        Parameters
        ----------
        adata_dict : Dict[str, anndata.AnnData]
            Dictionary mapping sample_id -> AnnData

        Returns
        -------
        all_cell_types : List[str]
            Sorted list of all unique cell types across samples
        """
        all_cell_types = set()
        total_nan_cells = 0
        total_cells = 0

        for sample_id, adata in adata_dict.items():
            # Get cell types, handling both categorical and non-categorical columns
            ct_col = adata.obs[self.cell_type_column]
            if hasattr(ct_col, "cat"):
                sample_cell_types = ct_col.cat.categories.tolist()
            else:
                sample_cell_types = ct_col.dropna().unique().tolist()

            # Count nan values
            n_nan = ct_col.isna().sum()
            total_nan_cells += n_nan
            total_cells += len(ct_col)

            # Filter out any nan from cell types list (shouldn't be in categories, but just in case)
            sample_cell_types = [ct for ct in sample_cell_types if pd.notna(ct)]

            all_cell_types.update(sample_cell_types)
            logger.debug(
                f"  {sample_id}: {len(sample_cell_types)} cell types"
                + (f" ({n_nan} cells with nan)" if n_nan > 0 else "")
            )

        if total_nan_cells > 0:
            logger.debug(
                f"  Warning: {total_nan_cells}/{total_cells} total cells have missing cell types"
            )

        return sorted(all_cell_types)

    def _run_patchify_multi(
        self, adata_dict: dict[str, anndata.AnnData], output_dir: str | None = None
    ) -> dict:
        """
        Run patch-based LRI analysis for multiple AnnData objects.

        Parameters
        ----------
        adata_dict : Dict[str, anndata.AnnData]
            Dictionary mapping sample_id -> AnnData
        output_dir : str, optional
            Directory to save results

        Returns
        -------
        results : dict
            Dictionary containing combined results from all samples
        """
        logger.debug("=" * 60)
        logger.debug("Starting patch-based LRI analysis (multi-sample mode)")
        logger.debug("=" * 60)
        logger.debug(f"Number of samples: {len(adata_dict)}")
        logger.debug(f"Patch size: {self.patch_size} μm")
        for sample_id, adata in adata_dict.items():
            logger.debug(
                f"  {sample_id}: {adata.shape[0]} cells, {adata.shape[1]} genes"
            )

        # Step 1: Get shared genes (intersection)
        logger.debug("[Step 1/6] Computing gene intersection...")
        shared_genes = self._get_shared_genes(adata_dict)

        # Step 2: Get all cell types (union) and set as shared
        logger.debug("[Step 2/6] Computing cell type union...")
        shared_cell_types = self._get_shared_cell_types(adata_dict)
        self.cell_types = shared_cell_types
        logger.debug(f"Total unique cell types: {len(shared_cell_types)}")

        # Step 3: Prepare LRI database using shared genes
        logger.debug("[Step 3/6] Preparing LRI database with shared genes...")
        lr_pairs, ligand_genes_list, receptor_genes_list, signaling_types = (
            self.prepare_lri_database(gene_names=shared_genes)
        )

        # Step 4: Create unified column structure
        # We need a reference adata to get cell types, but we'll use shared_cell_types
        logger.debug("[Step 4/6] Creating unified column structure...")
        column_names = self._create_column_structure_from_cell_types(
            shared_cell_types, signaling_types
        )
        self.column_names = column_names
        logger.debug(f"Total LRI columns: {len(column_names)}")

        # Step 5: Process each sample
        logger.debug("[Step 5/6] Processing each sample...")
        all_matrices = []
        all_patch_metadata = []
        sample_info = {}
        global_patch_idx = 0

        for sample_id, adata in adata_dict.items():
            logger.debug(f"--- Processing sample: {sample_id} ---")

            # Create spatial patches for this sample
            patch_assignments, patch_info = self.create_spatial_patches(adata)

            # Build patch-LRI matrix using the pre-defined column structure
            sample_matrix = self.build_patch_lri_matrix(adata, signaling_types)

            # Get unique patches in order (matching matrix rows)
            unique_patches = np.array(
                [p for p in np.unique(patch_assignments) if p in self.patch_coords]
            )
            n_patches = len(unique_patches)

            # Count cells per patch
            patch_cell_counts = {}
            for local_pid in patch_assignments:
                if local_pid in self.patch_coords:
                    patch_cell_counts[local_pid] = (
                        patch_cell_counts.get(local_pid, 0) + 1
                    )

            # Create patch metadata for this sample
            for i, local_patch_id in enumerate(unique_patches):
                center_x, center_y = self.patch_coords[local_patch_id]
                all_patch_metadata.append(
                    {
                        "sample_id": sample_id,
                        "local_patch_id": local_patch_id,
                        "global_patch_idx": global_patch_idx + i,
                        "center_x": center_x,
                        "center_y": center_y,
                        "n_cells": patch_cell_counts.get(local_patch_id, 0),
                    }
                )

            # Add patch_idx to adata.obs (global index)
            patch_id_to_global_idx = {
                pid: global_patch_idx + i for i, pid in enumerate(unique_patches)
            }
            patch_idx_array = np.array(
                [
                    patch_id_to_global_idx.get(patch_assignments[i], -1)
                    for i in range(adata.n_obs)
                ]
            )
            adata.obs["patch_idx"] = patch_idx_array
            logger.debug(
                f"  Added 'patch_idx' to adata.obs ({(patch_idx_array >= 0).sum()} cells assigned)"
            )

            # Store sample info
            sample_info[sample_id] = {
                "n_cells": adata.n_obs,
                "n_patches": n_patches,
                "global_patch_idx_start": global_patch_idx,
                "global_patch_idx_end": global_patch_idx + n_patches - 1,
            }

            all_matrices.append(sample_matrix)
            global_patch_idx += n_patches

        # Step 6: Combine results
        logger.debug("[Step 6/6] Combining results...")

        # Vertically stack all matrices
        combined_matrix = sparse_vstack(all_matrices, format="csr")
        logger.debug(f"Combined matrix shape: {combined_matrix.shape}")

        # Create metadata DataFrame
        patch_metadata_df = pd.DataFrame(all_patch_metadata)

        # Remove zero columns from combined matrix
        combined_matrix, column_names = self.remove_zero_columns(
            combined_matrix, column_names
        )

        # Save results (optional)
        if output_dir is not None:
            logger.debug("Saving results...")
            os.makedirs(output_dir, exist_ok=True)

            # Save sparse matrix
            matrix_file = os.path.join(output_dir, "patch_lri_matrix.npz")
            sparse.save_npz(matrix_file, combined_matrix)

            # Save column names
            columns_file = os.path.join(output_dir, "patch_lri_columns.csv")
            pd.DataFrame({"column_name": column_names}).to_csv(
                columns_file, index=False
            )

            # Save patch metadata
            patch_metadata_file = os.path.join(output_dir, "patch_metadata.csv")
            patch_metadata_df.to_csv(patch_metadata_file, index=False)

            # Save sample info
            sample_info_file = os.path.join(output_dir, "sample_info.csv")
            sample_info_df = pd.DataFrame(
                [{"sample_id": sid, **info} for sid, info in sample_info.items()]
            )
            sample_info_df.to_csv(sample_info_file, index=False)

            # Save analysis parameters
            params_file = os.path.join(output_dir, "analysis_parameters.csv")
            params_df = pd.DataFrame(
                {
                    "parameter": [
                        "patch_size",
                        "resource_name",
                        "n_samples",
                        "total_patches",
                        "n_lri_combinations",
                        "n_shared_genes",
                        "matrix_sparsity",
                    ],
                    "value": [
                        self.patch_size,
                        self.resource_name,
                        len(adata_dict),
                        combined_matrix.shape[0],
                        len(column_names),
                        len(shared_genes),
                        f"{(1 - combined_matrix.nnz / np.prod(combined_matrix.shape)) * 100:.2f}%",
                    ],
                }
            )
            params_df.to_csv(params_file, index=False)

            logger.debug(f"Results saved to: {output_dir}")
            logger.debug(f"- Patch-LRI matrix: {matrix_file}")
            logger.debug(f"- Column names: {columns_file}")
            logger.debug(f"- Patch metadata: {patch_metadata_file}")
            logger.debug(f"- Sample info: {sample_info_file}")
            logger.debug(f"- Analysis parameters: {params_file}")

        logger.debug("=" * 60)
        logger.debug("Multi-sample analysis complete!")
        logger.debug(f"Total patches: {combined_matrix.shape[0]}")
        logger.debug(f"Total LRI columns: {combined_matrix.shape[1]}")
        logger.debug("=" * 60)

        return {
            "patch_lri_matrix": combined_matrix,
            "column_names": column_names,
            "sample_info": sample_info,
        }

    def _create_column_structure_from_cell_types(
        self, cell_types: list[str], signaling_types: list[str]
    ) -> list[str]:
        """
        Create column names using a predefined list of cell types.
        This is used in multi-sample mode where we need consistent columns across samples.

        Parameters
        ----------
        cell_types : List[str]
            List of all cell types to include
        signaling_types : List[str]
            Signaling type for each LR pair

        Returns
        -------
        column_names : List[str]
            List of column names
        """
        column_names = []

        for idx, (lig, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]

            for lig_ct in cell_types:
                for rec_ct in cell_types:
                    if sig_type == "Cell-Cell Contact":
                        column_names.append(
                            f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}juxtacrine"
                        )
                    else:
                        if lig_ct == rec_ct:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}autocrine"
                            )
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )
                        else:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )

        return column_names


# ============================================================================
# Cell Neighborhood-based LRI Analysis
# (From single_cell.py)
# ============================================================================


class NeighborhoodLRIAnalyzer(BaseLRIAnalyzer):
    """
    Cell Neighborhood-based Ligand-Receptor Interaction Analyzer

    For each cell, defines a square neighborhood centered on that cell and counts
    ligand-receptor interactions within the neighborhood using all-to-all strategy.
    """

    def __init__(
        self,
        neighborhood_size: float = 50.0,
        resource_name: str = "cellchatdb",
        spliter: str = "|",
        cellchatdb_path: str | None = None,
        cellphonedb_path: str | None = None,
        cell_type_column: str = "cell_type",
        include_gene_expression: bool = False,
        raw_count_location: str = "X",
    ):
        """
        Initialize the neighborhood-based LRI analyzer.

        Parameters
        ----------
        neighborhood_size : float, default 50.0
            Size of square neighborhood in micrometers (edge length)
        resource_name : str, default 'cellchatdb'
            LRI database to use from liana
        spliter : str, default '|'
            Separator for column names
        cellchatdb_path : str, optional
            Path to local CellChatDB CSV file. If None, uses bundled database.
        cellphonedb_path : str, optional
            Path to local CellPhoneDB CSV file. If None, uses bundled database.
        cell_type_column : str
            Column name for cell types in adata.obs
        include_gene_expression : bool
            Whether to append gene expression to the matrix
        raw_count_location : str
            Where to get raw counts from ('X', 'raw', or 'layers:NAME')
        """
        # Initialize base class
        super().__init__(
            resource_name=resource_name,
            spliter=spliter,
            cellchatdb_path=cellchatdb_path,
            cellphonedb_path=cellphonedb_path,
            cell_type_column=cell_type_column,
        )

        # Neighborhood-specific attributes
        self.neighborhood_size = neighborhood_size
        self.neighborhoods = None
        self.cell_lri_matrix = None
        self.include_gene_expression = include_gene_expression
        self.raw_count_location = raw_count_location

    def get_raw_counts(self, adata: anndata.AnnData) -> np.ndarray:
        """
        Get raw count matrix from specified location.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data

        Returns
        -------
        raw_counts : np.ndarray or sparse matrix
            Raw count matrix
        """
        if self.raw_count_location == "X":
            return adata.X
        elif self.raw_count_location == "raw":
            if adata.raw is None:
                raise ValueError("adata.raw is None but raw_count_location='raw'")
            return adata.raw.X
        elif self.raw_count_location.startswith("layers:"):
            layer_name = self.raw_count_location.split(":", 1)[1]
            if layer_name not in adata.layers:
                raise ValueError(f"Layer '{layer_name}' not found in adata.layers")
            return adata.layers[layer_name]
        else:
            raise ValueError(f"Invalid raw_count_location: {self.raw_count_location}")

    def build_neighborhoods(self, adata: anndata.AnnData) -> dict[int, np.ndarray]:
        """
        For each cell, find all cells within its square neighborhood.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data

        Returns
        -------
        neighborhoods : Dict[int, np.ndarray]
            Dictionary mapping cell_index -> array of neighbor cell indices
        """
        logger.debug(f"Building neighborhoods (size={self.neighborhood_size}µm)...")
        coords = adata.obsm["spatial"][:, :2]
        n_cells = len(coords)

        # Build KD-Tree for efficient spatial queries
        tree = KDTree(coords)

        # Define half-width of square neighborhood
        half_size = self.neighborhood_size / 2.0

        neighborhoods = {}
        for i in tqdm(range(n_cells), desc="Building neighborhoods", unit="cell"):
            center = coords[i]

            # Define square boundaries
            x_min, x_max = center[0] - half_size, center[0] + half_size
            y_min, y_max = center[1] - half_size, center[1] + half_size

            # Query all points within the bounding box using KD-Tree
            # First get points within a circle, then filter to square
            radius = half_size * np.sqrt(2)  # Diagonal of square
            candidate_indices = tree.query_radius([center], r=radius)[0]

            # Filter to cells strictly within square
            candidate_coords = coords[candidate_indices]
            in_square = (
                (candidate_coords[:, 0] >= x_min)
                & (candidate_coords[:, 0] <= x_max)
                & (candidate_coords[:, 1] >= y_min)
                & (candidate_coords[:, 1] <= y_max)
            )

            neighborhoods[i] = candidate_indices[in_square]

        self.neighborhoods = neighborhoods
        logger.debug(
            f"Average neighborhood size: {np.mean([len(n) for n in neighborhoods.values()]):.1f} cells"
        )
        return neighborhoods

    def build_cell_lri_matrix(
        self, adata: anndata.AnnData, signaling_types: list[str]
    ) -> csr_matrix:
        """
        Build the cell-LRI interaction matrix with full vectorization.
        Supports both ligand and receptor complexes (using AND logic for co-expression).
        """
        logger.debug(
            "Building cell-LRI matrix with vectorized neighborhood interactions..."
        )

        n_cells = adata.n_obs

        # ─── 1) Index mappings ────────────────────────────────────────────────────
        ct_to_idx = {ct: i for i, ct in enumerate(self.cell_types)}

        # Check for nan/missing cell types and map them
        cell_type_values = adata.obs[self.cell_type_column].values
        cell_types_idx = []
        for ct in cell_type_values:
            if pd.isna(ct):
                # Assign -1 for nan cell types (will be excluded from analysis)
                cell_types_idx.append(-1)
            elif ct not in ct_to_idx:
                # Cell type not in shared cell types
                cell_types_idx.append(-1)
            else:
                cell_types_idx.append(ct_to_idx[ct])
        cell_types_idx = np.array(cell_types_idx, dtype=int)

        # Report cells with missing/invalid cell types
        n_invalid = np.sum(cell_types_idx == -1)
        if n_invalid > 0:
            logger.debug(
                f"  Warning: {n_invalid} cells ({n_invalid / len(cell_types_idx) * 100:.1f}%) have missing/invalid cell types and will be excluded"
            )

        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # ─── 2) Parse column metadata - NOW HANDLES LIGAND COMPLEXES ─────────────
        col_meta = []
        for idx, (lig_str, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]

            # Check if all ligand genes exist (支持配体复合体)
            lig_genes = _split_gene_complex(lig_str)
            if any(lg not in gene_to_idx for lg in lig_genes):
                continue

            # Check if all receptor genes exist
            rec_genes = _split_gene_complex(rec_str)
            if any(rg not in gene_to_idx for rg in rec_genes):
                continue

            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == "Cell-Cell Contact":
                        # 一列：juxtacrine
                        col_meta.append(
                            (
                                len(col_meta),  # column j
                                ct_to_idx[lig_ct],  # lig_ct_idx
                                ct_to_idx[rec_ct],  # rec_ct_idx
                                lig_str,  # ligand string (可能是复合体)
                                rec_str,  # receptor string (可能是复合体)
                                "juxtacrine",  # mode
                                sig_type,  # signaling type
                            )
                        )
                    else:
                        # 非接触
                        if lig_ct == rec_ct:
                            # 两列：autocrine + paracrine
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "autocrine",
                                    sig_type,
                                )
                            )
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "paracrine",
                                    sig_type,
                                )
                            )
                        else:
                            # 一列：paracrine
                            col_meta.append(
                                (
                                    len(col_meta),
                                    ct_to_idx[lig_ct],
                                    ct_to_idx[rec_ct],
                                    lig_str,
                                    rec_str,
                                    "paracrine",
                                    sig_type,
                                )
                            )

        n_columns = len(col_meta)
        logger.debug(f"Processing {n_cells} cells × {n_columns} LRI combinations")

        # ─── 3) Binarize expression ───────────────────────────────────────────────
        X = adata.X
        if sp.issparse(X):
            expr_bool = (X > 0).astype(int).tocsc()
        else:
            expr_bool = csr_matrix((X > 0).astype(int)).tocsc()

        expr_coo = expr_bool.tocoo()

        # ─── 4) Build cell-neighborhood adjacency matrix ──────────────────────────
        logger.debug("Building cell-neighborhood adjacency matrix...")
        rows, cols = [], []
        for cell_idx, neighbors in self.neighborhoods.items():
            rows.extend([cell_idx] * len(neighbors))
            cols.extend(neighbors)
        cell_nbr_matrix = coo_matrix(
            (np.ones(len(rows), dtype=int), (rows, cols)), shape=(n_cells, n_cells)
        ).tocsr()

        # ─── 5) Build neighborhood_by_lig for INDIVIDUAL ligand genes ─────────────
        neighborhood_by_lig = {}
        logger.debug("Building neighborhood-by-ligand matrices (individual genes)...")

        for ct_idx in range(len(self.cell_types)):
            mask_cells = cell_types_idx == ct_idx
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (
                    expr_coo.data[entry_mask],
                    (expr_coo.row[entry_mask], expr_coo.col[entry_mask]),
                ),
                shape=expr_bool.shape,
            )

            # Collect all individual ligand genes for this cell type (展开复合体)
            all_lig_genes_ct = set()
            for _, lct, _, lig_str, _, _, _ in col_meta:
                if lct == ct_idx:
                    all_lig_genes_ct.update(_split_gene_complex(lig_str))
            all_lig_genes_ct = sorted([g for g in all_lig_genes_ct if g in gene_to_idx])

            # Build matrix for each individual ligand gene
            lig_gene_matrices = {}
            for lig_gene in all_lig_genes_ct:
                lig_gene_idx = gene_to_idx[lig_gene]
                lig_mask = sub.col == lig_gene_idx
                rows_cells = sub.row[lig_mask]
                data_vals_ct = sub.data[lig_mask]

                cell_lig = coo_matrix(
                    (data_vals_ct, (rows_cells, np.zeros(len(rows_cells), dtype=int))),
                    shape=(n_cells, 1),
                )
                cell_lig.sum_duplicates()
                # 邻域聚合
                lig_gene_matrices[lig_gene] = cell_nbr_matrix.dot(cell_lig.tocsr())

            neighborhood_by_lig[ct_idx] = lig_gene_matrices

        # ─── 6) Build neighborhood_by_rec for INDIVIDUAL receptor genes ───────────
        neighborhood_by_rec = {}
        logger.debug("Building neighborhood-by-receptor matrices (individual genes)...")

        for ct_idx in range(len(self.cell_types)):
            mask_cells = cell_types_idx == ct_idx
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (
                    expr_coo.data[entry_mask],
                    (expr_coo.row[entry_mask], expr_coo.col[entry_mask]),
                ),
                shape=expr_bool.shape,
            )

            # Collect all individual receptor genes for this cell type
            all_rec_genes_ct = set()
            for _, _, rct, _, rec_str, _, _ in col_meta:
                if rct == ct_idx:
                    all_rec_genes_ct.update(_split_gene_complex(rec_str))
            all_rec_genes_ct = sorted([g for g in all_rec_genes_ct if g in gene_to_idx])

            # Build matrix for each individual receptor gene
            rec_gene_matrices = {}
            for rec_gene in all_rec_genes_ct:
                rec_gene_idx = gene_to_idx[rec_gene]
                rec_mask = sub.col == rec_gene_idx
                rows_cells = sub.row[rec_mask]
                data_vals_ct = sub.data[rec_mask]

                cell_rec = coo_matrix(
                    (data_vals_ct, (rows_cells, np.zeros(len(rows_cells), dtype=int))),
                    shape=(n_cells, 1),
                )
                cell_rec.sum_duplicates()
                # 邻域聚合
                rec_gene_matrices[rec_gene] = cell_nbr_matrix.dot(cell_rec.tocsr())

            neighborhood_by_rec[ct_idx] = rec_gene_matrices

        # ─── 7) Compute interactions with BOTH ligand and receptor complexes ──────
        logger.debug("Computing LRI interactions...")
        row_inds, col_inds, data_vals = [], [], []

        for j, lig_ct_idx, rec_ct_idx, lig_str, rec_str, mode, sig_type in tqdm(
            col_meta, desc="Processing LRI interactions", unit="interaction"
        ):
            # Ligand count (支持配体复合体，AND逻辑)
            lig_genes = _split_gene_complex(lig_str)
            if lig_genes[0] in neighborhood_by_lig[lig_ct_idx]:
                count_lig = np.array(
                    neighborhood_by_lig[lig_ct_idx][lig_genes[0]].toarray()
                ).ravel()
                for lig_gene in lig_genes[1:]:  # 如果是复合体
                    if lig_gene in neighborhood_by_lig[lig_ct_idx]:
                        count_lig = np.minimum(
                            count_lig,
                            np.array(
                                neighborhood_by_lig[lig_ct_idx][lig_gene].toarray()
                            ).ravel(),
                        )
                    else:
                        count_lig = np.zeros(n_cells, dtype=int)
                        break
            else:
                count_lig = np.zeros(n_cells, dtype=int)

            # Receptor count (支持受体复合体，AND逻辑)
            rec_genes = _split_gene_complex(rec_str)
            if rec_genes[0] in neighborhood_by_rec[rec_ct_idx]:
                count_rec = np.array(
                    neighborhood_by_rec[rec_ct_idx][rec_genes[0]].toarray()
                ).ravel()
                for rec_gene in rec_genes[1:]:
                    if rec_gene in neighborhood_by_rec[rec_ct_idx]:
                        count_rec = np.minimum(
                            count_rec,
                            np.array(
                                neighborhood_by_rec[rec_ct_idx][rec_gene].toarray()
                            ).ravel(),
                        )
                    else:
                        count_rec = np.zeros(n_cells, dtype=int)
                        break
            else:
                count_rec = np.zeros(n_cells, dtype=int)

            # Autocrine co-expression (支持配体和受体复合体)
            if lig_ct_idx == rec_ct_idx:
                coexpr = np.ones(n_cells, dtype=int)

                # All ligand genes must be co-expressed
                for lig_gene in lig_genes:
                    lig_gene_idx = gene_to_idx[lig_gene]
                    coexpr = coexpr * expr_bool[
                        :, lig_gene_idx
                    ].toarray().ravel().astype(int)

                # All receptor genes must be co-expressed
                for rec_gene in rec_genes:
                    rec_gene_idx = gene_to_idx[rec_gene]
                    coexpr = coexpr * expr_bool[
                        :, rec_gene_idx
                    ].toarray().ravel().astype(int)

                # Only cells of the correct type
                coexpr = coexpr * (cell_types_idx == lig_ct_idx).astype(int)
                auto = np.array(cell_nbr_matrix.dot(coexpr)).ravel()
            else:
                auto = np.zeros(n_cells, dtype=int)

            # Write interactions based on mode
            if mode == "juxtacrine":
                # Cell-Cell Contact: exclude autocrine (same-cell interactions)
                juxta = count_lig * count_rec - auto
                rows_nz = np.nonzero(juxta)[0]
                row_inds.extend(rows_nz.tolist())
                col_inds.extend([j] * len(rows_nz))
                data_vals.extend(juxta[rows_nz].tolist())

            elif mode == "autocrine":
                rows_nz = np.nonzero(auto)[0]
                row_inds.extend(rows_nz.tolist())
                col_inds.extend([j] * len(rows_nz))
                data_vals.extend(auto[rows_nz].tolist())

            else:  # paracrine
                if lig_ct_idx == rec_ct_idx:
                    para = count_lig * count_rec - auto  # 去掉同细胞的自分泌部分
                else:
                    para = count_lig * count_rec
                rows_nz = np.nonzero(para)[0]
                row_inds.extend(rows_nz.tolist())
                col_inds.extend([j] * len(rows_nz))
                data_vals.extend(para[rows_nz].tolist())

        # ─── 8) Assemble final sparse matrix ──────────────────────────────────────
        cell_lri_matrix = csr_matrix(
            (data_vals, (row_inds, col_inds)), shape=(n_cells, n_columns), dtype=int
        )
        self.cell_lri_matrix = cell_lri_matrix
        logger.debug(
            f"Matrix density: {cell_lri_matrix.nnz / (n_cells * n_columns) * 100:.2f}%"
        )
        return cell_lri_matrix

    def create_metadata_dataframe(self, adata: anndata.AnnData) -> pd.DataFrame:
        """
        Create metadata dataframe with cell information.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data

        Returns
        -------
        cell_metadata_df : pd.DataFrame
            DataFrame with cell metadata
        """
        logger.debug("Creating metadata dataframe...")

        coords = adata.obsm["spatial"][:, :2]
        neighborhood_sizes = [len(self.neighborhoods[i]) for i in range(adata.n_obs)]

        cell_metadata_df = pd.DataFrame(
            {
                "cell_id": adata.obs.index,
                # 'tma_id': adata.obs['tma_id'],
                "cell_type": adata.obs["cell_type"],
                "x_coord": coords[:, 0],
                "y_coord": coords[:, 1],
                "neighborhood_size": neighborhood_sizes,
            }
        )

        return cell_metadata_df

    def run_neighborhood(
        self,
        adata: anndata.AnnData | dict[str, anndata.AnnData],
        output_dir: str | None = None,
        required_columns: list[str] | None = None,
        multi_sample: bool = False,
        sample_column: str | None = None,
    ) -> dict:
        """
        Run the complete neighborhood-based LRI analysis.

        Supports three input modes:
        1. Single AnnData with multi_sample=False: standard single-sample analysis
        2. Dict of AnnData objects: automatically treated as multi-sample
        3. Single merged AnnData with multi_sample=True: split by sample_column

        For multiple samples, cell indices are made globally unique and all samples
        share the same column structure based on the intersection of genes.

        Parameters
        ----------
        adata : anndata.AnnData or Dict[str, anndata.AnnData]
            Single AnnData object, or a dictionary mapping sample_id -> AnnData.
            For multi-sample mode, all samples will be processed with a shared
            gene set (intersection) and unified column structure.
        output_dir : str, optional
            Directory to save results. If None, results are not saved to disk.
        required_columns : list of str, optional
            If provided, subset the cell-LRI matrix to only these columns (e.g., from patch analysis).
            This ensures column alignment with a reference matrix (like patch-LRI matrix).
            If None, removes zero columns as usual.
        multi_sample : bool, default False
            Whether to treat a single AnnData as containing multiple samples.
            Automatically set to True if adata is a dict.
            If True and adata is a single AnnData, sample_column must be provided.
        sample_column : str, optional
            Column name in adata.obs that identifies different samples.
            Required when multi_sample=True and adata is a single AnnData.
            Ignored when adata is a dict.

        Returns
        -------
        results : dict
            For single sample:
                - cell_lri_matrix: sparse matrix (n_cells x n_columns)
                - column_names: list of column names
                - neighborhoods: dict of cell neighborhoods
            For multiple samples:
                - cell_lri_matrix: sparse matrix (total_cells x n_columns)
                - column_names: list of column names
                - sample_info: dict mapping sample_id -> {n_cells, global_cell_idx_start, global_cell_idx_end, avg_neighborhood_size}

        Examples
        --------
        # Mode 1: Single sample
        >>> results = analyzer.run_neighborhood(adata)

        # Mode 2: Dict of samples
        >>> results = analyzer.run_neighborhood({'sample_A': adata_a, 'sample_B': adata_b})

        # Mode 3: Merged AnnData with sample column
        >>> results = analyzer.run_neighborhood(
        ...     adata_merged,
        ...     multi_sample=True,
        ...     sample_column='patient_id'
        ... )
        """
        # Route to appropriate method based on input type
        if isinstance(adata, dict):
            # Dict input: automatically multi-sample mode
            if len(adata) == 0:
                raise ValueError("adata_dict cannot be empty")
            elif len(adata) == 1:
                # Single sample in dict format - extract and run single mode
                sample_id, single_adata = next(iter(adata.items()))
                logger.debug(
                    f"Single sample detected ('{sample_id}'), running in single-sample mode"
                )
                return self._run_neighborhood_single(
                    single_adata, output_dir, required_columns
                )
            else:
                # Multiple samples
                return self._run_neighborhood_multi(adata, output_dir, required_columns)
        else:
            # Single AnnData object
            if multi_sample:
                # Split by sample_column and run multi-sample mode
                if sample_column is None:
                    raise ValueError(
                        "sample_column must be provided when multi_sample=True "
                        "and adata is a single AnnData object"
                    )
                if sample_column not in adata.obs.columns:
                    raise ValueError(
                        f"sample_column '{sample_column}' not found in adata.obs. "
                        f"Available columns: {list(adata.obs.columns)}"
                    )

                # Split merged AnnData into dict
                adata_dict = self._split_adata_by_sample(adata, sample_column)
                return self._run_neighborhood_multi(
                    adata_dict, output_dir, required_columns
                )
            else:
                # Standard single-sample mode
                return self._run_neighborhood_single(
                    adata, output_dir, required_columns
                )

    def _run_neighborhood_single(
        self,
        adata: anndata.AnnData,
        output_dir: str | None = None,
        required_columns: list[str] | None = None,
    ) -> dict:
        """
        Run neighborhood-based LRI analysis for a single AnnData object.
        This is the original run_neighborhood logic.
        """
        logger.debug("Starting cell neighborhood-based LRI analysis (single sample)...")
        logger.debug(f"Neighborhood size: {self.neighborhood_size} µm")
        logger.debug(f"Data shape: {adata.shape}")

        # Step 1: Build neighborhoods
        neighborhoods = self.build_neighborhoods(adata)

        # Step 2: Prepare LRI database
        lr_pairs, ligand_genes_list, receptor_genes_list, signaling_types = (
            self.prepare_lri_database(adata=adata)
        )

        # Step 3: Create column structure
        column_names = self.create_column_structure(adata, signaling_types)

        # Step 4: Build cell-LRI matrix
        cell_lri_matrix = self.build_cell_lri_matrix(adata, signaling_types)

        # Step 5: Optionally append gene expression
        if self.include_gene_expression:
            logger.debug("Appending gene expression counts to matrix...")

            # Get raw counts
            raw_counts = self.get_raw_counts(adata)
            logger.debug("Using raw counts")

            # Convert to sparse CSR if needed
            if not sp.issparse(raw_counts):
                raw_counts = sp.csr_matrix(raw_counts)
            else:
                raw_counts = raw_counts.tocsr()

            # Concatenate horizontally: [LRI features | gene expression]
            combined_matrix = sp.hstack([cell_lri_matrix, raw_counts], format="csr")

            # Update column names
            gene_column_names = [
                f"GENE{self.spliter}{gene}" for gene in adata.var_names
            ]
            all_column_names = column_names + gene_column_names

            logger.debug(f"Combined matrix shape: {combined_matrix.shape}")
            logger.debug(f"  - LRI features: {len(column_names)}")
            logger.debug(f"  - Gene features: {len(gene_column_names)}")
        else:
            combined_matrix = cell_lri_matrix
            all_column_names = column_names

        # Step 6: Column filtering
        if required_columns is not None:
            # Subset to required columns (for alignment with patch matrix)
            logger.debug(
                f"Subsetting to {len(required_columns)} required columns (from reference matrix)..."
            )

            # Build column name to index mapping
            col_to_idx = {name: i for i, name in enumerate(all_column_names)}

            # Find common columns (preserve required_columns order)
            common_cols = [name for name in required_columns if name in col_to_idx]
            col_indices = np.array(
                [col_to_idx[name] for name in common_cols], dtype=int
            )

            logger.debug(f"  Cell matrix columns: {len(all_column_names)}")
            logger.debug(f"  Required columns: {len(required_columns)}")
            logger.debug(f"  Common columns (kept): {len(common_cols)}")

            if len(common_cols) < len(required_columns):
                missing = len(required_columns) - len(common_cols)
                logger.debug(
                    f"  Warning: {missing} required columns not found in cell matrix (will be absent)"
                )

            # Subset matrix
            combined_matrix = combined_matrix[:, col_indices]
            all_column_names = common_cols
        else:
            # Remove zero columns (default behavior)
            combined_matrix, all_column_names = self.remove_zero_columns(
                combined_matrix, all_column_names
            )

        # Step 7: Save results (optional)
        if output_dir is not None:
            logger.debug("Saving results...")
            os.makedirs(output_dir, exist_ok=True)

            # Save sparse matrix
            matrix_file = os.path.join(output_dir, "cell_lri_matrix.npz")
            sparse.save_npz(matrix_file, combined_matrix)

            # Save column names
            columns_file = os.path.join(output_dir, "cell_lri_columns.csv")
            pd.DataFrame({"column_name": all_column_names}).to_csv(
                columns_file, index=False
            )

            # Save analysis parameters
            params_file = os.path.join(output_dir, "analysis_parameters.csv")
            params_df = pd.DataFrame(
                {
                    "parameter": [
                        "neighborhood_size",
                        "resource_name",
                        "n_cells",
                        "n_lri_combinations",
                        "n_genes",
                        "n_total_features",
                        "matrix_sparsity",
                        "avg_neighborhood_size",
                        "include_gene_expression",
                        "raw_count_location",
                    ],
                    "value": [
                        self.neighborhood_size,
                        self.resource_name,
                        adata.n_obs,
                        len(column_names),
                        len(adata.var_names) if self.include_gene_expression else 0,
                        len(all_column_names),
                        f"{(1 - combined_matrix.nnz / np.prod(combined_matrix.shape)) * 100:.2f}%",
                        f"{np.mean([len(n) for n in neighborhoods.values()]):.1f}",
                        self.include_gene_expression,
                        self.raw_count_location
                        if self.include_gene_expression
                        else "N/A",
                    ],
                }
            )
            params_df.to_csv(params_file, index=False)

            logger.debug(f"Results saved to: {output_dir}")
            logger.debug(f"- Cell-LRI matrix: {matrix_file}")
            logger.debug(f"- Column names: {columns_file}")
            logger.debug(f"- Analysis parameters: {params_file}")

        return {
            "cell_lri_matrix": combined_matrix,
            "column_names": all_column_names,
            "neighborhoods": neighborhoods,
        }

    def _get_shared_genes(self, adata_dict: dict[str, anndata.AnnData]) -> list[str]:
        """
        Get the intersection of genes across all AnnData objects.
        """
        gene_sets = [set(adata.var_names) for adata in adata_dict.values()]
        shared_genes = set.intersection(*gene_sets)
        logger.debug(
            f"Gene intersection across {len(adata_dict)} samples: {len(shared_genes)} genes"
        )
        for sample_id, adata in adata_dict.items():
            logger.debug(f"  {sample_id}: {len(adata.var_names)} genes")
        return sorted(shared_genes)

    def _get_shared_cell_types(
        self, adata_dict: dict[str, anndata.AnnData]
    ) -> list[str]:
        """
        Get the union of cell types across all AnnData objects.
        """
        all_cell_types = set()
        total_nan_cells = 0
        total_cells = 0

        for sample_id, adata in adata_dict.items():
            ct_col = adata.obs[self.cell_type_column]
            if hasattr(ct_col, "cat"):
                sample_cell_types = ct_col.cat.categories.tolist()
            else:
                sample_cell_types = ct_col.dropna().unique().tolist()

            n_nan = ct_col.isna().sum()
            total_nan_cells += n_nan
            total_cells += len(ct_col)

            sample_cell_types = [ct for ct in sample_cell_types if pd.notna(ct)]
            all_cell_types.update(sample_cell_types)
            logger.debug(
                f"  {sample_id}: {len(sample_cell_types)} cell types"
                + (f" ({n_nan} cells with nan)" if n_nan > 0 else "")
            )

        if total_nan_cells > 0:
            logger.debug(
                f"  Warning: {total_nan_cells}/{total_cells} total cells have missing cell types"
            )

        return sorted(all_cell_types)

    def _run_neighborhood_multi(
        self,
        adata_dict: dict[str, anndata.AnnData],
        output_dir: str | None = None,
        required_columns: list[str] | None = None,
    ) -> dict:
        """
        Run neighborhood-based LRI analysis for multiple AnnData objects.

        Parameters
        ----------
        adata_dict : Dict[str, anndata.AnnData]
            Dictionary mapping sample_id -> AnnData
        output_dir : str, optional
            Directory to save results
        required_columns : list of str, optional
            If provided, subset to these columns

        Returns
        -------
        results : dict
            Dictionary containing combined results from all samples
        """
        logger.debug("=" * 60)
        logger.debug("Starting neighborhood-based LRI analysis (multi-sample mode)")
        logger.debug("=" * 60)
        logger.debug(f"Number of samples: {len(adata_dict)}")
        logger.debug(f"Neighborhood size: {self.neighborhood_size} µm")
        for sample_id, adata in adata_dict.items():
            logger.debug(
                f"  {sample_id}: {adata.shape[0]} cells, {adata.shape[1]} genes"
            )

        # Step 1: Get shared genes (intersection)
        logger.debug("[Step 1/7] Computing gene intersection...")
        shared_genes = self._get_shared_genes(adata_dict)

        # Step 2: Get all cell types (union)
        logger.debug("[Step 2/7] Computing cell type union...")
        shared_cell_types = self._get_shared_cell_types(adata_dict)
        self.cell_types = shared_cell_types
        logger.debug(f"Total unique cell types: {len(shared_cell_types)}")

        # Step 3: Prepare LRI database using shared genes
        logger.debug("[Step 3/7] Preparing LRI database with shared genes...")
        lr_pairs, ligand_genes_list, receptor_genes_list, signaling_types = (
            self.prepare_lri_database(gene_names=shared_genes)
        )

        # Step 4: Create unified column structure
        logger.debug("[Step 4/7] Creating unified column structure...")
        column_names = self._create_column_structure_from_cell_types(
            shared_cell_types, signaling_types
        )
        self.column_names = column_names
        logger.debug(f"Total LRI columns: {len(column_names)}")

        # Step 5: Process each sample
        logger.debug("[Step 5/7] Processing each sample...")
        all_matrices = []
        sample_info = {}
        global_cell_idx = 0

        for sample_id, adata in adata_dict.items():
            logger.debug(f"--- Processing sample: {sample_id} ---")

            # Build neighborhoods for this sample
            neighborhoods = self.build_neighborhoods(adata)

            # Build cell-LRI matrix using the pre-defined column structure
            sample_matrix = self.build_cell_lri_matrix(adata, signaling_types)

            n_cells = adata.n_obs

            # Store sample info
            sample_info[sample_id] = {
                "n_cells": n_cells,
                "global_cell_idx_start": global_cell_idx,
                "global_cell_idx_end": global_cell_idx + n_cells - 1,
                "avg_neighborhood_size": np.mean(
                    [len(n) for n in neighborhoods.values()]
                ),
            }

            all_matrices.append(sample_matrix)
            global_cell_idx += n_cells

        # Step 6: Combine results
        logger.debug("[Step 6/7] Combining results...")

        # Vertically stack all matrices
        combined_matrix = sparse_vstack(all_matrices, format="csr")
        logger.debug(f"Combined matrix shape: {combined_matrix.shape}")

        # Handle gene expression (if enabled)
        if self.include_gene_expression:
            logger.debug("Appending gene expression counts to matrix...")
            all_gene_matrices = []
            for sample_id, adata in adata_dict.items():
                raw_counts = self.get_raw_counts(adata)
                if not sp.issparse(raw_counts):
                    raw_counts = sp.csr_matrix(raw_counts)
                else:
                    raw_counts = raw_counts.tocsr()
                all_gene_matrices.append(raw_counts)

            combined_gene_matrix = sparse_vstack(all_gene_matrices, format="csr")
            combined_matrix = sp.hstack(
                [combined_matrix, combined_gene_matrix], format="csr"
            )

            # Use shared genes for column names
            gene_column_names = [f"GENE{self.spliter}{gene}" for gene in shared_genes]
            all_column_names = column_names + gene_column_names
            logger.debug(f"  - LRI features: {len(column_names)}")
            logger.debug(f"  - Gene features: {len(gene_column_names)}")
        else:
            all_column_names = column_names

        # Step 7: Column filtering
        logger.debug("[Step 7/7] Filtering columns...")
        if required_columns is not None:
            logger.debug(f"Subsetting to {len(required_columns)} required columns...")
            col_to_idx = {name: i for i, name in enumerate(all_column_names)}
            common_cols = [name for name in required_columns if name in col_to_idx]
            col_indices = np.array(
                [col_to_idx[name] for name in common_cols], dtype=int
            )

            logger.debug(f"  Cell matrix columns: {len(all_column_names)}")
            logger.debug(f"  Required columns: {len(required_columns)}")
            logger.debug(f"  Common columns (kept): {len(common_cols)}")

            if len(common_cols) < len(required_columns):
                missing = len(required_columns) - len(common_cols)
                logger.debug(f"  Warning: {missing} required columns not found")

            combined_matrix = combined_matrix[:, col_indices]
            all_column_names = common_cols
        else:
            combined_matrix, all_column_names = self.remove_zero_columns(
                combined_matrix, all_column_names
            )

        # Save results (optional)
        if output_dir is not None:
            logger.debug("Saving results...")
            os.makedirs(output_dir, exist_ok=True)

            # Save sparse matrix
            matrix_file = os.path.join(output_dir, "cell_lri_matrix.npz")
            sparse.save_npz(matrix_file, combined_matrix)

            # Save column names
            columns_file = os.path.join(output_dir, "cell_lri_columns.csv")
            pd.DataFrame({"column_name": all_column_names}).to_csv(
                columns_file, index=False
            )

            # Save sample info
            sample_info_file = os.path.join(output_dir, "sample_info.csv")
            sample_info_df = pd.DataFrame(
                [{"sample_id": sid, **info} for sid, info in sample_info.items()]
            )
            sample_info_df.to_csv(sample_info_file, index=False)

            # Save analysis parameters
            params_file = os.path.join(output_dir, "analysis_parameters.csv")
            params_df = pd.DataFrame(
                {
                    "parameter": [
                        "neighborhood_size",
                        "resource_name",
                        "n_samples",
                        "total_cells",
                        "n_lri_combinations",
                        "n_shared_genes",
                        "matrix_sparsity",
                        "include_gene_expression",
                    ],
                    "value": [
                        self.neighborhood_size,
                        self.resource_name,
                        len(adata_dict),
                        combined_matrix.shape[0],
                        len(column_names),
                        len(shared_genes),
                        f"{(1 - combined_matrix.nnz / np.prod(combined_matrix.shape)) * 100:.2f}%",
                        self.include_gene_expression,
                    ],
                }
            )
            params_df.to_csv(params_file, index=False)

            logger.debug(f"Results saved to: {output_dir}")
            logger.debug(f"- Cell-LRI matrix: {matrix_file}")
            logger.debug(f"- Column names: {columns_file}")
            logger.debug(f"- Sample info: {sample_info_file}")
            logger.debug(f"- Analysis parameters: {params_file}")

        logger.debug("=" * 60)
        logger.debug("Multi-sample analysis complete!")
        logger.debug(f"Total cells: {combined_matrix.shape[0]}")
        logger.debug(f"Total LRI columns: {combined_matrix.shape[1]}")
        logger.debug("=" * 60)

        return {
            "cell_lri_matrix": combined_matrix,
            "column_names": all_column_names,
            "sample_info": sample_info,
        }

    def _create_column_structure_from_cell_types(
        self, cell_types: list[str], signaling_types: list[str]
    ) -> list[str]:
        """
        Create column names using a predefined list of cell types.
        """
        column_names = []

        for idx, (lig, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]

            for lig_ct in cell_types:
                for rec_ct in cell_types:
                    if sig_type == "Cell-Cell Contact":
                        column_names.append(
                            f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}juxtacrine"
                        )
                    else:
                        if lig_ct == rec_ct:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}autocrine"
                            )
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )
                        else:
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )

        return column_names


# Alias for backward compatibility
SingleCellLRIAnalyzer = NeighborhoodLRIAnalyzer
