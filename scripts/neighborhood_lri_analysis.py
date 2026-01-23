"""
Cell Neighborhood-based Ligand-Receptor Interaction Analysis

This module implements spatial neighborhood-based LRI analysis where each cell
is the center of a square neighborhood (neighborhood_size x neighborhood_size).
It counts all-to-all interactions within each cell's neighborhood.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse import coo_matrix, csr_matrix
from liana.resource import select_resource
import anndata
from pathlib import Path
import os
from typing import Tuple, List, Dict, Optional
import scipy.sparse as sp
from sklearn.neighbors import KDTree


class NeighborhoodLRIAnalyzer:
    """
    Cell Neighborhood-based Ligand-Receptor Interaction Analyzer
    
    For each cell, defines a square neighborhood centered on that cell and counts
    ligand-receptor interactions within the neighborhood using all-to-all strategy.
    """
    
    def __init__(self, neighborhood_size: float = 50.0, resource_name: str = 'cellchatdb', 
                 spliter: str = '|',cellchatdb_path: str = 'data/LRdatabase/CellChatDBv2.0.human.csv', 
                 cellphonedb_path: str = 'data/LRdatabase/CellPhoneDBv5.0.human.csv', 
                 cell_type_column: str = 'cell_type',
                 include_gene_expression: bool = True,
                 raw_count_location: str = 'X'):
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
        """
        self.neighborhood_size = neighborhood_size
        self.resource_name = resource_name
        self.lr_pairs = None
        self.cell_types = None
        self.column_names = None
        self.neighborhoods = None
        self.cell_lri_matrix = None
        self.spliter = spliter
        self.receptor_genes_list = None
        self.cellphonedb_path = cellphonedb_path
        self.cellchatdb_path = cellchatdb_path
        self.cell_type_column = cell_type_column
        self.include_gene_expression = include_gene_expression  # 新增
        self.raw_count_location = raw_count_location  # 新增
    
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
        if self.raw_count_location == 'X':
            return adata.X
        elif self.raw_count_location == 'raw':
            if adata.raw is None:
                raise ValueError("adata.raw is None but raw_count_location='raw'")
            return adata.raw.X
        elif self.raw_count_location.startswith('layers:'):
            layer_name = self.raw_count_location.split(':', 1)[1]
            if layer_name not in adata.layers:
                raise ValueError(f"Layer '{layer_name}' not found in adata.layers")
            return adata.layers[layer_name]
        else:
            raise ValueError(f"Invalid raw_count_location: {self.raw_count_location}")
    
    def build_neighborhoods(self, adata: anndata.AnnData) -> Dict[int, np.ndarray]:
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
        print(f"Building neighborhoods (size={self.neighborhood_size}µm)...")
        coords = adata.obsm['spatial'][:, :2]
        n_cells = len(coords)
        
        # Build KD-Tree for efficient spatial queries
        tree = KDTree(coords)
        
        # Define half-width of square neighborhood
        half_size = self.neighborhood_size / 2.0
        
        neighborhoods = {}
        for i in range(n_cells):
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
                (candidate_coords[:, 0] >= x_min) & 
                (candidate_coords[:, 0] <= x_max) &
                (candidate_coords[:, 1] >= y_min) & 
                (candidate_coords[:, 1] <= y_max)
            )
            
            neighborhoods[i] = candidate_indices[in_square]
            
            if (i + 1) % 1000 == 0:
                print(f"  Processed {i + 1}/{n_cells} cells")
        
        self.neighborhoods = neighborhoods
        print(f"Average neighborhood size: {np.mean([len(n) for n in neighborhoods.values()]):.1f} cells")
        return neighborhoods
    
    def prepare_lri_database(self, adata: anndata.AnnData) -> Tuple[List[Tuple], List[Tuple], List[List[str]], List[str]]:
        """
        Prepare ligand-receptor pairs from database.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
            
        Returns
        -------
        lr_pairs : List[Tuple]
            List of (ligand, receptor_str) tuples
        receptor_genes_list : List[List[str]]
            List of receptor gene lists for each LR pair
        signaling_types : List[str]
            List of signaling types for each LR pair
        """
        print(f"Loading {self.resource_name} database...")
        
        # Load from local CSV if cellchatdb or cellphonedb
        if self.resource_name.lower() == 'cellchatdb':
            resource = pd.read_csv(self.cellchatdb_path)
            # Check required columns
            required_cols = ['ligand', 'receptor', 'signaling_type']
            if not all(col in resource.columns for col in required_cols):
                raise ValueError(
                    f"CellChatDB CSV must contain {required_cols}. "
                    f"Found columns: {resource.columns.tolist()}"
                )
        elif self.resource_name.lower() == 'cellphonedb':
            resource = pd.read_csv(self.cellphonedb_path)
            # Check required columns
            required_cols = ['ligand', 'receptor', 'signaling_type']
            if not all(col in resource.columns for col in required_cols):
                raise ValueError(
                    f"CellPhoneDB CSV must contain {required_cols}. "
                    f"Found columns: {resource.columns.tolist()}"
                )
        else:
            # Use liana's select_resource for other databases
            resource = select_resource(self.resource_name)
            # LIANA doesn't have signaling_type, add it as 'Unknown'
            if 'signaling_type' not in resource.columns:
                resource['signaling_type'] = 'Unknown'
        
        # Filter pairs where ligand exists and ALL receptor genes exist
        lr_pairs = []
        ligand_genes_list = []  #new
        receptor_genes_list = []
        signaling_types = []
        
        for idx in range(len(resource)):
            ligand = resource.iloc[idx]['ligand']
            receptor_str = resource.iloc[idx]['receptor']
            signaling_type = resource.iloc[idx]['signaling_type']
            
            # Skip if ligand or receptor is NaN
            if pd.isna(ligand) or pd.isna(receptor_str):
                continue
            
            # Convert to string
            ligand_str = str(ligand)
            receptor_str = str(receptor_str)
            signaling_type = str(signaling_type) if not pd.isna(signaling_type) else 'Unknown'
            
            # Parse both ligand and receptor genes (split by '_')
            ligand_genes = ligand_str.split('_')
            receptor_genes = receptor_str.split('_')
            
            # Check if ligand and ALL receptor genes exist in data
            if (all(lig_gene in adata.var_names for lig_gene in ligand_genes) and 
                all(rec_gene in adata.var_names for rec_gene in receptor_genes)):
                lr_pairs.append((ligand, receptor_str))
                ligand_genes_list.append(ligand_genes)
                receptor_genes_list.append(receptor_genes)
                signaling_types.append(signaling_type)
        
        print(f"Initial L-R pairs in data: {len(lr_pairs)}")
        print(f"  Single receptor: {sum(1 for rg in receptor_genes_list if len(rg) == 1)}")
        print(f"  Multi receptor: {sum(1 for rg in receptor_genes_list if len(rg) > 1)}")
        
        # Print signaling type distribution
        from collections import Counter
        sig_type_counts = Counter(signaling_types)
        print(f"\nSignaling type distribution:")
        for sig_type, count in sig_type_counts.items():
            print(f"  {sig_type}: {count}")
        
        self.lr_pairs = lr_pairs
        self.ligand_genes_list = ligand_genes_list
        self.receptor_genes_list = receptor_genes_list
        self.signaling_types = signaling_types  # Store as instance variable
        
        return lr_pairs,ligand_genes_list, receptor_genes_list, signaling_types
    
    def create_column_structure(self, adata: anndata.AnnData, signaling_types: List[str]) -> List[str]:
        """
        Create column names for the cell-LRI matrix.
        Must match the exact logic in build_cell_lri_matrix().
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
        signaling_types : List[str]
            List of signaling types for each LR pair
            
        Returns
        -------
        column_names : List[str]
            List of column names in format: ligand_ct|receptor_ct|ligand|receptor|mode
        """
        self.cell_types = adata.obs['cell_type'].cat.categories.tolist()
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}
        column_names = []
        
        # Must match the exact filtering and logic in build_cell_lri_matrix
        for idx, (lig_str, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]
            
            # Check if all ligand genes exist (支持配体复合体)
            lig_genes = lig_str.split('_')
            if any(lg not in gene_to_idx for lg in lig_genes):
                continue  # Skip this LR pair entirely
                
            # Check if all receptor genes exist
            rec_genes = rec_str.split('_')
            if any(rg not in gene_to_idx for rg in rec_genes):
                continue  # Skip this LR pair entirely
            
            # Now create columns based on signaling type
            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == 'Cell-Cell Contact':
                        # Only juxtacrine for contact-dependent
                        column_names.append(
                            f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}"
                            f"{lig_str}{self.spliter}{rec_str}{self.spliter}juxtacrine"
                        )
                    else:
                        # Non-contact signaling
                        if lig_ct == rec_ct:
                            # Both autocrine and paracrine for same cell type
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}"
                                f"{lig_str}{self.spliter}{rec_str}{self.spliter}autocrine"
                            )
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}"
                                f"{lig_str}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )
                        else:
                            # Only paracrine for different cell types
                            column_names.append(
                                f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}"
                                f"{lig_str}{self.spliter}{rec_str}{self.spliter}paracrine"
                            )
        
        self.column_names = column_names
        return column_names
    
    def build_cell_lri_matrix(self, adata: anndata.AnnData, signaling_types: List[str]) -> csr_matrix:
        """
        Build the cell-LRI interaction matrix with full vectorization.
        Supports both ligand and receptor complexes (using AND logic for co-expression).
        """
        print("Building cell-LRI matrix with vectorized neighborhood interactions...")

        n_cells = adata.n_obs
        
        # ─── 1) Index mappings ────────────────────────────────────────────────────
        ct_to_idx = {ct: i for i, ct in enumerate(self.cell_types)}
        cell_types_idx = np.array(
            [ct_to_idx[ct] for ct in adata.obs[self.cell_type_column]],
            dtype=int
        )
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # ─── 2) Parse column metadata - NOW HANDLES LIGAND COMPLEXES ─────────────
        col_meta = []
        for idx, (lig_str, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]
            
            # Check if all ligand genes exist (支持配体复合体)
            lig_genes = lig_str.split('_')
            if any(lg not in gene_to_idx for lg in lig_genes):
                continue
                
            # Check if all receptor genes exist
            rec_genes = rec_str.split('_')
            if any(rg not in gene_to_idx for rg in rec_genes):
                continue

            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == 'Cell-Cell Contact':
                        # 一列：juxtacrine
                        col_meta.append((
                            len(col_meta),             # column j
                            ct_to_idx[lig_ct],        # lig_ct_idx
                            ct_to_idx[rec_ct],        # rec_ct_idx
                            lig_str,                  # ligand string (可能是复合体)
                            rec_str,                  # receptor string (可能是复合体)
                            'juxtacrine',             # mode
                            sig_type                  # signaling type
                        ))
                    else:
                        # 非接触
                        if lig_ct == rec_ct:
                            # 两列：autocrine + paracrine
                            col_meta.append((
                                len(col_meta),
                                ct_to_idx[lig_ct],
                                ct_to_idx[rec_ct],
                                lig_str,
                                rec_str,
                                'autocrine',
                                sig_type
                            ))
                            col_meta.append((
                                len(col_meta),
                                ct_to_idx[lig_ct],
                                ct_to_idx[rec_ct],
                                lig_str,
                                rec_str,
                                'paracrine',
                                sig_type
                            ))
                        else:
                            # 一列：paracrine
                            col_meta.append((
                                len(col_meta),
                                ct_to_idx[lig_ct],
                                ct_to_idx[rec_ct],
                                lig_str,
                                rec_str,
                                'paracrine',
                                sig_type
                            ))

        n_columns = len(col_meta)
        print(f"Processing {n_cells} cells × {n_columns} LRI combinations")

        # ─── 3) Binarize expression ───────────────────────────────────────────────
        X = adata.X
        if sp.issparse(X):
            expr_bool = (X > 0).astype(int).tocsc()
        else:
            expr_bool = csr_matrix((X > 0).astype(int)).tocsc()

        expr_coo = expr_bool.tocoo()

        # ─── 4) Build cell-neighborhood adjacency matrix ──────────────────────────
        print("Building cell-neighborhood adjacency matrix...")
        rows, cols = [], []
        for cell_idx, neighbors in self.neighborhoods.items():
            rows.extend([cell_idx] * len(neighbors))
            cols.extend(neighbors)
        cell_nbr_matrix = coo_matrix(
            (np.ones(len(rows), dtype=int), (rows, cols)),
            shape=(n_cells, n_cells)
        ).tocsr()

        # ─── 5) Build neighborhood_by_lig for INDIVIDUAL ligand genes ─────────────
        neighborhood_by_lig = {}
        print("Building neighborhood-by-ligand matrices (individual genes)...")

        for ct_idx in range(len(self.cell_types)):
            mask_cells = (cell_types_idx == ct_idx)
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (expr_coo.data[entry_mask],
                (expr_coo.row[entry_mask], expr_coo.col[entry_mask])),
                shape=expr_bool.shape
            )

            # Collect all individual ligand genes for this cell type (展开复合体)
            all_lig_genes_ct = set()
            for (_, lct, _, lig_str, _, _, _) in col_meta:
                if lct == ct_idx:
                    all_lig_genes_ct.update(lig_str.split('_'))
            all_lig_genes_ct = sorted([g for g in all_lig_genes_ct if g in gene_to_idx])

            # Build matrix for each individual ligand gene
            lig_gene_matrices = {}
            for lig_gene in all_lig_genes_ct:
                lig_gene_idx = gene_to_idx[lig_gene]
                lig_mask = (sub.col == lig_gene_idx)
                rows_cells = sub.row[lig_mask]
                data_vals_ct = sub.data[lig_mask]

                cell_lig = coo_matrix(
                    (data_vals_ct, (rows_cells, np.zeros(len(rows_cells), dtype=int))),
                    shape=(n_cells, 1)
                )
                cell_lig.sum_duplicates()
                # 邻域聚合
                lig_gene_matrices[lig_gene] = cell_nbr_matrix.dot(cell_lig.tocsr())

            neighborhood_by_lig[ct_idx] = lig_gene_matrices

        # ─── 6) Build neighborhood_by_rec for INDIVIDUAL receptor genes ───────────
        neighborhood_by_rec = {}
        print("Building neighborhood-by-receptor matrices (individual genes)...")

        for ct_idx in range(len(self.cell_types)):
            mask_cells = (cell_types_idx == ct_idx)
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (expr_coo.data[entry_mask],
                (expr_coo.row[entry_mask], expr_coo.col[entry_mask])),
                shape=expr_bool.shape
            )

            # Collect all individual receptor genes for this cell type
            all_rec_genes_ct = set()
            for (_, _, rct, _, rec_str, _, _) in col_meta:
                if rct == ct_idx:
                    all_rec_genes_ct.update(rec_str.split('_'))
            all_rec_genes_ct = sorted([g for g in all_rec_genes_ct if g in gene_to_idx])

            # Build matrix for each individual receptor gene
            rec_gene_matrices = {}
            for rec_gene in all_rec_genes_ct:
                rec_gene_idx = gene_to_idx[rec_gene]
                rec_mask = (sub.col == rec_gene_idx)
                rows_cells = sub.row[rec_mask]
                data_vals_ct = sub.data[rec_mask]

                cell_rec = coo_matrix(
                    (data_vals_ct, (rows_cells, np.zeros(len(rows_cells), dtype=int))),
                    shape=(n_cells, 1)
                )
                cell_rec.sum_duplicates()
                # 邻域聚合
                rec_gene_matrices[rec_gene] = cell_nbr_matrix.dot(cell_rec.tocsr())

            neighborhood_by_rec[ct_idx] = rec_gene_matrices

        # ─── 7) Compute interactions with BOTH ligand and receptor complexes ──────
        print("Computing LRI interactions...")
        row_inds, col_inds, data_vals = [], [], []

        for j, lig_ct_idx, rec_ct_idx, lig_str, rec_str, mode, sig_type in col_meta:
            if j % 500 == 0:
                print(f"  Progress: {j}/{n_columns}")

            # Ligand count (支持配体复合体，AND逻辑)
            lig_genes = lig_str.split('_')
            if lig_genes[0] in neighborhood_by_lig[lig_ct_idx]:
                count_lig = np.array(
                    neighborhood_by_lig[lig_ct_idx][lig_genes[0]].toarray()
                ).ravel()
                for lig_gene in lig_genes[1:]:  # 如果是复合体
                    if lig_gene in neighborhood_by_lig[lig_ct_idx]:
                        count_lig = np.minimum(
                            count_lig,
                            np.array(neighborhood_by_lig[lig_ct_idx][lig_gene].toarray()).ravel()
                        )
                    else:
                        count_lig = np.zeros(n_cells, dtype=int)
                        break
            else:
                count_lig = np.zeros(n_cells, dtype=int)

            # Receptor count (支持受体复合体，AND逻辑)
            rec_genes = rec_str.split('_')
            if rec_genes[0] in neighborhood_by_rec[rec_ct_idx]:
                count_rec = np.array(
                    neighborhood_by_rec[rec_ct_idx][rec_genes[0]].toarray()
                ).ravel()
                for rec_gene in rec_genes[1:]:
                    if rec_gene in neighborhood_by_rec[rec_ct_idx]:
                        count_rec = np.minimum(
                            count_rec,
                            np.array(neighborhood_by_rec[rec_ct_idx][rec_gene].toarray()).ravel()
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
                    coexpr = coexpr * expr_bool[:, lig_gene_idx].toarray().ravel().astype(int)
                
                # All receptor genes must be co-expressed
                for rec_gene in rec_genes:
                    rec_gene_idx = gene_to_idx[rec_gene]
                    coexpr = coexpr * expr_bool[:, rec_gene_idx].toarray().ravel().astype(int)
                
                # Only cells of the correct type
                coexpr = coexpr * (cell_types_idx == lig_ct_idx).astype(int)
                auto = np.array(cell_nbr_matrix.dot(coexpr)).ravel()
            else:
                auto = np.zeros(n_cells, dtype=int)

            # Write interactions based on mode
            if mode == "juxtacrine":
                total = count_lig * count_rec
                rows_nz = np.nonzero(total)[0]
                row_inds.extend(rows_nz.tolist())
                col_inds.extend([j] * len(rows_nz))
                data_vals.extend(total[rows_nz].tolist())

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
            (data_vals, (row_inds, col_inds)),
            shape=(n_cells, n_columns),
            dtype=int
        )
        self.cell_lri_matrix = cell_lri_matrix
        print(f"Matrix density: {cell_lri_matrix.nnz / (n_cells * n_columns) * 100:.2f}%")
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
        print("Creating metadata dataframe...")
        
        coords = adata.obsm['spatial'][:, :2]
        neighborhood_sizes = [len(self.neighborhoods[i]) for i in range(adata.n_obs)]
        
        cell_metadata_df = pd.DataFrame({
            'cell_id': adata.obs.index,
            # 'tma_id': adata.obs['tma_id'],
            'cell_type': adata.obs['cell_type'],
            'x_coord': coords[:, 0],
            'y_coord': coords[:, 1],
            'neighborhood_size': neighborhood_sizes
        })
        
        return cell_metadata_df
    
    def run_analysis(self, adata: anndata.AnnData, output_dir: str) -> Dict:
        """
        Run the complete neighborhood-based LRI analysis.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
        output_dir : str
            Directory to save results
            
        Returns
        -------
        results : dict
            Dictionary containing all analysis results
        """
        print("Starting cell neighborhood-based LRI analysis...")
        print(f"Neighborhood size: {self.neighborhood_size} µm")
        print(f"Data shape: {adata.shape}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Step 1: Build neighborhoods
        neighborhoods = self.build_neighborhoods(adata)
        
        # Step 2: Prepare LRI database
        lr_pairs,ligand_genes_list, receptor_genes_list, signaling_types = self.prepare_lri_database(adata)
        
        # Step 3: Create column structure
        column_names = self.create_column_structure(adata, signaling_types)
        
        # Step 4: Build cell-LRI matrix
        cell_lri_matrix = self.build_cell_lri_matrix(adata, signaling_types)

        # Step 5: Optionally append gene expression (新增部分)
        if self.include_gene_expression:
            print("Appending gene expression counts to matrix...")
            
            # Get raw counts
            raw_counts = self.get_raw_counts(adata)
            print('Using raw counts')
            
            # Convert to sparse CSR if needed
            if not sp.issparse(raw_counts):
                raw_counts = sp.csr_matrix(raw_counts)
            else:
                raw_counts = raw_counts.tocsr()
            
            # Concatenate horizontally: [LRI features | gene expression]
            combined_matrix = sp.hstack([cell_lri_matrix, raw_counts], format='csr')
            
            # Update column names
            gene_column_names = [f"GENE{self.spliter}{gene}" for gene in adata.var_names]
            all_column_names = column_names + gene_column_names
            
            print(f"Combined matrix shape: {combined_matrix.shape}")
            print(f"  - LRI features: {len(column_names)}")
            print(f"  - Gene features: {len(gene_column_names)}")
        else:
            combined_matrix = cell_lri_matrix
            all_column_names = column_names
        
        # Step 6: Create metadata dataframe
        cell_metadata_df = self.create_metadata_dataframe(adata)
        
        # Step 7: Save results
        print("Saving results...")
        
        # Save sparse matrix
        matrix_file = os.path.join(output_dir, 'cell_lri_matrix.npz')
        sparse.save_npz(matrix_file, combined_matrix)
        
        # Save column names
        columns_file = os.path.join(output_dir, 'cell_lri_columns.csv')
        pd.DataFrame({'column_name': all_column_names}).to_csv(columns_file, index=False)
        
        # Save metadata
        metadata_file = os.path.join(output_dir, 'cell_metadata.csv')
        cell_metadata_df.to_csv(metadata_file, index=False)
        
        # Save analysis parameters
        params_file = os.path.join(output_dir, 'analysis_parameters.csv')
        params_df = pd.DataFrame({
            'parameter': ['neighborhood_size', 'resource_name', 'n_cells', 
                         'n_lri_combinations', 'n_genes', 'n_total_features',
                         'matrix_sparsity', 'avg_neighborhood_size', 
                         'include_gene_expression', 'raw_count_location'],
            'value': [
                self.neighborhood_size,
                self.resource_name,
                adata.n_obs,
                len(column_names),
                len(adata.var_names) if self.include_gene_expression else 0,
                len(all_column_names),
                f"{(1 - combined_matrix.nnz / np.prod(combined_matrix.shape)) * 100:.2f}%",
                f"{cell_metadata_df['neighborhood_size'].mean():.1f}",
                self.include_gene_expression,
                self.raw_count_location if self.include_gene_expression else 'N/A'
            ]
        })
        params_df.to_csv(params_file, index=False)
        
        print(f"Results saved to: {output_dir}")
        print(f"- Cell-LRI matrix: {matrix_file}")
        print(f"- Column names: {columns_file}")
        print(f"- Cell metadata: {metadata_file}")
        print(f"- Analysis parameters: {params_file}")
        
        return {
            'cell_lri_matrix': combined_matrix,
            'column_names': all_column_names,
            'cell_metadata_df': cell_metadata_df,
            'neighborhoods': neighborhoods,
            'lr_pairs': lr_pairs
        }

    
    # def create_metadata_dataframe(self, adata: anndata.AnnData) -> pd.DataFrame:
    #     """
    #     Create metadata dataframe with cell information.
        
    #     Parameters
    #     ----------
    #     adata : anndata.AnnData
    #         Spatial transcriptomics data
            
    #     Returns
    #     -------
    #     cell_metadata_df : pd.DataFrame
    #         DataFrame with cell metadata
    #     """
    #     print("Creating metadata dataframe...")
        
    #     coords = adata.obsm['spatial'][:, :2]
    #     neighborhood_sizes = [len(self.neighborhoods[i]) for i in range(adata.n_obs)]
        
    #     cell_metadata_df = pd.DataFrame({
    #         'cell_id': adata.obs['cell_id'],
    #         'tma_id': adata.obs['tma_id'],
    #         'cell_type': adata.obs['cell_type'],
    #         'x_coord': coords[:, 0],
    #         'y_coord': coords[:, 1],
    #         'neighborhood_size': neighborhood_sizes
    #     })
        
    #     return cell_metadata_df
    
    # def run_analysis(self, adata: anndata.AnnData, output_dir: str) -> Dict:
    #     """
    #     Run the complete neighborhood-based LRI analysis.
        
    #     Parameters
    #     ----------
    #     adata : anndata.AnnData
    #         Spatial transcriptomics data
    #     output_dir : str
    #         Directory to save results
            
    #     Returns
    #     -------
    #     results : dict
    #         Dictionary containing all analysis results
    #     """
    #     print("Starting cell neighborhood-based LRI analysis...")
    #     print(f"Neighborhood size: {self.neighborhood_size} µm")
    #     print(f"Data shape: {adata.shape}")
        
    #     # Create output directory
    #     os.makedirs(output_dir, exist_ok=True)
        
    #     # Step 1: Build neighborhoods
    #     neighborhoods = self.build_neighborhoods(adata)
        
    #     # Step 2: Prepare LRI database
    #     lr_pairs, receptor_genes_list, signaling_types = self.prepare_lri_database(adata)
        
    #     # Step 3: Create column structure
    #     column_names = self.create_column_structure(adata)
        
    #     # Step 4: Build cell-LRI matrix
    #     cell_lri_matrix = self.build_cell_lri_matrix(adata, signaling_types)
        
    #     # Step 5: Create metadata dataframe
    #     cell_metadata_df = self.create_metadata_dataframe(adata)
        
    #     # Step 6: Save results
    #     print("Saving results...")
        
    #     # Save sparse matrix
    #     matrix_file = os.path.join(output_dir, 'cell_lri_matrix.npz')
    #     sparse.save_npz(matrix_file, cell_lri_matrix)
        
    #     # Save column names
    #     columns_file = os.path.join(output_dir, 'cell_lri_columns.csv')
    #     pd.DataFrame({'column_name': self.column_names}).to_csv(columns_file, index=False)
        
    #     # Save metadata
    #     metadata_file = os.path.join(output_dir, 'cell_metadata.csv')
    #     cell_metadata_df.to_csv(metadata_file, index=False)
        
    #     # Save analysis parameters
    #     params_file = os.path.join(output_dir, 'analysis_parameters.csv')
    #     params_df = pd.DataFrame({
    #         'parameter': ['neighborhood_size', 'resource_name', 'n_cells', 'n_lri_combinations', 'matrix_sparsity', 'avg_neighborhood_size'],
    #         'value': [
    #             self.neighborhood_size,
    #             self.resource_name,
    #             adata.n_obs,
    #             len(column_names),
    #             f"{(1 - cell_lri_matrix.nnz / np.prod(cell_lri_matrix.shape)) * 100:.2f}%",
    #             f"{cell_metadata_df['neighborhood_size'].mean():.1f}"
    #         ]
    #     })
    #     params_df.to_csv(params_file, index=False)
        
    #     print(f"Results saved to: {output_dir}")
    #     print(f"- Cell-LRI matrix: {matrix_file}")
    #     print(f"- Column names: {columns_file}")
    #     print(f"- Cell metadata: {metadata_file}")
    #     print(f"- Analysis parameters: {params_file}")
        
    #     return {
    #         'cell_lri_matrix': cell_lri_matrix,
    #         'column_names': column_names,
    #         'cell_metadata_df': cell_metadata_df,
    #         'neighborhoods': neighborhoods,
    #         'lr_pairs': lr_pairs
    #     }


def load_cell_lri_results(output_dir: str) -> Dict:
    """
    Load previously saved cell-LRI analysis results.
    
    Parameters
    ----------
    output_dir : str
        Directory containing saved results
        
    Returns
    -------
    results : dict
        Dictionary containing loaded results
    """
    print(f"Loading cell-LRI results from: {output_dir}")
    
    # Load sparse matrix
    matrix_file = os.path.join(output_dir, 'cell_lri_matrix.npz')
    cell_lri_matrix = sparse.load_npz(matrix_file)
    
    # Load column names
    columns_file = os.path.join(output_dir, 'cell_lri_columns.csv')
    column_names = pd.read_csv(columns_file)['column_name'].tolist()
    
    # Load metadata
    metadata_file = os.path.join(output_dir, 'cell_metadata.csv')
    cell_metadata_df = pd.read_csv(metadata_file)
    
    # Load parameters
    params_file = os.path.join(output_dir, 'analysis_parameters.csv')
    params_df = pd.read_csv(params_file)
    
    print(f"Loaded matrix shape: {cell_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")
    
    return {
        'cell_lri_matrix': cell_lri_matrix,
        'column_names': column_names,
        'cell_metadata_df': cell_metadata_df,
        'parameters': params_df
    }