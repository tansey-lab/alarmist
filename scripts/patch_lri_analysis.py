"""
Patch-based Ligand-Receptor Interaction Analysis

This module implements spatial patch-based LRI analysis for matrix factorization approaches.
It divides tissue into regular grid patches and counts all-to-all interactions within each patch.
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
# from numba import njit, prange
import scipy.sparse as sp
from scipy.sparse import kron, hstack


class PatchLRIAnalyzer:
    """
    Patch-based Ligand-Receptor Interaction Analyzer
    
    Divides spatial transcriptomics data into regular grid patches and counts
    ligand-receptor interactions within each patch using all-to-all strategy.
    """
    
    def __init__(self, patch_size: float = 50.0, resource_name: str = 'cellchatdb', spliter: str = '|'):
        """
        Initialize the patch-based LRI analyzer.
        
        Parameters
        ----------
        patch_size : float, default 50.0
            Size of spatial patches in micrometers
        resource_name : str, default 'cellchatdb'
            LRI database to use from liana
        """
        self.patch_size = patch_size
        self.resource_name = resource_name
        self.lr_pairs = None
        self.cell_types = None
        self.column_names = None
        self.patch_assignments = None
        self.patch_coords = None
        self.patch_lri_matrix = None
        self.spliter = spliter
        
    def create_spatial_patches(self, adata: anndata.AnnData) -> Tuple[np.ndarray, Dict]:
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
        coords = adata.obsm['spatial'][:, :2]
        
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
            'x_bins': x_bins,
            'y_bins': y_bins,
            'patch_coords': patch_coords,
            'n_patches': len(patch_coords)
        }
        
        self.patch_assignments = patch_assignments
        self.patch_coords = patch_coords
        
        return patch_assignments, patch_info
    
    def prepare_lri_database(self, adata: anndata.AnnData) -> List[Tuple]:
        """
        Prepare ligand-receptor pairs from CellChat database.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
            
        Returns
        -------
        lr_pairs : List[Tuple]
            List of (ligand, receptor) tuples present in the data
        """
        print(f"Loading {self.resource_name} database...")
        resource = select_resource(self.resource_name)
        
        # Filter L-R pairs present in the data
        lr_df = resource[
            (np.isin(resource.ligand, adata.var_names)) &
            (np.isin(resource.receptor, adata.var_names))
        ]
        
        print(f"Initial L-R pairs in data: {len(lr_df)}")
        
        lr_pairs = [(row.ligand, row.receptor) for _, row in lr_df.iterrows()]
        self.lr_pairs = lr_pairs
        return lr_pairs
    
    def create_column_structure(self, adata: anndata.AnnData) -> List[str]:
        """
        Create column names for the patch-LRI matrix.
        # ligand_ct==receptor_ct → two columns (auto + para); else one (para).
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
            
        Returns
        -------
        column_names : List[str]
            List of column names in format: ligand_ct-receptor_ct-ligand-receptor
        """
        self.cell_types = adata.obs['cell_type'].cat.categories.tolist()
        column_names = []
        
        for lig, rec in self.lr_pairs:
            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                   if lig_ct == rec_ct:
                       column_names.append(f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec}{self.spliter}autocrine")
                       column_names.append(f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec}{self.spliter}paracrine")
                   else:
                       column_names.append(f"{lig_ct}{self.spliter}{rec_ct}{self.spliter}{lig}{self.spliter}{rec}{self.spliter}paracrine")
        
        self.column_names = column_names
        return column_names
    
    
    def build_patch_lri_matrix(self, adata: anndata.AnnData, use_batch_processing: bool = True, batch_size: int = 500) -> sparse.csr_matrix:
        """
        Build the patch-LRI interaction matrix using efficient sparse operations.
        """
        print("Building patch-LRI matrix with sparse optimization...")

        # 1. Prepare basics
        unique_patches = np.array([p for p in np.unique(self.patch_assignments) if p in self.patch_coords])
        # map each global patch_id → local row index 0..n_patches-1
        patch_idx_map = {pid: idx for idx, pid in enumerate(unique_patches)}
        n_patches = len(unique_patches)
        n_columns = len(self.column_names)
        print(f"Processing {n_patches} patches × {n_columns} LRI combinations")

        # Map cell types and genes to indices
        ct_to_idx = {ct: i for i, ct in enumerate(self.cell_types)}
        cell_types_idx = np.array([ct_to_idx[ct] for ct in adata.obs['cell_type']], dtype=np.int64)
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # Prepare LRI indices
        lig_idxs, rec_idxs, lig_ct_idxs, rec_ct_idxs = [], [], [], []
        for col in self.column_names:
            lig_ct, rec_ct, lig, rec = col.split(self.spliter)
            lig_idxs.append(gene_to_idx[lig])
            rec_idxs.append(gene_to_idx[rec])
            lig_ct_idxs.append(ct_to_idx[lig_ct])
            rec_ct_idxs.append(ct_to_idx[rec_ct])
        lig_idxs = np.array(lig_idxs); rec_idxs = np.array(rec_idxs)
        lig_ct_idxs = np.array(lig_ct_idxs); rec_ct_idxs = np.array(rec_ct_idxs)
        n_valid = len(self.column_names)

        # 2. Binarize expression into COO
        expr_coo = sparse.coo_matrix(adata.X)
        expr_coo.data = np.ones_like(expr_coo.data)

        # 3. Build Patch×Ligand and Patch×Receptor per cell type
        patch_by_lig = {}
        patch_by_rec = {}
        print("Building patch-by-ligand and patch-by-receptor matrices...")
        for ct_idx in range(len(self.cell_types)):
            # Filter by cell type
            # 3.1 按细胞类型筛行：先得到每个非零条目是否属于该类型细胞
            mask_cells = (cell_types_idx == ct_idx)            # 长度 = n_cells
            entry_mask = mask_cells[expr_coo.row]              # 长度 = expr_coo.nnz

            # 3.2 构建只含该 ct_idx 行的子 COO
            sub = sparse.coo_matrix(
                (expr_coo.data[entry_mask],
                (expr_coo.row[entry_mask],
                expr_coo.col[entry_mask])),
                shape=expr_coo.shape
            )

            # 1) 先取出 sub 中这一 cell type 的 ligand 条目
            lig_mask = np.isin(sub.col, lig_idxs[lig_ct_idxs == ct_idx])
            row_cells = sub.row[lig_mask]           # 原 cell 索引
            col_genes = sub.col[lig_mask]           # 原基因全局索引
            data_vals = sub.data[lig_mask]
            
            # 2) map cell→local patch row index
            patch_rows = np.array([patch_idx_map[self.patch_assignments[c]] for c in row_cells], dtype=np.int64)
            
            # 3) map全局基因索引→ligand子集的 local col index
            lig_genes_ct = sorted(set(lig_idxs[lig_ct_idxs == ct_idx]))
            lig_to_local = {g:i for i,g in enumerate(lig_genes_ct)}
            local_cols = np.array([lig_to_local[g] for g in col_genes], dtype=np.int64)
            
            # 4) 构造 COO，并显式指定 shape=(n_patches, n_ligands_this_ct)
            lig_coo = sparse.coo_matrix(
                (data_vals, (patch_rows, local_cols)),
                shape=(n_patches, len(lig_genes_ct))
            )
            lig_coo.sum_duplicates()
            patch_by_lig[ct_idx] = lig_coo.tocsr()

            rec_mask = np.isin(sub.col, rec_idxs[rec_ct_idxs == ct_idx])
            row_cells = sub.row[rec_mask]
            col_genes = sub.col[rec_mask]
            data_vals = sub.data[rec_mask]
            
            patch_rows = np.array([patch_idx_map[self.patch_assignments[c]] for c in row_cells], dtype=np.int64)
            
            rec_genes_ct = sorted(set(rec_idxs[rec_ct_idxs == ct_idx]))
            rec_to_local = {g:i for i,g in enumerate(rec_genes_ct)}
            local_cols = np.array([rec_to_local[g] for g in col_genes], dtype=np.int64)
            
            rec_coo = sparse.coo_matrix(
                (data_vals, (patch_rows, local_cols)),
                shape=(n_patches, len(rec_genes_ct))
            )
            rec_coo.sum_duplicates()
            patch_by_rec[ct_idx] = rec_coo.tocsr()

        # 预缓存
        unique_ct = sorted(set(lig_ct_idxs))
        lig_ct2genes  = {ct: sorted(set(lig_idxs[lig_ct_idxs==ct]))  for ct in unique_ct}
        lig_ct2local  = {ct: {g:i for i,g in enumerate(lig_ct2genes[ct])} for ct in unique_ct}
        rec_ct2genes  = {ct: sorted(set(rec_idxs[rec_ct_idxs==ct]))  for ct in unique_ct}
        rec_ct2local  = {ct: {g:i for i,g in enumerate(rec_ct2genes[ct])} for ct in unique_ct}

        # 4. Compute interactions by Cartesian product per cell-type pair
        print("Computing predefined LRI combinations...")
        row_indices, col_indices, data_values = [], [], []

        for lri_idx in range(n_valid):
            if lri_idx % 500 == 0:  # 每500个组合打印一次进度
                print(f"  Progress: {lri_idx}/{n_valid}")
            
            lig_ct_idx = lig_ct_idxs[lri_idx]
            rec_ct_idx = rec_ct_idxs[lri_idx]

            Lm = patch_by_lig.get(lig_ct_idx)
            Rm = patch_by_rec.get(rec_ct_idx)
            if Lm is None or Rm is None:
                continue

            # 查本次组合在各自子矩阵里的列号
            lig_col = lig_ct2local[lig_ct_idx][ lig_idxs[lri_idx] ]
            rec_col = rec_ct2local[rec_ct_idx][ rec_idxs[lri_idx] ]

            # 取稀疏列并相乘
            lig_vec = Lm[:, lig_col]      # CSR slice
            rec_vec = Rm[:, rec_col]
            inter  = lig_vec.multiply(rec_vec)  # CSR

            nz_rows = inter.nonzero()[0]
            if nz_rows.size:
                row_indices.extend(nz_rows)
                col_indices.extend([lri_idx] * len(nz_rows))
                data_values.extend(inter.data)

        patch_lri_matrix = sparse.csr_matrix(
            (data_values, (row_indices, col_indices)),
            shape=(n_patches, n_valid),
            dtype=np.int64
        )
        self.patch_lri_matrix = patch_lri_matrix
        return patch_lri_matrix
    
    def build_patch_lri_matrix_with_mode(self,
                            adata: anndata.AnnData,
                            use_batch_processing: bool = True,
                            batch_size: int = 500) -> csr_matrix:
        """
        Build the patch-LRI interaction matrix, distinguishing true single-cell
        autocrine from paracrine interactions—while preserving your original
        vectorized speed in patch_by_lig/patch_by_rec.
        """
        print("Building patch-LRI matrix with autocrine/paracrine distinction...")

        # ─── 1) Prepare basics ────────────────────────────────────────────────────────
        unique_patches = np.array([
            p for p in np.unique(self.patch_assignments)
            if p in self.patch_coords
        ])
        patch_idx_map = {pid: i for i, pid in enumerate(unique_patches)}
        n_patches = len(unique_patches)
        n_columns = len(self.column_names)
        print(f"Processing {n_patches} patches × {n_columns} LRI combinations")

        # ─── 2) Index mappings ───────────────────────────────────────────────────────
        ct_to_idx   = {ct: i for i, ct in enumerate(self.cell_types)}
        cell_types_idx = np.array(
            [ct_to_idx[ct] for ct in adata.obs['cell_type']],
            dtype=int
        )
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # Parse column metadata: (col_idx, lig_ct_idx, rec_ct_idx, lig_gene_idx, rec_gene_idx, mode)
        col_meta = []
        for j, col in enumerate(self.column_names):
            lig_ct, rec_ct, lig, rec, mode = col.split(self.spliter)
            col_meta.append((
                j,
                ct_to_idx[lig_ct],
                ct_to_idx[rec_ct],
                gene_to_idx[lig],
                gene_to_idx[rec],
                mode
            ))

        # ─── 3) Binarize expression ───────────────────────────────────────────────────
        X = adata.X
        if sp.issparse(X):
            expr_bool = (X > 0).astype(int).tocsc()
        else:
            expr_bool = csr_matrix((X > 0).astype(int)).tocsc()

        expr_coo = expr_bool.tocoo()

        # ─── 4) Build patch_by_lig and patch_by_rec (vectorized) ───
        patch_by_lig = {}
        patch_by_rec = {}
        print("Building patch-by-ligand and patch-by-receptor matrices...")
        for ct_idx in range(len(self.cell_types)):
            # 3.1 filter according to cell type
            mask_cells = (cell_types_idx == ct_idx)
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix(
                (expr_coo.data[entry_mask],
                (expr_coo.row[entry_mask], expr_coo.col[entry_mask])),
                shape=expr_bool.shape
            )
            # Ligands
            lig_genes_ct = sorted({
                lig for (_, lct, _, lig, _, _) in col_meta if lct == ct_idx
            })
            lig_to_local = {g: i for i, g in enumerate(lig_genes_ct)}
            lig_mask = np.isin(sub.col, lig_genes_ct)
            rows_cells = sub.row[lig_mask]
            cols_genes = sub.col[lig_mask]
            data_vals_ct = sub.data[lig_mask]
            patch_rows = np.array([
                patch_idx_map[self.patch_assignments[c]]
                for c in rows_cells
            ], dtype=int)
            local_cols = np.array([lig_to_local[g] for g in cols_genes], dtype=int)
            lig_coo = coo_matrix(
                (data_vals_ct, (patch_rows, local_cols)),
                shape=(n_patches, len(lig_genes_ct))
            )
            lig_coo.sum_duplicates()
            patch_by_lig[ct_idx] = lig_coo.tocsr()

            # Receptors
            rec_genes_ct = sorted({
                rec for (_, _, rct, _, rec, _) in col_meta if rct == ct_idx
            })
            rec_to_local = {g: i for i, g in enumerate(rec_genes_ct)}
            rec_mask = np.isin(sub.col, rec_genes_ct)
            rows_cells = sub.row[rec_mask]
            cols_genes = sub.col[rec_mask]
            data_vals_ct = sub.data[rec_mask]
            patch_rows = np.array([
                patch_idx_map[self.patch_assignments[c]]
                for c in rows_cells
            ], dtype=int)
            local_cols = np.array([rec_to_local[g] for g in cols_genes], dtype=int)
            rec_coo = coo_matrix(
                (data_vals_ct, (patch_rows, local_cols)),
                shape=(n_patches, len(rec_genes_ct))
            )
            rec_coo.sum_duplicates()
            patch_by_rec[ct_idx] = rec_coo.tocsr()

        # ─── 5) Local-index caches ───────────────────────────────────────────────────
        lig_ct2local = {
            ct: {g: i for i, g in enumerate(sorted({
                lig for (_, lct, _, lig, _, _) in col_meta if lct == ct
            }))}
            for ct in range(len(self.cell_types))
        }
        rec_ct2local = {
            ct: {g: i for i, g in enumerate(sorted({
                rec for (_, _, rct, _, rec, _) in col_meta if rct == ct
            }))}
            for ct in range(len(self.cell_types))
        }

        # ─── 6) 构造 patch_by_cell：一次性稀疏矩阵，shape=(n_patches, n_cells) ───────────
        n_cells = adata.n_obs
        cell_patches = np.array(self.patch_assignments, dtype=int)
        pb_rows = [patch_idx_map[p] for p in cell_patches]
        pb_cols = list(range(n_cells))
        patch_by_cell = coo_matrix(
            (np.ones(n_cells, int), (pb_rows, pb_cols)),
            shape=(n_patches, n_cells)
        ).tocsr()

        # ─── 7) Compute autocrine/paracrine ───────────────────────────────────────────
        row_inds, col_inds, data_vals = [], [], []
        for j, lig_ct_idx, rec_ct_idx, lig_gene_idx, rec_gene_idx, mode in col_meta:
            # 1) 取原来向量化的 count_lig, count_rec
            lig_local = lig_ct2local[lig_ct_idx][lig_gene_idx]
            rec_local = rec_ct2local[rec_ct_idx][rec_gene_idx]
            count_lig = np.array(patch_by_lig[lig_ct_idx][:, lig_local].toarray()).ravel()
            count_rec = np.array(patch_by_rec[rec_ct_idx][:, rec_local].toarray()).ravel()

            # 2) 只要是 autocrine，都使用相同的 co-expression 向量
            if lig_ct_idx == rec_ct_idx:
                # coexpr 为同一细胞既表达 L 又表达 R，长度 = n_cells
                coexpr = (
                    expr_bool[:, lig_gene_idx].multiply(expr_bool[:, rec_gene_idx])
                    .toarray().ravel().astype(int)
                    *
                    (cell_types_idx == lig_ct_idx).astype(int)
                )
                auto = patch_by_cell.dot(coexpr)  # vectorized: shape=(n_patches,)
            else:
                auto = np.zeros(n_patches, int)

            # 3) 根据 mode 填值到对应列
            if mode == "autocrine":
                rows = np.nonzero(auto)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(auto[rows].tolist())
            else:  # paracrine
                if lig_ct_idx == rec_ct_idx:
                    total = count_lig * count_rec
                    para = total - auto
                else:
                    para = count_lig * count_rec

                rows = np.nonzero(para)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(para[rows].tolist())

        # ─── 8) 组装成最终稀疏矩阵 ───────────────────────────────────────────────────
        patch_lri_matrix = csr_matrix(
            (data_vals, (row_inds, col_inds)),
            shape=(n_patches, n_columns),
            dtype=int
        )
        self.patch_lri_matrix = patch_lri_matrix
        return patch_lri_matrix

    def create_metadata_dataframes(self, adata: anndata.AnnData) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Create metadata dataframes for patch-tma correspondence and cell-patch correspondence.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics data
            
        Returns
        -------
        patch_tma_df : pd.DataFrame
            DataFrame with patch_id and tma_id correspondence
        cell_patch_df : pd.DataFrame
            DataFrame with cell_id and patch_id correspondence
        """
        print("Creating metadata dataframes...")
        
        # Create cell-patch correspondence
        cell_patch_df = pd.DataFrame({
            'cell_id': adata.obs['cell_id'],
            'patch_id': self.patch_assignments,
            'tma_id': adata.obs['tma_id'],
            'cell_type': adata.obs['cell_type']
        })
        
        # Create patch-tma correspondence
        # For each patch, find the most common tma_id
        patch_tma_data = []
        unique_patches = np.unique(self.patch_assignments)
        unique_patches = [p for p in unique_patches if p in self.patch_coords]
        
        for patch_id in unique_patches:
            patch_cells = cell_patch_df[cell_patch_df['patch_id'] == patch_id]
            
            if len(patch_cells) > 0:
                # Get most common tma_id in this patch
                tma_counts = patch_cells['tma_id'].value_counts()
                most_common_tma = tma_counts.index[0]
                
                # Get patch center coordinates
                center_x, center_y = self.patch_coords[patch_id]
                
                patch_tma_data.append({
                    'patch_id': patch_id,
                    'tma_id': most_common_tma,
                    'center_x': center_x,
                    'center_y': center_y,
                    'n_cells': len(patch_cells),
                    'n_cell_types': patch_cells['cell_type'].nunique()
                })
        
        patch_tma_df = pd.DataFrame(patch_tma_data)
        
        return patch_tma_df, cell_patch_df
    
    def run_analysis(self, adata: anndata.AnnData, output_dir: str, 
                     use_batch_processing: bool = True, batch_size: int = 500) -> Dict:
        """
        Run the complete patch-based LRI analysis.
        
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
        print("Starting patch-based LRI analysis...")
        print(f"Patch size: {self.patch_size} μm")
        print(f"Data shape: {adata.shape}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Step 1: Create spatial patches
        patch_assignments, patch_info = self.create_spatial_patches(adata)
        
        # Step 2: Prepare LRI database
        lr_pairs = self.prepare_lri_database(adata)
        
        # Step 3: Create column structure
        column_names = self.create_column_structure(adata)
        
        # Step 4: Build patch-LRI matrix
        patch_lri_matrix = self.build_patch_lri_matrix_with_mode(adata, use_batch_processing, batch_size)
        
        # Step 5: Create metadata dataframes
        patch_tma_df, cell_patch_df = self.create_metadata_dataframes(adata)
        
        # Step 6: Save results
        print("Saving results...")
        
        # Save sparse matrix
        matrix_file = os.path.join(output_dir, 'patch_lri_matrix.npz')
        sparse.save_npz(matrix_file, patch_lri_matrix)
        
        # Save column names
        columns_file = os.path.join(output_dir, 'patch_lri_columns.csv')
        pd.DataFrame({'column_name': self.column_names}).to_csv(columns_file, index=False)
        
        # Save metadata dataframes
        patch_tma_file = os.path.join(output_dir, 'patch_tma_correspondence.csv')
        cell_patch_file = os.path.join(output_dir, 'cell_patch_correspondence.csv')
        
        patch_tma_df.to_csv(patch_tma_file, index=False)
        cell_patch_df.to_csv(cell_patch_file, index=False)
        
        # Save analysis parameters
        params_file = os.path.join(output_dir, 'analysis_parameters.csv')
        params_df = pd.DataFrame({
            'parameter': ['patch_size', 'resource_name', 'n_patches', 'n_lri_combinations', 'matrix_sparsity'],
            'value': [
                self.patch_size,
                self.resource_name,
                patch_info['n_patches'],
                len(column_names),
                f"{(1 - patch_lri_matrix.nnz / np.prod(patch_lri_matrix.shape)) * 100:.2f}%"
            ]
        })
        params_df.to_csv(params_file, index=False)
        
        print(f"Results saved to: {output_dir}")
        print(f"- Patch-LRI matrix: {matrix_file}")
        print(f"- Column names: {columns_file}")
        print(f"- Patch-TMA correspondence: {patch_tma_file}")
        print(f"- Cell-patch correspondence: {cell_patch_file}")
        print(f"- Analysis parameters: {params_file}")
        
        return {
            'patch_lri_matrix': patch_lri_matrix,
            'column_names': column_names,
            'patch_tma_df': patch_tma_df,
            'cell_patch_df': cell_patch_df,
            'patch_info': patch_info,
            'lr_pairs': lr_pairs
        }


def load_patch_lri_results(output_dir: str, sparse_matrix_name: str = 'patch_lri_matrix.npz') -> Dict:
    """
    Load previously saved patch-LRI analysis results.
    
    Parameters
    ----------
    output_dir : str
        Directory containing saved results
        
    Returns
    -------
    results : dict
        Dictionary containing loaded results
    """
    print(f"Loading patch-LRI results from: {output_dir}")
    
    # Load sparse matrix
    matrix_file = os.path.join(output_dir, sparse_matrix_name)
    patch_lri_matrix = sparse.load_npz(matrix_file)
    
    # Load column names
    columns_file = os.path.join(output_dir, 'patch_lri_columns.csv')
    column_names = pd.read_csv(columns_file)['column_name'].tolist()
    
    # Load metadata dataframes
    patch_tma_file = os.path.join(output_dir, 'patch_tma_correspondence.csv')
    cell_patch_file = os.path.join(output_dir, 'cell_patch_correspondence.csv')
    
    patch_tma_df = pd.read_csv(patch_tma_file)
    cell_patch_df = pd.read_csv(cell_patch_file)
    
    # Load parameters
    params_file = os.path.join(output_dir, 'analysis_parameters.csv')
    params_df = pd.read_csv(params_file)
    
    print(f"Loaded matrix shape: {patch_lri_matrix.shape}")
    print(f"Matrix sparsity: {params_df[params_df['parameter'] == 'matrix_sparsity']['value'].iloc[0]}")
    
    return {
        'patch_lri_matrix': patch_lri_matrix,
        'column_names': column_names,
        'patch_tma_df': patch_tma_df,
        'cell_patch_df': cell_patch_df,
        'parameters': params_df
    }


def permute_spatial_obsm(adata, seed=42):
    rng = np.random.RandomState(seed)
    coords = adata.obsm['spatial'][:, :2].copy()
    perm = rng.permutation(len(coords))

    adata.obsm['spatial'][:, :2] = coords[perm]

    return adata

def permute_expression(adata, seed=42):
    """
    Permute the cell×gene expression matrix by shuffling each column independently.
    """
    rng = np.random.RandomState(seed)
    X = adata.X
    # 保证稀疏矩阵格式
    if not sp.issparse(X):
        X = sparse.csr_matrix(X)
    # 转成 CSC 方便按列切片
    X_csc = X.tocsc()
    n_cells, n_genes = X_csc.shape

    rows, cols, data = [], [], []
    # 对每个基因（列）单独打乱
    for j in range(n_genes):
        col = X_csc[:, j]                  # 稀疏列
        dense_col = col.toarray().ravel()  # 转成一维密集向量
        rng.shuffle(dense_col)             # 随机打乱
        nz = np.nonzero(dense_col)[0]      # 找到非零位置
        rows.append(nz)                    
        cols.append(np.full_like(nz, j))   # 列索引全部为 j
        data.append(dense_col[nz])         # 对应的打乱后值

    # 合并所有列的结果，重建稀疏矩阵
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)
    adata.X = sparse.csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes))
    return adata

def permute_expression_optimized(adata, seed=42):
    """优化版本的表达式置换"""
    rng = np.random.RandomState(seed)
    X = adata.X
    if not sp.issparse(X):
        X = sparse.csr_matrix(X)
    
    # 转换为 COO 格式进行高效操作
    X_coo = X.tocoo()
    n_cells, n_genes = X_coo.shape
    
    # 为每个基因生成置换索引
    new_rows = []
    for j in range(n_genes):
        # 找到该基因的所有非零位置
        gene_mask = X_coo.col == j
        gene_rows = X_coo.row[gene_mask]
        
        # 生成置换
        perm_indices = rng.permutation(n_cells)
        # 将原行索引映射到新位置
        new_gene_rows = perm_indices[gene_rows]
        new_rows.append(new_gene_rows)
    
    # 重建稀疏矩阵
    all_new_rows = np.concatenate(new_rows)
    adata.X = sparse.csr_matrix(
        (X_coo.data, (all_new_rows, X_coo.col)), 
        shape=(n_cells, n_genes)
    )
    return adata