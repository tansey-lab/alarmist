#!/usr/bin/env python3
"""
Neighborhood-based LRI Analysis (adapter over PatchLRIAnalyzer)

Reuses your ALARMIST PatchLRIAnalyzer to compute a Neighborhood x LRI matrix,
given an external cell→neighborhood map.

Changes vs your original:
- We *do not* create spatial patches. Instead, we set `self.patch_assignments`
  to your neighborhood IDs (one integer per cell).
- We override build_* to subtract self-pairs for juxtacrine (Cell-Cell Contact)
  when ligand_ct == receptor_ct.  (# CHANGED: marked below)
- We save CSV (not TSV) and store the big matrix as .npz.

Author: (kept minimal & faithful to ALARMIST)
"""
import os
import numpy as np
import pandas as pd
from scipy import sparse as sp
from typing import Dict, List, Tuple, Optional

import anndata as ad
from patch_lri_analysis import PatchLRIAnalyzer

class SCNeighborhoodLRIAnalyzer(PatchLRIAnalyzer):
    """
    Adapter that reuses PatchLRIAnalyzer to work with arbitrary neighborhoods.
    Updated to use adata.obs.index consistently.
    """

    def use_precomputed_neighborhoods(
        self,
        adata: ad.AnnData,
        membership: pd.DataFrame,
        cell_col: str = "cell",
        nb_col: str = "neighborhood_id"
    ):
        """
        Set self.patch_assignments to neighborhood ids (one per cell).
        
        Parameters
        ----------
        adata : AnnData object
        membership : DataFrame with columns [cell_col, nb_col]
            - cell_col should match adata.obs.index
        cell_col : str
            Column name in membership for cell IDs
        nb_col : str
            Column name in membership for neighborhood IDs
        """
        # Always use obs.index as the reference
        adata_ids = pd.Index(adata.obs.index.astype(str).values, name="cell")
        
        # Normalize membership cell IDs to string
        mem = membership.copy()
        mem[cell_col] = mem[cell_col].astype(str)
        
        # Create alignment in adata order
        df = pd.DataFrame({"cell": adata_ids})
        df = df.merge(
            mem[[cell_col, nb_col]].rename(columns={cell_col: "cell"}), 
            on="cell", 
            how="left"
        )
        
        # Check for missing assignments
        n_missing = df[nb_col].isna().sum()
        if n_missing > 0:
            missing_cells = df[df[nb_col].isna()]["cell"].head(10).tolist()
            print(f"WARNING: {n_missing} cells have no neighborhood assignment")
            print(f"  Examples: {missing_cells}")
            
            # Option 1: Assign to singleton neighborhoods
            max_nb_id = df[nb_col].max()
            if pd.isna(max_nb_id):
                max_nb_id = -1
            for idx in df[df[nb_col].isna()].index:
                max_nb_id += 1
                df.loc[idx, nb_col] = max_nb_id
            print(f"  Assigned unmatched cells to singleton neighborhoods")
            
            # Option 2: Raise error (uncomment if preferred)
            # raise ValueError(f"{n_missing} cells missing neighborhood_id")
        
        # Set patch assignments (neighborhoods)
        self.patch_assignments = df[nb_col].astype(int).values
        
        # Create minimal patch coords (not used for neighborhoods but required by parent)
        unique_nbs = np.unique(self.patch_assignments)
        self.patch_coords = {int(pid): (np.nan, np.nan) for pid in unique_nbs}
        
        print(f"Loaded {len(unique_nbs)} neighborhoods for {len(adata)} cells")
        print(f"  Neighborhood sizes: min={np.bincount(self.patch_assignments).min()}, "
              f"max={np.bincount(self.patch_assignments).max()}")

    def build_neighborhood_lri_matrix_with_mode(
        self,
        adata: ad.AnnData,
        signaling_types: List[str],
        use_batch_processing: bool = True,
        batch_size: int = 500
    ) -> sp.csr_matrix:
        """
        Build neighborhood×LRI matrix with proper handling of juxtacrine interactions.
        For juxtacrine: exclude self-pairs when ligand_ct == receptor_ct.
        """
        import scipy.sparse as sparse
        from scipy.sparse import coo_matrix, csr_matrix

        print("Building NEIGHBORHOOD-LRI matrix with autocrine/paracrine distinction...")

        # Prepare basics
        unique_patches = np.array([p for p in np.unique(self.patch_assignments) 
                                  if p in self.patch_coords])
        patch_idx_map = {int(pid): i for i, pid in enumerate(unique_patches)}
        n_patches = len(unique_patches)
        n_columns = len(self.column_names)
        print(f"Processing {n_patches} neighborhoods × {n_columns} LRI combinations")

        # Index mappings
        ct_to_idx = {ct: i for i, ct in enumerate(self.cell_types)}
        cell_types_idx = np.array([ct_to_idx[ct] for ct in adata.obs[self.cell_type_column]], 
                                  dtype=int)
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

        # Parse column metadata
        col_meta = []
        for idx, (lig, rec_str) in enumerate(self.lr_pairs):
            sig_type = signaling_types[idx]
            for lig_ct in self.cell_types:
                for rec_ct in self.cell_types:
                    if sig_type == 'Cell-Cell Contact':
                        col_meta.append((
                            len(col_meta),
                            ct_to_idx[lig_ct],
                            ct_to_idx[rec_ct],
                            gene_to_idx[lig],
                            rec_str,
                            'juxtacrine',
                            sig_type
                        ))
                    else:
                        if lig_ct == rec_ct:
                            col_meta.append((len(col_meta), ct_to_idx[lig_ct], 
                                           ct_to_idx[rec_ct], gene_to_idx[lig], 
                                           rec_str, 'autocrine', sig_type))
                            col_meta.append((len(col_meta), ct_to_idx[lig_ct], 
                                           ct_to_idx[rec_ct], gene_to_idx[lig], 
                                           rec_str, 'paracrine', sig_type))
                        else:
                            col_meta.append((len(col_meta), ct_to_idx[lig_ct], 
                                           ct_to_idx[rec_ct], gene_to_idx[lig], 
                                           rec_str, 'paracrine', sig_type))

        # Binarize expression
        X = adata.X
        if sp.issparse(X):
            expr_bool = (X > 0).astype(int).tocsc()
        else:
            expr_bool = sp.csr_matrix((X > 0).astype(int)).tocsc()
        expr_coo = expr_bool.tocoo()

        # Build patch_by_lig matrices
        patch_by_lig = {}
        print("Building neighborhood-by-ligand matrices...")
        for ct_idx in range(len(self.cell_types)):
            mask_cells = (cell_types_idx == ct_idx)
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix((expr_coo.data[entry_mask],
                             (expr_coo.row[entry_mask], expr_coo.col[entry_mask])),
                            shape=expr_bool.shape)

            lig_genes_ct = sorted({lig for (_, lct, _, lig, _, _, _) in col_meta 
                                  if lct == ct_idx})
            lig_to_local = {g: i for i, g in enumerate(lig_genes_ct)}
            lig_mask = np.isin(sub.col, lig_genes_ct)
            rows_cells = sub.row[lig_mask]
            cols_genes = sub.col[lig_mask]
            data_vals_ct = sub.data[lig_mask]
            
            patch_rows = np.array([patch_idx_map[int(self.patch_assignments[c])] 
                                  for c in rows_cells], dtype=int)
            local_cols = np.array([lig_to_local[g] for g in cols_genes], dtype=int)

            lig_coo = coo_matrix((data_vals_ct, (patch_rows, local_cols)),
                                shape=(n_patches, len(lig_genes_ct)))
            lig_coo.sum_duplicates()
            patch_by_lig[ct_idx] = lig_coo.tocsr()

        # Build patch_by_rec matrices
        patch_by_rec = {}
        print("Building neighborhood-by-receptor matrices...")
        for ct_idx in range(len(self.cell_types)):
            mask_cells = (cell_types_idx == ct_idx)
            entry_mask = mask_cells[expr_coo.row]
            sub = coo_matrix((expr_coo.data[entry_mask],
                             (expr_coo.row[entry_mask], expr_coo.col[entry_mask])),
                            shape=expr_bool.shape)

            all_rec_genes_ct = set()
            for (_, _, rct, _, rec_str, _, _) in col_meta:
                if rct == ct_idx:
                    all_rec_genes_ct.update(rec_str.split('_'))
            all_rec_genes_ct = sorted(all_rec_genes_ct)

            rec_gene_matrices = {}
            for rec_gene in all_rec_genes_ct:
                rec_gene_idx = gene_to_idx[rec_gene]
                rec_mask = (sub.col == rec_gene_idx)
                rows_cells = sub.row[rec_mask]
                data_vals_ct = sub.data[rec_mask]
                patch_rows = np.array([patch_idx_map[int(self.patch_assignments[c])] 
                                      for c in rows_cells], dtype=int)

                rec_coo = coo_matrix((data_vals_ct, (patch_rows, np.zeros(len(patch_rows), dtype=int))),
                                    shape=(n_patches, 1))
                rec_coo.sum_duplicates()
                rec_gene_matrices[rec_gene] = rec_coo.tocsr()
            patch_by_rec[ct_idx] = rec_gene_matrices

        lig_ct2local = {
            ct: {g: i for i, g in enumerate(sorted({lig for (_, lct, _, lig, _, _, _) in col_meta 
                                                    if lct == ct}))}
            for ct in range(len(self.cell_types))
        }

        # Build neighborhood-by-cell matrix
        print("Building neighborhood-by-cell mapping...")
        n_cells = adata.n_obs
        cell_patches = np.array(self.patch_assignments, dtype=int)
        pb_rows = [patch_idx_map[int(p)] for p in cell_patches]
        pb_cols = list(range(n_cells))

        patch_by_cell = sp.coo_matrix((np.ones(n_cells, int), (pb_rows, pb_cols)),
                                      shape=(n_patches, n_cells)).tocsr()

        # Compute interactions
        print("Computing LRI interactions...")
        row_inds, col_inds, data_vals = [], [], []

        for j, lig_ct_idx, rec_ct_idx, lig_gene_idx, rec_str, mode, sig_type in col_meta:
            if j % 500 == 0:
                print(f"  Progress: {j}/{n_columns}")

            # Ligand counts
            lig_local = lig_ct2local[lig_ct_idx][lig_gene_idx]
            count_lig = np.array(patch_by_lig[lig_ct_idx][:, lig_local].toarray()).ravel()

            # Receptor counts (AND across subunits)
            rec_genes = rec_str.split('_')
            if len(rec_genes) == 1:
                count_rec = np.array(patch_by_rec[rec_ct_idx][rec_genes[0]].toarray()).ravel()
            else:
                count_rec = np.array(patch_by_rec[rec_ct_idx][rec_genes[0]].toarray()).ravel()
                for rec_gene in rec_genes[1:]:
                    count_rec = np.minimum(count_rec, 
                                          np.array(patch_by_rec[rec_ct_idx][rec_gene].toarray()).ravel())

            # Autocrine same-cell overlaps
            if lig_ct_idx == rec_ct_idx:
                coexpr = expr_bool[:, lig_gene_idx].toarray().ravel().astype(int)
                for rec_gene in rec_genes:
                    rec_gene_idx = gene_to_idx[rec_gene]
                    coexpr = coexpr * expr_bool[:, rec_gene_idx].toarray().ravel().astype(int)
                coexpr = coexpr * (cell_types_idx == lig_ct_idx).astype(int)
                auto = np.array(patch_by_cell.dot(coexpr)).ravel()
            else:
                auto = np.zeros(n_patches, dtype=int)

            if mode == "juxtacrine":
                total = count_lig * count_rec
                # Exclude self-pairs when ligand_ct == receptor_ct
                if lig_ct_idx == rec_ct_idx:
                    total = total - auto
                    total[total < 0] = 0
                rows = np.nonzero(total)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(total[rows].tolist())

            elif mode == "autocrine":
                rows = np.nonzero(auto)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(auto[rows].tolist())

            else:  # paracrine
                if lig_ct_idx == rec_ct_idx:
                    para = count_lig * count_rec - auto
                else:
                    para = count_lig * count_rec
                rows = np.nonzero(para)[0]
                row_inds.extend(rows.tolist())
                col_inds.extend([j] * len(rows))
                data_vals.extend(para[rows].tolist())

        # Assemble matrix
        mat = sp.csr_matrix((data_vals, (row_inds, col_inds)), 
                           shape=(n_patches, n_columns), dtype=int)
        self.patch_lri_matrix = mat
        print(f"Matrix shape: {mat.shape}")
        print(f"Matrix density: {mat.nnz / (n_patches * n_columns) * 100:.2f}%")
        print(f"Non-zero entries: {mat.nnz:,}")
        return mat

    def run_neighborhood_analysis(
        self,
        adata: ad.AnnData,
        output_dir: str,
        use_batch_processing: bool = True,
        batch_size: int = 500
    ) -> Dict:
        """
        Run the full pipeline using precomputed neighborhoods.
        Always uses adata.obs.index for cell identification.
        """
        os.makedirs(output_dir, exist_ok=True)

        # Prepare LR database
        lr_pairs, receptor_genes_list, signaling_types = self.prepare_lri_database(adata)

        # Create columns
        column_names = self.create_column_structure(adata, signaling_types, self.cell_type_column)

        # Build matrix
        nb_lri_matrix = self.build_neighborhood_lri_matrix_with_mode(
            adata, signaling_types, use_batch_processing, batch_size
        )

        # Create cell-neighborhood correspondence using obs.index
        cell_ids = adata.obs.index.astype(str).values
        cell_nb_df = pd.DataFrame({
            "cell": cell_ids,
            "patch_id": self.patch_assignments,
            "cell_type": adata.obs[self.cell_type_column].astype(str).values,
            "tma_id": adata.obs['tma_id'].values,
        })

        # Save outputs
        matrix_file = os.path.join(output_dir, "neighborhood_lri_matrix.npz")
        sp.save_npz(matrix_file, nb_lri_matrix)
        print(f"Saved matrix: {matrix_file}")

        columns_file = os.path.join(output_dir, "neighborhood_lri_columns.csv")
        pd.DataFrame({"column_name": column_names}).to_csv(columns_file, index=False)
        print(f"Saved columns: {columns_file}")

        cell_nb_file = os.path.join(output_dir, "cell_patch_correspondence.csv")
        cell_nb_df.to_csv(cell_nb_file, index=False)
        print(f"Saved correspondence: {cell_nb_file}")

        # Save parameters
        params_file = os.path.join(output_dir, "analysis_parameters.csv")
        sparsity = 1 - (nb_lri_matrix.nnz / np.prod(nb_lri_matrix.shape))
        unique_nbs = len(np.unique(self.patch_assignments))
        pd.DataFrame({
            "parameter": [
                "n_neighborhoods", 
                "n_lri_columns", 
                "matrix_sparsity",
                "matrix_nonzero_entries",
                "n_cells",
                "n_cell_types",
                "n_lr_pairs"
            ],
            "value": [
                unique_nbs, 
                len(column_names), 
                f"{sparsity*100:.2f}%",
                nb_lri_matrix.nnz,
                adata.n_obs,
                len(self.cell_types),
                len(lr_pairs)
            ]
        }).to_csv(params_file, index=False)
        print(f"Saved parameters: {params_file}")

        return {
            "neighborhood_lri_matrix": nb_lri_matrix,
            "column_names": column_names,
            "cell_neighborhood_df": cell_nb_df,
        }