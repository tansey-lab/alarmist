#!/usr/bin/env python3
"""
Create CSV files for Kenny's analysis

Generates two CSV files:
1. LRI factor loadings: rows = LRI tuples (sender, receiver, ligand, receptor), columns = motifs
2. Differential expression results: rows = (cell type, gene) pairs, columns = motifs (coeffs + p-values)

Usage:
    python scripts/create_kenny_csvs.py --bptf-dir results/bptf --glm-dir results/glm --output-dir results/kenny_csvs
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def load_bptf_results(bptf_dir):
    """Load BPTF results from directory"""
    print(f"Loading BPTF results from {bptf_dir}...")
    
    # Load LRI factors (motifs x LRIs)
    lri_factors = np.load(os.path.join(bptf_dir, 'lri_factors.npy'))
    
    # Load LRI motifs CSV for column names
    lri_motifs_df = pd.read_csv(os.path.join(bptf_dir, 'lri_motifs.csv'))
    
    print(f"Loaded LRI factors: {lri_factors.shape}")
    print(f"Loaded LRI motifs data: {len(lri_motifs_df)} entries")
    
    return lri_factors, lri_motifs_df


def create_lri_factor_csv(lri_factors, lri_motifs_df, output_path, splitter='|'):
    """
    Create CSV with LRI tuples as rows and motifs as columns
    
    Format:
    - Column A: Sender cell type
    - Column B: Receiver cell type  
    - Column C: Ligand gene
    - Column D: Receptor gene
    - Column E: Signaling type (autocrine/paracrine)
    - Columns F+: Motif factor loadings (one column per motif)
    """
    print("Creating LRI factor loadings CSV...")
    
    # Get unique LRI names
    unique_lris = lri_motifs_df['lri_name'].unique()
    n_motifs = lri_factors.shape[0]
    
    # Parse LRI names into components
    lri_data = []
    for lri_name in unique_lris:
        parts = lri_name.split(splitter)
        if len(parts) >= 4:
            sender = parts[0]
            receiver = parts[1] 
            ligand = parts[2]
            receptor = parts[3]
            signaling_type = parts[4] if len(parts) > 4 else 'unknown'
            
            # Get factor loadings for this LRI across all motifs
            lri_idx = lri_motifs_df[lri_motifs_df['lri_name'] == lri_name]['lri_idx'].iloc[0]
            
            # Create row with LRI info + factor loadings
            row = {
                'Sender': sender,
                'Receiver': receiver,
                'Ligand': ligand,
                'Receptor': receptor,
                'Signaling_Type': signaling_type
            }
            
            # Add factor loadings for each motif
            for motif_idx in range(n_motifs):
                row[f'Motif_{motif_idx}'] = lri_factors[motif_idx, lri_idx]
                
            lri_data.append(row)
    
    # Create DataFrame
    lri_df = pd.DataFrame(lri_data)
    
    # Save CSV
    lri_df.to_csv(output_path, index=False)
    print(f"LRI factor loadings CSV saved: {output_path}")
    print(f"Shape: {lri_df.shape}")
    
    return lri_df


def load_glm_results(glm_dir):
    """Load all GLM differential expression results"""
    print(f"Loading GLM results from {glm_dir}...")
    
    # Find all GLM result files
    glm_files = glob.glob(os.path.join(glm_dir, "motif_*_celltype_*_de_results.csv"))
    
    if not glm_files:
        print(f"No GLM results found in {glm_dir}")
        return None
    
    print(f"Found {len(glm_files)} GLM result files")
    
    # Parse file names to get motif and cell type info
    glm_data = {}
    
    for file_path in glm_files:
        filename = os.path.basename(file_path)
        # Parse: motif_{X}_celltype_{Y}_de_results.csv
        parts = filename.replace('.csv', '').split('_')
        
        try:
            motif_idx = int(parts[1])  # motif_X -> X
            celltype_start = parts.index('celltype') + 1
            celltype_parts = parts[celltype_start:-2]  # everything between celltype and 'de results'
            celltype = '_'.join(celltype_parts)
            
            # Load the results
            df = pd.read_csv(file_path)
            
            if motif_idx not in glm_data:
                glm_data[motif_idx] = {}
            glm_data[motif_idx][celltype] = df
            
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse filename {filename}: {e}")
            continue
    
    print(f"Loaded GLM data for {len(glm_data)} motifs")
    return glm_data


def create_de_results_csv_by_celltype(glm_data, output_dir):
    """
    Create separate CSV files for each cell type with genes as rows and motifs as columns
    
    Format for each cell type CSV:
    - Rows: Genes
    - Columns: Motifs with coefficient and p-value pairs
    """
    print("Creating differential expression results CSVs by cell type...")
    
    if glm_data is None:
        print("No GLM data available")
        return None
    
    motif_indices = sorted(glm_data.keys())
    
    # Collect all cell types
    all_celltypes = set()
    for motif_idx in motif_indices:
        all_celltypes.update(glm_data[motif_idx].keys())
    
    print(f"Found {len(all_celltypes)} cell types")
    print(f"Found {len(motif_indices)} motifs")
    
    celltype_dfs = {}
    
    for celltype in all_celltypes:
        print(f"Processing cell type: {celltype}")
        
        # Collect all genes for this cell type across all motifs
        all_genes = set()
        for motif_idx in motif_indices:
            if celltype in glm_data[motif_idx]:
                df = glm_data[motif_idx][celltype]
                all_genes.update(df['gene'])
        
        if not all_genes:
            print(f"No genes found for cell type {celltype}")
            continue
            
        # Create DataFrame with genes as rows
        genes_list = sorted(list(all_genes))
        celltype_df = pd.DataFrame({'Gene': genes_list})
        celltype_df.set_index('Gene', inplace=True)
        
        # Add columns for each motif (coefficient and p-value)
        for motif_idx in motif_indices:
            coeff_col = f'Motif_{motif_idx}_Coefficient'
            pval_col = f'Motif_{motif_idx}_PValue'
            
            # Initialize with NaN
            celltype_df[coeff_col] = np.nan
            celltype_df[pval_col] = np.nan
            
            # Fill in data if available for this motif
            if celltype in glm_data[motif_idx]:
                df = glm_data[motif_idx][celltype]
                
                for _, row in df.iterrows():
                    gene = row['gene']
                    if gene in celltype_df.index:
                        celltype_df.loc[gene, coeff_col] = row.get('logFC', np.nan)
                        celltype_df.loc[gene, pval_col] = row.get('pval', row.get('qval', np.nan))
        
        # Save CSV for this cell type
        celltype_clean = celltype.replace('/', '_').replace(' ', '_')
        csv_path = os.path.join(output_dir, f"de_results_{celltype_clean}.csv")
        celltype_df.to_csv(csv_path)
        
        celltype_dfs[celltype] = celltype_df
        print(f"Saved {celltype} results: {csv_path} (shape: {celltype_df.shape})")
    
    return celltype_dfs


def main():
    """Main function to create Kenny's CSV files"""
    parser = argparse.ArgumentParser(description='Create CSV files for Kenny analysis')
    parser.add_argument('--bptf-dir', required=True,
                       help='Directory containing BPTF results')
    parser.add_argument('--glm-dir', required=True,
                       help='Directory containing GLM results')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for CSV files')
    parser.add_argument('--splitter', default='|',
                       help='Separator used in LRI names')
    
    args = parser.parse_args()
    
    print("="*60)
    print("CREATING KENNY'S CSV FILES")
    print("="*60)
    print(f"BPTF directory: {args.bptf_dir}")
    print(f"GLM directory: {args.glm_dir}")
    print(f"Output directory: {args.output_dir}")
    print("="*60)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check input directories exist
    if not os.path.exists(args.bptf_dir):
        print(f"Error: BPTF directory not found: {args.bptf_dir}")
        return
        
    if not os.path.exists(args.glm_dir):
        print(f"Error: GLM directory not found: {args.glm_dir}")
        return
    
    try:
        # 1. Create LRI factor loadings CSV
        print("\n1. Processing BPTF results...")
        lri_factors, lri_motifs_df = load_bptf_results(args.bptf_dir)
        
        lri_csv_path = os.path.join(args.output_dir, "lri_factor_loadings.csv")
        lri_df = create_lri_factor_csv(lri_factors, lri_motifs_df, lri_csv_path, args.splitter)
        
        # 2. Create differential expression results CSV (separate by cell type)
        print("\n2. Processing GLM results...")
        glm_data = load_glm_results(args.glm_dir)
        
        if glm_data:
            celltype_dfs = create_de_results_csv_by_celltype(glm_data, args.output_dir)
        else:
            print("Warning: No GLM data found - skipping DE results CSV")
            celltype_dfs = None
        
        print("\n" + "="*60)
        print("CSV CREATION COMPLETED SUCCESSFULLY!")
        print("="*60)
        print(f"Files created in: {args.output_dir}")
        print("Files:")
        print(f"  1. lri_factor_loadings.csv")
        if celltype_dfs:
            print(f"  2. de_results_[celltype].csv (separate file for each cell type)")
            for celltype in sorted(celltype_dfs.keys()):
                celltype_clean = celltype.replace('/', '_').replace(' ', '_')
                print(f"     - de_results_{celltype_clean}.csv")
        
        # Print summary stats
        print(f"\nSummary:")
        print(f"  - LRI tuples: {len(lri_df)}")
        print(f"  - Motifs: {lri_factors.shape[0]}")
        if celltype_dfs:
            print(f"  - Cell types: {len(celltype_dfs)}")
            print(f"  - Motifs with DE data: {len(glm_data)}")
            for celltype, df in celltype_dfs.items():
                print(f"    * {celltype}: {df.shape[0]} genes")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()