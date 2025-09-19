#!/usr/bin/env python3
"""
Patch-LRI-BPTF Analysis Pipeline Runner

Simple script to run the complete pipeline with default parameters.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def run_command_with_conda(cmd, description, conda_env=None):
    """Run a command with optional conda environment activation"""
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    
    if conda_env:
        # Create the full command with conda activation
        full_cmd = f"conda run -n {conda_env} {cmd}"
        print(f"CONDA ENV: {conda_env}")
    else:
        full_cmd = cmd
        
    print(f"COMMAND: {full_cmd}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(full_cmd, shell=True, check=True, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {e}")
        if e.stdout:
            print(f"STDOUT: {e.stdout}")
        if e.stderr:
            print(f"STDERR: {e.stderr}")
        return False


def run_command(cmd, description):
    """Run a command without conda (for backward compatibility)"""
    return run_command_with_conda(cmd, description, conda_env=None)


def main():
    parser = argparse.ArgumentParser(description='Run Patch-LRI-BPTF pipeline')
    parser.add_argument('--data-file', required=True, help='Input h5ad file')
    parser.add_argument('--output-dir', default='results', help='Output directory')
    parser.add_argument('--n-components', type=int, default=15, help='BPTF components')
    parser.add_argument('--patch-size', type=float, default=50.0, help='Patch size in μm')
    parser.add_argument('--skip-glm', action='store_true', help='Skip GLM analysis')
    parser.add_argument('--skip-visualization', action='store_true', help='Skip GLM visualization')
    parser.add_argument('--bptf-conda-env', default='bptf', help='Conda environment for BPTF step (default: bptf)')
    parser.add_argument('--main-conda-env', default='tf2.10', help='Conda environment for other steps (default: tf2.10)')
    parser.add_argument('--no-conda', action='store_true', help='Run without conda environment switching')
    
    args = parser.parse_args()
    
    # Setup paths
    output_dir = Path(args.output_dir)
    patch_lri_dir = output_dir / "patch_lri"
    bptf_dir = output_dir / "bptf"
    glm_dir = output_dir / "glm_results"
    plots_dir = output_dir / "plots"
    
    # Create output directories
    output_dir.mkdir(exist_ok=True)
    
    script_dir = Path(__file__).parent / "scripts"
    
    print("="*60)
    print("PATCH-LRI-BPTF ANALYSIS PIPELINE")
    print("="*60)
    print(f"Input data: {args.data_file}")
    print(f"Output directory: {output_dir}")
    print(f"Patch size: {args.patch_size} μm")
    print(f"BPTF components: {args.n_components}")
    
    if not args.no_conda:
        print(f"BPTF conda env: {args.bptf_conda_env}")
        print(f"Main conda env: {args.main_conda_env}")
    else:
        print("Running without conda environment switching")
    
    # Determine conda environments to use
    main_env = None if args.no_conda else args.main_conda_env
    bptf_env = None if args.no_conda else args.bptf_conda_env
    
    # Step 1: Patch-LRI Analysis (use main environment)
    cmd1 = (f"python {script_dir}/01_run_patch_lri_analysis.py "
            f"--data-file {args.data_file} "
            f"--output-dir {patch_lri_dir} "
            f"--patch-size {args.patch_size}")
    
    if not run_command_with_conda(cmd1, "Step 1: Patch-LRI Analysis", main_env):
        sys.exit(1)
    
    # Step 2: BPTF Matrix Factorization (use BPTF environment)
    cmd2 = (f"python {script_dir}/02_bptf_matrix_factorization.py "
            f"--input-dir {patch_lri_dir} "
            f"--output-dir {bptf_dir} "
            f"--n-components {args.n_components}")
    
    if not run_command_with_conda(cmd2, "Step 2: BPTF Matrix Factorization", bptf_env):
        sys.exit(1)
    
    # Step 3: BPTF Visualization (use main environment)
    cmd3 = (f"python {script_dir}/03_bptf_visualization.py "
            f"--bptf-dir {bptf_dir} "
            f"--patch-dir {patch_lri_dir} "
            f"--data-file {args.data_file} "
            f"--output-dir {plots_dir}/bptf_plots")
    
    if not run_command_with_conda(cmd3, "Step 3: BPTF Visualization", main_env):
        print("Warning: BPTF visualization failed, but continuing...")
    
    # Step 4: GLM Analysis (use main environment)
    if not args.skip_glm:
        cmd4 = (f"python {script_dir}/04_poisson_glm.py "
                f"--data-file {args.data_file} "
                f"--results-dir {bptf_dir} "
                f"--patch-lri-dir {patch_lri_dir} "
                f"--output-dir {glm_dir}")
        
        if not run_command_with_conda(cmd4, "Step 4: GLM Analysis", main_env):
            print("Warning: GLM analysis failed, but continuing...")
    
    # Step 5: GLM Results Visualization (use main environment)
    if not args.skip_visualization and not args.skip_glm:
        cmd5 = (f"python {script_dir}/05_glm_results.py "
                f"--data-file {args.data_file} "
                f"--results-dir {glm_dir} "
                f"--output-dir {plots_dir}/glm_plots")
        
        if not run_command_with_conda(cmd5, "Step 5: GLM Results Visualization", main_env):
            print("Warning: GLM visualization failed, but continuing...")
    
    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Results saved in: {output_dir}")
    print("\nNext steps:")
    print("- Check results in the output directory")
    print("- Use notebooks/bptf_analysis.ipynb for detailed analysis")
    if args.skip_glm:
        print("- Run GLM analysis: python scripts/04_poisson_glm.py")
    if args.skip_visualization:
        print("- Generate GLM plots: python scripts/05_glm_results.py")


if __name__ == '__main__':
    main()