import os
import numpy as np
import pandas as pd

import scipy.sparse as sp
import sparse
from pathlib import Path
import matplotlib.pyplot as plt
from bptf import BPTF, save_bptf, load_bptf

import sys
import os


# Add scripts directory to path
scripts_path = os.path.join(os.path.dirname(os.getcwd()), "scripts")
if scripts_path not in sys.path:
    sys.path.insert(0, scripts_path)

from neighborhood_lri_analysis import load_cell_lri_results

INPUT_DIR = 'results/ES_cellphone_15/single_cell'
OUTPUT_DIR  = 'results/ES_cellphone_15/single_cell/bptf'
N_COMPONENTS= 15
MAX_ITER    = 200

random_state = 0

np.random.seed(random_state)

results = load_cell_lri_results(INPUT_DIR)

mat = results['cell_lri_matrix']
cols = results['column_names']
metadata_df = results['cell_metadata_df']

if isinstance(mat, sp.spmatrix):
    cell_lri_coo = sparse.COO(mat)


model = load_bptf('results/ES_cellphone_15/bptf/bptf_model.npz/bptf_1.npz')

G_DK_M = model.G_DK_M
patch_loadings = G_DK_M[0]         # shape = (num_patches, K)
lri_factors    = G_DK_M[1].T       # Turn into (K, num_lris)
print(f"Patch loadings shape: {patch_loadings.shape}")
print(f"LRI factors shape: {lri_factors.shape}")

n_cells, n_lri = cell_lri_coo.shape
K = model.n_components
alpha = model.alpha
print(f"Number of cells: {n_cells}, Number of LRIs: {n_lri}, Number of components: {K}, Alpha: {alpha}")

# 1) 新建一个与 cell×LRI 尺寸匹配的 BPTF
proj = BPTF(data_shape=(n_cells, n_lri), n_components=K, alpha=alpha)
proj._init(modes=[0, 1]) 

# 2) 把列 mode=1 的变分参数/期望 从 patch 模型拷贝过来（保持完全一致）
proj.shp_DK_M[1][:, :] = model.shp_DK_M[1]
proj.rte_DK_M[1][:, :] = model.rte_DK_M[1]
proj._update_cache(1)  # 刷新 E/G/_sumE/beta
proj.beta_M[1] = model.beta_M[1]

# 3) 固定列 mode（不更新）
proj._clamp_component(1)

# 4) 只更新行 mode=0
proj._update(cell_lri_coo, modes=[0], max_iter=MAX_ITER, verbose=True)

U_cell = proj.G_DK_M[0]
V_fixed = proj.G_DK_M[1]

# Save model
output_dir = OUTPUT_DIR
model_path = Path(output_dir) / 'bptf_model.npz'
# save_bptf(model, model_path)
print(f"BPTF model saved to: {model_path}")

# Save factor matrices as numpy arrays
np.save(os.path.join(output_dir, 'cell_loadings.npy'), U_cell)
np.save(os.path.join(output_dir, 'lri_factors.npy'), lri_factors)