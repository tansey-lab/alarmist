import os
import gc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import sparse
from pathlib import Path
from bptf import BPTF, load_bptf
from numpy.lib.format import open_memmap
import sys

# Optional: limit BLAS threads (reduces extra memory)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Add scripts directory to path
scripts_path = os.path.join(os.path.dirname(os.getcwd()), "scripts")
if scripts_path not in sys.path:
    sys.path.insert(0, scripts_path)

from neighborhood_lri_analysis import load_cell_lri_results

INPUT_DIR   = 'results/ES_immune/single_cell'
OUTPUT_DIR  = 'results/ES_immune/single_cell/bptf'
MAX_ITER    = 200
random_state = 0
np.random.seed(random_state)

# ---------- Load the patch-level model (for fixed mode=1) ----------
model = load_bptf('results/ES_immune/bptf/bptf_model.npz/bptf_5.npz')
K     = model.n_components
alpha = model.alpha
model_lri_names = pd.read_csv('/Users/jiayifan/tansey_lab/alarmist/results/ES_immune/patch_lri_columns.csv')
model_lri_names = np.array(model_lri_names['column_name'])

# ---------- Load data as CSR ----------
results = load_cell_lri_results(INPUT_DIR)
mat = results['cell_lri_matrix']
lri_names = np.array(results['column_names'])
# subset mat to same lris as the model has
# del results

# ---------- Subset mat to same LRIs as the model ----------
# Find intersection and preserve the model's order
# _, idx_in_mat, idx_in_model = np.intersect1d(
#     lri_names, model_lri_names, return_indices=True
# )

import numpy as np
from collections import defaultdict, deque

# 交集基于这列
a = lri_names
b = model_lri_names

# 哪些 b 的元素出现在 a 中（按 b 的原始顺序与重复都保留）
mask = np.isin(b, a)
idx_in_model = np.flatnonzero(mask)        # b 中匹配到的行号（保留重复与顺序）

# 为 a 中每个取值建一个“待分配索引队列”，用于一一配对（避免把同一个 a 位置配多次）
pos = defaultdict(deque)
for i, v in enumerate(a):
    pos[v].append(i)

# 将 b[mask] 中每个值，配到 a 中该值的“下一个尚未用过”的位置
idx_in_mat = np.array([pos[v].popleft() if pos[v] else -1 for v in b[mask]], dtype=int)
# 如果 b 的某个值在 a 中出现次数不足，会得到 -1，你也可以改成 np.nan

mat = mat[:, idx_in_mat]      # reorder / subset columns
lri_names = lri_names[idx_in_mat]
print(mat.shape[1], len(model_lri_names))
assert mat.shape[1] == len(model_lri_names), "Mismatch between LRI sets!"

# ---------- Convert matrix to efficient CSR format ----------
if not sp.isspmatrix_csr(mat):
    mat = mat.tocsr(copy=False)

mat.data = mat.data.astype(np.float64, copy=False)
mat.indices = mat.indices.astype(np.int32, copy=False)
mat.indptr = mat.indptr.astype(np.int64, copy=False)

n_cells, n_lri = mat.shape
print(f"Cell-LRI matrix: {n_cells:,} cells × {n_lri:,} LRIs")
print(f"BPTF model: K={K}, alpha={alpha}")

# Fixed column parameters (mode=1)
shp1 = model.shp_DK_M[1]
rte1 = model.rte_DK_M[1]
bet1 = model.beta_M[1]

# Save the fixed LRI factors once
V_fixed = model.G_DK_M[1]
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
np.save(output_dir / 'lri_factors.npy', V_fixed)

del model
gc.collect()

# ---------- Prepare memmap for U_cell ----------
u_mm = open_memmap(output_dir / 'cell_loadings.npy',
                   mode='w+', dtype='float64', shape=(n_cells, K))

# ---------- Chunked row-wise updates ----------
CHUNK = 50000  # adjust if still near memory limit

for start in range(0, n_cells, CHUNK):
    end = min(start + CHUNK, n_cells)
    mat_chunk_csr = mat[start:end, :]

    coo = mat_chunk_csr.tocoo(copy=False)
    coo.data = coo.data.astype(np.float64, copy=False)
    coo.row  = coo.row.astype(np.int32, copy=False)
    coo.col  = coo.col.astype(np.int32, copy=False)
    X_chunk = sparse.COO.from_scipy_sparse(coo)
    # Manually ensure dtype is consistent
    X_chunk.coords = X_chunk.coords.astype(np.int32, copy=False)
    if getattr(X_chunk, "has_duplicates", False):
        X_chunk = X_chunk.sum_duplicates()

    # Initialize a small BPTF for this chunk
    proj = BPTF(data_shape=(end - start, n_lri), n_components=K, alpha=alpha)
    proj._init(modes=[0, 1])

    # Copy fixed column parameters
    proj.shp_DK_M[1][...] = shp1
    proj.rte_DK_M[1][...] = rte1
    proj.beta_M[1]        = bet1
    proj._update_cache(1)
    proj._clamp_component(1)

    # Update only mode=0 for this chunk
    proj._update(X_chunk, modes=[0], max_iter=MAX_ITER, verbose=False)

    # Write results to disk
    u_mm[start:end, :] = proj.G_DK_M[0]

    del proj, X_chunk, coo, mat_chunk_csr
    gc.collect()
    print(f"Updated rows {start}–{end} / {n_cells}")

# Flush memmap
del u_mm
gc.collect()

# Optional: also write the transposed version
mm = open_memmap(output_dir / 'lri_factors_KxLRI.npy',
                 mode='w+', dtype='float64', shape=(K, n_lri))
for k in range(K):
    mm[k, :] = V_fixed[:, k]
del mm
gc.collect()

print(f"Streaming BPTF projection finished. Outputs saved under: {output_dir}")
