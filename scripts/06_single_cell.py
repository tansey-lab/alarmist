import os
import gc
import numpy as np
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

INPUT_DIR   = 'results/ES_cellphone_15/single_cell'
OUTPUT_DIR  = 'results/ES_cellphone_15/single_cell/bptf'
MAX_ITER    = 200
random_state = 0
np.random.seed(random_state)

# ---------- Load data as CSR ----------
results = load_cell_lri_results(INPUT_DIR)
mat = results['cell_lri_matrix']
del results

# Keep CSR to reduce index memory
if not sp.isspmatrix_csr(mat):
    mat = mat.tocsr(copy=False)
mat.data = mat.data.astype(np.float64, copy=False)
mat.indices = mat.indices.astype(np.int32, copy=False)
mat.indptr = mat.indptr.astype(np.int64, copy=False)

n_cells, n_lri = mat.shape

# ---------- Load the patch-level model (for fixed mode=1) ----------
model = load_bptf('results/ES_cellphone_15/bptf/bptf_model.npz/bptf_1.npz')
K     = model.n_components
alpha = model.alpha

print(f"Number of cells: {n_cells}, Number of LRIs: {n_lri}, K: {K}, Alpha: {alpha}")

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
