"""
Column name constants for alarmist.

All DataFrame and adata.obs column names should be defined here to ensure
consistency across the codebase. Use lowercase names only (with rare exceptions
that match upstream library conventions, e.g. W_max).

Usage:
    from alarmist.constants import COLUMN_NAME_GENE, COLUMN_NAME_LOGFC
    df[COLUMN_NAME_LOGFC] = coefficients
"""

# =============================================================================
# GLM Differential Expression Result Columns
# =============================================================================

COLUMN_NAME_GENE = "gene"
COLUMN_NAME_LOGFC = "logfc"
COLUMN_NAME_SE = "se"
COLUMN_NAME_PVAL = "pval"
COLUMN_NAME_QVAL = "qval"
COLUMN_NAME_SIGNIFICANT = "significant"
COLUMN_NAME_MOTIF = "motif"

# scanpy/rank_genes_groups output column names
COLUMN_NAME_P_ADJ = "p_adj"
COLUMN_NAME_LOGFOLDCHANGES = "logfoldchanges"


# =============================================================================
# LRI (Ligand-Receptor Interaction) Motif Columns
# =============================================================================

COLUMN_NAME_LRI_IDX = "lri_idx"
COLUMN_NAME_MOTIF_IDX = "motif_idx"
COLUMN_NAME_LRI_NAME = "lri_name"
COLUMN_NAME_FACTOR = "factor"
COLUMN_NAME_MEAN = "mean"
COLUMN_NAME_CELLTYPE1 = "celltype1"
COLUMN_NAME_CELLTYPE2 = "celltype2"
COLUMN_NAME_LIGAND = "ligand"
COLUMN_NAME_RECEPTOR = "receptor"
COLUMN_NAME_SIGNALING_TYPE = "signaling_type"
COLUMN_NAME_PATHWAY = "pathway"
COLUMN_NAME_MOTIF_NAME = "motif_name"


# =============================================================================
# Cell and Patch Metadata Columns
# =============================================================================

COLUMN_NAME_CELL_TYPE = "cell_type"
COLUMN_NAME_CELL_ID = "cell_id"
COLUMN_NAME_SAMPLE_ID = "sample_id"
COLUMN_NAME_PATCH_IDX = "patch_idx"
COLUMN_NAME_TMA_ID = "tma_id"


# =============================================================================
# Patch/Sample Metadata Table Columns
# =============================================================================

COLUMN_NAME_N_CELLS = "n_cells"
COLUMN_NAME_N_PATCHES = "n_patches"
COLUMN_NAME_GLOBAL_PATCH_IDX_START = "global_patch_idx_start"
COLUMN_NAME_GLOBAL_PATCH_IDX_END = "global_patch_idx_end"
COLUMN_NAME_GLOBAL_CELL_IDX_START = "global_cell_idx_start"
COLUMN_NAME_GLOBAL_CELL_IDX_END = "global_cell_idx_end"
COLUMN_NAME_AVG_NEIGHBORHOOD_SIZE = "avg_neighborhood_size"


# =============================================================================
# Generic Columns Used Across Tables
# =============================================================================

COLUMN_NAME_COLUMN_NAME = "column_name"
COLUMN_NAME_PARAMETER = "parameter"
COLUMN_NAME_VALUE = "value"
COLUMN_NAME_WEIGHT = "weight"
COLUMN_NAME_SOURCE = "source"
COLUMN_NAME_TARGET = "target"
COLUMN_NAME_MODE = "mode"


# =============================================================================
# Factorization Iteration History Columns
# =============================================================================

COLUMN_NAME_ITERATION = "iteration"
COLUMN_NAME_ELBO = "elbo"
COLUMN_NAME_DELTA = "delta"


# =============================================================================
# Derived/Computed Columns (used in plotting and scoring)
# =============================================================================

COLUMN_NAME_NEG_LOG10_Q = "neg_log10_q"
COLUMN_NAME_SCORE = "score"
COLUMN_NAME_W_MAX = "W_max"
COLUMN_NAME_FACTOR_RESCALED = "factor_rescaled"
COLUMN_NAME_FACTOR_LRNORM = "factor_lrnorm"
COLUMN_NAME_LR_GLOBAL_MEAN = "lr_global_mean"
COLUMN_NAME_CELL_PAIR = "cell_pair"
COLUMN_NAME_IN_TOP20 = "in_top20"
