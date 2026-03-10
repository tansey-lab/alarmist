"""
Column name constants for alarmist.

All DataFrame column names should be defined here to ensure consistency
across the codebase. Use lowercase names only.

Usage:
    from alarmist.constants import COLUMN_NAME_GENE, COLUMN_NAME_LOGFC
    df[COLUMN_NAME_LOGFC] = coefficients
"""

# =============================================================================
# GLM Differential Expression Result Columns
# =============================================================================

# Gene identifier
COLUMN_NAME_GENE = "gene"

# Log2 fold change (effect size from GLM)
COLUMN_NAME_LOGFC = "logfc"

# Standard error of the coefficient estimate
COLUMN_NAME_SE = "se"

# Raw p-value from statistical test
COLUMN_NAME_PVAL = "pval"

# FDR-adjusted p-value (q-value)
COLUMN_NAME_QVAL = "qval"

# Boolean flag indicating statistical significance
COLUMN_NAME_SIGNIFICANT = "significant"

# Motif identifier (0-indexed)
COLUMN_NAME_MOTIF = "motif"


# =============================================================================
# LRI (Ligand-Receptor Interaction) Motif Columns
# =============================================================================

# Index of the LRI in the original matrix
COLUMN_NAME_LRI_IDX = "lri_idx"

# Index of the motif from factorization
COLUMN_NAME_MOTIF_IDX = "motif_idx"

# Full LRI name string (e.g., "B_cells|T_cells|TGFB1|TGFBR1|paracrine")
COLUMN_NAME_LRI_NAME = "lri_name"

# Factor loading value from BPTF
COLUMN_NAME_FACTOR = "factor"

# Mean expression/activity value
COLUMN_NAME_MEAN = "mean"

# Sender cell type
COLUMN_NAME_CELLTYPE1 = "celltype1"

# Receiver cell type
COLUMN_NAME_CELLTYPE2 = "celltype2"

# Ligand gene name
COLUMN_NAME_LIGAND = "ligand"

# Receptor gene name (may include complex notation like "TGFBR2_TGFBR1")
COLUMN_NAME_RECEPTOR = "receptor"

# Signaling type: "paracrine", "autocrine", or "juxtacrine"
COLUMN_NAME_SIGNALING_TYPE = "signaling_type"


# =============================================================================
# Cell and Patch Metadata Columns
# =============================================================================

# Cell type annotation
COLUMN_NAME_CELL_TYPE = "cell_type"

# Sample identifier for multi-sample analyses
COLUMN_NAME_SAMPLE_ID = "sample_id"

# Patch index for spatial analyses
COLUMN_NAME_PATCH_IDX = "patch_idx"


# =============================================================================
# Derived/Computed Columns (used in plotting)
# =============================================================================

# Negative log10 of q-value for volcano plots
COLUMN_NAME_NEG_LOG10_Q = "neg_log10_q"

# Generic score column for ranking
COLUMN_NAME_SCORE = "score"
