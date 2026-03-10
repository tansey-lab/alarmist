"""Integration tests using ovarian cancer sample data."""

import pytest
from pathlib import Path

# Path to test fixture
FIXTURE_PATH = Path(__file__).parent / "fixtures" / "ovarian_cancer_sample.h5ad"


@pytest.fixture
def sample_adata():
    """Load the test AnnData fixture."""
    import scanpy as sc
    return sc.read_h5ad(FIXTURE_PATH)


@pytest.mark.skipif(not FIXTURE_PATH.exists(), reason="Test fixture not available")
class TestPatchLRIAnalyzer:
    """Integration tests for PatchLRIAnalyzer."""

    def test_analyzer_init(self):
        """Test analyzer initialization."""
        from alarmist import PatchLRIAnalyzer

        analyzer = PatchLRIAnalyzer(
            patch_size=50.0,
            resource_name='cellchatdb'
        )
        assert analyzer.patch_size == 50.0

    def test_prepare_lri_database(self, sample_adata):
        """Test LRI database preparation with real data."""
        from alarmist import PatchLRIAnalyzer

        analyzer = PatchLRIAnalyzer(patch_size=50.0)
        lr_pairs, ligands, receptors, signaling = analyzer.prepare_lri_database(
            adata=sample_adata
        )

        assert len(lr_pairs) > 0
        assert len(ligands) == len(lr_pairs)
        assert len(receptors) == len(lr_pairs)


@pytest.mark.skipif(not FIXTURE_PATH.exists(), reason="Test fixture not available")
class TestDataLoading:
    """Test data loading utilities."""

    def test_adata_structure(self, sample_adata):
        """Verify test data has required structure."""
        # Check spatial coordinates
        assert 'spatial' in sample_adata.obsm, "Missing spatial coordinates"

        # Check cell type/label column exists (various naming conventions)
        cell_type_cols = ['cell_type', 'cell_labels', 'celltype', 'annotation']
        has_cell_type = any(col in sample_adata.obs.columns for col in cell_type_cols)
        assert has_cell_type, f"Missing cell type column. Found: {list(sample_adata.obs.columns)}"

        # Check we have cells and genes
        assert sample_adata.n_obs > 0, "No cells in data"
        assert sample_adata.n_vars > 0, "No genes in data"
