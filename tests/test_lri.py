"""Tests for ligand-receptor interaction module."""

import pandas as pd


def test_import():
    """Test that alarmist can be imported."""
    import alarmist

    assert alarmist is not None


def test_lri_import():
    """Test that core LRI module can be imported."""
    from alarmist.core import lri

    assert lri is not None


def test_extract_lri_genes():
    """Test extract_lri_genes function."""
    from alarmist.core.glm import extract_lri_genes

    # Test with Series input
    # Expected format: sender_celltype | receiver_celltype | ligand | receptor | mode
    lri_names = pd.Series(
        [
            "CellTypeA|CellTypeB|LIGAND1|RECEPTOR1|secreted",
            "CellTypeA|CellTypeC|LIGAND2|RECEPTOR2|secreted",
            "CellTypeB|CellTypeC|LIGAND3|RECEPTOR3|contact",
        ]
    )

    genes = extract_lri_genes(lri_names)
    assert "LIGAND1" in genes
    assert "RECEPTOR1" in genes
    assert len(genes) == 6

    # Test with list input
    genes_list = extract_lri_genes(lri_names.tolist())
    assert len(genes_list) == 6
