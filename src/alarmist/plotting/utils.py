"""
Plotting utility functions
"""

import matplotlib.pyplot as plt
import numpy as np


def parse_lri_full(lri_name):
    """
    Parse LRI name into components

    Parameters
    ----------
    lri_name : str
        LRI name in format: cell1|cell2|ligand|receptor|signaling_type

    Returns
    -------
    tuple
        (cell1, cell2, ligand, receptor, signaling_type)
    """
    parts = lri_name.split('|')
    if len(parts) >= 5:
        return parts[0], parts[1], parts[2], parts[3], parts[4]
    elif len(parts) == 4:
        c1, c2, ligand, receptor = parts
        mode = 'autocrine' if c1 == c2 else 'paracrine'
        return c1, c2, ligand, receptor, mode
    elif len(parts) == 2:
        return 'unknown', 'unknown', parts[0], parts[1], 'unknown'
    else:
        return 'unknown', 'unknown', lri_name, lri_name, 'unknown'


def get_cell_type_colors(unique_ct):
    """
    Generate color map for cell types

    Parameters
    ----------
    unique_ct : list
        List of unique cell types

    Returns
    -------
    dict
        Mapping from cell type to color
    """
    ct_cmap = plt.get_cmap('tab20', len(unique_ct))
    return {ct: ct_cmap(i) for i, ct in enumerate(unique_ct)}
