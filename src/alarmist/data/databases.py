"""
LRI database management functions
"""

import pandas as pd
from liana.resource import select_resource


def load_lr_database(
    resource_name: str,
    cellchatdb_path: str | None = None,
    cellphonedb_path: str | None = None,
) -> pd.DataFrame:
    """
    Load ligand-receptor database

    Parameters
    ----------
    resource_name : str
        Name of LRI database ('cellchatdb', 'cellphonedb', or liana resource name)
    cellchatdb_path : str, optional
        Path to local CellChatDB CSV file
    cellphonedb_path : str, optional
        Path to local CellPhoneDB CSV file

    Returns
    -------
    pd.DataFrame
        Database with at least 'ligand', 'receptor', 'signaling_type' columns
    """
    print(f"Loading {resource_name} database...")

    # Load from local CSV if cellchatdb or cellphonedb
    if resource_name.lower() == "cellchatdb" and cellchatdb_path:
        resource = pd.read_csv(cellchatdb_path)
        # Check required columns
        required_cols = ["ligand", "receptor", "signaling_type"]
        if not all(col in resource.columns for col in required_cols):
            raise ValueError(
                f"CellChatDB CSV must contain {required_cols}. "
                f"Found columns: {resource.columns.tolist()}"
            )
    elif resource_name.lower() == "cellphonedb" and cellphonedb_path:
        resource = pd.read_csv(cellphonedb_path)
        # Check required columns
        required_cols = ["ligand", "receptor", "signaling_type"]
        if not all(col in resource.columns for col in required_cols):
            raise ValueError(
                f"CellPhoneDB CSV must contain {required_cols}. "
                f"Found columns: {resource.columns.tolist()}"
            )
    else:
        # Use liana's select_resource for other databases
        resource = select_resource(resource_name)
        # LIANA doesn't have signaling_type, add it as 'Unknown'
        if "signaling_type" not in resource.columns:
            resource["signaling_type"] = "Unknown"

    return resource
