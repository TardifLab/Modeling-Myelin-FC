# myelinfccoupling/io.py

import pandas as pd


def load_inputs(edges_path, nodes_path=None) -> pd.DataFrame:
    """
    Load edges (and optionally nodes) tables, and ensure required columns.

    Required edge columns after processing:
    - i, j          (node indices)
    - FC            (functional connectivity)
    - caliber       (edge caliber)
    - myelin        (edge myelin)
    - length        (edge length)
    - rsn_i, rsn_j  (RSN labels for each endpoint)
    """
    edges = pd.read_csv(edges_path)
    nodes = pd.read_csv(nodes_path) if nodes_path else None

    # If rsn_i / rsn_j are missing but we have nodes, merge RSN labels
    if (("rsn_i" not in edges.columns) or ("rsn_j" not in edges.columns)) and nodes is not None:
        # Expect nodes to have node_id and rsn
        if not {"node_id", "rsn"} <= set(nodes.columns):
            raise ValueError(
                "Nodes table must contain columns 'node_id' and 'rsn' "
                "to derive rsn_i / rsn_j."
            )
        nodes_i = nodes[["node_id", "rsn"]].rename(
            columns={"node_id": "i", "rsn": "rsn_i"}
        )
        edges = edges.merge(nodes_i, on="i", how="left")

        nodes_j = nodes[["node_id", "rsn"]].rename(
            columns={"node_id": "j", "rsn": "rsn_j"}
        )
        edges = edges.merge(nodes_j, on="j", how="left")

    need = [
        "i",
        "j",
        "FC",
        "caliber",
        "myelin",
        "length",
        "rsn_i",
        "rsn_j",
    ]
    missing = [c for c in need if c not in edges.columns]
    if missing:
        raise ValueError(f"Edges table missing columns: {missing}")

    return edges
