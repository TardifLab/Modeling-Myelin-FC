import pandas as pd

def load_inputs(edges_path, nodes_path=None):
    edges = pd.read_csv(edges_path)
    nodes = pd.read_csv(nodes_path) if nodes_path else None
    if ('rsn_i' not in edges.columns or 'rsn_j' not in edges.columns) and nodes is not None:
        nodes = nodes[['node_id','rsn']].rename(columns={'node_id':'i','rsn':'rsn_i'})
        edges = edges.merge(nodes, on='i', how='left')
        nodes = nodes.rename(columns={'i':'j','rsn_i':'rsn_j'})
        edges = edges.merge(nodes[['j','rsn_j']], on='j', how='left')
    need = ['i','j','FC','caliber','myelin','length','rsn_i','rsn_j']
    missing = [c for c in need if c not in edges.columns]
    if missing:
        raise ValueError(f"Edges table missing columns: {missing}")
    return edges
