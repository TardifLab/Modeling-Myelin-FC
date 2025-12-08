import numpy as np, pandas as pd
from . import model

def bin_by_myelin(edges: pd.DataFrame, nbins=10):
    q = np.linspace(0,1,nbins+1)
    edges = edges.copy()
    edges['my_bin'] = pd.qcut(edges['myelin'], q, duplicates='drop')
    return edges

def run_binned(edges: pd.DataFrame, level: str, nbins=10, with_interactions=True):
    edges_b = bin_by_myelin(edges, nbins)
    out = []
    for b, sub in edges_b.groupby('my_bin'):
        M = model.run(sub, level, with_interactions)
        M['bin'] = str(b)
        out.append(M)
    return pd.concat(out, ignore_index=True)

def bin_correlations(edges: pd.DataFrame, nbins=10):
    edges_b = bin_by_myelin(edges, nbins)
    rows = []
    for b, sub in edges_b.groupby('my_bin'):
        R = {'bin': str(b)}
        for v in ['caliber','myelin','length']:
            R[f'{v}_r'] = sub['FC'].corr(sub[v])
            R[f'{v}_n'] = sub[['FC',v]].dropna().shape[0]
        rows.append(R)
    return pd.DataFrame(rows)
