import click, pandas as pd
from .io import load_inputs
from . import model, binning

@click.command()
@click.option('--edges', type=click.Path(exists=True, readable=True), required=True)
@click.option('--nodes', type=click.Path(exists=True, readable=True), required=False, default=None)
@click.option('--out',   type=click.Path(), required=True)
def main(edges, nodes, out):
    import os
    os.makedirs(out, exist_ok=True)
    E = load_inputs(edges, nodes)
    # Models
    M_full = {
        'global':   model.run(E, 'global', True),
        'rsn_pairs':model.run(E, 'rsn_pairs', True),
        'nodewise': model.run(E, 'nodewise', True),
    }
    M_red = {
        'global':   model.run(E, 'global', False),
        'rsn_pairs':model.run(E, 'rsn_pairs', False),
        'nodewise': model.run(E, 'nodewise', False),
    }
    # Bins
    B_mod = {
        'global':   binning.run_binned(E, 'global', 10, True),
        'rsn_pairs':binning.run_binned(E, 'rsn_pairs', 10, True),
        'nodewise': binning.run_binned(E, 'nodewise', 10, True),
    }
    B_cor = {
        'global':   binning.bin_correlations(E, 10),
        'rsn_pairs':None,
        'nodewise': None,
    }
    # Write
    for k,v in M_full.items(): v.to_csv(os.path.join(out, f'{k}_full_py.csv'), index=False)
    for k,v in M_red.items():  v.to_csv(os.path.join(out, f'{k}_reduced_py.csv'), index=False)
    for k,v in B_mod.items():  v.to_csv(os.path.join(out, f'{k}_bins_py.csv'), index=False)
    B_cor['global'].to_csv(os.path.join(out, 'global_bin_corr_py.csv'), index=False)
    click.echo(f'Wrote results to {out}')
