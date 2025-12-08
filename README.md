# Myelin–FC Coupling (MATLAB + Python)

Pipeline to assess how white-matter **myelin** relates to **structural–functional coupling**.

## What’s here
- **DATA** derivative connectivity data in 200 & 400 node parcellation — see `data/`
- **MATLAB** main analysis code (used in paper) — see `matlab/`
- **Python** companion code mirroring logic of Matlab package — see `python/`

## Quick start (MATLAB)
```matlab
addpath(genpath('matlab'));
mk_edgescsv_from_mat()            % edit data in here
cfg = config.default_config();    % edit paths in here
main;                             % runs full pipeline
```

## Quick start (Python)
```bash
cd python
pip install -e .
myelin-fc-run --edges ../data/edges_fc_BOLDin.csv --out ../out --fc-label BOLD
```

## Inputs (edges table)
Columns required in `edges.csv`:
`i,j,FC,caliber,myelin,length` and optional `rsn_i,rsn_j`.
If RSNs are not on edges, provide `nodes.csv` with `node_id,rsn`.

## Outputs
All results saved to `out/` as:
- subdirectories by FC modality 
- CSVs by modeling resolution / analysis stage
+ a single `results_all.mat` (MATLAB only).

---

2025 Mark C Nelson – MNI — GNU License (code), see `LICENSE`.
