#!/usr/bin/env python
import pyarrow.parquet as pq
import pandas as pd
import os
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
PROJ = str(PROJECT_ROOT)
EMB = f'{PROJ}/figures/fig2/_intermediate/spatial_umap_expr_harmony.csv'
SCT = f'{PROJ}/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet'
OUT = f'{PROJ}/figures/fig2/_intermediate/spatial_umap_feat.csv'
renumber = pd.read_csv(f'{PROJ}/results/crossregion_v1/program_renumber_map.tsv', sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
old_to_new = dict(zip(renumber['old_P'], renumber['new_P']))
if len(old_to_new) != 54 or sorted(old_to_new.values()) != list(range(1, 55)):
    raise ValueError()
FEAT_OLD = [37, 34, 7]
FEAT = [f'program_{i}' for i in FEAT_OLD]
FEATURE_NAMES = {f'program_{i}': f'P{old_to_new[i]}' for i in FEAT_OLD}
emb_bins = set(pd.read_csv(EMB, usecols=['bin'])['bin'].astype(str))
pf = pq.ParquetFile(SCT)
parts = []
for b in pf.iter_batches(batch_size=400000, columns=['bin'] + FEAT):
    d = b.to_pandas()
    d['bin'] = d['bin'].astype(str)
    d = d[d['bin'].isin(emb_bins)]
    if len(d):
        parts.append(d)
out = pd.concat(parts)
out = out.rename(columns=FEATURE_NAMES)
out.to_csv(OUT, index=False)
