import pyarrow.parquet as pq
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
base = str(PROJECT_ROOT / 'results/crossregion_v1/')
work = str(PROJECT_ROOT / 'scripts/fig2/')
renumber = pd.read_csv(base + 'program_renumber_map.tsv', sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
OLD_IDS = renumber['old_P'].tolist()
if len(OLD_IDS) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
OLD_TO_NEW = {f'program_{old}': f'P{new}' for (old, new) in zip(renumber['old_P'], renumber['new_P'])}
NEW_TO_OLD = {new: old for (old, new) in OLD_TO_NEW.items()}
PROGS = list(OLD_TO_NEW.values())
CELLTYPES = ['AST', 'CHANDELIER', 'ENDO', 'ET', 'L2-L3 IT LINC00507', 'L3-L4 IT RORB', 'L4-L5 IT RORB', 'L6 CAR3', 'L6 CT', 'L6 IT', 'L6B', 'LAMP5', 'MICRO', 'NDNF', 'NP', 'OLIGO', 'OPC', 'PAX6', 'PVALB', 'SST', 'VIP', 'VLMC']
LAYER_ORDER = ['ARACHNOID', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM']
g = pd.read_csv(work + 'prog_x_layer_global.tsv', sep='\t').set_index('majorDomain')[PROGS]
l6score = g.loc['L6'] - g.drop('L6').max()
l6prog = l6score.idxmax()
EXEMPLARS = {'P5': 'L2/3 (P5)', l6prog: f'L6 ({l6prog})', 'P33': 'OLIGO / WM (P33)', 'P51': 'ENDO / vascular (P51)', 'P26': 'Inhibitory (P26)'}
meta = pd.read_parquet(work + '_meta_cache.parquet').set_index('bin')
cov = meta.groupby('chip')['majorDomain'].nunique()
full = cov[cov == 8].index
sizes = meta[meta.chip.isin(full)].groupby('chip').size().sort_values(ascending=False)
REP = sizes.index[0]
rep_bins = meta.index[meta.chip == REP]
rep_meta = meta.loc[rep_bins, ['chip', 'x', 'y', 'region', 'majorDomain']].reset_index()
rep_meta.to_csv(work + 'repchip_meta.tsv', sep='\t', index=False)
repset = set(rep_bins)
ex_cols = [NEW_TO_OLD[p] for p in EXEMPLARS]
chunks = []
pf = pq.ParquetFile(base + 'spatial_bin50_program_score_SCT.parquet')
for b in pf.iter_batches(batch_size=400000, columns=['bin'] + ex_cols):
    d = b.to_pandas()
    d = d[d['bin'].isin(repset)]
    if len(d):
        chunks.append(d)
repprog = pd.concat(chunks)
repprog = repprog.rename(columns=OLD_TO_NEW)
repprog.to_csv(work + 'repchip_progscores.tsv', sep='\t', index=False)
chunks = []
pf = pq.ParquetFile(base + 'spatial_bin50_rctd_weights.parquet')
for b in pf.iter_batches(batch_size=400000, columns=['bin'] + CELLTYPES + ['rctd_pass_mask']):
    d = b.to_pandas()
    d = d[d['bin'].isin(repset)]
    if len(d):
        chunks.append(d)
reprctd = pd.concat(chunks)
reprctd.to_csv(work + 'repchip_rctd.tsv', sep='\t', index=False)
