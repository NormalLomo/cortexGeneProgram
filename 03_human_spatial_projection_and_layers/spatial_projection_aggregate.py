import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
base = str(PROJECT_ROOT / 'results/crossregion_v1/')
work = str(PROJECT_ROOT / 'scripts/fig2/')
LAYER_ORDER = ['ARACHNOID', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM']
renumber = pd.read_csv(base + 'program_renumber_map.tsv', sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
OLD_IDS = renumber['old_P'].tolist()
if len(OLD_IDS) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
OLD_PROGS = [f'program_{i}' for i in OLD_IDS]
PROGS = [f'P{i}' for i in renumber['new_P']]
meta = pd.read_parquet(work + '_meta_cache.parquet', columns=['bin', 'chip', 'x', 'y', 'region', 'majorDomain'])
meta = meta.set_index('bin')
pf = pq.ParquetFile(base + 'spatial_bin50_program_score_SCT.parquet')
from collections import defaultdict
chunks_cl = defaultdict(list)
chunks_g = defaultdict(list)
n_done = 0
for batch in pf.iter_batches(batch_size=400000, columns=['bin'] + OLD_PROGS):
    df = batch.to_pandas()
    df = df.set_index('bin')
    sub = meta.loc[df.index, ['chip', 'majorDomain']]
    mat = df[OLD_PROGS].to_numpy(dtype=np.float32)
    chips = sub['chip'].to_numpy()
    layers = sub['majorDomain'].to_numpy()
    key_cl = pd.Series(list(zip(chips, layers)))
    for (k, idx) in key_cl.groupby(key_cl).groups.items():
        chunks_cl[k].append(mat[np.asarray(idx)])
    ser_ly = pd.Series(layers)
    for (ly, idx) in ser_ly.groupby(ser_ly).groups.items():
        chunks_g[ly].append(mat[np.asarray(idx)])
    n_done += len(df)
rows = []
for ly in LAYER_ORDER:
    if ly in chunks_g and len(chunks_g[ly]):
        arr = np.concatenate(chunks_g[ly], axis=0)
        med = np.median(arr, axis=0)
        rows.append([ly] + list(med))
        del arr
gdf = pd.DataFrame(rows, columns=['majorDomain'] + PROGS)
gdf.to_csv(work + 'prog_x_layer_global.tsv', sep='\t', index=False)
recs = []
for ((ch, ly), clist) in chunks_cl.items():
    arr = np.concatenate(clist, axis=0)
    n = arr.shape[0]
    if n > 0:
        med = np.median(arr, axis=0)
        for (i, p) in enumerate(PROGS):
            recs.append((ch, ly, p, med[i], n))
    del arr
pcl = pd.DataFrame(recs, columns=['chip', 'majorDomain', 'program', 'mean_z', 'n'])
pcl.to_csv(work + 'prog_x_layer_per_chip.tsv', sep='\t', index=False)
G = gdf.set_index('majorDomain')[PROGS]
exemplar = []
for p in PROGS:
    col = G[p]
    exemplar.append((p, col.idxmax(), col.max() - col.min(), col.max()))
ex = pd.DataFrame(exemplar, columns=['program', 'peak_layer', 'contrast', 'peak_z']).sort_values('contrast', ascending=False)
ex.to_csv(work + 'program_layer_specificity.tsv', sep='\t', index=False)
