#!/usr/bin/env python
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as mt
from sklearn.decomposition import PCA
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
OUT = str(PROJECT_ROOT / 'results/crossregion_v1')
df = pd.read_parquet(f'{OUT}/cell_program_region_subclass.parquet')
rp = pd.read_csv(f'{OUT}/region_program_mean.tsv', sep='\t', index_col=0)
z = pd.read_csv(f'{OUT}/program_region_zscore.tsv', sep='\t', index_col=0)
renumber = pd.read_csv(f'{OUT}/program_renumber_map.tsv', sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = pd.to_numeric(renumber['new_P'].astype(str).str.removeprefix('P'))
renumber = renumber.sort_values('new_P')
if renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
progs = renumber['old_P'].astype(str).tolist()
if not set(progs).issubset(set(map(str, rp.columns))):
    raise ValueError()
regs = list(rp.index)
rows = []
for p in progs:
    groups = [df.loc[df.region == r, p].dropna().values for r in regs]
    (F, pv) = stats.f_oneway(*groups)
    allvals = np.concatenate(groups)
    gm = allvals.mean()
    ssb = sum((len(g) * (g.mean() - gm) ** 2 for g in groups))
    sst = ((allvals - gm) ** 2).sum()
    eta2 = ssb / sst if sst > 0 else np.nan
    rows.append((p, F, pv, eta2))
v = pd.DataFrame(rows, columns=['program', 'F', 'p', 'eta2_region'])
v['fdr'] = mt.multipletests(v['p'].values, method='fdr_bh')[1]
cv = rp.std(axis=0, ddof=1) / rp.mean(axis=0)
cv.index = [str(i) for i in cv.index]
v['cv'] = cv.reindex(v['program']).values
v['class'] = np.where((v['fdr'] < 0.05) & (v['eta2_region'] > 0.05), 'variable', 'stable')
v.to_csv(f'{OUT}/program_variability.tsv', sep='\t', index=False)
zmat = z.loc[regs, [str(c) for c in z.columns]] if all((str(c) in [str(x) for x in z.columns] for c in progs)) else z.loc[regs]
zmat = z.loc[regs]
zmat.columns = [str(c) for c in zmat.columns]
zmat = zmat[progs]
pca = PCA(n_components=4).fit(zmat.values)
coords = pca.transform(zmat.values)
evr = pca.explained_variance_ratio_
rc = pd.DataFrame({'region': regs, 'PC1': coords[:, 0], 'PC2': coords[:, 1]})
rc.to_csv(f'{OUT}/region_pc_coords.tsv', sep='\t', index=False)
pc1 = coords[:, 0]
g = []
for p in progs:
    y = rp.loc[regs, p].values
    (r, pp) = stats.pearsonr(y, pc1)
    sl = np.polyfit(pc1, y, 1)[0]
    g.append((p, sl, r, pp))
grad = pd.DataFrame(g, columns=['program', 'axis_slope', 'axis_r', 'axis_p'])
grad.to_csv(f'{OUT}/program_gradient.tsv', sep='\t', index=False)
rc.to_csv(f'{OUT}/region_axis_gradient.tsv', sep='\t', index=False)
nvar = int((v['class'] == 'variable').sum())
nsta = int((v['class'] == 'stable').sum())
order = rc.sort_values('PC1')['region'].tolist()
