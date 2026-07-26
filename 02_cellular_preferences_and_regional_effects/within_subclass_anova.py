#!/usr/bin/env python
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as mt
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
OUT = str(PROJECT_ROOT / 'results/crossregion_v1')
df = pd.read_parquet(f'{OUT}/cell_program_region_subclass.parquet')
renumber = pd.read_csv(f'{OUT}/program_renumber_map.tsv', sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
old_to_new = dict(zip(renumber['old_P'].astype(str), 'P' + renumber['new_P'].astype(str)))
progs = list(old_to_new.values())
df = df.rename(columns=old_to_new)
subs = sorted(df.subclass.unique())
rows = []
for s in subs:
    d = df[df.subclass == s]
    rc = d.region.value_counts()
    use = rc[rc >= 20].index.tolist()
    if len(use) < 3:
        continue
    for p in progs:
        groups = [d.loc[d.region == r, p].values for r in use]
        (F, pv) = stats.f_oneway(*groups)
        gm = d[p].mean()
        ssb = sum((len(g) * (g.mean() - gm) ** 2 for g in groups))
        sst = ((d[p] - gm) ** 2).sum()
        eta2 = ssb / sst if sst > 0 else np.nan
        rows.append((s, p, F, pv, eta2, len(use), len(d)))
w = pd.DataFrame(rows, columns=['subclass', 'program', 'F', 'p', 'eta2', 'n_regions', 'n_cells'])
w['fdr'] = mt.multipletests(w.p.fillna(1), method='fdr_bh')[1]
w = w[['subclass', 'program', 'F', 'p', 'fdr', 'eta2', 'n_regions', 'n_cells']]
w.to_csv(f'{OUT}/within_subclass_region_eta2.tsv', sep='\t', index=False)
agg = w.groupby('subclass').agg(median_eta2=('eta2', 'median'), mean_eta2=('eta2', 'mean'), n_sig_programs=('fdr', lambda x: ((x < 0.05) & (w.loc[x.index, 'eta2'] > 0.05)).sum()), total_cells=('n_cells', 'first')).sort_values('median_eta2', ascending=False)
agg.to_csv(f'{OUT}/subclass_driver_rank.tsv', sep='\t')
top = w.sort_values('eta2', ascending=False).head(10)
panel_specs = [('14', 'P13', 'L3-L4 IT RORB'), ('6', 'P6', 'L2-L3 IT LINC00507'), ('1', 'P1', 'L6 IT')]
panel_rows = []
for (raw_program, program, subclass) in panel_specs:
    subset = df.loc[df['subclass'].eq(subclass), ['region', 'subclass', program]].copy()
    subset = pd.concat([x.sample(n=min(len(x), 4000), random_state=0) for (_, x) in subset.groupby('region', sort=False)], ignore_index=True)
    subset = subset.rename(columns={program: 'act'})
    subset['program'] = raw_program
    panel_rows.append(subset[['program', 'subclass', 'region', 'act']])
pd.concat(panel_rows, ignore_index=True).to_csv(f'{OUT}/panel_g_subsample.tsv', sep='\t', index=False)
