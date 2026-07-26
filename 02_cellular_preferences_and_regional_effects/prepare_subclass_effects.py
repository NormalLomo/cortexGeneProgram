#!/usr/bin/env python
import os, sys, numpy as np, pandas as pd
from pathlib import Path
np.random.seed(0)
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
RES = f'{PROJECT_ROOT}/results/crossregion_v1'
OUT = RES
PARQ = os.path.join(RES, 'cell_program_region_subclass.parquet')
renumber = pd.read_csv(os.path.join(RES, 'program_renumber_map.tsv'), sep='\t')
renumber = renumber[renumber['status'].eq('kept')].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
OLD_PROGS = renumber['old_P'].astype(str).tolist()
PROGS = [f'P{i}' for i in renumber['new_P']]
OLD_TO_NEW = dict(zip(OLD_PROGS, PROGS))
df = pd.read_parquet(PARQ)
df = df.rename(columns=OLD_TO_NEW)
if [p for p in PROGS if p not in df.columns]:
    raise ValueError()
progvar = pd.read_csv(os.path.join(RES, 'program_variability.tsv'), sep='\t')
progvar['program'] = progvar['program'].astype(str).map(lambda x: x if x.startswith('P') else OLD_TO_NEW.get(x.removeprefix('program_')))
progvar = progvar[progvar['program'].isin(PROGS)].copy()

def region_eta2(values, region):
    grand = values.mean()
    ss_tot = ((values - grand) ** 2).sum()
    if ss_tot <= 0:
        return 0.0
    ss_between = 0.0
    for (r, idx) in region.groupby(region).groups.items():
        v = values.loc[idx]
        ss_between += len(v) * (v.mean() - grand) ** 2
    return ss_between / ss_tot
top_progs = progvar.sort_values('eta2_region', ascending=False).head(15)['program'].astype(str).tolist()
rows = []
region_all = df['region']
for p in top_progs:
    v = df[p]
    grand = v.mean()
    ss_tot = ((v - grand) ** 2).sum()
    ss_between_region = 0.0
    for (r, idx) in region_all.groupby(region_all).groups.items():
        vr = v.loc[idx]
        ss_between_region += len(vr) * (vr.mean() - grand) ** 2
    eta_global = ss_between_region / ss_tot
    ss_within_subclass_region = 0.0
    for (sc, idx_sc) in df.groupby('subclass').groups.items():
        vsc = v.loc[idx_sc]
        rsc = region_all.loc[idx_sc]
        gmean_sc = vsc.mean()
        for (r, idx_r) in rsc.groupby(rsc).groups.items():
            vr = vsc.loc[idx_r]
            ss_within_subclass_region += len(vr) * (vr.mean() - gmean_sc) ** 2
    eta_within = ss_within_subclass_region / ss_tot
    eta_comp = max(eta_global - eta_within, 0.0)
    rows.append(dict(program=p, eta_global=eta_global, cell_autonomous=eta_within, compositional=eta_comp))
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'panel_c_partition.tsv'), sep='\t', index=False)
SC = 'L3-L4 IT RORB'
sub = df[df['subclass'] == SC].copy()
X = sub[PROGS].values.astype(np.float32)
from sklearn.preprocessing import StandardScaler
Xs = StandardScaler().fit_transform(X)
import scanpy as sc
ad = sc.AnnData(Xs)
ad.obs['region'] = sub['region'].values
ad.obs['P13'] = sub['P13'].values
sc.pp.pca(ad, n_comps=30)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30)
sc.tl.umap(ad, random_state=0)
um = ad.obsm['X_umap']
outf = pd.DataFrame({'UMAP1': um[:, 0], 'UMAP2': um[:, 1], 'region': sub['region'].values, 'P13': sub['P13'].values})
if len(outf) > 60000:
    outf = outf.sample(60000, random_state=0)
outf.to_csv(os.path.join(OUT, 'panel_f_umap.tsv'), sep='\t', index=False)
TOP_SC = ['L6 CT', 'L6 IT', 'ET', 'NP', 'L6B', 'AST', 'OPC', 'L3-L4 IT RORB']
N_BOOT = 20
boot_rows = []
for sc_name in TOP_SC:
    sub_g = df[df['subclass'] == sc_name]
    n = len(sub_g)
    vals = sub_g[PROGS].values
    reg = sub_g['region'].values
    full_eta = {}
    reg_ser = pd.Series(reg)
    for (j, p) in enumerate(PROGS):
        full_eta[p] = region_eta2(pd.Series(vals[:, j]), reg_ser)
    top10 = sorted(full_eta, key=full_eta.get, reverse=True)[:10]
    idx_top = [PROGS.index(p) for p in top10]
    for b in range(N_BOOT):
        bidx = np.random.randint(0, n, n)
        vb = vals[bidx][:, idx_top]
        rb = pd.Series(reg[bidx])
        for (k, p) in enumerate(top10):
            e = region_eta2(pd.Series(vb[:, k]), rb)
            boot_rows.append(dict(subclass=sc_name, program=p, boot=b, eta2=e, full_eta2=full_eta[p]))
pd.DataFrame(boot_rows).to_csv(os.path.join(OUT, 'panel_g_bootstrap.tsv'), sep='\t', index=False)
