#!/usr/bin/env python
import os, json, gc
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
import statsmodels.stats.multitest as mt
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
BASE = str(PROJECT_ROOT)
CNMF = f'{BASE}/results/cnmf_snrna_joint_full1M_v1'
CRX = f'{BASE}/results/crossregion_v1'
OUTD = f'{CRX}/k_robustness'
os.makedirs(OUTD, exist_ok=True)
OBS = f'{BASE}/inputs/snRNA_1M_obs.csv'
CS = {40: f'{CNMF}/snrna_joint_full1M_v1_k40_cell_scores.tsv', 50: f'{CNMF}/snrna_joint_full1M_v1_k50_cell_scores.tsv', 60: f'{CNMF}/snrna_joint_full1M_v1_k60_cell_scores.tsv', 65: f'{CNMF}/snrna_joint_full1M_v1_k65_k65_cell_scores.tsv', 70: f'{CNMF}/snrna_joint_full1M_v1_k70_cell_scores.tsv'}
FL = {40: f'{CNMF}/snrna_joint_full1M_v1_k40_factor_loadings.tsv', 50: f'{CNMF}/snrna_joint_full1M_v1_k50_factor_loadings.tsv', 60: f'{CNMF}/snrna_joint_full1M_v1_k60_factor_loadings.tsv', 65: f'{CNMF}/snrna_joint_full1M_v1_k65_k65_factor_loadings.tsv', 70: f'{CNMF}/snrna_joint_full1M_v1_k70_factor_loadings.tsv'}
ALTK = [40, 50, 65, 70]
DEEP_EX = {'L6 CT', 'L6 IT', 'ET', 'NP', 'L6B', 'L6 CAR3'}
GLIA = {'AST', 'OLIGO', 'OPC', 'MICRO'}
DRIVER_CLASSES = DEEP_EX | GLIA
LAMINAR_IT = {'L3-L4 IT RORB', 'L4-L5 IT RORB', 'L2-L3 IT LINC00507', 'L6 IT', 'L6 CAR3'}
RORB_IT = {'L3-L4 IT RORB', 'L4-L5 IT RORB'}
obs = pd.read_csv(OBS, index_col=0)[['region', 'subclass']]

def per_program_variability(df, progs, regs):
    rp = df.groupby('region')[progs].mean()
    rows = []
    for p in progs:
        groups = [df.loc[df.region == r, p].dropna().values for r in regs]
        (F, pv) = stats.f_oneway(*groups)
        allv = np.concatenate(groups)
        gm = allv.mean()
        ssb = sum((len(g) * (g.mean() - gm) ** 2 for g in groups))
        sst = ((allv - gm) ** 2).sum()
        eta2 = ssb / sst if sst > 0 else np.nan
        rows.append((p, F, pv, eta2))
    v = pd.DataFrame(rows, columns=['program', 'F', 'p', 'eta2_region'])
    v['fdr'] = mt.multipletests(v['p'].values, method='fdr_bh')[1]
    cv = rp.std(axis=0, ddof=1) / rp.mean(axis=0)
    cv.index = [str(i) for i in cv.index]
    v['cv'] = cv.reindex(v['program']).values
    v['class'] = np.where((v['fdr'] < 0.05) & (v['eta2_region'] > 0.05), 'variable', 'stable')
    return (v, rp)

def within_subclass_eta2(df, progs):
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
    agg = w.groupby('subclass').agg(median_eta2=('eta2', 'median'), mean_eta2=('eta2', 'mean'), n_sig_programs=('fdr', lambda x: ((x < 0.05) & (w.loc[x.index, 'eta2'] > 0.05)).sum()), total_cells=('n_cells', 'first')).sort_values('median_eta2', ascending=False)
    return (w, agg)

def program_subclass_assoc(df, progs):
    sp = df.groupby('subclass')[progs].mean()
    spz = (sp - sp.mean(axis=0)) / sp.std(axis=0)
    top = spz.idxmax(axis=0)
    return top
results = {}
for K in [60] + ALTK:
    cs = pd.read_csv(CS[K], sep='\t', index_col=0)
    progs = list(cs.columns)
    df = cs.join(obs, how='inner')
    del cs
    gc.collect()
    ov = df['region'].notna().sum()
    regs = sorted(df.region.unique())
    (v, rp) = per_program_variability(df, progs, regs)
    (w, agg) = within_subclass_eta2(df, progs)
    assoc = program_subclass_assoc(df, progs)
    results[K] = dict(var=v, agg=agg, assoc=assoc, n_cells=int(len(df)), n_prog=len(progs))
    v.to_csv(f'{OUTD}/program_variability_k{K}.tsv', sep='\t', index=False)
    agg.to_csv(f'{OUTD}/subclass_driver_rank_k{K}.tsv', sep='\t')
    assoc.rename('top_subclass').to_csv(f'{OUTD}/program_top_subclass_k{K}.tsv', sep='\t')
    nvar = int((v['class'] == 'variable').sum())
    del df, w
    gc.collect()
fl60 = pd.read_csv(FL[60], sep='\t', index_col=0)
match = {}
for K in ALTK:
    flk = pd.read_csv(FL[K], sep='\t', index_col=0)
    genes = fl60.columns.intersection(flk.columns)
    A = flk[genes].values
    B = fl60[genes].values
    Az = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
    Bz = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
    C = Az @ Bz.T / Az.shape[1]
    (ri, ci) = linear_sum_assignment(-C)
    rows = []
    for (i, j) in zip(ri, ci):
        rows.append((str(flk.index[i]), str(fl60.index[j]), float(C[i, j])))
    m = pd.DataFrame(rows, columns=['prog_K', 'match_k60', 'corr']).sort_values('corr', ascending=False)
    m.to_csv(f'{OUTD}/match_k{K}_to_k60.tsv', sep='\t', index=False)
    match[K] = m
    del flk
    gc.collect()
v60 = results[60]['var'].set_index('program')
eta60 = v60['eta2_region']
conc_rows = []
scatter = {}
for K in ALTK:
    m = match[K]
    vk = results[K]['var'].set_index('program')
    sub = []
    for (_, r) in m.iterrows():
        (pk, p60, cc) = (r['prog_K'], r['match_k60'], r['corr'])
        if pk in vk.index and p60 in eta60.index:
            sub.append((pk, vk.loc[pk, 'eta2_region'], p60, eta60.loc[p60], cc))
    sc = pd.DataFrame(sub, columns=['prog_K', 'eta_K', 'prog_60', 'eta_60', 'corr'])
    scatter[K] = sc
    (rho_all, p_all) = stats.spearmanr(sc['eta_K'], sc['eta_60'])
    scc = sc[sc['corr'] > 0.5]
    (rho_c, p_c) = stats.spearmanr(scc['eta_K'], scc['eta_60']) if len(scc) >= 4 else (np.nan, np.nan)
    (pear, _) = stats.pearsonr(sc['eta_K'], sc['eta_60'])
    conc_rows.append((K, len(sc), rho_all, p_all, len(scc), rho_c, p_c, pear, m['corr'].median()))
conc = pd.DataFrame(conc_rows, columns=['K', 'n_matched', 'spearman_all', 'p_all', 'n_conf', 'spearman_conf', 'p_conf', 'pearson_all', 'match_corr_median'])
conc.to_csv(f'{OUTD}/concordance_summary.tsv', sep='\t', index=False)
TOPN = 6
drv_rows = []
for K in [60] + ALTK:
    agg = results[K]['agg']
    topN = list(agg.index[:TOPN])
    n_deep = sum((s in DEEP_EX for s in topN))
    n_glia = sum((s in GLIA for s in topN))
    n_drv = sum((s in DRIVER_CLASSES for s in topN))
    drv_rows.append((K, TOPN, n_deep, n_glia, n_drv, n_drv / TOPN, ';'.join(topN)))
drv = pd.DataFrame(drv_rows, columns=['K', 'topN', 'n_deepEX', 'n_glia', 'n_driver', 'driver_frac', 'top_subclasses'])
drv.to_csv(f'{OUTD}/topdriver_recurrence.tsv', sep='\t', index=False)
tv_rows = []
for K in [60] + ALTK:
    v = results[K]['var']
    assoc = results[K]['assoc']
    tvp = v.sort_values('eta2_region', ascending=False).iloc[0]['program']
    tv_eta = v['eta2_region'].max()
    sc_assoc = assoc.get(tvp, 'NA')
    is_laminar = sc_assoc in LAMINAR_IT
    is_rorb = sc_assoc in RORB_IT
    top3 = v.sort_values('eta2_region', ascending=False).head(3)['program'].tolist()
    top3_sc = [assoc.get(p, 'NA') for p in top3]
    tv_rows.append((K, tvp, round(float(tv_eta), 4), sc_assoc, is_laminar, is_rorb, ';'.join((f'{p}:{s}' for (p, s) in zip(top3, top3_sc)))))
tv = pd.DataFrame(tv_rows, columns=['K', 'top_var_prog', 'top_var_eta2', 'assoc_subclass', 'is_laminar_IT', 'is_RORB', 'top3_prog_subclass'])
tv.to_csv(f'{OUTD}/top_variable_program_identity.tsv', sep='\t', index=False)
for K in ALTK:
    scatter[K].to_csv(f'{OUTD}/scatter_eta_k{K}_vs_k60.tsv', sep='\t', index=False)
summary = {'n_cells': results[60]['n_cells'], 'K_list': [60] + ALTK, 'n_programs': {str(K): results[K]['n_prog'] for K in [60] + ALTK}, 'concordance': conc.to_dict(orient='records'), 'topdriver_recurrence': drv.to_dict(orient='records'), 'top_variable_identity': tv.to_dict(orient='records')}
with open(f'{OUTD}/k_robustness_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
