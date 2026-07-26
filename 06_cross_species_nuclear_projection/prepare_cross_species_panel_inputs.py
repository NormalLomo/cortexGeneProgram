#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os, json, re
import numpy as np
import pandas as pd
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BASE = str(PROJECT_ROOT / 'results/xspecies_humanmap_v1')
CROSS = str(PROJECT_ROOT / 'results/crossregion_v1')
OUT = os.environ.get('CORTEX_FIGURE_OUTPUT', os.path.join(BASE, 'figures', 'cross_species_nuclear_projection', 'data'))
os.makedirs(OUT, exist_ok=True)
renumber = _annotation_table(pd.read_csv(os.path.join(CROSS, 'program_renumber_map.tsv'), sep='\t'))
renumber = renumber[renumber['status'].astype(str).str.lower().eq('kept') & renumber['new_P'].notna()].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
if renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
old_to_new = dict(zip('P' + renumber['old_P'].astype(str), 'P' + renumber['new_P'].astype(str)))
new_to_old = {new: old for (old, new) in old_to_new.items()}

def retain_programs(table):
    table = table[table['program'].isin(old_to_new)].copy()
    table['source_program'] = table['program']
    table['program'] = table['program'].map(old_to_new)
    if table['program'].nunique() != 54:
        raise ValueError()
    return table
nm = _annotation_table(pd.read_csv(os.path.join(CROSS, 'program_names.tsv'), sep='\t'))
source_col = 'cnmf_component' if 'cnmf_component' in nm.columns else 'program'
nm['source_program'] = 'P' + nm[source_col].astype(str)
nm = nm[nm['source_program'].isin(old_to_new)].copy()
nm['program'] = nm['source_program'].map(old_to_new)
name_map = dict(zip(nm['program'], nm['program'] + ' ' + nm['name_short'].astype(str)))

def short_name(p, raw):
    nm = name_map.get(p, raw)
    if not isinstance(nm, str) or nm.strip() == '':
        nm = raw
    return nm
sig = _annotation_table(pd.read_csv(os.path.join(BASE, 'conservation_significance_per_program.csv')))
sig = retain_programs(sig)
sig['func_short'] = [short_name(p, f) for (p, f) in zip(sig['program'], sig['func_name'])]
strat = _annotation_table(pd.read_csv(os.path.join(BASE, 'conservation_function_stratification.csv')))
strat = retain_programs(strat)
decay = _annotation_table(pd.read_csv(os.path.join(BASE, 'decay_per_program_full.csv')))
decay = retain_programs(decay)
a = decay.melt(id_vars='program', value_vars=['human_macaque', 'human_mouse'], var_name='pair', value_name='cosine')
pair_label = {'human_macaque': 'human-macaque', 'human_mouse': 'human-mouse', 'mouse_macaque': 'mouse-macaque'}
a['pair'] = a['pair'].map(pair_label)
a.to_csv(os.path.join(OUT, 'panelA_pairs_long.csv'), index=False)
med = a.groupby('pair')['cosine'].median().reset_index().rename(columns={'cosine': 'median'})
med.to_csv(os.path.join(OUT, 'panelA_medians.csv'), index=False)
b = sig[['program', 'func_short', 'h_mac_cosine', 'h_mac_sig', 'h_mac_fdr']].copy()
b = b.sort_values('h_mac_cosine', ascending=True).reset_index(drop=True)
b['order'] = np.arange(len(b))

def stars(fdr):
    if fdr < 0.001:
        return '***'
    if fdr < 0.01:
        return '**'
    if fdr < 0.05:
        return '*'
    return ''
b['star'] = b['h_mac_fdr'].apply(stars)
b['label'] = b['func_short'] + ' ' + b['star']
b['label'] = b['label'].str.strip()
b.to_csv(os.path.join(OUT, 'panelB_lollipop.csv'), index=False)
c = strat[['program', 'func_class', 'brain_relevance', 'h_mac_cosine', 'conservation_tier']].copy()
c['tier'] = c['conservation_tier'].where(c['conservation_tier'].isin(['class_a', 'class_b']), 'unknown')
c.to_csv(os.path.join(OUT, 'panelC_strat.csv'), index=False)
with open(os.path.join(BASE, 'conservation_null_summary.json')) as f:
    nullj = json.load(f)
pd.DataFrame([{'class_a_median': nullj['class_a_mac_median'], 'class_b_median': nullj['class_b_mac_median'], 'MWU_p': nullj['class_a_vs_b_mwu_p'], 'n_class_a': nullj['n_class_a'], 'n_class_b': nullj['n_class_b']}]).to_csv(os.path.join(OUT, 'panelC_stat.csv'), index=False)
bottom_prog = sig[sig['confidence'] == 'class_a'].sort_values('h_mac_cosine').iloc[0]['program']
cases = ['P8', 'P1', 'P14', bottom_prog]
mou = np.load(os.path.join(BASE, 'loadings_mouse.npz'), allow_pickle=True)
mac = np.load(os.path.join(BASE, 'loadings_macaque.npz'), allow_pickle=True)
progs_mou = list(mou['programs'])
progs_mac = list(mac['programs'])
rows = []
for (sp, dat, progs) in [('macaque', mac, progs_mac), ('mouse', mou, progs_mou)]:
    H = dat['Hload']
    R = dat['refit']
    for p in cases:
        i = progs.index(new_to_old[p])
        h = H[i]
        r = R[i]
        for (hv, rv) in zip(h, r):
            rows.append({'program': p, 'species': sp, 'human_loading': float(hv), 'species_refit': float(rv)})
d = pd.DataFrame(rows)

def subsample(g, n=1500):
    if len(g) <= n:
        return g
    top = g.reindex(g['human_loading'].abs().sort_values(ascending=False).index).head(n // 2)
    rest = g.drop(top.index).sample(n - len(top), random_state=0)
    return pd.concat([top, rest])
d = d.groupby(['program', 'species'], group_keys=False).apply(subsample)
d.to_csv(os.path.join(OUT, 'panelD_gene_scatter.csv'), index=False)
cosrows = []
for p in cases:
    cosrows.append({'program': p, 'func_short': sig.loc[sig.program == p, 'func_short'].values[0], 'cos_mac': float(sig.loc[sig.program == p, 'h_mac_cosine'].values[0]), 'cos_mou': float(sig.loc[sig.program == p, 'h_mou_cosine'].values[0])})
pd.DataFrame(cosrows).to_csv(os.path.join(OUT, 'panelD_case_cos.csv'), index=False)
e = sig[['program', 'func_short', 'h_mac_cosine', 'h_mac_fdr', 'h_mac_sig', 'h_mac_diag_rank', 'h_mou_cosine', 'h_mou_fdr', 'h_mou_sig', 'h_mou_diag_rank']].copy()
e = e.sort_values('h_mac_cosine', ascending=False).reset_index(drop=True)
e.to_csv(os.path.join(OUT, 'panelE_significance.csv'), index=False)
sc = _annotation_table(pd.read_csv(os.path.join(CROSS, 'region_subclass_program_mean.tsv'), sep='\t'))
sc['program'] = 'P' + sc['program'].astype(str)
sc = retain_programs(sc)
psc = sc.groupby(['program', 'subclass'])['mean'].mean().reset_index()

def zscore(g):
    m = g['mean'].mean()
    s = g['mean'].std(ddof=0)
    g['z'] = (g['mean'] - m) / (s if s > 0 else 1)
    return g
psc = psc.groupby('program', group_keys=False).apply(zscore)
psc = psc.merge(sig[['program', 'h_mac_cosine']], on='program', how='left')
psc.to_csv(os.path.join(OUT, 'panelF_program_subclass.csv'), index=False)
dom = _annotation_table(pd.read_csv(os.path.join(CROSS, 'cluster_confusion', 'view2_domain_mean.tsv'), sep='\t'))
dom_long = dom.melt(id_vars='majorDomain', var_name='program', value_name='score')
dom_long['program'] = 'P' + dom_long['program'].astype(str)
dom_long = retain_programs(dom_long)
dom_long = dom_long.merge(sig[['program', 'h_mac_cosine']], on='program', how='left')
dom_long.to_csv(os.path.join(OUT, 'panelG_program_laminar.csv'), index=False)
unc = strat[['program', 'func_name', 'func_class', 'brain_relevance', 'conservation_tier']].copy()
unc['func_short'] = [short_name(p, f) for (p, f) in zip(unc['program'], unc['func_name'])]
unc['other'] = (unc['brain_relevance'] == 'other') | unc['conservation_tier'].isin(['True']) | ~unc['func_class'].isin(['glia_vascular_cytoskeleton_immune', 'neuron_synapse_axon'])
other_list = unc[unc['other']][['program', 'func_short', 'func_class', 'brain_relevance']]
other_list.to_csv(os.path.join(OUT, 'func_attribution_other.csv'), index=False)
