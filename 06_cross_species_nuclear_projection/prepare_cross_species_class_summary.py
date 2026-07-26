#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BASE = str(PROJECT_ROOT / 'results/xspecies_humanmap_v1')
D = os.environ.get('CORTEX_FIGURE_DATA', os.path.join(BASE, 'figures', 'cross_species_nuclear_projection', 'data'))
OUT = os.environ.get('CORTEX_FIGURE_OUTPUT', D)
AUTH = str(PROJECT_ROOT / 'results/crossregion_v1/program_names.tsv')
renumber = _annotation_table(pd.read_csv(PROJECT_ROOT / 'results/crossregion_v1/program_renumber_map.tsv', sep='\t'))
renumber = renumber[renumber['status'].astype(str).str.lower().eq('kept') & renumber['new_P'].notna()].copy()
renumber['old_P'] = renumber['old_P'].astype(int)
renumber['new_P'] = renumber['new_P'].astype(int)
renumber = renumber.sort_values('new_P')
if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
old_to_new = dict(zip('P' + renumber['old_P'].astype(str), 'P' + renumber['new_P'].astype(str)))
auth = _annotation_table(pd.read_csv(AUTH, sep='\t'))
source_col = 'cnmf_component' if 'cnmf_component' in auth.columns else 'program'
auth = auth[[source_col, 'name_short', 'confidence']]
auth['program'] = 'P' + auth[source_col].astype(str)
auth = auth[auth['program'].isin(old_to_new)].copy()
auth['program'] = auth['program'].map(old_to_new)
auth = auth.rename(columns={'confidence': 'brain_class'})
auth['star'] = np.where(auth['brain_class'] == 'class_b', '*', '')
auth['name_short_star'] = auth['name_short'] + auth['star']
fc = _annotation_table(pd.read_csv(os.path.join(D, 'func_class_curated.csv')))
fc = fc[['program', 'function_class', 'confidence', 'h_mac_cosine', 'h_mou_cosine']]
source_programs = set(fc['program'].astype(str))
if source_programs == {f'P{i}' for i in range(1, 61)}:
    fc = fc[fc['program'].isin(old_to_new)].copy()
    fc['program'] = fc['program'].map(old_to_new)
elif source_programs == set(old_to_new):
    fc['program'] = fc['program'].map(old_to_new)
elif source_programs != {f'P{i}' for i in range(1, 55)}:
    raise ValueError()
fc.to_csv(os.path.join(OUT, 'func_class_curated_current.csv'), index=False)
df = auth.merge(fc, on='program', how='left')
if len(df) != 54 or df['program'].nunique() != 54:
    raise ValueError()
inner_map = {'glia_oligo_myelin': 'glia (oligo/astro)', 'glia_astrocyte': 'glia (oligo/astro)', 'glia_microglia_immune': 'microglia / immune', 'vascular': 'vascular', 'cytoskeleton': 'cytoskeleton', 'neuron_synapse': 'neuronal synaptic', 'neuron_ion_channel': 'neuronal synaptic', 'neuron_axon_guidance': 'neuronal synaptic', 'neuron_neuropeptide': 'neuropeptide / IEG', 'neuron_activity_IEG': 'neuropeptide / IEG', 'other_unresolved': 'unresolved', 'metabolic_housekeeping': 'unresolved'}
df['inner'] = df['function_class'].map(inner_map)
df['outer'] = np.where(df['function_class'].str.startswith('neuron'), 'neuron', 'non-neuron')
df.loc[df['function_class'].isin(['other_unresolved', 'metabolic_housekeeping']), 'outer'] = 'unresolved'
df['low_conf'] = df['confidence'] == 'low'
inner_order = df.groupby('inner')['h_mou_cosine'].median().sort_values().index.tolist()
rows = []
for (_, r) in df.iterrows():
    for (sp_lab, col) in [('猴', 'h_mac_cosine'), ('鼠', 'h_mou_cosine')]:
        rows.append({'program': r['program'], 'name_short': r['name_short_star'], 'inner': r['inner'], 'outer': r['outer'], 'confidence': r['confidence'], 'low_conf': r['low_conf'], 'species': sp_lab, 'cosine': r[col]})
long = pd.DataFrame(rows)
long['inner'] = pd.Categorical(long['inner'], categories=inner_order, ordered=True)
long.to_csv(os.path.join(OUT, 'panelC_v2_strat.csv'), index=False)
core = df[df['outer'] != 'unresolved']
stat = {}
for (axis, key) in [('h_mou_cosine', 'mou'), ('h_mac_cosine', 'mac')]:
    n = core[core['outer'] == 'neuron'][axis]
    nn = core[core['outer'] == 'non-neuron'][axis]
    (u, p) = mannwhitneyu(n, nn, alternative='two-sided')
    stat[f'{key}_p'] = p
    stat[f'{key}_neuron_med'] = float(n.median())
    stat[f'{key}_nonneuron_med'] = float(nn.median())
stat['n_neuron'] = int((core['outer'] == 'neuron').sum())
stat['n_nonneuron'] = int((core['outer'] == 'non-neuron').sum())
g = df.groupby('inner').agg(n=('program', 'size'), mou_med=('h_mou_cosine', 'median'), mac_med=('h_mac_cosine', 'median')).reindex(inner_order)
g.to_csv(os.path.join(OUT, 'panelC_v2_inner_med.csv'))
pd.DataFrame([stat]).to_csv(os.path.join(OUT, 'panelC_v2_stat.csv'), index=False)
