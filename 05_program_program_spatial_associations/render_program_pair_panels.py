#!/usr/bin/env python3
import os
from pathlib import Path
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
RESULTS = ROOT / 'results/crossregion_v1'
FINAL = RESULTS / 'markcorr_v2/final'
INPUT_DIR = Path(os.environ['MARKCORR_INPUT_DIR'])
BETWEEN = RESULTS / 'markcorr_betweenchip_v1'
PROGRAM_MAP = RESULTS / 'program_renumber_map.tsv'
OUTPUT = Path(os.environ.get('CORTEX_FIGURE_OUTPUT', ROOT / 'figures/fig5'))
mapping = pd.read_csv(PROGRAM_MAP, sep='\t')
mapping = mapping[mapping['status'].eq('kept')].copy()
mapping['old_P'] = mapping['old_P'].astype(int)
mapping['new_P'] = mapping['new_P'].astype(int)
old_to_new = dict(zip(mapping['old_P'], mapping['new_P']))
new_to_old = dict(zip(mapping['new_P'], mapping['old_P']))
name_by_new = dict(zip(mapping['new_P'], mapping['functional_name']))
if len(old_to_new) != 54 or sorted(old_to_new.values()) != list(range(1, 55)):
    raise ValueError()

def old_program(values):
    return values.astype(str).str.replace('program_', '', regex=False).astype(int)

def pair_label(new_a, new_b):
    return f'P{new_a} {name_by_new[new_a]} × P{new_b} {name_by_new[new_b]}'

def save_figure(fig, stem):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT / f'{stem}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
between = pd.read_csv(BETWEEN / 'betweenchip_progprog_stouffer_q.tsv', sep='\t')
between['new_a'] = old_program(between['A_name'])
between['new_b'] = old_program(between['B_name'])
between = between[between['new_a'].between(1, 54) & between['new_b'].between(1, 54)].copy()
between[['new_a', 'new_b']] = np.sort(between[['new_a', 'new_b']].to_numpy(), axis=1)
between['headline'] = between['is_headline'].astype(str).str.lower().eq('true')
loo = pd.read_csv(INPUT_DIR / 'progprog_donor_loo.tsv', sep='\t')
loo['old_a'] = loo['A_idx'].astype(int) + 1
loo['old_b'] = loo['B_idx'].astype(int) + 1
loo = loo[loo['old_a'].isin(old_to_new) & loo['old_b'].isin(old_to_new)].copy()
loo['new_a'] = loo['old_a'].map(old_to_new)
loo['new_b'] = loo['old_b'].map(old_to_new)
loo[['new_a', 'new_b']] = np.sort(loo[['new_a', 'new_b']].to_numpy(), axis=1)
loo = loo.sort_values('full_log2', key=lambda x: x.abs(), ascending=False).head(22).copy()
loo['label'] = [pair_label(a, b) for (a, b) in zip(loo['new_a'], loo['new_b'])]
loo['lo'] = np.log2(loo['loo_min_median_g'].clip(lower=np.finfo(float).tiny))
loo['hi'] = np.log2(loo['loo_max_median_g'].clip(lower=np.finfo(float).tiny))
loo = loo.sort_values('full_log2')
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 6, 'axes.titlesize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
(fig, ax) = plt.subplots(figsize=(5.4, 4.5), dpi=300)
y = np.arange(len(loo))
stable = loo['stable_90'].astype(str).str.lower().eq('true')
ax.hlines(y, loo['lo'], loo['hi'], color='#7f8c8d', linewidth=1)
ax.scatter(loo['full_log2'], y, c=np.where(stable, '#2e933c', '#777777'), s=15, zorder=3)
ax.axvline(0, color='#888888', linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(loo['label'], fontsize=4.2)
ax.set_xlabel('log2 g, with donor omission range')
ax.set_title('Section and donor stability of selected program pairs', loc='left')
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, 'Fig5i')
byarea = pd.read_csv(FINAL / 'progprog_byarea_median_iqr.tsv', sep='\t')
byarea['old_a'] = old_program(byarea['A'])
byarea['old_b'] = old_program(byarea['B'])
byarea = byarea[byarea['old_a'].isin(old_to_new) & byarea['old_b'].isin(old_to_new)].copy()
byarea['new_a'] = byarea['old_a'].map(old_to_new)
byarea['new_b'] = byarea['old_b'].map(old_to_new)
byarea[['new_a', 'new_b']] = np.sort(byarea[['new_a', 'new_b']].to_numpy(), axis=1)
byarea = byarea[byarea['ring_um'].eq(25)].copy()
display_pairs = [(33, 34), (8, 32), (36, 26), (16, 43)]
area_order = ['DLPFC', 'SMG', 'M1', 'V1', 'AG', 'VLPFC', 'SPL', 'FPPFC', 'S1']
(fig, axes) = plt.subplots(1, 4, figsize=(7.1, 2.65), dpi=300, sharey=True)
for (ax, (new_a, new_b)) in zip(axes, display_pairs):
    (a, b) = sorted((new_a, new_b))
    data = byarea[byarea['new_a'].eq(a) & byarea['new_b'].eq(b) & byarea['area'].isin(area_order)].copy()
    data['area'] = pd.Categorical(data['area'], categories=area_order, ordered=True)
    data = data.sort_values('area')
    colors = np.where(data['log2_median_g'].to_numpy() >= 0, '#b2182b', '#2166ac')
    ax.bar(np.arange(len(data)), data['log2_median_g'], color=colors, width=0.7)
    ax.axhline(0, color='#888888', linewidth=0.4)
    ax.set_xticks(np.arange(len(data)))
    ax.set_xticklabels(data['area'], rotation=60, ha='right', fontsize=4.2)
    ax.set_title(f'P{a} × P{b}')
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('log2 median g')
fig.tight_layout()
save_figure(fig, 'Fig5j')
headline = between[between['headline']].copy()
headline = headline.sort_values('median_log2g', key=lambda x: x.abs(), ascending=False).head(40)
headline_keys = set(zip(headline['new_a'], headline['new_b']))
regional = byarea[byarea['unstable'].eq(0)].copy()
regional = regional[[pair in headline_keys for pair in zip(regional['new_a'], regional['new_b'])]]
ranges = regional.groupby(['new_a', 'new_b'], as_index=False)['log2_median_g'].agg(lo='min', hi='max')
ranges['span'] = ranges['hi'] - ranges['lo']
ranges = ranges.sort_values('span', ascending=False).head(22).sort_values('span')
ranges['label'] = [pair_label(a, b) for (a, b) in zip(ranges['new_a'], ranges['new_b'])]
(fig, ax) = plt.subplots(figsize=(5.4, 4.7), dpi=300)
y = np.arange(len(ranges))
ax.hlines(y, ranges['lo'], ranges['hi'], color='#777777', linewidth=1)
ax.scatter(ranges['lo'], y, c='#2166ac', s=12)
ax.scatter(ranges['hi'], y, c='#b2182b', s=12)
ax.axvline(0, color='#888888', linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(ranges['label'], fontsize=4.1)
ax.set_xlabel('regional range of log2 median g')
ax.set_title('Program pairs with the widest regional variation', loc='left')
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, 'Fig5k')
