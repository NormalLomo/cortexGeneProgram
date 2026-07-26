#!/usr/bin/env python3
import os
from pathlib import Path
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Wedge
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
AGING = ROOT / 'results/crossregion_v1/program_aging'
OUTPUT = ROOT / 'figures/fig8'
FDR = pd.read_csv(AGING / 'fdr_matrix_N150.tsv', sep='\t').set_index('new_P')
LONG = pd.read_csv(AGING / 'enrichment_long_N150.tsv', sep='\t')
TOP = pd.read_csv(AGING / 'top30_pairs_N150.tsv', sep='\t')
SIG = pd.read_csv(AGING / 'significant_pairs_N150.tsv', sep='\t')
COUNTS = pd.read_csv(AGING / 'category_program_counts_N150.tsv', sep='\t')
SUMMARY = pd.read_csv(AGING / 'program_aging_summary_N150.tsv', sep='\t')
expected = {f'P{i}' for i in range(1, 55)}
if set(FDR.index) != expected:
    raise ValueError()
LONG = LONG[LONG['new_P'].isin(expected)].copy()
TOP = TOP[TOP['new_P'].isin(expected)].copy()
SIG = SIG[SIG['new_P'].isin(expected)].copy()
SUMMARY = SUMMARY[SUMMARY['new_P'].isin(expected)].copy()
for frame in [LONG, TOP, SIG, SUMMARY]:
    if not set(frame['new_P'].dropna()).issubset(expected):
        raise ValueError()
SET_LABEL = {'HAGR_GenAge_human': 'GenAge', 'HAGR_CellAge_induces_senescence': 'CellAge+', 'HAGR_CellAge_inhibits_senescence': 'CellAge-', 'SAUL_SEN_MAYO': 'SenMayo', 'LIU2026_CELL_ACCELERATED_AGING_PROTEINS': 'AA proteins', 'LIU2026_CELL_DECELERATED_AGING_PROTEINS': 'DA proteins'}
SET_ORDER = ['LIU2026_CELL_ACCELERATED_AGING_PROTEINS', 'LIU2026_CELL_DECELERATED_AGING_PROTEINS', 'SAUL_SEN_MAYO', 'HAGR_CellAge_induces_senescence', 'HAGR_CellAge_inhibits_senescence', 'HAGR_GenAge_human']
CATEGORY_COLORS = {'Accelerated aging proteome': '#b4493f', 'Decelerated aging proteome': '#2b6cb0', 'SASP/SenMayo': '#8064a2', 'Cellular senescence': '#c9823b', 'Aging curated': '#4e7e59'}
CLASS_COLORS = {'nonneuron': '#a65e5e', 'vascular': '#c9823b', 'glia': '#4e9a72', 'exc': '#4a6fa5', 'inh': '#8a65a8'}

def short_name(value, length=24):
    value = str(value)
    return value if len(value) <= length else value[:length - 1] + '…'
sig_order = SIG.sort_values('fdr_retained54')['new_P'].drop_duplicates().tolist()
row_names = SIG.sort_values('fdr_retained54').drop_duplicates('new_P').set_index('new_P').reindex(sig_order)['name_short']
matrix = FDR.reindex(index=sig_order, columns=SET_ORDER)
odds = LONG.pivot(index='new_P', columns='aging_set', values='odds_ratio').reindex(index=sig_order, columns=SET_ORDER)
sig_keys = set(zip(SIG['new_P'], SIG['aging_set']))
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 6, 'axes.titlesize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
fig = plt.figure(figsize=(7.1, 5.25), dpi=300)
axes = {'f': fig.add_axes([0.08, 0.54, 0.42, 0.4]), 'g': fig.add_axes([0.64, 0.6, 0.3, 0.34]), 'h': fig.add_axes([0.04, 0.08, 0.28, 0.34]), 'i': fig.add_axes([0.36, 0.07, 0.34, 0.34]), 'j': fig.add_axes([0.75, 0.08, 0.22, 0.34])}
cmap = LinearSegmentedColormap.from_list('aging', ['#f5f5f2', '#f1d2b9', '#d97854', '#a9363e', '#5d183a'])
neglog = -np.log10(matrix.clip(lower=1e-300))
vmax = max(4.0, float(np.nanpercentile(neglog.to_numpy(), 98)))
norm = Normalize(0, vmax)
ax = axes['f']
for x in range(len(SET_ORDER) + 1):
    ax.axvline(x - 0.5, color='#e7e7e4', linewidth=0.3, zorder=0)
for y in range(len(sig_order) + 1):
    ax.axhline(y - 0.5, color='#e7e7e4', linewidth=0.3, zorder=0)
(xs, ys, cs, sizes) = ([], [], [], [])
for (y, program) in enumerate(sig_order):
    for (x, aging_set) in enumerate(SET_ORDER):
        xs.append(x)
        ys.append(y)
        cs.append(neglog.loc[program, aging_set])
        odds_value = odds.loc[program, aging_set]
        sizes.append(8 + min(float(odds_value), 10) * 4 if (program, aging_set) in sig_keys else 0)
points = ax.scatter(xs, ys, c=cs, s=sizes, cmap=cmap, norm=norm, edgecolors='#333333', linewidths=0.15)
ax.set_xlim(-0.5, len(SET_ORDER) - 0.5)
ax.set_ylim(len(sig_order) - 0.5, -0.5)
ax.set_xticks(range(len(SET_ORDER)))
ax.set_xticklabels([SET_LABEL[s] for s in SET_ORDER], rotation=40, ha='right')
ax.set_yticks(range(len(sig_order)))
ax.set_yticklabels([f'{p} {short_name(row_names[p])}' for p in sig_order], fontsize=4.6)
ax.set_title('Aging gene set enrichment', loc='left')
fig.colorbar(points, ax=ax, fraction=0.035, pad=0.02, label='−log10 FDR')
top = TOP.head(10).iloc[::-1].reset_index(drop=True)
top['neglog'] = -np.log10(top['fdr_retained54'].clip(lower=1e-300))
ax = axes['g']
y = np.arange(len(top))
colors = [CATEGORY_COLORS.get(value, '#777777') for value in top['category']]
ax.hlines(y, 0, top['neglog'], color=colors, linewidth=1.3)
ax.scatter(top['neglog'], y, s=18, c=colors, edgecolors='#222222', linewidths=0.2)
ax.set_yticks(y)
ax.set_yticklabels([f'{row.new_P} {short_name(row.name_short, 16)} / {SET_LABEL.get(row.aging_set, row.aging_set)}' for row in top.itertuples()], fontsize=4.2)
ax.set_xlabel('−log10 FDR')
ax.set_title('Top aging linked programs', loc='left')
ax.spines[['top', 'right']].set_visible(False)
relations = SIG.sort_values('fdr_retained54').head(12).copy()
sets = [value for value in SET_ORDER if value in set(relations['aging_set'])]
programs = relations['new_P'].drop_duplicates().head(9).tolist()
relations = relations[relations['new_P'].isin(programs)]
set_y = dict(zip(sets, np.linspace(0.88, 0.12, len(sets))))
program_y = dict(zip(programs, np.linspace(0.88, 0.12, len(programs))))
ax = axes['h']
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
for row in relations.itertuples():
    ax.plot([0.26, 0.72], [set_y[row.aging_set], program_y[row.new_P]], color=CATEGORY_COLORS.get(row.category, '#777777'), linewidth=0.4 + min(float(row.odds_ratio), 10) / 5, alpha=0.45)
for (aging_set, y_value) in set_y.items():
    ax.text(0.23, y_value, SET_LABEL.get(aging_set, aging_set), ha='right', va='center', fontsize=4.5)
for (program, y_value) in program_y.items():
    ax.text(0.75, y_value, program, ha='left', va='center', fontsize=4.8)
ax.set_title('Shared aging relationships', loc='left')
summary_class = SUMMARY.set_index('new_P')['dominant_class'].to_dict()
donut_data = SIG.copy()
donut_data['dominant_class'] = donut_data['new_P'].map(summary_class).fillna(donut_data['dominant_class'])
category_order = ['Accelerated aging proteome', 'Decelerated aging proteome', 'SASP/SenMayo', 'Cellular senescence', 'Aging curated']
classes = ['nonneuron', 'vascular', 'glia', 'exc', 'inh']
ax = axes['i']
ax.set_xlim(-0.6, 4.6)
ax.set_ylim(-1.0, 0.8)
ax.set_aspect('equal')
ax.axis('off')
for (x, category) in enumerate(category_order):
    subset = donut_data[donut_data['category'].eq(category)]
    counts = subset.groupby('dominant_class')['new_P'].nunique()
    total = int(counts.sum())
    start = 90
    for cell_class in classes:
        value = int(counts.get(cell_class, 0))
        if value == 0 or total == 0:
            continue
        angle = 360 * value / total
        ax.add_patch(Wedge((x, 0), 0.42, start, start + angle, width=0.19, facecolor=CLASS_COLORS[cell_class], edgecolor='white', linewidth=0.4))
        start += angle
    ax.text(x, 0, str(total), ha='center', va='center', fontsize=7, fontweight='bold')
    ax.text(x, -0.6, category.replace(' proteome', ''), ha='center', va='top', fontsize=4.2, rotation=25)
ax.set_title('Cellular identity', loc='left')
ax.legend(handles=[Patch(facecolor=CLASS_COLORS[value], label=value) for value in classes], frameon=False, fontsize=4, ncol=3, loc='upper center')
chosen = [('P54', 'LIU2026_CELL_ACCELERATED_AGING_PROTEINS'), ('P49', 'LIU2026_CELL_ACCELERATED_AGING_PROTEINS'), ('P8', 'LIU2026_CELL_DECELERATED_AGING_PROTEINS'), ('P53', 'LIU2026_CELL_ACCELERATED_AGING_PROTEINS')]
ax = axes['j']
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.set_title('Aging overlap genes', loc='left')
for (index, (program, aging_set)) in enumerate(chosen):
    row = SIG[SIG['new_P'].eq(program) & SIG['aging_set'].eq(aging_set)]
    if row.empty:
        continue
    row = row.iloc[0]
    y = 0.88 - index * 0.24
    genes = str(row['overlap_genes']).split(';')[:3]
    ax.text(0, y, f'{program} / {SET_LABEL[aging_set]}', fontsize=5, fontweight='bold', va='top')
    ax.text(0, y - 0.07, f"OR {row['odds_ratio']:.1f}; FDR {row['fdr_retained54']:.1e}", fontsize=4.2, va='top')
    ax.text(0, y - 0.14, ', '.join(genes), fontsize=4.2, va='top')
for (label, ax) in axes.items():
    ax.text(-0.08, 1.02, label, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')
OUTPUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT / 'Fig8f_j.pdf', bbox_inches='tight')
fig.savefig(OUTPUT / 'Fig8f_j.png', dpi=300, bbox_inches='tight')
plt.close(fig)
