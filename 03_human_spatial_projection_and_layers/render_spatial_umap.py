#!/usr/bin/env python3
import os
from pathlib import Path
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
EMBEDDING = ROOT / 'figures/fig2/_intermediate/spatial_umap_expr_harmony.csv'
OVERLAYS = ROOT / 'figures/fig2/_intermediate/spatial_umap_feat.csv'
PROGRAM_MAP = ROOT / 'results/crossregion_v1/program_renumber_map.tsv'
OUTPUT = ROOT / 'figures/fig3'
mapping = pd.read_csv(PROGRAM_MAP, sep='\t')
mapping = mapping[mapping['status'].eq('kept')].copy()
mapping['old_P'] = mapping['old_P'].astype(int)
mapping['new_P'] = mapping['new_P'].astype(int)
old_to_new = dict(zip(mapping['old_P'], mapping['new_P']))
if len(old_to_new) != 54 or sorted(old_to_new.values()) != list(range(1, 55)):
    raise ValueError()
embedding = pd.read_csv(EMBEDDING)
overlays = pd.read_csv(OVERLAYS)
display_old = [37, 34, 7]
display_columns = [f'program_{i}' for i in display_old]
required = {'bin', 'UMAP1', 'UMAP2', 'majorDomain'}
if not required.issubset(embedding.columns) or not {'bin', *display_columns}.issubset(overlays.columns):
    raise ValueError()
data = embedding.merge(overlays[['bin', *display_columns]], on='bin', how='inner', validate='one_to_one')
domain_order = ['ARACHNOID', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM']
domain_colors = {'ARACHNOID': '#8c8c8c', 'L1': '#d7b5d8', 'L2': '#c994c7', 'L3': '#df65b0', 'L4': '#e7298a', 'L5': '#ce1256', 'L6': '#980043', 'WM': '#4d4d4d'}
program_names = {37: 'Oligodendrocyte/myelin', 34: 'Activity dependent IEG', 7: 'Cation channel (interneuron)'}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 6, 'axes.titlesize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
(fig, axes) = plt.subplots(1, 4, figsize=(7.1, 1.95), dpi=300)
for domain in domain_order:
    part = data[data['majorDomain'].eq(domain)]
    axes[0].scatter(part['UMAP1'], part['UMAP2'], s=0.35, c=domain_colors[domain], linewidths=0, rasterized=True, label=domain)
axes[0].set_title('Tissue domain')
axes[0].legend(frameon=False, markerscale=5, fontsize=4.6, ncol=2, handletextpad=0.2, columnspacing=0.5)
for (ax, old_id, column) in zip(axes[1:], display_old, display_columns):
    values = data[column].to_numpy(dtype=float)
    (lo, hi) = np.nanpercentile(values, [2, 98])
    image = ax.scatter(data['UMAP1'], data['UMAP2'], s=0.35, c=values, cmap='viridis', vmin=lo, vmax=hi, linewidths=0, rasterized=True)
    ax.set_title(f'P{old_to_new[old_id]} {program_names[old_id]}')
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    colorbar.ax.tick_params(labelsize=4, length=1)
for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    for spine in ax.spines.values():
        spine.set_visible(False)
OUTPUT.mkdir(parents=True, exist_ok=True)
fig.tight_layout(pad=0.5, w_pad=0.5)
fig.savefig(OUTPUT / 'Fig3i.pdf', bbox_inches='tight')
fig.savefig(OUTPUT / 'Fig3i.png', dpi=300, bbox_inches='tight')
plt.close(fig)
