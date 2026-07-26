#!/usr/bin/env python3
import os
from pathlib import Path
import matplotlib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
RESULTS = ROOT / 'results/crossregion_v1'
FINAL = RESULTS / 'markcorr_v2/final'
INPUT_DIR = Path(os.environ['MARKCORR_INPUT_DIR'])
BETWEEN = RESULTS / 'markcorr_betweenchip_v1'
PROGRAM_MAP = RESULTS / 'program_renumber_map.tsv'
OUTPUT = ROOT / 'figures/fig4'
SUPPLEMENT = ROOT / 'figures/figS6'
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

def save_figure(fig, directory, stem):
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(directory / f'{stem}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
pairs = pd.read_csv(FINAL / 'cellprog_pairs_median_iqr.tsv', sep='\t')
pairs['old_P'] = old_program(pairs['B'])
pairs = pairs[pairs['old_P'].isin(old_to_new)].copy()
pairs['new_P'] = pairs['old_P'].map(old_to_new)
pairs['program'] = 'P' + pairs['new_P'].astype(str)
r25 = pairs[pairs['ring_um'].eq(25)].copy()
between = pd.read_csv(BETWEEN / 'betweenchip_cellprog_stouffer_q.tsv', sep='\t')
between['new_P'] = old_program(between['B_name'])
between = between[between['new_P'].between(1, 54)].copy()
between['headline'] = between['is_headline'].astype(str).str.lower().eq('true')
headline = between[['A_name', 'new_P', 'q_bh', 'headline']].rename(columns={'A_name': 'A'})
r25 = r25.merge(headline, on=['A', 'new_P'], how='left')
logo = pd.read_csv(INPUT_DIR / 'cellprog_logo.tsv', sep='\t')
logo['old_P'] = old_program(logo['B'])
logo = logo[logo['old_P'].isin(old_to_new)].copy()
logo['new_P'] = logo['old_P'].map(old_to_new)
logo['survives'] = logo['survives'].astype(str).str.lower().eq('true')
logo['label'] = [f'{cell} · P{program} {name_by_new[program]}' for (cell, program) in zip(logo['A'], logo['new_P'])]
forest_pairs = [('L2-L3 IT LINC00507', 41), ('L2-L3 IT LINC00507', 32), ('L2-L3 IT LINC00507', 33), ('L6 IT', 49), ('VIP', 30), ('L6 CT', 19), ('L6B', 28), ('OPC', 39), ('OLIGO', 12), ('OLIGO', 33), ('OLIGO', 23), ('OLIGO', 32), ('OLIGO', 41)]
forest = pd.DataFrame(forest_pairs, columns=['A', 'new_P'])
forest['old_P'] = forest['new_P'].map(new_to_old)
forest = forest.merge(r25, on=['A', 'new_P', 'old_P'], how='left', validate='one_to_one')
forest = forest.merge(logo[['A', 'old_P', 'survives']], on=['A', 'old_P'], how='left')
if forest['log2_median_g'].isna().any():
    raise ValueError()
forest['label'] = forest['A'] + ' · P' + forest['new_P'].astype(str) + ' ' + forest['B_name']
forest['lo'] = np.log2(forest['q1_g'].clip(lower=np.finfo(float).tiny))
forest['hi'] = np.log2(forest['q3_g'].clip(lower=np.finfo(float).tiny))
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 6, 'axes.titlesize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
(fig, ax) = plt.subplots(figsize=(4.8, 3.5), dpi=300)
y = np.arange(len(forest))
colors = np.where(forest['log2_median_g'].to_numpy() >= 0, '#b2182b', '#2166ac')
ax.hlines(y, forest['lo'], forest['hi'], color=colors, linewidth=1)
ax.scatter(forest['log2_median_g'], y, c=colors, s=18, edgecolors='white', linewidths=0.3)
surviving = forest['survives'].fillna(False).to_numpy()
ax.scatter(forest.loc[surviving, 'log2_median_g'], y[surviving], facecolors='none', edgecolors='#222222', s=34, linewidths=0.6)
ax.axvline(0, color='#888888', linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(forest['label'])
ax.invert_yaxis()
ax.set_xlabel('log2 median g at 25 µm')
ax.set_title('Selected local cell and program relationships', loc='left')
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, OUTPUT, 'Fig4b')
logo_display = logo[logo['orig_log2g'].notna() & logo['logo_log2g'].notna() & logo['orig_log2g'].abs().ge(0.32)].copy()
logo_display = logo_display.sort_values('orig_log2g', key=lambda values: values.abs(), ascending=False).head(10)
logo_display = logo_display.sort_values('orig_log2g')
(fig, ax) = plt.subplots(figsize=(5.2, 3.4), dpi=300)
y = np.arange(len(logo_display))
ax.hlines(y, logo_display['orig_log2g'], logo_display['logo_log2g'], color='#bdbdbd', linewidth=1)
ax.scatter(logo_display['orig_log2g'], y, color='#333333', s=18, label='full')
colors = np.where(logo_display['survives'], '#b2182b', '#e08214')
ax.scatter(logo_display['logo_log2g'], y, color=colors, marker='^', s=22, label='gene removed')
for (index, row) in enumerate(logo_display.itertuples()):
    ax.text(max(row.orig_log2g, row.logo_log2g) + 0.04, index, f'J={row.jaccard:.2f}', va='center', fontsize=5, color='#555555')
ax.axvline(0, color='#888888', linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(logo_display['label'])
ax.set_xlabel('log2 median g')
ax.set_title('Gene removal retention', loc='left')
ax.legend(frameon=False, fontsize=5, loc='lower right')
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, OUTPUT, 'Fig4d')
region_pairs = [('OLIGO', 32), ('OLIGO', 41), ('OLIGO', 33), ('OLIGO', 23), ('OLIGO', 12), ('OPC', 39), ('L6B', 28), ('L6 CT', 19), ('VIP', 30), ('L6 IT', 49)]
region_order = ['DLPFC', 'S1', 'AG', 'M1', 'V1', 'SMG']
byarea = pd.read_csv(FINAL / 'cellprog_byarea_median_iqr.tsv', sep='\t')
byarea['old_P'] = old_program(byarea['B'])
byarea = byarea[byarea['old_P'].isin(old_to_new)].copy()
byarea['new_P'] = byarea['old_P'].map(old_to_new)
selected_rows = []
for (cell, new_id) in region_pairs:
    part = byarea[byarea['A'].eq(cell) & byarea['new_P'].eq(new_id) & byarea['ring_um'].eq(25) & byarea['area'].isin(region_order)].copy()
    part['row'] = f"{cell} · P{new_id} {part['B_name']}"
    selected_rows.append(part)
region_data = pd.concat(selected_rows, ignore_index=True)
region_matrix = region_data.pivot(index='row', columns='area', values='log2_median_g').reindex(columns=region_order).reindex([f"{cell} · P{new_id} {region_data[region_data['A'].eq(cell) & region_data['new_P'].eq(new_id)]['B_name'].iloc[0]}" for (cell, new_id) in region_pairs])
(fig, ax) = plt.subplots(figsize=(5.2, 3.9), dpi=300)
limit = np.nanpercentile(np.abs(region_matrix.to_numpy()), 98)
image = ax.imshow(region_matrix, aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit)
ax.set_xticks(np.arange(len(region_order)))
ax.set_xticklabels(region_order, rotation=45, ha='right')
ax.set_yticks(np.arange(len(region_matrix.index)))
ax.set_yticklabels(region_matrix.index)
ax.set_title('Regional variation in selected local relationships', loc='left')
fig.colorbar(image, ax=ax, label='log2 median g')
fig.tight_layout()
save_figure(fig, OUTPUT, 'Fig4k')
cell_order = ['L2-L3 IT LINC00507', 'L4-L5 IT RORB', 'L3-L4 IT RORB', 'SST', 'VIP', 'L6 IT', 'L6 CAR3', 'OPC', 'L6 CT', 'L6B', 'ET', 'CHANDELIER', 'PVALB', 'NDNF', 'PAX6', 'NP', 'LAMP5', 'AST', 'MICRO', 'ENDO', 'VLMC', 'OLIGO']
program_order = [42, 35, 27, 20, 30, 25, 7, 46, 26, 43, 18, 17, 29, 52, 37, 19, 28, 44, 1, 24, 11, 10, 22, 15, 6, 9, 13, 16, 3, 5, 2, 38, 50, 36, 45, 54, 51, 53, 39, 34, 47, 14, 40, 48, 31, 49, 4, 21, 8, 12, 23, 33, 32, 41]
full_matrix = r25.pivot(index='A', columns='new_P', values='log2_median_g').reindex(index=cell_order, columns=program_order)
if full_matrix.shape != (22, 54):
    raise ValueError()
(fig, ax) = plt.subplots(figsize=(7.1, 3.4), dpi=300)
limit = np.nanpercentile(np.abs(full_matrix.to_numpy()), 98)
image = ax.imshow(full_matrix, aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit)
ax.set_xticks(np.arange(54))
ax.set_xticklabels([f'P{i}' for i in program_order], rotation=90, fontsize=4.2)
ax.set_yticks(np.arange(22))
ax.set_yticklabels(cell_order, fontsize=4.8)
fig.colorbar(image, ax=ax, label='log2 median g')
fig.tight_layout()
save_figure(fig, SUPPLEMENT, 'FigS6a')
target_chip = 'B02111F5'
meta_parts = []
meta_file = pq.ParquetFile(RESULTS / 'spatial_bin50_meta.parquet')
for batch in meta_file.iter_batches(batch_size=250000, columns=['bin', 'chip', 'x', 'y']):
    part = batch.to_pandas()
    part = part[part['chip'].eq(target_chip)]
    if not part.empty:
        meta_parts.append(part[['bin', 'x', 'y']])
coordinates = pd.concat(meta_parts, ignore_index=True)
target_bins = set(coordinates['bin'].astype(str))
score_parts = []
score_file = pq.ParquetFile(RESULTS / 'spatial_bin_program_score.parquet')
for batch in score_file.iter_batches(batch_size=250000, columns=['bin', 'program_37']):
    part = batch.to_pandas()
    part['bin'] = part['bin'].astype(str)
    part = part[part['bin'].isin(target_bins)]
    if not part.empty:
        score_parts.append(part)
scores = pd.concat(score_parts, ignore_index=True)
weight_parts = []
weight_file = pq.ParquetFile(RESULTS / 'spatial_bin_rctd_weights.parquet')
for batch in weight_file.iter_batches(batch_size=250000, columns=['bin', 'OLIGO', 'rctd_chip']):
    part = batch.to_pandas()
    part = part[part['rctd_chip'].eq(target_chip)]
    if not part.empty:
        weight_parts.append(part[['bin', 'OLIGO']])
weights = pd.concat(weight_parts, ignore_index=True)
maps = coordinates.merge(scores, on='bin', how='inner').merge(weights, on='bin', how='inner')
if len(maps) > 120000:
    maps = maps.sample(120000, random_state=42)
(fig, axes) = plt.subplots(1, 2, figsize=(6.8, 3.2), dpi=300)
image1 = axes[0].scatter(maps['x'], maps['y'], c=maps['program_37'], s=0.2, cmap='viridis', vmin=-0.5, vmax=1.5, linewidths=0, rasterized=True)
image2 = axes[1].scatter(maps['x'], maps['y'], c=maps['OLIGO'], s=0.2, cmap='magma', vmin=0, vmax=0.5, linewidths=0, rasterized=True)
axes[0].set_title('P33 Oligodendrocyte/myelin')
axes[1].set_title('OLIGO weight')
for ax in axes:
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_axis_off()
fig.colorbar(image1, ax=axes[0], fraction=0.035, pad=0.01)
fig.colorbar(image2, ax=axes[1], fraction=0.035, pad=0.01)
fig.suptitle(target_chip, y=0.98)
fig.tight_layout()
save_figure(fig, SUPPLEMENT, 'FigS6c')
profile_programs = [(33, 'OLIGO'), (51, 'ENDO'), (5, 'L2-L3 IT LINC00507')]
(fig, axes) = plt.subplots(1, 3, figsize=(7.1, 3.4), dpi=300)
for (ax, (new_id, peak_cell)) in zip(axes, profile_programs):
    profile = r25[r25['new_P'].eq(new_id)].sort_values('log2_median_g')
    colors = np.where(profile['A'].eq(peak_cell), '#b2182b', '#7f8c8d')
    ax.barh(np.arange(len(profile)), profile['log2_median_g'], color=colors)
    ax.set_yticks(np.arange(len(profile)))
    ax.set_yticklabels(profile['A'], fontsize=4.1)
    ax.axvline(0, color='#888888', linewidth=0.4)
    ax.set_title(f"P{new_id} {profile['B_name'].iloc[0]}")
    ax.set_xlabel('log2 median g')
    ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, SUPPLEMENT, 'FigS6d')
distance_programs = [32, 33, 41]
(fig, ax) = plt.subplots(figsize=(4.8, 3.4), dpi=300)
colors = ['#4e79a7', '#59a14f', '#e15759']
for (new_id, color) in zip(distance_programs, colors):
    curve = pairs[pairs['A'].eq('OLIGO') & pairs['new_P'].eq(new_id)].sort_values('ring_um')
    y = curve['log2_median_g'].to_numpy()
    lo = np.log2(curve['q1_g'].clip(lower=np.finfo(float).tiny).to_numpy())
    hi = np.log2(curve['q3_g'].clip(lower=np.finfo(float).tiny).to_numpy())
    ax.fill_between(curve['ring_um'], lo, hi, color=color, alpha=0.15)
    ax.plot(curve['ring_um'], y, color=color, linewidth=1.2, label=f"P{new_id} {curve['B_name'].iloc[0]}")
ax.axhline(0, color='#888888', linewidth=0.4)
ax.set_xlabel('distance (µm)')
ax.set_ylabel('log2 median g')
ax.set_title('Spatial decay from OLIGO neighborhoods', loc='left')
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
save_figure(fig, SUPPLEMENT, 'FigS6e')
