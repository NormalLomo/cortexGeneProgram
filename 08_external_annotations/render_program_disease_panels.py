#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os, json
import numpy as np
import pandas as pd
import matplotlib
from pathlib import Path
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Liberation Sans', 'Nimbus Sans', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 6, 'axes.titlesize': 7, 'axes.labelsize': 6.2, 'xtick.labelsize': 5.4, 'ytick.labelsize': 5.4, 'legend.fontsize': 5.4, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 1.8, 'ytick.major.size': 1.8, 'axes.edgecolor': '#333333', 'pdf.use14corefonts': False})
FS_MIN = 5.0
ROOT = str(PROJECT_ROOT)
PD_DIR = f'{ROOT}/results/crossregion_v1/program_disease'
LOAD = f'{ROOT}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv'
NAMES = f'{ROOT}/results/crossregion_v1/program_names.tsv'
RENUM = f'{ROOT}/results/crossregion_v1/program_renumber_map.tsv'
RZ = f'{ROOT}/results/crossregion_v1/program_region_zscore.tsv'
BB = f'{ROOT}/data/brainbase/disease_gene_associations.txt'
FIGDIR = f'{ROOT}/figures/fig8'
os.makedirs(FIGDIR, exist_ok=True)
PRIMARY_N = 150
FDR_SIG = 0.05
renum_df = _annotation_table(pd.read_csv(RENUM, sep='\t'))
renum_df = renum_df[renum_df['status'].astype(str).str.lower().eq('kept') & renum_df['new_P'].notna()].copy()
renum_df['old_P'] = renum_df['old_P'].astype(int)
renum_df['new_P'] = renum_df['new_P'].astype(int)
renum_df = renum_df.sort_values('new_P')
old_to_new = dict(zip(renum_df['old_P'], renum_df['new_P']))
meta = json.load(open(f'{PD_DIR}/meta.json'))
diseases_all = meta['diseases']
cat_of = meta['category']
dn = meta['disease_n_genes']
long_df = _annotation_table(pd.read_csv(f'{PD_DIR}/enrichment_long_N{PRIMARY_N}.tsv', sep='\t'))
program_ids = set(long_df['program'].astype(int))
retained_old = set(renum_df['old_P'])
if program_ids != set(range(1, 55)):
    raise ValueError()
fdr_mat = long_df.pivot(index='program', columns='disease', values='fdr')
or_mat = long_df.pivot(index='program', columns='disease', values='odds_ratio')
_sig54 = long_df[long_df['fdr'] < FDR_SIG]
catsum = _sig54.groupby('category')['program'].nunique().rename('n_programs').to_frame()
dom = _annotation_table(pd.read_csv(f'{PD_DIR}/program_dom_subclass.tsv', sep='\t'))
dom['program'] = dom['program'].astype(int)
if retained_old.issubset(set(dom['program'])):
    dom = dom[dom['program'].isin(retained_old)].copy()
    dom['program'] = dom['program'].map(old_to_new)
elif set(dom['program']) != set(range(1, 55)):
    raise ValueError()
dom = dom.set_index('program')
names_df = _annotation_table(pd.read_csv(NAMES, sep='\t'))
source_col = 'cnmf_component' if 'cnmf_component' in names_df.columns else 'program'
names_df = names_df.merge(renum_df[['old_P', 'new_P']], left_on=source_col, right_on='old_P', suffixes=('_source', ''))
names_df = names_df.set_index('new_P')
name_short = names_df['name_short'].to_dict()
star_map = {idx: '*' if row['confidence'] == 'class_b' else '' for (idx, row) in names_df.iterrows()}

def disp_label(program):
    return f"P{program}{star_map.get(program, '')} {name_short.get(program, '')}"
or_mat = or_mat[fdr_mat.columns]
CAT_ORDER = ['Neurodegenerative', 'Psychiatric/Neurodev', 'Tumor', 'Developmental', 'Other/Vascular']
CAT_COL = {'Neurodegenerative': '#C0392B', 'Psychiatric/Neurodev': '#2471A3', 'Tumor': '#7D3C98', 'Developmental': '#1E8449', 'Other/Vascular': '#B9770E'}

def lineage(sc):
    if sc in ('MICRO',):
        return 'Microglia'
    if sc in ('AST',):
        return 'Astrocyte'
    if sc in ('OLIGO', 'OPC'):
        return 'Oligo/OPC'
    if sc in ('ENDO', 'VLMC'):
        return 'Vascular'
    if sc in ('PVALB', 'SST', 'VIP', 'LAMP5', 'NDNF', 'PAX6', 'CHANDELIER'):
        return 'Inhibitory'
    return 'Excitatory'
LIN_COL = {'Microglia': '#8E44AD', 'Astrocyte': '#16A085', 'Oligo/OPC': '#D68910', 'Vascular': '#CB4335', 'Inhibitory': '#2E86C1', 'Excitatory': '#34495E'}
neglog = -np.log10(fdr_mat.clip(lower=1e-300))
sig_progs = sorted(set(long_df.loc[long_df['fdr'] < FDR_SIG, 'program']))
A_rows = sig_progs
sub_nl = neglog.loc[A_rows]
sub_or = or_mat.loc[A_rows]
if len(A_rows) > 2:
    Z = linkage(pdist(sub_nl.values, metric='euclidean'), method='average')
    order = [A_rows[i] for i in leaves_list(Z)]
else:
    order = A_rows
sub_nl = sub_nl.loc[order]
sub_or = sub_or.loc[order]
sig_dis = [d for d in fdr_mat.columns if (fdr_mat[d] < FDR_SIG).any()]
sig_dis = sorted(sig_dis, key=lambda d: (CAT_ORDER.index(cat_of[d]), -dn[d]))
sub_nl = sub_nl[sig_dis]
sub_or = sub_or[sig_dis]
prog_labels = [disp_label(p) for p in order]
best = long_df[long_df['fdr'] < FDR_SIG].sort_values('fdr').groupby('program').first().reset_index()
best['neglogfdr'] = -np.log10(best['fdr'].clip(lower=1e-300))
best = best.sort_values('neglogfdr', ascending=False).head(15)
best['dom'] = best['program'].map(dom['dom_subclass'])
best['lineage'] = best['dom'].map(lineage)
EXEMPLARS = [(36, "Alzheimer's Disease"), (1, 'Autism Spectrum Disorder'), (1, 'Epilepsy'), (53, 'Glioma')]
load = _annotation_table(pd.read_csv(LOAD, sep='\t', index_col=0))
load.index = [int(i) for i in load.index]
load = load.loc[renum_df['old_P']].copy()
load.index = renum_df['new_P']
bb = _annotation_table(pd.read_csv(BB, sep='\t'))
bb.columns = [c.strip() for c in bb.columns]
gcol = 'Gene symbol' if 'Gene symbol' in bb.columns else 'Gene'
mrna = bb[bb['Type'] == 'mRNA']
dgenes = {d: set(mrna.loc[mrna['Disease'] == d, gcol].astype(str)) for (_, d) in EXEMPLARS}
rz = _annotation_table(pd.read_csv(RZ, sep='\t', index_col=0))
rz.columns = [int(c) for c in rz.columns]
rz = rz[[old for old in renum_df['old_P'] if old in rz.columns]].copy()
rz.columns = [old_to_new[int(column)] for column in rz.columns]
top_for_region = list(best['program'].head(12))
top_for_region = [p for p in top_for_region if p in rz.columns]
rz_sub = rz[top_for_region].T
REGION_ORDER = ['FPPFC', 'VLPFC', 'DLPFC', 'ACC', 'M1', 'S1', 'S1E', 'PoCG', 'SMG', 'AG', 'SPL', 'STG', 'ITG', 'V1']
rz_sub = rz_sub[[r for r in REGION_ORDER if r in rz_sub.columns]]
FW = 180.0 / 25.4
FH = 224.0 / 25.4
fig = plt.figure(figsize=(FW, FH), dpi=300)
gs = GridSpec(2, 2, figure=fig, height_ratios=[3.05, 1.0], width_ratios=[1.18, 1.0], hspace=0.3, wspace=0.34, left=0.115, right=0.965, top=0.952, bottom=0.058)
gs_right = gs[0, 1].subgridspec(2, 1, height_ratios=[1.16, 1.0], hspace=0.42)
gs_bot = gs[1, :].subgridspec(1, 2, width_ratios=[2.95, 1.0], wspace=0.3)

def panel_tag(ax, s, dx=-0.052, dy=1.022):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='right')
axa = fig.add_subplot(gs[0, 0])
(nR, nC) = sub_nl.shape
cmap_a = LinearSegmentedColormap.from_list('fdr', ['#F2F2F2', '#FBD9B6', '#F0894E', '#C0392B', '#7B1B12'])
vmax = float(np.nanpercentile(sub_nl.values[sub_nl.values > 0], 98)) if (sub_nl.values > 0).any() else 3
vmax = max(vmax, 4)
norm_a = Normalize(vmin=0, vmax=vmax)
or_vals = sub_or.values
(smin, smax) = (1.0, np.nanpercentile(or_vals, 97))

def or_to_size(o):
    o = np.clip(o, 1.0, smax)
    return 5 + 46 * (o - smin) / (smax - smin + 1e-09)
(xs, ys, cs, ss) = ([], [], [], [])
for i in range(nR):
    for j in range(nC):
        v = sub_nl.values[i, j]
        xs.append(j)
        ys.append(i)
        cs.append(v if v > 0 else 0)
        o = sub_or.values[i, j]
        ss.append(or_to_size(o) if v > -np.log10(FDR_SIG) else 0.0)
axa.set_xlim(-0.5, nC - 0.5)
axa.set_ylim(-0.5, nR - 0.5)
axa.set_aspect('equal')
for j in range(nC + 1):
    axa.axvline(j - 0.5, color='#E8E8E8', lw=0.3, zorder=0)
for i in range(nR + 1):
    axa.axhline(i - 0.5, color='#E8E8E8', lw=0.3, zorder=0)
sc = axa.scatter(xs, ys, c=cs, s=ss, cmap=cmap_a, norm=norm_a, edgecolors='#333333', linewidths=0.18, zorder=3)
axa.set_xticks(range(nC))
axa.set_xticklabels(sig_dis, rotation=55, ha='right', va='top', fontsize=FS_MIN)
axa.set_yticks(range(nR))
axa.set_yticklabels(prog_labels, fontsize=FS_MIN)
axa.invert_yaxis()
axa.tick_params(length=1.2)
for (j, d) in enumerate(sig_dis):
    axa.add_patch(plt.Rectangle((j - 0.5, -1.15), 1, 0.7, color=CAT_COL[cat_of[d]], clip_on=False, lw=0))
axa.set_ylim(nR - 0.5, -1.3)
axa.set_title('Program × brain-disease enrichment\n(BrainBase mRNA, top-150 loading genes)', fontsize=7, pad=8)
panel_tag(axa, 'a', dx=-0.3)
cax = axa.inset_axes([1.045, 0.74, 0.03, 0.24])
cb = fig.colorbar(sc, cax=cax)
cb.set_label('-log10 FDR', fontsize=FS_MIN, labelpad=2)
cb.ax.tick_params(labelsize=FS_MIN, length=1.2)
leg_or = [2, 4, 8]
hs = [axa.scatter([], [], s=or_to_size(o), color='#888888', edgecolors='#333333', linewidths=0.18, label=f'{o}') for o in leg_or]
size_leg = axa.legend(handles=hs, title='odds ratio', loc='upper left', bbox_to_anchor=(1.04, 0.66), frameon=False, labelspacing=0.7, handletextpad=0.3, fontsize=FS_MIN, title_fontsize=FS_MIN + 0.4)
axa.add_artist(size_leg)
catleg = [Patch(fc=CAT_COL[c], ec='none', label=c) for c in CAT_ORDER]
axa.legend(handles=catleg, loc='upper left', bbox_to_anchor=(1.04, 0.38), frameon=False, fontsize=FS_MIN, title='disease category', title_fontsize=FS_MIN + 0.4, handlelength=1.0, handleheight=1.0)
axb = fig.add_subplot(gs_right[0, 0])
b = best.iloc[::-1].reset_index(drop=True)
yb = np.arange(len(b))
cols = [LIN_COL[l] for l in b['lineage']]
axb.hlines(yb, 0, b['neglogfdr'], color=cols, lw=1.4, zorder=1)
axb.scatter(b['neglogfdr'], yb, s=22, c=cols, edgecolors='#222222', linewidths=0.25, zorder=2)
axb.set_yticks(yb)
axb.set_yticklabels([disp_label(p) for p in b['program']], fontsize=FS_MIN)
for (i, row) in b.iterrows():
    axb.text(row['neglogfdr'] + 0.4, i, f"{row['disease']}", va='center', ha='left', fontsize=FS_MIN, color='#444444')
axb.set_xlabel('-log10 FDR (best disease hit)', fontsize=FS_MIN + 0.8)
axb.set_xlim(0, b['neglogfdr'].max() * 1.55)
axb.axvline(-np.log10(FDR_SIG), color='#999999', ls='--', lw=0.5)
axb.set_title('Top disease-enriched programs', fontsize=7, pad=4)
axb.spines[['top', 'right']].set_visible(False)
linleg = [Patch(fc=LIN_COL[l], ec='none', label=l) for l in ['Excitatory', 'Inhibitory', 'Microglia', 'Astrocyte', 'Oligo/OPC', 'Vascular']]
axb.legend(handles=linleg, loc='lower right', frameon=False, fontsize=FS_MIN, title='dominant subclass', title_fontsize=FS_MIN + 0.4, ncol=1, handlelength=1.0)
panel_tag(axb, 'b', dx=-0.14)
axc = fig.add_subplot(gs_right[1, 0])
M = rz_sub.values
vlim = np.nanpercentile(np.abs(M), 98)
cmap_c = LinearSegmentedColormap.from_list('rz', ['#2166AC', '#92C5DE', '#F7F7F7', '#F4A582', '#B2182B'])
im = axc.imshow(M, aspect='equal', cmap=cmap_c, vmin=-vlim, vmax=vlim)
axc.set_xticks(range(rz_sub.shape[1]))
axc.set_xticklabels(rz_sub.columns, rotation=55, ha='right', fontsize=FS_MIN)
axc.set_yticks(range(rz_sub.shape[0]))
axc.set_yticklabels([disp_label(p) for p in rz_sub.index], fontsize=FS_MIN)
axc.set_title('Region enrichment of disease programs', fontsize=7, pad=4)
axc.tick_params(length=1.2)
ccx = axc.inset_axes([1.03, 0.18, 0.045, 0.64])
cbc = fig.colorbar(im, cax=ccx)
cbc.set_label('region z', fontsize=FS_MIN, labelpad=2)
cbc.ax.tick_params(labelsize=FS_MIN, length=1.2)
axc.text(1.0, -0.3, 'anterior → posterior', transform=axc.transAxes, fontsize=FS_MIN, ha='right', color='#666666')
panel_tag(axc, 'c', dx=-0.14)
gse = gs_bot[0, 0].subgridspec(1, len(EXEMPLARS), wspace=1.05)
for (k, (program, dis)) in enumerate(EXEMPLARS):
    axe = fig.add_subplot(gse[0, k])
    row = load.loc[program]
    top = row.sort_values(ascending=False).head(PRIMARY_N)
    le = [g for g in top.index if g in dgenes[dis]]
    le = sorted(le, key=lambda g: row[g], reverse=True)[:10]
    vals = [row[g] for g in le][::-1]
    gg = le[::-1]
    col = LIN_COL[lineage(dom.loc[program, 'dom_subclass'])]
    axe.barh(range(len(gg)), vals, color=col, edgecolor='#222', linewidth=0.2, height=0.74)
    axe.set_yticks(range(len(gg)))
    axe.set_yticklabels(gg, fontsize=FS_MIN, style='italic')
    axe.set_xlabel('loading', fontsize=FS_MIN, labelpad=1)
    short = {"Alzheimer's Disease": 'AD', 'Autism Spectrum Disorder': 'ASD', 'Epilepsy': 'Epilepsy', 'Glioma': 'Glioma'}[dis]
    axe.set_title(f'P{program}→{short}', fontsize=FS_MIN + 0.5, pad=2)
    axe.tick_params(length=1.0, labelsize=FS_MIN)
    axe.spines[['top', 'right']].set_visible(False)
    if k == 0:
        panel_tag(axe, 'e', dx=-0.46, dy=1.05)
axg = fig.add_subplot(gs_bot[0, 1])
cs2 = catsum.reindex(CAT_ORDER).fillna(0)
yc = np.arange(len(CAT_ORDER))
axg.barh(yc, cs2['n_programs'].values, color=[CAT_COL[c] for c in CAT_ORDER], edgecolor='#222', linewidth=0.2, height=0.7)
axg.set_yticks(yc)
axg.set_yticklabels(CAT_ORDER, fontsize=FS_MIN)
for (i, v) in enumerate(cs2['n_programs'].values):
    axg.text(v + 0.3, i, int(v), va='center', fontsize=FS_MIN)
axg.set_xlabel('n programs enriched', fontsize=FS_MIN + 0.8, labelpad=1)
axg.set_xlim(0, cs2['n_programs'].max() * 1.28)
axg.invert_yaxis()
axg.tick_params(length=1.0, labelsize=FS_MIN)
axg.spines[['top', 'right']].set_visible(False)
axg.set_title('Programs per disease class', fontsize=7, pad=4)
panel_tag(axg, 'd', dx=-0.5)
out_pdf = f'{FIGDIR}/fig8_program_disease.pdf'
out_png = f'{FIGDIR}/fig8_program_disease.png'
full_pdf = f'{FIGDIR}/fig8_full.pdf'
full_png = f'{FIGDIR}/fig8_full.png'
fig.savefig(out_pdf, dpi=300)
fig.savefig(out_png, dpi=320)
fig.savefig(full_pdf, dpi=300)
fig.savefig(full_png, dpi=320)
