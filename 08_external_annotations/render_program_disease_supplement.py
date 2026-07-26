#!/usr/bin/env python
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
import networkx as nx
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 6, 'axes.titlesize': 7, 'axes.labelsize': 6.2, 'xtick.labelsize': 5.4, 'ytick.labelsize': 5.4, 'legend.fontsize': 5.4, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 1.8, 'ytick.major.size': 1.8, 'axes.edgecolor': '#333333', 'pdf.use14corefonts': False})
FS_MIN = 5.0
ROOT = str(PROJECT_ROOT)
PD_DIR = f'{ROOT}/results/crossregion_v1/program_disease'
EDDIR = f'{ROOT}/figures/extended'
os.makedirs(EDDIR, exist_ok=True)
PRIMARY_N = 150
FDR_SIG = 0.05
meta = json.load(open(f'{PD_DIR}/meta.json'))
cat_of = meta['category']
long_df = pd.read_csv(f'{PD_DIR}/enrichment_long_N{PRIMARY_N}.tsv', sep='\t')
if set(long_df['program'].astype(int)) != set(range(1, 55)):
    raise ValueError()
dom = pd.read_csv(f'{PD_DIR}/program_dom_subclass.tsv', sep='\t')

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    ranked = np.empty_like(p)
    prev = 1.0
    n = len(p)
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = p[idx] * n / rank
        prev = min(prev, val)
        ranked[idx] = min(prev, 1.0)
    return ranked
_RMAP_PATH = f'{ROOT}/results/crossregion_v1/program_renumber_map.tsv'
_rmap = pd.read_csv(_RMAP_PATH, sep='\t')
_rmap = _rmap[_rmap['status'].astype(str).str.lower().eq('kept') & _rmap['new_P'].notna()].copy()
_rmap['old_P'] = _rmap['old_P'].astype(int)
_rmap['new_P'] = _rmap['new_P'].astype(int)
_rmap = _rmap.sort_values('new_P')
_old2new = dict(zip(_rmap['old_P'], _rmap['new_P']))
dom['program'] = dom['program'].astype(int)
if set(dom['program']) == set(_rmap['old_P']):
    dom = dom[dom['program'].isin(_old2new)].copy()
    dom['program'] = dom['program'].map(_old2new)
elif set(dom['program']) != set(range(1, 55)):
    raise ValueError()
dom = dom.set_index('program')
sens_rows = []
for _top_n in (100, 150, 200):
    _df = pd.read_csv(f'{PD_DIR}/enrichment_long_N{_top_n}.tsv', sep='\t')
    if set(_df['program'].astype(int)) != set(range(1, 55)):
        raise ValueError()
    _sig = _df[_df['fdr'] < FDR_SIG]
    sens_rows.append({'topN': _top_n, 'n_sig_pairs': int(len(_sig)), 'n_programs_hit': int(_sig['program'].nunique()), 'n_diseases_hit': int(_sig['disease'].nunique())})
sens = pd.DataFrame(sens_rows)
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
FW = 180.0 / 25.4
FH = 150.0 / 25.4
fig = plt.figure(figsize=(FW, FH), dpi=300)
gs = GridSpec(1, 2, figure=fig, width_ratios=[2.35, 1.0], wspace=0.18, left=0.045, right=0.965, top=0.918, bottom=0.105)

def panel_tag(ax, s, dx=-0.02, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='right')
axd = fig.add_subplot(gs[0, 0])
sigpairs = long_df[long_df['fdr'] < FDR_SIG].copy()
G = nx.Graph()
progs_d = sorted(sigpairs['program'].unique())
dis_d = sorted(sigpairs['disease'].unique(), key=lambda d: (CAT_ORDER.index(cat_of[d]), d))
for p in progs_d:
    G.add_node(('P', p))
for d in dis_d:
    G.add_node(('D', d))
for (_, r) in sigpairs.iterrows():
    G.add_edge(('P', r['program']), ('D', r['disease']), w=-np.log10(max(r['fdr'], 1e-300)))
posd = {}
for (i, p) in enumerate(progs_d):
    posd['P', p] = (0.0, 1 - i / max(len(progs_d) - 1, 1))
for (i, d) in enumerate(dis_d):
    posd['D', d] = (1.0, 1 - i / max(len(dis_d) - 1, 1))
ews = np.array([G[u][v]['w'] for (u, v) in G.edges()])
ew_n = 0.2 + 1.6 * (ews - ews.min()) / (np.ptp(ews) + 1e-09)
for ((u, v), w) in zip(G.edges(), ew_n):
    (x0, y0) = posd[u]
    (x1, y1) = posd[v]
    dnode = v if v[0] == 'D' else u
    axd.plot([x0, x1], [y0, y1], lw=w, color=CAT_COL[cat_of[dnode[1]]], alpha=0.45, zorder=1)
for p in progs_d:
    (x, y) = posd['P', p]
    axd.scatter([x], [y], s=16, color=LIN_COL[lineage(dom.loc[p, 'dom_subclass'])], edgecolors='#222', linewidths=0.2, zorder=3)
    axd.text(x - 0.025, y, f'P{p}', ha='right', va='center', fontsize=FS_MIN, color='#333')
for d in dis_d:
    (x, y) = posd['D', d]
    axd.scatter([x], [y], s=20, marker='s', color=CAT_COL[cat_of[d]], edgecolors='#222', linewidths=0.2, zorder=3)
    axd.text(x + 0.025, y, d, ha='left', va='center', fontsize=FS_MIN, color='#333')
axd.set_xlim(-0.32, 1.42)
axd.set_ylim(-0.18, 1.1)
axd.axis('off')
axd.text(0.0, 1.07, 'programs', ha='center', fontsize=FS_MIN + 0.6, fontweight='bold')
axd.text(1.0, 1.07, 'diseases', ha='center', fontsize=FS_MIN + 0.6, fontweight='bold')
axd.set_title('Program–disease significant network (FDR < 0.05)', fontsize=7, pad=4)
linleg = [Line2D([], [], marker='o', ls='', mfc=LIN_COL[l], mec='#222', mew=0.2, ms=4, label=l) for l in ['Excitatory', 'Inhibitory', 'Microglia', 'Astrocyte', 'Oligo/OPC', 'Vascular']]
l1 = axd.legend(handles=linleg, loc='upper left', bbox_to_anchor=(-0.02, 0.045), frameon=False, fontsize=FS_MIN, title='program dominant subclass', title_fontsize=FS_MIN + 0.4, handletextpad=0.3, borderaxespad=0.0, ncol=3, columnspacing=1.1)
axd.add_artist(l1)
catleg = [Line2D([], [], marker='s', ls='', mfc=CAT_COL[c], mec='#222', mew=0.2, ms=4, label=c) for c in CAT_ORDER]
axd.legend(handles=catleg, loc='upper left', bbox_to_anchor=(0.52, 0.045), frameon=False, fontsize=FS_MIN, title='disease category', title_fontsize=FS_MIN + 0.4, handletextpad=0.3, borderaxespad=0.0, ncol=2, columnspacing=1.1)
panel_tag(axd, 'a', dx=0.0)
axf = fig.add_subplot(gs[0, 1])
axf.plot(sens['topN'], sens['n_sig_pairs'], '-o', color='#C0392B', lw=1.2, ms=4, mec='#222', mew=0.3)
axf.fill_between(sens['topN'], sens['n_sig_pairs'], color='#C0392B', alpha=0.12)
for (_, r) in sens.iterrows():
    axf.text(r['topN'], r['n_sig_pairs'] + 2, int(r['n_sig_pairs']), ha='center', fontsize=FS_MIN)
axf.set_xticks(sens['topN'])
axf.set_xlabel('top-N loading genes', fontsize=FS_MIN + 1.0, labelpad=2)
axf.set_ylabel('n significant program-disease pairs', fontsize=FS_MIN + 1.0, labelpad=2)
axf.set_ylim(0, sens['n_sig_pairs'].max() * 1.3)
axf.tick_params(length=1.4, labelsize=FS_MIN + 0.5)
axf.spines[['top', 'right']].set_visible(False)
axf.set_title('N-sensitivity of enrichment', fontsize=7, pad=4)
axf.axvline(150, color='#999999', ls='--', lw=0.5)
axf.text(150, axf.get_ylim()[1] * 0.97, 'primary N=150', ha='center', va='top', fontsize=FS_MIN, color='#666666', rotation=0)
panel_tag(axf, 'b', dx=-0.18)
ed_pdf = f'{EDDIR}/ed_program_disease_supp.pdf'
ed_png = f'{EDDIR}/ed_program_disease_supp.png'
fig.savefig(ed_pdf, dpi=300)
fig.savefig(ed_png, dpi=320)
