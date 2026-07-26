#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os
import re
import numpy as np
import pandas as pd
import matplotlib
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as _mpl_font
_mpl_font.rcParams['font.family'] = 'sans-serif'
_mpl_font.rcParams['font.sans-serif'] = ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans']
_mpl_font.rcParams['pdf.fonttype'] = 42
_mpl_font.rcParams['ps.fonttype'] = 42
_mpl_font.rcParams['svg.fonttype'] = 'none'
_mpl_font.rcParams['mathtext.fontset'] = 'dejavusans'
_mpl_font.rcParams['mathtext.default'] = 'regular'
from matplotlib.patches import Rectangle
import matplotlib.font_manager as fm
ROOT = PROJECT_ROOT
CR = os.path.join(ROOT, 'results', 'crossregion_v1')
FIGDIR = os.path.join(ROOT, 'figures', 'extended')
os.makedirs(FIGDIR, exist_ok=True)
names = _annotation_table(pd.read_csv(os.path.join(CR, 'program_names.tsv'), sep='\t'))
SUPP_DIR = os.path.join(ROOT, 'results', 'crossregion_v1')
ann = _annotation_table(pd.read_csv(os.path.join(SUPP_DIR, 'TableS1_program_annotation.tsv'), sep='\t'))
gobp = _annotation_table(pd.read_csv(os.path.join(CR, 'program_annotation_gobp.tsv'), sep='\t'))
names_kept = names[names['new_P'] != 'EXCLUDED'].copy().reset_index(drop=True)
excluded_old = names[names['new_P'] == 'EXCLUDED']['cnmf_component'].tolist()

def old_pid(x):
    return int(re.sub('[^0-9]', '', str(x)))
gobp['cnmf_component'] = gobp['program'].astype(int)
ann_kept = ann[~ann['cnmf_component'].isin(excluded_old)].copy()
gobp_kept = gobp[~gobp['cnmf_component'].isin(excluded_old)].copy()

def split_term(s):
    s = str(s)
    m = re.search('\\(GO:(\\d+)\\)', s)
    go = 'GO:' + m.group(1) if m else ''
    txt = re.sub('\\s*\\(GO:\\d+\\)\\s*$', '', s).strip()
    return (txt, go)
_split = [split_term(s) for s in gobp_kept['top_BP_term']]
gobp_kept = gobp_kept.copy()
gobp_kept['gobp_term'] = [t for (t, g) in _split]
gobp_kept['gobp_id'] = [g for (t, g) in _split]

def top8_genes(s):
    parts = re.split('[;,]', str(s))
    parts = [p.strip() for p in parts if p.strip()]
    return ', '.join(parts[:8])
gobp_kept['top8_loading_genes'] = gobp_kept['top10_loading_genes'].map(top8_genes)
m = names_kept.merge(ann_kept[['cnmf_component', 'dominant_subclass', 'dominant_class']].rename(columns={'dominant_class': 'class'}), on='cnmf_component', how='left').merge(gobp_kept[['cnmf_component', 'gobp_term', 'gobp_id', 'top8_loading_genes']], on='cnmf_component', how='left')
m['new_pid_int'] = m['new_P'].map(lambda x: int(re.sub('[^0-9]', '', str(x))))
m = m.sort_values('new_pid_int').reset_index(drop=True)
supp = pd.DataFrame({'program_id': m['new_P'], 'functional_name': m['name_short'], 'dominant_class': m['class'], 'dominant_subclass': m['dominant_subclass'], 'top_GOBP_term': m['gobp_term'], 'GOBP_id': m['gobp_id'], 'NES': m['brain_term_NES'].round(3), 'FDR': m['fdr'], 'confidence': m['confidence'], 'top8_loading_genes': m['top8_loading_genes']})
CLASS_COL = {'exc': '#3B6FB6', 'inh': '#C25E5E', 'glia': '#5E9E7A', 'nonneuron': '#8A6FB0', 'vascular': '#C99A3B'}
CLASS_LABEL = {'exc': 'Excitatory', 'inh': 'Inhibitory', 'glia': 'Glia', 'nonneuron': 'Immune/Other', 'vascular': 'Vascular'}

def fdr_fmt(v):
    try:
        v = float(v)
    except Exception:
        return str(v)
    if v < 0.001:
        return f'{v:.1e}'
    return f'{v:.3f}'

def truncate(s, n):
    s = '' if pd.isna(s) else str(s)
    return s if len(s) <= n else s[:n - 1] + '…'
COLS = [('P', 0.0, 0.034, 'left'), ('Program', 0.034, 0.15, 'left'), ('Class', 0.15, 0.215, 'left'), ('Subclass', 0.215, 0.305, 'left'), ('Top GO:BP term', 0.305, 0.56, 'left'), ('NES', 0.56, 0.605, 'right'), ('FDR', 0.605, 0.67, 'right'), ('Top-8 loading genes', 0.67, 1.0, 'left')]
TRUNC = {'Program': 26, 'Subclass': 17, 'Top GO:BP term': 40, 'Top-8 loading genes': 52}
N = len(supp)
HALF = 27
blocks = [supp.iloc[:HALF].reset_index(drop=True), supp.iloc[HALF:].reset_index(drop=True)]
(FIG_W, FIG_H) = (17.0, 11.0)
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
HEADER_FS = 6.4
CELL_FS = 5.4
TITLE_FS = 11.0
SUB_FS = 7.0
fig.text(0.5, 0.975, 'Supplementary Data 3  |  Cross-region cNMF gene programs: GO:BP functional annotation (P1–P54)', ha='center', va='top', fontsize=TITLE_FS, fontweight='bold')
fig.text(0.5, 0.952, '', ha='center', va='top', fontsize=SUB_FS - 0.5, color='#222222')
lx = 0.3
for cls in ['exc', 'inh', 'glia', 'nonneuron', 'vascular']:
    fig.patches.append(Rectangle((lx, 0.928), 0.011, 0.011, transform=fig.transFigure, facecolor=CLASS_COL[cls], edgecolor='none', zorder=5))
    fig.text(lx + 0.014, 0.9335, CLASS_LABEL[cls], ha='left', va='center', fontsize=SUB_FS - 1.0, color='#222222')
    lx += 0.014 + 0.011 + 0.009 * len(CLASS_LABEL[cls]) + 0.01
PANEL_TOP = 0.915
PANEL_BOT = 0.018
PANEL_GAP = 0.03
PANEL_W = (1.0 - 2 * 0.012 - PANEL_GAP) / 2.0
PANEL_LEFTS = [0.012, 0.012 + PANEL_W + PANEL_GAP]
for (bi, (df, ax_left)) in enumerate(zip(blocks, PANEL_LEFTS)):
    ax = fig.add_axes([ax_left, PANEL_BOT, PANEL_W, PANEL_TOP - PANEL_BOT])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    nrows = len(df)
    header_h = 1.0 / (nrows + 1.6)
    row_h = (1.0 - header_h) / nrows

    def ytop_of_row(i):
        return 1.0 - header_h - i * row_h
    ax.add_patch(Rectangle((0, 1.0 - header_h), 1.0, header_h, facecolor='#2A2A2A', edgecolor='none', zorder=1))
    for (label, x0, x1, align) in COLS:
        if align == 'right':
            (tx, ha) = (x1 - 0.004, 'right')
        else:
            (tx, ha) = (x0 + 0.004, 'left')
        ax.text(tx, 1.0 - header_h / 2.0, label, ha=ha, va='center', fontsize=HEADER_FS, color='white', fontweight='bold', zorder=3)
    for i in range(nrows):
        r = df.iloc[i]
        y1 = ytop_of_row(i)
        y0 = y1 - row_h
        cls = r['dominant_class']
        ccol = CLASS_COL.get(cls, '#999999')
        sig = r['confidence'] == 'class_a'
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y0), 1.0, row_h, facecolor='#F2F2F2', edgecolor='none', zorder=0.5))
        ax.add_patch(Rectangle((0, y0), 0.006, row_h, facecolor=ccol, edgecolor='none', zorder=2))
        txt_col = '#111111' if sig else '#8A8A8A'
        weight = 'bold' if sig else 'normal'
        yc = (y0 + y1) / 2.0
        cells = {'P': r['program_id'], 'Program': truncate(r['functional_name'], TRUNC['Program']), 'Class': CLASS_LABEL.get(cls, str(cls)), 'Subclass': truncate(r['dominant_subclass'], TRUNC['Subclass']), 'Top GO:BP term': truncate(r['top_GOBP_term'], TRUNC['Top GO:BP term']), 'NES': f"{float(r['NES']):.2f}", 'FDR': fdr_fmt(r['FDR']), 'Top-8 loading genes': truncate(r['top8_loading_genes'], TRUNC['Top-8 loading genes'])}
        for (label, x0c, x1c, align) in COLS:
            val = cells[label]
            if label == 'Class':
                cc = ccol if sig else '#9AA0A6'
                ax.text(x0c + 0.01, yc, val, ha='left', va='center', fontsize=CELL_FS, color=cc, fontweight=weight, fontstyle='italic', zorder=3)
                continue
            if align == 'right':
                (tx, ha) = (x1c - 0.004, 'right')
            else:
                (tx, ha) = (x0c + 0.01 if label == 'P' else x0c + 0.004, 'left')
            fst = 'italic' if label == 'Top-8 loading genes' else 'normal'
            ax.text(tx, yc, val, ha=ha, va='center', fontsize=CELL_FS, color=txt_col, fontweight=weight, fontstyle=fst, zorder=3)
    ax.add_patch(Rectangle((0, 0), 1.0, 1.0, fill=False, edgecolor='#BBBBBB', linewidth=0.4, zorder=4))
pdf_path = os.path.join(FIGDIR, 'ed_fig4_program_table_54.pdf')
png_path = os.path.join(FIGDIR, 'ed_fig4_program_table_54.png')
fig.savefig(pdf_path)
fig.savefig(png_path, dpi=300)
