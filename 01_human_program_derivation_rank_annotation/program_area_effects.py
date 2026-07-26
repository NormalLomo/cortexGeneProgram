#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os
import zipfile
import numpy as np
import pandas as pd
import matplotlib as mpl
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as _mpl_font
_mpl_font.rcParams['font.family'] = 'sans-serif'
_mpl_font.rcParams['font.sans-serif'] = ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans']
_mpl_font.rcParams['pdf.fonttype'] = 42
_mpl_font.rcParams['ps.fonttype'] = 42
_mpl_font.rcParams['svg.fonttype'] = 'none'
_mpl_font.rcParams['mathtext.fontset'] = 'dejavusans'
_mpl_font.rcParams['mathtext.default'] = 'regular'
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
mpl.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 6.0, 'axes.titlesize': 7.0, 'axes.labelsize': 6.0, 'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5, 'legend.fontsize': 5.5, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.0, 'ytick.major.size': 2.0, 'axes.edgecolor': '#333333', 'savefig.dpi': 600})
PROJ = str(PROJECT_ROOT)
SCORES = os.path.join(PROJ, 'results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv')
OBS = os.path.join(PROJ, 'inputs/snRNA_1M_obs.csv')
NAMES = os.path.join(PROJ, 'results/crossregion_v1/program_names.tsv')
VARI = os.path.join(PROJ, 'results/crossregion_v1/program_variability.tsv')
RETAIN_MAP = os.environ.get('RETAIN_MAP', os.path.join(PROJ, 'tables/TableS3_program_annotation.tsv'))
SUMMARY_ZIP = os.environ.get('PROGRAM_SUMMARY_ZIP')
OUTDIR = os.environ.get('FIGS3_OUTPUT_DIR', os.path.join(PROJ, 'figures/extended'))
RESDIR = os.path.join(PROJ, 'results/crossregion_v1')
os.makedirs(OUTDIR, exist_ok=True)
CCOL = {'us': '#D7642C', 'edlein': '#3C6E9C'}
SIGCOL = {'class_a': '#2E6E4E', 'class_b': '#B0883B'}
_map_df = _annotation_table(pd.read_csv(RETAIN_MAP, sep='\t'))
_map_df['new_int'] = _map_df['new_P'].astype(str).str.removeprefix('P').astype(int)
_map_df['old_int'] = _map_df['cnmf_component'].astype(int)
_OLD_TO_LABEL = dict(zip(_map_df['old_int'], 'P' + _map_df['new_int'].astype(str)))
EXCLUDED_OLD = set(range(1, 61)) - set(_OLD_TO_LABEL)
_names_df = _map_df[['new_int', 'confidence']].rename(columns={'new_int': 'new_P'})
_CONF_MAP = {}
for (_, row) in _names_df.iterrows():
    _CONF_MAP[int(row['new_P'])] = str(row['confidence'])

def _star_for_old(old_p):
    label = _OLD_TO_LABEL.get(old_p, '')
    if label.startswith('P'):
        try:
            new_int = int(label[1:])
            return '*' if _CONF_MAP.get(new_int, '') == 'class_b' else ''
        except ValueError:
            return ''
    return ''

def prog_label(old_p):
    base = _OLD_TO_LABEL.get(old_p, 'P%d' % old_p)
    return base + _star_for_old(old_p)
REPLOT_ONLY = os.environ.get('REPLOT_ONLY', '0') == '1'

def _build_tab_from_summary():
    if not SUMMARY_ZIP:
        raise RuntimeError()
    with zipfile.ZipFile(SUMMARY_ZIP) as zf:
        with zf.open('programs_master.tsv') as fh:
            master = _annotation_table(pd.read_csv(fh, sep='\t'))
    required = {'program', 'name_short', 'confidence', 'cohort_partial_eta2', 'region_partial_eta2', 'subclass_eta2', 'cohort_gap_us_minus_edlein_SD', 'area_variable'}
    absent = required - set(master.columns)
    master['program'] = master['program'].astype(int)
    tab = master[master['program'].isin(_OLD_TO_LABEL)].copy()
    tab['depth_corr_log10nCount'] = np.nan
    tab['area_variable'] = tab['area_variable'].astype(str).str.lower().eq('true')
    tab = tab.set_index('program').sort_values('cohort_partial_eta2', ascending=False)
    return tab

def _build_tab_from_data():
    S = _annotation_table(pd.read_csv(SCORES, sep='\t', index_col=0))
    S.columns = [int(c) for c in S.columns]
    S = S.sort_index(axis=1)
    S = S.loc[:, sorted(_OLD_TO_LABEL)]
    PROGS = list(S.columns)
    obs = _annotation_table(pd.read_csv(OBS, index_col=0, usecols=['Unnamed: 0', 'batch', 'region', 'subclass', 'class', 'nCount_RNA', 'nFeature_RNA']))
    obs = obs.reindex(S.index)
    for c in ['batch', 'region', 'subclass', 'class']:
        obs[c] = obs[c].astype(str)
    log_nc = np.log10(obs['nCount_RNA'].values.astype(float) + 1.0)

    def design_matrix(factors):
        n = len(factors[0])
        cols = [np.ones(n)]
        for f in factors:
            cats = pd.Categorical(f)
            codes = cats.codes
            K = len(cats.categories)
            for k in range(1, K):
                cols.append((codes == k).astype(float))
        return np.column_stack(cols)

    def rss_all(X, Y):
        XtX = X.T @ X
        XtY = X.T @ Y
        beta = np.linalg.solve(XtX, XtY)
        Yhat = X @ beta
        resid = Y - Yhat
        return (resid * resid).sum(axis=0)
    Y = S.values.astype(np.float64)
    Ymean = Y.mean(axis=0, keepdims=True)
    Ystd = Y.std(axis=0, ddof=0, keepdims=True)
    Ystd[Ystd == 0] = 1.0
    Yz = (Y - Ymean) / Ystd
    N = Yz.shape[0]
    TSS = (Yz * Yz).sum(axis=0)
    sub = obs['subclass']
    reg = obs['region']
    bat = obs['batch']
    X_sub = design_matrix([sub])
    X_sc = design_matrix([sub, bat])
    X_sr = design_matrix([sub, reg])
    X_bat = design_matrix([bat])
    X_reg = design_matrix([reg])
    rss_sub = rss_all(X_sub, Yz)
    rss_sc = rss_all(X_sc, Yz)
    rss_sr = rss_all(X_sr, Yz)
    rss_bat = rss_all(X_bat, Yz)
    rss_reg = rss_all(X_reg, Yz)
    ss_cohort = rss_sub - rss_sc
    ss_region = rss_sub - rss_sr
    cohort_partial_eta2 = ss_cohort / (ss_cohort + rss_sc)
    region_partial_eta2 = ss_region / (ss_region + rss_sr)
    subclass_eta2 = 1.0 - rss_sub / TSS
    cohort_marg_eta2 = 1.0 - rss_bat / TSS
    region_marg_eta2 = 1.0 - rss_reg / TSS
    lnc = (log_nc - log_nc.mean()) / log_nc.std()
    depth_corr = (Yz * lnc[:, None]).sum(axis=0) / N
    is_us = bat.values == 'us'
    gap = Yz[is_us].mean(axis=0) - Yz[~is_us].mean(axis=0)
    names_raw = _annotation_table(pd.read_csv(NAMES, sep='\t'))
    names = names_raw[names_raw['new_P'].astype(str).str.startswith('P')].set_index('new_P')
    vari_raw = _annotation_table(pd.read_csv(VARI, sep='\t'))
    vari = vari_raw[vari_raw['new_P'].astype(str).str.startswith('P')].set_index('new_P') if 'new_P' in vari_raw.columns else vari_raw.set_index('program')
    AREA_VARIABLE = set(vari.index[vari['class'] == 'variable'].tolist())
    tab = pd.DataFrame({'program': PROGS, 'name_short': [names.loc[p, 'name_short'] if p in names.index else str(p) for p in PROGS], 'confidence': [names.loc[p, 'confidence'] if p in names.index else 'NA' for p in PROGS], 'cohort_partial_eta2': cohort_partial_eta2, 'region_partial_eta2': region_partial_eta2, 'subclass_eta2': subclass_eta2, 'cohort_marginal_eta2': cohort_marg_eta2, 'region_marginal_eta2': region_marg_eta2, 'depth_corr_log10nCount': depth_corr, 'cohort_gap_us_minus_edlein_SD': gap, 'area_variable': [p in AREA_VARIABLE for p in PROGS]}).set_index('program')
    tab = tab.sort_values('cohort_partial_eta2', ascending=False)
    OUTTAB = os.path.join(RESDIR, 'program_cohort_eta2.tsv')
    tab.to_csv(OUTTAB, sep='\t')
    return tab
if SUMMARY_ZIP:
    tab = _build_tab_from_summary()
elif REPLOT_ONLY:
    OUTTAB = os.path.join(RESDIR, 'program_cohort_eta2.tsv')
    tab = _annotation_table(pd.read_csv(OUTTAB, sep='\t')).set_index('program')
    tab.index = tab.index.astype(int)
    tab = tab.loc[sorted(_OLD_TO_LABEL)]
else:
    tab = _build_tab_from_data()
ce = tab['cohort_partial_eta2']
for p in ce.head(8).index:
    r = tab.loc[p]
sub14 = tab[tab['area_variable']].sort_values('region_partial_eta2', ascending=False)
for p in sub14.index:
    r = sub14.loc[p]
    ratio = r['region_partial_eta2'] / max(r['cohort_partial_eta2'], 1e-09)
dflag = tab[tab['depth_corr_log10nCount'].abs() > 0.3].sort_values('depth_corr_log10nCount')
for p in dflag.index:
    r = dflag.loc[p]
fig = plt.figure(figsize=(7.2, 7.6))
gs = GridSpec(3, 1, figure=fig, height_ratios=[1.35, 1.0, 1.05], hspace=0.62, left=0.1, right=0.965, top=0.945, bottom=0.175)

def sigc(p):
    return SIGCOL.get(tab.loc[p, 'confidence'], '#888888')
axa = fig.add_subplot(gs[0, 0])
xr = tab['region_partial_eta2'].values
yc = tab['cohort_partial_eta2'].values
cols = [sigc(p) for p in tab.index]
is_area_variable = tab['area_variable'].values
axa.scatter(xr[~is_area_variable], yc[~is_area_variable], s=18, c=[cols[i] for i in range(len(cols)) if not is_area_variable[i]], alpha=0.78, linewidths=0.3, edgecolors='white', zorder=3)
axa.scatter(xr[is_area_variable], yc[is_area_variable], s=34, facecolors=[cols[i] for i in range(len(cols)) if is_area_variable[i]], edgecolors='#111111', linewidths=0.9, alpha=0.95, zorder=4)
mx = max(xr.max(), yc.max()) * 1.08
axa.plot([0, mx], [0, mx], ls='--', lw=0.6, c='#999999', zorder=1)
for (thr, lab) in [(0.05, '0.05'), (0.1, '0.10')]:
    axa.axhline(thr, ls=':', lw=0.5, c='#C24C4C', zorder=1)
    axa.text(mx * 0.995, thr, ' cohort=%s' % lab, fontsize=4.6, color='#C24C4C', va='bottom', ha='right')
to_label = set(ce.head(6).index) | set(tab.index[is_area_variable])
try:
    from adjustText import adjust_text
    texts = []
    for p in to_label:
        texts.append(axa.text(tab.loc[p, 'region_partial_eta2'], tab.loc[p, 'cohort_partial_eta2'], prog_label(p), fontsize=4.7, color='#222222', zorder=5))
    adjust_text(texts, ax=axa, only_move={'text': 'xy'}, expand_text=(1.05, 1.2), expand_points=(1.05, 1.2), force_text=(0.4, 0.55), force_points=(0.25, 0.3), arrowprops=dict(arrowstyle='-', color='#888888', lw=0.25, alpha=0.7), max_move=8.0)
except Exception as _e:
    for p in to_label:
        axa.annotate(prog_label(p), (tab.loc[p, 'region_partial_eta2'], tab.loc[p, 'cohort_partial_eta2']), fontsize=4.7, xytext=(2.2, 2.2), textcoords='offset points', color='#222222', zorder=5)
axa.set_xlabel('Region partial-$\\eta^2$  (usage ~ subclass + region)')
axa.set_ylabel('Cohort partial-$\\eta^2$\n(usage ~ subclass + cohort)')
axa.set_xlim(-mx * 0.02, mx)
axa.set_ylim(-mx * 0.02, mx)
axa.set_title('(a)  Per-program cohort vs region effect, controlling for subclass identity', loc='left', pad=4)
leg = [Patch(fc=SIGCOL['class_a'], label='class_a'), Patch(fc=SIGCOL['class_b'], label='class_b'), plt.Line2D([0], [0], marker='o', ls='', mfc='#ccc', mec='#111', mew=0.9, ms=5, label='area region-variable (n=%d)' % int(is_area_variable.sum()))]
axa.legend(handles=leg, loc='upper left', frameon=False, fontsize=5.0, handletextpad=0.4, borderaxespad=0.2)
axb = fig.add_subplot(gs[1, 0])
order = tab.index.tolist()
vals = tab['cohort_partial_eta2'].values
barcols = [sigc(p) for p in order]
xpos = np.arange(len(order))
bars = axb.bar(xpos, vals, color=barcols, width=0.78, linewidth=0.2, edgecolor='white')
for (i, p) in enumerate(order):
    if tab.loc[p, 'area_variable']:
        bars[i].set_edgecolor('#111111')
        bars[i].set_linewidth(0.8)
axb.axhline(0.05, ls=':', lw=0.6, c='#C24C4C')
axb.axhline(0.1, ls='--', lw=0.6, c='#8B2E2E')
axb.text(len(order) - 0.5, 0.05, ' 0.05', fontsize=4.8, color='#C24C4C', va='bottom', ha='right')
axb.text(len(order) - 0.5, 0.1, ' 0.10', fontsize=4.8, color='#8B2E2E', va='bottom', ha='right')
axb.set_xticks(xpos)
axb.set_xticklabels([prog_label(p) for p in order], rotation=90, fontsize=4.2)
axb.set_xlim(-0.7, len(order) - 0.3)
axb.set_ylabel('Cohort partial-$\\eta^2$')
axb.set_title('(b)  All 54 programs ranked by cohort partial-$\\eta^2$  (outlined = area variable; color = annotation class)', loc='left', pad=4)
axb.margins(y=0.12)
axc = fig.add_subplot(gs[2, 0])
s14 = tab[tab['area_variable']].sort_values('region_partial_eta2', ascending=False)
p14 = s14.index.tolist()
xp = np.arange(len(p14))
w = 0.4
axc.bar(xp - w / 2, s14['region_partial_eta2'].values, width=w, color='#3C6E9C', label='region partial-$\\eta^2$', linewidth=0.2, edgecolor='white')
axc.bar(xp + w / 2, s14['cohort_partial_eta2'].values, width=w, color='#D7642C', label='cohort partial-$\\eta^2$', linewidth=0.2, edgecolor='white')
axc.axhline(0.05, ls=':', lw=0.5, c='#555555')
axc.set_xticks(xp)
lbls = ['%s %s' % (prog_label(p), s14.loc[p, 'name_short']) for p in p14]
axc.set_xticklabels(lbls, rotation=40, ha='right', fontsize=4.6)
axc.set_ylabel('partial-$\\eta^2$')
n_regdom = int((s14['region_partial_eta2'] > s14['cohort_partial_eta2']).sum())
n_cohdom = len(s14) - n_regdom
axc.set_title('(c)  Programs with area variation: %d/%d region-dominated, %d/%d cohort-dominated' % (n_regdom, len(s14), n_cohdom, len(s14)), loc='left', pad=4)
ymax_c = max(s14['region_partial_eta2'].max(), s14['cohort_partial_eta2'].max())
for (i, p) in enumerate(p14):
    if s14.loc[p, 'cohort_partial_eta2'] > s14.loc[p, 'region_partial_eta2']:
        axc.plot(i, ymax_c * 1.04, marker='v', ms=3.2, color='#8B2E2E', zorder=5)
axc.legend(loc='upper right', frameon=False, fontsize=5.2, handletextpad=0.4)
axc.margins(y=0.12)
OUTPDF = os.path.join(OUTDIR, 'ed_program_area_effects.pdf')
OUTPNG = os.path.join(OUTDIR, 'ed_program_area_effects.png')
OUTSVG = os.path.join(OUTDIR, 'ed_program_area_effects.svg')
fig.savefig(OUTPDF)
fig.savefig(OUTPNG, dpi=300)
fig.savefig(OUTSVG)
