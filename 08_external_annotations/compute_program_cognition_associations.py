#!/usr/bin/env python
import warnings
warnings.filterwarnings('ignore')
import os, json, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BASE = str(PROJECT_ROOT)
OUT = f'{BASE}/results/crossregion_v1/program_cognition'
FIGD = f'{BASE}/figures/fig7'
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIGD, exist_ok=True)
os.makedirs(f'{OUT}/ns_data', exist_ok=True)
RNG = np.random.default_rng(0)
REGION_MAP = {'FPPFC': ([1], 'high'), 'DLPFC': ([4], 'high'), 'VLPFC': ([5, 6], 'high'), 'M1': ([7], 'high'), 'ACC': ([28, 29], 'high'), 'S1': ([17], 'high'), 'PoCG': ([17], 'high-collinear'), 'S1E': ([17], 'low'), 'STG': ([9, 10], 'high'), 'ITG': ([14, 15, 16], 'high'), 'SMG': ([19, 20], 'high'), 'SPL': ([18], 'high'), 'AG': ([21], 'high'), 'V1': ([24, 47, 48], 'high')}
TERMS = ['default mode', 'working memory', 'language', 'attention', 'semantic', 'episodic memory', 'visual', 'motor', 'pain', 'salience', 'executive', 'reward']
RENUM_PATH = f'{BASE}/results/crossregion_v1/program_renumber_map.tsv'
renum_df = pd.read_csv(RENUM_PATH, sep='\t')
renum_df = renum_df[renum_df['status'].astype(str).str.lower().eq('kept') & renum_df['new_P'].notna()].copy()
renum_df['old_P'] = renum_df['old_P'].astype(int)
renum_df['new_P'] = renum_df['new_P'].astype(int)
renum_df = renum_df.sort_values('new_P')
if len(renum_df) != 54 or renum_df['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
if renum_df['new_P'].tolist() != list(range(1, 55)):
    raise ValueError()
old_to_new = dict(zip(renum_df['old_P'].astype(str), 'P' + renum_df['new_P'].astype(str)))
from nilearn import datasets, image
from nilearn.maskers import NiftiLabelsMasker
ho = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
ho_img = ho.maps if hasattr(ho, 'maps') else image.load_img(ho['maps'])
ho_data = np.asarray(ho_img.dataobj)
ho_aff = ho_img.affine
from nimare.extract import fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.meta.cbma.mkda import MKDAChi2
from nimare.dataset import Dataset
ns_files = fetch_neurosynth(data_dir=f'{OUT}/ns_data', version='7', overwrite=False, source='abstract', vocab='terms')
ns = ns_files[0] if isinstance(ns_files, list) else ns_files
dset_path = f'{OUT}/ns_data/neurosynth_dataset.pkl.gz'
if os.path.exists(dset_path):
    dset = Dataset.load(dset_path)
else:
    dset = convert_neurosynth_to_dataset(coordinates_file=ns['coordinates'], metadata_file=ns['metadata'], annotations_files=ns['features'])
    dset.save(dset_path)

def term_map(term):
    feat = 'terms_abstract_tfidf__' + term
    if feat not in dset.annotations.columns:
        cands = [c for c in dset.annotations.columns if c.endswith('__' + term)]
        if not cands:
            cands = [c for c in dset.annotations.columns if term in c]
        if not cands:
            raise KeyError()
        feat = cands[0]
    ids = dset.get_studies_by_label(feat, label_threshold=0.001)
    ids_other = list(set(dset.ids) - set(ids))
    d1 = dset.slice(ids)
    d2 = dset.slice(ids_other)
    mkda = MKDAChi2()
    res = mkda.fit(d1, d2)
    img = res.get_map('z_desc-association')
    return (img, len(ids), feat)
region_term = pd.DataFrame(index=list(REGION_MAP.keys()), columns=TERMS, dtype=float)
term_info = {}
for t in TERMS:
    try:
        (zimg, nstud, feat) = term_map(t)
    except Exception as e:
        raise
    ho_rs = image.resample_to_img(ho_img, zimg, interpolation='nearest')
    ho_rs_data = np.asarray(ho_rs.dataobj)
    zdata = np.asarray(zimg.dataobj)
    for (reg, (idxs, conf)) in REGION_MAP.items():
        mask = np.isin(ho_rs_data, idxs)
        vals = zdata[mask]
        vals = vals[np.isfinite(vals)]
        region_term.loc[reg, t] = float(np.nanmean(vals)) if vals.size else np.nan
    term_info[t] = {'n_studies': int(nstud), 'feature': feat}
region_term.to_csv(f'{OUT}/region_term_neurosynth.tsv', sep='\t')
json.dump(term_info, open(f'{OUT}/term_info.json', 'w'), indent=2)
try:
    sch = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
    sch_img = image.load_img(sch['maps'])
    labs = [l.decode() if isinstance(l, bytes) else l for l in sch.labels]
    dmn_idx = [i + 1 for (i, l) in enumerate(labs) if 'Default' in l]
    sch_data = np.asarray(sch_img.dataobj)
    dmn_mask_img = image.new_img_like(sch_img, np.isin(sch_data, dmn_idx).astype('int8'))
    json.dump({'n_dmn_parcels': len(dmn_idx)}, open(f'{OUT}/dmn_meta.json', 'w'))
except Exception as e:
    raise
pz = pd.read_csv(f'{BASE}/results/crossregion_v1/program_region_zscore.tsv', sep='\t', index_col=0)
pz.columns = [str(c) for c in pz.columns]
names_df = pd.read_csv(f'{BASE}/results/crossregion_v1/program_names.tsv', sep='\t')
names_df = names_df.merge(renum_df[['old_P', 'new_P']], left_on='cnmf_component', right_on='old_P')
names_df['program'] = 'P' + names_df['new_P'].astype(str)
names_df = names_df.set_index('program')
name_short = names_df['name_short'].to_dict()
conf_map = names_df['confidence'].to_dict()
regions = [r for r in REGION_MAP if r in pz.index and r in region_term.index]
RT = region_term.loc[regions].astype(float)
source_columns = [str(old) for old in renum_df['old_P']]
if not set(source_columns).issubset(pz.columns):
    raise ValueError()
PZ = pz.loc[regions, source_columns].astype(float)
PZ.columns = [old_to_new[column] for column in PZ.columns]
if list(PZ.columns) != [f'P{i}' for i in range(1, 55)]:
    raise ValueError()
progs = list(PZ.columns)
good_terms = [t for t in TERMS if RT[t].notna().sum() >= 6]
RT = RT[good_terms]

def spear_matrix(PZmat, RTmat):
    (nP, nT) = (PZmat.shape[1], RTmat.shape[1])
    M = np.full((nP, nT), np.nan)
    for (j, p) in enumerate(PZmat.columns):
        x = PZmat[p].values
        for (k, t) in enumerate(RTmat.columns):
            y = RTmat[t].values
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 6:
                M[j, k] = spearmanr(x[ok], y[ok]).correlation
    return M
R = spear_matrix(PZ, RT)
Rdf = pd.DataFrame(R, index=progs, columns=good_terms)
Rdf.to_csv(f'{OUT}/program_term_spearman.tsv', sep='\t')
NPERM = 10000
perm_max = np.zeros((NPERM,))
exceed = np.zeros_like(R)
idx = np.arange(len(regions))
PZv = PZ.values
RTv = RT.values
for b in range(NPERM):
    perm = RNG.permutation(idx)
    PZp = PZv[perm, :]
    for k in range(RTv.shape[1]):
        y = RTv[:, k]
        oky = np.isfinite(y)
        for j in range(PZp.shape[1]):
            x = PZp[:, j]
            ok = oky & np.isfinite(x)
            if ok.sum() >= 6:
                rp = spearmanr(x[ok], y[ok]).correlation
                if abs(rp) >= abs(R[j, k]):
                    exceed[j, k] += 1
pmat = (exceed + 1) / (NPERM + 1)
Pdf = pd.DataFrame(pmat, index=progs, columns=good_terms)
Pdf.to_csv(f'{OUT}/program_term_permp.tsv', sep='\t')
flat = pmat.flatten()
(rej, q, _, _) = multipletests(flat, method='fdr_bh')
Qdf = pd.DataFrame(q.reshape(pmat.shape), index=progs, columns=good_terms)
Qdf.to_csv(f'{OUT}/program_term_fdr.tsv', sep='\t')
recs = []
for p in progs:
    ns = name_short.get(p, '')
    cf = conf_map.get(p, '')
    for t in good_terms:
        recs.append((p, ns, cf, t, Rdf.loc[p, t], Pdf.loc[p, t], Qdf.loc[p, t]))
top = pd.DataFrame(recs, columns=['program', 'name_short', 'confidence', 'term', 'spearman_r', 'perm_p', 'fdr_q'])
top = top.reindex(top.spearman_r.abs().sort_values(ascending=False).index)
top.to_csv(f'{OUT}/top_program_term_pairs.tsv', sep='\t', index=False)
dmn_term = 'default mode' if 'default mode' in good_terms else good_terms[0]
dmn_rank = Rdf[dmn_term].sort_values(ascending=False)
dmn_rank.to_csv(f'{OUT}/program_DMN_ranking.tsv', sep='\t', header=[dmn_term])
import matplotlib
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
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
plt.rcParams.update({'font.size': 6, 'axes.linewidth': 0.5, 'font.family': 'sans-serif', 'font.sans-serif': ['Nimbus Sans', 'Liberation Sans', 'DejaVu Sans'], 'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.major.width': 0.4, 'ytick.major.width': 0.4, 'pdf.fonttype': 42, 'ps.fonttype': 42})
Rc = Rdf.fillna(0)
row_order = leaves_list(linkage(pdist(Rc.values), method='average')) if Rc.shape[0] > 2 else range(Rc.shape[0])
col_order = leaves_list(linkage(pdist(Rc.values.T), method='average')) if Rc.shape[1] > 2 else range(Rc.shape[1])
Rc = Rc.iloc[row_order, col_order]

def make_prog_lab(program, maxlen=34):
    return f"{program} {name_short.get(program, '')}"[:maxlen]
prog_lab = [make_prog_lab(p) for p in Rc.index]
term_lab = list(Rc.columns)
fig = plt.figure(figsize=(180 / 25.4, 220 / 25.4), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.1, wspace=0.04, hspace=0.07)
outer = fig.add_gridspec(3, 1, height_ratios=[3.2, 1.25, 1.55])
gtop = outer[0].subgridspec(1, 2, width_ratios=[1.55, 0.95], wspace=0.0)
gmid = outer[1].subgridspec(1, 3, wspace=0.04)
gbot = outer[2].subgridspec(1, 1)

def panel_tag(ax, s, dx=-30):
    ax.annotate(s, xy=(0, 1), xycoords='axes fraction', xytext=(dx, 8), textcoords='offset points', fontsize=10, fontweight='bold', va='bottom', ha='left')
axa = fig.add_subplot(gtop[0, 0])
vmax = np.nanmax(np.abs(Rc.values))
vmax = max(vmax, 0.3)
im = axa.imshow(Rc.values, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
axa.set_xticks(range(len(term_lab)))
axa.set_xticklabels(term_lab, rotation=45, ha='right', fontsize=5.5)
axa.set_yticks(range(len(prog_lab)))
axa.set_yticklabels(prog_lab, fontsize=4.4)
axa.tick_params(axis='y', pad=1, length=1.5)
Qc = Qdf.reindex(index=Rc.index, columns=Rc.columns)
for i in range(Rc.shape[0]):
    for j in range(Rc.shape[1]):
        if Qc.values[i, j] < 0.1:
            axa.text(j, i, '\\N{DAGGER}', ha='center', va='center', fontsize=6, color='k')
axa.set_title('Program x cognitive-term (Spearman; dagger FDR<0.1)', fontsize=6, pad=3)
cb = fig.colorbar(im, ax=axa, fraction=0.035, pad=0.015)
cb.ax.tick_params(labelsize=5)
cb.set_label('Spearman r', fontsize=5)
panel_tag(axa, 'a', dx=-100)
axb = fig.add_subplot(gtop[0, 1])
dr = Rdf[dmn_term].dropna().sort_values()
sel = pd.concat([dr.head(8), dr.tail(8)])
ylab = [make_prog_lab(p, maxlen=30) for p in sel.index]
cols = ['#b2182b' if v > 0 else '#2166ac' for v in sel.values]
axb.barh(range(len(sel)), sel.values, color=cols, height=0.78)
axb.set_yticks(range(len(sel)))
axb.set_yticklabels(ylab, fontsize=4.8)
axb.set_ylim(-0.6, len(sel) - 0.4)
axb.axvline(0, color='k', lw=0.4)
axb.set_xlabel(f"Spearman r vs '{dmn_term}'", fontsize=5.5)
axb.set_title('Programs ranked by\nDMN correlation', fontsize=6, pad=3)
panel_tag(axb, 'b')
ex = top.dropna(subset=['spearman_r']).copy()
chosen = []
seen = set()
for (_, r) in ex.iterrows():
    if r['term'] in seen:
        continue
    chosen.append(r)
    seen.add(r['term'])
    if len(chosen) == 3:
        break
for (ci, r) in enumerate(chosen):
    axx = fig.add_subplot(gmid[0, ci])
    program = r['program']
    t = r['term']
    x = PZ[program].values
    y = RT[t].values
    ok = np.isfinite(x) & np.isfinite(y)
    axx.scatter(x[ok], y[ok], s=12, c='#444', zorder=3)
    for (rg, xi, yi) in zip(np.array(regions)[ok], x[ok], y[ok]):
        axx.annotate(rg, (xi, yi), fontsize=4.3, xytext=(2, 2), textcoords='offset points')
    if ok.sum() > 2:
        b = np.polyfit(x[ok], y[ok], 1)
        xs = np.linspace(x[ok].min(), x[ok].max(), 20)
        axx.plot(xs, np.polyval(b, xs), color='#b2182b', lw=0.8)
    ns_disp = name_short.get(program, '')
    axx.set_xlabel(f'{program} {ns_disp}'[:24] + ' (region z)', fontsize=5)
    axx.set_ylabel(f"NS '{t}' z", fontsize=5)
    axx.set_title(f"r={r['spearman_r']:.2f}, q={r['fdr_q']:.2f}", fontsize=5.6)
    axx.tick_params(labelsize=4.5)
    if ci == 0:
        panel_tag(axx, 'c')
axd = fig.add_subplot(gbot[0, 0])
RTplot = region_term.loc[regions, good_terms].astype(float)
rr = leaves_list(linkage(pdist(np.nan_to_num(RTplot.values)), 'average')) if RTplot.shape[0] > 2 else range(RTplot.shape[0])
RTp = RTplot.iloc[rr]
imd = axd.imshow(RTp.values, aspect='auto', cmap='viridis')
axd.set_xticks(range(len(good_terms)))
axd.set_xticklabels(good_terms, rotation=45, ha='right', fontsize=5.5)
ylabs = [f"{r}{('*' if REGION_MAP[r][1] == 'low' else '')}" for r in RTp.index]
axd.set_yticks(range(len(RTp.index)))
axd.set_yticklabels(ylabs, fontsize=5.5)
axd.set_title('Neurosynth region x term cognitive atlas (mean z in ROI; *S1E low-conf)', fontsize=6, pad=3)
cb2 = fig.colorbar(imd, ax=axd, fraction=0.02, pad=0.01)
cb2.ax.tick_params(labelsize=5)
cb2.set_label('mean z', fontsize=5)
panel_tag(axd, 'd')
fig.suptitle('Fig. S10 | Program-cognition imaging-transcriptomics (Neurosynth DMN), exploratory analysis (n=14 regions)', fontsize=7)
fig.savefig(f'{FIGD}/figS10_cognition.pdf', dpi=400)
fig.savefig(f'{FIGD}/figS10_cognition.png', dpi=300)
