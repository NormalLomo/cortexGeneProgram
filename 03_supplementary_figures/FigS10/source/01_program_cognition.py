#!/usr/bin/env python
"""
Fig 7 = gene-program x cognitive-term/functional-network association via Neurosynth.

Method mirrors the reference NC paper (Nat Commun 2025, s41467-025-62793-9), which
"employed the Neurosynth framework for meta-analysis, which provides term-region
associations" to relate cortical biology to cognitive functions (e.g. working memory
~ prefrontal; pain/motor ~ central sulcus). Same 14 cortical regions.

Pipeline:
  1. Region -> Harvard-Oxford cortical ROI mask (MNI152 2mm). Mapping locked in REGION_MAP.
  2. Neurosynth v7 (nimare) -> per-term meta-analytic z map (association test) for ~12 terms.
  3. Sample each term map's mean z within each region ROI -> region x term matrix.
  4. program_region_zscore.tsv (14x60) x region x term -> program x term Spearman matrix.
  5. FDR (BH) across 60 x n_term; null = region-label permutation (10k) -> empirical p (weak surrogate).
  6. DMN cross-check via Schaefer-400/Yeo-7 Default-network mask.
  7. Panels a-d, native vector PDF+PNG.

CAVEATS (printed + in report): n=14 regions (low power); region-label permutation null is a
weak surrogate (does NOT respect spatial autocorrelation -> liberal); ROI alignment is gyral-
level (S1E has no standard atlas correlate, mapped to Postcentral ~ S1/PoCG, flagged low-conf);
program z-scores are cohort-uncorrected.
"""
import warnings; warnings.filterwarnings("ignore")
import os, json, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests

BASE = "CORTEX_PROGRAM_ROOT"
OUT  = f"{BASE}/results/crossregion_v1/program_cognition"
FIGD = f"{BASE}/figures/fig7"
os.makedirs(OUT, exist_ok=True); os.makedirs(FIGD, exist_ok=True)
os.makedirs(f"{OUT}/ns_data", exist_ok=True)
RNG = np.random.default_rng(0)

# ------------------------------------------------------------------ region map
# Harvard-Oxford cort-maxprob-thr25-2mm label indices (see labels list).
REGION_MAP = {
    "FPPFC": ([1], "high"),                 # Frontal Pole
    "DLPFC": ([4], "high"),                 # Middle Frontal Gyrus
    "VLPFC": ([5,6], "high"),               # IFG tri + oper
    "M1":    ([7], "high"),                 # Precentral Gyrus
    "ACC":   ([28,29], "high"),             # Paracingulate + Cingulate ant
    "S1":    ([17], "high"),                # Postcentral Gyrus
    "PoCG":  ([17], "high-collinear"),      # Postcentral (= S1 territory)
    "S1E":   ([17], "low"),                 # no std correlate -> Postcentral (~S1)
    "STG":   ([9,10], "high"),              # STG ant+post
    "ITG":   ([14,15,16], "high"),          # ITG ant+post+TO
    "SMG":   ([19,20], "high"),             # SMG ant+post
    "SPL":   ([18], "high"),                # Superior Parietal Lobule
    "AG":    ([21], "high"),                # Angular Gyrus
    "V1":    ([24,47,48], "high"),          # Intracalcarine + Supracalcarine + Occipital Pole
}

# Neurosynth terms (~12) spanning canonical cognitive domains + DMN
TERMS = ["default mode","working memory","language","attention","semantic",
         "episodic memory","visual","motor","pain","salience","executive","reward"]

# ------------------------------------------------------------------ PROGRAM MAP
RENUM_PATH = f"{BASE}/results/crossregion_v1/program_renumber_map.tsv"
renum_df = pd.read_csv(RENUM_PATH, sep="\t")
# old_P (int) -> new_P (int or "EXCLUDED")
EXCLUDED_OLD = set(renum_df.loc[renum_df.new_P == "EXCLUDED", "old_P"].astype(int))
old_to_new = {}
for _, row in renum_df.iterrows():
    op = int(row.old_P)
    if row.new_P == "EXCLUDED":
        old_to_new[op] = None
    else:
        old_to_new[op] = int(row.new_P)
print(f"[renumber] excluded old_P: {sorted(EXCLUDED_OLD)}", flush=True)
print(f"[renumber] kept: {sum(v is not None for v in old_to_new.values())} programs", flush=True)

def disp_label(old_p_str):
    """Convert old_P string column name to display string 'P{new_n}'."""
    op = int(old_p_str)
    nn = old_to_new.get(op)
    if nn is None:
        return f"[FLAG-EXCLUDED-oldP{op}]"
    return str(nn)

# ------------------------------------------------------------------ atlases
print("[1] fetching atlases ...", flush=True)
from nilearn import datasets, image
from nilearn.maskers import NiftiLabelsMasker
ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
ho_img = ho.maps if hasattr(ho,"maps") else image.load_img(ho["maps"])
ho_data = np.asarray(ho_img.dataobj)
ho_aff = ho_img.affine
print("    HO atlas shape", ho_data.shape, flush=True)

# ------------------------------------------------------------------ neurosynth maps
print("[2] building Neurosynth term maps (this is the slow step) ...", flush=True)
from nimare.extract import fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.meta.cbma.mkda import MKDAChi2
from nimare.dataset import Dataset

ns_files = fetch_neurosynth(data_dir=f"{OUT}/ns_data", version="7", overwrite=False,
                            source="abstract", vocab="terms")
ns = ns_files[0] if isinstance(ns_files, list) else ns_files
dset_path = f"{OUT}/ns_data/neurosynth_dataset.pkl.gz"
if os.path.exists(dset_path):
    dset = Dataset.load(dset_path)
else:
    dset = convert_neurosynth_to_dataset(
        coordinates_file=ns["coordinates"], metadata_file=ns["metadata"],
        annotations_files=ns["features"])
    dset.save(dset_path)
print("    NS dataset: %d studies" % len(dset.ids), flush=True)

# resample HO atlas to the NS mask space once, for ROI sampling on term maps
def term_map(term):
    """MKDA Chi2 association-test z map for a term (studies w/ term vs without)."""
    feat = "terms_abstract_tfidf__" + term
    if feat not in dset.annotations.columns:
        # fallback: find a column containing the term
        cands = [c for c in dset.annotations.columns if c.endswith("__"+term)]
        if not cands:
            cands = [c for c in dset.annotations.columns if term in c]
        if not cands:
            raise KeyError(term)
        feat = cands[0]
    ids = dset.get_studies_by_label(feat, label_threshold=0.001)
    ids_other = list(set(dset.ids) - set(ids))
    d1 = dset.slice(ids); d2 = dset.slice(ids_other)
    mkda = MKDAChi2()
    res = mkda.fit(d1, d2)
    # association test z (zForward = P(term|activation)); use the consistency z map
    img = res.get_map("z_desc-association")
    return img, len(ids), feat

region_term = pd.DataFrame(index=list(REGION_MAP.keys()), columns=TERMS, dtype=float)
term_info = {}
for t in TERMS:
    try:
        zimg, nstud, feat = term_map(t)
    except Exception as e:
        print("    term FAIL %s: %r" % (t, e), flush=True); continue
    # resample atlas to this z map grid
    ho_rs = image.resample_to_img(ho_img, zimg, interpolation="nearest")
    ho_rs_data = np.asarray(ho_rs.dataobj)
    zdata = np.asarray(zimg.dataobj)
    for reg,(idxs,conf) in REGION_MAP.items():
        mask = np.isin(ho_rs_data, idxs)
        vals = zdata[mask]
        vals = vals[np.isfinite(vals)]
        region_term.loc[reg, t] = float(np.nanmean(vals)) if vals.size else np.nan
    term_info[t] = {"n_studies": int(nstud), "feature": feat}
    print("    term %-16s n=%-6d done" % (t, nstud), flush=True)

region_term.to_csv(f"{OUT}/region_term_neurosynth.tsv", sep="\t")
json.dump(term_info, open(f"{OUT}/term_info.json","w"), indent=2)
print("    saved region_term_neurosynth.tsv", flush=True)

# ------------------------------------------------------------------ DMN cross-check (Yeo-7 Default via Schaefer-400)
print("[3] DMN cross-check (Schaefer-400/Yeo-7 Default) ...", flush=True)
try:
    sch = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
    sch_img = image.load_img(sch["maps"])   # sch["maps"] is a path string
    labs = [l.decode() if isinstance(l,bytes) else l for l in sch.labels]
    dmn_idx = [i+1 for i,l in enumerate(labs) if "Default" in l]  # parcel ids 1-based
    sch_data = np.asarray(sch_img.dataobj)
    dmn_mask_img = image.new_img_like(sch_img, np.isin(sch_data, dmn_idx).astype("int8"))
    # mean neurosynth "default mode" z within Yeo-DMN vs HO-region overlap -> sanity only
    json.dump({"n_dmn_parcels": len(dmn_idx)}, open(f"{OUT}/dmn_meta.json","w"))
    print("    Yeo-DMN parcels:", len(dmn_idx), flush=True)
except Exception as e:
    print("    DMN cross-check skipped:", repr(e)[:120], flush=True)

# ------------------------------------------------------------------ program x term correlation
print("[4] program x term Spearman + permutation null ...", flush=True)
pz = pd.read_csv(f"{BASE}/results/crossregion_v1/program_region_zscore.tsv", sep="\t", index_col=0)
pz.columns = [str(c) for c in pz.columns]            # 14 regions x 60 programs

# Load names keyed by cnmf_component (= old_P int)
names_df = pd.read_csv(f"{BASE}/results/crossregion_v1/program_names.tsv", sep="\t")
names_df = names_df.set_index("cnmf_component")      # index = old_P int
name_short = names_df["name_short"].to_dict()         # {old_P_int: name_short}
conf_map   = names_df["confidence"].to_dict()          # {old_P_int: confidence}

regions = [r for r in REGION_MAP if r in pz.index and r in region_term.index]
RT = region_term.loc[regions].astype(float)          # 14 x term

# Filter excluded programs from PZ (data layer: use old_P columns as-is; drop excluded)
all_progs_old = [c for c in pz.columns]
kept_progs_old = [c for c in all_progs_old if int(c) not in EXCLUDED_OLD]
PZ = pz.loc[regions, kept_progs_old].astype(float)   # 14 x 54 (old_P cols, excluded removed)
progs = list(PZ.columns)  # old_P strings, kept only
print(f"[renumber] PZ programs kept: {len(progs)} (excluded {len(all_progs_old)-len(progs)})", flush=True)

# Flag check: verify no excluded_old in progs
flag_in_progs = [p for p in progs if int(p) in EXCLUDED_OLD]
if flag_in_progs:
    print(f"[FLAG] WARNING: excluded old_P still in progs: {flag_in_progs}", flush=True)
else:
    print(f"[renumber] FLAG check PASS: no excluded old_P in progs", flush=True)

# drop terms that are all-NaN
good_terms = [t for t in TERMS if RT[t].notna().sum() >= 6]
RT = RT[good_terms]

def spear_matrix(PZmat, RTmat):
    nP, nT = PZmat.shape[1], RTmat.shape[1]
    M = np.full((nP, nT), np.nan)
    for j,p in enumerate(PZmat.columns):
        x = PZmat[p].values
        for k,t in enumerate(RTmat.columns):
            y = RTmat[t].values
            ok = np.isfinite(x)&np.isfinite(y)
            if ok.sum()>=6:
                M[j,k] = spearmanr(x[ok], y[ok]).correlation
    return M

R = spear_matrix(PZ, RT)                              # 54 x nT observed
Rdf = pd.DataFrame(R, index=progs, columns=good_terms)
Rdf.to_csv(f"{OUT}/program_term_spearman.tsv", sep="\t")

# region-label permutation null: shuffle region rows of PZ, recompute, 10k
NPERM = 10000
perm_max = np.zeros((NPERM,))                         # for context
exceed = np.zeros_like(R)                             # |r_perm| >= |r_obs|
idx = np.arange(len(regions))
PZv = PZ.values; RTv = RT.values
# precompute ranks of RT columns (spearman = pearson on ranks); handle nan per pair simply by recompute
for b in range(NPERM):
    perm = RNG.permutation(idx)
    PZp = PZv[perm,:]
    # vectorized-ish spearman via per-column; acceptable at 54xnT, 10k
    for k in range(RTv.shape[1]):
        y = RTv[:,k]
        oky = np.isfinite(y)
        for j in range(PZp.shape[1]):
            x = PZp[:,j]
            ok = oky & np.isfinite(x)
            if ok.sum()>=6:
                rp = spearmanr(x[ok], y[ok]).correlation
                if abs(rp) >= abs(R[j,k]):
                    exceed[j,k]+=1
    if (b+1)%2000==0: print("    perm %d/%d" % (b+1,NPERM), flush=True)
pmat = (exceed+1)/(NPERM+1)
Pdf = pd.DataFrame(pmat, index=progs, columns=good_terms)
Pdf.to_csv(f"{OUT}/program_term_permp.tsv", sep="\t")

# BH-FDR across all 54 x nT
flat = pmat.flatten()
rej, q, _, _ = multipletests(flat, method="fdr_bh")
Qdf = pd.DataFrame(q.reshape(pmat.shape), index=progs, columns=good_terms)
Qdf.to_csv(f"{OUT}/program_term_fdr.tsv", sep="\t")

# top pairs table (display new_P)
recs=[]
for p in progs:
    op = int(p)
    nn = old_to_new[op]
    ns = name_short.get(op, "")
    cf = conf_map.get(op, "")
    for t in good_terms:
        recs.append((f"P{nn}", ns, cf, t,
                     Rdf.loc[p,t], Pdf.loc[p,t], Qdf.loc[p,t]))
top = pd.DataFrame(recs, columns=["program","name_short","confidence","term","spearman_r","perm_p","fdr_q"])
top = top.reindex(top.spearman_r.abs().sort_values(ascending=False).index)
top.to_csv(f"{OUT}/top_program_term_pairs.tsv", sep="\t", index=False)
print("    TOP pairs:\n", top.head(15).to_string(index=False), flush=True)

# DMN ranking
dmn_term = "default mode" if "default mode" in good_terms else good_terms[0]
dmn_rank = Rdf[dmn_term].sort_values(ascending=False)
dmn_rank.to_csv(f"{OUT}/program_DMN_ranking.tsv", sep="\t", header=[dmn_term])
print("    DMN top programs:\n", dmn_rank.head(8).to_string(), flush=True)

# ================================================================== FIGURE
print("[5] rendering figure ...", flush=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Use a shared sans-serif font stack across renderers.
import matplotlib as _mpl_font
_mpl_font.rcParams["font.family"] = "sans-serif"
_mpl_font.rcParams["font.sans-serif"] = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
_mpl_font.rcParams["pdf.fonttype"] = 42
_mpl_font.rcParams["ps.fonttype"] = 42
_mpl_font.rcParams["svg.fonttype"] = "none"
_mpl_font.rcParams["mathtext.fontset"] = "dejavusans"  # sans math
_mpl_font.rcParams["mathtext.default"] = "regular"  # math uses font.family
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
plt.rcParams.update({"font.size":6,"axes.linewidth":0.5,"font.family":"sans-serif","font.sans-serif":["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "xtick.major.size":2,"ytick.major.size":2,"xtick.major.width":0.4,"ytick.major.width":0.4,
    "pdf.fonttype":42,"ps.fonttype":42})

# cluster programs (rows) and terms (cols) on observed R
Rc = Rdf.fillna(0)
row_order = leaves_list(linkage(pdist(Rc.values), method="average")) if Rc.shape[0]>2 else range(Rc.shape[0])
col_order = leaves_list(linkage(pdist(Rc.values.T), method="average")) if Rc.shape[1]>2 else range(Rc.shape[1])
Rc = Rc.iloc[row_order, col_order]

# Build display labels using new_P numbering
def make_prog_lab(old_p_str, maxlen=34):
    op = int(old_p_str)
    nn = old_to_new.get(op)
    ns = name_short.get(op, "")
    if nn is None:
        return f"[FLAG-EXCLUDED-oldP{op}]"
    return f"P{nn} {ns}"[:maxlen]

prog_lab = [make_prog_lab(p) for p in Rc.index]
term_lab = list(Rc.columns)

fig = plt.figure(figsize=(180/25.4, 220/25.4), constrained_layout=True)  # 180mm wide composed
fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.10, wspace=0.04, hspace=0.07)
outer = fig.add_gridspec(3, 1, height_ratios=[3.2, 1.25, 1.55])
gtop = outer[0].subgridspec(1, 2, width_ratios=[1.55, 0.95], wspace=0.0)
gmid = outer[1].subgridspec(1, 3, wspace=0.04)
gbot = outer[2].subgridspec(1, 1)

def panel_tag(ax, s, dx=-30):
    ax.annotate(s, xy=(0,1), xycoords="axes fraction", xytext=(dx,8),
                textcoords="offset points", fontsize=10, fontweight="bold", va="bottom", ha="left")

# (a) program x term heatmap (clustered) — constrained_layout reserves room for row labels
axa = fig.add_subplot(gtop[0,0])
vmax = np.nanmax(np.abs(Rc.values)); vmax = max(vmax, 0.3)
im = axa.imshow(Rc.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
axa.set_xticks(range(len(term_lab))); axa.set_xticklabels(term_lab, rotation=45, ha="right", fontsize=5.5)
axa.set_yticks(range(len(prog_lab))); axa.set_yticklabels(prog_lab, fontsize=4.4)
axa.tick_params(axis="y", pad=1, length=1.5)
# mark FDR<0.1 cells
Qc = Qdf.reindex(index=Rc.index, columns=Rc.columns)
for i in range(Rc.shape[0]):
    for j in range(Rc.shape[1]):
        if Qc.values[i,j] < 0.1:
            axa.text(j,i,"*",ha="center",va="center",fontsize=6,color="k")
axa.set_title("Program x cognitive-term (Spearman; * FDR<0.1)", fontsize=6, pad=3)
cb=fig.colorbar(im, ax=axa, fraction=0.035, pad=0.015); cb.ax.tick_params(labelsize=5); cb.set_label("Spearman r", fontsize=5)
panel_tag(axa,"a",dx=-100)

# (b) DMN-highlight ranking (top + bottom programs by DMN r)
axb = fig.add_subplot(gtop[0,1])
dr = Rdf[dmn_term].dropna().sort_values()
sel = pd.concat([dr.head(8), dr.tail(8)])
ylab = [make_prog_lab(p, maxlen=30) for p in sel.index]
cols = ["#b2182b" if v>0 else "#2166ac" for v in sel.values]
axb.barh(range(len(sel)), sel.values, color=cols, height=0.78)
axb.set_yticks(range(len(sel))); axb.set_yticklabels(ylab, fontsize=4.8)
axb.set_ylim(-0.6, len(sel)-0.4)
axb.axvline(0,color="k",lw=0.4)
axb.set_xlabel(f"Spearman r vs '{dmn_term}'", fontsize=5.5)
axb.set_title("Programs ranked by\nDMN correlation", fontsize=6, pad=3)
panel_tag(axb,"b")

# (c) exemplar scatters: pick top-|r| pairs across distinct terms
ex = top.dropna(subset=["spearman_r"]).copy()
chosen=[]; seen=set()
for _,r in ex.iterrows():
    if r["term"] in seen: continue
    chosen.append(r); seen.add(r["term"])
    if len(chosen)==3: break
for ci,r in enumerate(chosen):
    axx = fig.add_subplot(gmid[0,ci])
    p_new = r["program"]  # already "P{nn}" string
    # find old_P for this new label to retrieve PZ column
    nn = int(p_new.lstrip("P"))
    op_match = [op for op, nv in old_to_new.items() if nv == nn]
    if not op_match:
        print(f"[FLAG] exemplar {p_new}: no old_P match, skipping", flush=True)
        continue
    op = op_match[0]
    old_p_str = str(op)
    t = r["term"]
    x=PZ[old_p_str].values; y=RT[t].values; ok=np.isfinite(x)&np.isfinite(y)
    axx.scatter(x[ok],y[ok],s=12,c="#444",zorder=3)
    for rg,xi,yi in zip(np.array(regions)[ok], x[ok], y[ok]):
        axx.annotate(rg,(xi,yi),fontsize=4.3,xytext=(2,2),textcoords="offset points")
    if ok.sum()>2:
        b=np.polyfit(x[ok],y[ok],1); xs=np.linspace(x[ok].min(),x[ok].max(),20)
        axx.plot(xs,np.polyval(b,xs),color="#b2182b",lw=0.8)
    ns_disp = name_short.get(op, "")
    axx.set_xlabel(f"P{nn} {ns_disp}"[:24]+" (region z)", fontsize=5)
    axx.set_ylabel(f"NS '{t}' z", fontsize=5)
    axx.set_title(f"r={r['spearman_r']:.2f}, q={r['fdr_q']:.2f}", fontsize=5.6)
    axx.tick_params(labelsize=4.5)
    if ci==0: panel_tag(axx,"c")

# (d) region x term Neurosynth atlas heatmap
axd = fig.add_subplot(gbot[0,0])
RTplot = region_term.loc[regions, good_terms].astype(float)
# cluster regions
rr = leaves_list(linkage(pdist(np.nan_to_num(RTplot.values)),"average")) if RTplot.shape[0]>2 else range(RTplot.shape[0])
RTp = RTplot.iloc[rr]
imd = axd.imshow(RTp.values, aspect="auto", cmap="viridis")
axd.set_xticks(range(len(good_terms))); axd.set_xticklabels(good_terms, rotation=45, ha="right", fontsize=5.5)
ylabs=[f"{r}{'*' if REGION_MAP[r][1]=='low' else ''}" for r in RTp.index]
axd.set_yticks(range(len(RTp.index))); axd.set_yticklabels(ylabs, fontsize=5.5)
axd.set_title("Neurosynth region x term cognitive atlas (mean z in ROI; *S1E low-conf)", fontsize=6, pad=3)
cb2=fig.colorbar(imd, ax=axd, fraction=0.02, pad=0.01); cb2.ax.tick_params(labelsize=5); cb2.set_label("mean z", fontsize=5)
panel_tag(axd,"d")

fig.suptitle("Fig. S10 | Program-cognition imaging-transcriptomics (Neurosynth DMN), exploratory analysis (n=14 regions)",
             fontsize=7)
fig.savefig(f"{FIGD}/figS10_cognition.pdf", dpi=400)
fig.savefig(f"{FIGD}/figS10_cognition.png", dpi=300)
fig.savefig("/tmp/figS10_cognition.png", dpi=300)
print("    saved figure", flush=True)
print("DONE", flush=True)
