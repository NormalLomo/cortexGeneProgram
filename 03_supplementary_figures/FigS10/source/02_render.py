#!/usr/bin/env python
"""Render Fig. S10 from canonical cognition matrices with the retained-program layout.

Program labels map from old cNMF component numbers (1-60) to new_P numbers
(P1-P54) through program_renumber_map.tsv. Excluded programs
(old 9/18/19/35/52/57) are omitted from all panels, and name_short is read
through cnmf_component.

Program map format:
  program_renumber_map.tsv: old_P (int), new_P (int or "EXCLUDED")
  program_names.tsv: new_P in format "P1"..."P54" or "EXCLUDED"
  Mapping: old cNMF int -> new int -> "P{new_int}" for name lookup
"""
import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd
BASE="CORTEX_PROGRAM_ROOT"
OUT=f"{BASE}/results/crossregion_v1/program_cognition"; FIGD=f"{BASE}/figures/fig7"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
# Use a 180 mm portrait page with Liberation Sans and a 5 pt minimum font size.
FS_MIN = 5.0
matplotlib.rcParams["pdf.fonttype"] = 42      # embed TrueType, no Type3
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
plt.rcParams.update({"font.size":6,"axes.titlesize":7,"axes.labelsize":6.2,
    "xtick.labelsize":5.4,"ytick.labelsize":5.4,"legend.fontsize":5.4,
    "axes.linewidth":0.5,"axes.edgecolor":"#333333",
    "xtick.major.size":1.8,"ytick.major.size":1.8,"xtick.major.width":0.5,"ytick.major.width":0.5})

# ================================================================== RENUMBER MAP
# Load renumber map: old_P (cNMF component int) -> new_P (int or "EXCLUDED")
renum = pd.read_csv(f"{BASE}/results/crossregion_v1/program_renumber_map.tsv", sep="\t")
# old2new: str(old_cNMF_int) -> "P{new_int}" or "EXCLUDED"
old2new = {}
for _, row in renum.iterrows():
    old_str = str(int(row["old_P"]))
    new_val = str(row["new_P"])
    if new_val.upper() == "EXCLUDED":
        old2new[old_str] = "EXCLUDED"
    else:
        # new_P is integer (e.g. 36) -> convert to "P36"
        try:
            old2new[old_str] = f"P{int(float(new_val))}"
        except ValueError:
            old2new[old_str] = "EXCLUDED"

excluded_old = {k for k, v in old2new.items() if v == "EXCLUDED"}
print(f"[renumber] excluded old cNMF: {sorted(excluded_old, key=int)}", flush=True)

# name_short / confidence lookup via program_names.tsv
# Columns: new_P ("P1"..."P54" or "EXCLUDED"), cnmf_component (int), name_short, confidence
names = pd.read_csv(f"{BASE}/results/crossregion_v1/program_names.tsv", sep="\t")
name_short = {}   # "P{n}" -> name_short string
conf_map   = {}   # "P{n}" -> confidence string
for _, row in names.iterrows():
    np_str = str(row["new_P"])   # e.g. "P1", "P36", "EXCLUDED"
    name_short[np_str] = str(row["name_short"])
    conf_map[np_str]   = str(row["confidence"])
print(f"[renumber] name_short loaded for {len(name_short)} entries (sample P36: {name_short.get('P36','?')})", flush=True)

def to_newP(old_str):
    """Return 'P{new_int}' or None if excluded."""
    v = old2new.get(str(old_str))
    if v is None or v == "EXCLUDED":
        return None
    return v

def newP_label(np_str, maxlen=34):
    """Format display label: 'P36 Microglial immune activation' (+ '*' if brain-weak)"""
    sname = name_short.get(np_str, "")
    star = "*" if conf_map.get(np_str, "") == "brain-weak" else ""
    return f"{np_str} {sname}{star}"[:maxlen]

# ================================================================== LOAD CACHED MATRICES
Rdf = pd.read_csv(f"{OUT}/program_term_spearman.tsv", sep="\t", index_col=0); Rdf.index=Rdf.index.astype(str)
Qdf = pd.read_csv(f"{OUT}/program_term_fdr.tsv", sep="\t", index_col=0); Qdf.index=Qdf.index.astype(str)
region_term = pd.read_csv(f"{OUT}/region_term_neurosynth.tsv", sep="\t", index_col=0)
top = pd.read_csv(f"{OUT}/top_program_term_pairs.tsv", sep="\t"); top["program"]=top["program"].astype(str)
pz = pd.read_csv(f"{BASE}/results/crossregion_v1/program_region_zscore.tsv", sep="\t", index_col=0)
pz.columns=[str(c) for c in pz.columns]

# ================================================================== APPLY RENUMBER
# 1. Drop EXCLUDED from Rdf / Qdf (index = old cNMF str)
keep_idx = [i for i in Rdf.index if to_newP(i) is not None]
Rdf = Rdf.loc[keep_idx].copy()
Qdf = Qdf.loc[keep_idx].copy()

# 2. Rename index old str -> "P{new_int}"
Rdf.index = [to_newP(i) for i in Rdf.index]
Qdf.index = [to_newP(i) for i in Qdf.index]

# 3. Drop EXCLUDED from pz columns, rename
keep_cols = [c for c in pz.columns if to_newP(c) is not None]
pz = pz[keep_cols].copy()
pz.columns = [to_newP(c) for c in pz.columns]

# 4. Remap top_program_term_pairs: drop excluded, remap program col
top["_new_P"] = top["program"].apply(lambda x: to_newP(str(x)))
top = top[top["_new_P"].notna()].copy()
top["program"] = top["_new_P"]
top = top.drop(columns=["_new_P"])
# refresh name_short and confidence from new_P
top["name_short"] = top["program"].apply(lambda p: name_short.get(p, ""))
top["confidence"]  = top["program"].apply(lambda p: conf_map.get(p, ""))

print(f"[renumber] Rdf shape after drop+remap: {Rdf.shape}", flush=True)
print(f"[renumber] Rdf index sample (first 5): {list(Rdf.index[:5])}", flush=True)
print(f"[renumber] pz columns count: {len(pz.columns)}", flush=True)
print(f"[renumber] top rows: {len(top)}", flush=True)
print(f"[renumber] label check P36: {newP_label('P36')}", flush=True)

# ================================================================== FIGURE SETUP
REGION_MAP_CONF = {"FPPFC":"high","DLPFC":"high","VLPFC":"high","M1":"high","ACC":"high","S1":"high",
 "PoCG":"high-collinear","S1E":"low","STG":"high","ITG":"high","SMG":"high","SPL":"high","AG":"high","V1":"high"}
good_terms = list(Rdf.columns)
regions = [r for r in REGION_MAP_CONF if r in pz.index and r in region_term.index]
RT = region_term.loc[regions, good_terms].astype(float)
PZ = pz.loc[regions].astype(float)
dmn_term = "default mode" if "default mode" in good_terms else good_terms[0]

Rc = Rdf.fillna(0)
row_order = leaves_list(linkage(pdist(Rc.values),"average")); col_order = leaves_list(linkage(pdist(Rc.values.T),"average"))
Rc = Rc.iloc[row_order, col_order]
# prog_lab now uses new_P index (e.g. "P36 Microglial immune activation")
prog_lab = [newP_label(p, 34) for p in Rc.index]
term_lab = list(Rc.columns)

FW = 180.0/25.4                 # 7.087 in == 180 mm wide
FH = 224.0/25.4                 # 8.819 in == 224 mm tall  -> AR(w/h)=0.804 portrait
fig = plt.figure(figsize=(FW, FH), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.10, wspace=0.04, hspace=0.07)
outer = fig.add_gridspec(3,1,height_ratios=[3.35,1.25,1.55])
gtop=outer[0].subgridspec(1,2,width_ratios=[1.5,1.0],wspace=0.06)
gmid=outer[1].subgridspec(1,3,wspace=0.30); gbot=outer[2].subgridspec(1,1)
def tag(ax,s,dx=-0.06,dy=1.015): ax.text(dx,dy,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom",ha="right")

# ================================================================== PANEL a: program x term heatmap
axa=fig.add_subplot(gtop[0,0])
vmax=max(np.nanmax(np.abs(Rc.values)),0.3)
im=axa.imshow(Rc.values,aspect="equal",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
axa.set_xticks(range(len(term_lab))); axa.set_xticklabels(term_lab,rotation=45,ha="right",fontsize=5.4)
axa.set_yticks(range(len(prog_lab))); axa.set_yticklabels(prog_lab,fontsize=5.0); axa.tick_params(axis="y",pad=1,length=1.5)
Qc=Qdf.reindex(index=Rc.index,columns=Rc.columns)
for i in range(Rc.shape[0]):
    for j in range(Rc.shape[1]):
        if Qc.values[i,j]<0.1: axa.text(j,i,"*",ha="center",va="center",fontsize=5.5,color="k")
axa.set_title("Program x cognitive-term (Spearman; * FDR<0.1)",fontsize=6,pad=3)
cb=fig.colorbar(im,ax=axa,fraction=0.030,pad=0.02); cb.ax.tick_params(labelsize=5); cb.set_label("Spearman r",fontsize=5); tag(axa,"a",dx=-0.30,dy=1.085)

# ================================================================== PANEL b: DMN ranking
axb=fig.add_subplot(gtop[0,1])
dr=Rdf[dmn_term].dropna().sort_values(); sel=pd.concat([dr.head(8),dr.tail(8)])
ylab=[newP_label(p, 30) for p in sel.index]
cols=["#b2182b" if v>0 else "#2166ac" for v in sel.values]
axb.barh(range(len(sel)),sel.values,color=cols,height=0.78); axb.set_yticks(range(len(sel))); axb.set_yticklabels(ylab,fontsize=5.0)
axb.set_ylim(-0.6,len(sel)-0.4); axb.axvline(0,color="k",lw=0.4); axb.set_xlabel(f"Spearman r vs '{dmn_term}'",fontsize=5.5)
axb.tick_params(axis="y",pad=1,length=1.5); axb.tick_params(axis="x",labelsize=5.0)
axb.set_title("Programs ranked by\nDMN correlation",fontsize=6,pad=3); tag(axb,"b",dx=-0.42)

# ================================================================== PANEL c: exemplar scatters
ex=top.dropna(subset=["spearman_r"]).copy(); chosen=[]; seen=set()
for _,r in ex.iterrows():
    if r["term"] in seen: continue
    chosen.append(r); seen.add(r["term"])
    if len(chosen)==3: break
for ci,r in enumerate(chosen):
    axx=fig.add_subplot(gmid[0,ci]); p=r["program"]; t=r["term"]
    x=PZ[p].values; y=RT[t].values; ok=np.isfinite(x)&np.isfinite(y)
    axx.scatter(x[ok],y[ok],s=12,c="#444",zorder=3)
    for rg,xi,yi in zip(np.array(regions)[ok],x[ok],y[ok]): axx.annotate(rg,(xi,yi),fontsize=5.0,xytext=(2,2),textcoords="offset points")
    if ok.sum()>2:
        b=np.polyfit(x[ok],y[ok],1); xs=np.linspace(x[ok].min(),x[ok].max(),20); axx.plot(xs,np.polyval(b,xs),color="#b2182b",lw=0.8)
    axx.set_xlabel(f"{newP_label(p,24)} (region z)",fontsize=5.2); axx.set_ylabel(f"NS '{t}' z",fontsize=5.2)
    axx.set_title(f"r={r['spearman_r']:.2f}, q={r['fdr_q']:.2f}",fontsize=5.8); axx.tick_params(labelsize=5.0)
    axx.set_box_aspect(1.0)
    if ci==0: tag(axx,"c",dx=-0.22)

# ================================================================== PANEL d: region x term atlas heatmap
axd=fig.add_subplot(gbot[0,0])
RTp=RT.copy(); rr=leaves_list(linkage(pdist(np.nan_to_num(RTp.values)),"average")); RTp=RTp.iloc[rr]
imd=axd.imshow(RTp.values,aspect="equal",cmap="viridis")
axd.set_xticks(range(len(good_terms))); axd.set_xticklabels(good_terms,rotation=45,ha="right",fontsize=5.5)
ylabs=[f"{r}{'*' if REGION_MAP_CONF[r]=='low' else ''}" for r in RTp.index]
axd.set_yticks(range(len(RTp.index))); axd.set_yticklabels(ylabs,fontsize=5.5); axd.tick_params(axis="y",pad=1,length=1.5)
axd.set_title("Neurosynth region x term cognitive atlas\n(mean z in ROI; *S1E low-conf)",fontsize=6,pad=3)
cb2=fig.colorbar(imd,ax=axd,fraction=0.018,pad=0.01); cb2.ax.tick_params(labelsize=5); cb2.set_label("mean z",fontsize=5); tag(axd,"d",dx=-0.085,dy=1.16)

fig.suptitle("Fig. S10 | Gene-program association with cognitive functional architecture (Neurosynth, n=14 regions)", fontsize=7)
fig.savefig(f"{FIGD}/figS10_cognition.pdf",dpi=400); fig.savefig(f"{FIGD}/figS10_cognition.png",dpi=300)
fig.savefig(f"{FIGD}/figS10_full.pdf",dpi=400); fig.savefig(f"{FIGD}/figS10_full.png",dpi=320)
fig.savefig("/tmp/figS10_full.png",dpi=320)
fig.savefig("/tmp/figS10_full.svg")
print("RENDER_DONE")
