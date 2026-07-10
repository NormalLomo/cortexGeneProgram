"""Python (matplotlib) panels:
 a) rep-chip tissue scatter colored by majorDomain
 b) 6-program small-multiples z-activity on rep chip
 f) program vs RCTD weight hexbin (3 pairs)
All vector PDF, points rasterized."""
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import PathCollection
import matplotlib.gridspec as gridspec

work="CORTEX_PROGRAM_ROOT/scripts/fig2/"
outd="CORTEX_PROGRAM_ROOT/figures/fig2/"

mpl.rcParams.update({
 "pdf.fonttype":42, "ps.fonttype":42, "font.family":"sans-serif",
 "font.sans-serif":["Arial","DejaVu Sans"], "font.size":7,
 "axes.linewidth":0.5, "axes.edgecolor":"#222222",
 "xtick.major.width":0.5,"ytick.major.width":0.5,
 "xtick.major.size":2.5,"ytick.major.size":2.5,
 "savefig.bbox":"tight","savefig.pad_inches":0.02,
})

LAYER_ORDER=["ARACHNOID","L1","L2","L3","L4","L5","L6","WM"]
# perceptually-ordered pia->WM palette (viridis-like but distinct, anatomical feel)
LAYER_COLORS={
 "ARACHNOID":"#9E9E9E","L1":"#3B4CC0","L2":"#5A78D6","L3":"#7DA0E0",
 "L4":"#36A66B","L5":"#E8A93B","L6":"#D65A3B","WM":"#7B3294"}

# z-activity diverging colormap (CNS style): blue-white-red, perceptual
ZCMAP=LinearSegmentedColormap.from_list("zdiv",
  ["#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B"])

choices=open(work+"_choices.txt").read()
REP=choices.split("REP=")[1].split("\n")[0].strip()
EX={}
for line in choices.splitlines():
    if line.startswith("EX\t"):
        _,k,v=line.split("\t"); EX[k]=v

# ---- program -> functional GO:BP label lookup (from program_names.tsv) ----
# label = "P{n} {name_short}"  (trailing " P{n}" stripped; "*" suffix if brain-weak)
_nm=pd.read_csv("CORTEX_PROGRAM_ROOT/results/crossregion_v1/program_names.tsv", sep="\t")
def _mk_label(n):
    r=_nm[_nm["program"]==int(n)]
    if len(r)==0: return f"P{n}"
    ns=str(r["name_short"].iloc[0]).strip()
    import re as _re
    ns=_re.sub(r"\s+P\d+$","",ns).strip()       # strip any trailing " P{n}"
    star="*" if str(r["confidence"].iloc[0])=="brain-weak" else ""
    return f"P{n} {ns}{star}"
PLABEL={int(r["program"]):_mk_label(r["program"]) for _,r in _nm.iterrows()}
def plabel(prog_key):
    """accept 'program_37' or '37' or 37 -> functional label."""
    n=int(str(prog_key).replace("program_",""))
    return PLABEL.get(n, f"P{n}")
# anatomical/cell-type context tag per exemplar (kept as spatial anchor)
EX_CONTEXT={"program_5":"L2/3","program_35":"L4","program_53":"L6",
            "program_37":"OLIGO / WM","program_56":"ENDO / vascular",
            "program_29":"Inhibitory"}
def ex_title(prog_key, sep=" · "):
    ctx=EX_CONTEXT.get(prog_key)
    lab=plabel(prog_key)
    return f"{ctx}{sep}{lab}" if ctx else lab

meta=pd.read_csv(work+"repchip_meta.tsv", sep="\t")
prog=pd.read_csv(work+"repchip_progscores.tsv", sep="\t")
df=meta.merge(prog, on="bin")
# flip y so pia/top reads naturally (keep aspect equal)
df["yy"]=df["y"].max()-df["y"]
print("rep chip bins for plot:", len(df))

# ---------------- PANEL A ----------------
fig,ax=plt.subplots(figsize=(3.2,3.2))
for ly in LAYER_ORDER:
    s=df[df.majorDomain==ly]
    ax.scatter(s["x"], s["yy"], s=1.2, c=LAYER_COLORS[ly], marker="s",
               linewidths=0, label=ly, rasterized=True)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title(f"Cortical layers (chip {REP})", fontsize=8, pad=4)
leg=ax.legend(loc="center left", bbox_to_anchor=(1.0,0.5), frameon=False,
              markerscale=4, handletextpad=0.3, labelspacing=0.35, fontsize=6.5,
              title="majorDomain", title_fontsize=7)
# scale bar 1 mm = 1000 um
x0=df["x"].min()+ (df["x"].max()-df["x"].min())*0.05
y0=df["yy"].min()+ (df["yy"].max()-df["yy"].min())*0.02
ax.plot([x0,x0+1000],[y0,y0], color="k", lw=1.5)
ax.text(x0+500,y0+ (df["yy"].max()-df["yy"].min())*0.02,"1 mm",ha="center",va="bottom",fontsize=6)
fig.savefig(outd+"fig2_a.pdf", dpi=450)
plt.close(fig); print("panel a done")

# ---------------- PANEL B ----------------
exprogs=list(EX.keys())
fig=plt.figure(figsize=(7.2,4.8))
gs=gridspec.GridSpec(2,3, figure=fig, wspace=0.16, hspace=0.26)
for i,p in enumerate(exprogs):
    ax=fig.add_subplot(gs[i//3, i%3])
    z=df[p].to_numpy()
    vlim=np.nanpercentile(np.abs(z),98)
    sc=ax.scatter(df["x"], df["yy"], c=z, s=1.0, marker="s", linewidths=0,
                  cmap=ZCMAP, norm=Normalize(-vlim,vlim), rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    # two-line title: anatomical context on top, functional GO:BP label below
    ax.set_title(ex_title(p, sep="\n"), fontsize=6.5, pad=2, linespacing=1.05)
    cb=fig.colorbar(sc, ax=ax, fraction=0.040, pad=0.02, shrink=0.7)
    cb.ax.tick_params(labelsize=5.5, width=0.4, length=1.8)
    cb.outline.set_linewidth(0.4)
    cb.set_label("z-activity", fontsize=5.5)
fig.suptitle(f"Spatial program activity (chip {REP})", fontsize=9, y=0.98)
fig.savefig(outd+"fig2_b.pdf", dpi=450)
plt.close(fig); print("panel b done")

# ---------------- PANEL F ----------------
fdat=pd.read_csv(work+"panelf_pairs.tsv", sep="\t")
PAIRS=[("program_37","OLIGO","OLIGO weight", plabel("program_37")),
       ("program_56","ENDO","ENDO weight", plabel("program_56")),
       ("program_5","L2-L3 IT LINC00507","L2/3 IT weight", plabel("program_5"))]
fig,axes=plt.subplots(1,3, figsize=(6.8,2.5))
for ax,(pp,cc,xl,tt) in zip(axes,PAIRS):
    x=fdat[cc].to_numpy(); y=fdat[pp].to_numpy()
    ok=np.isfinite(x)&np.isfinite(y)
    hb=ax.hexbin(x[ok], y[ok], gridsize=42, cmap="magma", bins="log",
                 mincnt=1, linewidths=0, rasterized=True)
    r=np.corrcoef(x[ok],y[ok])[0,1]
    ax.set_xlabel(xl, fontsize=7); ax.set_ylabel(f"{plabel(pp)} z", fontsize=6)
    ax.set_title(tt, fontsize=6.8)
    ax.text(0.04,0.93,f"r = {r:.2f}", transform=ax.transAxes, fontsize=7,
            va="top", bbox=dict(fc="white",ec="none",alpha=0.7,pad=1))
    ax.spines[["top","right"]].set_visible(False)
    cb=fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=5.5,width=0.4,length=1.8); cb.outline.set_linewidth(0.4)
    cb.set_label("log$_{10}$ bins", fontsize=5.5)
fig.suptitle("Program activity vs RCTD cell-type weight", fontsize=9, y=1.02)
fig.tight_layout()
fig.savefig(outd+"fig2_f.pdf", dpi=450)
plt.close(fig); print("panel f done")
print("PYTHON PANELS DONE")
