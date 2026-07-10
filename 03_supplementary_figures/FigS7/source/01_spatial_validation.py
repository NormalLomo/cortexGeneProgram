#!/usr/bin/env python
"""
Extended Data Fig 6 : spatial program-scoring validation (SCT).

Quality-control / method-validation figure for the SCT-corrected spatial program
projection. Demonstrates that the per-bin program scores (results/crossregion_v1/
spatial_bin50_program_score_SCT.parquet, depth-equalized SCT counts, zero per-bin
library normalization) land where the underlying biology predicts, are not a
sequencing-depth readout, and reproduce across chips. SCT data only.

Panels:
 (a) GAD ground-truth      : GAD1+GAD2 mean expression (counts / 10k UMI) per majorDomain
                             (GM L1-L6 vs WM vs ARACHNOID) -> validates the inhibitory
                             prior is gray-matter-high.
 (b) cell-class validation : per-majorDomain mean program-z (SCT) for a panel of programs
                             grouped by cell class -- inhibitory (P7/P26/P29) peak GM,
                             oligodendrocyte (P33/P41) peak WM, excitatory laminar (P8/P31)
                             peak their layer, microglia (P36) / astrocyte (P14) distributed.
 (c) depth control         : per-program correlation(per-domain mean z, per-domain mean UMI)
                             across the 60 programs (SCT) -- distribution near 0, scores are
                             not a depth artifact.
 (d) cross-chip repro      : per-program correlation of the 8-domain SCT profile between every
                             chip pair (14 cortical regions) -- distribution, median ~high.

SCT canonical = results/crossregion_v1/spatial_bin50_program_score_SCT.parquet
Writes figures/extended/ed_fig6_spatial_validation.{pdf,png}
"""
import os, glob, sys, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# FONT UNIFY (W-figfont-unify 2026-06-26): Nimbus Sans cross-engine
import matplotlib as _mpl_font
_mpl_font.rcParams["font.family"] = "sans-serif"
_mpl_font.rcParams["font.sans-serif"] = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
_mpl_font.rcParams["pdf.fonttype"] = 42
_mpl_font.rcParams["ps.fonttype"] = 42
_mpl_font.rcParams["svg.fonttype"] = "none"
_mpl_font.rcParams["mathtext.fontset"] = "dejavusans"  # sans math
_mpl_font.rcParams["mathtext.default"] = "regular"  # math uses font.family
from matplotlib.patches import Patch

BASE = "CORTEX_PROGRAM_ROOT"
RES  = f"{BASE}/results/crossregion_v1"
FIGD = f"{BASE}/figures/extended"
CACHE = f"{BASE}/scripts/extended/_ed6_cache"
SPATIAL = "CORTEX_PROGRAM_DATA_ROOT/neuropeptide_cortex/data/human/spatial/bin50"
os.makedirs(FIGD, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

SCT_PARQ = f"{RES}/spatial_bin50_program_score_SCT.parquet"

# majorDomain display order: gray matter laminae -> WM -> arachnoid
MD_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM", "ARACHNOID"]
GM = ["L1", "L2", "L3", "L4", "L5", "L6"]

# ---- program picks, grouped by cell class (verified from program_celltype_coloc.tsv) ----
# renumbered 2026-06-20: old_P -> new_P via results/crossregion_v1/program_renumber_map.tsv
# inhibitory interneuron programs (GM-high)
INH = [7, 26, 29]          # old[7,29,32]->new[7,26,29] | P7 Cation channel (interneuron), P26 Glutamate-receptor (postsyn.), P29 Ionotropic glutamate receptors
# oligodendrocyte / myelin programs (WM-high)
OLIGO = [33, 41]           # old[37,45]->new[33,41] | coloc OLIGO 0.74 / 0.77
# excitatory laminar programs (peak their layer)
EX = [8, 31]               # old[8,34]->new[8,31] | P8 Neurofilament cytoskeleton (pan-neuronal), P31 Activity-dependent IEG
MICRO = [36]               # old[40]->new[36] | Microglial immune activation (coloc MICRO top)
ASTRO = [14]               # old[15]->new[14] | Astrocyte glutamate transport (coloc AST 0.61)

# class -> (programs, expected-localization label)
CLASS_BLOCKS = [
    ("Inhibitory", INH,   "gray matter"),
    ("Oligo/myelin", OLIGO, "white matter"),
    ("Excitatory", EX,    "its layer"),
    ("Microglia", MICRO,  "distributed"),
    ("Astrocyte", ASTRO,  "distributed"),
]


def progcol(p): return f"program_{p}"


# ======================================================================
# STAGE 1 : build small per-majorDomain summaries from the heavy parquet
# ======================================================================
def domain_means(parquet, tag):
    """per-majorDomain mean of every program-z + mean UMI. cached."""
    out = f"{CACHE}/domain_mean_{tag}.tsv"
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t", index_col=0)
    import pyarrow.parquet as pq
    print(f"[stage1] domain_means {tag}: reading parquet ...", flush=True)
    pcols = [progcol(i) for i in range(1, 61)]
    df = pq.read_table(parquet, columns=["majorDomain", "bin_total_umi"] + pcols).to_pandas()
    g = df.groupby("majorDomain", observed=True)
    m = g[pcols].mean()
    m["bin_total_umi"] = g["bin_total_umi"].mean()
    m["n_bins"] = g.size()
    m = m.reindex(MD_ORDER)
    m.to_csv(out, sep="\t")
    print(f"[stage1] wrote {out}  shape={m.shape}", flush=True)
    return m


def domain_means_byregion(parquet, tag):
    """per-(region, majorDomain) mean program-z, long table. cached.
    Used for the cross-chip reproducibility panel."""
    out = f"{CACHE}/domain_mean_byregion_{tag}.tsv"
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t")
    import pyarrow.parquet as pq
    print(f"[stage1] domain_means_byregion {tag}: reading parquet ...", flush=True)
    pcols = [progcol(i) for i in range(1, 61)]
    df = pq.read_table(parquet, columns=["region", "majorDomain"] + pcols).to_pandas()
    g = df.groupby(["region", "majorDomain"], observed=True)[pcols].mean().reset_index()
    g.to_csv(out, sep="\t", index=False)
    print(f"[stage1] wrote {out}  shape={g.shape}", flush=True)
    return g


def _h5_obs(om, k):
    """read an obs column from an open h5ad obs group via h5py (handles categorical)."""
    import h5py
    o = om[k]
    if isinstance(o, h5py.Group):  # categorical
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in o["categories"][:]])
        return cats[o["codes"][:]]
    return o[:]


def gad_ground_truth():
    """GAD1+GAD2 CP10k mean per majorDomain, pooled across all 44 bin50 chips. cached.
    Pure-h5py reader (anndata backed open was ~110s/chip on these CSR layers; h5py ~4s/chip)."""
    out = f"{CACHE}/gad_by_majordomain.tsv"
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t", index_col=0)
    import h5py
    import scipy.sparse as sp
    files = sorted(glob.glob(f"{SPATIAL}/*_bin50.h5ad"))
    sums = {md: 0.0 for md in MD_ORDER}
    cnts = {md: 0 for md in MD_ORDER}
    for i, fp in enumerate(files):
        f = h5py.File(fp, "r")
        vg = f["var"]
        idxkey = vg.attrs.get("_index", "_index")
        vnames = np.array([x.decode() if isinstance(x, bytes) else x for x in vg[idxkey][:]])
        if not ({"GAD1", "GAD2"} <= set(vnames.tolist())):
            f.close(); continue
        gi = [int(np.where(vnames == g)[0][0]) for g in ["GAD1", "GAD2"]]
        grp = f["layers/counts"]
        M = sp.csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]),
                          shape=tuple(grp.attrs["shape"]))
        gad = np.asarray(M[:, gi].sum(axis=1)).ravel()        # GAD1+GAD2 raw count per bin
        om = f["obs"]
        md = np.array([x.decode() if isinstance(x, bytes) else x for x in _h5_obs(om, "majorDomain")])
        tot = np.asarray(_h5_obs(om, "bin_total_umi"), dtype=float)
        tot[tot == 0] = np.nan
        cp10k = gad / tot * 1e4                                # per 10k UMI per bin
        for k in MD_ORDER:
            m = md == k
            if m.any():
                v = cp10k[m]; v = v[np.isfinite(v)]
                sums[k] += float(v.sum()); cnts[k] += int(v.size)
        f.close()
        if (i + 1) % 10 == 0:
            print(f"[stage1 gad] {i+1}/{len(files)} chips", flush=True)
    res = pd.DataFrame({
        "gad_cp10k_mean": [sums[k] / cnts[k] if cnts[k] else np.nan for k in MD_ORDER],
        "n_bins": [cnts[k] for k in MD_ORDER],
    }, index=MD_ORDER)
    res.to_csv(out, sep="\t")
    print(f"[stage1] wrote {out}")
    return res


print("=== STAGE 1: building summaries (cached) ===")
gad = gad_ground_truth()
dm_sct = domain_means(SCT_PARQ, "SCT")
dmr_sct = domain_means_byregion(SCT_PARQ, "SCT")


# per-BIN corr(program score, bin_total_umi) for all 60 programs (SCT, panel c).
# This is the correct granularity for the depth-artifact test: a per-program score that
# merely read out sequencing depth would correlate strongly with per-bin UMI. Computed
# over all ~5.7M bins (heavy) and cached.
def binlevel_umi_corr(parquet, tag):
    out = f"{CACHE}/binlevel_umi_corr_{tag}.tsv"
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t")
    import pyarrow.parquet as pq
    from scipy.stats import spearmanr
    print(f"[stage1] binlevel_umi_corr {tag}: reading parquet (per-bin) ...", flush=True)
    pcols = [progcol(i) for i in range(1, 61)]
    t = pq.read_table(parquet, columns=["bin_total_umi"] + pcols).to_pandas()
    u = t["bin_total_umi"].values.astype(float)
    rows = []
    for i in range(1, 61):
        z = t[progcol(i)].values.astype(float)
        ok = np.isfinite(z) & np.isfinite(u)
        pear = np.corrcoef(z[ok], u[ok])[0, 1]
        sp = spearmanr(z[ok], u[ok]).correlation
        rows.append((i, pear, sp))
    d = pd.DataFrame(rows, columns=["program", "pearson_umi", "spearman_umi"])
    d.to_csv(out, sep="\t", index=False)
    print(f"[stage1] wrote {out}", flush=True)
    return d


binc = binlevel_umi_corr(SCT_PARQ, "SCT")
corr_sct = binc["pearson_umi"].values.astype(float)   # per-bin Pearson, 60 programs
print(f"median per-bin corr( score , UMI ) across 60 programs  SCT={np.nanmedian(corr_sct):+.3f}")


# per-program cross-chip reproducibility (panel d):
# for each program, build region x domain matrix; correlate every region pair's
# 8-domain profile; take the mean pairwise correlation as that program's repro score.
def crosschip_repro(dmr):
    regions = sorted(dmr["region"].unique())
    perprog_mean = []
    all_pair = []  # flattened pairwise corrs across all programs, for histogram
    for i in range(1, 61):
        pc = progcol(i)
        # region x domain (ordered)
        mat = (dmr.pivot(index="region", columns="majorDomain", values=pc)
                  .reindex(index=regions, columns=MD_ORDER))
        vals = mat.values  # n_region x 8
        pcs = []
        for a, b in itertools.combinations(range(len(regions)), 2):
            va, vb = vals[a], vals[b]
            ok = np.isfinite(va) & np.isfinite(vb)
            if ok.sum() > 2 and np.std(va[ok]) > 0 and np.std(vb[ok]) > 0:
                r = np.corrcoef(va[ok], vb[ok])[0, 1]
                pcs.append(r); all_pair.append(r)
        perprog_mean.append(np.nanmean(pcs) if pcs else np.nan)
    return np.array(perprog_mean), np.array(all_pair), regions


repro_perprog, repro_pairs, REGIONS = crosschip_repro(dmr_sct)
print(f"cross-chip repro: {len(REGIONS)} regions, median per-program profile r="
      f"{np.nanmedian(repro_perprog):+.3f}  (n pairwise={len(repro_pairs)})")


# ======================================================================
# STAGE 2 : figure (native matplotlib, vector PDF, project style)
# ======================================================================
plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "figure.dpi": 200,
})

# palette
C_GM  = "#4C72B0"   # gray matter laminae
C_WM  = "#B2182B"   # white matter
C_AR  = "#8172B3"   # arachnoid
C_SCT = "#2166AC"   # SCT accent
PNAMES = pd.read_csv(f"{RES}/program_names.tsv", sep="\t", index_col=0)


def md_colors(idx):
    return [C_WM if k == "WM" else C_AR if k == "ARACHNOID" else C_GM for k in idx]


def pshort(p):
    # cln-N3 fix 2026-06-25: PNAMES index is 'P7' string (new_P format),
    # but caller passes integer p; must convert to 'P{p}' to avoid KeyError
    # which previously caused fallback f'P{p}' → double-label 'P7 P7'.
    key = f"P{p}"
    s = PNAMES.loc[key, "name_short"] if key in PNAMES.index else f"P{p}"
    return s if len(s) <= 22 else s[:20] + "…"   # truncate long names for the heatmap rows


fig = plt.figure(figsize=(7.6, 6.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.08],
                      hspace=0.58, wspace=0.62,
                      left=0.080, right=0.965, top=0.93, bottom=0.085)

# ---------------- (a) GAD ground-truth ----------------
axa = fig.add_subplot(gs[0, 0])
xa = np.arange(len(MD_ORDER))
va = gad["gad_cp10k_mean"].values
axa.bar(xa, va, 0.72, color=md_colors(MD_ORDER), edgecolor="black", linewidth=0.5)
axa.set_xticks(xa); axa.set_xticklabels(MD_ORDER, rotation=45, ha="right")
axa.set_ylabel("GAD1+GAD2  (counts / 10k UMI)")
axa.set_title("a  Inhibitory ground truth: GAD is gray-matter", loc="left", fontweight="bold")
axa.annotate("", xy=(-0.3, max(va)*1.02), xytext=(5.3, max(va)*1.02),
             arrowprops=dict(arrowstyle="-", lw=0.8, color="0.35"))
axa.text(2.5, max(va)*1.07, "gray matter", ha="center", va="bottom", fontsize=6, color="0.35")
axa.set_ylim(0, max(va)*1.20)
axa.legend(handles=[Patch(fc=C_GM, ec="k", lw=.4, label="GM laminae"),
                    Patch(fc=C_WM, ec="k", lw=.4, label="white matter"),
                    Patch(fc=C_AR, ec="k", lw=.4, label="arachnoid")],
           loc="upper right", handlelength=1.1)

# ---------------- (b) cell-class spatial validation (SCT, small heatmap) ----------------
axb = fig.add_subplot(gs[0, 1])
prog_rows, row_labels, block_spans = [], [], []
r0 = 0
for cname, progs, _exp in CLASS_BLOCKS:
    for p in progs:
        prog_rows.append(dm_sct[progcol(p)].reindex(MD_ORDER).values.astype(float))
        row_labels.append(f"P{p} {pshort(p)}")
    block_spans.append((cname, r0, r0 + len(progs)))
    r0 += len(progs)
H = np.vstack(prog_rows)                         # n_prog x 8, already mean-z (global z per program)
vmax = np.nanmax(np.abs(H))
im = axb.imshow(H, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
axb.set_xticks(np.arange(len(MD_ORDER))); axb.set_xticklabels(MD_ORDER, rotation=45, ha="right")
axb.set_yticks(np.arange(len(row_labels))); axb.set_yticklabels(row_labels, fontsize=5.6)
# class block dividers (cell-class identity is already given by the row names)
for cname, a, b in block_spans:
    if a > 0:
        axb.axhline(a - 0.5, color="k", lw=0.8)
axb.set_title("b  Programs land where biology expects (SCT)",
              loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=axb, fraction=0.045, pad=0.02)
cb.set_label("mean z", fontsize=6); cb.ax.tick_params(labelsize=5.5)
cb.outline.set_linewidth(0.5)

# ---------------- (c) depth control (SCT only, per-bin) ----------------
axc = fig.add_subplot(gs[1, 0])
med_c = np.nanmedian(corr_sct)
hi_c = float(np.nanmax(corr_sct))
bins = np.linspace(-0.1, max(0.6, hi_c * 1.05), 26)
axc.hist(corr_sct, bins=bins, color=C_SCT, alpha=0.82, edgecolor="black",
         linewidth=0.3)
axc.axvline(0, color="0.4", lw=0.6, ls="--")
axc.axvline(med_c, color="#A1352B", lw=1.0, label=f"median = {med_c:.2f}")
axc.set_xlabel("per-bin corr( program score , bin UMI )  (54 programs)")
axc.set_ylabel("# programs")
axc.set_title("c  Depth control: scores do not track sequencing depth",
              loc="left", fontweight="bold")
axc.set_xlim(-0.12, bins[-1] + 0.02)
axc.legend(loc="upper right", handlelength=1.1)

# ---------------- (d) cross-chip reproducibility (SCT) ----------------
axd = fig.add_subplot(gs[1, 1])
bins2 = np.linspace(-1, 1, 25)
med_d = np.nanmedian(repro_perprog)
axd.hist(repro_perprog, bins=bins2, color="#2A8A5B", alpha=0.85, edgecolor="black",
         linewidth=0.3)
axd.axvline(0, color="0.4", lw=0.6, ls="--")
axd.axvline(med_d, color="#0B5A37", lw=1.0, label=f"median = {med_d:+.2f}")
axd.set_xlabel("cross-chip corr of 8-domain profile  (mean over region pairs)")
axd.set_ylabel("# programs")
axd.set_title(f"d  Cross-chip reproducibility ({len(REGIONS)} regions)",
              loc="left", fontweight="bold")
axd.set_xlim(-1.05, 1.05)
axd.legend(loc="upper left", handlelength=1.1)

# super-caption
fig.text(0.5, 0.006,
         "SCT-corrected spatial program scores reproduce the expected laminar/white-matter biology, "
         "are decoupled from sequencing depth, and are reproducible across cortical chips.",
         ha="center", va="bottom", fontsize=5.8, color="0.30")

out_pdf = f"{FIGD}/ed_fig6_spatial_validation.pdf"
out_png = f"{FIGD}/ed_fig6_spatial_validation.png"
fig.savefig(out_pdf)
fig.savefig(out_png, dpi=300)
# local-pull copy for the user
fig.savefig("/tmp/ed_fig6_spatial_qc.png", dpi=300)
print("WROTE", out_pdf, out_png, "/tmp/ed_fig6_spatial_qc.png")

# ======================================================================
# REPORT NUMBERS
# ======================================================================
print("\n=== REPORT NUMBERS ===")
gm_gad = np.nanmean([gad.loc[k, "gad_cp10k_mean"] for k in GM])
print("GAD CP10k by domain:", dict(zip(MD_ORDER, np.round(gad["gad_cp10k_mean"].values, 3))))
print("GAD: GM-mean=%.3f  WM=%.3f  ARACHNOID=%.3f"
      % (gm_gad, gad.loc["WM", "gad_cp10k_mean"], gad.loc["ARACHNOID", "gad_cp10k_mean"]))


def peak_domain(p):
    v = dm_sct[progcol(p)].reindex(MD_ORDER).values.astype(float)
    return MD_ORDER[int(np.nanargmax(v))]


def peak_domain_gm(p):
    v = dm_sct[progcol(p)].reindex(GM).values.astype(float)
    return GM[int(np.nanargmax(v))]


print("INH peak domains:  " + ", ".join(f"P{p}->{peak_domain(p)}" for p in INH))
print("OLIGO peak domains:" + ", ".join(f"P{p}->{peak_domain(p)}" for p in OLIGO))
print("EX peak domains:   " + ", ".join(f"P{p}->{peak_domain(p)} (GM-arg {peak_domain_gm(p)})" for p in EX))
print("MICRO peak: " + ", ".join(f"P{p}->{peak_domain(p)}" for p in MICRO)
      + "  ASTRO peak: " + ", ".join(f"P{p}->{peak_domain(p)}" for p in ASTRO))
print(f"SCT per-bin score<->UMI median corr (60 programs) = {np.nanmedian(corr_sct):+.3f}  "
      f"[IQR {np.nanpercentile(corr_sct,25):+.3f},{np.nanpercentile(corr_sct,75):+.3f}]  "
      f"max={np.nanmax(corr_sct):+.3f}")
print(f"cross-chip repro median per-program r = {np.nanmedian(repro_perprog):+.3f}  "
      f"over {len(REGIONS)} regions ({len(repro_pairs)} pairwise corrs)")
print("fig size inches:", fig.get_size_inches())
