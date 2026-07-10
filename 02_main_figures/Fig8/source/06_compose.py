#!/usr/bin/env python3
"""Compact candidate Fig. 8: disease and aging relevance of cortical programs.

This composition writes only to the public figure output directory and does not
overwrite the current SUBMIT Fig. 8.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, Patch, PathPatch, Wedge
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from statsmodels.stats.multitest import multipletests


ROOT = Path("__PRIVATE_CANONICAL_ROOT__")
PD_DIR = ROOT / "results/crossregion_v1/program_disease"
AGE_DIR = ROOT / "results/crossregion_v1/program_aging"
LOAD = ROOT / "results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv"
NAMES = ROOT / "results/crossregion_v1/program_names.tsv"
RENUM = ROOT / "results/crossregion_v1/program_renumber_map.tsv"
RZ = ROOT / "results/crossregion_v1/program_region_zscore.tsv"
BB = ROOT / "data/brainbase/disease_gene_associations.txt"
OUT = ROOT / "figures/Fig8/outputs"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY_N = 150
FDR_SIG = 0.05
COHORT_TECHNICAL = {9, 18, 19, 35, 52, 57}

FS_MIN = 5.0

matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "Arial", "DejaVu Sans"],
        "font.size": 6,
        "axes.titlesize": 7,
        "axes.labelsize": 6.2,
        "xtick.labelsize": 5.2,
        "ytick.labelsize": 5.2,
        "legend.fontsize": 5.0,
        "axes.linewidth": 0.45,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
        "xtick.major.size": 1.5,
        "ytick.major.size": 1.5,
        "axes.edgecolor": "#333333",
        "pdf.use14corefonts": False,
    }
)

CAT_ORDER = ["Neurodegenerative", "Psychiatric/Neurodev", "Tumor", "Developmental", "Other/Vascular"]
CAT_COL = {
    "Neurodegenerative": "#B4493F",
    "Psychiatric/Neurodev": "#3B74A8",
    "Tumor": "#8064A2",
    "Developmental": "#4E9A72",
    "Other/Vascular": "#C9823B",
}

LIN_COL = {
    "Microglia": "#8A65A8",
    "Astrocyte": "#4E9A72",
    "Oligo/OPC": "#C9823B",
    "Vascular": "#B75A58",
    "Inhibitory": "#2B6CB0",
    "Excitatory": "#4A5568",
}

AGE_SET_LABEL = {
    "HAGR_GenAge_human": "GenAge",
    "HAGR_CellAge_induces_senescence": "CellAge+",
    "HAGR_CellAge_inhibits_senescence": "CellAge-",
    "SAUL_SEN_MAYO": "SenMayo",
    "LIU2026_CELL_ACCELERATED_AGING_PROTEINS": "AA proteins",
    "LIU2026_CELL_DECELERATED_AGING_PROTEINS": "DA proteins",
}

AGE_SET_ORDER = [
    "LIU2026_CELL_ACCELERATED_AGING_PROTEINS",
    "LIU2026_CELL_DECELERATED_AGING_PROTEINS",
    "SAUL_SEN_MAYO",
    "HAGR_CellAge_induces_senescence",
    "HAGR_CellAge_inhibits_senescence",
    "HAGR_GenAge_human",
]

AGE_CAT_COL = {
    "Accelerated aging proteome": "#B4493F",
    "Decelerated aging proteome": "#2B6CB0",
    "SASP/SenMayo": "#8064A2",
    "Cellular senescence": "#C9823B",
    "Aging curated": "#4E7E59",
}

CLASS_COL = {
    "exc": "#4A6FA5",
    "inh": "#8A65A8",
    "glia": "#4E9A72",
    "nonneuron": "#A65E5E",
    "vascular": "#C9823B",
}

DISEASE_SHORT = {
    "Alzheimer's Disease": "AD",
    "Parkinson's Disease": "PD",
    "Amyotrophic Lateral Sclerosis": "ALS",
    "Huntington's Disease": "HD",
    "Multiple Sclerosis": "MS",
    "Autism Spectrum Disorder": "ASD",
    "Schizophrenia": "SCZ",
    "Bipolar Disorder": "BD",
    "Major Depressive Disorder": "MDD",
    "Intellectual Disability": "ID",
    "Epilepsy": "Epilepsy",
    "Stroke": "Stroke",
    "Glioma": "Glioma",
    "Glioblastoma": "GBM",
    "Neuroblastoma": "NBL",
    "Down Syndrome": "Down",
}


def panel_tag(ax: plt.Axes, tag: str, dx: float = -0.05, dy: float = 1.02) -> None:
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=9.5, fontweight="bold", ha="right", va="bottom")


def lineage(sc: str) -> str:
    if sc in ("MICRO",):
        return "Microglia"
    if sc in ("AST",):
        return "Astrocyte"
    if sc in ("OLIGO", "OPC"):
        return "Oligo/OPC"
    if sc in ("ENDO", "VLMC"):
        return "Vascular"
    if sc in ("PVALB", "SST", "VIP", "LAMP5", "NDNF", "PAX6", "CHANDELIER"):
        return "Inhibitory"
    return "Excitatory"


def compact_name(text: str, max_chars: int = 34) -> str:
    text = str(text)
    repl = {
        "morphogenesis": "morphog.",
        "signaling": "sig.",
        "regulation": "reg.",
        "transport": "transp.",
        "cytoskeleton": "cytoskel.",
        "Neurofilament": "Nfl.",
        "Microglial": "Microgl.",
        "Postsynaptic": "Postsyn.",
        "Synaptic-vesicle": "Syn.-ves.",
        "chemokine": "chemok.",
        "complement": "compl.",
    }
    for a, b in repl.items():
        text = text.replace(a, b)
    if len(text) > max_chars:
        text = text[: max_chars - 1] + "."
    return text


def load_program_maps():
    renum = pd.read_csv(RENUM, sep="\t")
    old_to_new = {}
    for _, row in renum.iterrows():
        old = int(row["old_P"])
        new_str = str(row["new_P"])
        old_to_new[old] = None if new_str == "EXCLUDED" else int(new_str)

    names = pd.read_csv(NAMES, sep="\t").set_index("cnmf_component")
    star = {idx: ("*" if row["confidence"] == "brain-weak" else "") for idx, row in names.iterrows()}
    name_short = names["name_short"].to_dict()

    def label(old_p: int, max_chars: int = 34) -> str:
        new_n = old_to_new.get(int(old_p))
        if new_n is None:
            return f"[excluded old_P{old_p}]"
        return f"P{new_n}{star.get(int(old_p), '')} {compact_name(name_short.get(int(old_p), ''), max_chars=max_chars)}"

    return old_to_new, names, label


def load_disease_data(label_func):
    meta = json.load(open(PD_DIR / "meta.json"))
    cat_of = meta["category"]
    dn = meta["disease_n_genes"]
    long_df = pd.read_csv(PD_DIR / f"enrichment_long_N{PRIMARY_N}.tsv", sep="\t")
    long_df = long_df[~long_df["program"].isin(COHORT_TECHNICAL)].copy()
    long_df["fdr"] = multipletests(long_df["pval"].values, method="fdr_bh")[1]
    fdr_mat = long_df.pivot(index="program", columns="disease", values="fdr")
    or_mat = long_df.pivot(index="program", columns="disease", values="odds_ratio")
    neglog = -np.log10(fdr_mat.clip(lower=1e-300))
    sig_progs = sorted(set(long_df.loc[long_df["fdr"] < FDR_SIG, "program"]))
    sig_dis = [d for d in fdr_mat.columns if (fdr_mat[d] < FDR_SIG).any()]
    sig_dis = sorted(sig_dis, key=lambda d: (CAT_ORDER.index(cat_of[d]), -dn[d]))
    sub_nl = neglog.loc[sig_progs, sig_dis]
    sub_or = or_mat.loc[sig_progs, sig_dis]
    if len(sig_progs) > 2:
        z = linkage(pdist(sub_nl.fillna(0).values, metric="euclidean"), method="average")
        order = [sig_progs[i] for i in leaves_list(z)]
    else:
        order = sig_progs
    sub_nl = sub_nl.loc[order]
    sub_or = sub_or.loc[order]
    prog_labels = [label_func(p, max_chars=24) for p in order]

    dom = pd.read_csv(PD_DIR / "program_dom_subclass.tsv", sep="\t").set_index("program")
    best = (
        long_df[long_df["fdr"] < FDR_SIG]
        .sort_values("fdr")
        .groupby("program")
        .first()
        .reset_index()
    )
    best["neglogfdr"] = -np.log10(best["fdr"].clip(lower=1e-300))
    best = best.sort_values("neglogfdr", ascending=False).head(12)
    best["dom"] = best["program"].map(dom["dom_subclass"])
    best["lineage"] = best["dom"].map(lineage)

    rz = pd.read_csv(RZ, sep="\t", index_col=0)
    rz.columns = [int(c) for c in rz.columns]
    top_for_region = [p for p in best["program"].head(10) if p in rz.columns]
    region_order = ["FPPFC", "VLPFC", "DLPFC", "ACC", "M1", "S1", "S1E", "PoCG", "SMG", "AG", "SPL", "STG", "ITG", "V1"]
    rz_sub = rz[top_for_region].T
    rz_sub = rz_sub[[r for r in region_order if r in rz_sub.columns]]
    return long_df, sub_nl, sub_or, sig_dis, cat_of, prog_labels, best, dom, rz_sub


def load_aging_data():
    long_df = pd.read_csv(AGE_DIR / f"enrichment_long_N{PRIMARY_N}.tsv", sep="\t")
    sig = pd.read_csv(AGE_DIR / f"significant_pairs_N{PRIMARY_N}.tsv", sep="\t")
    sig_programs = list(sig.sort_values("fdr_retained54")["new_P"].drop_duplicates())
    rows = (
        long_df[long_df["new_P"].isin(sig_programs) & long_df["retained54"]]
        .drop_duplicates("new_P")
        .set_index("new_P")
        .loc[sig_programs]
        .reset_index()
    )
    row_labels = [f"{r.new_P} {compact_name(r.name_short, 24)}" for r in rows.itertuples()]
    row_index = dict(zip(rows["new_P"], range(len(rows))))
    set_order = [s for s in AGE_SET_ORDER if s in long_df["aging_set"].unique()]
    plot_df = long_df[
        long_df["retained54"] & long_df["new_P"].isin(row_index) & long_df["aging_set"].isin(set_order)
    ].copy()
    plot_df["neglog"] = -np.log10(plot_df["fdr_retained54"].clip(lower=1e-300))
    return long_df, sig, plot_df, set_order, row_labels, row_index


def draw_disease_dotplot(ax, sub_nl, sub_or, sig_dis, cat_of, prog_labels):
    cmap = LinearSegmentedColormap.from_list("disease_fdr", ["#F5F5F2", "#F2D2BC", "#D97854", "#A9363E", "#5D183A"])
    vals = sub_nl.values
    vmax = max(4.0, float(np.nanpercentile(vals[vals > 0], 98)) if (vals > 0).any() else 4.0)
    norm = Normalize(0, vmax)
    or_vals = sub_or.values
    smax = float(np.nanpercentile(or_vals[np.isfinite(or_vals)], 97))
    smax = max(smax, 4.0)

    xs, ys, cs, ss = [], [], [], []
    for i in range(sub_nl.shape[0]):
        for j in range(sub_nl.shape[1]):
            v = float(sub_nl.values[i, j])
            o = float(sub_or.values[i, j])
            xs.append(j)
            ys.append(i)
            cs.append(v)
            ss.append(4 + 34 * (min(o, smax) - 1) / max(smax - 1, 1e-9) if v > -np.log10(FDR_SIG) else 0)
    ax.set_facecolor("#FBFBFA")
    for j in range(sub_nl.shape[1] + 1):
        ax.axvline(j - 0.5, color="#E7E7E4", lw=0.28, zorder=0)
    for i in range(sub_nl.shape[0] + 1):
        ax.axhline(i - 0.5, color="#E7E7E4", lw=0.28, zorder=0)
    sc = ax.scatter(xs, ys, c=cs, s=ss, cmap=cmap, norm=norm, edgecolors="#333333", linewidths=0.15)
    for j, d in enumerate(sig_dis):
        ax.add_patch(plt.Rectangle((j - 0.5, -1.05), 1, 0.55, color=CAT_COL[cat_of[d]], clip_on=False, lw=0))
    ax.set_xlim(-0.5, sub_nl.shape[1] - 0.5)
    ax.set_ylim(sub_nl.shape[0] - 0.5, -1.15)
    ax.set_xticks(range(len(sig_dis)))
    ax.set_xticklabels([DISEASE_SHORT.get(d, d[:8]) for d in sig_dis], rotation=45, ha="right", fontsize=FS_MIN)
    ax.set_yticks(range(len(prog_labels)))
    ax.set_yticklabels(prog_labels, fontsize=4.15)
    ax.tick_params(length=1)
    ax.set_title("Brain-disease enrichment of cortical programs", loc="left", pad=4)
    cax = ax.inset_axes([0.972, 0.57, 0.018, 0.30])
    cb = plt.colorbar(sc, cax=cax)
    cb.set_label("-log10 FDR", fontsize=4.8, labelpad=1)
    cb.ax.tick_params(labelsize=4.8, length=1)
    handles = [ax.scatter([], [], s=4 + 34 * (o - 1) / max(smax - 1, 1e-9), color="#888888", edgecolors="#333333", linewidths=0.15, label=str(o)) for o in [2, 4, 8]]
    ax.legend(handles=handles, title="OR", loc="lower right", bbox_to_anchor=(0.995, 0.02), frameon=False, fontsize=4.8, title_fontsize=4.8, handletextpad=0.25)
    panel_tag(ax, "a", dx=-0.08)


def draw_disease_lollipop(ax, best, label_func):
    b = best.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(b))
    cols = [LIN_COL[l] for l in b["lineage"]]
    ax.hlines(y, 0, b["neglogfdr"], color=cols, lw=1.3)
    ax.scatter(b["neglogfdr"], y, s=20, c=cols, edgecolors="#222222", linewidths=0.22)
    labels = [label_func(p, max_chars=24) for p in b["program"]]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=4.8)
    for i, row in b.iterrows():
        ax.text(row["neglogfdr"] + 0.25, i, DISEASE_SHORT.get(row["disease"], row["disease"][:8]), va="center", fontsize=4.8, color="#444444")
    ax.set_xlim(0, b["neglogfdr"].max() * 1.42)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("Top disease-linked programs", loc="left", pad=3)
    ax.axvline(-np.log10(FDR_SIG), color="#888888", lw=0.45, ls="--")
    ax.spines[["top", "right"]].set_visible(False)
    panel_tag(ax, "b", dx=-0.09)


def draw_region_heatmap(ax, rz_sub, label_func):
    m = rz_sub.values
    vlim = max(1.0, float(np.nanpercentile(np.abs(m), 98)))
    cmap = LinearSegmentedColormap.from_list("rz", ["#2B6CB0", "#A5C9DE", "#F7F7F7", "#E8A582", "#B4493F"])
    im = ax.imshow(m, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim)
    ax.set_xticks(range(rz_sub.shape[1]))
    ax.set_xticklabels(rz_sub.columns, rotation=45, ha="right", fontsize=FS_MIN)
    ax.set_yticks(range(rz_sub.shape[0]))
    ax.set_yticklabels([label_func(p, max_chars=23) for p in rz_sub.index], fontsize=4.8)
    ax.tick_params(length=1)
    ax.set_title("Cortical-area bias of disease-linked programs", loc="left", pad=3)
    cb = plt.colorbar(im, cax=ax.inset_axes([1.02, 0.18, 0.028, 0.60]))
    cb.set_label("region z", fontsize=4.8, labelpad=1)
    cb.ax.tick_params(labelsize=4.8, length=1)
    panel_tag(ax, "c", dx=-0.09)


def draw_disease_category_counts(ax, long_df, cat_of):
    sig = long_df[long_df["fdr"] < FDR_SIG].copy()
    sig["category"] = sig["disease"].map(cat_of)
    counts = sig.groupby("category")["program"].nunique().reindex(CAT_ORDER).fillna(0).astype(int)
    cats = [c for c in CAT_ORDER if counts[c] > 0]
    y = np.arange(len(cats))[::-1]
    vals = counts.loc[cats].values
    ax.barh(y, vals, color=[CAT_COL[c] for c in cats], height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + 0.5, yi, str(v), va="center", ha="left", fontsize=5.0)
    short = {
        "Neurodegenerative": "ND",
        "Psychiatric/Neurodev": "Psych/ND",
        "Tumor": "Tumor",
        "Developmental": "Dev",
        "Other/Vascular": "Other",
    }
    ax.set_yticks(y)
    ax.set_yticklabels([short[c] for c in cats], fontsize=4.7)
    ax.set_xlim(0, max(vals) * 1.22)
    ax.set_xlabel("programs")
    ax.set_title("Disease classes", loc="left", pad=3)
    ax.grid(axis="x", color="#E6E6E6", lw=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    panel_tag(ax, "d", dx=-0.20)


def draw_disease_gene_cards(ax, old_to_new, dom):
    exemplars = [
        (40, "Alzheimer's Disease"),
        (1, "Autism Spectrum Disorder"),
        (1, "Epilepsy"),
        (59, "Glioma"),
    ]
    load = pd.read_csv(LOAD, sep="\t", index_col=0)
    load.index = [int(i) for i in load.index]
    bb = pd.read_csv(BB, sep="\t")
    bb.columns = [c.strip() for c in bb.columns]
    gcol = "Gene symbol" if "Gene symbol" in bb.columns else "Gene"
    mrna = bb[bb["Type"] == "mRNA"]
    dgenes = {d: set(mrna.loc[mrna["Disease"] == d, gcol].astype(str)) for _, d in exemplars}

    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0, 1.02, "Leading disease-overlap genes", transform=ax.transAxes, fontsize=7, fontweight="bold", ha="left")
    panel_tag(ax, "e", dx=-0.04, dy=1.02)
    for col, (old_p, dis) in enumerate(exemplars):
        row = load.loc[old_p]
        top = row.sort_values(ascending=False).head(PRIMARY_N)
        genes = [g for g in top.index if g in dgenes[dis]]
        genes = sorted(genes, key=lambda g: row[g], reverse=True)[:9]
        lineage_name = lineage(dom.loc[old_p, "dom_subclass"])
        color = LIN_COL[lineage_name]
        cx_i = col % 2
        cy_i = col // 2
        x0 = cx_i + 0.04
        y0 = 0.62 if cy_i == 0 else 0.03
        w = 0.88
        new_n = old_to_new.get(old_p, old_p)
        title = f"P{new_n} / {DISEASE_SHORT.get(dis, dis)}"
        rect = FancyBboxPatch((x0, y0 + 0.24), w, 0.15, boxstyle="round,pad=0.014,rounding_size=0.014", fc=color, ec="none", alpha=0.92)
        ax.add_patch(rect)
        ax.text(x0 + 0.03, y0 + 0.315, title, color="white", fontsize=5.7, fontweight="bold", va="center")
        body = "\n".join(genes[:4])
        ax.text(x0 + 0.03, y0 + 0.19, body, fontsize=4.8, ha="left", va="top", color="#222222", linespacing=1.03)


def draw_aging_dotplot(ax, plot_df, set_order, row_labels, row_index):
    col_index = {s: i for i, s in enumerate(set_order)}
    cmap = LinearSegmentedColormap.from_list("aging_fdr", ["#F5F5F2", "#F1D2B9", "#D97854", "#A9363E", "#5D183A"])
    vmax = max(4.0, float(np.nanpercentile(plot_df["neglog"], 98)))
    norm = Normalize(0, vmax)
    odds_hi = float(np.nanpercentile(plot_df["odds_ratio"].replace(np.inf, np.nan).dropna(), 97))
    odds_hi = max(odds_hi, 4.0)
    xs, ys, cs, ss = [], [], [], []
    for _, row in plot_df.iterrows():
        xs.append(col_index[row["aging_set"]])
        ys.append(row_index[row["new_P"]])
        cs.append(row["neglog"])
        if row["fdr_retained54"] < FDR_SIG:
            odds = min(float(row["odds_ratio"]), odds_hi)
            ss.append(7 + 42 * (odds - 1) / max(odds_hi - 1, 1e-9))
        else:
            ss.append(0)
    ax.set_facecolor("#FBFBFA")
    for j in range(len(set_order) + 1):
        ax.axvline(j - 0.5, color="#E7E7E4", lw=0.28, zorder=0)
    for i in range(len(row_labels) + 1):
        ax.axhline(i - 0.5, color="#E7E7E4", lw=0.28, zorder=0)
    sc = ax.scatter(xs, ys, c=cs, s=ss, cmap=cmap, norm=norm, edgecolors="#333333", linewidths=0.15)
    ax.set_xlim(-0.5, len(set_order) - 0.5)
    ax.set_ylim(len(row_labels) - 0.5, -0.5)
    ax.set_xticks(range(len(set_order)))
    ax.set_xticklabels([AGE_SET_LABEL.get(s, s) for s in set_order], rotation=35, ha="right", fontsize=FS_MIN)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=4.8)
    ax.tick_params(length=1)
    ax.set_title("Aging gene-set recovery", loc="left", pad=3)
    cb = plt.colorbar(sc, cax=ax.inset_axes([0.970, 0.57, 0.018, 0.30]))
    cb.set_label("-log10 FDR", fontsize=4.8, labelpad=1)
    cb.ax.tick_params(labelsize=4.8, length=1)
    handles = [ax.scatter([], [], s=7 + 42 * (o - 1) / max(odds_hi - 1, 1e-9), color="#888888", edgecolors="#333333", linewidths=0.15, label=str(o)) for o in [2, 4, 8]]
    ax.legend(handles=handles, title="OR", loc="lower right", bbox_to_anchor=(0.995, 0.02), frameon=False, fontsize=4.8, title_fontsize=4.8, handletextpad=0.25)
    panel_tag(ax, "f", dx=-0.08)


def draw_aging_top(ax, sig):
    top = sig.sort_values("fdr_retained54").head(10).copy().iloc[::-1].reset_index(drop=True)
    top["neglog"] = -np.log10(top["fdr_retained54"].clip(lower=1e-300))
    y = np.arange(len(top))
    colors = [AGE_CAT_COL.get(c, "#777777") for c in top["category"]]
    labels = [f"{r.new_P} {compact_name(r.name_short, 16)} / {AGE_SET_LABEL.get(r.aging_set, r.aging_set).replace(' proteins', '')}" for r in top.itertuples()]
    ax.hlines(y, 0, top["neglog"], color=colors, lw=1.3)
    ax.scatter(top["neglog"], y, s=20, c=colors, edgecolors="#222222", linewidths=0.22)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=4.35)
    ax.set_xlabel("-log10 FDR")
    ax.set_title("Top aging-linked programs", loc="left", pad=3)
    ax.grid(axis="x", color="#E5E5E5", lw=0.30)
    ax.spines[["top", "right"]].set_visible(False)
    panel_tag(ax, "g", dx=-0.09)


def draw_aging_bipartite(ax, sig):
    top = sig.sort_values("fdr_retained54").head(12).copy()
    set_order = [s for s in AGE_SET_ORDER if s in set(top["aging_set"])]
    prog_order = (
        top.sort_values("fdr_retained54")
        .drop_duplicates("new_P")
        .head(9)[["new_P", "name_short", "dominant_class"]]
        .copy()
    )
    prog_ids = prog_order["new_P"].tolist()
    top = top[top["new_P"].isin(prog_ids)].copy()

    left_y = {s: y for s, y in zip(set_order, np.linspace(0.84, 0.16, len(set_order)))}
    right_y = {p: y for p, y in zip(prog_ids, np.linspace(0.88, 0.12, len(prog_ids)))}
    lx, rx = 0.01, 0.73
    ax.set_xlim(-0.02, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_tag(ax, "h", dx=-0.035, dy=1.01)

    # edges first
    max_or = max(4.0, float(np.nanpercentile(top["odds_ratio"], 95)))
    for _, row in top.sort_values("fdr_retained54", ascending=False).iterrows():
        y0 = left_y[row["aging_set"]]
        y1 = right_y[row["new_P"]]
        color = AGE_CAT_COL.get(row["category"], "#777777")
        lw = 0.45 + 2.0 * (min(float(row["odds_ratio"]), max_or) - 1) / max(max_or - 1, 1e-9)
        path = MplPath(
            [(lx + 0.20, y0), (0.36, y0), (0.50, y1), (rx - 0.04, y1)],
            [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
        )
        ax.add_patch(PathPatch(path, fc="none", ec=color, lw=lw, alpha=0.42, capstyle="round"))

    # left gene-set nodes
    for s in set_order:
        y = left_y[s]
        cat = top.loc[top["aging_set"] == s, "category"].iloc[0]
        color = AGE_CAT_COL.get(cat, "#777777")
        ax.add_patch(FancyBboxPatch((lx, y - 0.035), 0.20, 0.070, boxstyle="round,pad=0.006,rounding_size=0.010", fc=color, ec="white", lw=0.4))
        ax.text(lx + 0.10, y, AGE_SET_LABEL.get(s, s), ha="center", va="center", fontsize=4.3, color="white", fontweight="bold")

    # right program nodes
    prog_cls = {r.new_P: r.dominant_class for r in prog_order.itertuples()}
    for p in prog_ids:
        y = right_y[p]
        color = CLASS_COL.get(prog_cls.get(p, ""), "#777777")
        ax.scatter([rx], [y], s=38, c=[color], edgecolors="#222222", linewidths=0.25, zorder=4)
        ax.text(rx + 0.025, y, f"{p}", ha="left", va="center", fontsize=4.6)

    ax.text(lx, 0.97, "aging set", fontsize=5.0, color="#555555", ha="left")
    ax.text(rx - 0.02, 0.97, "program", fontsize=5.0, color="#555555", ha="left")
    ax.text(0.03, 0.02, "line width = OR", transform=ax.transAxes, fontsize=4.4, color="#666666", ha="left", va="bottom")


def draw_aging_donuts(ax, sig):
    keep_cats = ["Accelerated aging proteome", "Decelerated aging proteome", "SASP/SenMayo", "Cellular senescence", "Aging curated"]
    count = (
        sig.assign(category=pd.Categorical(sig["category"], keep_cats, ordered=True))
        .groupby(["category", "dominant_class"], observed=True)["new_P"]
        .nunique()
        .reset_index(name="n_programs")
    )
    classes = ["nonneuron", "vascular", "glia", "exc", "inh"]
    piv = count.pivot(index="category", columns="dominant_class", values="n_programs").reindex(keep_cats).fillna(0)
    piv = piv[[c for c in classes if c in piv.columns]]
    x_pos = np.linspace(0.0, 3.35, len(piv))
    donut_r = 0.41
    ax.set_xlim(x_pos[0] - 0.56, x_pos[-1] + 0.56)
    ax.set_ylim(-0.96, 0.72)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.set_title("Cellular identity", loc="left", pad=2)
    panel_tag(ax, "i", dx=-0.07)
    for x, (cat, vals) in zip(x_pos, piv.iterrows()):
        total = float(vals.sum())
        start = 90
        for cls, value in vals.items():
            if value <= 0:
                continue
            theta = 360 * float(value) / total
            ax.add_patch(Wedge((x, 0), donut_r, start, start + theta, width=0.18, fc=CLASS_COL.get(cls, "#777777"), ec="white", lw=0.45))
            start += theta
        ax.text(x, 0, f"{int(total)}", fontsize=7.0, fontweight="bold", ha="center", va="center")
        ax.text(x, -0.66, str(cat).replace(" proteome", "").replace("Cellular ", "Cell "), fontsize=4.9, ha="center", va="top", rotation=22)
    handles = [Patch(fc=CLASS_COL.get(c, "#777777"), ec="none", label=c) for c in piv.columns]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.00, 1.03), ncol=5, frameon=False, fontsize=3.8, title=None, handlelength=0.50, columnspacing=0.22)


def draw_aging_gene_cards(ax, sig):
    chosen = [
        ("P54", "LIU2026_CELL_ACCELERATED_AGING_PROTEINS"),
        ("P49", "LIU2026_CELL_ACCELERATED_AGING_PROTEINS"),
        ("P8", "LIU2026_CELL_DECELERATED_AGING_PROTEINS"),
        ("P53", "LIU2026_CELL_ACCELERATED_AGING_PROTEINS"),
    ]
    records = []
    for p, term in chosen:
        row = sig[(sig["new_P"] == p) & (sig["aging_set"] == term)]
        if not row.empty:
            records.append(row.iloc[0])
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.axis("off")
    ax.text(0, 1.02, "Aging-overlap genes", transform=ax.transAxes, fontsize=7, fontweight="bold", ha="left")
    panel_tag(ax, "j", dx=-0.07, dy=1.02)
    for col, row in enumerate(records):
        genes = str(row["overlap_genes"]).split(";")[:9]
        x0 = (col % 2) + 0.05
        y0 = 1.04 if col < 2 else 0.08
        w = 0.88
        color = AGE_CAT_COL.get(row["category"], "#777777")
        title = f"{row['new_P']} / {AGE_SET_LABEL.get(row['aging_set'], row['aging_set'])}"
        rect = FancyBboxPatch((x0, y0 + 0.50), w, 0.24, boxstyle="round,pad=0.014,rounding_size=0.014", fc=color, ec="none", alpha=0.92)
        ax.add_patch(rect)
        ax.text(x0 + 0.03, y0 + 0.62, title, color="white", fontsize=4.9, fontweight="bold", va="center")
        ax.text(x0 + 0.03, y0 + 0.43, f"OR {row['odds_ratio']:.1f}; FDR {row['fdr_retained54']:.1e}", fontsize=3.7, color="#333333", ha="left", va="top")
        body = "\n".join([", ".join(genes[:2]), genes[2] if len(genes) > 2 else ""])
        ax.text(x0 + 0.03, y0 + 0.26, body, fontsize=4.0, ha="left", va="top", color="#222222", linespacing=1.0)


def main() -> None:
    old_to_new, names, label_func = load_program_maps()
    dis_long, sub_nl, sub_or, sig_dis, cat_of, prog_labels, best, dom, rz_sub = load_disease_data(label_func)
    _, age_sig, age_plot, age_set_order, age_row_labels, age_row_index = load_aging_data()

    fig = plt.figure(figsize=(180 / 25.4, 222 / 25.4), dpi=300)
    fig.patch.set_facecolor("white")

    axes = {
        "a": fig.add_axes([0.145, 0.620, 0.330, 0.350]),
        "b": fig.add_axes([0.655, 0.720, 0.280, 0.240]),
        "c": fig.add_axes([0.145, 0.500, 0.330, 0.075]),
        "d": fig.add_axes([0.555, 0.570, 0.110, 0.080]),
        "e": fig.add_axes([0.690, 0.510, 0.275, 0.140]),
        "f": fig.add_axes([0.145, 0.265, 0.330, 0.190]),
        "g": fig.add_axes([0.655, 0.285, 0.280, 0.170]),
        "h": fig.add_axes([0.035, 0.035, 0.340, 0.140]),
        "i": fig.add_axes([0.355, 0.020, 0.385, 0.170]),
        "j": fig.add_axes([0.790, 0.035, 0.180, 0.140]),
    }

    draw_disease_dotplot(axes["a"], sub_nl, sub_or, sig_dis, cat_of, prog_labels)
    draw_disease_lollipop(axes["b"], best, label_func)
    draw_region_heatmap(axes["c"], rz_sub, label_func)
    draw_disease_category_counts(axes["d"], dis_long, cat_of)
    draw_disease_gene_cards(axes["e"], old_to_new, dom)
    draw_aging_dotplot(axes["f"], age_plot, age_set_order, age_row_labels, age_row_index)
    draw_aging_top(axes["g"], age_sig)
    draw_aging_bipartite(axes["h"], age_sig)
    draw_aging_donuts(axes["i"], age_sig)
    draw_aging_gene_cards(axes["j"], age_sig)

    for ext in ["pdf", "png", "svg"]:
        path = OUT / f"Fig8_imagegen_layout_v17.{ext}"
        if ext == "png":
            fig.savefig(path, dpi=450)
        else:
            fig.savefig(path)
    print("SAVED", OUT / "Fig8_imagegen_layout_v17.pdf")
    print("page_mm", 180, 222)
    print("disease_matrix", sub_nl.shape, "aging_pairs", len(age_sig))


if __name__ == "__main__":
    main()
