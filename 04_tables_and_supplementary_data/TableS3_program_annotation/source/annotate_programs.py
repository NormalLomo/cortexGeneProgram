#!/usr/bin/env python
"""Annotate 60 cNMF gene programs with meaningful short names.
Integrates: top loading genes, program x subclass identity, spatial coloc,
variability, biotype flags. Writes program_annotation.tsv. Does NOT modify sources.
"""
import os
import numpy as np
import pandas as pd

BASE = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
LOAD = "CORTEX_PROGRAM_ROOT/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv"
RSCM = os.path.join(BASE, "region_subclass_program_mean.tsv")
COLOC = os.path.join(BASE, "program_celltype_coloc.tsv")
VAR = os.path.join(BASE, "program_variability.tsv")
GENEINFO = "CORTEX_PROGRAM_ROOT/inputs/geneInfo_snRNA.csv"
OUT = os.path.join(BASE, "program_annotation.tsv")

# ---- 1. top loading genes per program ----
load = pd.read_csv(LOAD, sep="\t", index_col=0)
# rows = programs (1..60), columns = genes
load.index = [str(i) for i in load.index]
TOPN = 30
top_genes = {}
for prog in load.index:
    row = load.loc[prog]
    tg = row.sort_values(ascending=False).head(TOPN).index.tolist()
    top_genes[prog] = tg

# ---- 2. program x subclass mean (average over regions) ----
rscm = pd.read_csv(RSCM, sep="\t")
rscm["program"] = rscm["program"].astype(str)
# mean program score per subclass, averaging equally across regions
psc = rscm.groupby(["program", "subclass"])["mean"].mean().unstack("subclass")
# specificity from positive part (program scores can be negative)
def specificity(vec):
    v = vec.values.astype(float)
    v = np.clip(v, 0, None)
    s = v.sum()
    if s <= 0:
        return 0.0, "(none)", 0.0
    p = v / s
    p_nz = p[p > 0]
    ent = -(p_nz * np.log(p_nz)).sum()
    ent_norm = ent / np.log(len(v))  # 0=specific, 1=uniform
    spec = 1.0 - ent_norm            # 1=specific, 0=broad
    dom = vec.idxmax()
    return spec, dom, vec.max()

spec_map, dom_map = {}, {}
for prog in psc.index:
    spec, dom, _ = specificity(psc.loc[prog])
    spec_map[prog] = round(spec, 4)
    dom_map[prog] = dom

# ---- 3. spatial coloc top cell type ----
coloc = pd.read_csv(COLOC, sep="\t", index_col=0)
coloc.index = [i.replace("program_", "") for i in coloc.index]
spatial_top = {p: coloc.loc[p].idxmax() for p in coloc.index}

# ---- 4. variability ----
var = pd.read_csv(VAR, sep="\t")
var["program"] = var["program"].astype(str)
var_map = dict(zip(var["program"], var["class"]))
eta2_map = dict(zip(var["program"], var["eta2_region"]))

# ---- 5. biotype ----
gi = pd.read_csv(GENEINFO)
biotype = dict(zip(gi["gene_name"], gi["gene_biotype"]))

# ---- naming knowledge ----
# subclass -> (short label, class)
SUB = {
    "AST": ("Astro", "glia"),
    "CHANDELIER": ("Chandelier", "inh"),
    "ENDO": ("Endothelial", "vascular"),
    "ET": ("L5-ET", "exc"),
    "L2-L3 IT LINC00507": ("L2/3-IT", "exc"),
    "L3-L4 IT RORB": ("L3/4-IT", "exc"),
    "L4-L5 IT RORB": ("L4/5-IT", "exc"),
    "L6 CAR3": ("L6-CAR3", "exc"),
    "L6 CT": ("L6-CT", "exc"),
    "L6 IT": ("L6-IT", "exc"),
    "L6B": ("L6b", "exc"),
    "LAMP5": ("LAMP5-IN", "inh"),
    "MICRO": ("Microglia", "nonneuron"),
    "NDNF": ("NDNF-IN", "inh"),
    "NP": ("L5/6-NP", "exc"),
    "OLIGO": ("Oligo", "glia"),
    "OPC": ("OPC", "glia"),
    "PAX6": ("PAX6-IN", "inh"),
    "PVALB": ("PVALB-IN", "inh"),
    "SST": ("SST-IN", "inh"),
    "VIP": ("VIP-IN", "inh"),
    "VLMC": ("VLMC", "vascular"),
}

# functional marker sets -> (tag, priority). Checked against top genes.
FUNC = [
    ("Myelin",       {"MOBP","PLP1","MBP","MAG","MOG","CLDN11","CNP","MYRF","OPALIN","ASPA"}),
    ("Astro",        {"AQP4","GFAP","SLC1A2","SLC1A3","GJA1","ALDH1L1","FGFR3"}),
    ("Microglia",    {"CX3CR1","P2RY12","CSF1R","C1QA","C1QB","C1QC","AIF1","CTSS","TYROBP","ITGAM"}),
    ("Vascular-BBB", {"CLDN5","FLT1","PECAM1","VWF","CLEC14A","A2M","ABCB1","EPAS1"}),
    ("Mural-VLMC",   {"DCN","COL1A1","COL1A2","PDGFRB","RGS5","ACTA2","TAGLN","COL3A1","LUM","PTGDS"}),
    ("OPC",          {"PDGFRA","CSPG4","OLIG1","OLIG2","SOX10","LHFPL3","PCDH15"}),
    ("Mito-ubiq",    {"MT-CO1","MT-CO2","MT-CO3","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ATP6","MT-CYB","MTRNR2L12"}),
    ("Ribosomal",    {"RPL13","RPS18","RPL10","RPS27","RPL41","RPS6","RPL3","RPLP1","RPS2","RPS3"}),
    ("Synaptic",     {"SYT1","SNAP25","GRIN1","GRIN2B","DLG2","NRXN1","NRXN3","NLGN1","CAMK2A","GRIA1","SYN1"}),
    ("HeatShock",    {"HSPA1A","HSPA1B","DNAJB1","HSPH1","FOS","JUN","EGR1","ARC","NR4A1"}),
]

def func_tag(genes):
    gset = set(genes)
    best = None; best_hits = 0
    for tag, mk in FUNC:
        hits = len(gset & mk)
        if hits > best_hits:
            best_hits = hits; best = tag
    return best, best_hits

def sig_gene(genes, dom_label):
    # pick a representative protein-coding gene from top genes, avoid generic mito/ribo
    skip_prefix = ("MT-","RPL","RPS","MTRNR")
    for g in genes:
        if g.startswith(skip_prefix):
            continue
        bt = biotype.get(g, "")
        return g
    return genes[0]

def lncrna_flag(genes):
    flags = []
    for g in genes[:10]:
        bt = biotype.get(g, "")
        if "lncRNA" in bt or "lincRNA" in bt or "antisense" in bt:
            flags.append(g)
    return flags

# ---- assemble ----
rows = []
SPEC_BROAD = 0.10  # below -> broadly shared (uniform across subclasses)
for prog in [str(i) for i in range(1, 61)]:
    tg = top_genes[prog]
    dom = dom_map[prog]
    spec = spec_map[prog]
    sub_label, klass = SUB[dom]
    ftag, fhits = func_tag(tg)
    lflags = lncrna_flag(tg)
    sg = sig_gene(tg, sub_label)

    # decide name
    # 1) strong functional signature (>=2 marker hits) overrides for glia/nonneuron/vascular
    name = None
    if ftag and fhits >= 2 and ftag in ("Myelin","Astro","Microglia","Vascular-BBB","Mural-VLMC","OPC"):
        name = f"{ftag} ({sg})"
    elif ftag in ("Mito-ubiq","Ribosomal","HeatShock") and fhits >= 2:
        name = f"{ftag} ({sg})"
    elif spec < SPEC_BROAD:
        # broadly shared -> functional / ubiquitous
        if ftag and fhits >= 1:
            name = f"{ftag} ({sg})"
        else:
            name = f"Ubiquitous ({sg})"
    else:
        # identity-driven name
        name = f"{sub_label} ({sg})"

    rows.append({
        "program": f"P{int(prog):02d}",
        "_pnum": int(prog),
        "dominant_subclass": dom,
        "subclass_specificity": spec,
        "class": klass,
        "top5_genes": ",".join(tg[:5]),
        "spatial_top_celltype": spatial_top.get(prog, "NA"),
        "variable_stable": var_map.get(prog, "NA"),
        "_base_name": name,
        "notable_lncRNA_or_flag": (",".join(lflags) if lflags else ""),
        "_eta2": round(float(eta2_map.get(prog, np.nan)), 4),
    })

df = pd.DataFrame(rows).sort_values("_pnum").reset_index(drop=True)

# ---- disambiguate duplicate names ----
# names currently like "P01 ..." not yet prefixed; build proposed_name = P## + base
def with_prefix(r):
    return f"{r['program']} {r['_base_name']}"
df["proposed_name"] = df.apply(with_prefix, axis=1)
# the P## prefix already makes every name unique, but the human-readable BASE
# (without prefix) may collide; disambiguate the base part using 2nd signature gene
base_counts = df["_base_name"].value_counts()
dup_bases = set(base_counts[base_counts > 1].index)
for i, r in df.iterrows():
    if r["_base_name"] in dup_bases:
        prog = str(r["_pnum"])
        tg = top_genes[prog]
        # find a 2nd distinguishing protein-coding gene different from one already used
        used = r["_base_name"].split("(")[-1].rstrip(")")
        skip_prefix = ("MT-","RPL","RPS","MTRNR")
        alt = None
        for g in tg:
            if g == used or g.startswith(skip_prefix):
                continue
            alt = g; break
        if alt:
            newbase = r["_base_name"].replace(f"({used})", f"({used}/{alt})")
            df.at[i, "_base_name"] = newbase
            df.at[i, "proposed_name"] = f"{r['program']} {newbase}"

# final columns in required order
out_cols = ["program","dominant_subclass","subclass_specificity","class",
            "top5_genes","spatial_top_celltype","variable_stable",
            "proposed_name","notable_lncRNA_or_flag"]
out = df[out_cols].copy()

# ---- validate ----
assert len(out) == 60, f"expected 60 rows got {len(out)}"
assert out["proposed_name"].notna().all(), "missing proposed_name"
assert out["proposed_name"].nunique() == 60, f"duplicate proposed_name: {out['proposed_name'].nunique()}"

out.to_csv(OUT, sep="\t", index=False)

# ---- report ----
print("WROTE", OUT, "rows", len(out))
print("=== 60 PROGRAM -> PROPOSED NAME ===")
print("; ".join(out["proposed_name"].tolist()))
print()
print("=== FULL TABLE ===")
for _, r in out.iterrows():
    print(f"{r['proposed_name']}\t[dom={r['dominant_subclass']} spec={r['subclass_specificity']} "
          f"spatial={r['spatial_top_celltype']} {r['variable_stable']}] top5={r['top5_genes']}"
          + (f" lnc={r['notable_lncRNA_or_flag']}" if r['notable_lncRNA_or_flag'] else ""))
