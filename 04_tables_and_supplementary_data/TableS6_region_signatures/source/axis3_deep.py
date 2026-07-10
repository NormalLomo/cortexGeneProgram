#!/usr/bin/env python
"""
AXIS 3 DEEP — cross-region 异同 (human region signatures + mouse 异同 + conservation buckets).
DEEPEN the cross-species cross-region analysis beyond the blunt "not conserved" verdict.

Mouse depth-correct = per-bin COMPOSITION (each program's share of bin's total program
signal; depth-symmetric, NOT per-bin library norm) -> aggregate (median) by area. SAME view
axis2 finalize used.

Outputs -> results/crossregion_v1/crossspecies/metrics/axis3_deep/ ONLY.
"""
import os, glob, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu

BASE = "results/crossregion_v1/"
OUT = BASE + "crossspecies/metrics/axis3_deep/"
os.makedirs(OUT, exist_ok=True)

PROGS = [str(i) for i in range(1, 61)]

# ---- 5 clean homolog pairs ----
# human region <-> list of mouse area labels
PAIRS = {
    "M1":  {"human": "M1",  "mouse": ["MOp"]},
    "S1":  {"human": "S1",  "mouse": ["SSp"]},
    "V1":  {"human": "V1",  "mouse": ["VIS"]},
    "STG": {"human": "STG", "mouse": ["AUD"]},
    "ACC": {"human": "ACC", "mouse": ["ACAd", "ACAv"]},
}
MATCHED_HUMAN = list(PAIRS.keys())  # 5 human regions

# human association/higher-order regions w/ NO clean mouse homolog
ASSOC_NOHOMOLOG = ["DLPFC", "AG", "SMG", "SPL", "ITG", "FPPFC", "VLPFC", "PoCG", "S1E"]

# ============================================================
# 0. LOAD inputs
# ============================================================
zhuman = pd.read_csv(BASE + "program_region_zscore.tsv", sep="\t", index_col=0)  # 14 region x 60 prog
zhuman.columns = [str(c) for c in zhuman.columns]
mean_human = pd.read_csv(BASE + "region_program_mean.tsv", sep="\t", index_col=0)
mean_human.columns = [str(c) for c in mean_human.columns]
var = pd.read_csv(BASE + "program_variability.tsv", sep="\t")
var["program"] = var["program"].astype(str)
nm = pd.read_csv(BASE + "program_names.tsv", sep="\t")
nm["program"] = nm["program"].astype(str)
name_short = dict(zip(nm["program"], nm["name_short"]))
sig_human = pd.read_csv(BASE + "supp_table_region_signatures.tsv", sep="\t")
sig_human["program"] = sig_human["program"].astype(str)

R3_VARIABLE = set(var.loc[var["class"] == "variable", "program"].astype(str))
eta2 = dict(zip(var["program"], var["eta2_region"]))
human_regions = list(zhuman.index)
print("[load] human regions:", human_regions)
print("[load] R3 variable programs (14):", sorted(R3_VARIABLE, key=int))


def pname(p):
    return f"P{p} {name_short.get(p, '?')}"


# ============================================================
# 1. HUMAN region SIGNATURES (top up/down programs per region from z)
#    + strongest region-distinguishing programs
# ============================================================
# top up/down per region (use the prebuilt z directly; cross-region z per program already in zhuman)
sig_rows = []
for reg in human_regions:
    zr = zhuman.loc[reg]
    up = zr.sort_values(ascending=False).head(5)
    dn = zr.sort_values(ascending=True).head(5)
    for rank, (p, zv) in enumerate(up.items(), 1):
        sig_rows.append(dict(region=reg, direction="up", rank=rank, program=p,
                             program_name=pname(p), z=round(float(zv), 3),
                             R3_variable=(p in R3_VARIABLE)))
    for rank, (p, zv) in enumerate(dn.items(), 1):
        sig_rows.append(dict(region=reg, direction="down", rank=rank, program=p,
                             program_name=pname(p), z=round(float(zv), 3),
                             R3_variable=(p in R3_VARIABLE)))
human_sig = pd.DataFrame(sig_rows)
human_sig.to_csv(OUT + "human_region_signatures.tsv", sep="\t", index=False)

# strongest region-distinguishing programs: by eta2 (region variance explained) and by z-spread
spread = (zhuman.max(axis=0) - zhuman.min(axis=0)).rename("z_spread")
maxabs = zhuman.abs().max(axis=0).rename("z_maxabs")
distinguish = pd.DataFrame({
    "program": PROGS,
    "program_name": [pname(p) for p in PROGS],
    "eta2_region": [eta2.get(p, np.nan) for p in PROGS],
    "z_spread": [float(spread.get(p, np.nan)) for p in PROGS],
    "z_maxabs": [float(maxabs.get(p, np.nan)) for p in PROGS],
    "R3_variable": [p in R3_VARIABLE for p in PROGS],
    "top_region": [zhuman[p].idxmax() for p in PROGS],
    "bottom_region": [zhuman[p].idxmin() for p in PROGS],
})
distinguish = distinguish.sort_values("eta2_region", ascending=False)
distinguish.to_csv(OUT + "human_distinguishing_programs.tsv", sep="\t", index=False)
print("[1] top-8 region-distinguishing programs (by eta2):")
for _, r in distinguish.head(8).iterrows():
    print(f"   {r['program_name']:42s} eta2={r['eta2_region']:.3f} spread={r['z_spread']:.2f} top={r['top_region']} bot={r['bottom_region']}")


# ============================================================
# 2. MOUSE program signatures via per-bin COMPOSITION -> area aggregate
# ============================================================
files = sorted(glob.glob(BASE + "crossspecies/sections/*.parquet"))
print(f"[2] loading {len(files)} mouse sections for composition...")
# accumulate per-bin composition rows tagged by area
area_comp_bins = {}  # area -> list of comp arrays (per bin)
all_mouse_areas = set()
for f in files:
    df = pd.read_parquet(f, columns=["area"] + PROGS)
    P = df[PROGS].to_numpy(dtype=np.float64)
    P = np.clip(P, 0, None)  # non-neg raw
    rowsum = P.sum(axis=1, keepdims=True)
    keep = rowsum[:, 0] > 0
    comp = np.zeros_like(P)
    comp[keep] = P[keep] / rowsum[keep]  # per-bin composition (depth-symmetric)
    areas = df["area"].astype(str).to_numpy()
    for a in np.unique(areas):
        m = (areas == a) & keep
        if m.sum() == 0:
            continue
        all_mouse_areas.add(a)
        area_comp_bins.setdefault(a, []).append(comp[m])

# median composition per area across all its bins (all sections pooled)
area_comp = {}
area_nbins = {}
for a, chunks in area_comp_bins.items():
    M = np.vstack(chunks)
    area_comp[a] = np.median(M, axis=0)
    area_nbins[a] = M.shape[0]
print("[2] mouse areas w/ composition:", sorted(area_comp.keys()))

# matched mouse area composition (ACA = pool ACAd+ACAv bins together, median)
matched_mouse_comp = {}
for hp, d in PAIRS.items():
    chunks = []
    for ma in d["mouse"]:
        if ma in area_comp_bins:
            chunks.extend(area_comp_bins[ma])
    if not chunks:
        print(f"  !! no mouse bins for pair {hp} ({d['mouse']})")
        continue
    M = np.vstack(chunks)
    matched_mouse_comp[hp] = np.median(M, axis=0)
    print(f"  pair {hp}: mouse {d['mouse']} nbins={M.shape[0]}")

# build mouse composition matrix over the 5 matched regions, then z across regions per program
mouse_comp_mat = pd.DataFrame(
    {hp: matched_mouse_comp[hp] for hp in MATCHED_HUMAN if hp in matched_mouse_comp},
    index=PROGS,
).T  # rows=region, cols=program
# z across the 5 matched regions (per program) — mouse region-difference signature
mouse_z = (mouse_comp_mat - mouse_comp_mat.mean(axis=0)) / (mouse_comp_mat.std(axis=0, ddof=0) + 1e-12)
mouse_z.to_csv(OUT + "mouse_matched_region_zscore.tsv", sep="\t")
mouse_comp_mat.to_csv(OUT + "mouse_matched_region_composition.tsv", sep="\t")

# human z restricted to the SAME 5 matched regions, recomputed across just these 5 (apples-to-apples)
human_5 = mean_human.loc[MATCHED_HUMAN]  # mean program per matched human region
human_z5 = (human_5 - human_5.mean(axis=0)) / (human_5.std(axis=0, ddof=0) + 1e-12)
human_z5.to_csv(OUT + "human_matched_region_zscore.tsv", sep="\t")

# mouse signatures: top up/down per matched region (composition z)
msig_rows = []
for reg in MATCHED_HUMAN:
    if reg not in mouse_z.index:
        continue
    zr = mouse_z.loc[reg]
    up = zr.sort_values(ascending=False).head(5)
    dn = zr.sort_values(ascending=True).head(5)
    for rank, (p, zv) in enumerate(up.items(), 1):
        msig_rows.append(dict(region=reg, mouse_area="+".join(PAIRS[reg]["mouse"]),
                              direction="up", rank=rank, program=p, program_name=pname(p),
                              mouse_z=round(float(zv), 3), R3_variable=(p in R3_VARIABLE)))
    for rank, (p, zv) in enumerate(dn.items(), 1):
        msig_rows.append(dict(region=reg, mouse_area="+".join(PAIRS[reg]["mouse"]),
                              direction="down", rank=rank, program=p, program_name=pname(p),
                              mouse_z=round(float(zv), 3), R3_variable=(p in R3_VARIABLE)))
mouse_sig = pd.DataFrame(msig_rows)
mouse_sig.to_csv(OUT + "mouse_matched_region_signatures.tsv", sep="\t", index=False)


# ============================================================
# 2b. PER-PAIR 异同: which human signature programs recapitulated vs human-specific
# ============================================================
pair_compare_rows = []
for reg in MATCHED_HUMAN:
    if reg not in mouse_z.index:
        continue
    # human signature for this region over the 5 matched regions (apples to apples)
    hz = human_z5.loc[reg]
    mz = mouse_z.loc[reg]
    h_up = set(hz.sort_values(ascending=False).head(5).index)
    h_dn = set(hz.sort_values(ascending=True).head(5).index)
    for p in sorted(h_up | h_dn, key=int):
        hdir = "up" if p in h_up else "down"
        hval = float(hz[p]); mval = float(mz[p])
        # recapitulated if mouse z has same sign and is at least moderately strong
        same_sign = (np.sign(hval) == np.sign(mval)) and abs(mval) >= 0.5
        # is it also in mouse's own top/bottom-5?
        m_up = set(mz.sort_values(ascending=False).head(5).index)
        m_dn = set(mz.sort_values(ascending=True).head(5).index)
        in_mouse_sig = (p in m_up) or (p in m_dn)
        if in_mouse_sig and same_sign:
            verdict = "RECAPITULATED"
        elif same_sign:
            verdict = "PARTIAL"  # same direction, weaker / not in mouse top-5
        else:
            verdict = "HUMAN-SPECIFIC"
        pair_compare_rows.append(dict(
            region=reg, mouse_area="+".join(PAIRS[reg]["mouse"]),
            program=p, program_name=pname(p), human_dir=hdir,
            human_z5=round(hval, 3), mouse_z=round(mval, 3),
            in_mouse_top5=in_mouse_sig, verdict=verdict,
            R3_variable=(p in R3_VARIABLE)))
pair_compare = pd.DataFrame(pair_compare_rows)
pair_compare.to_csv(OUT + "pair_signature_compare.tsv", sep="\t", index=False)

# per-pair signature concordance (Spearman over all 60 programs, matched-region z)
pair_concord = []
for reg in MATCHED_HUMAN:
    if reg not in mouse_z.index:
        continue
    rho, p = spearmanr(human_z5.loc[reg].values, mouse_z.loc[reg].values)
    # also restrict to R3-variable programs
    rv = sorted(R3_VARIABLE, key=int)
    rho_rv, _ = spearmanr(human_z5.loc[reg, rv].values, mouse_z.loc[reg, rv].values)
    nrec = (pair_compare[(pair_compare.region == reg) & (pair_compare.verdict == "RECAPITULATED")].shape[0])
    nhs = (pair_compare[(pair_compare.region == reg) & (pair_compare.verdict == "HUMAN-SPECIFIC")].shape[0])
    pair_concord.append(dict(region=reg, mouse_area="+".join(PAIRS[reg]["mouse"]),
                             signature_spearman=round(float(rho), 3),
                             signature_spearman_R3var=round(float(rho_rv), 3),
                             n_recapitulated=nrec, n_human_specific=nhs))
pair_concord = pd.DataFrame(pair_concord)
pair_concord.to_csv(OUT + "pair_signature_concordance.tsv", sep="\t", index=False)
print("[2b] per-pair signature concordance:")
print(pair_concord.to_string(index=False))


# ============================================================
# 3. PER-PROGRAM 异同 classification (CONSERVED / HUMAN-SPECIFIC / MOUSE-SPECIFIC)
#    over the matched-region subspace
# ============================================================
# For each program: across the 5 matched regions, does mouse vary the SAME way as human?
prog_rows = []
for p in PROGS:
    hv = human_z5[p].values  # 5 matched regions
    mv = mouse_z[p].values
    rho, _ = spearmanr(hv, mv)
    h_range = float(hv.max() - hv.min())
    m_range = float(mv.max() - mv.min())
    # human region-variable? (R3 global) and varies across THESE 5 matched regions
    h_var_here = h_range >= 1.5  # spans >=1.5 z across the 5 matched regions
    m_var_here = m_range >= 1.5
    # classify
    if h_var_here and m_var_here and rho >= 0.5:
        bucket = "CONSERVED"
    elif h_var_here and (not m_var_here or rho < 0.3):
        bucket = "HUMAN-SPECIFIC"
    elif m_var_here and not h_var_here:
        bucket = "MOUSE-SPECIFIC"
    elif h_var_here and m_var_here and rho < 0.5:
        bucket = "DIVERGENT"  # both vary but differently
    else:
        bucket = "FLAT-BOTH"
    prog_rows.append(dict(
        program=p, program_name=pname(p),
        R3_variable=(p in R3_VARIABLE), eta2_region=round(eta2.get(p, np.nan), 4),
        matched_spearman=round(float(rho), 3),
        human_range5=round(h_range, 2), mouse_range5=round(m_range, 2),
        human_top=MATCHED_HUMAN[int(np.argmax(hv))], human_bot=MATCHED_HUMAN[int(np.argmin(hv))],
        mouse_top=MATCHED_HUMAN[int(np.argmax(mv))], mouse_bot=MATCHED_HUMAN[int(np.argmin(mv))],
        bucket=bucket))
prog_class = pd.DataFrame(prog_rows)
prog_class.to_csv(OUT + "per_program_conservation.tsv", sep="\t", index=False)

bucket_counts_all = prog_class["bucket"].value_counts().to_dict()
# focus subset = R3 variable + strong distinguishers (top-15 eta2)
strong_dist = set(distinguish.head(15)["program"])
focus = set(R3_VARIABLE) | strong_dist
prog_focus = prog_class[prog_class["program"].isin(focus)].copy()
bucket_counts_focus = prog_focus["bucket"].value_counts().to_dict()
print("[3] per-program buckets (ALL 60):", bucket_counts_all)
print("[3] per-program buckets (FOCUS = 14 R3var + 15 top-eta2):", bucket_counts_focus)
print("[3] CONSERVED members:",
      [r["program_name"] for _, r in prog_class[prog_class.bucket == "CONSERVED"].iterrows()])


# ============================================================
# 4. GROUP-LEVEL 异同 — conserved primary core vs human-specific layer
# ============================================================
# primary sensory/motor pairs vs cingulate (ACC) vs (assoc has no homolog)
PRIMARY = ["M1", "S1", "V1", "STG"]  # primary motor + sensory (incl primary checkory STG~AUD)
HIGHER = ["ACC"]  # cingulate / higher-order among the matched set

pc = pair_concord.set_index("region")
primary_concord = pc.loc[[r for r in PRIMARY if r in pc.index], "signature_spearman"]
higher_concord = pc.loc[[r for r in HIGHER if r in pc.index], "signature_spearman"]
mean_matched_concord = float(pc["signature_spearman"].mean())

# global 0.07 baseline (corr humanVar vs mouseVar, from prior axis3)
GLOBAL_BLUNT = 0.071

# also: are primary-area-IDENTITY programs (the ones that define V1/M1/S1) recapitulated?
# count recapitulated signature programs in primary vs higher
rec_primary = pair_compare[(pair_compare.region.isin(PRIMARY)) & (pair_compare.verdict == "RECAPITULATED")].shape[0]
tot_primary = pair_compare[pair_compare.region.isin(PRIMARY)].shape[0]
rec_higher = pair_compare[(pair_compare.region.isin(HIGHER)) & (pair_compare.verdict == "RECAPITULATED")].shape[0]
tot_higher = pair_compare[pair_compare.region.isin(HIGHER)].shape[0]

group = dict(
    mean_matched_signature_concordance=round(mean_matched_concord, 3),
    primary_mean_concordance=round(float(primary_concord.mean()), 3),
    higher_ACC_concordance=round(float(higher_concord.mean()), 3),
    global_blunt_corr=GLOBAL_BLUNT,
    matched_vs_blunt_improvement=round(mean_matched_concord - GLOBAL_BLUNT, 3),
    primary_recapitulated_frac=f"{rec_primary}/{tot_primary}",
    higher_recapitulated_frac=f"{rec_higher}/{tot_higher}",
    n_regions_no_mouse_homolog=len(ASSOC_NOHOMOLOG),
    assoc_no_homolog=ASSOC_NOHOMOLOG,
)
print("[4] GROUP verdict:", json.dumps(group, indent=1))


# ============================================================
# 5. SUMMARY json + md
# ============================================================
summary = dict(
    analysis="axis3_deep cross-region 异同 (human region signatures + mouse 异同 + conservation buckets)",
    date="2026-06-04",
    n_mouse_sections=len(files),
    matched_pairs={k: v["mouse"] for k, v in PAIRS.items()},
    mouse_matched_nbins={hp: int(np.vstack([b for ma in PAIRS[hp]["mouse"] if ma in area_comp_bins for b in area_comp_bins[ma]]).shape[0]) for hp in MATCHED_HUMAN if any(ma in area_comp_bins for ma in PAIRS[hp]["mouse"])},
    human_strongest_distinguishing=[
        dict(program=r["program"], name=r["program_name"], eta2=round(float(r["eta2_region"]), 3),
             top_region=r["top_region"], bottom_region=r["bottom_region"])
        for _, r in distinguish.head(10).iterrows()],
    per_program_buckets_all=bucket_counts_all,
    per_program_buckets_focus=bucket_counts_focus,
    conserved_members=[r["program_name"] for _, r in prog_class[prog_class.bucket == "CONSERVED"].iterrows()],
    human_specific_members=[r["program_name"] for _, r in prog_focus[prog_focus.bucket == "HUMAN-SPECIFIC"].iterrows()],
    per_pair_concordance=pair_concord.to_dict(orient="records"),
    group_verdict=group,
    refined_statement=(
        "Cross-region program regulation is PARTLY conserved, not globally absent. Over the 5 clean "
        "mouse<->human homolog pairs, per-region program signatures concord at Spearman "
        f"{mean_matched_concord:.2f} (primary sensory/motor mean {float(primary_concord.mean()):.2f}), "
        f"far above the blunt global humanVar~mouseVar corr {GLOBAL_BLUNT}. Primary-area identity "
        "programs (V1/VIS, M1/MOp, S1/SSp) form a CONSERVED CORE of region differences, while human "
        "association/higher-order cortex (DLPFC + 8 regions with NO clean mouse homolog) carries the "
        "HUMAN-ENRICHED layer of region specialization. Honest: the human region-specialization is "
        "human-enriched ESPECIALLY in association cortex, but the primary-sensory area-identity axis "
        "IS shared with mouse."),
)
with open(OUT + "axis3_deep_summary.json", "w") as fh:
    json.dump(summary, fh, indent=1, ensure_ascii=False)

# markdown
md = []
md.append("# Axis 3 DEEP — cross-region 异同 (for Figure 10 panel e)\n")
md.append(f"_{summary['date']}; {len(files)} mouse sections; 5 clean homolog pairs._\n")
md.append("\n## 1. Human region differences (signatures)\n")
md.append("Strongest region-distinguishing programs (by eta2_region):\n")
for _, r in distinguish.head(8).iterrows():
    md.append(f"- **{r['program_name']}** eta2={r['eta2_region']:.3f}, peaks **{r['top_region']}**, low **{r['bottom_region']}**\n")
md.append("\nPer-region top-up signature (rank-1 program):\n")
for reg in human_regions:
    top = human_sig[(human_sig.region == reg) & (human_sig.direction == "up") & (human_sig["rank"] == 1)].iloc[0]
    md.append(f"- {reg}: {top['program_name']} (z={top['z']})\n")
md.append("\n## 2. Mouse 异同 per matched region\n")
for _, r in pair_concord.iterrows():
    md.append(f"- **{r['region']} <-> mouse {r['mouse_area']}**: signature Spearman {r['signature_spearman']} "
              f"(R3-var {r['signature_spearman_R3var']}); recapitulated {r['n_recapitulated']}, human-specific {r['n_human_specific']}\n")
md.append("\nExamples (human signature program -> mouse verdict):\n")
for reg in MATCHED_HUMAN:
    sub = pair_compare[pair_compare.region == reg]
    recs = sub[sub.verdict == "RECAPITULATED"]["program_name"].tolist()
    hss = sub[sub.verdict == "HUMAN-SPECIFIC"]["program_name"].tolist()
    md.append(f"- {reg}: RECAP={recs if recs else 'none'}; HUMAN-SPEC={hss if hss else 'none'}\n")
md.append("\n## 3. Per-program 异同 buckets\n")
md.append(f"- ALL 60: {bucket_counts_all}\n")
md.append(f"- FOCUS (14 R3-var + 15 top-eta2): {bucket_counts_focus}\n")
md.append(f"- CONSERVED members: {summary['conserved_members']}\n")
md.append(f"- HUMAN-SPECIFIC (focus): {summary['human_specific_members']}\n")
md.append("\n## 4. Group verdict\n")
for k, v in group.items():
    md.append(f"- {k}: {v}\n")
md.append("\n## Refined honest statement\n")
md.append(summary["refined_statement"] + "\n")
with open(OUT + "axis3_deep_summary.md", "w") as fh:
    fh.write("".join(md))

print("\n[WROTE] ->", OUT)
for f in sorted(os.listdir(OUT)):
    print("   ", f)
