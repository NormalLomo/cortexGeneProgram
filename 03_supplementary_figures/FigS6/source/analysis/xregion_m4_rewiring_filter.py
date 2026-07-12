#!/usr/bin/env python
"""
Metric 4: Multi-step filter for region-specific REWIRING hits.

A program is a "rewiring candidate" when its EXPRESSION composition is conserved
across regions but its CO-ACTIVATION neighborhood is NOT (partners turn over),
with spatial co-localization corroborating, and anchored to variability/signature
tables.

Gate A (expression conserved high): m1 expr_cons_auroc >= A_thr (top tertile / >=0.90)
Gate B (neighborhood conserved low): m2b neigh_cons_auroc <= B_thr (bottom tertile)
Gate C (co-activation + spatial agree on rewiring direction):
        program is among lower-conserved in BOTH m2b and m2a (rank agreement),
        AND has high partner turnover.
Gate D (anchor): program is in program_variability 'variable' set OR appears in
        supp_table_region_signatures (region-distinguishing).

For each hit, identify the target region(s) where its neighborhood deviates most
(lowest mean self-AUROC as regionB) and report partner turnover + anchors.

Inputs (xregion_auroc/): m1_*, m2b_*, m2a_*
  + program_variability.tsv, supp_table_region_signatures.tsv, program_names.tsv
Output (xregion_auroc/):
  m4_program_master_table.tsv   per-program all scores + gates
  m4_rewiring_hits.tsv          hit short table
"""
import os
import numpy as np
import pandas as pd

BASE = os.environ.get("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
RES  = f"{BASE}/results/crossregion_v1"
OUT  = os.environ.get("XREGION_OUTPUT_DIR", f"{RES}/xregion_auroc")
RETAIN_MAP = os.environ.get("RETAIN_MAP", f"{BASE}/tables/TableS3_program_annotation.tsv")
os.makedirs(OUT, exist_ok=True)

def main():
    mapping = pd.read_csv(RETAIN_MAP, sep="\t")
    assert len(mapping) == 54
    mapping["old_int"] = mapping["cnmf_component"].astype(int)
    mapping["new_int"] = mapping["new_P"].astype(str).str.removeprefix("P").astype(int)
    assert mapping["new_int"].tolist() == list(range(1, 55))
    old_to_new = dict(zip(mapping["old_int"], mapping["new_int"]))
    m1 = pd.read_csv(f"{OUT}/m1_expr_conservation_per_program.tsv", sep="\t")
    m2b= pd.read_csv(f"{OUT}/m2b_neighborhood_conservation_per_program.tsv", sep="\t")
    m2a= pd.read_csv(f"{OUT}/m2a_spatial_neighborhood_conservation_per_program.tsv", sep="\t")
    turn=pd.read_csv(f"{OUT}/m2b_partner_turnover.tsv", sep="\t")
    var = pd.read_csv(f"{RES}/program_variability.tsv", sep="\t")[["program","class","eta2_region","fdr"]]
    var["program"] = var["program"].map(old_to_new)
    var = var.dropna(subset=["program"]); var["program"] = var["program"].astype(int)
    names=pd.read_csv(f"{RES}/program_names.tsv", sep="\t")[["cnmf_component","name_short","name_full","confidence","fdr"]]
    names["program"] = names["cnmf_component"].map(old_to_new)
    names = names.dropna(subset=["program"]); names["program"] = names["program"].astype(int)
    names = names.drop(columns=["cnmf_component"])
    names=names.rename(columns={"fdr":"name_fdr"})
    sig = pd.read_csv(f"{RES}/supp_table_region_signatures.tsv", sep="\t")

    M = m1.merge(m2b, on="program", suffixes=("","_b")).merge(turn, on="program")
    M = M.merge(m2a[["program","spatial_neigh_cons_auroc"]], on="program", how="left")
    M = M.merge(var, on="program", how="left").merge(names, on="program", how="left")

    # signature membership
    sigset = {old_to_new[p] for p in sig["program"].unique() if p in old_to_new}
    M["in_region_signature"] = M["program"].isin(sigset)

    # gates
    A_thr = M["expr_cons_auroc"].quantile(2/3)      # top tertile
    B_thr = M["neigh_cons_auroc"].quantile(1/3)      # bottom tertile
    M["gateA_expr_conserved"] = M["expr_cons_auroc"] >= max(A_thr, 0.90)
    M["gateB_neigh_low"]      = M["neigh_cons_auroc"] <= B_thr
    # gate C: rank agreement low-conservation in both 2b and 2a + high turnover
    M["rank_neigh"]   = M["neigh_cons_auroc"].rank()           # low rank = low conservation
    M["rank_spatial"] = M["spatial_neigh_cons_auroc"].rank()
    turn_thr = M["partner_turnover"].quantile(2/3)
    n_programs = len(M)
    assert n_programs == 54
    neigh_rank_cut = int(np.ceil(0.50 * n_programs))
    spatial_rank_cut = int(np.ceil((35 / 60) * n_programs))
    M["gateC_coact_spatial_agree"] = (
        (M["rank_neigh"] <= neigh_rank_cut) &
        (M["rank_spatial"] <= spatial_rank_cut) &
        (M["partner_turnover"] >= turn_thr)
    )
    M["gateD_anchor"] = (M["class"]=="variable") | (M["in_region_signature"])

    M["rewiring_gap"] = M["expr_cons_auroc"] - M["neigh_cons_auroc"]
    M["n_gates_pass"] = M[["gateA_expr_conserved","gateB_neigh_low","gateC_coact_spatial_agree","gateD_anchor"]].sum(1)
    M = M.sort_values(["n_gates_pass","rewiring_gap"], ascending=False)
    M.to_csv(f"{OUT}/m4_program_master_table.tsv", sep="\t", index=False)

    # Strict sequential hits. A relaxed result is never substituted after a zero stage.
    hits = M[(M["gateA_expr_conserved"]) & (M["gateB_neigh_low"]) &
             (M["gateC_coact_spatial_agree"]) & (M["gateD_anchor"])].copy()
    fallback = M[(M["gateA_expr_conserved"]) & (M["gateB_neigh_low"]) &
                 (M["n_gates_pass"] >= 3)].copy()
    fallback.to_csv(f"{OUT}/m4_rewiring_fallback_hits.tsv", sep="\t", index=False)

    # target region(s): regions (as regionB) where neighborhood self-AUROC is lowest
    neighlong = pd.read_csv(f"{OUT}/m2b_neigh_program_region_self_auroc.tsv", sep="\t")
    tgt_rows=[]
    for p in hits["program"]:
        sub = neighlong[neighlong["program"]==p].groupby("regionB")["self_auroc"].mean().sort_values()
        tgt_rows.append((p, ";".join(sub.index[:3]), round(float(sub.iloc[0]),3)))
    tgtdf=pd.DataFrame(tgt_rows,columns=["program","rewiring_target_regions","min_region_auroc"])
    hits=hits.merge(tgtdf,on="program",how="left")

    cols=["program","name_short","name_full","confidence","name_fdr","class","eta2_region",
          "expr_cons_auroc","neigh_cons_auroc","spatial_neigh_cons_auroc","rewiring_gap",
          "partner_turnover","mean_partner_jaccard","in_region_signature",
          "rewiring_target_regions","min_region_auroc","n_gates_pass"]
    cols=[c for c in cols if c in hits.columns]
    hits[cols].to_csv(f"{OUT}/m4_rewiring_hits.tsv", sep="\t", index=False)

    print(f"[m4] A_thr(expr>=)={max(A_thr,0.90):.3f}  B_thr(neigh<=)={B_thr:.3f}  turn_thr={turn_thr:.3f}", flush=True)
    print(f"[m4] rank cutoffs: neigh<={neigh_rank_cut}, spatial<={spatial_rank_cut}", flush=True)
    print(f"[m4] strict hits={len(hits)} separate fallback candidates={len(fallback)}", flush=True)
    print(hits[cols].to_string(index=False))
    # funnel counts
    print("\n[m4] FUNNEL:")
    print(f"  all programs           : {len(M)}")
    print(f"  Gate A (expr conserved): {M['gateA_expr_conserved'].sum()}")
    print(f"  +Gate B (neigh low)    : {(M['gateA_expr_conserved']&M['gateB_neigh_low']).sum()}")
    print(f"  +Gate C (coact+spatial): {(M['gateA_expr_conserved']&M['gateB_neigh_low']&M['gateC_coact_spatial_agree']).sum()}")
    print(f"  +Gate D (anchor) = HITS: {len(hits)}")
    # save funnel
    fn=pd.DataFrame({
        "stage":["all_programs","gateA_expr_conserved","+gateB_neigh_low","+gateC_coact_spatial","+gateD_anchor_HITS"],
        "n":[len(M), int(M['gateA_expr_conserved'].sum()),
             int((M['gateA_expr_conserved']&M['gateB_neigh_low']).sum()),
             int((M['gateA_expr_conserved']&M['gateB_neigh_low']&M['gateC_coact_spatial_agree']).sum()),
             len(hits)]})
    fn.to_csv(f"{OUT}/m4_filter_funnel.tsv", sep="\t", index=False)
    print("[m4] DONE", flush=True)

if __name__=="__main__":
    main()
