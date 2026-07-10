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
import numpy as np
import pandas as pd

BASE = "CORTEX_PROGRAM_ROOT"
RES  = f"{BASE}/results/crossregion_v1"
OUT  = f"{RES}/xregion_auroc"

def main():
    m1 = pd.read_csv(f"{OUT}/m1_expr_conservation_per_program.tsv", sep="\t")
    m2b= pd.read_csv(f"{OUT}/m2b_neighborhood_conservation_per_program.tsv", sep="\t")
    m2a= pd.read_csv(f"{OUT}/m2a_spatial_neighborhood_conservation_per_program.tsv", sep="\t")
    turn=pd.read_csv(f"{OUT}/m2b_partner_turnover.tsv", sep="\t")
    var = pd.read_csv(f"{RES}/program_variability.tsv", sep="\t")[["program","class","eta2_region","fdr"]]
    names=pd.read_csv(f"{RES}/program_names.tsv", sep="\t")[["program","name_short","name_full","confidence","fdr"]]
    names=names.rename(columns={"fdr":"name_fdr"})
    sig = pd.read_csv(f"{RES}/supp_table_region_signatures.tsv", sep="\t")

    M = m1.merge(m2b, on="program", suffixes=("","_b")).merge(turn, on="program")
    M = M.merge(m2a[["program","spatial_neigh_cons_auroc"]], on="program", how="left")
    M = M.merge(var, on="program", how="left").merge(names, on="program", how="left")

    # signature membership
    sigset = set(sig["program"].unique())
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
    M["gateC_coact_spatial_agree"] = (
        (M["rank_neigh"] <= 30) & (M["rank_spatial"] <= 35) & (M["partner_turnover"] >= turn_thr)
    )
    M["gateD_anchor"] = (M["class"]=="variable") | (M["in_region_signature"])

    M["rewiring_gap"] = M["expr_cons_auroc"] - M["neigh_cons_auroc"]
    M["n_gates_pass"] = M[["gateA_expr_conserved","gateB_neigh_low","gateC_coact_spatial_agree","gateD_anchor"]].sum(1)
    M = M.sort_values(["n_gates_pass","rewiring_gap"], ascending=False)
    M.to_csv(f"{OUT}/m4_program_master_table.tsv", sep="\t", index=False)

    # hits = pass all 4 gates (or >=3 with A&B mandatory)
    hits = M[(M["gateA_expr_conserved"]) & (M["gateB_neigh_low"]) &
             (M["gateC_coact_spatial_agree"]) & (M["gateD_anchor"])].copy()
    if len(hits)==0:
        hits = M[(M["gateA_expr_conserved"]) & (M["gateB_neigh_low"]) & (M["n_gates_pass"]>=3)].copy()

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
    print(f"[m4] n hits = {len(hits)}", flush=True)
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
