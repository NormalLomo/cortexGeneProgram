#!/usr/bin/env python
"""03_bh_fdr.py — BH-FDR across all pairs within each hypothesis space.

Reads stouffer_prelim TSVs (from script 02), applies BH-FDR separately within:
  - cellprog: 22 × 54 = 1,188 pairs
  - progprog: C(54,2) = 1,431 upper-triangle pairs (A_name < B_name)

Outputs:
  betweenchip_{mode}_stouffer_q.tsv — full table with q_bh column
  betweenchip_{mode}_headline.tsv  — q<0.05 AND |log2g|>0.32 AND frac_same≥0.85

Effect-size threshold: |log2 g| > 0.32 (same as v10 production; spec §5)
Same-sign gate: frac_same_sign >= 0.85 (>=37/44; spec §5)

Usage:
  python 03_bh_fdr.py [--smoke]
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

RESDIR = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUTDIR = os.path.join(RESDIR, "markcorr_betweenchip_v1")

ES_THRESH = 0.32   # |log2g| headline threshold
FS_THRESH = 0.85   # frac_same_sign headline gate

def apply_bh(df_in: pd.DataFrame, mode: str, suffix: str) -> pd.DataFrame:
    df = df_in.copy()

    # For progprog: only upper triangle (A_name < B_name) to avoid double-counting
    if mode == 'progprog':
        df = df[df['A_name'] < df['B_name']].copy()

    # BH over p_stouffer within this hypothesis space
    pvals = df['p_stouffer'].values.clip(0, 1)
    _, q_bh, _, _ = multipletests(pvals, method='fdr_bh')
    df['q_bh'] = q_bh

    # Headline gate
    df['is_headline'] = (
        (df['q_bh'] < 0.05) &
        (df['median_log2g'].abs() > ES_THRESH) &
        (df['frac_same_sign'] >= FS_THRESH)
    )

    # Sort by q_bh
    df = df.sort_values('q_bh').reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    suffix = '_smoke' if args.smoke else ''

    for mode in ['cellprog', 'progprog']:
        in_f = os.path.join(OUTDIR, f"betweenchip_{mode}_stouffer_prelim{suffix}.tsv")
        if not os.path.exists(in_f):
            # smoke may have only cellprog
            print(f"Skipping {mode}: {in_f} not found")
            continue
        print(f"[{time.strftime('%H:%M:%S')}] BH-FDR for {mode} from {in_f}")
        df = pd.read_csv(in_f, sep='\t')
        df_out = apply_bh(df, mode, suffix)

        out_f = os.path.join(OUTDIR, f"betweenchip_{mode}_stouffer_q{suffix}.tsv")
        df_out.to_csv(out_f, sep='\t', index=False)

        head_f = os.path.join(OUTDIR, f"betweenchip_{mode}_headline{suffix}.tsv")
        df_out[df_out['is_headline']].to_csv(head_f, sep='\t', index=False)

        n_total = len(df_out)
        n_q05   = (df_out['q_bh'] < 0.05).sum()
        n_head  = df_out['is_headline'].sum()
        print(f"  {mode}: {n_total} pairs tested | q<0.05: {n_q05} | headline: {n_head}")
        print(f"  Saved -> {out_f}")
        print(f"  Saved -> {head_f}")

    print(f"[{time.strftime('%H:%M:%S')}] DONE")

if __name__ == '__main__':
    main()
