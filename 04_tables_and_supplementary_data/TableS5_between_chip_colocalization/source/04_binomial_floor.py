#!/usr/bin/env python
"""04_binomial_floor.py — Distribution-free binomial p floor.

For each pair, compute one-sided binomial test:
  H0: P(same sign) = 0.5, n = n_chips (44 typically)
  p_binom = P(K >= k_same | Binomial(n, 0.5))

Final reported p: p_final = min(p_stouffer, p_binom)

Also flags pairs where p_stouffer is at permutation floor (< 1/n_perm = 1e-3)
→ in those cases p_binom is the operative p.

Reads: betweenchip_{mode}_stouffer_q.tsv (from script 03)
Writes: betweenchip_{mode}_stouffer_q.tsv (overwrites with new columns added)
  new columns: p_binom, p_final, at_perm_floor

Usage:
  python 04_binomial_floor.py [--smoke]
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
from scipy.stats import binomtest

RESDIR = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUTDIR = os.path.join(RESDIR, "markcorr_betweenchip_v1")

PERM_FLOOR = 1e-3   # 1/n_perm; p_stouffer below this may be at floor

def compute_binom_p(k: int, n: int) -> float:
    """One-sided binomial: P(K >= k | Binom(n, 0.5))."""
    if k <= 0:
        return 1.0
    return float(binomtest(int(k), int(n), 0.5, alternative='greater').pvalue)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    suffix = '_smoke' if args.smoke else ''

    for mode in ['cellprog', 'progprog']:
        in_f = os.path.join(OUTDIR, f"betweenchip_{mode}_stouffer_q{suffix}.tsv")
        if not os.path.exists(in_f):
            print(f"Skipping {mode}: {in_f} not found")
            continue
        print(f"[{time.strftime('%H:%M:%S')}] Binomial floor for {mode}")
        df = pd.read_csv(in_f, sep='\t')

        # Compute binomial p
        df['p_binom'] = [compute_binom_p(row['n_same_sign'], row['n_chips'])
                         for _, row in df.iterrows()]
        df['at_perm_floor'] = df['p_stouffer'] < PERM_FLOOR
        df['p_final'] = np.minimum(df['p_stouffer'], df['p_binom'])

        # Re-apply headline gate using p_final (for q-cut we keep q_bh on p_stouffer;
        # p_final is reported in TableS3; headline uses q_bh from BH of p_stouffer)
        # Note: p_binom not FDR-corrected separately (spec §10.3 says "p_final = min(p_s, p_b)")
        df['p_final_note'] = np.where(
            df['at_perm_floor'],
            'binomial_operative',
            'stouffer_operative'
        )

        df.to_csv(in_f, sep='\t', index=False)

        n_floor = df['at_perm_floor'].sum()
        n_44_44 = (df['n_same_sign'] == df['n_chips']).sum()
        print(f"  {mode}: {n_floor} pairs at perm floor | {n_44_44} with n_same=n_chips")
        if n_44_44 > 0:
            ex = df[df['n_same_sign'] == df['n_chips']].iloc[0]
            print(f"  Example 44/44 pair: {ex['A_name']} x {ex['B_name']} "
                  f"p_binom={ex['p_binom']:.3e}")
        print(f"  Overwrote -> {in_f}")

    print(f"[{time.strftime('%H:%M:%S')}] DONE")

if __name__ == '__main__':
    main()
