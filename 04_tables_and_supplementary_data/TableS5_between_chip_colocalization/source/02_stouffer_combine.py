#!/usr/bin/env python
"""02_stouffer_combine.py — Between-chip Stouffer combination.

Reads per_chip_Z TSV from script 01, combines across 44 chips:
  Z_combined = sum(w_i * Z_i) / sqrt(sum(w_i^2))   with w_i = sqrt(n_in_tissue_i)
  p_stouffer = 2 * Phi^c(|Z_combined|)  (two-sided)

Also computes:
  - frac_same_sign: fraction of chips with same sign as Z_combined
  - n_same_sign: count
  - median_log2g: median log2 g(25µm) across chips
  - iqr_log2g: IQR
  - Z_unweighted: unweighted Stouffer (sensitivity check)

Output: betweenchip_{mode}_stouffer_prelim.tsv
  (p-values without BH; BH done in script 03)

Usage:
  python 02_stouffer_combine.py --mode cellprog [--smoke]
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
from scipy.stats import norm

RESDIR = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
OUTDIR = os.path.join(RESDIR, "markcorr_betweenchip_v1")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', required=True, choices=['cellprog', 'progprog'])
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    suffix = '_smoke' if args.smoke else ''
    in_f = os.path.join(OUTDIR, f"betweenchip_{args.mode}_per_chip_Z{suffix}.tsv")
    if not os.path.exists(in_f):
        # try without suffix (e.g. shard outputs concatenated)
        in_f = os.path.join(OUTDIR, f"betweenchip_{args.mode}_per_chip_Z.tsv")
    print(f"[{time.strftime('%H:%M:%S')}] Reading {in_f}")
    df = pd.read_csv(in_f, sep='\t')
    print(f"  {len(df)} rows, {df['chip_id'].nunique()} chips, "
          f"{df['A_name'].nunique()} A labels, {df['B_name'].nunique()} B labels")

    rows = []
    for (mode, An, Bn), grp in df.groupby(['mode', 'A_name', 'B_name']):
        Zi = grp['Z_i'].values.astype(np.float64)
        n_bin = grp['n_in_tissue'].values.astype(np.float64)
        log2g = grp['log2_g_obs'].values.astype(np.float64)
        chips = grp['chip_id'].values
        n_chips = len(grp)

        w = np.sqrt(np.maximum(n_bin, 1.0))

        # Weighted Stouffer
        Z_w = np.dot(w, Zi) / np.sqrt(np.dot(w, w))
        p_w = float(2 * norm.sf(np.abs(Z_w)))

        # Unweighted Stouffer (sensitivity)
        Z_uw = Zi.sum() / np.sqrt(n_chips)
        p_uw = float(2 * norm.sf(np.abs(Z_uw)))

        # Same-sign stats
        dominant_sign = np.sign(np.nanmedian(log2g))
        same_sign = (np.sign(Zi) == dominant_sign).sum()
        frac_same = float(same_sign / n_chips)

        # Effect size
        med_g = float(np.nanmedian(log2g))
        q25, q75 = np.nanpercentile(log2g, [25, 75])
        iqr_g = float(q75 - q25)

        rows.append({
            'mode': mode,
            'A_name': An,
            'B_name': Bn,
            'n_chips': n_chips,
            'Z_combined': float(Z_w),
            'p_stouffer': p_w,
            'Z_combined_unweighted': float(Z_uw),
            'p_stouffer_unweighted': p_uw,
            'median_log2g': med_g,
            'iqr_log2g': iqr_g,
            'frac_same_sign': frac_same,
            'n_same_sign': int(same_sign),
        })

    out_df = pd.DataFrame(rows)
    out_f = os.path.join(OUTDIR, f"betweenchip_{args.mode}_stouffer_prelim{suffix}.tsv")
    out_df.to_csv(out_f, sep='\t', index=False)
    print(f"Saved {len(out_df)} rows -> {out_f}")

    # Quick summary stats
    sig = (out_df['p_stouffer'] < 0.05).sum()
    print(f"Pairs with p_stouffer < 0.05 (before BH): {sig}/{len(out_df)}")
    print(f"Z_combined range: [{out_df['Z_combined'].min():.2f}, {out_df['Z_combined'].max():.2f}]")
    print(f"[{time.strftime('%H:%M:%S')}] DONE")

if __name__ == '__main__':
    main()
