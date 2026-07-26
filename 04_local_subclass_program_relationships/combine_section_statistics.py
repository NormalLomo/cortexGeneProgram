#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, norm
from statsmodels.stats.multitest import multipletests
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
OUTPUT = ROOT / 'results/crossregion_v1/markcorr_betweenchip_v1'
EFFECT_THRESHOLD = 0.32
SIGN_THRESHOLD = 0.85

def calculate_mode(mode):
    source = OUTPUT / f'betweenchip_{mode}_per_chip_Z.tsv'
    table = pd.read_csv(source, sep='\t')
    rows = []
    for ((label_a, label_b), group) in table.groupby(['A_name', 'B_name'], sort=False):
        z_values = group['Z_i'].to_numpy(dtype=np.float64)
        sizes = group['n_in_tissue'].to_numpy(dtype=np.float64)
        effects = group['log2_g_obs'].to_numpy(dtype=np.float64)
        weights = np.sqrt(np.maximum(sizes, 1.0))
        weighted_z = float(np.dot(weights, z_values) / np.sqrt(np.dot(weights, weights)))
        unweighted_z = float(z_values.sum() / np.sqrt(len(z_values)))
        dominant_sign = np.sign(np.nanmedian(effects))
        same_sign = int((np.sign(z_values) == dominant_sign).sum())
        (lower, upper) = np.nanpercentile(effects, [25, 75])
        rows.append({'mode': mode, 'A_name': label_a, 'B_name': label_b, 'n_chips': len(group), 'Z_combined': weighted_z, 'p_stouffer': float(2 * norm.sf(abs(weighted_z))), 'Z_combined_unweighted': unweighted_z, 'p_stouffer_unweighted': float(2 * norm.sf(abs(unweighted_z))), 'median_log2g': float(np.nanmedian(effects)), 'iqr_log2g': float(upper - lower), 'frac_same_sign': float(same_sign / len(group)), 'n_same_sign': same_sign})
    result = pd.DataFrame(rows)
    if mode == 'progprog':
        first = result['A_name'].str.replace('program_', '', regex=False).astype(int)
        second = result['B_name'].str.replace('program_', '', regex=False).astype(int)
        result = result[first < second].copy()
    (_, q_values, _, _) = multipletests(result['p_stouffer'].to_numpy(dtype=float).clip(0, 1), method='fdr_bh')
    result['q_bh'] = q_values
    result['is_headline'] = result['q_bh'].lt(0.05) & result['median_log2g'].abs().gt(EFFECT_THRESHOLD) & result['frac_same_sign'].ge(SIGN_THRESHOLD)
    result['p_binom'] = [float(binomtest(int(row.n_same_sign), int(row.n_chips), 0.5, alternative='greater').pvalue) for row in result.itertuples()]
    result['at_perm_floor'] = result['p_stouffer'].lt(0.001)
    result['p_final'] = np.minimum(result['p_stouffer'], result['p_binom'])
    result['p_final_note'] = np.where(result['p_stouffer'].le(result['p_binom']), 'stouffer_operative', 'binomial_operative')
    required_rows = 22 * 54 if mode == 'cellprog' else 54 * 53 // 2
    if len(result) != required_rows:
        raise ValueError()
    result = result.sort_values('q_bh').reset_index(drop=True)
    result.to_csv(OUTPUT / f'betweenchip_{mode}_stouffer_q.tsv', sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--mode', choices=['cellprog', 'progprog', 'both'], default='both')
    arguments = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    modes = ['cellprog', 'progprog'] if arguments.mode == 'both' else [arguments.mode]
    for mode in modes:
        calculate_mode(mode)
if __name__ == '__main__':
    main()
