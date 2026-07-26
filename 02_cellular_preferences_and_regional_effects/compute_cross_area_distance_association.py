#!/usr/bin/env python
import os, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
BASE = os.environ.get('CORTEX_NMF_ROOT', os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
OUT = os.environ.get('XREGION_OUTPUT_DIR', f'{BASE}/results/crossregion_v1/p37_current54')
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(0)

def perm_p_spearman(x, y, n=10000):
    (rho0, _) = spearmanr(x, y)
    cnt = 0
    for _ in range(n):
        yp = RNG.permutation(y)
        (r, _) = spearmanr(x, yp)
        if abs(r) >= abs(rho0):
            cnt += 1
    return (rho0, (cnt + 1) / (n + 1))

def boot_ci_spearman(x, y, n=5000):
    rs = []
    idx = np.arange(len(x))
    for _ in range(n):
        s = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(x[s])) < 3:
            continue
        (r, _) = spearmanr(x[s], y[s])
        rs.append(r)
    rs = np.array(rs)
    return (np.nanpercentile(rs, 2.5), np.nanpercentile(rs, 97.5))

def main():
    pc = pd.read_csv(f'{BASE}/results/crossregion_v1/region_pc_coords.tsv', sep='\t', index_col=0)
    regions = list(pc.index)
    try:
        rag = pd.read_csv(f'{BASE}/results/crossregion_v1/region_axis_gradient.tsv', sep='\t')
        ap_order = pc['PC1'].rank().to_dict()
    except Exception:
        ap_order = pc['PC1'].rank().to_dict()
    expr = pd.read_csv(f'{OUT}/m1_expr_pairwise_auroc_matrix.tsv', sep='\t', index_col=0)
    neigh = pd.read_csv(f'{OUT}/m2b_neigh_pairwise_auroc_matrix.tsv', sep='\t', index_col=0)
    rows = []
    for (rA, rB) in itertools.combinations(regions, 2):
        td = float(np.sqrt(((pc.loc[rA] - pc.loc[rB]) ** 2).sum()))
        ad = abs(ap_order[rA] - ap_order[rB])
        ea = np.nanmean([expr.loc[rA, rB], expr.loc[rB, rA]])
        na = np.nanmean([neigh.loc[rA, rB], neigh.loc[rB, rA]])
        rows.append((rA, rB, td, ad, ea, na))
    pts = pd.DataFrame(rows, columns=['regionA', 'regionB', 'trans_dist', 'ap_dist', 'expr_auroc', 'neigh_auroc'])
    pts.to_csv(f'{OUT}/m3_distance_curve_points.tsv', sep='\t', index=False)
    stat_rows = []
    for (ymet, ycol) in [('expr', 'expr_auroc'), ('neigh', 'neigh_auroc')]:
        for (dname, dcol) in [('transcriptomic', 'trans_dist'), ('anatomical_AP', 'ap_dist')]:
            x = pts[dcol].to_numpy(float)
            y = pts[ycol].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            (x, y) = (x[ok], y[ok])
            (rho, pp) = perm_p_spearman(x, y, n=10000)
            (lo, hi) = boot_ci_spearman(x, y, n=5000)
            stat_rows.append((ymet, dname, rho, pp, lo, hi, len(x)))
    pd.DataFrame(stat_rows, columns=['conservation_metric', 'distance', 'spearman_rho', 'perm_p', 'ci_lo', 'ci_hi', 'n_pairs']).to_csv(f'{OUT}/m3_distance_curve_stats.tsv', sep='\t', index=False)
if __name__ == '__main__':
    main()
