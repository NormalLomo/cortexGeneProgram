#!/usr/bin/env python
import os, sys, itertools
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import rankdata
BASE = os.environ.get('CORTEX_NMF_ROOT', os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
OUT = os.environ.get('XREGION_OUTPUT_DIR', f'{BASE}/results/crossregion_v1/p37_current54')
PARQ = f'{BASE}/results/crossregion_v1/cell_program_region_subclass.parquet'
RETAIN_MAP = os.environ.get('RETAIN_MAP', f'{BASE}/results/crossregion_v1/program_renumber_map.tsv')
os.makedirs(OUT, exist_ok=True)
K = 10

def spearman_corr_matrix(X):
    (n, P) = X.shape
    R = np.empty_like(X, dtype=float)
    for j in range(P):
        R[:, j] = rankdata(X[:, j])
    R -= R.mean(0)
    cov = R.T @ R
    d = np.sqrt(np.diag(cov))
    denom = np.outer(d, d)
    with np.errstate(invalid='ignore', divide='ignore'):
        C = cov / denom
    C[~np.isfinite(C)] = 0.0
    np.fill_diagonal(C, 1.0)
    return C

def main():
    mapping = pd.read_csv(RETAIN_MAP, sep='\t')
    mapping = mapping[mapping['status'].astype(str).str.lower().eq('kept') & mapping['new_P'].notna()].copy()
    mapping['old_int'] = mapping['old_P'].astype(int)
    mapping['new_int'] = mapping['new_P'].astype(str).str.removeprefix('P').astype(int)
    progs = mapping['old_int'].astype(str).tolist()
    tbl = pq.read_table(PARQ, columns=progs + ['region'])
    df = tbl.to_pandas()
    regions = sorted(df['region'].unique())
    P = len(progs)
    progint = mapping['new_int'].tolist()
    C = {}
    for r in regions:
        Xr = df.loc[df['region'] == r, progs].to_numpy(float)
        C[r] = spearman_corr_matrix(Xr)
        pd.DataFrame(C[r], index=progint, columns=progint).to_csv(f'{OUT}/m2b_coact_corr_{r}.tsv', sep='\t')
    topk = {}
    for r in regions:
        Cr = C[r]
        d = {}
        for pi in range(P):
            v = Cr[pi].copy()
            v[pi] = -np.inf
            d[pi] = set(np.argsort(-v)[:K].tolist())
        topk[r] = d
    self_rows = []
    ridx = {r: i for (i, r) in enumerate(regions)}
    Rn = len(regions)
    pair_sum = np.zeros((Rn, Rn))
    pair_n = np.zeros((Rn, Rn))
    for (rA, rB) in itertools.permutations(regions, 2):
        CB = C[rB]
        for pi in range(P):
            partners = topk[rA][pi]
            scores = CB[pi].copy()
            mask = np.ones(P, bool)
            mask[pi] = False
            qs = np.where(mask)[0]
            lab = np.array([1 if q in partners else 0 for q in qs])
            sc = scores[qs]
            npos = lab.sum()
            nneg = len(lab) - npos
            if npos == 0 or nneg == 0:
                auc = np.nan
            else:
                order = rankdata(sc)
                auc = (order[lab == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)
            self_rows.append((progint[pi], rA, rB, auc))
            if np.isfinite(auc):
                pair_sum[ridx[rA], ridx[rB]] += auc
                pair_n[ridx[rA], ridx[rB]] += 1
    selfdf = pd.DataFrame(self_rows, columns=['program', 'regionA', 'regionB', 'self_auroc'])
    selfdf.to_csv(f'{OUT}/m2b_neigh_program_region_self_auroc.tsv', sep='\t', index=False)
    perprog = selfdf.groupby('program')['self_auroc'].agg(['mean', 'count']).reset_index()
    perprog.columns = ['program', 'neigh_cons_auroc', 'n_pairs']
    perprog = perprog.sort_values('neigh_cons_auroc', ascending=False)
    perprog.to_csv(f'{OUT}/m2b_neighborhood_conservation_per_program.tsv', sep='\t', index=False)
    with np.errstate(invalid='ignore'):
        pairmat = pair_sum / pair_n
    pd.DataFrame(pairmat, index=regions, columns=regions).to_csv(f'{OUT}/m2b_neigh_pairwise_auroc_matrix.tsv', sep='\t')
    turn_rows = []
    for pi in range(P):
        js = []
        for (rA, rB) in itertools.combinations(regions, 2):
            a = topk[rA][pi]
            b = topk[rB][pi]
            j = len(a & b) / len(a | b) if a | b else np.nan
            js.append(j)
        turn_rows.append((progint[pi], float(np.nanmean(js)), 1.0 - float(np.nanmean(js))))
    pd.DataFrame(turn_rows, columns=['program', 'mean_partner_jaccard', 'partner_turnover']).to_csv(f'{OUT}/m2b_partner_turnover.tsv', sep='\t', index=False)
if __name__ == '__main__':
    main()
