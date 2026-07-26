#!/usr/bin/env python
import os, sys, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata
BASE = os.environ.get('CORTEX_NMF_ROOT', os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
OUT = os.environ.get('XREGION_OUTPUT_DIR', f'{BASE}/results/crossregion_v1/p37_current54')
SRC = f'{BASE}/results/crossregion_v1/region_subclass_program_mean.tsv'
RETAIN_MAP = os.environ.get('RETAIN_MAP', f'{BASE}/results/crossregion_v1/program_renumber_map.tsv')
os.makedirs(OUT, exist_ok=True)

def zshape(v):
    v = np.asarray(v, float)
    sd = v.std()
    if sd == 0 or not np.isfinite(sd):
        return v - v.mean()
    return (v - v.mean()) / sd

def spearman_vec(a, B):
    ar = rankdata(a)
    ar = ar - ar.mean()
    out = np.empty(B.shape[1])
    for j in range(B.shape[1]):
        br = rankdata(B[:, j])
        br = br - br.mean()
        denom = np.sqrt((ar ** 2).sum() * (br ** 2).sum())
        out[j] = (ar * br).sum() / denom if denom > 0 else 0.0
    return out

def main():
    mapping = pd.read_csv(RETAIN_MAP, sep='\t')
    mapping = mapping[mapping['status'].astype(str).str.lower().eq('kept') & mapping['new_P'].notna()].copy()
    mapping['old_int'] = mapping['old_P'].astype(int)
    mapping['new_int'] = mapping['new_P'].astype(str).str.removeprefix('P').astype(int)
    old_to_new = dict(zip(mapping['old_int'], mapping['new_int']))
    df = pd.read_csv(SRC, sep='\t')
    df['program'] = df['program'].astype(int)
    df = df[df['program'].isin(mapping['old_int'])].copy()
    regions = sorted(df['region'].unique())
    programs = mapping['old_int'].tolist()
    display_programs = mapping['new_int'].tolist()
    subclasses = sorted(df['subclass'].unique())
    P = len(programs)
    S = len(subclasses)
    R = len(regions)
    pidx = {p: i for (i, p) in enumerate(programs)}
    sidx = {s: i for (i, s) in enumerate(subclasses)}
    Mraw = {r: np.zeros((S, P)) for r in regions}
    for (r, sub, prog, m) in df[['region', 'subclass', 'program', 'mean']].itertuples(index=False):
        Mraw[r][sidx[sub], pidx[prog]] = m
    M = {}
    for r in regions:
        Z = np.apply_along_axis(zshape, 0, Mraw[r])
        M[r] = Z
    self_rows = []
    pair_self_sum = np.zeros((R, R))
    pair_self_n = np.zeros((R, R))
    ridx = {r: i for (i, r) in enumerate(regions)}
    for (rA, rB) in itertools.permutations(regions, 2):
        A = M[rA]
        Bm = M[rB]
        for pi in range(P):
            sims = spearman_vec(A[:, pi], Bm)
            s_self = sims[pi]
            others = np.delete(sims, pi)
            gt = np.sum(s_self > others)
            eq = np.sum(s_self == others)
            auroc = (gt + 0.5 * eq) / (P - 1)
            self_rows.append((display_programs[pi], rA, rB, auroc))
            pair_self_sum[ridx[rA], ridx[rB]] += auroc
            pair_self_n[ridx[rA], ridx[rB]] += 1
    selfdf = pd.DataFrame(self_rows, columns=['program', 'regionA', 'regionB', 'self_auroc'])
    selfdf.to_csv(f'{OUT}/m1_expr_program_region_self_auroc.tsv', sep='\t', index=False)
    perprog = selfdf.groupby('program')['self_auroc'].agg(['mean', 'count']).reset_index()
    perprog.columns = ['program', 'expr_cons_auroc', 'n_pairs']
    perprog = perprog.sort_values('expr_cons_auroc', ascending=False)
    perprog.to_csv(f'{OUT}/m1_expr_conservation_per_program.tsv', sep='\t', index=False)
    with np.errstate(invalid='ignore'):
        pairmat = pair_self_sum / pair_self_n
    pm = pd.DataFrame(pairmat, index=regions, columns=regions)
    pm.to_csv(f'{OUT}/m1_expr_pairwise_auroc_matrix.tsv', sep='\t')
if __name__ == '__main__':
    main()
