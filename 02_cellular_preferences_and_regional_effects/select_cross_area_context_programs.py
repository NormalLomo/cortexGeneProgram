#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
BASE = os.environ.get('CORTEX_NMF_ROOT', os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
RES = f'{BASE}/results/crossregion_v1'
OUT = os.environ.get('XREGION_OUTPUT_DIR', f'{BASE}/results/crossregion_v1/p37_current54')
RETAIN_MAP = os.environ.get('RETAIN_MAP', f'{BASE}/results/crossregion_v1/program_renumber_map.tsv')
os.makedirs(OUT, exist_ok=True)

def main():
    mapping = pd.read_csv(RETAIN_MAP, sep='\t')
    mapping = mapping[mapping['status'].astype(str).str.lower().eq('kept') & mapping['new_P'].notna()].copy()
    mapping['old_int'] = mapping['old_P'].astype(int)
    mapping['new_int'] = mapping['new_P'].astype(str).str.removeprefix('P').astype(int)
    old_to_new = dict(zip(mapping['old_int'], mapping['new_int']))
    m1 = pd.read_csv(f'{OUT}/m1_expr_conservation_per_program.tsv', sep='\t')
    m2b = pd.read_csv(f'{OUT}/m2b_neighborhood_conservation_per_program.tsv', sep='\t')
    m2a = pd.read_csv(f'{OUT}/m2a_spatial_neighborhood_conservation_per_program.tsv', sep='\t')
    turn = pd.read_csv(f'{OUT}/m2b_partner_turnover.tsv', sep='\t')
    var = pd.read_csv(f'{RES}/program_variability.tsv', sep='\t')[['program', 'class', 'eta2_region', 'fdr']]
    var['program'] = var['program'].map(old_to_new)
    var = var.dropna(subset=['program'])
    var['program'] = var['program'].astype(int)
    names = pd.read_csv(f'{RES}/program_names.tsv', sep='\t')[['cnmf_component', 'name_short', 'name_full', 'confidence', 'fdr']]
    names['program'] = names['cnmf_component'].map(old_to_new)
    names = names.dropna(subset=['program'])
    names['program'] = names['program'].astype(int)
    names = names.drop(columns=['cnmf_component'])
    names = names.rename(columns={'fdr': 'name_fdr'})
    sig = pd.read_csv(f'{RES}/supp_table_region_signatures.tsv', sep='\t')
    M = m1.merge(m2b, on='program', suffixes=('', '_b')).merge(turn, on='program')
    M = M.merge(m2a[['program', 'spatial_neigh_cons_auroc']], on='program', how='left')
    M = M.merge(var, on='program', how='left').merge(names, on='program', how='left')
    sigset = {old_to_new[p] for p in sig['program'].unique() if p in old_to_new}
    M['in_region_signature'] = M['program'].isin(sigset)
    A_thr = M['expr_cons_auroc'].quantile(2 / 3)
    B_thr = M['neigh_cons_auroc'].quantile(1 / 3)
    M['gateA_expr_conserved'] = M['expr_cons_auroc'] >= max(A_thr, 0.9)
    M['gateB_neigh_low'] = M['neigh_cons_auroc'] <= B_thr
    M['rank_neigh'] = M['neigh_cons_auroc'].rank()
    M['rank_spatial'] = M['spatial_neigh_cons_auroc'].rank()
    turn_thr = M['partner_turnover'].quantile(2 / 3)
    n_programs = len(M)
    neigh_rank_cut = int(np.ceil(0.5 * n_programs))
    spatial_rank_cut = int(np.ceil(0.5833333333333334 * n_programs))
    M['gateC_coact_spatial_agree'] = (M['rank_neigh'] <= neigh_rank_cut) & (M['rank_spatial'] <= spatial_rank_cut) & (M['partner_turnover'] >= turn_thr)
    M['gateD_anchor'] = (M['class'] == 'variable') | M['in_region_signature']
    M['rewiring_gap'] = M['expr_cons_auroc'] - M['neigh_cons_auroc']
    M['n_gates_pass'] = M[['gateA_expr_conserved', 'gateB_neigh_low', 'gateC_coact_spatial_agree', 'gateD_anchor']].sum(1)
    M = M.sort_values(['n_gates_pass', 'rewiring_gap'], ascending=False)
    M.to_csv(f'{OUT}/m4_program_master_table.tsv', sep='\t', index=False)
    hits = M[M['gateA_expr_conserved'] & M['gateB_neigh_low'] & M['gateC_coact_spatial_agree'] & M['gateD_anchor']].copy()
    neighlong = pd.read_csv(f'{OUT}/m2b_neigh_program_region_self_auroc.tsv', sep='\t')
    tgt_rows = []
    for p in hits['program']:
        sub = neighlong[neighlong['program'] == p].groupby('regionB')['self_auroc'].mean().sort_values()
        tgt_rows.append((p, ';'.join(sub.index[:3]), round(float(sub.iloc[0]), 3)))
    tgtdf = pd.DataFrame(tgt_rows, columns=['program', 'rewiring_target_regions', 'min_region_auroc'])
    hits = hits.merge(tgtdf, on='program', how='left')
    cols = ['program', 'name_short', 'name_full', 'confidence', 'name_fdr', 'class', 'eta2_region', 'expr_cons_auroc', 'neigh_cons_auroc', 'spatial_neigh_cons_auroc', 'rewiring_gap', 'partner_turnover', 'mean_partner_jaccard', 'in_region_signature', 'rewiring_target_regions', 'min_region_auroc', 'n_gates_pass']
    cols = [c for c in cols if c in hits.columns]
    hits[cols].to_csv(f'{OUT}/m4_rewiring_hits.tsv', sep='\t', index=False)
    fn = pd.DataFrame({'stage': ['all_programs', 'gateA_expr_conserved', '+gateB_neigh_low', '+gateC_coact_spatial', '+gateD_anchor_HITS'], 'n': [len(M), int(M['gateA_expr_conserved'].sum()), int((M['gateA_expr_conserved'] & M['gateB_neigh_low']).sum()), int((M['gateA_expr_conserved'] & M['gateB_neigh_low'] & M['gateC_coact_spatial_agree']).sum()), len(hits)]})
    fn.to_csv(f'{OUT}/m4_filter_funnel.tsv', sep='\t', index=False)
if __name__ == '__main__':
    main()
