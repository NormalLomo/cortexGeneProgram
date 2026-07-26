#!/usr/bin/env python
def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import os, json
import numpy as np
import pandas as pd
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BASE = '/home/luomeng/DATA/cortex_nmf_program'
OUT = f'{BASE}/results/xspecies_humanmap_v1'
NAMES = f'{BASE}/results/crossregion_v1/program_names.tsv'
ANNOT = f'{BASE}/results/crossregion_v1/program_annotation_gobp.tsv'
N_PERM = 2000
RNG = np.random.default_rng(12345)
SIG_FDR = 0.05

def cos_row(a, B):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(B, axis=1)
    out = B @ a
    denom = nb * na
    with np.errstate(invalid='ignore', divide='ignore'):
        out = out / denom
    out[~np.isfinite(out)] = np.nan
    return out

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty(n)
    out[order] = fdr
    return out

def species_null(npz_path, species):
    z = np.load(npz_path, allow_pickle=True)
    H = z['Hload'].astype(np.float64)
    R = z['refit'].astype(np.float64)
    programs = [str(x) for x in z['programs']]
    (n_prog, G) = H.shape
    true_cos = np.array([float(R[i] @ H[i] / (np.linalg.norm(R[i]) * np.linalg.norm(H[i]) + 1e-300)) for i in range(n_prog)])
    p_a = np.empty(n_prog)
    perc_a = np.empty(n_prog)
    diag_rank = np.empty(n_prog, dtype=int)
    for i in range(n_prog):
        allc = cos_row(R[i], H)
        null = np.delete(allc, i)
        null = null[np.isfinite(null)]
        t = true_cos[i]
        p_a[i] = (np.sum(null >= t) + 1) / (len(null) + 1)
        perc_a[i] = 100.0 * np.mean(null < t)
        diag_rank[i] = int(np.sum(allc >= allc[i]))
    p_b = np.empty(n_prog)
    for i in range(n_prog):
        hi = H[i]
        nh = np.linalg.norm(hi)
        ri = R[i]
        nr = np.linalg.norm(ri)
        t = true_cos[i]
        cnt = 0
        for _ in range(N_PERM):
            perm = RNG.permutation(G)
            c = ri[perm] @ hi / (nr * nh + 1e-300)
            if c >= t:
                cnt += 1
        p_b[i] = (cnt + 1) / (N_PERM + 1)
    fdr_a = bh_fdr(p_a)
    fdr_b = bh_fdr(p_b)
    identity_gate = diag_rank <= max(1, n_prog // 3)
    sig_primary = ((fdr_b < SIG_FDR) & identity_gate).astype(int)
    df = pd.DataFrame({'program': programs, f'{species}_cosine': true_cos, f'{species}_p_perm': p_b, f'{species}_fdr_perm': fdr_b, f'{species}_diag_rank': diag_rank, f'{species}_percentile_xprog': perc_a, f'{species}_p_xprog': p_a, f'{species}_fdr_xprog': fdr_a, f'{species}_identity_gate': identity_gate.astype(int), f'{species}_sig': sig_primary})
    return (df, true_cos, programs)
NEURON_KEYS = ['synap', 'axon', 'neuron', 'neuropeptide', 'glutamate', 'gaba', 'cholinergic', 'nerve impulse', 'neural crest', 'neurotransmit', 'retrograde axonal', 'myelinat']
GLIA_VASC_CYTO_IMMUNE_KEYS = ['lipoprotein', 'sterol', 'extracellular matrix', 'external encapsulating', 'vasoconstrict', 'blood vessel', 'angiogen', 'endothelial', 'cytoskelet', 'cilium', 'muscle contraction', 'actin', 'b cell', 't cell', 'leukocyte', 'lymphocyte', 'immune', 'antimicrobial', 'osteoclast', 'humoral', 'splic', 'iron-sulfur', 'proton transmembrane', 'secretion', 'kinase', 'mrna', 'mirna', 'transcript']
BRAIN_RELATED_KEYS = NEURON_KEYS + ['secretion', 'calcium ion transmembrane', 'anion transmembrane', 'amino acid transmembrane', 'vasoconstrict', 'blood vessel', 'angiogen', 'endothelial', 'myelinat', 'growth', 'chemotaxis', 'zinc ion']
NON_BRAIN_KEYS = ['kidney', 'eye', 'ear', 'taste', 'skeletal', 'cardiac', 'muscle cell fate', 'striated muscle', 'endoderm', 'embryonic skeletal', 'sensory organ', 'sensory perception of taste', 'osteoclast', 'epithelial', 'retinoic acid', 'nodal', 'neural crest', 'morphogenesis of', 'substantia nigra']

def classify(func_lower):
    big = 'neuron_synapse_axon' if any((k in func_lower for k in NEURON_KEYS)) else 'glia_vascular_cytoskeleton_immune'
    return big

def brain_relevance(func_lower):
    if any((k in func_lower for k in NEURON_KEYS)):
        return 'brain'
    if any((k in func_lower for k in NON_BRAIN_KEYS)) and (not any((k in func_lower for k in BRAIN_RELATED_KEYS))):
        return 'non_brain_generic'
    return 'other'

def main():
    (df_mac, tc_mac, progs) = species_null(f'{OUT}/loadings_macaque.npz', 'h_mac')
    (df_mou, tc_mou, _) = species_null(f'{OUT}/loadings_mouse.npz', 'h_mou')
    names = _annotation_table(pd.read_csv(NAMES, sep='\t'))
    names['program'] = ['P' + str(int(p)) for p in names['program']]
    func_map = dict(zip(names['program'], names['name_short']))
    conf_map = dict(zip(names['program'], names['confidence']))
    CLASS_B = {p for (p, c) in conf_map.items() if str(c) == 'class_b'}
    ann = _annotation_table(pd.read_csv(ANNOT, sep='\t'))
    ann['program'] = ['P' + str(int(p)) for p in ann['program']]
    top3_map = dict(zip(ann['program'], ann.get('top3_BP_terms', ann['program'])))
    m = df_mac.merge(df_mou, on='program', how='inner')
    m['func_name'] = m['program'].map(func_map)
    m['top3_BP'] = m['program'].map(top3_map)
    sig = pd.DataFrame({'program': m['program'], 'func_name': m['func_name'], 'confidence': m['program'].map(conf_map), 'h_mac_cosine': m['h_mac_cosine'], 'h_mac_p': m['h_mac_p_perm'], 'h_mac_fdr': m['h_mac_fdr_perm'], 'h_mac_sig': m['h_mac_sig'].astype(bool), 'h_mou_cosine': m['h_mou_cosine'], 'h_mou_p': m['h_mou_p_perm'], 'h_mou_fdr': m['h_mou_fdr_perm'], 'h_mou_sig': m['h_mou_sig'].astype(bool), 'h_mac_diag_rank': m['h_mac_diag_rank'], 'h_mac_identity_gate': m['h_mac_identity_gate'], 'h_mac_p_xprog': m['h_mac_p_xprog'], 'h_mac_fdr_xprog': m['h_mac_fdr_xprog'], 'h_mac_percentile_xprog': m['h_mac_percentile_xprog'], 'h_mou_diag_rank': m['h_mou_diag_rank'], 'h_mou_identity_gate': m['h_mou_identity_gate'], 'h_mou_p_xprog': m['h_mou_p_xprog'], 'h_mou_fdr_xprog': m['h_mou_fdr_xprog'], 'h_mou_percentile_xprog': m['h_mou_percentile_xprog']})
    sig['_n'] = sig['program'].str.lstrip('P').astype(int)
    sig = sig.sort_values('_n').drop(columns='_n').reset_index(drop=True)
    sig.to_csv(f'{OUT}/conservation_significance_per_program.csv', index=False)
    rows = []
    for (_, r) in sig.iterrows():
        fn = str(r['func_name']).lower()
        big = classify(fn)
        brain = brain_relevance(fn)
        mac_c = r['h_mac_cosine']
        mac_sig = bool(r['h_mac_sig'])
        mac_rank1 = int(r['h_mac_diag_rank']) == 1
        tier = 'class_b' if r['program'] in CLASS_B else 'class_a'
        rows.append({'program': r['program'], 'func_name': r['func_name'], 'confidence': str(conf_map.get(r['program'], '')), 'func_class': big, 'brain_relevance': brain, 'h_mac_cosine': mac_c, 'h_mac_sig': mac_sig, 'h_mac_diag_rank': int(r['h_mac_diag_rank']), 'h_mou_cosine': r['h_mou_cosine'], 'h_mou_sig': bool(r['h_mou_sig']), 'conservation_tier': tier})
    strat = pd.DataFrame(rows)
    strat.to_csv(f'{OUT}/conservation_function_stratification.csv', index=False)
    n_mac_sig = int(sig['h_mac_sig'].sum())
    n_mou_sig = int(sig['h_mou_sig'].sum())
    n_mac_rank1 = int((sig['h_mac_diag_rank'] == 1).sum())
    n_mou_rank1 = int((sig['h_mou_diag_rank'] == 1).sum())
    n_mac_gate = int(sig['h_mac_identity_gate'].sum())
    n_mou_gate = int(sig['h_mou_identity_gate'].sum())
    glia = strat[strat['func_class'] == 'glia_vascular_cytoskeleton_immune']
    neur = strat[strat['func_class'] == 'neuron_synapse_axon']
    bsig = strat[strat['conservation_tier'] == 'class_a']['h_mac_cosine']
    bweak = strat[strat['conservation_tier'] == 'class_b']['h_mac_cosine']
    pmw = float('nan')
    try:
        from scipy.stats import mannwhitneyu
        bs = bsig.dropna()
        bw = bweak.dropna()
        if len(bs) and len(bw):
            (u, pmw) = mannwhitneyu(bs, bw, alternative='greater')
    except Exception as ex:
        raise
    s = sig.sort_values('h_mac_cosine', ascending=False)
    summary = {'null_primary': 'gene-label permutation N=2000 (BH-FDR)', 'n_mac_sig': n_mac_sig, 'n_mou_sig': n_mou_sig, 'n_mac_rank1_xprog': n_mac_rank1, 'n_mou_rank1_xprog': n_mou_rank1, 'n_mac_gate_pass': n_mac_gate, 'n_mou_gate_pass': n_mou_gate, 'glia_mac_median': float(np.nanmedian(glia['h_mac_cosine'])), 'neuron_mac_median': float(np.nanmedian(neur['h_mac_cosine'])), 'class_a_mac_median': float(np.nanmedian(bsig)), 'class_b_mac_median': float(np.nanmedian(bweak)), 'n_class_a': int(len(bsig)), 'n_class_b': int(len(bweak)), 'class_a_vs_b_mwu_p': float(pmw)}
    with open(f'{OUT}/conservation_null_summary.json', 'w') as fo:
        json.dump(summary, fo, indent=2)
if __name__ == '__main__':
    main()
