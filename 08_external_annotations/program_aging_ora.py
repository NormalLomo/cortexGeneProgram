#!/usr/bin/env python3
from __future__ import annotations
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
ROOT = PROJECT_ROOT
LOAD = ROOT / 'results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv'
NAMES = ROOT / 'results/crossregion_v1/program_names.tsv'
RENUMBER = ROOT / 'results/crossregion_v1/program_renumber_map.tsv'
GENESETS = ROOT / 'data/aging_gene_sets/aging_gene_sets_main_panel.gmt'
GENESET_META = ROOT / 'data/aging_gene_sets/aging_gene_set_metadata.tsv'
OUT = ROOT / 'results/crossregion_v1/program_aging'
NS = [100, 150, 200]
PRIMARY_N = 150
MIN_GENES_IN_BG = 10
FDR_SIG = 0.05

def read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            sets[name] = genes
    return sets

def enrich(top_set: set[str], term_set: set[str], n_bg: int) -> tuple[float, float, int]:
    k = len(term_set)
    n = len(top_set)
    x = len(top_set & term_set)
    p = hypergeom.sf(x - 1, n_bg, k, n) if x > 0 else 1.0
    a = x
    b = n - x
    c = k - x
    d = n_bg - n - k + x
    (odds, _) = fisher_exact([[a, b], [c, d]], alternative='greater')
    return (float(p), float(odds), int(x))

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    load = pd.read_csv(LOAD, sep='\t', index_col=0)
    load.index = [int(i) for i in load.index]
    renumber = pd.read_csv(RENUMBER, sep='\t')
    renumber = renumber[renumber['status'].astype(str).str.lower().eq('kept') & renumber['new_P'].notna()].copy()
    renumber['old_P'] = renumber['old_P'].astype(int)
    renumber['new_P'] = renumber['new_P'].astype(int)
    renumber = renumber.sort_values('new_P')
    if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
        raise ValueError()
    if renumber['new_P'].tolist() != list(range(1, 55)):
        raise ValueError()
    load = load.loc[renumber['old_P']].copy()
    load.index = renumber['new_P']
    genes = [str(g).upper() for g in load.columns]
    load.columns = genes
    gene_bg = set(genes)
    n_bg = len(gene_bg)
    if n_bg != 18742:
        raise RuntimeError()
    names = pd.read_csv(NAMES, sep='\t')
    names = names.merge(renumber[['old_P', 'new_P']], left_on='cnmf_component', right_on='old_P', suffixes=('_source', ''))
    name_short = dict(zip(names['new_P'], names['name_short']))
    name_full = dict(zip(names['new_P'], names['name_full']))
    confidence = dict(zip(names['new_P'], names['confidence']))
    retained = pd.read_csv(SCRIPT_DIR / 'inputs/program_display_map.tsv', sep='\t')
    retained['program'] = retained['program'].astype(int)
    retained['retained_program'] = retained['new_P'].astype(str).str.removeprefix('P').astype(int)
    retained_programs = set(range(1, 55))
    if dict(zip(retained['program'], retained['retained_program'])) != dict(zip(renumber['old_P'], renumber['new_P'])):
        raise ValueError()
    display_program = dict(zip(retained['retained_program'], retained['display_program'].astype(int)))
    new_p = {program: f'P{program}' for program in retained_programs}
    retained_class = dict(zip(retained['retained_program'], retained['dominant_class'].astype(str)))
    dominant_subclass = dict(zip(retained['retained_program'], retained['dominant_subclass'].astype(str)))
    raw_sets = read_gmt(GENESETS)
    aging_sets: dict[str, set[str]] = {}
    for (term, geneset) in raw_sets.items():
        in_bg = geneset & gene_bg
        if len(in_bg) >= MIN_GENES_IN_BG:
            aging_sets[term] = in_bg
    if not aging_sets:
        raise RuntimeError()
    meta = pd.read_csv(GENESET_META, sep='\t')
    meta_map = meta.set_index('set_name').to_dict('index')
    category_map = {'HAGR_GenAge_human': 'Aging curated', 'HAGR_CellAge_induces_senescence': 'Cellular senescence', 'HAGR_CellAge_inhibits_senescence': 'Cellular senescence', 'HAGR_CellAge_overexpressed_signature': 'Senescence expression', 'HAGR_CellAge_underexpressed_signature': 'Senescence expression', 'SAUL_SEN_MAYO': 'SASP/SenMayo', 'LIU2026_CELL_ACCELERATED_AGING_PROTEINS': 'Accelerated aging proteome', 'LIU2026_CELL_DECELERATED_AGING_PROTEINS': 'Decelerated aging proteome', 'LIU2026_CELL_PROTEOME_CLOCK_FEATURES': 'Aging clock proteome', 'LIU2026_CELL_ORGAN_PROTEIN_CLOCK_FEATURES': 'Aging clock proteome', 'LIU2026_CELL_PROXY_PROTEIN_CLOCK_FEATURES': 'Aging clock proteome'}
    aging_terms = [t for t in raw_sets.keys() if t in aging_sets]
    topn_genes: dict[int, dict[int, set[str]]] = {}
    for n in NS:
        topn_genes[n] = {}
        for program in load.index:
            top = load.loc[program].sort_values(ascending=False).head(n).index
            topn_genes[n][program] = set(top)
    records: dict[int, pd.DataFrame] = {}
    sig_counts: dict[int, int] = {}
    for n in NS:
        rows = []
        for program in load.index:
            top_set = topn_genes[n][program]
            for term in aging_terms:
                term_set = aging_sets[term]
                (pval, odds, overlap) = enrich(top_set, term_set, n_bg)
                overlap_genes = sorted(top_set & term_set)
                rows.append({'program': program, 'display_program': display_program.get(program, np.nan), 'new_P': new_p[program], 'name_short': name_short.get(program, str(program)), 'name_full': name_full.get(program, str(program)), 'confidence': confidence.get(program, ''), 'dominant_class': retained_class.get(program, ''), 'dominant_subclass': dominant_subclass.get(program, ''), 'aging_set': term, 'category': category_map.get(term, 'Aging'), 'n_aging_genes': len(term_set), 'topN': n, 'overlap': overlap, 'pval': pval, 'odds_ratio': odds, 'overlap_genes': ';'.join(overlap_genes)})
        df = pd.DataFrame(rows)
        df['fdr'] = multipletests(df['pval'].values, method='fdr_bh')[1]
        df['fdr_retained54'] = df['fdr']
        records[n] = df
        sig_counts[n] = int((df['fdr'] < FDR_SIG).sum())
        df.to_csv(OUT / f'enrichment_long_N{n}.tsv', sep='\t', index=False)
    prim = records[PRIMARY_N]
    sig = prim[prim['fdr'] < FDR_SIG].sort_values('fdr').copy()
    sig.to_csv(OUT / f'significant_pairs_N{PRIMARY_N}.tsv', sep='\t', index=False)
    top30 = sig.head(30)[['program', 'display_program', 'new_P', 'name_short', 'dominant_class', 'dominant_subclass', 'aging_set', 'category', 'overlap', 'n_aging_genes', 'odds_ratio', 'fdr', 'fdr_retained54', 'overlap_genes']].copy()
    top30.to_csv(OUT / f'top30_pairs_N{PRIMARY_N}.tsv', sep='\t', index=False)
    fdr_mat = prim.pivot(index='new_P', columns='aging_set', values='fdr')
    or_mat = prim.pivot(index='new_P', columns='aging_set', values='odds_ratio')
    neglog = -np.log10(fdr_mat.clip(lower=1e-300))
    fdr_mat.to_csv(OUT / f'fdr_matrix_N{PRIMARY_N}.tsv', sep='\t')
    or_mat.to_csv(OUT / f'or_matrix_N{PRIMARY_N}.tsv', sep='\t')
    neglog.to_csv(OUT / f'neglog10fdr_matrix_N{PRIMARY_N}.tsv', sep='\t')
    sensitivity_rows = []
    for n in NS:
        df = records[n]
        selected = df[df['fdr'] < FDR_SIG]
        sensitivity_rows.append({'topN': n, 'n_sig_pairs': int(len(selected)), 'n_programs_hit': int(selected['new_P'].nunique()), 'n_aging_sets_hit': int(selected['aging_set'].nunique())})
    pd.DataFrame(sensitivity_rows).to_csv(OUT / 'n_sensitivity.tsv', sep='\t', index=False)
    category_counts = sig.groupby('category')['new_P'].nunique().sort_values(ascending=False).rename('n_programs').reset_index()
    category_counts.to_csv(OUT / f'category_program_counts_N{PRIMARY_N}.tsv', sep='\t', index=False)
    program_summary = sig.assign(neglog10_fdr=-np.log10(sig['fdr'].clip(lower=1e-300))).groupby(['new_P', 'program', 'name_short', 'dominant_class', 'dominant_subclass'], as_index=False).agg(n_aging_sets_hit=('aging_set', 'nunique'), best_aging_set=('aging_set', lambda s: s.iloc[0]), best_fdr=('fdr', 'min'), best_odds_ratio=('odds_ratio', 'max'), total_overlap=('overlap', 'sum')).sort_values(['n_aging_sets_hit', 'best_fdr'], ascending=[False, True])
    program_summary.to_csv(OUT / f'program_aging_summary_N{PRIMARY_N}.tsv', sep='\t', index=False)
    meta_out = {'analysis': 'program_aging_ora', 'n_background': n_bg, 'loadings': str(LOAD), 'retained_program_map': str(RENUMBER), 'programs_n': len(retained_programs), 'gene_sets': {term: {'n_genes_in_background': len(aging_sets[term]), 'category': category_map.get(term, 'Aging'), 'source': meta_map.get(term, {}).get('source', ''), 'description': meta_map.get(term, {}).get('description', '')} for term in aging_terms}, 'NS': NS, 'primary_N': PRIMARY_N, 'fdr_sig': FDR_SIG, 'n_sig_by_N': sig_counts}
    (OUT / 'meta.json').write_text(json.dumps(meta_out, indent=2))
if __name__ == '__main__':
    main()
