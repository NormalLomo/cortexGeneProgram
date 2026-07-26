#!/usr/bin/env python
import os, json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
ROOT = str(PROJECT_ROOT)
LOAD = f'{ROOT}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_factor_loadings.tsv'
NAMES = f'{ROOT}/results/crossregion_v1/program_names.tsv'
RENUM = f'{ROOT}/results/crossregion_v1/program_renumber_map.tsv'
BB = f'{ROOT}/data/brainbase/disease_gene_associations.txt'
OUT = f'{ROOT}/results/crossregion_v1/program_disease'
os.makedirs(OUT, exist_ok=True)
NS = [100, 150, 200]
PRIMARY_N = 150
MIN_GENES = 10
FDR_SIG = 0.05
load = pd.read_csv(LOAD, sep='\t', index_col=0)
load.index = [int(i) for i in load.index]
renumber = pd.read_csv(RENUM, sep='\t')
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
GENES = list(load.columns)
gene_set_bg = set(GENES)
N_BG = len(GENES)
names = pd.read_csv(NAMES, sep='\t')
source_col = 'cnmf_component' if 'cnmf_component' in names.columns else 'program'
names = names.merge(renumber[['old_P', 'new_P']], left_on=source_col, right_on='old_P', suffixes=('_source', ''))
name_short = dict(zip(names['new_P'], names['name_short']))
name_full = dict(zip(names['new_P'], names['name_full']))
conf = dict(zip(names['new_P'], names.get('confidence', pd.Series(dtype=str))))
bb = pd.read_csv(BB, sep='\t')
bb.columns = [c.strip() for c in bb.columns]
gene_col = 'Gene symbol' if 'Gene symbol' in bb.columns else 'Gene'
mrna = bb[bb['Type'] == 'mRNA'].copy()
mrna[gene_col] = mrna[gene_col].astype(str).str.strip()
disease_sets = {}
for (dis, sub) in mrna.groupby('Disease'):
    g = set(sub[gene_col]) & gene_set_bg
    if len(g) >= MIN_GENES:
        disease_sets[dis] = g
diseases = sorted(disease_sets, key=lambda d: -len(disease_sets[d]))
CATEGORY = {"Alzheimer's Disease": 'Neurodegenerative', "Parkinson's Disease": 'Neurodegenerative', 'Amyotrophic Lateral Sclerosis': 'Neurodegenerative', "Huntington's Disease": 'Neurodegenerative', 'Multiple Sclerosis': 'Neurodegenerative', 'Motor Neurone Disease': 'Neurodegenerative', 'Charcot-Marie-Tooth Disease': 'Neurodegenerative', 'Autism Spectrum Disorder': 'Psychiatric/Neurodev', 'Schizophrenia': 'Psychiatric/Neurodev', 'Bipolar Disorder': 'Psychiatric/Neurodev', 'Major Depressive Disorder': 'Psychiatric/Neurodev', 'Depression': 'Psychiatric/Neurodev', 'Intellectual Disability': 'Psychiatric/Neurodev', 'Attention Deficit Hyperactivity Disorder': 'Psychiatric/Neurodev', 'Rett Syndrome': 'Psychiatric/Neurodev', 'Narcolepsy': 'Psychiatric/Neurodev', 'Restless legs Syndrome': 'Psychiatric/Neurodev', 'Glioma': 'Tumor', 'Neuroblastoma': 'Tumor', 'Glioblastoma': 'Tumor', 'Meningioma': 'Tumor', 'Medulloblastoma': 'Tumor', 'Down Syndrome': 'Developmental', 'Cerebral Palsy': 'Developmental', 'Hydrocephalus': 'Developmental', 'Muscular Dystrophy': 'Developmental', 'Tuberous Sclerosis': 'Developmental', 'Agenesis Corpus Callosum': 'Developmental', 'Prader-Willi Syndrome': 'Developmental', 'Epilepsy': 'Other/Vascular', 'Stroke': 'Other/Vascular', 'Intracranial Aneurysm': 'Other/Vascular', 'Cerebral Aneurysm': 'Other/Vascular', 'Dystonia': 'Other/Vascular', 'Gaucher Disease': 'Other/Vascular', 'Dysautonomia': 'Other/Vascular', 'Coma': 'Other/Vascular', 'Encephalitis': 'Other/Vascular'}
cat_of = {d: CATEGORY.get(d, 'Other/Vascular') for d in diseases}
topN_genes = {}
for N in NS:
    topN_genes[N] = {}
    for p in load.index:
        row = load.loc[p]
        top = row.sort_values(ascending=False).head(N).index
        topN_genes[N][p] = set(top)

def enrich(top_set, dis_set):
    K = len(dis_set)
    n = len(top_set)
    x = len(top_set & dis_set)
    p = hypergeom.sf(x - 1, N_BG, K, n) if x > 0 else 1.0
    a = x
    b = n - x
    c = K - x
    d = N_BG - n - K + x
    (orr, _) = fisher_exact([[a, b], [c, d]], alternative='greater')
    return (p, orr, x)
records = {}
nsig_by_N = {}
for N in NS:
    rows = []
    for p in load.index:
        ts = topN_genes[N][p]
        for dis in diseases:
            ds = disease_sets[dis]
            (pv, orr, ov) = enrich(ts, ds)
            ovg = sorted(ts & ds)
            rows.append(dict(program=p, name_short=name_short.get(p, str(p)), disease=dis, category=cat_of[dis], n_disease_genes=len(ds), topN=N, overlap=ov, pval=pv, odds_ratio=orr, overlap_genes=';'.join(ovg)))
    df = pd.DataFrame(rows)
    df['fdr'] = multipletests(df['pval'].values, method='fdr_bh')[1]
    records[N] = df
    nsig = int((df['fdr'] < FDR_SIG).sum())
    nsig_by_N[N] = nsig
    df.to_csv(f'{OUT}/enrichment_long_N{N}.tsv', sep='\t', index=False)
prim = records[PRIMARY_N]
fdr_mat = prim.pivot(index='program', columns='disease', values='fdr')
or_mat = prim.pivot(index='program', columns='disease', values='odds_ratio')
ov_mat = prim.pivot(index='program', columns='disease', values='overlap')
neglog = -np.log10(fdr_mat.clip(lower=1e-300))
col_order = sorted(diseases, key=lambda d: (cat_of[d], -len(disease_sets[d])))
fdr_mat = fdr_mat[col_order]
or_mat = or_mat[col_order]
ov_mat = ov_mat[col_order]
neglog = neglog[col_order]
fdr_mat.to_csv(f'{OUT}/fdr_matrix_N{PRIMARY_N}.tsv', sep='\t')
or_mat.to_csv(f'{OUT}/or_matrix_N{PRIMARY_N}.tsv', sep='\t')
neglog.to_csv(f'{OUT}/neglog10fdr_matrix_N{PRIMARY_N}.tsv', sep='\t')
sig = prim[prim['fdr'] < FDR_SIG].sort_values('fdr').copy()
sig.to_csv(f'{OUT}/significant_pairs_N{PRIMARY_N}.tsv', sep='\t', index=False)
sens = pd.DataFrame({'topN': NS, 'n_sig_pairs': [nsig_by_N[N] for N in NS]})
extra = []
for N in NS:
    s = records[N][records[N]['fdr'] < FDR_SIG]
    extra.append(dict(topN=N, n_programs_hit=s['program'].nunique(), n_diseases_hit=s['disease'].nunique()))
sens = sens.merge(pd.DataFrame(extra), on='topN')
sens.to_csv(f'{OUT}/n_sensitivity.tsv', sep='\t', index=False)
catsum = sig.groupby('category')['program'].nunique().reindex(['Neurodegenerative', 'Psychiatric/Neurodev', 'Tumor', 'Developmental', 'Other/Vascular']).fillna(0).astype(int)
catsum.to_csv(f'{OUT}/category_program_counts_N{PRIMARY_N}.tsv', sep='\t', header=['n_programs'])
meta = dict(n_background=N_BG, n_diseases=len(diseases), diseases=diseases, disease_n_genes={d: len(disease_sets[d]) for d in diseases}, category=cat_of, NS=NS, primary_N=PRIMARY_N, n_sig_by_N=nsig_by_N, fdr_sig=FDR_SIG)
with open(f'{OUT}/meta.json', 'w') as f:
    json.dump(meta, f, indent=2)
top_pairs = sig.head(30)[['program', 'name_short', 'disease', 'category', 'overlap', 'n_disease_genes', 'odds_ratio', 'fdr', 'overlap_genes']].copy()
top_pairs.to_csv(f'{OUT}/top30_pairs_N{PRIMARY_N}.tsv', sep='\t', index=False)
