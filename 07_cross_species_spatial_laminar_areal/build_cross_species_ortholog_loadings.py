#!/usr/bin/env python
import os, json, time
import numpy as np, pandas as pd
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
t0 = time.time()
PROJ = os.environ.get('CORTEX_PROGRAM_ROOT', '/home/luomeng/DATA/cortex_nmf_program')
ORTHO = '/home/luomeng/DATA/pseudo_human_brain_spatial_reference_jiaxi/data/ortholog/biomart_mouse_to_human.tsv'
SPECTRA = f'{PROJ}/results/cnmf_snrna_joint_full1M_v1/cnmf_work/snrna_joint_full1M_v1/snrna_joint_full1M_v1.gene_spectra_tpm.k_60.dt_0_15.txt'
RETAIN_MAP = f'{PROJ}/results/crossregion_v1/program_renumber_map.tsv'
OUT = f'{PROJ}/results/crossregion_v1/crossspecies'
os.makedirs(OUT, exist_ok=True)
spec = pd.read_csv(SPECTRA, sep='\t', index_col=0)
spec.index = [str(i) for i in spec.index]
program_map = pd.read_csv(RETAIN_MAP, sep='\t')
program_map = program_map.loc[program_map['status'].astype(str).str.lower().eq('kept') & program_map['new_P'].notna()].copy()
program_map['old_id'] = program_map['old_P'].astype(int).astype(str)
program_map['new_order'] = program_map['new_P'].astype(str).str.removeprefix('P').astype(int)
program_map['new_id'] = 'P' + program_map['new_order'].astype(str)
program_map = program_map.sort_values('new_order')
spec = spec.loc[program_map['old_id']].copy()
spec.index = program_map['new_id'].tolist()
human_genes = list(spec.columns)
orth = pd.read_csv(ORTHO, sep='\t')
orth = orth[orth.orthology_type == 'ortholog_one2one'].copy()
orth = orth.dropna(subset=['human_symbol', 'mouse_symbol'])
hs = orth.human_symbol.value_counts()
ms = orth.mouse_symbol.value_counts()
orth = orth[orth.human_symbol.isin(hs[hs == 1].index) & orth.mouse_symbol.isin(ms[ms == 1].index)]
h2m = dict(zip(orth.human_symbol, orth.mouse_symbol))
retained = [g for g in human_genes if g in h2m]
pct = 100 * len(retained) / len(human_genes)
cov = []
for p in spec.index:
    top = spec.loc[p].sort_values(ascending=False).head(50).index
    cov.append((p, float(np.mean([g in h2m for g in top]))))
cov_df = pd.DataFrame(cov, columns=['program', 'top50_ortholog_frac']).set_index('program')
cov_df.to_csv(f'{OUT}/ortholog_coverage_per_program.tsv', sep='\t')
spec_sub = spec[retained].copy()
spec_sub.columns = [h2m[g] for g in retained]
mouse_load = spec_sub.T.copy()
mouse_load.index.name = 'mouse_symbol'
mouse_load.columns = program_map['new_id'].tolist()
mouse_load.to_parquet(f'{OUT}/mouse_loadings.parquet')
summary = {'n_human_spectra_genes': len(human_genes), 'n_clean_1to1_pairs': len(h2m), 'n_retained_ortholog_genes': len(retained), 'pct_retained': round(pct, 2), 'per_program_top50_coverage_median': round(float(cov_df.top50_ortholog_frac.median()), 3), 'per_program_top50_coverage_min': round(float(cov_df.top50_ortholog_frac.min()), 3), 'loadings_nonneg_min': round(float(mouse_load.values.min()), 6), 'runtime_sec': round(time.time() - t0, 1)}
with open(f'{OUT}/ortholog_build_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
