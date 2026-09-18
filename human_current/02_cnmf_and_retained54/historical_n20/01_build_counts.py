#!/usr/bin/env python3
"""Historical preparation plus 20-run cNMF branch, not the complete 100-run history.
Input is full-cohort SCT-corrected snRNA_1M.h5ad. Apply V3 blacklist.
K=30..90.
"""
import os, sys, time, shutil, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blacklist_v3 import is_blacklisted

SRC = os.path.join(os.environ["SNRNA_WORK_DIR"], "snRNA_1M.h5ad")
OUT_BASE = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/cnmf_snrna_joint_full1M_v1")
LOG_DIR = OUT_BASE
os.makedirs(OUT_BASE, exist_ok=True)

NAME = 'snrna_joint_full1M_v1'
K_RANGE = [30, 40, 50, 60, 70, 80, 90]
N_ITER = 20
N_WORKERS = 48
N_TOP_GENES = 3000
SEED = 42
DENSITY_T = 0.15

LOG = f'{LOG_DIR}/{NAME}.log'
logf = open(LOG, 'w')
def log(m):
    s = f'[{time.strftime("%H:%M:%S")}] {m}'
    print(s, flush=True); logf.write(s+'\n'); logf.flush()

t0 = time.time()
log('=== Historical n20 branch on full-cohort SCT-corrected counts ===')
log(f'K={K_RANGE}, n_iter={N_ITER}, top_genes={N_TOP_GENES}, density_t={DENSITY_T}')

log('load SCT-corrected snRNA matrix')
a = ad.read_h5ad(SRC)
log(f'  shape: {a.shape}')

# Convert ENSG to symbol if needed
GENE_INFO = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/inputs/geneInfo_snRNA.csv")
gi = pd.read_csv(GENE_INFO)
ensg2sym = {r.gene_id: r.gene_name for r in gi.itertuples() if pd.notna(r.gene_name) and isinstance(r.gene_name, str)}
new_var = [ensg2sym.get(g, None) for g in a.var_names]
kept = []; seen = set(); new_names = []
for i, n in enumerate(new_var):
    if n is None or n in seen or n == 'nan': continue
    seen.add(n); kept.append(i); new_names.append(n)
a = a[:, kept].copy()
a.var_names = pd.Index(new_names)
a.var_names_make_unique()
log(f'  ENSG → symbol: {a.shape}')

# V3 blacklist
n_before = a.shape[1]
keep = np.array([not is_blacklisted(g, cell_group=None) for g in a.var_names])
a = a[:, keep].copy()
log(f'  V3 blacklist dropped {n_before-a.shape[1]}, remain {a.shape[1]}')

if not sp.issparse(a.X):
    a.X = sp.csr_matrix(a.X)
a.X = a.X.astype(np.float32)
sc.pp.filter_genes(a, min_cells=20)
sc.pp.filter_cells(a, min_counts=500)
log(f'  after filter: {a.shape}')

counts_h5 = f'{OUT_BASE}/{NAME}_counts.h5ad'
a.write_h5ad(counts_h5)
log(f'  wrote {counts_h5}')
if "--counts-only" in sys.argv:
    logf.close()
    sys.exit(0)

work_dir = f'{OUT_BASE}/cnmf_work'
if os.path.exists(work_dir):
    raise FileExistsError("Choose a new output directory; existing cNMF results are not deleted.")

from cnmf import cNMF
cnmf_obj = cNMF(output_dir=work_dir, name=NAME)
cnmf_obj.prepare(
    counts_fn=counts_h5,
    components=K_RANGE,
    n_iter=N_ITER,
    seed=SEED,
    num_highvar_genes=N_TOP_GENES,
    densify=False,
)
log('prepared, factorize ...')

import multiprocessing as mp
def factor_worker(wid):
    obj = cNMF(output_dir=work_dir, name=NAME)
    obj.factorize(worker_i=wid, total_workers=N_WORKERS)

t_fact = time.time()
procs = []
for w in range(N_WORKERS):
    p = mp.Process(target=factor_worker, args=(w,))
    p.start(); procs.append(p)
for p in procs:
    p.join()
log(f'factorize done | {time.time()-t_fact:.1f}s')

cnmf_obj.combine()
log('combine done')

import matplotlib; matplotlib.use('Agg')
try:
    cnmf_obj.k_selection_plot(close_fig=True)
except Exception as e:
    log(f'k_sel warn: {e}')

for k in [40, 50, 60, 70]:
    log(f'consensus K={k}')
    try:
        cnmf_obj.consensus(k=k, density_threshold=DENSITY_T,
                           show_clustering=True, close_clustergram_fig=True)
        gep_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.gene_spectra_score.k_{k}.dt_*.txt')
        usage_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.usages.k_{k}.dt_*.consensus.txt')
        if gep_cand:
            gep = pd.read_csv(gep_cand[0], sep='\t', index_col=0)
            gep.to_csv(f'{OUT_BASE}/{NAME}_k{k}_factor_loadings.tsv', sep='\t')
            log(f'  k={k} GEP exported {gep.shape}')
        if usage_cand:
            u = pd.read_csv(usage_cand[0], sep='\t', index_col=0)
            u.to_csv(f'{OUT_BASE}/{NAME}_k{k}_cell_scores.tsv', sep='\t')
    except Exception as e:
        log(f'  K={k} err: {e}')

with open(f'{OUT_BASE}/{NAME}_done.flag', 'w') as fh:
    fh.write(f'done at {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
log(f'\nDONE | {(time.time()-t0)/60:.1f} min')
logf.close()
