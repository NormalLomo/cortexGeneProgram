#!/usr/bin/env python3
"""Historical n20 K=65 add-on; this is not a recovered 80-run batch or 80+20 fusion."""
import os, time, glob
import pandas as pd

OUT_BASE = os.environ["CORTEX_PROGRAM_ROOT"] + "/results/cnmf_snrna_joint_full1M_v1"
COUNTS = os.environ["CNMF_COUNTS_H5AD"]
NAME = 'snrna_joint_full1M_v1_k65'
K_RANGE = [65]
N_ITER = 20
N_WORKERS = 60
N_TOP_GENES = 3000
SEED = 42
DENSITY_T = 0.15

LOG = f'{OUT_BASE}/{NAME}.log'
logf = open(LOG, 'w')
def log(m):
    s = f'[{time.strftime("%H:%M:%S")}] {m}'
    print(s, flush=True); logf.write(s+'\n'); logf.flush()

t0 = time.time()
log('=== cNMF K=65 add-on (60 workers) from counts.h5ad ===')

work_dir = f'{OUT_BASE}/cnmf_work_k65'
if os.path.exists(work_dir):
    raise FileExistsError("Choose a new output directory; existing cNMF results are not deleted.")

from cnmf import cNMF
cnmf_obj = cNMF(output_dir=work_dir, name=NAME)
cnmf_obj.prepare(counts_fn=COUNTS, components=K_RANGE, n_iter=N_ITER, seed=SEED, num_highvar_genes=N_TOP_GENES, densify=False)
log('prepared, factorize 60 workers ...')

import multiprocessing as mp
def factor_worker(wid):
    obj = cNMF(output_dir=work_dir, name=NAME)
    obj.factorize(worker_i=wid, total_workers=N_WORKERS)

t_fact = time.time()
procs = []
for w in range(N_WORKERS):
    p = mp.Process(target=factor_worker, args=(w,)); p.start(); procs.append(p)
for p in procs: p.join()
log(f'factorize done | {(time.time()-t_fact)/60:.1f}min')

cnmf_obj.combine()
log('combine done')

import matplotlib; matplotlib.use('Agg')
k = 65
log(f'consensus K={k}')
cnmf_obj.consensus(k=k, density_threshold=DENSITY_T, show_clustering=True, close_clustergram_fig=True)
gep_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.gene_spectra_score.k_{k}.dt_*.txt')
usage_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.usages.k_{k}.dt_*.consensus.txt')
if gep_cand:
    pd.read_csv(gep_cand[0], sep='\t', index_col=0).to_csv(f'{OUT_BASE}/{NAME}_k{k}_factor_loadings.tsv', sep='\t')
    log(f'  k={k} GEP exported')
if usage_cand:
    pd.read_csv(usage_cand[0], sep='\t', index_col=0).to_csv(f'{OUT_BASE}/{NAME}_k{k}_cell_scores.tsv', sep='\t')

log(f'DONE | {(time.time()-t0)/60:.1f} min')
logf.close()
