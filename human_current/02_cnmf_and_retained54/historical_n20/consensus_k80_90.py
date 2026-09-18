#!/usr/bin/env python3
"""Historical n20 consensus for ranks K=80 and K=90, reusing factorization outputs.
This does not recover an 80-run batch or recreate the adopted 80+20 fusion.
"""
import os, time, glob
import pandas as pd
import matplotlib; matplotlib.use('Agg')
from cnmf import cNMF

OUT_BASE = os.environ["CORTEX_PROGRAM_ROOT"] + "/results/cnmf_snrna_joint_full1M_v1"
NAME = 'snrna_joint_full1M_v1'
work_dir = f'{OUT_BASE}/cnmf_work'
DENSITY_T = 0.15
KS = [80, 90]

LOG = f'{OUT_BASE}/{NAME}_consensus_k80_90.log'
logf = open(LOG, 'w')
def log(m):
    s = f'[{time.strftime("%H:%M:%S")}] {m}'
    print(s, flush=True); logf.write(s+'\n'); logf.flush()

assert os.path.exists(f'{work_dir}/{NAME}'), f'missing {work_dir}/{NAME}'
log(f'=== supplement consensus K={KS} (reuse existing factorize, no rmtree) ===')
cnmf_obj = cNMF(output_dir=work_dir, name=NAME)
for k in KS:
    log(f'consensus K={k}')
    try:
        cnmf_obj.consensus(k=k, density_threshold=DENSITY_T, show_clustering=True, close_clustergram_fig=True)
        gep_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.gene_spectra_score.k_{k}.dt_*.txt')
        usage_cand = glob.glob(f'{work_dir}/{NAME}/{NAME}.usages.k_{k}.dt_*.consensus.txt')
        if gep_cand:
            pd.read_csv(gep_cand[0], sep='\t', index_col=0).to_csv(f'{OUT_BASE}/{NAME}_k{k}_factor_loadings.tsv', sep='\t')
            log(f'  k={k} GEP exported')
        if usage_cand:
            pd.read_csv(usage_cand[0], sep='\t', index_col=0).to_csv(f'{OUT_BASE}/{NAME}_k{k}_cell_scores.tsv', sep='\t')
            log(f'  k={k} usages exported')
    except Exception as e:
        log(f'  K={k} err: {e}')
log('DONE')
logf.close()
