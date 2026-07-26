#!/usr/bin/env python3
import os, sys, math
import numpy as np
import scipy.sparse as sp
PROJECT_ROOT = os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')
RESDIR = os.path.join(PROJECT_ROOT, 'results')
WORK = os.path.join(RESDIR, 'ic_kselection')
COUNTS = os.path.join(RESDIR, 'cnmf_snrna_joint_full1M_v1', 'counts.h5ad')
HVG = os.path.join(RESDIR, 'cnmf_snrna_joint_full1M_v1', 'cnmf_work', 'snrna_joint_full1M_v1', 'snrna_joint_full1M_v1.overdispersed_genes.txt')
GPU_PATCH_DIR = os.environ.get('CORTEX_GPU_CNMF_PATCH', os.path.join(RESDIR, 'cnmf_joint_xspecies_v1', 'gpu_cnmf'))
N_SUB = 150000
SEED_SUB = 42
K_GRID = [30, 40, 50, 55, 60, 65, 70, 80, 90]
RESTARTS = [0, 1, 2]
MAX_ITER = 1500
TOL = 1e-05
EPS = 1e-09

def main():
    sys.path.insert(0, GPU_PATCH_DIR)
    import torch
    import gpu_cnmf_patch as gp
    if not torch.cuda.is_available():
        raise RuntimeError()
    dev = 'cuda'
    with open(HVG) as f:
        hvg = [l.strip() for l in f if l.strip()]
    import h5py
    with h5py.File(COUNTS, 'r') as h:
        Xg = h['X']
        shape = tuple((int(s) for s in Xg.attrs['shape']))
        (n_total, n_var) = shape
        vidx_key = h['var'].attrs['_index']
        var_names = h['var'][vidx_key].asstr()[:] if hasattr(h['var'][vidx_key], 'asstr') else np.array([x.decode() if isinstance(x, bytes) else x for x in h['var'][vidx_key][:]])
        var_pos = {g: i for (i, g) in enumerate(var_names)}
        gcols = np.array([var_pos[g] for g in hvg if g in var_pos], dtype=np.int64)
        n_genes = len(gcols)
        remap = np.full(n_var, -1, dtype=np.int32)
        remap[gcols] = np.arange(n_genes, dtype=np.int32)
        rng = np.random.default_rng(SEED_SUB)
        n_sub = min(N_SUB, n_total)
        idx = np.sort(rng.choice(n_total, size=n_sub, replace=False))
        indptr_ds = Xg['indptr']
        data_ds = Xg['data']
        indices_ds = Xg['indices']
        indptr = indptr_ds[:]
        out_indptr = np.zeros(n_sub + 1, dtype=np.int64)
        data_chunks = []
        col_chunks = []
        CHUNK = 5000
        ridx = 0
        for c0 in range(0, n_sub, CHUNK):
            c1 = min(c0 + CHUNK, n_sub)
            rows = idx[c0:c1]
            lo = int(indptr[rows[0]])
            hi = int(indptr[rows[-1] + 1])
            blk_data = data_ds[lo:hi]
            blk_idx = indices_ds[lo:hi]
            for r in rows:
                (s, e) = (int(indptr[r]), int(indptr[r + 1]))
                cols = blk_idx[s - lo:e - lo]
                vals = blk_data[s - lo:e - lo]
                mapped = remap[cols]
                keep = mapped >= 0
                kc = mapped[keep]
                kv = vals[keep]
                if kc.size:
                    col_chunks.append(kc.astype(np.int32))
                    data_chunks.append(kv.astype(np.float32))
                out_indptr[ridx + 1] = out_indptr[ridx] + int(kc.size)
                ridx += 1
        all_cols = np.concatenate(col_chunks) if col_chunks else np.zeros(0, np.int32)
        all_data = np.concatenate(data_chunks) if data_chunks else np.zeros(0, np.float32)
        X = sp.csr_matrix((all_data, all_cols, out_indptr), shape=(n_sub, n_genes))
    X = X.astype(np.float32)
    (n_cells, n_genes) = X.shape
    n_obs = n_cells * n_genes
    Vd = np.ascontiguousarray(X.toarray(), dtype=np.float32)
    V = torch.from_numpy(Vd).to(dev)
    del Vd
    Vnorm2 = float((V * V).sum().item())
    lgam_const = float(torch.lgamma(V + 1.0).sum().item())
    log_n_obs = math.log(n_obs)
    rows = []
    for K in K_GRID:
        best = None
        for s in RESTARTS:
            (H, W, n_it, relF) = gp.gpu_munmf(X, K=K, max_iter=MAX_ITER, seed=s, tol=TOL, device=dev)
            if best is None or relF < best['relF']:
                best = {'H': H, 'W': W, 'n_it': n_it, 'relF': relF, 'seed': s}
        gp.clear_V_cache()
        H = best['H']
        W = best['W']
        relF = best['relF']
        Ht = torch.from_numpy(H).to(dev)
        Wt = torch.from_numpy(W).to(dev)
        WH = Wt @ Ht
        resid = V - WH
        RSS = float((resid * resid).sum().item())
        sigma2 = RSS / n_obs
        sigma2 = max(sigma2, EPS)
        loglik_gauss = -0.5 * n_obs * (math.log(2.0 * math.pi * sigma2) + 1.0)
        Lam = torch.clamp(WH, min=EPS)
        loglik_pois = float((V * torch.log(Lam) - Lam).sum().item()) - lgam_const
        k_params = (n_cells + n_genes) * K
        Wpos = torch.clamp(Wt, min=0.0)
        rowsum = Wpos.sum(dim=1, keepdim=True) + EPS
        tau = Wpos / rowsum
        EN = float(-(tau * torch.log(tau + EPS)).sum().item())

        def crit(loglik):
            bic = -2.0 * loglik + k_params * log_n_obs
            aic = -2.0 * loglik + 2.0 * k_params
            mdl = -loglik + 0.5 * k_params * log_n_obs
            icl = bic + 2.0 * EN
            return (bic, aic, mdl, icl)
        (bic_g, aic_g, mdl_g, icl_g) = crit(loglik_gauss)
        (bic_p, aic_p, mdl_p, icl_p) = crit(loglik_pois)
        relF_recomp = math.sqrt(RSS) / math.sqrt(max(Vnorm2, EPS))
        row = dict(K=K, best_seed=best['seed'], n_iter=best['n_it'], relFrob=relF, relFrob_recomp=relF_recomp, RSS=RSS, k_params=k_params, n_obs=n_obs, entropy_EN=EN, loglik_gauss=loglik_gauss, loglik_pois=loglik_pois, BIC_gauss=bic_g, AIC_gauss=aic_g, MDL_gauss=mdl_g, ICL_gauss=icl_g, BIC_pois=bic_p, AIC_pois=aic_p, MDL_pois=mdl_p, ICL_pois=icl_p)
        rows.append(row)
        del Ht, Wt, WH, resid, Lam, Wpos, tau
        torch.cuda.empty_cache()
    import pandas as pd
    df = pd.DataFrame(rows).sort_values('K').reset_index(drop=True)
    tsv = os.path.join(WORK, 'ic_kselection_results.tsv')
    df.to_csv(tsv, sep='\t', index=False)

if __name__ == '__main__':
    main()
