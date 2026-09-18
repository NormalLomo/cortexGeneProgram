#!/usr/bin/env python
"""
04 (I/O helper) — bridges raw rich bin50 h5ad <-> R sctransform, and FOLDS the SCT result
back INTO the rich h5ad so <chip>_bin50.h5ad becomes the single all-in-one "完全体". Modes:

  export : <chip>_bin50.h5ad RAW counts (layers['counts']) --> genes x cells mtx + genes/bcs
           (input for 04_sct.R). NOTE reads layers['counts'] (raw), NOT X — after folding X=SCT.

  build  : preserve the rich input and write SCT-corrected counts to --output:
             X                  = SCT-corrected depth-equalized counts, on the FULL rich gene
                                  axis (genes sctransform dropped -> implicit 0)
             layers['counts']   = raw counts (preserved)
             layers['sct_log1p']= log1p(X)
             var['sct_modeled'] = bool (True for genes sctransform kept)
             uns['sct']         = {scale_factor, n_modeled, ...}
           obs / obsm / other uns (RCTD weights, cell_count_proxy, celltype_count_frac,
           domain, rctd_cell_type_names, ...) are left untouched.
           SCT-corrected source: --corrected <mtx> (fresh R run) OR --sct <existing _sct.h5ad> (backfill).

Env: cellist.
"""
import os, argparse
import numpy as np, scipy.sparse as sp, scipy.io as sio
import anndata as ad


def _raw_from_rich(a):
    X = a.layers['counts'] if 'counts' in a.layers else a.X
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def do_export(chip, rich, tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    a = ad.read_h5ad(rich)
    X = _raw_from_rich(a)                                # RAW counts (cells x G)
    gxc = X.T.tocsc()
    gxc.data = np.rint(gxc.data)
    gxc = gxc.astype(np.int32)                           # genes x cells, integer
    sio.mmwrite(os.path.join(tmpdir, f'{chip}_counts.mtx'), gxc)
    with open(os.path.join(tmpdir, f'{chip}_genes.txt'), 'w') as fh:
        fh.write('\n'.join(map(str, a.var_names)) + '\n')
    with open(os.path.join(tmpdir, f'{chip}_bcs.txt'), 'w') as fh:
        fh.write('\n'.join(map(str, a.obs_names)) + '\n')
    print(f'EXPORT\t{chip}\t{X.shape[0]} cells x {X.shape[1]} genes (raw)', flush=True)


def _corrected_from_mtx(corrected_mtx):
    cg = sio.mmread(corrected_mtx).tocsr()               # genes(kept) x cells
    genes = [l.strip() for l in open(corrected_mtx + '.genes')]
    bcs   = [l.strip() for l in open(corrected_mtx + '.bcs')]
    return cg.T.tocsr().astype(np.float32), genes, bcs   # cells x genes


def _corrected_from_sct_h5ad(sct_h5ad):
    s = ad.read_h5ad(sct_h5ad)
    X = s.X.tocsr().astype(np.float32) if sp.issparse(s.X) else sp.csr_matrix(s.X, dtype=np.float32)
    return X, list(s.var_names), list(s.obs_names)       # cells x genes (SCT subset)


def do_build(chip, rich, corrected_mtx, sct_h5ad, residuals_h5, scale_factor, output):
    if os.path.realpath(output) == os.path.realpath(rich):
        raise ValueError("--output must differ from the read-only --rich input")
    a = ad.read_h5ad(rich)
    raw = _raw_from_rich(a)                               # cells x G (full gene axis)
    a.layers['counts'] = raw                             # preserve raw

    if corrected_mtx:
        Xc, sgenes, sbcs = _corrected_from_mtx(corrected_mtx)
    else:
        Xc, sgenes, sbcs = _corrected_from_sct_h5ad(sct_h5ad)

    # align corrected cells -> rich obs order
    pos = {b: i for i, b in enumerate(sbcs)}
    order = np.array([pos[b] for b in a.obs_names])
    Xc = Xc[order]

    # map corrected genes -> rich var positions; project onto full gene axis (dropped -> 0)
    vpos = {g: i for i, g in enumerate(a.var_names)}
    keep = [j for j, g in enumerate(sgenes) if g in vpos]
    gidx = np.array([vpos[sgenes[j]] for j in keep], dtype=np.int64)
    Xc = Xc[:, keep].tocoo()
    G = a.n_vars
    full = sp.csr_matrix((Xc.data, (Xc.row, gidx[Xc.col])), shape=(a.n_obs, G), dtype=np.float32)

    a.X = full
    log = full.copy(); log.data = np.log1p(log.data)
    a.layers['sct_log1p'] = log
    modeled = np.zeros(G, dtype=bool); modeled[gidx] = True

    # --- optional: fold Pearson residuals from 04_sct.R (dense, full gene axis) ---
    # Robust to either disk layout (hdf5r/rhdf5 may store as (genes,cells) or (cells,genes));
    # chunked read converts float64->float32 on the fly to keep memory bounded.
    n_res = 0
    if residuals_h5 and os.path.exists(residuals_h5):
        import h5py
        with h5py.File(residuals_h5, 'r') as fh:
            rgenes_raw = fh['genes'][:]
            rbcs_raw   = fh['bcs'][:]
            ds = fh['residuals']
            n_g_disk, n_c_disk = len(rgenes_raw), len(rbcs_raw)
            R = np.empty((n_c_disk, n_g_disk), dtype=np.float32)
            if ds.shape == (n_g_disk, n_c_disk):           # disk: genes x cells -> transpose into R
                CHUNK = 1024
                for s in range(0, n_g_disk, CHUNK):
                    e = min(s + CHUNK, n_g_disk)
                    R[:, s:e] = ds[s:e, :].astype(np.float32).T
            elif ds.shape == (n_c_disk, n_g_disk):          # disk already cells x genes
                CHUNK = 8192
                for s in range(0, n_c_disk, CHUNK):
                    e = min(s + CHUNK, n_c_disk)
                    R[s:e, :] = ds[s:e, :].astype(np.float32)
            else:
                raise ValueError(f'residuals shape {ds.shape} != ({n_g_disk},{n_c_disk}) or ({n_c_disk},{n_g_disk})')
        rgenes = [s.decode() if isinstance(s, bytes) else str(s) for s in rgenes_raw]
        rbcs   = [s.decode() if isinstance(s, bytes) else str(s) for s in rbcs_raw]
        rpos = {b: i for i, b in enumerate(rbcs)}
        cell_order = np.array([rpos[b] for b in a.obs_names])
        R = R[cell_order]
        rgidx = np.array([vpos[g] for g in rgenes if g in vpos], dtype=np.int64)
        full_res = np.zeros((a.n_obs, G), dtype=np.float32)
        full_res[:, rgidx] = R
        a.layers['sct_residuals'] = full_res
        n_res = len(rgenes)

    a.var['sct_modeled'] = modeled
    a.uns['sct'] = {'scale_factor': int(scale_factor), 'n_modeled': int(modeled.sum()),
                    'X': 'SCT corrected depth-equalized counts (full gene axis; dropped genes=0)',
                    'raw_layer': 'counts', 'log1p_layer': 'sct_log1p',
                    'residuals_layer': 'sct_residuals' if n_res else 'NONE'}

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    tmp = output + '.tmp'
    a.write_h5ad(tmp, compression='gzip'); os.replace(tmp, output)
    extras = f' + sct_residuals({n_res} modeled genes, dense float32)' if n_res else ''
    print(f'BUILD\t{chip}\twrote SCT to {output}: X=SCT {a.n_obs}x{G} '
          f'({int(modeled.sum())} modeled genes) + layers[counts,sct_log1p]{extras}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['export', 'build'])
    ap.add_argument('--chip', required=True)
    ap.add_argument('--rich', required=True)
    ap.add_argument('--tmpdir')
    ap.add_argument('--corrected')                       # mtx from a fresh R run
    ap.add_argument('--sct')                             # existing _sct.h5ad (backfill)
    ap.add_argument('--residuals')                       # optional residuals HDF5 from 04_sct.R
    ap.add_argument('--scale_factor', default='1352')
    ap.add_argument('--output', help='Independent SCT h5ad destination; required for build')
    a = ap.parse_args()
    if a.mode == 'export':
        do_export(a.chip, a.rich, a.tmpdir)
    else:
        if not a.output:
            ap.error("build requires --output; the rich input is never overwritten")
        do_build(a.chip, a.rich, a.corrected, a.sct, a.residuals, a.scale_factor, a.output)
