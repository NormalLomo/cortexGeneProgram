#!/usr/bin/env python
"""
PER-CHIP bin50 spatial program scoring on the neighbor's SCT-CORRECTED data.

This is the CANONICAL fix replacing the buggy per-bin-CPM scoring AND it is
distinct from the earlier dense Pearson-residual attempt (04b). Here we consume
the neighbor's *pre-computed* SCT-corrected counts directly from the sct/ h5ad
files -- no residual computation, no normalization step is performed by us.

Run:  python 04_spatial_score_sct.py <CHIP_ID>

For ONE chip:
  - Load {chip}_sct.h5ad from the SCT directory.
  - Use .X  = SCT depth-equalized corrected counts (positive, sparse CSR float32,
    all chips equalized to scale_factor 1352). We do NOT use layers['log1p'] and
    we apply NO further per-bin / library normalization. SCT already depth-equalized;
    adding per-bin scaling would re-create the CPM bug. This is intentional.
  - Load NON-NEGATIVE gene_spectra_tpm (60 programs x 18742 genes, HGNC columns, min=0).
  - shared = intersection(adata.var_names, tpm columns), computed PER CHIP because SCT
    drops low-count genes so the gene set varies per chip. Spectra gene order kept.
  - score[bin, prog] = X[:, shared] . tpm[prog, shared].T   (plain non-negative dot
    product, sparse @ dense; NO per-bin normalization beyond using SCT .X as-is).
  - Write _SCT_score_perchip/{chip}.parquet : bin, majorDomain, domain(guard),
    region(guard), bin_total_umi, x, y, + program_1..program_60 (raw dot scores,
    column order = tpm.index = program_names.tsv 'program' order = P1..P60).

Atomic write (.tmp then os.replace) so a partial file never looks complete.
Resume-safe: caller skips chips whose output parquet already exists.
Thread caps (OMP/OPENBLAS/MKL) are set by the caller env to avoid BLAS oversubscription.
"""
import os, sys, gc, time
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

t0 = time.time()
def log(*a):
    print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

OUTDIR   = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")
TPM_PATH = (os.environ[(__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "")] + "/results/cnmf_snrna_joint_full1M_v1/"
            "cnmf_work/snrna_joint_full1M_v1/"
            "snrna_joint_full1M_v1.gene_spectra_tpm.k_60.dt_0_15.txt")
SCT_DIR  = os.environ["SPATIAL_SCT_DIR"]
TMPDIR   = os.path.join(OUTDIR, "_SCT_score_perchip")
CHUNK    = 50000        # bins per matmul chunk (sparse rows are light)

if len(sys.argv) != 2:
    sys.exit("usage: 04_spatial_score_sct.py <CHIP_ID>")
chip = sys.argv[1]
os.makedirs(TMPDIR, exist_ok=True)
outp = os.path.join(TMPDIR, f"{chip}.parquet")
if os.path.exists(outp):
    log(f"{chip}: output exists, skip"); sys.exit(0)

log(f"chip={chip}  threads OMP={os.environ.get('OMP_NUM_THREADS')}")

# --- program loadings (programs x genes, HGNC symbol columns, non-negative) ---
tpm = pd.read_csv(TPM_PATH, sep="\t", index_col=0)
tpm.index = tpm.index.astype(str)
program_ids = list(tpm.index)                          # '1'..'60'
tpm_genes   = list(tpm.columns)
prog_cols   = [f"program_{p}" for p in program_ids]    # program_1..program_60 (== P1..P60)
assert (tpm.to_numpy() >= 0).all(), "gene_spectra_tpm must be non-negative"

# --- load chip SCT-corrected counts ---
fp = os.path.join(SCT_DIR, f"{chip}_sct.h5ad")
A = ad.read_h5ad(fp)
n = A.n_obs
obs = A.obs
X = A.X
if not sparse.issparse(X):
    X = sparse.csr_matrix(X)
X = X.tocsr().astype(np.float32)
# sanity: SCT .X must be non-negative and finite
assert X.data.size == 0 or (np.isfinite(X.data).all() and (X.data >= 0).all()), \
    f"{chip}: SCT .X must be non-negative finite"
st_var = list(A.var_names)
del A
gc.collect()

# --- gene intersection (HGNC symbols, same space), PER CHIP ---
pos = {g: i for i, g in enumerate(st_var)}
shared = [g for g in tpm_genes if g in pos]            # keep spectra gene order
st_idx = np.array([pos[g] for g in shared], dtype=np.int64)
W  = tpm.loc[:, shared].to_numpy(dtype=np.float32)     # (60, n_shared) >= 0
WT = np.ascontiguousarray(W.T)                         # (n_shared, 60)
log(f"{chip}: {n} bins, {len(st_var)} genes; shared = {len(shared)} of {len(tpm_genes)} tpm genes")

# --- plain non-negative dot product (NO per-bin normalization; SCT .X used as-is) ---
Xsub = X[:, st_idx].tocsr()
del X
gc.collect()
sc_out = np.empty((n, len(program_ids)), dtype=np.float32)
s = 0
while s < n:
    e = min(s + CHUNK, n)
    sc_out[s:e, :] = Xsub[s:e].dot(WT)
    s = e
del Xsub
gc.collect()

# --- assemble per-chip parquet ---
df = pd.DataFrame(sc_out, columns=prog_cols)
df.insert(0, "bin", np.asarray(obs.index))
df["majorDomain"]   = obs["majorDomain"].astype(str).values
df["domain"]        = obs["domain"].astype(str).values if "domain" in obs.columns else ""
df["region"]        = obs["region"].astype(str).values if "region" in obs.columns else ""
df["bin_total_umi"] = np.asarray(obs["bin_total_umi"], dtype=np.float64)
df["x"]             = np.asarray(obs["x"])
df["y"]             = np.asarray(obs["y"])

tmp_out = outp + ".tmp"
df.to_parquet(tmp_out, index=False)
os.replace(tmp_out, outp)   # atomic
log(f"{chip}: wrote {outp}  ({n} bins x {len(prog_cols)} programs)  DONE")
