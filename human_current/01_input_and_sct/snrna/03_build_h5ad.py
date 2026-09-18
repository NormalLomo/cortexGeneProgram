# =============================================================================
# STEP 03 / 03  —  assemble per-donor chunks into the final AnnData .h5ad
# IN : chunks_1M/ (per-donor .mtx) + snRNA_1M_obs.csv + snRNA_1M_var.csv  (this dir, from 02)
# DOES: mmread each chunk -> scipy hstack -> transpose to cells x genes -> AnnData -> write
# OUT: snRNA_1M.h5ad  (this dir)
# PREV: 02_prep_full_1M.R
# NOTE: assembly is in Python (NOT merged into the R step 02) because the full matrix has
#   ~5.18e9 nonzeros > 2^31, exceeding R dgCMatrix's int32 index limit.
# =============================================================================
"""Assemble the 1M snRNA matrix from per-donor mtx chunks into AnnData .h5ad."""
import numpy as _np_patch, pandas as _pd_patch
_pd_patch.Series.nonzero = lambda self: (_np_patch.asarray(self).nonzero()[0],)
import os, time, numpy as np, pandas as pd, anndata as ad
import scipy.io as sio, scipy.sparse as sp, glob

T0 = time.time()
INP = os.environ["SNRNA_WORK_DIR"]
LOG = open(os.path.join(INP, "build_h5ad.log"), "w")
def log(m): line=f"[{time.time()-T0:7.1f}s] {m}"; print(line); LOG.write(line+"\n"); LOG.flush()

chunk_dir = os.path.join(INP, "chunks_1M")
chunks = sorted(glob.glob(os.path.join(chunk_dir, "counts_*.mtx")))

log(f"=== build snRNA_1M.h5ad")
log(f"Found {len(chunks)} chunk mtx files")

# Assemble cells × genes matrix from chunks
mats = []
barcodes_all = []
for mtx_fn in chunks:
    donor = os.path.basename(mtx_fn).replace("counts_", "").replace(".mtx", "")
    m = sio.mmread(mtx_fn).tocsc()
    bc = [l.strip() for l in open(mtx_fn.replace("counts_", "barcodes_").replace(".mtx", ".txt"))]
    mats.append(m)
    barcodes_all.extend(bc)
    log(f"  {donor}: {m.shape}, bc={len(bc)}")

# Stack along cells axis (column = cells in genes×cells), then transpose to cells × genes
log("hstack chunks")
m_full = sp.hstack(mats).tocsr()
m_full = m_full.T.tocsr()
log(f"Assembled cells×genes: {m_full.shape}")

obs = pd.read_csv(os.path.join(INP, "snRNA_1M_obs.csv"), index_col=0)
var = pd.read_csv(os.path.join(INP, "snRNA_1M_var.csv"))
# Reorder obs by barcodes_all (chunks are donor-ordered)
obs = obs.loc[barcodes_all]
adata = ad.AnnData(X=m_full, obs=obs, var=pd.DataFrame(index=var["gene"]))
log(f"AnnData: {adata.shape}")
adata.obs["cohort"] = adata.obs["donor"].apply(lambda d: "S" if str(d).startswith("S") else "H")
log(f"Cohort: {dict(adata.obs['cohort'].value_counts())}")

adata.write(os.path.join(INP, "snRNA_1M.h5ad"))
log("Wrote snRNA_1M.h5ad")
LOG.close()
