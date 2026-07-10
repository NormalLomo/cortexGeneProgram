# =============================================================================
# STEP 00 (mouse) — stream dense h5ad -> Isocortex subset CSR h5ad (BPCells-readable)
# REDO v3: 不寫巨型 .mtx（前版 6.2GB .mtx -> R readMM 整塊 -> 爆 int32 崩潰）。
#          改：stream dense X 的 Isocortex 行，分塊用 h5py 增量寫 CSR h5ad，
#          indptr 強制 int64（nnz=2.30e9 > 2^31）；全程不在記憶體持有整塊 scipy CSR。
# IN : Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad  (dense uint32, 4.4M x 21899)
#      snRNA_cellMeta.csv  (per-cell meta keyed by barcode; brain_struct=='Isocortex' ~761k)
# DOES: stream dense X (10k-row blocks by default) only for Isocortex rows -> per-block CSR ->
#       h5py append into ONE CSR h5ad (cells x genes); X group = AnnData csr_matrix spec.
# OUT: iso_subset/counts_csr.h5ad  (X only, cells x genes, indptr int64)
#      iso_subset/{barcodes.txt, genes.csv, meta.csv}
# NEXT: 01_prepare_mouse_raw_counts.py validates the emitted raw-UMI subset.
# =============================================================================
import argparse
import os, time, h5py, numpy as np, pandas as pd
import scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

parser = argparse.ArgumentParser()
parser.add_argument("--raw-h5ad", required=True)
parser.add_argument("--metadata-csv", required=True)
parser.add_argument("--output-dir", required=True)
parser.add_argument("--expected-nuclei", type=int, default=761378)
parser.add_argument("--block-cells", type=int, default=10_000)
args = parser.parse_args()
RAW = args.raw_h5ad
META = args.metadata_csv
OUT = args.output_dir
os.makedirs(OUT, exist_ok=True)
H5OUT = os.path.join(OUT, "counts_csr.h5ad")
if args.block_cells <= 0:
    raise ValueError("--block-cells must be positive")

# --- 1. meta: keep Isocortex, key on barcode (first unnamed col) ---
log("read cellMeta (selected cols)")
hdr = pd.read_csv(META, nrows=0)
bc_col = hdr.columns[0]
usecols = [bc_col, "derived_cell_libs", "library", "region", "sub_region",
           "brain_struct", "donor_id", "sex", "Annotation", "cell_class", "num_cells_postQC"]
usecols = [c for c in usecols if c in hdr.columns]
meta = pd.read_csv(META, usecols=usecols)
meta = meta.rename(columns={bc_col: "barcode"})
meta = meta[meta["brain_struct"] == "Isocortex"].copy()
meta = meta.set_index("barcode")
log(f"Isocortex meta rows: {len(meta)}; libs={meta['library'].nunique()}")

# --- 2. h5ad obs index (barcodes) + var ---
log("open raw h5ad (backed)")
f = h5py.File(RAW, "r")
obs_bc = np.array([b.decode() if isinstance(b, bytes) else b for b in f["obs/_index"][:]])
var_idx = np.array([b.decode() if isinstance(b, bytes) else b for b in f["var/_index"][:]])
gn_cats = np.array([b.decode() if isinstance(b, bytes) else b for b in f["var/gene_name/categories"][:]])
gn_codes = f["var/gene_name/codes"][:]
gene_name = np.where(gn_codes >= 0, gn_cats[gn_codes.clip(min=0)], "")
n_cells, n_genes = f["X"].shape
log(f"raw X: {n_cells} x {n_genes}; obs_bc={len(obs_bc)} var={len(var_idx)}")

# --- 3. row index of Isocortex cells (preserve h5ad order) ---
bc_to_row = pd.Series(np.arange(len(obs_bc)), index=obs_bc)
present = meta.index[meta.index.isin(bc_to_row.index)]
rows = bc_to_row.loc[present].values
order = np.argsort(rows)
rows_sorted = rows[order]
bc_sorted = present.values[order]
n_keep = len(rows_sorted)
log(f"Isocortex cells present in h5ad: {n_keep} / {len(meta)}")
if n_keep != args.expected_nuclei:
    raise ValueError(f"Isocortex subset contains {n_keep} nuclei; expected {args.expected_nuclei}")

# --- 4. create CSR h5ad skeleton (AnnData csr_matrix spec), append per block ---
# X: cells x genes ; indptr int64 (nnz may exceed 2^31), indices int32, data uint32
log("create CSR h5ad skeleton (resizable datasets)")
hf = h5py.File(H5OUT, "w")
hf.attrs["encoding-type"] = "anndata"
hf.attrs["encoding-version"] = "0.1.0"
X = hf.create_group("X")
X.attrs["encoding-type"] = "csr_matrix"
X.attrs["encoding-version"] = "0.1.0"
X.attrs["shape"] = np.array([n_keep, n_genes], dtype="int64")
d_data = X.create_dataset("data", shape=(0,), maxshape=(None,), dtype="uint32",
                          chunks=(2**20,), compression="gzip", compression_opts=2)
d_ind  = X.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32",
                          chunks=(2**20,), compression="gzip", compression_opts=2)
# indptr length = n_keep+1, int64
d_ptr  = X.create_dataset("indptr", shape=(n_keep + 1,), dtype="int64")
d_ptr[0] = 0

# minimal obs/var so AnnData/BPCells can read group="X"
obs_g = hf.create_group("obs")
obs_g.attrs["encoding-type"] = "dataframe"; obs_g.attrs["encoding-version"] = "0.2.0"
obs_g.attrs["_index"] = "_index"; obs_g.attrs["column-order"] = np.array([], dtype="S1")
obs_g.create_dataset("_index", data=np.array([b.encode() for b in bc_sorted]))
var_g = hf.create_group("var")
var_g.attrs["encoding-type"] = "dataframe"; var_g.attrs["encoding-version"] = "0.2.0"
var_g.attrs["_index"] = "_index"; var_g.attrs["column-order"] = np.array([], dtype="S1")
var_g.create_dataset("_index", data=np.array([str(g).encode() for g in var_idx]))

# --- 5. stream dense X rows -> per-block CSR -> append (never hold full CSR) ---
log("stream dense X (Isocortex rows) -> append CSR h5ad")
Xraw = f["X"]
# 10,000 dense uint32 cells x 21,899 genes is about 0.82 GiB before temporary
# arrays and sparse conversion; increase only after measuring available memory.
BLK = args.block_cells
sel_mask = np.zeros(n_cells, dtype=bool); sel_mask[rows_sorted] = True
nnz_total = 0
ptr_offset = 0          # running cumulative nnz for indptr
written_cells = 0
# map: global h5ad row order of kept cells is already increasing (rows_sorted sorted)
for start in range(0, n_cells, BLK):
    end = min(start + BLK, n_cells)
    local = np.where(sel_mask[start:end])[0]
    if local.size == 0:
        continue
    sub = Xraw[start:end, :][local, :]            # dense uint32 block, kept rows in increasing order
    m = sp.csr_matrix(sub)                         # cells x genes, block-local (small, int32 fine)
    bn = m.data.shape[0]
    # append data + indices
    d_data.resize((nnz_total + bn,)); d_data[nnz_total:nnz_total + bn] = m.data.astype("uint32")
    d_ind.resize((nnz_total + bn,));  d_ind[nnz_total:nnz_total + bn]  = m.indices.astype("int32")
    # indptr: m.indptr is length (rows+1), local cumulative; shift by ptr_offset, write rows 1..rows
    nrows_blk = m.shape[0]
    d_ptr[written_cells + 1: written_cells + 1 + nrows_blk] = (m.indptr[1:].astype("int64") + ptr_offset)
    nnz_total += bn
    ptr_offset += bn
    written_cells += nrows_blk
    log(f"  block {start}:{end} kept {local.size} rows; nnz_total={nnz_total} cells={written_cells}")
    del sub, m
f.close()
hf.flush()
log(f"assembled CSR h5ad: shape=({written_cells},{n_genes}) nnz={nnz_total}")
assert written_cells == n_keep, f"row mismatch {written_cells} != {n_keep}"
assert int(d_ptr[-1]) == nnz_total, f"indptr tail {int(d_ptr[-1])} != nnz {nnz_total}"
hf.close()

# --- 6. sidecar files for R 01 ---
log("write barcodes/genes/meta sidecars")
with open(os.path.join(OUT, "barcodes.txt"), "w") as fh:
    fh.write("\n".join(bc_sorted) + "\n")
pd.DataFrame({"gene": var_idx, "gene_name": gene_name}).to_csv(
    os.path.join(OUT, "genes.csv"), index=False)
meta.loc[bc_sorted].to_csv(os.path.join(OUT, "meta.csv"))
log("DONE 00_fix_subset_mouse (CSR h5ad, BPCells-readable)")
