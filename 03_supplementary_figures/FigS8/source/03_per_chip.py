#!/usr/bin/env python
"""
Per-chip lattice-ring program x program spatial MARK CROSS-CORRELATION (Stoyan k_mm(r)).

For ONE chip (one parquet row group):
  points = bins at lattice (ix,iy)=(x-xmin)/50 ; marks = 60 program-z vectors (SCT).
  For each ring r-bin k in {1..R} lattice units (k*50 um), gather all bin PAIRS (i,j)
  separated by Euclidean lattice distance ~= k, and compute the full 60x60 cross-product
  matrix in ONE matmul:  S_k = sum_pairs  z_i (outer) z_j  ,  count n_k = #pairs.
  We accumulate BOTH directions (offset and its negative) so S_k is symmetrized over
  ordered pairs; thus k_mm is symmetric in A,B at the pair level.

Mark cross-correlation (normalized):
  kbar_AB(r) = ( S_k[A,B] / n_k ) / ( mu_A * mu_B )
where mu_A = global mean of program A z over this chip's bins.
(Stoyan k_mm normalization by product of mean marks. With z-scored marks mu can be ~0,
 so we ALSO store the raw mean cross-product  C_k[A,B]=S_k[A,B]/n_k  which is the
 distance-binned cross-COVARIANCE of the marks -- the robust quantity used downstream.)

Neighbor pairing is done by a vectorized lattice hash-join (searchsorted), NO O(N^2),
NO python per-pair loop. One matmul per ring offset.

Writes: <outdir>/<chip>.npz with
  C   : (R,60,60) mean cross-product per ring (cross-covariance of marks)
  Sii : (R,60)    auto term sum_i z_i[A]*z_i[A]/n  per ring (diag self, for ref)
  npairs : (R,) number of ordered pairs per ring
  mu  : (60,) mean program z
  sd  : (60,) std program z
  nbins, chip, region
Usage: spatial_crosscorr_perchip.py <chip_id> <region> <row_group_index>
"""
import sys, os, time
import numpy as np
import pyarrow.parquet as pq

PARQUET = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/spatial_bin50_program_score_SCT.parquet"
OUTDIR  = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/spatial_crosscorr/_perchip"
R = 10          # max ring in lattice units -> 500 um
NPROG = 60
LAT = 50        # lattice spacing (um)

def ring_offsets(R):
    import collections
    rm = collections.defaultdict(list)
    for dx in range(-R, R+1):
        for dy in range(-R, R+1):
            if dx == 0 and dy == 0: continue
            d = np.hypot(dx, dy)
            k = int(round(d))
            if 1 <= k <= R:
                rm[k].append((dx, dy))
    return rm

def main():
    chip = sys.argv[1]; region = sys.argv[2]; rg = int(sys.argv[3])
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"{chip}.npz")
    if os.path.exists(out):
        print(f"[{chip}] already done, skip"); return
    t0 = time.time()
    pf = pq.ParquetFile(PARQUET)
    cols = ['bin','x','y'] + [f'program_{i}' for i in range(1, NPROG+1)]
    tbl = pf.read_row_group(rg, columns=cols)
    bins = tbl.column('bin').to_pylist()
    chips_here = set(b.split('_')[0] for b in bins)
    assert chips_here == {chip}, f"row group {rg} chip mismatch {chips_here} vs {chip}"
    x = np.asarray(tbl.column('x'), dtype=np.int64)
    y = np.asarray(tbl.column('y'), dtype=np.int64)
    Z = np.empty((len(bins), NPROG), dtype=np.float64)
    for i in range(NPROG):
        Z[:, i] = np.asarray(tbl.column(f'program_{i+1}'), dtype=np.float64)
    del tbl
    # NaN guard
    nan_mask = np.isnan(Z).any(axis=1)
    if nan_mask.any():
        keep = ~nan_mask
        x=x[keep]; y=y[keep]; Z=Z[keep]
        print(f"[{chip}] dropped {nan_mask.sum()} NaN bins")
    n = len(x)
    mu = Z.mean(axis=0); sd = Z.std(axis=0)

    # lattice integer coords
    ix = ((x - x.min()) // LAT).astype(np.int64)
    iy = ((y - y.min()) // LAT).astype(np.int64)
    W = int(iy.max()) + 1               # encode key = ix*(H) + iy with H>max(iy)
    H = W + 2*R + 5
    key = ix * H + iy
    order = np.argsort(key, kind='stable')
    key_s = key[order]
    Z_s = Z[order]                       # marks sorted by key for searchsorted matching

    rm = ring_offsets(R)
    C = np.zeros((R, NPROG, NPROG), dtype=np.float64)     # sum of outer products per ring
    npairs = np.zeros(R, dtype=np.int64)

    for k in range(1, R+1):
        Sk = np.zeros((NPROG, NPROG), dtype=np.float64)
        cnt = 0
        for (dx, dy) in rm[k]:
            shift = dx * H + dy
            tgt = key_s + shift                      # neighbor key for each source bin
            pos = np.searchsorted(key_s, tgt)
            valid = pos < len(key_s)
            pos2 = pos[valid]
            hit = key_s[pos2] == tgt[valid]
            src_idx = np.nonzero(valid)[0][hit]      # index into sorted arrays (source)
            dst_idx = pos2[hit]                      # matched neighbor index (sorted)
            if src_idx.size == 0:
                continue
            # one matmul: sum_pairs z_src (outer) z_dst  -> (60,60)
            Sk += Z_s[src_idx].T @ Z_s[dst_idx]
            cnt += src_idx.size
        if cnt > 0:
            C[k-1] = Sk / cnt
        npairs[k-1] = cnt

    np.savez_compressed(out,
                        C=C.astype(np.float32),
                        npairs=npairs,
                        mu=mu.astype(np.float32),
                        sd=sd.astype(np.float32),
                        nbins=np.int64(n),
                        chip=chip, region=region)
    print(f"[{chip}] region={region} nbins={n} npairs(r=50um)={npairs[0]} "
          f"done in {time.time()-t0:.1f}s -> {out}", flush=True)

if __name__ == "__main__":
    main()
