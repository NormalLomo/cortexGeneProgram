#!/usr/bin/env python
"""01_torus_shift_perchip.py — Per-chip torus-shift null distribution.

For each chip and each (A, B) pair at r=25µm, compute:
  - T_obs = log2 g_AB(25µm) [from existing per_chip NPZ]
  - null_T[k] for k=1..N_PERM via torus-shift of mark A
  - Z_i = (T_obs - mu_null) / sd_null  (Welford FP64)
  - n_in_tissue = N_bins for that chip

Output: betweenchip_{mode}_per_chip_Z.tsv
  columns: mode, A_name, B_name, chip_id, T_obs, Z_i, n_in_tissue, mu_null, sd_null

Usage:
  python 01_torus_shift_perchip.py --mode cellprog [--smoke] [--chip-start i] [--chip-end j]
  python 01_torus_shift_perchip.py --mode progprog [--smoke]
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    _HAS_CUPY = True
except Exception:
    cp = None; csp = None; _HAS_CUPY = False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = "CORTEX_PROGRAM_ROOT"
RESDIR = os.path.join(ROOT, "results", "crossregion_v1")
PERCHIP_DIR = os.path.join(RESDIR, "markcorr", "per_chip")
OUTDIR = os.path.join(RESDIR, "markcorr_betweenchip_v1")
F_META = os.path.join(RESDIR, "spatial_bin50_meta.parquet")
F_SCT  = os.path.join(RESDIR, "spatial_bin50_program_score_SCT.parquet")
F_RCTD = os.path.join(RESDIR, "spatial_bin50_rctd_weights.parquet")

# ---------------------------------------------------------------------------
# Constants (match markcorr_runner_lib.py)
# ---------------------------------------------------------------------------
CELLTYPES = ['AST','CHANDELIER','ENDO','ET','L2-L3 IT LINC00507','L3-L4 IT RORB',
             'L4-L5 IT RORB','L6 CAR3','L6 CT','L6 IT','L6B','LAMP5','MICRO','NDNF',
             'NP','OLIGO','OPC','PAX6','PVALB','SST','VIP','VLMC']  # 22
# The public analysis retains 54 mapped programs from the 60-column NPZ input.
N_PROG_SPEC = 54
PROGRAMS = [f'program_{i}' for i in range(1, N_PROG_SPEC + 1)]  # 54
PROGRAMS_60 = [f'program_{i}' for i in range(1, 61)]            # for parquet load
# The retained program set is cnmf {1..60}\{9,18,19,35,52,57}.
EXCLUDED_CNMF = {9, 18, 19, 35, 52, 57}
KEPT_CNMF = [i for i in range(1, 61) if i not in EXCLUDED_CNMF]
KEPT_PROGRAMS = [f'program_{i}' for i in KEPT_CNMF]  # 54 columns aligned with TableS1.new_P P1..P54
assert len(KEPT_PROGRAMS) == 54, f'KEPT_PROGRAMS expected 54, got {len(KEPT_PROGRAMS)}'

UMI_FLOOR = 100
# bin50 grid step = 50 DNB px = 25µm
GRID_STEP_PX = 50.0   # DNB px per grid unit
R_FOCAL_UM   = 25.0   # focal ring lower bound (µm)
R_FOCAL_PX   = R_FOCAL_UM * 2  # 50 DNB px
# Ring edges in px (DNB; 1µm = 2 px)
RING_EDGES_PX = tuple(e * 2 for e in (0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500))
# Ring index 0 = [0, 50px) = [0, 25µm) = within-bin + nearest neighbours
R_IDX = 0

N_PERM = 1000
SEED = 42
# Minimum shift in grid units (must skip at least 1 grid unit = 50px = 25µm)
MIN_SHIFT_GU = 1

# ---------------------------------------------------------------------------
# GPU helpers
# ---------------------------------------------------------------------------
def _to_numpy(a):
    if _HAS_CUPY and isinstance(a, cp.ndarray):
        return cp.asnumpy(a)
    return np.asarray(a)

def _free_gpu():
    if _HAS_CUPY:
        try:
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except: pass

# ---------------------------------------------------------------------------
# Ring adjacency (single ring at r_focal) — half-open annulus [r_lo, r_hi)
# coords are in PIXEL units (DNB px)
# ---------------------------------------------------------------------------
def build_focal_ring_cpu(coords_px: np.ndarray, r_lo: float, r_hi: float) -> sp.csr_matrix:
    """Build sparse CSR adjacency for the half-open annulus [r_lo, r_hi).
    coords_px: (N, 2) float, in DNB px units.
    Returns CSR matrix on CPU."""
    N = len(coords_px)
    tree = cKDTree(coords_px)
    eps = 1e-9 * max(1.0, r_hi)
    outer_pairs = tree.query_pairs(r=r_hi - eps, output_type='ndarray')  # dist < r_hi
    if r_lo > 0:
        inner_pairs = tree.query_pairs(r=r_lo - eps, output_type='ndarray')  # dist < r_lo
        outer_set = set(map(tuple, outer_pairs)) if len(outer_pairs) else set()
        inner_set = set(map(tuple, inner_pairs)) if len(inner_pairs) else set()
        ring_pairs = np.array(sorted(outer_set - inner_set), dtype=np.int64)
    else:
        # [0, r_hi): include self-pairs + nearest neighbours
        self_idx = np.arange(N)
        if len(outer_pairs):
            rows = np.concatenate([outer_pairs[:, 0], outer_pairs[:, 1], self_idx])
            cols = np.concatenate([outer_pairs[:, 1], outer_pairs[:, 0], self_idx])
        else:
            rows = self_idx; cols = self_idx
        data = np.ones(len(rows), dtype=np.float32)
        return sp.csr_matrix((data, (rows, cols)), shape=(N, N))

    if len(ring_pairs) == 0:
        return sp.csr_matrix((N, N), dtype=np.float32)
    rows = np.concatenate([ring_pairs[:, 0], ring_pairs[:, 1]])
    cols = np.concatenate([ring_pairs[:, 1], ring_pairs[:, 0]])
    data = np.ones(len(rows), dtype=np.float32)
    return sp.csr_matrix((data, (rows, cols)), shape=(N, N))

# ---------------------------------------------------------------------------
# g_AB for given W_A, W_B and adjacency matrix
# ---------------------------------------------------------------------------
def compute_gAB(M_cpu: sp.csr_matrix, W_A: np.ndarray, W_B: np.ndarray,
                E_A: np.ndarray, E_B: np.ndarray, use_gpu: bool) -> np.ndarray:
    """Compute g_AB (nA, nB) for given marks.
    g_AB = cross / (nnz * E_A[:,None] * E_B[None,:])
    Returns numpy (nA, nB)."""
    nnz = M_cpu.nnz
    if nnz == 0:
        return np.ones((W_A.shape[1], W_B.shape[1]), dtype=np.float64)

    if use_gpu and _HAS_CUPY:
        M = csp.csr_matrix(M_cpu.astype(np.float32))
        wA = cp.asarray(W_A.astype(np.float32))
        wB = cp.asarray(W_B.astype(np.float32))
        eA = cp.asarray(E_A.astype(np.float64))
        eB = cp.asarray(E_B.astype(np.float64))
        cross = (wA.T @ (M @ wB)).astype(cp.float64)  # (nA, nB)
        denom = nnz * eA[:, None] * eB[None, :]
        g = cp.where(denom > 0, cross / denom, cp.ones_like(cross))
        return cp.asnumpy(g)
    else:
        M = M_cpu.astype(np.float32)
        wA = W_A.astype(np.float32)
        wB = W_B.astype(np.float32)
        cross = (wA.T @ (M @ wB)).astype(np.float64)
        denom = nnz * E_A[:, None] * E_B[None, :]
        return np.where(denom > 0, cross / denom, np.ones_like(cross, dtype=np.float64))

# ---------------------------------------------------------------------------
# Torus-shift of mark A on 2D bin50 grid
#
# The bin50 grid has GRID_STEP_PX=50 DNB px spacing.
# coords_px: (N, 2) float — actual pixel coords (4900, 4550, ...).
# All coords are multiples of GRID_STEP_PX.
# Shifting by (dx_gu, dy_gu) grid units wraps x on torus, reflects-clamps y.
# ---------------------------------------------------------------------------
def torus_shift_grid(W_A: np.ndarray, coords_px: np.ndarray,
                     dx_gu: int, dy_gu: int) -> np.ndarray:
    """Torus-shift of mark A by (dx_gu, dy_gu) grid units.

    x: full torus wrap (modular arithmetic on grid index).
    y: reflect-clamp — if new_y_idx falls outside [0, Ny-1], the bin gets 0.

    Strategy:
      1. Convert coords to grid-index space: idx = (coords - min) / step
      2. Apply shift: new_x_idx = (x_idx + dx_gu) % Lx; new_y_idx = y_idx + dy_gu
      3. For each dest bin i, its shifted SOURCE is the bin that was at
         (x_idx[i] - dx_gu) % Lx, y_idx[i] - dy_gu) — i.e. we pull from source.
      4. Look up source bin index in a sorted coord-key array.
      5. Dest bins where source y is out of [0, Ny-1] get W_A=0.
    """
    coords_min = coords_px.min(axis=0)
    coords_idx = np.round((coords_px - coords_min) / GRID_STEP_PX).astype(np.int64)
    xi = coords_idx[:, 0]
    yi = coords_idx[:, 1]
    Lx = int(xi.max()) + 1
    Ly = int(yi.max()) + 1

    # Source index for each dest bin: where does dest bin i pull from?
    src_xi = ((xi - dx_gu) % Lx).astype(np.int64)
    src_yi = (yi - dy_gu).astype(np.int64)
    valid_y = (src_yi >= 0) & (src_yi < Ly)

    # Build lookup: linear key = xi * Ly + yi → row index
    key_all = xi * Ly + yi                          # (N,) keys for all bins
    src_key = src_xi * Ly + src_yi                  # (N,) source keys

    sort_order = np.argsort(key_all)                # sort by key
    sorted_keys = key_all[sort_order]

    N, nA = W_A.shape
    W_shifted = np.zeros_like(W_A)

    # For each dest bin i, find pos of src_key[i] in sorted_keys
    pos = np.searchsorted(sorted_keys, src_key)
    in_range = valid_y & (pos < N)
    # Verify key actually matches (bin may not exist if tissue has gaps)
    matched = in_range & (sorted_keys[np.minimum(pos, N - 1)] == src_key)

    dest_idx = np.where(matched)[0]
    src_idx  = sort_order[pos[matched]]
    W_shifted[dest_idx] = W_A[src_idx]
    return W_shifted

# ---------------------------------------------------------------------------
# Per-chip null computation
# ---------------------------------------------------------------------------
def compute_chip_null(chip_id: str, W_A: np.ndarray, W_B: np.ndarray,
                      coords_px: np.ndarray, n_perm: int, seed: int,
                      use_gpu: bool, verbose: bool = False) -> dict:
    """Run torus-shift null for one chip. Returns per-pair Z, n_tissue."""
    N = len(coords_px)
    nA, nB = W_A.shape[1], W_B.shape[1]

    # Grid dimensions (in grid units)
    coords_min = coords_px.min(axis=0)
    coords_idx = np.round((coords_px - coords_min) / GRID_STEP_PX).astype(np.int64)
    Lx = int(coords_idx[:, 0].max()) + 1
    Ly = int(coords_idx[:, 1].max()) + 1

    # y-shift limit: ±Ly/3 grid units
    half_Ly3 = max(MIN_SHIFT_GU + 1, Ly // 3)

    # Build focal ring adjacency in pixel coords
    r_lo = RING_EDGES_PX[R_IDX]
    r_hi = RING_EDGES_PX[R_IDX + 1]
    M_cpu = build_focal_ring_cpu(coords_px, r_lo, r_hi)

    if verbose:
        print(f"  Chip {chip_id}: N={N}, grid {Lx}x{Ly}, ring nnz={M_cpu.nnz}")

    # Per-chip mark means over all in-tissue bins
    E_A = W_A.mean(axis=0).astype(np.float64)
    E_B = W_B.mean(axis=0).astype(np.float64)

    # Observed g
    g_obs = compute_gAB(M_cpu, W_A, W_B, E_A, E_B, use_gpu)  # (nA, nB)
    with np.errstate(divide='ignore', invalid='ignore'):
        log2_g_obs = np.where(g_obs > 0, np.log2(np.maximum(g_obs, 1e-30)), np.nan)

    # Welford FP64 accumulators for null
    mean_g = np.zeros((nA, nB), dtype=np.float64)
    M2_g   = np.zeros((nA, nB), dtype=np.float64)

    rng = np.random.default_rng(seed + hash(chip_id) % (2**31))

    for k in range(n_perm):
        # x-shift: random grid unit from MIN_SHIFT_GU to Lx (full torus)
        # ensure dx != 0 to avoid identity
        dx_gu = int(rng.integers(MIN_SHIFT_GU, Lx))
        # y-shift: symmetric ±[1, half_Ly3]
        dy_gu = int(rng.integers(-half_Ly3, half_Ly3 + 1))
        if dy_gu == 0:
            dy_gu = int(rng.choice([-1, 1]))

        W_A_shifted = torus_shift_grid(W_A, coords_px, dx_gu, dy_gu)

        # Recompute E_A over shifted (non-zero) bins
        E_A_shift = W_A_shifted.mean(axis=0).astype(np.float64)

        g_null = compute_gAB(M_cpu, W_A_shifted, W_B, E_A_shift, E_B, use_gpu)

        # Welford update
        n = k + 1
        delta = g_null - mean_g
        mean_g += delta / n
        M2_g  += delta * (g_null - mean_g)

        if verbose and (k + 1) % max(1, n_perm // 5) == 0:
            print(f"    perm {k+1}/{n_perm} done")

    var_g = M2_g / max(n_perm - 1, 1)
    sd_g = np.sqrt(np.maximum(var_g, 0)) + 1e-12
    Z_i = (g_obs - mean_g) / sd_g   # (nA, nB)

    _free_gpu()
    return {
        'log2_g_obs': log2_g_obs,
        'Z_i': Z_i,
        'n_in_tissue': N,
        'mu_null': mean_g,
        'sd_null': sd_g,
    }

# ---------------------------------------------------------------------------
# Load mark data for all chips
# ---------------------------------------------------------------------------
def load_marks(mode: str, logfn=print):
    """Load W_A, W_B, coords for each chip. Same mask as markcorr_runner_lib."""
    logfn("Loading meta/rctd/sct parquets...")
    meta = pd.read_parquet(F_META, columns=["bin", "chip", "x", "y", "region", "majorDomain"])
    sct_cols = ["bin", "bin_total_umi"] + PROGRAMS_60
    sct  = pd.read_parquet(F_SCT, columns=sct_cols)
    rctd = pd.read_parquet(F_RCTD, columns=["bin", "rctd_pass_mask"] + CELLTYPES)
    assert (meta.bin.values == sct.bin.values).all()
    assert (meta.bin.values == rctd.bin.values).all()

    keep = (rctd.rctd_pass_mask.values.astype(bool)
            & (meta.majorDomain.values != "ARACHNOID")
            & (sct.bin_total_umi.values >= UMI_FLOOR))
    logfn(f"Mask: {keep.sum()}/{len(keep)} bins kept")

    chip_arr = meta.chip.values
    xy = np.column_stack([meta.x.values, meta.y.values]).astype(np.float64)
    # Load the retained 54 programs by component ID rather than positional slicing.
    Wct   = np.clip(np.nan_to_num(rctd[CELLTYPES].values.astype(np.float32), nan=0.0), 0, None)
    Wprog54 = np.clip(sct[KEPT_PROGRAMS].values.astype(np.float32), 0, None)
    assert Wprog54.shape[1] == N_PROG_SPEC, f'Wprog54 shape mismatch: {Wprog54.shape}'

    chip_data = {}
    for c in pd.unique(chip_arr):
        sl = (chip_arr == c) & keep
        if sl.sum() < 50:
            continue
        if mode == 'cellprog':
            WA = Wct[sl].copy()
            WB = Wprog54[sl].copy()
        else:  # progprog
            WA = Wprog54[sl].copy()
            WB = Wprog54[sl].copy()
        chip_data[c] = (xy[sl].copy(), WA, WB)
    logfn(f"Loaded {len(chip_data)} chips for mode={mode}")
    return chip_data

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', required=True, choices=['cellprog', 'progprog'])
    parser.add_argument('--smoke', action='store_true', help='Smoke: 1 chip, n_perm=10')
    parser.add_argument('--n-perm', type=int, default=N_PERM)
    parser.add_argument('--chip-start', type=int, default=None)
    parser.add_argument('--chip-end',   type=int, default=None)
    parser.add_argument('--no-gpu', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    use_gpu = _HAS_CUPY and not args.no_gpu
    n_perm = 10 if args.smoke else args.n_perm
    os.makedirs(OUTDIR, exist_ok=True)

    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] mode={args.mode} smoke={args.smoke} "
          f"n_perm={n_perm} gpu={use_gpu} cupy={_HAS_CUPY}")

    chip_data = load_marks(args.mode)
    chips = sorted(chip_data.keys())

    if args.smoke:
        chips = chips[:1]
    elif args.chip_start is not None:
        end = args.chip_end if args.chip_end else len(chips)
        chips = chips[args.chip_start:end]

    A_names = CELLTYPES if args.mode == 'cellprog' else PROGRAMS
    B_names = PROGRAMS

    rows = []
    for ci, chip_id in enumerate(chips):
        coords_px, W_A, W_B = chip_data[chip_id]
        print(f"[{time.strftime('%H:%M:%S')}] chip {ci+1}/{len(chips)} {chip_id} "
              f"N={len(coords_px)}")
        res = compute_chip_null(chip_id, W_A, W_B, coords_px, n_perm, SEED, use_gpu,
                                verbose=args.verbose)
        Z_i    = res['Z_i']           # (nA, nB)
        log2g  = res['log2_g_obs']    # (nA, nB)
        n_tiss = res['n_in_tissue']
        mu     = res['mu_null']
        sd     = res['sd_null']

        for ai, an in enumerate(A_names):
            for bi, bn in enumerate(B_names):
                rows.append({
                    'mode':        args.mode,
                    'A_name':      an,
                    'B_name':      bn,
                    'chip_id':     chip_id,
                    'log2_g_obs':  float(log2g[ai, bi]),
                    'Z_i':         float(Z_i[ai, bi]),
                    'n_in_tissue': int(n_tiss),
                    'mu_null':     float(mu[ai, bi]),
                    'sd_null':     float(sd[ai, bi]),
                })

        elapsed = time.time() - t0
        print(f"  elapsed {elapsed:.1f}s total | ~{elapsed/(ci+1):.1f}s/chip")

    df = pd.DataFrame(rows)
    suffix = '_smoke' if args.smoke else ''
    if args.chip_start is not None:
        suffix += f'_shard{args.chip_start}_{args.chip_end or len(chips)}'
    out_f = os.path.join(OUTDIR, f"betweenchip_{args.mode}_per_chip_Z{suffix}.tsv")
    df.to_csv(out_f, sep='\t', index=False)
    print(f"Saved {len(df)} rows -> {out_f}")
    print(f"[{time.strftime('%H:%M:%S')}] DONE mode={args.mode} total_elapsed={time.time()-t0:.1f}s")

if __name__ == '__main__':
    main()
