"""Spatial mark cross-correlation core; statistic/null adapted from internal
spatial point-pattern code (PROXIMA snapshot core_v2*).

Half-open annulus ring builder on a regular lattice, normalized mark
cross-correlation g_AB(r), asymmetric anchor-on-A rho_{B|A}, and a one-sided
shuffle_A permutation null with FP64 Welford accumulators.

Marks (w_A, w_B) are continuous, non-negative per-bin weights (e.g. RCTD cell-type
weights, program scores). Distances/ring edges are in the SAME units as `coords`.

GPU via cupy when available; falls back to numpy/scipy otherwise.
"""
from typing import List, Sequence
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    _HAS_CUPY = True
except Exception:  # pragma: no cover
    cp = None
    csp = None
    _HAS_CUPY = False


# ---------------------------------------------------------------------------
# 1. Ring adjacency: half-open annuli [edge_k, edge_{k+1}) per chip
# ---------------------------------------------------------------------------
def build_ring_adjacency(coords: np.ndarray, ring_edges: Sequence[float]) -> List[sp.csr_matrix]:
    """Half-open annulus adjacency matrices on a single chip's lattice.

    Given K+1 monotone `ring_edges` = [e0, e1, ..., eK], returns K sparse CSR
    matrices. Ring k covers the half-open band [e_k, e_{k+1}) measured in the
    units of `coords`.

    M_k[i, j] = 1 iff e_k <= dist(i, j) < e_{k+1}.

    Special case: if ring 0 starts at e0 == 0, its band [0, e1) includes the
    self-distance 0. On a lattice with minimum inter-bin spacing >= e1, the only
    pairs in [0, e1) are the self-pairs, so ring 0 is the identity (within-bin
    self co-occurrence). We set it to identity explicitly. If e1 is large enough
    to also capture true neighbour pairs, those are included too (KDTree-based).
    """
    coords = np.asarray(coords, dtype=np.float64)
    N = coords.shape[0]
    edges = list(map(float, ring_edges))
    assert all(edges[i] < edges[i + 1] for i in range(len(edges) - 1)), "ring_edges must be strictly increasing"
    n_rings = len(edges) - 1

    tree = cKDTree(coords)
    # ball_pairs(r): set of (i<j) with dist <= r, as ndarray
    cache = {}

    def ball(r):
        if r not in cache:
            cache[r] = tree.query_pairs(r=r, output_type='ndarray')
        return cache[r]

    def pairset(arr):
        return set(map(tuple, arr)) if len(arr) else set()

    adjs: List[sp.csr_matrix] = []
    for k in range(n_rings):
        e_in, e_out = edges[k], edges[k + 1]
        # half-open [e_in, e_out): dist < e_out  AND  dist >= e_in
        # query_pairs uses dist <= r (closed). For half-open upper bound we take
        # pairs with dist < e_out = (pairs <= e_out) minus (pairs == e_out).
        # On a discrete lattice exact ties at e_out matter; subtract the closed
        # ball at the largest distance strictly below e_out is hard, so we use a
        # tiny epsilon below e_out for the open upper bound.
        eps = 1e-9 * max(1.0, e_out)
        outer = pairset(ball(e_out - eps))           # dist <= e_out - eps  ->  dist < e_out
        if e_in > 0:
            inner = pairset(ball(e_in - eps))         # dist < e_in
            ring_pairs = np.array(sorted(outer - inner), dtype=np.int64)
        else:
            ring_pairs = np.array(sorted(outer), dtype=np.int64)

        if k == 0 and e_in == 0.0:
            # within-bin self-pairs (self distance 0 is in [0, e1)); add identity
            self_rows = np.arange(N)
            if len(ring_pairs):
                rows = np.concatenate([ring_pairs[:, 0], ring_pairs[:, 1], self_rows])
                cols = np.concatenate([ring_pairs[:, 1], ring_pairs[:, 0], self_rows])
            else:
                rows = self_rows
                cols = self_rows
            data = np.ones(len(rows), dtype=np.float32)
            adjs.append(sp.csr_matrix((data, (rows, cols)), shape=(N, N)))
            continue

        if len(ring_pairs) == 0:
            adjs.append(sp.csr_matrix((N, N), dtype=np.float32))
            continue
        rows = np.concatenate([ring_pairs[:, 0], ring_pairs[:, 1]])
        cols = np.concatenate([ring_pairs[:, 1], ring_pairs[:, 0]])
        data = np.ones(len(rows), dtype=np.float32)
        adjs.append(sp.csr_matrix((data, (rows, cols)), shape=(N, N)))
    return adjs


# ---------------------------------------------------------------------------
# backend helpers
# ---------------------------------------------------------------------------
def _xp(use_gpu):
    return (cp, csp) if (use_gpu and _HAS_CUPY) else (np, sp)


def _to_backend_sparse(M: sp.csr_matrix, use_gpu):
    if use_gpu and _HAS_CUPY:
        return csp.csr_matrix(M.astype(np.float32))
    return M.astype(np.float32)


def _asarray(a, use_gpu):
    if use_gpu and _HAS_CUPY:
        return cp.asarray(a)
    return np.asarray(a)


def _to_numpy(a, use_gpu):
    if use_gpu and _HAS_CUPY:
        return cp.asnumpy(a)
    return np.asarray(a)


# ---------------------------------------------------------------------------
# core statistic for one permutation/observation
# ---------------------------------------------------------------------------
def _curves(adjs_b, degs, nzs, nnzs, W_A, W_B, E_A, E_B, use_gpu, collect_cross=False):
    """Compute (rho_B_given_A, g_AB) arrays of shape (n_A, n_B, n_R).

    g_AB(r)   = <w_A . w_B>_ring / (<w_A> <w_B>)               [symmetric]
    rho_{B|A} = (sum_i w_A[i] * mean_B_neighbour(i)) / (sum_i w_A[i] * E_B)  [anchor on A]

    If collect_cross=True, also returns the raw numerator tensor
    cross_all (n_A, n_B, n_R) = sum_{(i,j) in ring} wA_i wB_j, so the unratio'd
    ingredient can be dumped for offline re-aggregation. The returned (rBA, g)
    are byte-identical whether or not collect_cross is set.
    """
    xp, _ = _xp(use_gpu)
    n_A = W_A.shape[1]
    n_B = W_B.shape[1]
    n_R = len(adjs_b)
    rBA = xp.ones((n_A, n_B, n_R))
    g = xp.ones((n_A, n_B, n_R))
    cross_all = xp.zeros((n_A, n_B, n_R)) if collect_cross else None
    for r_idx in range(n_R):
        M = adjs_b[r_idx]
        nnz = nnzs[r_idx]
        if nnz == 0:
            continue
        deg = degs[r_idx]
        nz = nzs[r_idx]
        sum_B_nb = M @ W_B                         # (N, n_B)
        # mean over neighbours (degree-normalized), zero where no neighbours
        mean_B_nb = xp.where(nz[:, None], sum_B_nb / deg[:, None], xp.zeros_like(sum_B_nb))
        W_A_nz = W_A[nz]
        # rho_{B|A}: anchor on A
        num_BA = W_A_nz.T @ mean_B_nb[nz]          # (n_A, n_B)
        sum_A = W_A_nz.sum(axis=0)                 # (n_A,)
        den_BA = sum_A[:, None] * E_B[None, :]     # (n_A, n_B)
        rBA[:, :, r_idx] = xp.where(den_BA > 0, num_BA / den_BA, xp.ones_like(num_BA))
        # g_AB: <w_A w_B>_ring / (<w_A><w_B>); <.>_ring = (sum over directed pairs)/nnz
        cross = W_A.T @ sum_B_nb                    # (n_A, n_B) = sum_{(i,j) in ring} wA_i wB_j
        denom_g = nnz * E_A[:, None] * E_B[None, :]
        g[:, :, r_idx] = xp.where(denom_g > 0, cross / denom_g, xp.ones_like(cross))
        if collect_cross:
            cross_all[:, :, r_idx] = cross
    if collect_cross:
        return rBA, g, cross_all
    return rBA, g


def compute_all_pairs(
    coords: np.ndarray,
    W_A: np.ndarray,
    W_B: np.ndarray,
    ring_edges: Sequence[float] = (0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500),
    n_perm: int = 1000,
    seed: int = 42,
    use_gpu: bool = True,
    return_log2: bool = True,
    verbose: bool = False,
    return_ingredients: bool = False,
):
    """One-sided (shuffle_A) mark cross-correlation with permutation null.

    Returns dict:
      ring_edges, g_AB (n_A,n_B,n_R), log2_g_AB, rho_B_given_A,
      Z_g, Z_rho_BA, n_perm.

    Null: permute mark A only (each A column an independent randperm of bins),
    keep B fixed. Welford FP64 mean/var over permutations -> one-sided Z.
    """
    xp, _ = _xp(use_gpu)
    on_gpu = use_gpu and _HAS_CUPY
    coords = np.asarray(coords, dtype=np.float64)
    W_A = np.ascontiguousarray(W_A, dtype=np.float32)
    W_B = np.ascontiguousarray(W_B, dtype=np.float32)
    N, n_A = W_A.shape
    n_B = W_B.shape[1]
    edges = list(map(float, ring_edges))
    n_R = len(edges) - 1

    E_A_np = W_A.mean(axis=0).astype(np.float64)
    E_B_np = W_B.mean(axis=0).astype(np.float64)

    adjs_cpu = build_ring_adjacency(coords, edges)
    adjs_b = [_to_backend_sparse(M, on_gpu) for M in adjs_cpu]
    degs, nzs, nnzs = [], [], []
    for M in adjs_b:
        d = xp.asarray(np.asarray(adjs_cpu[len(degs)].sum(axis=1)).ravel(), dtype=xp.float64) if False else None
        # compute degree on backend
        dd = xp.asarray(M.sum(axis=1)).ravel()
        degs.append(dd)
        nzs.append(dd > 0)
        nnzs.append(int(M.nnz))

    W_A_b = _asarray(W_A, on_gpu)
    W_B_b = _asarray(W_B, on_gpu)
    E_A_b = _asarray(E_A_np, on_gpu)
    E_B_b = _asarray(E_B_np, on_gpu)

    if return_ingredients:
        obs_BA, obs_g, obs_cross = _curves(
            adjs_b, degs, nzs, nnzs, W_A_b, W_B_b, E_A_b, E_B_b, on_gpu,
            collect_cross=True)
    else:
        obs_BA, obs_g = _curves(adjs_b, degs, nzs, nnzs, W_A_b, W_B_b, E_A_b, E_B_b, on_gpu)

    # Welford FP64 over permutations
    mean_BA = xp.zeros((n_A, n_B, n_R), dtype=xp.float64)
    M2_BA = xp.zeros((n_A, n_B, n_R), dtype=xp.float64)
    mean_g = xp.zeros((n_A, n_B, n_R), dtype=xp.float64)
    M2_g = xp.zeros((n_A, n_B, n_R), dtype=xp.float64)

    # cupy 13.x default_rng Generator lacks .permutation; use module-level
    # permutation for the GPU path (seeded), numpy Generator for the CPU path.
    if on_gpu:
        cp.random.seed(seed)
        _perm1 = lambda: cp.random.permutation(N)
    else:
        _rng = np.random.default_rng(seed)
        _perm1 = lambda: _rng.permutation(N)
    for p in range(n_perm):
        # permute mark A only; each A column independent randperm of bins
        perm = xp.stack([_perm1() for _ in range(n_A)], axis=1)  # (N, n_A)
        W_A_p = xp.take_along_axis(W_A_b, perm, axis=0)
        rBA, g = _curves(adjs_b, degs, nzs, nnzs, W_A_p, W_B_b, E_A_b, E_B_b, on_gpu)
        k = p + 1
        d1 = rBA.astype(xp.float64) - mean_BA
        mean_BA += d1 / k
        M2_BA += d1 * (rBA.astype(xp.float64) - mean_BA)
        d3 = g.astype(xp.float64) - mean_g
        mean_g += d3 / k
        M2_g += d3 * (g.astype(xp.float64) - mean_g)
        if verbose and (p + 1) % max(1, n_perm // 10) == 0:
            print(f"  perm {p+1}/{n_perm}")

    def zscore(obs, mean, M2):
        var = M2 / max(n_perm - 1, 1)
        std = xp.sqrt(xp.clip(var, 0, None)) + 1e-10
        return (obs.astype(xp.float64) - mean) / std

    Z_g = zscore(obs_g, mean_g, M2_g)
    Z_BA = zscore(obs_BA, mean_BA, M2_BA)

    g_np = _to_numpy(obs_g, on_gpu)
    out = {
        'ring_edges': np.array(edges),
        'g_AB': g_np,
        'rho_B_given_A': _to_numpy(obs_BA, on_gpu),
        'Z_g': _to_numpy(Z_g, on_gpu),
        'Z_rho_BA': _to_numpy(Z_BA, on_gpu),
        'n_perm': n_perm,
        'null_type': 'shuffle_A',
    }
    if return_ingredients:
        # Raw observed ingredients for offline re-aggregation. Per-chip g is
        # reproduced exactly by cross / (nnz[None,None,:] * E_A[:,None,None] *
        # E_B[None,:,None]); pooled-g across chips = sum(cross) / sum(nnz*E_A*E_B).
        out['cross'] = _to_numpy(obs_cross, on_gpu)          # (n_A, n_B, n_R)
        out['nnz'] = np.array(nnzs, dtype=np.int64)          # (n_R,)
        out['E_A_vec'] = E_A_np.copy()                       # (n_A,)
        out['E_B_vec'] = E_B_np.copy()                       # (n_B,)
    if return_log2:
        with np.errstate(divide='ignore', invalid='ignore'):
            out['log2_g_AB'] = np.log2(np.where(g_np > 0, g_np, np.nan))
    return out
