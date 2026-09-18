"""Engine shared by 07_markcorr_cellprog.py and 08_markcorr_progprog.py.

Per chip: build ring adjacency, compute g_AB + Z_rho_BA via the VALIDATED
markcorr_core.compute_all_pairs (one-sided shuffle_A FP64-Welford null).
Aggregate per-chip -> GLOBAL and per-area (14 areas) using equal-weight FP64
Welford over per-chip tensors (matches PROXIMA source; see markcorr_aggregate).
Then BH + Bonferroni across pairs at the peak ring; tier shortlist/headline.

Coords are DNB px; ring edges in um -> px = um * 2.
"""
import os, sys, time, traceback
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from markcorr_core import compute_all_pairs
from markcorr_aggregate import Welford

try:  # optional: free cupy pools between chips in a long sharded run
    import cupy as _cp
except Exception:  # pragma: no cover
    _cp = None


def _free_gpu_pools():
    """Release cupy device + pinned memory pools so a shard processing many
    chips serially does not accumulate GPU memory across chips. No-op on CPU."""
    if _cp is not None:
        try:
            _cp.get_default_memory_pool().free_all_blocks()
            _cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass

RESDIR = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")
# v2: identical METHOD/outputs as v1 markcorr/, plus per_chip/ ingredient dumps.
# v1 markcorr/ is left untouched.
OUTDIR = os.path.join(RESDIR, "markcorr_v2")
PERCHIP_DIR = os.path.join(OUTDIR, "per_chip")
F_META = os.path.join(RESDIR, "spatial_bin50_meta.parquet")
F_SCT  = os.path.join(RESDIR, "spatial_bin50_program_score_SCT.parquet")
F_RCTD = os.path.join(RESDIR, "spatial_bin50_rctd_weights.parquet")

RING_EDGES_UM = (0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500)
RING_EDGES_PX = tuple(e * 2 for e in RING_EDGES_UM)  # DNB px = um*2
N_RINGS = len(RING_EDGES_UM) - 1

CELLTYPES = ['AST','CHANDELIER','ENDO','ET','L2-L3 IT LINC00507','L3-L4 IT RORB',
             'L4-L5 IT RORB','L6 CAR3','L6 CT','L6 IT','L6B','LAMP5','MICRO','NDNF',
             'NP','OLIGO','OPC','PAX6','PVALB','SST','VIP','VLMC']  # 22
PROGRAMS = [f'program_{i}' for i in range(1, 61)]  # 60

UMI_FLOOR = 100  # low-UMI mask floor (1% pct=58, deep-trunc region)
TIER_SHORT = 0.14   # |log2 g| shortlist
TIER_HEAD  = 0.32   # |log2 g| headline


def log(msg, logf):
    logf.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n"); logf.flush()


_PERCHIP_README = """per_chip/ — per-chip mark cross-correlation ingredient dumps (v2)
================================================================
One compressed .npz per (mode, chip): <mode>_<chipID>.npz
  mode in {cellprog (A=22 cell-types, B=60 programs),
           progprog (A=B=60 programs)}.

Purpose: the production cross-chip aggregate (markcorr_v2/<mode>_gr.npz) uses
EQUAL-WEIGHT mean of per-chip-normalized g (byte-identical method to v1
markcorr/). These dumps store the RAW unratio'd ingredients so ANY other
aggregation (pooled, median, trimmed, weighted) can be recomputed OFFLINE
without re-running the 20h pipeline.

Arrays per file:
  cross        (nA, nB, nR) float64  Σ_{(i,j) in ring} wA_i · wB_j  (numerator)
  nnz          (nR,)        int64    directed pair count per ring (denominator count)
  E_A_vec      (nA,)        float64  per-CHIP per-mark mean of A (<w_A> over masked bins)
  E_B_vec      (nB,)        float64  per-CHIP per-mark mean of B (<w_B> over masked bins)
  g            (nA, nB, nR) float64  per-chip g this chip fed into the aggregate
  rho_BA       (nA, nB, nR) float64  per-chip rho_{B|A} (anchor-on-A) fed into Z aggregate
  ring_edges_um(nR+1,)               ring edges in microns
  A_names, B_names                   mark labels
  chip, area, N_bins, mode           metadata (0-d arrays)

Reconstruction identities (let d = nR ring index):
  per-chip g[a,b,r]  ==  cross[a,b,r] / (nnz[r] * E_A_vec[a] * E_B_vec[b])
                        (where nnz[r]*E_A*E_B > 0; else g==1 by convention)

  EQUAL-WEIGHT aggregate (== v1/v2 production):
      g_eq[a,b,r] = mean_chip( g[a,b,r] )   # nan/inf -> 1 before averaging
  POOLED aggregate:
      g_pool[a,b,r] = sum_chip(cross[a,b,r]) / sum_chip(nnz[r]*E_A[a]*E_B[b])
  MEDIAN aggregate:
      g_med[a,b,r] = median_chip( g[a,b,r] )

NOTE: Z / null-SD logic is UNCHANGED from v1; these dumps do NOT carry null
draws — they are observed-statistic ingredients only.
"""


def _write_perchip_readme():
    p = os.path.join(PERCHIP_DIR, "README.txt")
    with open(p, "w") as f:
        f.write(_PERCHIP_README)


def load_marks(b_is_program, logf):
    """Return dict chip -> (coords px, W_A, W_B, area) on masked bins.
    A,B per mode: cellprog A=22 RCTD, B=60 prog. progprog A=B=60 prog."""
    log("loading meta/rctd/sct ...", logf)
    meta = pd.read_parquet(F_META, columns=["bin", "chip", "x", "y", "region", "majorDomain"])
    sct_cols = ["bin", "bin_total_umi"] + PROGRAMS
    sct = pd.read_parquet(F_SCT, columns=sct_cols)
    rctd = pd.read_parquet(F_RCTD, columns=["bin", "rctd_pass_mask"] + CELLTYPES)
    assert (meta.bin.values == sct.bin.values).all() and (meta.bin.values == rctd.bin.values).all()

    keep = (rctd.rctd_pass_mask.values.astype(bool)
            & (meta.majorDomain.values != "ARACHNOID")
            & (sct.bin_total_umi.values >= UMI_FLOOR))
    log(f"mask keeps {keep.sum()}/{len(keep)} bins ({100*keep.mean():.1f}%)", logf)

    chip = meta.chip.values
    region = meta.region.values
    xy = np.column_stack([meta.x.values, meta.y.values]).astype(np.float64)
    Wprog = sct[PROGRAMS].values.astype(np.float32)
    Wct = rctd[CELLTYPES].values.astype(np.float32)
    # clip negatives to 0 (marks must be non-negative; SCT z can be <0)
    Wprog = np.clip(Wprog, 0, None)
    Wct = np.clip(np.nan_to_num(Wct, nan=0.0), 0, None)

    out = {}
    for c in pd.unique(chip):
        sl = (chip == c) & keep
        if sl.sum() < 50:
            continue
        WA = Wct[sl] if not b_is_program else Wprog[sl]
        WB = Wprog[sl]
        # area = modal region of the chip
        rr = region[sl]
        area = pd.Series(rr).mode().iloc[0]
        out[c] = (xy[sl], WA.copy(), WB.copy(), area)
    log(f"{len(out)} chips usable", logf)
    return out


def parse_runner_args(argv):
    """Shared CLI parse for 07/08 (Piece 1). Returns kwargs dict for run().
      --chips N               first-N usable chips (chip_limit)
      --shard i --nshards K   process chips where index%K==i (disjoint subset)
      --chip-start a [--chip-end b]   process chips[a:b]
    --shard/--nshards and --chip-start/--chip-end are mutually exclusive."""
    def _val(flag, cast=int):
        if flag in argv:
            return cast(argv[argv.index(flag) + 1])
        return None
    kw = dict(
        chip_limit=_val("--chips"),
        shard=_val("--shard"),
        nshards=_val("--nshards"),
        chip_start=_val("--chip-start"),
        chip_end=_val("--chip-end"),
    )
    if (kw["shard"] is not None) != (kw["nshards"] is not None):
        raise SystemExit("--shard and --nshards must be given together")
    if kw["shard"] is not None and (kw["chip_start"] is not None
                                    or kw["chip_end"] is not None):
        raise SystemExit("use EITHER --shard/--nshards OR --chip-start/--chip-end")
    return kw


def run(mode, n_perm, chip_limit=None, shard=None, nshards=None,
        chip_start=None, chip_end=None):
    """mode: 'cellprog' (A=22ct,B=60prog) | 'progprog' (A=B=60prog).

    chip_limit: if set (int), only process the first N usable chips. Non-invasive
    subset limiter for smoke/dev; None = all chips (production).

    SHARDING (Piece 1; for K-way concurrent execution on the single GPU). Two
    mutually-exclusive ways to select a DISJOINT subset of the usable chips so K
    processes can run concurrently, each writing its own per_chip/*.npz:
      shard=i, nshards=K  -> this process handles chips where (index % K == i).
      chip_start=a, chip_end=b -> this process handles chips[a:b] (half-open).
    Selection order: chip_limit (first-N) is applied FIRST, then shard/range.
    When EITHER sharding mode is active, NO single process sees all chips, so the
    cross-chip aggregate (_gr/_byarea/_niches) is NOT written here — the per_chip
    dumps are the source of truth, aggregated offline by
    markcorr_aggregate_from_dumps.py (Piece 2). Default (both None) = unchanged
    v1/v2 behavior: process all chips AND write the in-line equal-weight aggregate.
    """
    sharding = (shard is not None and nshards is not None) or \
               (chip_start is not None or chip_end is not None)
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PERCHIP_DIR, exist_ok=True)
    _write_perchip_readme()
    done_f = os.path.join(OUTDIR, "_DONE")
    fail_f = os.path.join(OUTDIR, "_FAILED")
    prefix = ""
    logname = f"{prefix}{mode}.runlog"
    logf = open(os.path.join(OUTDIR, logname), "a")
    try:
        b_is_program = True  # B always programs
        a_is_program = (mode == "progprog")
        nA = 60 if a_is_program else 22
        nB = 60
        A_names = PROGRAMS if a_is_program else CELLTYPES
        B_names = PROGRAMS

        chips = load_marks(b_is_program=a_is_program, logf=logf)
        chip_ids = list(chips.keys())
        if chip_limit is not None:
            chip_ids = chip_ids[:int(chip_limit)]
        # ---- shard / range selection (Piece 1) ----
        all_ids = list(chip_ids)  # post-limit ordering = the canonical shard index space
        if shard is not None and nshards is not None:
            shard, nshards = int(shard), int(nshards)
            assert 0 <= shard < nshards, f"need 0<=shard<nshards, got {shard}/{nshards}"
            chip_ids = [c for j, c in enumerate(all_ids) if j % nshards == shard]
            log(f"SHARD {shard}/{nshards}: {len(chip_ids)}/{len(all_ids)} chips "
                f"-> {chip_ids}", logf)
        elif chip_start is not None or chip_end is not None:
            a = 0 if chip_start is None else int(chip_start)
            b = len(all_ids) if chip_end is None else int(chip_end)
            chip_ids = all_ids[a:b]
            log(f"RANGE [{a}:{b}): {len(chip_ids)}/{len(all_ids)} chips "
                f"-> {chip_ids}", logf)
        log(f"mode={mode} n_perm={n_perm} chips={len(chip_ids)} "
            f"sharding={sharding} shape={nA}x{nB}x{N_RINGS}", logf)

        glob_g = Welford((nA, nB, N_RINGS))
        glob_z = Welford((nA, nB, N_RINGS))
        area_g, area_z = {}, {}
        area_of = {}

        for i, c in enumerate(chip_ids):
            coords, WA, WB, area = chips[c]
            t0 = time.time()
            dump_path = os.path.join(PERCHIP_DIR, f"{prefix}{mode}_{c}.npz")
            # ---- RESUMABLE skip-existing (non-invasive): if this chip's per_chip
            # dump already exists, DO NOT recompute. Load g/rho_BA from the dump
            # and still feed the in-line aggregate (so non-sharded default stays
            # correct on resume). When no dump exists, behavior is unchanged.
            if os.path.exists(dump_path):
                d = np.load(dump_path, allow_pickle=True)
                g = d["g"].astype(np.float64)
                # reconstruct Z_rho_BA is not stored; for aggregate we need z.
                # The dumps carry rho_BA (observed) but NOT null SD, so Z cannot
                # be rebuilt from a dump. For the SHARDED production run the
                # in-line aggregate is skipped anyway (sharding=True), so z is
                # unused. For the (rare) non-sharded resume we recompute z only.
                if not sharding:
                    res = compute_all_pairs(
                        coords, WA, WB, ring_edges=RING_EDGES_PX,
                        n_perm=n_perm, seed=42, use_gpu=True, return_log2=False,
                        verbose=False, return_ingredients=False)
                    z = res["Z_rho_BA"]
                else:
                    z = None
                g = np.nan_to_num(g, nan=1.0, posinf=1.0, neginf=1.0)
                if not sharding:
                    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
                    glob_g.update(g); glob_z.update(z)
                    if area not in area_g:
                        area_g[area] = Welford((nA, nB, N_RINGS)); area_z[area] = Welford((nA, nB, N_RINGS))
                    area_g[area].update(g); area_z[area].update(z)
                    area_of[c] = area
                log(f"chip {i+1}/{len(chip_ids)} {c} SKIP (dump exists) "
                    f"{time.time()-t0:.1f}s", logf)
                _free_gpu_pools()
                continue
            res = compute_all_pairs(
                coords, WA, WB, ring_edges=RING_EDGES_PX,
                n_perm=n_perm, seed=42, use_gpu=True, return_log2=False, verbose=False,
                return_ingredients=True)
            g = res["g_AB"]; z = res["Z_rho_BA"]
            # ---- per-chip ingredient dump (ADD-ONLY; does not affect aggregation) ----
            # Dump RAW values (pre nan-sanitize) so any aggregation is reconstructable.
            np.savez_compressed(
                dump_path,
                cross=res["cross"].astype(np.float64),        # (nA,nB,nR) Σ wA·wB over ring pairs
                nnz=res["nnz"].astype(np.int64),              # (nR,) directed pair count per ring
                E_A_vec=res["E_A_vec"].astype(np.float64),    # (nA,) per-chip per-mark mean of A
                E_B_vec=res["E_B_vec"].astype(np.float64),    # (nB,) per-chip per-mark mean of B
                g=g.astype(np.float64),                       # (nA,nB,nR) per-chip g produced by pipeline
                rho_BA=res["rho_B_given_A"].astype(np.float64),  # (nA,nB,nR) per-chip rho_{B|A}
                ring_edges_um=np.array(RING_EDGES_UM),
                A_names=np.array(A_names), B_names=np.array(B_names),
                chip=np.array(c), area=np.array(area), N_bins=np.array(len(coords)),
                mode=np.array(mode))
            g = np.nan_to_num(g, nan=1.0, posinf=1.0, neginf=1.0)
            z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
            glob_g.update(g); glob_z.update(z)
            if area not in area_g:
                area_g[area] = Welford((nA, nB, N_RINGS)); area_z[area] = Welford((nA, nB, N_RINGS))
            area_g[area].update(g); area_z[area].update(z)
            area_of[c] = area
            log(f"chip {i+1}/{len(chip_ids)} {c} area={area} N={len(coords)} "
                f"{time.time()-t0:.1f}s g[min/max]={g.min():.3f}/{g.max():.3f}", logf)
            # ---- free cupy pools so memory does not accumulate across chips ----
            _free_gpu_pools()

        if sharding:
            # This process saw only a disjoint subset of chips; the cross-chip
            # aggregate is meaningless here. per_chip dumps already written above
            # are the source of truth -> aggregate offline with
            # markcorr_aggregate_from_dumps.py. Do NOT write _gr/_byarea/_niches.
            log(f"SHARD/RANGE done {mode}: dumped {len(chip_ids)} per_chip npz; "
                f"skipping cross-chip aggregate (use markcorr_aggregate_from_dumps.py)",
                logf)
            logf.close()
            return True

        # ---- GLOBAL finalize ----
        g_mean, g_sd, g_n = glob_g.finalize()
        z_mean, z_sd, _ = glob_z.finalize()
        log2g = np.log2(np.clip(g_mean, 1e-9, None))
        np.savez_compressed(
            os.path.join(OUTDIR, f"{prefix}{mode}_gr.npz"),
            ring_edges_um=np.array(RING_EDGES_UM), g=g_mean, g_sd=g_sd,
            z_rho_BA=z_mean, log2g=log2g, n_chips=g_n,
            A_names=np.array(A_names), B_names=np.array(B_names))

        # ---- per-area finalize ----
        area_list = sorted(area_g.keys())
        area_g_arr = np.stack([area_g[a].finalize()[0] for a in area_list])     # (nArea,nA,nB,nR)
        area_z_arr = np.stack([area_z[a].finalize()[0] for a in area_list])
        area_n = np.array([area_g[a].n for a in area_list])
        np.savez_compressed(
            os.path.join(OUTDIR, f"{prefix}{mode}_byarea.npz"),
            ring_edges_um=np.array(RING_EDGES_UM), areas=np.array(area_list),
            area_n_chips=area_n, g=area_g_arr,
            z_rho_BA=area_z_arr, log2g=np.log2(np.clip(area_g_arr, 1e-9, None)),
            A_names=np.array(A_names), B_names=np.array(B_names))

        # ---- niche tables (GLOBAL): per (pair,ring) p from |Z|, two-sided, BH+Bonf ----
        rows = []
        for ia in range(nA):
            for ib in range(nB):
                if a_is_program and ia == ib:
                    continue  # skip self program-program diagonal
                for ir in range(N_RINGS):
                    z = z_mean[ia, ib, ir]
                    l2 = log2g[ia, ib, ir]
                    p = 2 * norm.sf(abs(z))
                    rows.append((A_names[ia], B_names[ib], RING_EDGES_UM[ir + 1],
                                 l2, z, p))
        df = pd.DataFrame(rows, columns=["A", "B", "ring_um", "log2g", "Z_rho_BA", "p"])
        df["padj_bh"] = multipletests(df.p.values, method="fdr_bh")[1]
        df["padj_bonf"] = np.clip(df.p.values * len(df), 0, 1)
        df["tier"] = "none"
        sh = (df.padj_bh < 0.05) & (df.log2g.abs() > TIER_SHORT)
        df.loc[sh, "tier"] = "shortlist"
        hd = (df.padj_bh < 0.05) & (df.log2g.abs() > TIER_HEAD)
        df.loc[hd, "tier"] = "headline"
        df = df.rename(columns={"A": "A", "B": "B"})
        df.to_csv(os.path.join(OUTDIR, f"{prefix}{mode}_niches.tsv"),
                  sep="\t", index=False)

        # ---- per-area niche table ----
        arows = []
        for ai, a in enumerate(area_list):
            zA = area_z_arr[ai]; l2A = np.log2(np.clip(area_g_arr[ai], 1e-9, None))
            for ia in range(nA):
                for ib in range(nB):
                    if a_is_program and ia == ib:
                        continue
                    for ir in range(N_RINGS):
                        z = zA[ia, ib, ir]
                        arows.append((a, A_names[ia], B_names[ib], RING_EDGES_UM[ir + 1],
                                      l2A[ia, ib, ir], z, 2 * norm.sf(abs(z))))
        adf = pd.DataFrame(arows, columns=["area", "A", "B", "ring_um", "log2g", "Z_rho_BA", "p"])
        adf["padj_bh"] = multipletests(adf.p.values, method="fdr_bh")[1]
        adf["padj_bonf"] = np.clip(adf.p.values * len(adf), 0, 1)
        adf["tier"] = "none"
        sh = (adf.padj_bh < 0.05) & (adf.log2g.abs() > TIER_SHORT)
        adf.loc[sh, "tier"] = "shortlist"
        hd = (adf.padj_bh < 0.05) & (adf.log2g.abs() > TIER_HEAD)
        adf.loc[hd, "tier"] = "headline"
        adf.to_csv(os.path.join(OUTDIR, f"{prefix}{mode}_niches_byarea.tsv"),
                   sep="\t", index=False)

        log(f"DONE {mode}: global g {g_mean.shape}, areas={len(area_list)}, "
            f"niches={len(df)}, headline={int((df.tier=='headline').sum())}", logf)
        logf.close()
        return True
    except Exception:
        tb = traceback.format_exc()
        log("FAILED:\n" + tb, logf)
        with open(os.path.join(OUTDIR, "_FAILED"), "w") as ff:
            ff.write(f"{mode}\n{tb}\n")
        logf.close()
        raise
