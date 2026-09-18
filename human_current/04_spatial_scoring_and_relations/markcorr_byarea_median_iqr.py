#!/usr/bin/env python
"""Per-AREA cross-chip g product: median (primary) + IQR error bar, stratified
by cortical region (the `area` field in each per_chip dump).

Companion to markcorr_median_iqr.py (global, median over ALL 44 chips) and
markcorr_aggregate_from_dumps.py. Reuses the SAME loaders/logic
(_load_dumps, dumped `g` sanitized nan/inf->1.0), but GROUPS chips by `area`
and computes, for EACH area and EACH (A,B,ring) element across that area's chips:

  median_g        primary per-area cross-chip aggregate (median of per-chip g)
  log2_median_g   log2 of median_g
  q1_g, q3_g      25th / 75th percentile of per-chip g within the area (IQR bar)
  iqr_g           q3_g - q1_g
  frac_same_sign  fraction of that area's per-chip log2g sharing sign of log2(median_g)
  n_chips_area    # chips in that area (sanity)
  unstable        1 if n_chips_area < 3 (IQR/median unstable), else 0

Sanitization identical to the global median aggregator (guarantees consistency).
Z / null-SD / significance NOT touched.

Outputs (under markcorr_v2/final/) — ADD-ONLY, global files untouched:
  <mode>_byarea_median_iqr.npz   arrays shaped (nArea, nA, nB, nR):
        median_g, log2_median_g, q1_g, q3_g, iqr_g, frac_same_sign
      + n_chips_area (nArea,), unstable (nArea,), area_names (nArea,),
        A_names, B_names, ring_edges_um.
  <mode>_byarea_median_iqr.tsv   long table (one row per area,A,B,ring_um),
        RESTRICTED to headline pairs (|global log2_median_g| > 0.32 at ANY ring)
        to keep size sane; FULL tensor lives in the npz.

Ring label convention: ring index i -> outer edge ring_edges_um[i+1]
(ring 0 -> 25um), matching the global product.

Usage:
  python markcorr_byarea_median_iqr.py --mode cellprog
  python markcorr_byarea_median_iqr.py --mode progprog
  python markcorr_byarea_median_iqr.py --mode both
"""
import os
import argparse
import numpy as np

from markcorr_aggregate_from_dumps import _load_dumps, OUTDIR

RESDIR = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")
FINALDIR = os.path.join(OUTDIR, "final")
PROGRAM_NAMES_TSV = os.path.join(RESDIR, "program_names.tsv")
HEADLINE_THR = 0.32  # |global log2_median_g| threshold for TSV headline pairs
MIN_CHIPS_STABLE = 3


def _load_program_names():
    m = {}
    if not os.path.exists(PROGRAM_NAMES_TSV):
        return m
    with open(PROGRAM_NAMES_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci = header.index("name_short")
        except ValueError:
            ci = 2
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0].strip():
                continue
            m[f"program_{parts[0].strip()}"] = parts[ci] if len(parts) > ci else parts[0]
    return m


def _name_of(label, prog_map):
    return prog_map.get(str(label), str(label))


def _headline_mask(mode, A_names, B_names):
    """Boolean (nA,nB): True if global |log2_median_g| > thr at ANY ring.
    Falls back to all-True if the global npz is missing."""
    gpath = os.path.join(FINALDIR, f"{mode}_median_iqr.npz")
    nA, nB = len(A_names), len(B_names)
    if not os.path.exists(gpath):
        print(f"[{mode}] WARN: global {gpath} missing -> TSV keeps ALL pairs")
        return np.ones((nA, nB), dtype=bool)
    g = np.load(gpath, allow_pickle=True)
    l2 = g["log2_median_g"]                      # (nA,nB,nR)
    return (np.abs(l2) > HEADLINE_THR).any(axis=2)


def compute(mode):
    paths = _load_dumps(mode)
    dumps = [np.load(p, allow_pickle=True) for p in paths]
    chips = [str(d["chip"]) for d in dumps]
    areas = [str(d["area"]) for d in dumps]
    A_names = dumps[0]["A_names"]
    B_names = dumps[0]["B_names"]
    ring_edges_um = dumps[0]["ring_edges_um"].astype(int)

    # per-chip g, sanitized identically to the global median aggregator
    gs = np.stack([np.nan_to_num(d["g"].astype(np.float64),
                                 nan=1.0, posinf=1.0, neginf=1.0)
                   for d in dumps], axis=0)            # (nChip,nA,nB,nR)
    eps = 1e-9
    log2gs = np.log2(np.clip(gs, eps, None))

    area_names = sorted(set(areas))
    nArea = len(area_names)
    nA, nB, nR = gs.shape[1:]

    median_g = np.empty((nArea, nA, nB, nR), dtype=np.float64)
    q1_g = np.empty_like(median_g)
    q3_g = np.empty_like(median_g)
    frac_same_sign = np.empty_like(median_g)
    n_chips_area = np.empty(nArea, dtype=np.int32)
    unstable = np.empty(nArea, dtype=np.int32)

    areas_arr = np.array(areas)
    for ai, area in enumerate(area_names):
        idx = np.where(areas_arr == area)[0]
        n_chips_area[ai] = len(idx)
        unstable[ai] = 1 if len(idx) < MIN_CHIPS_STABLE else 0
        sub = gs[idx]                                  # (nSub,nA,nB,nR)
        med = np.median(sub, axis=0)
        median_g[ai] = med
        q1_g[ai] = np.percentile(sub, 25, axis=0)
        q3_g[ai] = np.percentile(sub, 75, axis=0)
        med_sign = np.sign(np.log2(np.clip(med, eps, None)))[None, ...]
        frac_same_sign[ai] = (np.sign(log2gs[idx]) == med_sign).mean(axis=0)

    log2_median_g = np.log2(np.clip(median_g, eps, None))
    iqr_g = q3_g - q1_g

    os.makedirs(FINALDIR, exist_ok=True)
    npz_path = os.path.join(FINALDIR, f"{mode}_byarea_median_iqr.npz")
    np.savez_compressed(
        npz_path,
        median_g=median_g, log2_median_g=log2_median_g,
        q1_g=q1_g, q3_g=q3_g, iqr_g=iqr_g,
        frac_same_sign=frac_same_sign,
        n_chips_area=n_chips_area, unstable=unstable,
        area_names=np.array(area_names),
        A_names=A_names, B_names=B_names, ring_edges_um=ring_edges_um)

    # long TSV restricted to headline pairs (global |log2_median_g|>thr at any ring)
    prog_map = _load_program_names()
    hmask = _headline_mask(mode, A_names, B_names)
    n_pairs_kept = int(hmask.sum())
    tsv_path = os.path.join(FINALDIR, f"{mode}_byarea_median_iqr.tsv")
    cols = ["area", "n_chips_area", "unstable", "A", "B", "A_name", "B_name",
            "ring_um", "median_g", "log2_median_g", "q1_g", "q3_g", "iqr_g",
            "frac_same_sign"]
    nrow = 0
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for ai, area in enumerate(area_names):
            nc = int(n_chips_area[ai]); uns = int(unstable[ai])
            for ia in range(nA):
                A = str(A_names[ia]); An = _name_of(A, prog_map)
                for ib in range(nB):
                    if not hmask[ia, ib]:
                        continue
                    B = str(B_names[ib]); Bn = _name_of(B, prog_map)
                    for ir in range(nR):
                        ring_um = int(ring_edges_um[ir + 1])
                        fh.write("\t".join([
                            area, str(nc), str(uns), A, B, An, Bn, str(ring_um),
                            f"{median_g[ai,ia,ib,ir]:.6g}",
                            f"{log2_median_g[ai,ia,ib,ir]:.6g}",
                            f"{q1_g[ai,ia,ib,ir]:.6g}",
                            f"{q3_g[ai,ia,ib,ir]:.6g}",
                            f"{iqr_g[ai,ia,ib,ir]:.6g}",
                            f"{frac_same_sign[ai,ia,ib,ir]:.4f}",
                        ]) + "\n")
                        nrow += 1

    print(f"[{mode}] {len(dumps)} chips, {nArea} areas -> {npz_path}")
    print(f"[{mode}] areas: " + ", ".join(
        f"{a}={n}{'*' if u else ''}"
        for a, n, u in zip(area_names, n_chips_area, unstable)))
    print(f"[{mode}] headline pairs kept={n_pairs_kept}/{nA*nB} "
          f"-> TSV {tsv_path} rows={nrow}")
    print(f"[{mode}] median_g min/max={median_g.min():.4f}/{median_g.max():.4f}")
    return npz_path, tsv_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["cellprog", "progprog", "both"])
    a = ap.parse_args()
    modes = ("cellprog", "progprog") if a.mode == "both" else (a.mode,)
    os.makedirs(FINALDIR, exist_ok=True)
    for m in modes:
        compute(m)


if __name__ == "__main__":
    main()
