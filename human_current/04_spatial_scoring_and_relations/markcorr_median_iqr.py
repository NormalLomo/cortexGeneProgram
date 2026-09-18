#!/usr/bin/env python
"""Final cross-chip g product: median (primary) + IQR error bar.

Companion to markcorr_aggregate_from_dumps.py. Reuses that module's loaders
(_load_dumps, _g_of via the dumped `g`) to build, for EACH (A,B,ring) element
across the 44 chips:

  median_g        primary cross-chip aggregate (== agg_<mode>_median.npz g)
  log2_median_g   log2 of median_g
  q1_g, q3_g      25th / 75th percentile of per-chip g across chips (IQR bar)
  log2_q1,log2_q3 log2 of q1_g / q3_g
  iqr_g           q3_g - q1_g
  n_chips_finite  # chips with a finite g for that element (sanity ~44)
  frac_same_sign  fraction of per-chip log2g sharing the sign of log2(median_g)
                  (direction reproducibility across chips; Z deferred)

Per-chip g is read from the dumped `g` array, sanitized exactly as the inline
aggregate / existing median aggregator does (nan/inf -> 1.0). This guarantees
median_g is byte-identical to agg_<mode>_median.npz.

Z / null-SD / significance are NOT touched here.

Outputs (under markcorr_v2/final/):
  <mode>_median_iqr.npz            all arrays + ring_edges_um + A_names/B_names
  <mode>_pairs_median_iqr.tsv      long table, one row per (A,B,ring_um)
  README.md                        column docs (written once)

Ring label convention: ring index i spans [ring_edges_um[i], ring_edges_um[i+1]);
ring_um column = OUTER edge = ring_edges_um[i+1] (so ring index 0 -> 25um,
matching the r=25um anchors).

Usage:
  python markcorr_median_iqr.py --mode cellprog
  python markcorr_median_iqr.py --mode progprog
  python markcorr_median_iqr.py --mode both
"""
import os
import argparse
import numpy as np

# reuse the existing aggregator's loaders (same dir)
from markcorr_aggregate_from_dumps import _load_dumps, OUTDIR

RESDIR = (__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "/results/crossregion_v1")
FINALDIR = os.path.join(OUTDIR, "final")
PROGRAM_NAMES_TSV = os.path.join(RESDIR, "program_names.tsv")

README = """\
# markcorr_v2/final — cross-chip g final product (median + IQR)

Primary estimator (locked): **median of per-chip g** across the 44 chips, per
(A, B, ring). Error bar (locked): **cross-chip IQR = [Q1, Q3]** (25th–75th
percentile of per-chip g). Z / significance are DEFERRED — not computed here.

Per-chip g is the dumped `g` from per_chip/<mode>_<chip>.npz, sanitized
(nan/inf -> 1.0) exactly as the inline / existing median aggregator. Therefore
`median_g` here is byte-identical to agg_<mode>_median.npz.

## Files
- `<mode>_median_iqr.npz` — arrays of shape (nA, nB, nR):
  median_g, log2_median_g, q1_g, q3_g, log2_q1, log2_q3, iqr_g,
  n_chips_finite, frac_same_sign; plus ring_edges_um, A_names, B_names, n_chips.
- `<mode>_pairs_median_iqr.tsv` — long table, one row per (A, B, ring_um).

## TSV columns
- A, B            — raw labels (cell type / program_<k> for cellprog; program_<k>
                    x program_<k> for progprog).
- A_name, B_name  — human-readable name (programs -> name_short from
                    program_names.tsv; cell types kept as-is).
- ring_um         — OUTER edge of the ring (ring_edges_um[i+1]); ring index 0 = 25um.
- median_g        — PRIMARY cross-chip aggregate (median of per-chip g).
- log2_median_g   — log2(median_g).
- q1_g, q3_g      — 25th / 75th percentile of per-chip g (the IQR error bar).
- iqr_g           — q3_g - q1_g.
- n_chips         — # chips with a finite g for that (A,B,ring) (sanity ~44).
- frac_same_sign  — fraction of per-chip log2g with the SAME sign as
                    log2(median_g); direction reproducibility across chips.

median is primary; [Q1, Q3] is the error bar; Z is deferred.
"""


def _load_program_names():
    """program_<k> -> name_short, from program_names.tsv (col 0 = int, col 2)."""
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


def compute(mode):
    paths = _load_dumps(mode)
    dumps = [np.load(p, allow_pickle=True) for p in paths]
    chips = [str(d["chip"]) for d in dumps]
    A_names = dumps[0]["A_names"]
    B_names = dumps[0]["B_names"]
    ring_edges_um = dumps[0]["ring_edges_um"].astype(int)

    # stack per-chip g, sanitized identically to the median aggregator
    gs = np.stack([np.nan_to_num(d["g"].astype(np.float64),
                                 nan=1.0, posinf=1.0, neginf=1.0)
                   for d in dumps], axis=0)            # (nChip, nA, nB, nR)
    nChip = gs.shape[0]

    median_g = np.median(gs, axis=0)
    q1_g = np.percentile(gs, 25, axis=0)
    q3_g = np.percentile(gs, 75, axis=0)
    iqr_g = q3_g - q1_g

    eps = 1e-9
    log2_median_g = np.log2(np.clip(median_g, eps, None))
    log2_q1 = np.log2(np.clip(q1_g, eps, None))
    log2_q3 = np.log2(np.clip(q3_g, eps, None))

    # n_chips_finite: count chips whose RAW dumped g is finite for that element
    raw = np.stack([d["g"].astype(np.float64) for d in dumps], axis=0)
    n_chips_finite = np.isfinite(raw).sum(axis=0).astype(np.int32)

    # frac_same_sign: fraction of per-chip log2g sharing sign with log2(median_g)
    log2gs = np.log2(np.clip(gs, eps, None))           # (nChip,nA,nB,nR)
    med_sign = np.sign(log2_median_g)[None, ...]
    chip_sign = np.sign(log2gs)
    same = (chip_sign == med_sign)
    # when median sign == 0 (g==1 exactly), define agreement as chip log2g==0
    frac_same_sign = same.mean(axis=0)

    os.makedirs(FINALDIR, exist_ok=True)
    npz_path = os.path.join(FINALDIR, f"{mode}_median_iqr.npz")
    np.savez_compressed(
        npz_path,
        median_g=median_g, log2_median_g=log2_median_g,
        q1_g=q1_g, q3_g=q3_g, log2_q1=log2_q1, log2_q3=log2_q3,
        iqr_g=iqr_g, n_chips_finite=n_chips_finite,
        frac_same_sign=frac_same_sign,
        ring_edges_um=ring_edges_um, A_names=A_names, B_names=B_names,
        n_chips=nChip, chips=np.array(chips))

    # long TSV
    prog_map = _load_program_names()
    nA, nB, nR = median_g.shape
    tsv_path = os.path.join(FINALDIR, f"{mode}_pairs_median_iqr.tsv")
    cols = ["A", "B", "A_name", "B_name", "ring_um", "median_g",
            "log2_median_g", "q1_g", "q3_g", "iqr_g", "n_chips", "frac_same_sign"]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for ia in range(nA):
            A = str(A_names[ia]); An = _name_of(A, prog_map)
            for ib in range(nB):
                B = str(B_names[ib]); Bn = _name_of(B, prog_map)
                for ir in range(nR):
                    ring_um = int(ring_edges_um[ir + 1])  # outer edge
                    fh.write("\t".join([
                        A, B, An, Bn, str(ring_um),
                        f"{median_g[ia,ib,ir]:.6g}",
                        f"{log2_median_g[ia,ib,ir]:.6g}",
                        f"{q1_g[ia,ib,ir]:.6g}",
                        f"{q3_g[ia,ib,ir]:.6g}",
                        f"{iqr_g[ia,ib,ir]:.6g}",
                        str(int(n_chips_finite[ia,ib,ir])),
                        f"{frac_same_sign[ia,ib,ir]:.4f}",
                    ]) + "\n")

    print(f"[{mode}] {nChip} chips -> {npz_path}")
    print(f"[{mode}] long table -> {tsv_path}  rows={nA*nB*nR}")
    print(f"[{mode}] median_g min/max={median_g.min():.4f}/{median_g.max():.4f} "
          f"n_chips_finite min/max={n_chips_finite.min()}/{n_chips_finite.max()}")
    return npz_path, tsv_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["cellprog", "progprog", "both"])
    a = ap.parse_args()
    modes = ("cellprog", "progprog") if a.mode == "both" else (a.mode,)
    os.makedirs(FINALDIR, exist_ok=True)
    with open(os.path.join(FINALDIR, "README.md"), "w") as fh:
        fh.write(README)
    for m in modes:
        compute(m)


if __name__ == "__main__":
    main()
