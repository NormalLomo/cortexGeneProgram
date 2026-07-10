#!/usr/bin/env python
"""08 — mark cross-correlation: program (A=60 z) x program (B=60 z), directional.

Per chip: ring adjacency + g_AB / Z_rho_BA + n_perm=1000 shuffle_A null
(validated markcorr_core). Aggregate per-chip -> GLOBAL (60x60x10 directional)
and per-area (14 areas) via equal-weight FP64 Welford (matches PROXIMA source).
BH + Bonferroni across pairs (self-diagonal excluded); tier shortlist/headline.
Mask: exclude ARACHNOID + rctd_pass_mask + low-UMI (bin_total_umi>=100).

Outputs: progprog_gr.npz, progprog_byarea.npz, progprog_niches.tsv,
progprog_niches_byarea.tsv.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from markcorr_runner_lib import run, parse_runner_args

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    n_perm = 10 if smoke else 1000
    kw = parse_runner_args(sys.argv)
    run(mode="progprog", smoke=smoke, n_perm=n_perm, **kw)
