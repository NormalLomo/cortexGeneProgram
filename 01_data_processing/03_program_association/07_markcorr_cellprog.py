#!/usr/bin/env python
"""07 — mark cross-correlation: cell-type (A=22 RCTD weights) x program (B=60 z).

Per chip: ring adjacency + g_AB / Z_rho_BA + n_perm=1000 shuffle_A null
(validated markcorr_core). Aggregate per-chip -> GLOBAL and per-area (14 areas)
via equal-weight FP64 Welford (matches PROXIMA source pooling). BH + Bonferroni
across pairs; tier shortlist(padj<0.05 & |log2g|>0.14)/headline(>0.32).
Mask: exclude ARACHNOID + rctd_pass_mask + low-UMI (bin_total_umi>=100).

Outputs (results/crossregion_v1/markcorr/):
  cellprog_gr.npz (22x60x10 global g/z/log2g), cellprog_byarea.npz,
  cellprog_niches.tsv, cellprog_niches_byarea.tsv.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from markcorr_runner_lib import run, parse_runner_args

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    n_perm = 10 if smoke else 1000
    kw = parse_runner_args(sys.argv)
    run(mode="cellprog", smoke=smoke, n_perm=n_perm, **kw)
