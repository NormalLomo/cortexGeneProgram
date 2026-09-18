"""Shared cross-chip aggregation for mark cross-correlation runners (07/08).

AGGREGATION LOGIC — matched to source PROXIMA multi-chip drivers
(proxima_cross_tissue_runner.py / proxima_v1_1.proxima_inhom_multichip):
  Every source driver pools chips by EQUAL-WEIGHT per-chip averaging
  (`Z_sum += z; Z_count += 1; out = Z_sum/Z_count`) — each chip = ONE
  observation, NOT pair-count weighted. The per-bin means <w_A>,<w_B>
  (E_A,E_B) are PER-CHIP (recomputed inside each chip via W.mean(axis=0)),
  never pooled across chips. We match this exactly: FP64 Welford over the
  per-chip g_AB (and Z_rho_BA) tensors, equal weight, per-chip means.
"""
import numpy as np


class Welford:
    """FP64 Welford accumulator over equal-weight per-chip observation tensors."""
    def __init__(self, shape):
        self.n = 0
        self.mean = np.zeros(shape, dtype=np.float64)
        self.M2 = np.zeros(shape, dtype=np.float64)

    def update(self, x):
        x = np.asarray(x, dtype=np.float64)
        self.n += 1
        d = x - self.mean
        self.mean += d / self.n
        self.M2 += d * (x - self.mean)

    def finalize(self):
        mean = self.mean.copy()
        if self.n > 1:
            var = self.M2 / (self.n - 1)
        else:
            var = np.zeros_like(self.M2)
        sd = np.sqrt(np.clip(var, 0, None))
        return mean, sd, self.n
