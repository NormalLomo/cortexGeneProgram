#!/usr/bin/env python
"""Panel C variance partition (correct, additive decomposition).
Between-region variance of program activity decomposed into:
  - compositional: differential subclass composition across regions
  - cell-autonomous: within-subclass activity shifts across regions
via counterfactual region means holding the other factor at its global value.
"""
import os, numpy as np, pandas as pd
RES = "CORTEX_PROGRAM_ROOT/results/crossregion_v1"
df = pd.read_parquet(os.path.join(RES, "cell_program_region_subclass.parquet"))
progvar = pd.read_csv(os.path.join(RES, "program_variability.tsv"), sep="\t")
PROGS = [str(i) for i in range(1, 61)]

top_progs = progvar.sort_values("eta2_region", ascending=False).head(15)["program"].astype(str).tolist()
regions = sorted(df["region"].unique())
subclasses = sorted(df["subclass"].unique())

# region cell counts (weights for variance over regions) and global subclass fraction
n_r = df.groupby("region").size()
w_region = (n_r / n_r.sum()).reindex(regions).values            # P(region)
# global subclass fraction
n_s = df.groupby("subclass").size()
w_s_global = (n_s / n_s.sum()).reindex(subclasses).values       # P(subclass)
# composition within region: w_{r,s}
comp = (df.groupby(["region", "subclass"]).size()
          .unstack(fill_value=0))
comp = comp.div(comp.sum(axis=1), axis=0).reindex(index=regions, columns=subclasses).fillna(0).values  # r x s

def wvar(mu, w):
    mbar = np.sum(w * mu)
    return np.sum(w * (mu - mbar) ** 2)

rows = []
for p in top_progs:
    # mean of subclass s in region r: m_{r,s}
    m_rs = (df.groupby(["region", "subclass"])[p].mean()
              .unstack()).reindex(index=regions, columns=subclasses).values  # r x s, may have NaN
    # global subclass mean m_s (cell-weighted)
    m_s = df.groupby("subclass")[p].mean().reindex(subclasses).values
    # fill NaN in m_rs with global subclass mean (subclass absent in region)
    fill = np.broadcast_to(m_s, m_rs.shape)
    m_rs_f = np.where(np.isnan(m_rs), fill, m_rs)
    # observed region mean
    mu_obs = np.nansum(comp * m_rs_f, axis=1)
    # cell-autonomous: fix composition to global, vary activity
    mu_ca = np.nansum(w_s_global[None, :] * m_rs_f, axis=1)
    # compositional: fix activity to global subclass mean, vary composition
    mu_comp = np.nansum(comp * m_s[None, :], axis=1)
    V_obs = wvar(mu_obs, w_region)
    V_ca = wvar(mu_ca, w_region)
    V_comp = wvar(mu_comp, w_region)
    # rescale CA+COMP to sum to observed between-region variance (folds interaction proportionally)
    s = V_ca + V_comp
    if s <= 0:
        frac_ca, frac_comp = 0.5, 0.5
    else:
        frac_ca, frac_comp = V_ca / s, V_comp / s
    rows.append(dict(program=p,
                     V_between_region=V_obs,
                     cell_autonomous_frac=frac_ca,
                     compositional_frac=frac_comp,
                     eta2_region=float(progvar.loc[progvar["program"] == int(p), "eta2_region"].iloc[0])))
    print(f"p{p}: V_obs={V_obs:.4f} CA={frac_ca:.2f} COMP={frac_comp:.2f}", flush=True)

pd.DataFrame(rows).to_csv(os.path.join(RES, "panel_c_partition.tsv"), sep="\t", index=False)
print("DONE panel C", flush=True)
