#!/usr/bin/env python3
"""
Data-prep for the ggplot port of Figure B spatial tissue panels.
Exports small filtered long-format CSVs (per-chip, masked, per-panel vmin/vmax) so the
R patchwork script can plot geom_point fields natively (no arrow dep in R).

Faithful to scripts/figmarkcorr_spatial/build_markcorr_spatial.py:
  - SCT program score (NOT CPM)
  - mask: rctd_pass_mask & bin_total_umi>=200  -> else greyed
  - per-chip robust scaling: vmin = pct1 (if <0 else 0), vmax = pct99 over valid bins
  - aspect equal handled in R (coord_fixed); y inverted via flip in R
Outputs to results/crossregion_v1/markcorr_v2/spatial_for_R/
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from program_id_contract import load_program_contract

D   = 'CORTEX_PROGRAM_ROOT/results/crossregion_v1'
OUT = f'{D}/markcorr_v2/spatial_for_R'
os.makedirs(OUT, exist_ok=True)

CHIPS = {'B01012B2': 'DLPFC (prefrontal)', 'D00865B3': 'V1 (occipital)', 'A01186A4': 'M1 (motor)'}
CHIP_LIST = list(CHIPS.keys())
UMI_MIN = 200

P_MYELIN, P_OLIGODEV, P_PROTON, P_MICRO, P_SYN = 45, 13, 8, 40, 29
SOURCE_PROGRAMS = (P_MYELIN, P_OLIGODEV, P_PROTON, P_MICRO, P_SYN)
contract = load_program_contract(Path(D) / 'program_renumber_map.tsv')
if not set(SOURCE_PROGRAMS).issubset(contract.old_to_new):
    raise ValueError('Fig. 5 spatial exemplars must all be retained programs')

# SINGLE-SOURCE NAMING 20260615: authoritative display name = program_names.tsv
# name_short (NOT name_full). Strip any trailing '*' so the weak-FDR star is added
# ONLY dynamically by figB_large.R star_pid() (avoid double star / hard-coded star).
pn = pd.read_csv(f'{D}/program_names.tsv', sep='\t')
NAME = {int(str(p).removeprefix('P')): str(s) for p, s in zip(pn['program'], pn['name_short'])}


def retained_id(source_id):
    return contract.old_to_new[int(source_id)]


def name(source_id):
    new_id = retained_id(source_id)
    return NAME.get(new_id, f'P{new_id}').rstrip('*').rstrip()


pd.DataFrame([
    {
        'source_old_id': source_id,
        'retained_id': retained_id(source_id),
        'source_column': f'program_{source_id}',
        'display_name': name(source_id),
    }
    for source_id in SOURCE_PROGRAMS
]).to_csv(f'{OUT}/program_id_provenance.tsv', sep='\t', index=False)

print('meta...')
meta = pq.read_table(f'{D}/spatial_bin50_meta.parquet',
                     columns=['bin','chip','x','y','region']).to_pandas()
meta = meta[meta.chip.isin(CHIP_LIST)].reset_index(drop=True)
keep = set(meta['bin'])

print('rctd...')
rctd = pq.read_table(f'{D}/spatial_bin50_rctd_weights.parquet',
                     columns=['bin','rctd_pass_mask']).to_pandas()
rctd = rctd[rctd['bin'].isin(keep)].reset_index(drop=True)

print('sct...')
progs = list(SOURCE_PROGRAMS)
cols  = [f'program_{i}' for i in progs]
sct = pq.read_table(f'{D}/spatial_bin50_program_score_SCT.parquet',
                    columns=['bin','bin_total_umi']+cols).to_pandas()
sct = sct[sct['bin'].isin(keep)].reset_index(drop=True)

df = meta.merge(rctd, on='bin', how='left').merge(sct, on='bin', how='left')
df['valid'] = (df['rctd_pass_mask'].fillna(False)) & (df['bin_total_umi'].fillna(0) >= UMI_MIN)
print('valid frac per chip:\n', df.groupby('chip')['valid'].mean())

# ---- export per (chip, program) panel long table with per-panel scaling ----
def panel_table(sub, prog):
    col = f'program_{prog}'
    out = pd.DataFrame({'x': sub['x'].values, 'y': sub['y'].values,
                        'val': sub[col].values.astype(float),
                        'valid': sub['valid'].values})
    vv = out.loc[out['valid'] & np.isfinite(out['val']), 'val'].values
    if len(vv):
        vmax = np.percentile(vv, 99); lo = np.percentile(vv, 1)
        vmin = lo if lo < 0 else 0.0
        if vmax <= vmin: vmax = vmin + 1e-6
    else:
        vmin, vmax = 0.0, 1.0
    out['vmin'] = vmin; out['vmax'] = vmax
    # clip plotting value to [vmin,vmax] for invalid->NA
    out.loc[~out['valid'], 'val'] = np.nan
    return out, vmin, vmax

rows = []
# spatial1: 3 chips x [Myelin P45 | OligoDev P13 | Neurofilament P8]
# ptitle is built from authoritative name_short (no hard-coded names, no star;
# R adds weak-FDR star dynamically).
SP1 = [(P_MYELIN, f'{name(P_MYELIN)} (P{retained_id(P_MYELIN)})', 'viridis'),
       (P_OLIGODEV, f'{name(P_OLIGODEV)} (P{retained_id(P_OLIGODEV)})', 'viridis'),
       (P_PROTON, f'{name(P_PROTON)} (P{retained_id(P_PROTON)})', 'magma')]
for ci, chip in enumerate(CHIP_LIST):
    sub = df[df.chip == chip]
    for pj, (prog, ttl, cmap) in enumerate(SP1):
        t, vmin, vmax = panel_table(sub, prog)
        t['chip'] = chip; t['region_lab'] = CHIPS[chip]
        t['prog'] = retained_id(prog); t['ptitle'] = ttl; t['cmap'] = cmap
        t['row'] = ci; t['col'] = pj
        rows.append(t)
sp1 = pd.concat(rows, ignore_index=True)
# round to shrink CSV; keep only needed cols
for c in ['x','y','val','vmin','vmax']:
    sp1[c] = sp1[c].astype('float32')
sp1[['x','y','val','valid','vmin','vmax','chip','region_lab','prog','ptitle','cmap','row','col']]\
    .to_csv(f'{OUT}/sp1_long.csv', index=False)
print('sp1', sp1.shape, 'panels', sp1.groupby(['row','col']).ngroups)

# spatial2: V1 chip, P40 micro + P29 syn maps  + scatter
# ptitle from authoritative name_short (no hard-coded names, no star).
mchip = 'D00865B3'
subm = df[df.chip == mchip]
rows2 = []
for prog, ttl, cmap in [(P_MICRO, f'{name(P_MICRO)} (P{retained_id(P_MICRO)})', 'viridis'),
                        (P_SYN, f'{name(P_SYN)} (P{retained_id(P_SYN)})', 'viridis')]:
    t, vmin, vmax = panel_table(subm, prog)
    t['prog'] = retained_id(prog); t['ptitle'] = ttl; t['cmap'] = cmap
    rows2.append(t)
sp2 = pd.concat(rows2, ignore_index=True)
for c in ['x','y','val','vmin','vmax']:
    sp2[c] = sp2[c].astype('float32')
sp2[['x','y','val','valid','vmin','vmax','prog','ptitle','cmap']]\
    .to_csv(f'{OUT}/sp2_maps.csv', index=False)

# scatter: micro vs syn over valid bins
v = subm['valid'].values
a = subm[f'program_{P_MICRO}'].values[v]
b = subm[f'program_{P_SYN}'].values[v]
ok = np.isfinite(a) & np.isfinite(b)
a, b = a[ok], b[ok]
r = np.corrcoef(a, b)[0, 1]
sc = pd.DataFrame({'micro': a, 'syn': b})
# subsample scatter for vector size sanity (keep <=20k pts)
if len(sc) > 20000:
    sc = sc.sample(20000, random_state=7).reset_index(drop=True)
sc.astype('float32').to_csv(f'{OUT}/sp2_scatter.csv', index=False)
with open(f'{OUT}/sp2_meta.txt', 'w') as f:
    f.write(f'mchip\t{mchip}\nregion\t{CHIPS[mchip]}\nr\t{r:.4f}\nP_SYN\t{retained_id(P_SYN)}\nsyn_name\t{name(P_SYN)}\n')
print('sp2 maps', sp2.shape, 'scatter', sc.shape, 'r=%.3f' % r)
print('DONE ->', OUT)
