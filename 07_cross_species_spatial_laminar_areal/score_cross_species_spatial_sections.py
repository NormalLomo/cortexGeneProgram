#!/usr/bin/env python
import os, sys, time, re
import numpy as np, pandas as pd
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BIN = 50
CHUNK = 25000000
PROJ = '/home/luomeng/DATA/cortex_nmf_program'
LOAD = f'{PROJ}/results/crossregion_v1/crossspecies/mouse_loadings.parquet'
t0 = time.time()

def parse_region_table(region_tsv):
    rt = pd.read_csv(region_tsv, sep='\t', dtype={'gene_area': 'int64'})
    out = {}
    for (ga, name) in zip(rt.gene_area, rt.area_name.astype(str)):
        if 'Isocortex' not in name:
            continue
        toks = name.split('-')
        if len(toks) < 2:
            continue
        area = toks[1]
        layer = 'NA'
        if len(toks) >= 3:
            m = re.search('(\\d+[ab]?)$', toks[2])
            if m:
                raw = m.group(1)
                layer = '2/3' if raw == '23' else raw
        out[int(ga)] = (area, layer)
    return out

def main():
    (gem, region_tsv, mouse_id, section_id, out_dir) = sys.argv[1:6]
    os.makedirs(out_dir, exist_ok=True)
    load = pd.read_parquet(LOAD)
    prog_cols = list(load.columns)
    gene_pos = {g: i for (i, g) in enumerate(load.index)}
    G = load.shape[0]
    L = load.values.astype(np.float64)
    iso = parse_region_table(region_tsv)
    iso_ids = np.array(sorted(iso.keys()), dtype=np.int64)
    iso_id_set = set((int(x) for x in iso_ids))
    BIG = np.int64(1000000)
    gene_parts = []
    meta_parts = []
    reader = pd.read_csv(gem, sep='\t', usecols=['gene', 'x', 'y', 'umi_count', 'gene_area'], dtype={'gene': 'category', 'x': 'int32', 'y': 'int32', 'umi_count': 'int32', 'gene_area': 'int64'}, chunksize=CHUNK)
    total_rows = 0
    iso_rows = 0
    for (ci, ch) in enumerate(reader):
        total_rows += len(ch)
        ga = ch.gene_area.values
        mask = np.isin(ga, iso_ids)
        if not mask.any():
            continue
        ch = ch[mask]
        iso_rows += len(ch)
        bx = (ch.x.values // BIN).astype(np.int64)
        by = (ch.y.values // BIN).astype(np.int64)
        binid = bx * BIG + by
        gidx = ch.gene.map(gene_pos).values
        umi = ch.umi_count.values.astype(np.float64)
        mdf = pd.DataFrame({'binid': binid, 'gene_area': ch.gene_area.values, 'umi': umi})
        meta_parts.append(mdf.groupby(['binid', 'gene_area'], observed=True).agg(n_dnb=('umi', 'size'), n_umi=('umi', 'sum')).reset_index())
        gm = ~pd.isna(gidx)
        gdf = pd.DataFrame({'binid': binid[gm], 'gidx': gidx[gm].astype(np.int64), 'cnt': umi[gm]})
        gene_parts.append(gdf.groupby(['binid', 'gidx'], observed=True)['cnt'].sum().reset_index())
    if not meta_parts:
        open(os.path.join(out_dir, f'{mouse_id}__{section_id}.EMPTY'), 'w').close()
        return
    meta = pd.concat(meta_parts, ignore_index=True)
    meta = meta.groupby(['binid', 'gene_area'], as_index=False).sum()
    bin_tot = meta.groupby('binid', as_index=False).agg(n_dnb=('n_dnb', 'sum'), n_umi=('n_umi', 'sum'))
    meta = meta.sort_values(['binid', 'n_dnb'], ascending=[True, False])
    top_ga = meta.drop_duplicates('binid', keep='first')[['binid', 'gene_area']]
    bin_tot = bin_tot.merge(top_ga, on='binid', how='left')
    bin_tot['area'] = bin_tot.gene_area.map(lambda g: iso[int(g)][0])
    bin_tot['layer'] = bin_tot.gene_area.map(lambda g: iso[int(g)][1])
    gene = pd.concat(gene_parts, ignore_index=True)
    gene = gene.groupby(['binid', 'gidx'], as_index=False)['cnt'].sum()
    bins = bin_tot.binid.values
    binpos = {int(b): i for (i, b) in enumerate(bins)}
    M = np.zeros((len(bins), G), dtype=np.float64)
    bi = gene.binid.map(binpos).values
    M[bi, gene.gidx.values.astype(np.int64)] = gene.cnt.values
    n_genes = (M > 0).sum(axis=1)
    scores = M @ L
    out = pd.DataFrame({'bin_x': (bins // BIG).astype(np.int32), 'bin_y': (bins % BIG).astype(np.int32), 'area': bin_tot.area.values, 'layer': bin_tot.layer.values, 'n_dnb': bin_tot.n_dnb.values.astype(np.int32), 'n_umi': bin_tot.n_umi.values.astype(np.int64), 'n_genes': n_genes.astype(np.int32)})
    sc = pd.DataFrame(scores, columns=prog_cols)
    out = pd.concat([out, sc], axis=1)
    out['mouse'] = mouse_id
    out['section'] = section_id
    out_fp = os.path.join(out_dir, f'{mouse_id}__{section_id}.parquet')
    out.to_parquet(out_fp)
if __name__ == '__main__':
    main()
