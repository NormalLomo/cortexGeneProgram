#!/usr/bin/env python
import os, sys, time, json
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
BASE = str(PROJECT_ROOT)
NP = '/home/luomeng/DATA/neuropeptide_cortex'
SPECTRA = f'{BASE}/results/cnmf_snrna_joint_full1M_v1/cnmf_work/snrna_joint_full1M_v1/snrna_joint_full1M_v1.gene_spectra_score.k_60.dt_0_15.txt'
RETAIN_MAP = f'{BASE}/results/crossregion_v1/program_renumber_map.tsv'
ORTHO = f'{BASE}/data/ortholog/ortholog_human_mouse_cynomolgus_compara91.tsv'
OUT = f'{BASE}/results/xspecies_humanmap_v1'
H5 = {'mouse': f'{NP}/data/mouse/snrna/snRNA_mouse.h5ad', 'macaque': f'{NP}/data/monkey/snrna/snRNA_monkey.h5ad'}
SP_ORTHCOL = {'mouse': 'mouse_ensembl', 'macaque': 'cyno_symbol'}
CHUNK = 50000

def load_human_spectra():
    df = pd.read_csv(SPECTRA, sep='\t', index_col=0)
    df.index = [str(int(i)) if str(i).isdigit() else str(i).removeprefix('P') for i in df.index]
    renumber = pd.read_csv(RETAIN_MAP, sep='\t')
    renumber = renumber[renumber['status'].eq('kept')].copy()
    renumber['old_P'] = renumber['old_P'].astype(int)
    renumber['new_P'] = renumber['new_P'].astype(int)
    renumber = renumber.sort_values('new_P')
    if len(renumber) != 54 or renumber['new_P'].tolist() != list(range(1, 55)):
        raise ValueError()
    df = df.loc[renumber['old_P'].astype(str)].copy()
    df.index = [f'P{i}' for i in renumber['new_P']]
    return df

def load_ortho():
    o = pd.read_csv(ORTHO, sep='\t')
    return o

def species_common(spectra, ortho, species):
    orth_col = SP_ORTHCOL[species]
    om = ortho[['human_symbol', orth_col]].dropna().drop_duplicates()
    om = om.drop_duplicates('human_symbol').drop_duplicates(orth_col)
    h2s = dict(zip(om['human_symbol'], om[orth_col]))
    with h5py.File(H5[species], 'r') as f:
        sp_genes = [g.decode() if isinstance(g, bytes) else str(g) for g in f['var']['gene'][:]]
    sp_pos = {g: i for (i, g) in enumerate(sp_genes)}
    (common_human, common_sp_idx) = ([], [])
    for hg in spectra.columns:
        sg = h2s.get(hg)
        if sg is not None and sg in sp_pos:
            common_human.append(hg)
            common_sp_idx.append(sp_pos[sg])
    return (common_human, np.array(common_sp_idx, dtype=np.int64))

def iter_chunks(f, common_sp_idx, n_obs):
    indptr_all = f['X']['indptr'][:]
    data_d = f['X']['data']
    ind_d = f['X']['indices']
    ncol = int(f['X'].attrs['shape'][1])
    for s in range(0, n_obs, CHUNK):
        e = min(s + CHUNK, n_obs)
        (d0, d1) = (indptr_all[s], indptr_all[e])
        dat = data_d[d0:d1]
        ind = ind_d[d0:d1]
        iptr = indptr_all[s:e + 1] - d0
        Xchunk = sparse.csr_matrix((dat, ind, iptr), shape=(e - s, ncol))
        Xc = Xchunk[:, common_sp_idx].astype(np.float32)
        yield (s, e, Xc)

def lognorm_cp10k(Xc):
    tot = np.asarray(Xc.sum(1)).ravel()
    tot[tot == 0] = 1.0
    Xn = Xc.multiply(10000.0 / tot[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn

def run_species(species, spectra, ortho):
    t0 = time.time()
    (common_human, common_sp_idx) = species_common(spectra, ortho, species)
    n_common = len(common_human)
    Hload = spectra[common_human].values.astype(np.float64)
    n_prog = Hload.shape[0]
    with h5py.File(H5[species], 'r') as f:
        n_obs = int(f['X'].attrs['shape'][0])
        gsum = np.zeros(n_common, dtype=np.float64)
        gsq = np.zeros(n_common, dtype=np.float64)
        for (s, e, Xc) in iter_chunks(f, common_sp_idx, n_obs):
            Xn = lognorm_cp10k(Xc)
            Xd = np.asarray(Xn.todense(), dtype=np.float64)
            gsum += Xd.sum(0)
            gsq += (Xd * Xd).sum(0)
        mu = gsum / n_obs
        var = gsq / n_obs - mu * mu
        var[var < 1e-12] = 1e-12
        sd = np.sqrt(var)
        scores = np.empty((n_obs, n_prog), dtype=np.float64)
        for (s, e, Xc) in iter_chunks(f, common_sp_idx, n_obs):
            Xn = lognorm_cp10k(Xc)
            Xd = np.asarray(Xn.todense(), dtype=np.float64)
            Z = (Xd - mu) / sd
            scores[s:e] = Z @ Hload.T
        smu = scores.mean(0)
        ssd = scores.std(0)
        ssd[ssd == 0] = 1.0
        scores_z = (scores - smu) / ssd
        del scores
        acc = np.zeros((n_common, n_prog), dtype=np.float64)
        for (s, e, Xc) in iter_chunks(f, common_sp_idx, n_obs):
            Xn = lognorm_cp10k(Xc)
            Xd = np.asarray(Xn.todense(), dtype=np.float64)
            Z = (Xd - mu) / sd
            acc += Z.T @ scores_z[s:e]
        refit = (acc / n_obs).T

    def cos(a, b):
        (na, nb) = (np.linalg.norm(a), np.linalg.norm(b))
        if na == 0 or nb == 0:
            return np.nan
        return float(a @ b / (na * nb))
    rec = np.array([cos(Hload[p], refit[p]) for p in range(n_prog)])
    secs = round(time.time() - t0, 1)
    return {'species': species, 'common_human': common_human, 'common_sp_idx': common_sp_idx, 'Hload': Hload, 'refit': refit, 'recovery': rec, 'n_common': n_common, 'n_obs': n_obs, 'seconds': secs, 'programs': list(spectra.index)}

def pair_decay(res_a, res_b, spectra):
    set_a = {g: i for (i, g) in enumerate(res_a['common_human'])}
    set_b = {g: i for (i, g) in enumerate(res_b['common_human'])}
    shared = [g for g in res_a['common_human'] if g in set_b]
    ia = [set_a[g] for g in shared]
    ib = [set_b[g] for g in shared]
    n_prog = res_a['refit'].shape[0]
    Ra = res_a['refit'][:, ia]
    Rb = res_b['refit'][:, ib]

    def cos(a, b):
        (na, nb) = (np.linalg.norm(a), np.linalg.norm(b))
        if na == 0 or nb == 0:
            return np.nan
        return float(a @ b / (na * nb))
    c = np.array([cos(Ra[p], Rb[p]) for p in range(n_prog)])
    return (c, len(shared))

def main():
    os.makedirs(OUT, exist_ok=True)
    spectra = load_human_spectra()
    ortho = load_ortho()
    programs = list(spectra.index)
    results = {}
    for sp in ['mouse', 'macaque']:
        results[sp] = run_species(sp, spectra, ortho)
    for sp in ['mouse', 'macaque']:
        r = results[sp]
        np.savez_compressed(f'{OUT}/loadings_{sp}.npz', Hload=r['Hload'].astype(np.float32), refit=r['refit'].astype(np.float32), common_human=np.array(r['common_human'], dtype=object), programs=np.array(r['programs'], dtype=object))
    rec_df = pd.DataFrame({'program': programs})
    for sp in ['mouse', 'macaque']:
        rec_df[f'recovery_cosine_{sp}'] = results[sp]['recovery']
    rec_df.to_csv(f'{OUT}/recovery_cosine_per_program_full.csv', index=False)
    hm_mac = results['macaque']['recovery']
    hm_mou = results['mouse']['recovery']
    (mm, n_shared) = pair_decay(results['mouse'], results['macaque'], spectra)
    decay_per = pd.DataFrame({'program': programs, 'human_macaque': hm_mac, 'human_mouse': hm_mou, 'mouse_macaque': mm})
    decay_per.to_csv(f'{OUT}/decay_per_program_full.csv', index=False)

    def stats(name, arr, ncommon):
        return {'pair': name, 'mean_cosine': float(np.nanmean(arr)), 'median_cosine': float(np.nanmedian(arr)), 'n_pass_0.3': int(np.nansum(arr >= 0.3)), 'n_programs': int(len(arr)), 'n_common_genes': int(ncommon)}
    summary = pd.DataFrame([stats('human-macaque', hm_mac, results['macaque']['n_common']), stats('human-mouse', hm_mou, results['mouse']['n_common']), stats('mouse-macaque', mm, n_shared)])
    summary.to_csv(f'{OUT}/decay_summary_three_pairs_full.csv', index=False)
    meta = {'chunk': CHUNK, 'n_programs': len(programs), 'mouse_n_cells': results['mouse']['n_obs'], 'macaque_n_cells': results['macaque']['n_obs'], 'mouse_n_common': results['mouse']['n_common'], 'macaque_n_common': results['macaque']['n_common'], 'mouse_macaque_n_shared': int(n_shared), 'recovery_macaque_median': float(np.nanmedian(hm_mac)), 'recovery_mouse_median': float(np.nanmedian(hm_mou)), 'mouse_macaque_median': float(np.nanmedian(mm)), 'mouse_seconds': results['mouse']['seconds'], 'macaque_seconds': results['macaque']['seconds']}
    with open(f'{OUT}/meta_full.json', 'w') as fo:
        json.dump(meta, fo, indent=2)
if __name__ == '__main__':
    main()
