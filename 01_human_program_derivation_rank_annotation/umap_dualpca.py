#!/usr/bin/env python
import importlib.util
import numpy as np
import pandas as pd
import scanpy as sc
import os
import site
import sys
import types
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
SEED = 0
np.random.seed(SEED)
PROJ = str(PROJECT_ROOT)
UMAP_CSV = f'{PROJ}/figures/fig1/_intermediate/umap_embedding.csv'
CNMF_SCORES = f'{PROJ}/results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_k60_cell_scores.tsv'
OBS_CSV = f'{PROJ}/inputs/snRNA_1M_obs.csv'
RENUMBER_MAP = f'{PROJ}/results/crossregion_v1/program_renumber_map.tsv'
DPCA_DIM = 50
scslat_candidates = [Path(location) / 'scSLAT' for location in [*getattr(site, 'getsitepackages', lambda : [])(), site.getusersitepackages(), *sys.path] if location]
scslat_dir = next((candidate for candidate in scslat_candidates if (candidate / 'model' / 'batch.py').is_file()), None)
if scslat_dir is None:
    raise ImportError()
scslat_package = types.ModuleType('scSLAT')
scslat_package.__path__ = [str(scslat_dir)]
sys.modules['scSLAT'] = scslat_package
scslat_model = types.ModuleType('scSLAT.model')
scslat_model.__path__ = [str(scslat_dir / 'model')]
sys.modules['scSLAT.model'] = scslat_model
scslat_utils_spec = importlib.util.spec_from_file_location('scSLAT.utils', scslat_dir / 'utils.py')
if scslat_utils_spec is None or scslat_utils_spec.loader is None:
    raise ImportError()
scslat_utils = importlib.util.module_from_spec(scslat_utils_spec)
sys.modules['scSLAT.utils'] = scslat_utils
scslat_utils_spec.loader.exec_module(scslat_utils)
scslat_batch_spec = importlib.util.spec_from_file_location('scSLAT.model.batch', scslat_dir / 'model' / 'batch.py')
if scslat_batch_spec is None or scslat_batch_spec.loader is None:
    raise ImportError()
scslat_batch = importlib.util.module_from_spec(scslat_batch_spec)
sys.modules['scSLAT.model.batch'] = scslat_batch
scslat_batch_spec.loader.exec_module(scslat_batch)
dual_pca = scslat_batch.dual_pca
umap_old = pd.read_csv(UMAP_CSV, index_col=0)
barcodes = umap_old.index
scores = pd.read_csv(CNMF_SCORES, sep='\t', index_col=0)
scores.columns = [str(c) for c in scores.columns]
renumber = pd.read_csv(RENUMBER_MAP, sep='\t')
retained_raw = renumber.loc[renumber['status'].eq('kept'), 'old_P'].astype(str).tolist()
usage = scores.loc[barcodes, retained_raw]
obs = pd.read_csv(OBS_CSV, index_col=0)
batch = obs.loc[barcodes, 'batch'].astype(str)
mask_edlein = (batch == 'edlein').values
mask_us = (batch == 'us').values
X_edlein = usage.values[mask_edlein].astype('float32')
Y_us = usage.values[mask_us].astype('float32')
ax = sc.AnnData(X_edlein.copy())
ay = sc.AnnData(Y_us.copy())
sc.pp.scale(ax)
sc.pp.scale(ay)
Xs = np.asarray(ax.X, dtype='float32')
Ys = np.asarray(ay.X, dtype='float32')
(Z_x, Z_y) = dual_pca(Xs, Ys, dim=DPCA_DIM, singular=True, backend='sklearn', use_gpu=False)
Z_x = np.asarray(Z_x)
Z_y = np.asarray(Z_y)
dpca = np.zeros((usage.shape[0], DPCA_DIM), dtype='float32')
dpca[mask_edlein] = Z_x
dpca[mask_us] = Z_y
ad = sc.AnnData(np.zeros((usage.shape[0], 1), dtype='float32'))
ad.obs_names = list(barcodes)
ad.obsm['X_dpca'] = dpca
sc.pp.neighbors(ad, use_rep='X_dpca', n_neighbors=30, metric='euclidean', random_state=SEED)
sc.tl.umap(ad, min_dist=0.3, random_state=SEED)
emb = ad.obsm['X_umap']
out = pd.DataFrame(index=barcodes)
out.index.name = umap_old.index.name
out['UMAP1'] = emb[:, 0]
out['UMAP2'] = emb[:, 1]
out['subclass'] = umap_old['subclass'].values
out['region'] = umap_old['region'].values
out['class'] = umap_old['class'].values
out['dominant_program'] = usage.idxmax(axis=1).astype(int).values
out['batch'] = batch.values
out.to_csv(UMAP_CSV)
