#!/usr/bin/env python3
import argparse
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
try:
    import cupy as cp
    import cupyx.scipy.sparse as csp
    HAS_CUPY = True
except Exception:
    cp = None
    csp = None
    HAS_CUPY = False
ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
RESULTS = ROOT / 'results/crossregion_v1'
OUTPUT = RESULTS / 'markcorr_betweenchip_v1'
META_FILE = RESULTS / 'spatial_bin50_meta.parquet'
SCORE_FILE = RESULTS / 'spatial_bin50_program_score_SCT.parquet'
SUBCLASS_FILE = RESULTS / 'spatial_bin50_rctd_weights.parquet'
MAP_FILE = RESULTS / 'program_renumber_map.tsv'
CELLTYPES = ['AST', 'CHANDELIER', 'ENDO', 'ET', 'L2-L3 IT LINC00507', 'L3-L4 IT RORB', 'L4-L5 IT RORB', 'L6 CAR3', 'L6 CT', 'L6 IT', 'L6B', 'LAMP5', 'MICRO', 'NDNF', 'NP', 'OLIGO', 'OPC', 'PAX6', 'PVALB', 'SST', 'VIP', 'VLMC']
UMI_FLOOR = 100
GRID_STEP_PX = 50.0
RING_UPPER_PX = 50.0
N_PERM = 1000
SEED = 42

def load_program_map():
    table = pd.read_csv(MAP_FILE, sep='\t')
    table = table[table['status'].eq('kept')].copy()
    table['old_P'] = table['old_P'].astype(int)
    table['new_P'] = table['new_P'].astype(int)
    table = table.sort_values('new_P')
    if table['new_P'].tolist() != list(range(1, 55)):
        raise ValueError()
    source = [f'program_{value}' for value in table['old_P']]
    current = [f'program_{value}' for value in table['new_P']]
    return (source, current)

def build_focal_ring(coords):
    count = len(coords)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=RING_UPPER_PX - 1e-09, output_type='ndarray')
    self_index = np.arange(count)
    if len(pairs):
        rows = np.concatenate([pairs[:, 0], pairs[:, 1], self_index])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0], self_index])
    else:
        rows = self_index
        cols = self_index
    data = np.ones(len(rows), dtype=np.float32)
    return sp.csr_matrix((data, (rows, cols)), shape=(count, count))

def calculate_g(adjacency, values_a, values_b, mean_a, mean_b, use_gpu):
    links = adjacency.nnz
    if links == 0:
        return np.ones((values_a.shape[1], values_b.shape[1]), dtype=np.float64)
    if use_gpu and HAS_CUPY:
        matrix = csp.csr_matrix(adjacency.astype(np.float32))
        array_a = cp.asarray(values_a.astype(np.float32))
        array_b = cp.asarray(values_b.astype(np.float32))
        average_a = cp.asarray(mean_a.astype(np.float64))
        average_b = cp.asarray(mean_b.astype(np.float64))
        cross = (array_a.T @ (matrix @ array_b)).astype(cp.float64)
        denominator = links * average_a[:, None] * average_b[None, :]
        result = cp.where(denominator > 0, cross / denominator, cp.ones_like(cross))
        return cp.asnumpy(result)
    matrix = adjacency.astype(np.float32)
    cross = (values_a.astype(np.float32).T @ (matrix @ values_b.astype(np.float32))).astype(np.float64)
    denominator = links * mean_a[:, None] * mean_b[None, :]
    return np.where(denominator > 0, cross / denominator, np.ones_like(cross))

def translate_mark(values, coords, shift_x, shift_y):
    origin = coords.min(axis=0)
    grid = np.round((coords - origin) / GRID_STEP_PX).astype(np.int64)
    x_index = grid[:, 0]
    y_index = grid[:, 1]
    x_size = int(x_index.max()) + 1
    y_size = int(y_index.max()) + 1
    source_x = ((x_index - shift_x) % x_size).astype(np.int64)
    source_y = (y_index - shift_y).astype(np.int64)
    valid_y = (source_y >= 0) & (source_y < y_size)
    keys = x_index * y_size + y_index
    source_keys = source_x * y_size + source_y
    order = np.argsort(keys)
    sorted_keys = keys[order]
    positions = np.searchsorted(sorted_keys, source_keys)
    in_range = valid_y & (positions < len(values))
    matched = in_range & (sorted_keys[np.minimum(positions, len(values) - 1)] == source_keys)
    translated = np.zeros_like(values)
    destination = np.where(matched)[0]
    source = order[positions[matched]]
    translated[destination] = values[source]
    return translated

def calculate_section(section, values_a, values_b, coords, permutations, seed, use_gpu):
    adjacency = build_focal_ring(coords)
    mean_a = values_a.mean(axis=0).astype(np.float64)
    mean_b = values_b.mean(axis=0).astype(np.float64)
    observed = calculate_g(adjacency, values_a, values_b, mean_a, mean_b, use_gpu)
    log_observed = np.log2(np.maximum(observed, 1e-30))
    null_mean = np.zeros_like(observed, dtype=np.float64)
    null_m2 = np.zeros_like(observed, dtype=np.float64)
    origin = coords.min(axis=0)
    grid = np.round((coords - origin) / GRID_STEP_PX).astype(np.int64)
    x_size = int(grid[:, 0].max()) + 1
    y_size = int(grid[:, 1].max()) + 1
    y_limit = max(2, y_size // 3)
    generator = np.random.default_rng(seed + hash(section) % 2 ** 31)
    for index in range(permutations):
        shift_x = int(generator.integers(1, x_size))
        shift_y = int(generator.integers(-y_limit, y_limit + 1))
        if shift_y == 0:
            shift_y = int(generator.choice([-1, 1]))
        shifted_a = translate_mark(values_a, coords, shift_x, shift_y)
        shifted_mean = shifted_a.mean(axis=0).astype(np.float64)
        null_value = calculate_g(adjacency, shifted_a, values_b, shifted_mean, mean_b, use_gpu)
        count = index + 1
        delta = null_value - null_mean
        null_mean += delta / count
        null_m2 += delta * (null_value - null_mean)
    null_sd = np.sqrt(np.maximum(null_m2 / max(permutations - 1, 1), 0)) + 1e-12
    z_value = (observed - null_mean) / null_sd
    if use_gpu and HAS_CUPY:
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    return (log_observed, z_value, null_mean, null_sd)

def load_sections(mode):
    (source_programs, current_programs) = load_program_map()
    meta = pd.read_parquet(META_FILE, columns=['bin', 'chip', 'x', 'y', 'majorDomain'])
    score = pd.read_parquet(SCORE_FILE, columns=['bin', 'bin_total_umi'] + source_programs)
    subclass = pd.read_parquet(SUBCLASS_FILE, columns=['bin', 'rctd_pass_mask'] + CELLTYPES)
    if not np.array_equal(meta['bin'].to_numpy(), score['bin'].to_numpy()):
        raise ValueError()
    if not np.array_equal(meta['bin'].to_numpy(), subclass['bin'].to_numpy()):
        raise ValueError()
    keep = subclass['rctd_pass_mask'].to_numpy(dtype=bool) & meta['majorDomain'].ne('ARACHNOID').to_numpy() & score['bin_total_umi'].ge(UMI_FLOOR).to_numpy()
    coordinates = meta[['x', 'y']].to_numpy(dtype=np.float64)
    subclass_values = np.clip(np.nan_to_num(subclass[CELLTYPES].to_numpy(dtype=np.float32), nan=0.0), 0, None)
    program_values = np.clip(score[source_programs].to_numpy(dtype=np.float32), 0, None)
    sections = {}
    chips = meta['chip'].to_numpy()
    for chip in pd.unique(chips):
        selected = (chips == chip) & keep
        if selected.sum() < 50:
            continue
        if mode == 'cellprog':
            values_a = subclass_values[selected].copy()
        else:
            values_a = program_values[selected].copy()
        sections[str(chip)] = (coordinates[selected].copy(), values_a, program_values[selected].copy())
    names_a = CELLTYPES if mode == 'cellprog' else current_programs
    return (sections, names_a, current_programs)

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--mode', required=True, choices=['cellprog', 'progprog'])
    parser.add_argument('--n-perm', type=int, default=N_PERM)
    parser.add_argument('--chip-start', type=int)
    parser.add_argument('--chip-end', type=int)
    parser.add_argument('--no-gpu', action='store_true')
    arguments = parser.parse_args()
    use_gpu = HAS_CUPY and (not arguments.no_gpu)
    (sections, names_a, names_b) = load_sections(arguments.mode)
    section_names = sorted(sections)
    if arguments.chip_start is not None:
        section_names = section_names[arguments.chip_start:arguments.chip_end]
    rows = []
    for section in section_names:
        (coords, values_a, values_b) = sections[section]
        (observed, z_value, null_mean, null_sd) = calculate_section(section, values_a, values_b, coords, arguments.n_perm, SEED, use_gpu)
        for (index_a, name_a) in enumerate(names_a):
            for (index_b, name_b) in enumerate(names_b):
                rows.append({'mode': arguments.mode, 'A_name': name_a, 'B_name': name_b, 'chip_id': section, 'log2_g_obs': float(observed[index_a, index_b]), 'Z_i': float(z_value[index_a, index_b]), 'n_in_tissue': int(len(coords)), 'mu_null': float(null_mean[index_a, index_b]), 'sd_null': float(null_sd[index_a, index_b])})
    OUTPUT.mkdir(parents=True, exist_ok=True)
    suffix = ''
    if arguments.chip_start is not None:
        end = arguments.chip_end if arguments.chip_end is not None else len(section_names)
        suffix = f'_shard{arguments.chip_start}_{end}'
    output = OUTPUT / f'betweenchip_{arguments.mode}_per_chip_Z{suffix}.tsv'
    pd.DataFrame(rows).to_csv(output, sep='\t', index=False)
if __name__ == '__main__':
    main()
