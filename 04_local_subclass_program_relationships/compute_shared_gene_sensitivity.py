#!/usr/bin/env python3
from __future__ import annotations
import argparse
import base64
import math
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--candidate-pairs', required=True)
    parser.add_argument('--marker-top50', required=True)
    parser.add_argument('--tpm', required=True)
    parser.add_argument('--program-map', required=True)
    parser.add_argument('--meta', required=True)
    parser.add_argument('--rctd', required=True)
    parser.add_argument('--score-meta', required=True)
    parser.add_argument('--raw-score-dir', required=True)
    parser.add_argument('--h5ad-dir', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--pair-keys', default='')
    parser.add_argument('--pair-keys-b64', default='')
    parser.add_argument('--mapping-output', default='')
    parser.add_argument('--top-n', type=int, default=50)
    parser.add_argument('--expected-chips', type=int, default=44)
    parser.add_argument('--headline-threshold', type=float, default=0.32)
    parser.add_argument('--retention-threshold', type=float, default=0.5)
    return parser.parse_args()

def decode_strings(values: np.ndarray) -> list[str]:
    return [value.decode('utf-8', errors='replace') if isinstance(value, (bytes, np.bytes_)) else str(value) for value in values]

def h5_index(group: h5py.Group) -> list[str]:
    key = group.attrs.get('_index', '_index')
    if isinstance(key, bytes):
        key = key.decode()
    return decode_strings(group[str(key)][:])

def log(message: str) -> None:
    pass

def read_x_columns(path: Path, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    with h5py.File(path, 'r', rdcc_nbytes=256 * 1024 * 1024, rdcc_nslots=1000003) as handle:
        obs = h5_index(handle['obs'])
        if not genes:
            return (np.empty((len(obs), 0), dtype=np.float32), obs)
        var = h5_index(handle['var'])
        position = {gene: index for (index, gene) in enumerate(var)}
        encoded = handle['X']
        if encoded.attrs.get('encoding-type') != 'csr_matrix':
            raise AssertionError()
        shape = tuple((int(value) for value in encoded.attrs['shape']))
        if shape != (len(obs), len(var)):
            raise AssertionError()
        output = np.empty((shape[0], len(genes)), dtype=np.float32)
        output.fill(0.0)
        present = [(index, position[gene]) for (index, gene) in enumerate(genes) if gene in position]
        if not present:
            return (output, obs)
        output_columns = np.asarray([index for (index, _) in present], dtype=np.int64)
        wanted = np.asarray([column for (_, column) in present], dtype=np.int64)
        data = encoded['data']
        indices = encoded['indices']
        indptr = np.asarray(encoded['indptr'][:], dtype=np.int64)
        for start in range(0, shape[0], 4096):
            end = min(start + 4096, shape[0])
            (lo, hi) = (int(indptr[start]), int(indptr[end]))
            values = np.asarray(data[lo:hi], dtype=np.float32)
            columns = np.asarray(indices[lo:hi], dtype=np.int32)
            if values.size and (not np.isfinite(values).all() or np.any(values < 0)):
                raise AssertionError()
            block = sparse.csr_matrix((values, columns, indptr[start:end + 1] - lo), shape=(end - start, shape[1]), dtype=np.float32)
            output[start:end, output_columns] = block[:, wanted].toarray()
        return (output, obs)

def finalize_moments(count: int, total: float, total_sq: float) -> tuple[float, float]:
    mean = total / count
    variance = max(total_sq / count - mean * mean, 0.0)
    sd = math.sqrt(variance)
    if not np.isfinite(sd) or sd <= 0:
        raise AssertionError()
    return (mean, sd)

def mean32(values: np.ndarray) -> float:
    return float(np.asarray(values, dtype=np.float32)[:, None].mean(axis=0)[0])

def ring0_g(cell_mark: np.ndarray, program_mark: np.ndarray) -> float:
    cell_mark = np.ascontiguousarray(cell_mark, dtype=np.float32)
    program_mark = np.ascontiguousarray(program_mark, dtype=np.float32)
    denominator = np.float64(len(cell_mark)) * mean32(cell_mark) * mean32(program_mark)
    if denominator <= 0 or not np.isfinite(denominator):
        return 1.0
    return float(np.float64(np.dot(cell_mark, program_mark)) / denominator)

def parse_pair_keys(value: str, encoded: str) -> set[tuple[str, str]]:
    if value and encoded:
        raise AssertionError()
    if encoded:
        value = base64.b64decode(encoded.encode('ascii'), validate=True).decode('utf-8')
    if not value.strip():
        return set()
    keys = set()
    for item in value.split(','):
        parts = item.split('|', 1)
        if len(parts) != 2 or not all((part.strip() for part in parts)):
            raise AssertionError()
        keys.add((parts[0].strip(), parts[1].strip()))
    return keys

def prepare_pairs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = pd.read_csv(args.candidate_pairs, sep='\t')
    required = {'A', 'B', 'A_idx', 'B_idx'}
    missing = sorted(required - set(pairs.columns))
    if missing:
        raise AssertionError()
    pairs['A'] = pairs['A'].astype(str)
    pairs['B'] = pairs['B'].astype(str)
    pairs['A_idx'] = pairs['A_idx'].astype(int)
    pairs['B_idx'] = pairs['B_idx'].astype(int)
    if pairs[['A', 'B']].duplicated().any():
        raise AssertionError()
    requested = parse_pair_keys(args.pair_keys, args.pair_keys_b64)
    if requested:
        pairs = pairs[pairs.apply(lambda row: (row.A, row.B) in requested, axis=1)].copy()
        observed = set(zip(pairs['A'], pairs['B']))
        if observed != requested:
            raise AssertionError()
    if pairs.empty:
        raise AssertionError()
    markers = pd.read_csv(args.marker_top50, sep='\t')
    if not {'cell_type', 'top_genes'}.issubset(markers.columns):
        raise AssertionError()
    marker_sets = {str(row.cell_type): str(row.top_genes).split(';')[:args.top_n] for row in markers.itertuples(index=False)}
    tpm = pd.read_csv(args.tpm, sep='\t', index_col=0)
    tpm.index = tpm.index.astype(str)
    tpm.columns = tpm.columns.astype(str)
    program_rows = {program: program.removeprefix('program_') for program in pairs['B'].unique()}
    absent = [program for (program, row) in program_rows.items() if row not in tpm.index]
    if absent:
        raise AssertionError()
    map_frame = pd.read_csv(args.program_map, sep='\t')
    if not {'new_P', 'cnmf_component'}.issubset(map_frame.columns):
        raise AssertionError()
    map_frame['cnmf_component'] = map_frame['cnmf_component'].astype(int)
    if int((map_frame['new_P'].astype(str) != 'EXCLUDED').sum()) != 54:
        raise AssertionError()
    component_to_label = dict(zip(map_frame['cnmf_component'], map_frame['new_P'].astype(str)))
    records = []
    for row in pairs.itertuples(index=False):
        if row.A not in marker_sets:
            raise AssertionError()
        component = int(row.B.removeprefix('program_'))
        if component - 1 != int(row.B_idx):
            raise AssertionError()
        if component not in component_to_label:
            raise AssertionError()
        top_program = list(tpm.loc[program_rows[row.B]].sort_values(ascending=False, kind='mergesort').head(args.top_n).index)
        shared = sorted(set(marker_sets[row.A]) & set(top_program))
        n_shared = len(shared)
        jaccard = n_shared / len(set(marker_sets[row.A]) | set(top_program))
        if 'n_shared' in pairs.columns and int(row.n_shared) != n_shared:
            raise AssertionError()
        if 'jaccard' in pairs.columns and (not np.isclose(float(row.jaccard), jaccard, atol=1e-12, rtol=0)):
            raise AssertionError()
        records.append({'A': row.A, 'B': row.B, 'A_idx': int(row.A_idx), 'B_idx': int(row.B_idx), 'n_shared': n_shared, 'jaccard': jaccard, 'shared_genes': shared, 'program_label': component_to_label[component]})
    output = pd.DataFrame(records)
    return (output, map_frame)

def load_global_inputs(args: argparse.Namespace, pairs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, slice], dict[str, np.ndarray]]:
    cells = list(dict.fromkeys(pairs['A'].astype(str)))
    meta = pd.read_parquet(args.meta, columns=['bin', 'chip', 'majorDomain'])
    bins = meta['bin'].astype(str).to_numpy()
    chips = meta['chip'].astype(str).to_numpy()
    domains = meta['majorDomain'].astype(str).to_numpy()
    starts = np.flatnonzero(np.r_[True, chips[1:] != chips[:-1]])
    ends = np.r_[starts[1:], len(chips)]
    chip_order = [str(chips[start]) for start in starts]
    if len(chip_order) != args.expected_chips or len(set(chip_order)) != args.expected_chips:
        raise AssertionError()
    chip_slices = {chip: slice(int(start), int(end)) for (chip, start, end) in zip(chip_order, starts, ends)}
    rctd = pd.read_parquet(args.rctd, columns=['bin', 'rctd_pass_mask'] + cells)
    if not np.array_equal(rctd['bin'].astype(str).to_numpy(), bins):
        raise AssertionError()
    weights = {cell: np.clip(np.nan_to_num(rctd[cell].to_numpy(dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0), 0, None) for cell in cells}
    score_meta = pd.read_parquet(args.score_meta, columns=['bin', 'bin_total_umi'])
    if not np.array_equal(score_meta['bin'].astype(str).to_numpy(), bins):
        raise AssertionError()
    keep = rctd['rctd_pass_mask'].to_numpy(dtype=bool) & (domains != 'ARACHNOID') & (score_meta['bin_total_umi'].to_numpy(dtype=np.float64) >= 100)
    return (bins, keep, chip_order, chip_slices, weights)

def collect_stats(args: argparse.Namespace, pairs: pd.DataFrame, bins: np.ndarray, chip_order: list[str], chip_slices: dict[str, slice], tpm: pd.DataFrame) -> tuple[dict[str, tuple[np.float32, np.float32]], dict[int, tuple[float, float]]]:
    programs = list(dict.fromkeys(pairs['B'].astype(str)))
    p_index = {program: index for (index, program) in enumerate(programs)}
    gene_union = sorted({gene for genes in pairs['shared_genes'] for gene in genes})
    gene_index = {gene: index for (index, gene) in enumerate(gene_union)}
    raw_sums = np.zeros(len(programs), dtype=np.float64)
    raw_sumsq = np.zeros(len(programs), dtype=np.float64)
    logo_sums = np.zeros(len(pairs), dtype=np.float64)
    logo_sumsq = np.zeros(len(pairs), dtype=np.float64)
    n_total = 0
    h5_dir = Path(args.h5ad_dir)
    raw_dir = Path(args.raw_score_dir)
    for (number, chip) in enumerate(chip_order, start=1):
        log(f'statistics {number}/{len(chip_order)} {chip}')
        expected = bins[chip_slices[chip]]
        (x, h5_bins) = read_x_columns(h5_dir / f'{chip}_bin50.h5ad', gene_union)
        if not np.array_equal(np.asarray(h5_bins, dtype=object), expected):
            raise AssertionError()
        raw_frame = pd.read_parquet(raw_dir / f'{chip}.parquet', columns=['bin'] + programs)
        if not np.array_equal(raw_frame['bin'].astype(str).to_numpy(), expected):
            raise AssertionError()
        raw = np.ascontiguousarray(raw_frame[programs].to_numpy(dtype=np.float32), dtype=np.float32)
        raw64 = raw.astype(np.float64, copy=False)
        raw_sums += raw64.sum(axis=0)
        raw_sumsq += (raw64 * raw64).sum(axis=0)
        n_total += len(raw)
        for (pair_index, pair) in enumerate(pairs.itertuples(index=False)):
            if pair.n_shared == 0:
                candidate = raw[:, p_index[pair.B]]
            else:
                positions = [gene_index[gene] for gene in pair.shared_genes]
                loading = tpm.loc[pair.B.removeprefix('program_'), pair.shared_genes].to_numpy(dtype=np.float32)
                contribution = np.asarray((x[:, positions] * loading[None, :]).sum(axis=1, dtype=np.float32), dtype=np.float32)
                candidate = np.asarray(raw[:, p_index[pair.B]] - contribution, dtype=np.float32)
            values = candidate.astype(np.float64, copy=False)
            logo_sums[pair_index] += values.sum()
            logo_sumsq[pair_index] += np.dot(values, values)
    raw_parameters = {}
    for (program, index) in p_index.items():
        (mean, sd) = finalize_moments(n_total, raw_sums[index], raw_sumsq[index])
        raw_parameters[program] = (np.float32(mean), np.float32(sd))
    logo_parameters = {index: finalize_moments(n_total, logo_sums[index], logo_sumsq[index]) for index in range(len(pairs))}
    return (raw_parameters, logo_parameters)

def compute_results(args: argparse.Namespace, pairs: pd.DataFrame, bins: np.ndarray, keep: np.ndarray, chip_order: list[str], chip_slices: dict[str, slice], weights: dict[str, np.ndarray], tpm: pd.DataFrame, raw_parameters: dict[str, tuple[np.float32, np.float32]], logo_parameters: dict[int, tuple[float, float]]) -> pd.DataFrame:
    programs = list(dict.fromkeys(pairs['B'].astype(str)))
    p_index = {program: index for (index, program) in enumerate(programs)}
    gene_union = sorted({gene for genes in pairs['shared_genes'] for gene in genes})
    gene_index = {gene: index for (index, gene) in enumerate(gene_union)}
    original_g = [[] for _ in range(len(pairs))]
    logo_g = [[] for _ in range(len(pairs))]
    h5_dir = Path(args.h5ad_dir)
    raw_dir = Path(args.raw_score_dir)
    for (number, chip) in enumerate(chip_order, start=1):
        log(f'statistics {number}/{len(chip_order)} {chip}')
        expected = bins[chip_slices[chip]]
        local_keep = keep[chip_slices[chip]]
        (x, h5_bins) = read_x_columns(h5_dir / f'{chip}_bin50.h5ad', gene_union)
        if not np.array_equal(np.asarray(h5_bins, dtype=object), expected):
            raise AssertionError()
        raw_frame = pd.read_parquet(raw_dir / f'{chip}.parquet', columns=['bin'] + programs)
        if not np.array_equal(raw_frame['bin'].astype(str).to_numpy(), expected):
            raise AssertionError()
        raw = np.ascontiguousarray(raw_frame[programs].to_numpy(dtype=np.float32), dtype=np.float32)
        for (pair_index, pair) in enumerate(pairs.itertuples(index=False)):
            raw_score = raw[:, p_index[pair.B]]
            (raw_mean, raw_sd) = raw_parameters[pair.B]
            baseline = np.clip((raw_score - raw_mean) / raw_sd, 0, None).astype(np.float32)
            cell_mark = weights[pair.A][chip_slices[chip]][local_keep]
            original_g[pair_index].append(ring0_g(cell_mark, baseline[local_keep]))
            if pair.n_shared == 0:
                logo_mark = baseline
            else:
                positions = [gene_index[gene] for gene in pair.shared_genes]
                loading = tpm.loc[pair.B.removeprefix('program_'), pair.shared_genes].to_numpy(dtype=np.float32)
                contribution = np.asarray((x[:, positions] * loading[None, :]).sum(axis=1, dtype=np.float32), dtype=np.float32)
                raw_logo = np.asarray(raw_score - contribution, dtype=np.float32)
                (logo_mean, logo_sd) = logo_parameters[pair_index]
                logo_mark = np.clip((raw_logo - logo_mean) / logo_sd, 0, None).astype(np.float32)
            logo_g[pair_index].append(ring0_g(cell_mark, logo_mark[local_keep]))
    rows = []
    for (pair_index, pair) in enumerate(pairs.itertuples(index=False)):
        original = float(np.log2(np.median(np.asarray(original_g[pair_index], dtype=np.float64))))
        logo = float(np.log2(np.median(np.asarray(logo_g[pair_index], dtype=np.float64))))
        delta = logo - original
        retention = abs(logo) / abs(original) if original != 0 else np.nan
        survives = bool(np.sign(original) == np.sign(logo) and abs(logo) > args.headline_threshold and np.isfinite(retention) and (retention >= args.retention_threshold))
        rows.append({'A': pair.A, 'B': pair.B, 'A_idx': pair.A_idx, 'B_idx': pair.B_idx, 'n_shared': pair.n_shared, 'jaccard': pair.jaccard, 'orig_log2g': original, 'logo_log2g': logo, 'delta_log2g': delta, 'retention': retention, 'survives': survives})
    return pd.DataFrame(rows)

def main() -> int:
    args = parse_args()
    (pairs, map_frame) = prepare_pairs(args)
    tpm = pd.read_csv(args.tpm, sep='\t', index_col=0)
    tpm.index = tpm.index.astype(str)
    tpm.columns = tpm.columns.astype(str)
    (bins, keep, chip_order, chip_slices, weights) = load_global_inputs(args, pairs)
    (raw_parameters, logo_parameters) = collect_stats(args, pairs, bins, chip_order, chip_slices, tpm)
    output = compute_results(args, pairs, bins, keep, chip_order, chip_slices, weights, tpm, raw_parameters, logo_parameters)
    columns = ['A', 'B', 'A_idx', 'B_idx', 'n_shared', 'jaccard', 'orig_log2g', 'logo_log2g', 'delta_log2g', 'retention', 'survives']
    output = output.loc[:, columns]
    if args.output == '-':
        sys.stdout.write(output.to_csv(sep='\t', index=False))
    else:
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        output.to_csv(path, sep='\t', index=False)
    if args.mapping_output:
        mapping = map_frame.loc[:, ['cnmf_component', 'new_P']].copy()
        mapping.to_csv(args.mapping_output, sep='\t', index=False)
    return 0
if __name__ == '__main__':
    raise SystemExit()
