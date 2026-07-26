#!/usr/bin/env python3
import argparse
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
import torch.nn.functional as F

def retained_programs(map_path):
    frame = pd.read_csv(map_path, sep='\t')
    frame = frame[frame['new_P'].astype(str).str.fullmatch('\\d+')].copy()
    frame['old_P'] = frame['old_P'].astype(int)
    frame['new_P'] = frame['new_P'].astype(int)
    frame = frame.sort_values('new_P')
    if frame['new_P'].tolist() != list(range(1, 55)):
        raise ValueError()
    return frame

def physical_ring_kernels(step_um, upper_edges_um, dtype, device):
    max_step = int(np.ceil(max(upper_edges_um) / step_um))
    axis = torch.arange(-max_step, max_step + 1, dtype=dtype, device=device)
    (yy, xx) = torch.meshgrid(axis, axis, indexing='ij')
    distance = torch.sqrt(xx * xx + yy * yy) * step_um
    kernels = []
    lower = 0.0
    for upper in upper_edges_um:
        kernel = ((distance > lower) & (distance <= upper)).to(dtype)
        kernel[max_step, max_step] = 0
        kernels.append(kernel)
        lower = float(upper)
    return (torch.stack(kernels), max_step)

def process_section(parquet_file, row_group, chip, region, raw_programs, new_programs, output_path, kernels, max_step, device, dtype, raw_step, step_um, upper_edges_um):
    columns = ['bin', 'x', 'y'] + [f'program_{p}' for p in raw_programs]
    table = parquet_file.read_row_group(row_group, columns=columns)
    bins = table['bin'].to_pylist()
    chips_here = {str(value).split('_')[0] for value in bins}
    if chips_here != {chip}:
        raise ValueError()
    x = np.asarray(table['x'], dtype=np.int64)
    y = np.asarray(table['y'], dtype=np.int64)
    marks = np.column_stack([np.asarray(table[f'program_{p}'], dtype=np.float64) for p in raw_programs])
    finite = np.isfinite(marks).all(axis=1)
    x = x[finite]
    y = y[finite]
    marks = marks[finite]
    x0 = int(x.min())
    y0 = int(y.min())
    if np.any((x - x0) % raw_step) or np.any((y - y0) % raw_step):
        raise ValueError()
    ix = ((x - x0) // raw_step).astype(np.int64)
    iy = ((y - y0) // raw_step).astype(np.int64)
    height = int(iy.max()) + 1
    width = int(ix.max()) + 1
    if len(np.unique(iy * width + ix)) != len(ix):
        raise ValueError()
    n_program = len(raw_programs)
    grid = torch.zeros((1, n_program, height, width), dtype=dtype, device=device)
    mask = torch.zeros((1, 1, height, width), dtype=dtype, device=device)
    tx = torch.as_tensor(ix, dtype=torch.long, device=device)
    ty = torch.as_tensor(iy, dtype=torch.long, device=device)
    values = torch.as_tensor(marks, dtype=dtype, device=device)
    grid[0, :, ty, tx] = values.T
    mask[0, 0, ty, tx] = 1
    ring_count = len(upper_edges_um)
    grouped_weights = kernels[:, None, :, :].repeat(n_program, 1, 1, 1)
    with torch.no_grad():
        neighbor_sum = F.conv2d(grid, grouped_weights, padding=max_step, groups=n_program)
        neighbor_sum = neighbor_sum.reshape(1, n_program, ring_count, height, width)[0]
        neighbor_count = F.conv2d(mask, kernels[:, None, :, :], padding=max_step)[0]
        source = grid[0, :, ty, tx].T.contiguous()
        covariance = torch.empty((ring_count, n_program, n_program), dtype=dtype, device=device)
        pair_count = torch.empty(ring_count, dtype=dtype, device=device)
        for ring_index in range(ring_count):
            target = neighbor_sum[:, ring_index, ty, tx].T.contiguous()
            count = neighbor_count[ring_index, ty, tx].sum()
            pair_count[ring_index] = count
            matrix = source.T @ target / count
            covariance[ring_index] = (matrix + matrix.T) / 2
    np.savez_compressed(output_path, C=covariance.cpu().numpy().astype(np.float32), npairs=np.rint(pair_count.cpu().numpy()).astype(np.int64), chip=chip, region=region, row_group=np.int64(row_group), nbins=np.int64(len(ix)), raw_program=np.asarray(raw_programs, dtype=np.int16), new_program=np.asarray(new_programs, dtype=np.int16), raw_coordinate_step_dnb=np.float64(raw_step), physical_step_um=np.float64(step_um), ring_lower_um=np.asarray([0] + list(upper_edges_um[:-1]), dtype=np.float64), ring_upper_um=np.asarray(upper_edges_um, dtype=np.float64), interval_convention='lower_open_upper_closed')

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--root', default='/mnt/storage/home/luomeng/DATA/cortex_nmf_program')
    parser.add_argument('--output', required=True)
    parser.add_argument('--sections', default='')
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    root = Path(args.root)
    output = Path(args.output)
    per_section = output / 'per_section'
    per_section.mkdir(parents=True, exist_ok=True)
    parquet_path = root / 'results/crossregion_v1/spatial_bin50_program_score_SCT.parquet'
    map_path = root / 'results/crossregion_v1/program_renumber_map.tsv'
    chipmap_path = root / 'results/crossregion_v1/spatial_crosscorr/_chipmap.json'
    mapping = retained_programs(map_path)
    raw_programs = mapping['old_P'].tolist()
    new_programs = mapping['new_P'].tolist()
    chipmap = json.loads(chipmap_path.read_text())
    order = chipmap['order']
    selected = set(filter(None, args.sections.split(',')))
    if args.device == 'cuda' and (not torch.cuda.is_available()):
        raise RuntimeError()
    device = torch.device(args.device)
    dtype = torch.float64
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cuda.matmul.allow_tf32 = False
    raw_step = 50.0
    step_um = raw_step / 2.0
    upper_edges_um = np.arange(50.0, 501.0, 50.0)
    (kernels, max_step) = physical_ring_kernels(step_um, upper_edges_um, dtype, device)
    parquet_file = pq.ParquetFile(parquet_path)
    if parquet_file.num_row_groups != len(order):
        raise ValueError()
    run_log = []
    for (chip, region, row_group) in order:
        if selected and chip not in selected:
            continue
        destination = per_section / f'{chip}.npz'
        if destination.exists():
            run_log.append({'chip': chip, 'region': region, 'row_group': row_group, 'status': 'existing', 'seconds': 0})
            continue
        started = time.time()
        process_section(parquet_file, int(row_group), chip, region, raw_programs, new_programs, destination, kernels, max_step, device, dtype, raw_step, step_um, upper_edges_um)
        if device.type == 'cuda':
            torch.cuda.synchronize()
        elapsed = time.time() - started
        run_log.append({'chip': chip, 'region': region, 'row_group': row_group, 'status': 'computed', 'seconds': elapsed})
    pd.DataFrame(run_log).to_csv(output / 'per_section_run_log.tsv', sep='\t', index=False)
    metadata = {'source_parquet': str(parquet_path), 'program_map': str(map_path), 'chip_map': str(chipmap_path), 'coordinate_conversion': '50 DNB pixels = 25 micrometres', 'annuli_um': [{'lower': float(lower), 'upper': float(upper), 'bounds': '(lower,upper]'} for (lower, upper) in zip([0] + list(upper_edges_um[:-1]), upper_edges_um)], 'program_order': [{'raw_program': int(raw), 'retained_program': int(new)} for (raw, new) in zip(raw_programs, new_programs)], 'device': str(device), 'dtype': str(dtype)}
    (output / 'covariance_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
if __name__ == '__main__':
    main()
