#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any
import yaml
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXPORTER = Path(__file__).with_name('export_cnmf_outputs.py')
FORMAL_COMPONENTS = (30, 40, 50, 60, 70, 80, 100, 120, 130, 140, 150, 175, 200)
FORMAL_N_ITER = 100

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--data-root', type=Path, default=Path(os.environ.get('CORTEX_PROGRAM_DATA_ROOT', REPOSITORY_ROOT / 'data')))
    parser.add_argument('--results-root', type=Path, default=Path(os.environ.get('CORTEX_PROGRAM_RESULTS_ROOT', REPOSITORY_ROOT / 'results')))
    return parser.parse_args()

def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding='utf-8'))
    config.setdefault('components', list(FORMAL_COMPONENTS))
    config.setdefault('n_iter', FORMAL_N_ITER)
    required = {'counts_h5ad', 'output_dir', 'selected_components', 'seed', 'num_highvar_genes', 'local_density_threshold'}
    absent = sorted(required.difference(config))
    if absent:
        raise ValueError()
    if tuple(config['components']) != FORMAL_COMPONENTS:
        raise ValueError()
    if config['n_iter'] != FORMAL_N_ITER:
        raise ValueError()
    if config['selected_components'] not in config['components']:
        raise ValueError()
    return config

def resolve(path_text: str, root: Path) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else root / path

def build_plan(config: dict[str, Any], data_root: Path, results_root: Path) -> list[list[str]]:
    counts = resolve(config['counts_h5ad'], data_root)
    output_dir = resolve(config['output_dir'], results_root)
    name = 'human_cnmf'
    components = [str(component) for component in config['components']]
    selected = str(config['selected_components'])
    threshold = str(config['local_density_threshold'])
    return [['cnmf', 'prepare', '--output-dir', str(output_dir), '--name', name, '-c', str(counts), '-k', *components, '--n-iter', str(config['n_iter']), '--seed', str(config['seed']), '--numgenes', str(config['num_highvar_genes'])], ['cnmf', 'factorize', '--output-dir', str(output_dir), '--name', name, '--worker-index', '0', '--total-workers', '1'], ['cnmf', 'combine', '--output-dir', str(output_dir), '--name', name], ['cnmf', 'consensus', '--output-dir', str(output_dir), '--name', name, '--components', selected, '--local-density-threshold', threshold], [sys.executable, str(EXPORTER), '--cnmf-run-dir', str(output_dir / name), '--output-dir', str(output_dir), '--name', name, '--components', selected, '--local-density-threshold', threshold]]

def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    commands = build_plan(config, args.data_root.resolve(), args.results_root.resolve())
    counts = resolve(config['counts_h5ad'], args.data_root)
    if not counts.is_file():
        raise FileNotFoundError()
    for command in commands:
        subprocess.run(command, check=True)
if __name__ == '__main__':
    main()
