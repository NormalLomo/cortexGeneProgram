#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


LAYERS = ['L1', 'L2/3', 'L4', 'L5', 'L6']
SPECIES = ['human', 'macaque', 'mouse']
PROGRAMS = [f'program_{value}' for value in range(1, 61)]


def zscore(values: np.ndarray) -> np.ndarray:
    mean = values.mean(axis=0)
    sd = values.std(axis=0, ddof=0)
    sd[sd == 0] = 1.0
    return (values - mean) / sd


def map_layer(species: str, value: object) -> str | None:
    text = str(value).lower()
    if species == 'human':
        return {'l1': 'L1', 'l2': 'L2/3', 'l3': 'L2/3', 'l3a': 'L2/3', 'l3b': 'L2/3', 'l4': 'L4', 'l5': 'L5', 'l5a': 'L5', 'l5b': 'L5', 'l6': 'L6', 'l6b': 'L6'}.get(text)
    if species == 'macaque':
        return {'l1': 'L1', 'l2': 'L2/3', 'l3': 'L2/3', 'l4': 'L4', 'l5': 'L5', 'l6': 'L6'}.get(text)
    return {'1': 'L1', '2/3': 'L2/3', '4': 'L4', '5': 'L5', '6': 'L6', '6a': 'L6'}.get(text)


def preference(profile: pd.Series) -> float:
    return float(profile[['L1', 'L2/3']].mean() - profile[['L5', 'L6']].mean())


def classify(values: list[float]) -> tuple[str, str]:
    signs = np.sign(values).astype(int)
    same = bool(np.all(signs == signs[0]) and signs[0] != 0)
    if same and np.all(np.abs(values) >= 0.2):
        return ('A', 'shared')
    if same:
        return ('B', 'weak_common')
    return ('C', 'interspecific_variation')


def section_rows(species: str, section: str, frame: pd.DataFrame, layer_column: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    mapped = frame[layer_column].map(lambda value: map_layer(species, value))
    scores = zscore(frame[PROGRAMS].to_numpy(np.float64))
    rows = []
    coverage = []
    for layer in LAYERS:
        mask = mapped.eq(layer).to_numpy()
        count = int(mask.sum())
        coverage.append({'species': species, 'section': section, 'layer': layer, 'n_bins': count, 'included': count >= 20})
        if count < 20:
            continue
        means = scores[mask].mean(axis=0)
        rows.extend({'species': species, 'section': section, 'layer': layer, 'source_program': program, 'section_mean': float(value)} for program, value in zip(PROGRAMS, means))
    return rows, coverage


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--project-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root
    os.chdir(root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = pd.read_csv(root / 'results/crossregion_v1/program_names.tsv', sep='\t')
    names['source_program'] = 'program_' + names['cnmf_component'].astype(int).astype(str)
    names = names.loc[names['new_P'].ne('EXCLUDED'), ['source_program', 'new_P', 'cnmf_component', 'name_short']].copy()
    names['program_number'] = names['new_P'].str.removeprefix('P').astype(int)
    names = names.sort_values('program_number').drop(columns='program_number')
    observations = []
    coverage = []
    human = sorted(glob.glob(str(root / 'results/crossregion_v1/_spatial_score_perchip/*.parquet')))
    for path in human:
        frame = pd.read_parquet(path, columns=PROGRAMS + ['domain', 'rctd_pass_mask'])
        frame = frame.loc[frame['rctd_pass_mask'].eq(True)]
        rows, counts = section_rows('human', Path(path).stem, frame, 'domain')
        observations.extend(rows)
        coverage.extend(counts)
    macaque = sorted(glob.glob(str(root / 'results/xspecies_humanmap_v1/spatial_xspecies/monkey/results_sct/sections_bin50/*.parquet')))
    for path in macaque:
        frame = pd.read_parquet(path, columns=PROGRAMS + ['layer'])
        rows, counts = section_rows('macaque', Path(path).stem, frame, 'layer')
        observations.extend(rows)
        coverage.extend(counts)
    mouse = pd.read_parquet(root / 'results/xspecies_humanmap_v1/spatial_xspecies/mouse/sections_sct/mouse_bin50_program_score_SCT.parquet', columns=PROGRAMS + ['layer', 'section'])
    for section, frame in mouse.groupby('section', sort=True):
        rows, counts = section_rows('mouse', str(section), frame, 'layer')
        observations.extend(rows)
        coverage.extend(counts)
    observations = pd.DataFrame(observations)
    primary = observations.groupby(['species', 'source_program', 'layer'], observed=True)['section_mean'].median().unstack('layer').reindex(columns=LAYERS)
    display = primary.groupby(level='species', group_keys=False).apply(lambda frame: frame.subtract(frame.median(axis=0), axis=1), include_groups=False)
    profile_z = primary.groupby(level='species', group_keys=False).apply(lambda frame: frame.sub(frame.mean(axis=1), axis=0).div(frame.std(axis=1, ddof=0).replace(0, 1), axis=0), include_groups=False)
    records = []
    profile_records = []
    for entry in names.itertuples(index=False):
        values = []
        for species in SPECIES:
            raw_profile = profile_z.loc[(species, entry.source_program)]
            values.append(preference(raw_profile))
            for layer in LAYERS:
                profile_records.append({'new_P': entry.new_P, 'source_program': entry.source_program, 'name_short': entry.name_short, 'species': species, 'layer': layer, 'display_score': float(display.loc[(species, entry.source_program), layer]), 'within_program_z': float(raw_profile[layer])})
        code, label = classify(values)
        records.append({'new_P': entry.new_P, 'source_program': entry.source_program, 'cnmf_component': int(entry.cnmf_component), 'name_short': entry.name_short, 'human_preference': values[0], 'macaque_preference': values[1], 'mouse_preference': values[2], 'class_code': code, 'class_label': label, 'human_macaque_spearman': float(profile_z.loc[('human', entry.source_program)].corr(profile_z.loc[('macaque', entry.source_program)], method='spearman')), 'human_mouse_spearman': float(profile_z.loc[('human', entry.source_program)].corr(profile_z.loc[('mouse', entry.source_program)], method='spearman'))})
    classification = pd.DataFrame(records).sort_values('new_P', key=lambda value: value.str.removeprefix('P').astype(int))
    classification.to_csv(args.output_dir / 'fig7a_laminar_classification.tsv', sep='\t', index=False)
    pd.DataFrame(profile_records).to_csv(args.output_dir / 'fig7b_laminar_profiles.tsv', sep='\t', index=False)
    pd.DataFrame(coverage).to_csv(args.output_dir / 'fig7_laminar_coverage.tsv', sep='\t', index=False)
    classification['class_code'].value_counts().reindex(['A', 'B', 'C'], fill_value=0).rename_axis('class_code').reset_index(name='program_count').to_csv(args.output_dir / 'fig7_laminar_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
