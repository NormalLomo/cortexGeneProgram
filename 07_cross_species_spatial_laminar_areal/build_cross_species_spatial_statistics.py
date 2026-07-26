#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import f_oneway
from cross_species_spatial_statistics import COMMON_LAYERS, benjamini_hochberg, build_pair_index, classify_laminar, normalize_layer, partial_spearman_matrix, retained_program_map, wilson_interval
import os
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
PROJECT = Path('/mnt/storage/home/luomeng/DATA/cortex_nmf_program')
SECTION_DIRS = {'human': PROJECT / 'results/crossregion_v1/_SCT_score_perchip', 'macaque': PROJECT / 'results/xspecies_humanmap_v1/spatial_xspecies/monkey/results_sct/sections_bin50', 'mouse': PROJECT / 'results/xspecies_humanmap_v1/spatial_xspecies/mouse/sections_sct'}
LAYER_COLUMNS = {'human': 'majorDomain', 'macaque': 'layer', 'mouse': 'layer'}
CATEGORY_MEMBERS = {'microglia': ['P36', 'P45', 'P49', 'P54'], 'synaptic': ['P16', 'P26', 'P43'], 'neuropil': ['P8', 'P21'], 'myelin': ['P12', 'P23', 'P32', 'P33']}
CATEGORY_CONTRASTS = {'microglia_synaptic': ('microglia', 'synaptic'), 'neuropil_myelin': ('neuropil', 'myelin')}
MIN_SECTION_BINS = 50
MIN_AREA_BINS = 20

def section_files(species: str) -> list[Path]:
    files = sorted(SECTION_DIRS[species].glob('*.parquet'))
    if species == 'mouse':
        files = [p for p in files if 'program_score' not in p.name]
    expected = 49 if species == 'mouse' else 44
    if len(files) != expected:
        raise RuntimeError()
    return files

def section_scale(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    means = np.nanmean(values, axis=0)
    sds = np.nanstd(values, axis=0, ddof=0)
    sds[sds == 0] = np.nan
    return (values - means) / sds

def category_membership(program_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_program = program_map.set_index('new_P')
    membership_rows = []
    pair_rows = []
    for (category, programs) in CATEGORY_MEMBERS.items():
        for program in programs:
            row = by_program.loc[program]
            membership_rows.append({'category': category, 'new_P': program, 'cnmf_component': int(row['cnmf_component']), 'score_column': row['score_column'], 'name_short': row['name_short']})
    for (contrast, (left, right)) in CATEGORY_CONTRASTS.items():
        for a in CATEGORY_MEMBERS[left]:
            for b in CATEGORY_MEMBERS[right]:
                (program_a, program_b) = (a, b)
                (a_num, b_num) = (int(program_a[1:]), int(program_b[1:]))
                if a_num > b_num:
                    (program_a, program_b) = (program_b, program_a)
                pair_rows.append({'contrast': contrast, 'category_a': left, 'category_b': right, 'program_a': program_a, 'program_b': program_b, 'pair_id': f'{program_a}__{program_b}'})
    return (pd.DataFrame(membership_rows), pd.DataFrame(pair_rows))

def load_section(path: Path, species: str, score_columns: list[str]) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    metadata = [LAYER_COLUMNS[species], 'bin_total_umi']
    if species == 'mouse':
        metadata.append('area')
    table = pq.read_table(path, columns=metadata + score_columns)
    frame = table.to_pandas()
    layers = frame[LAYER_COLUMNS[species]].map(lambda value: normalize_layer(species, value))
    common = layers.notna().to_numpy()
    scores = frame.loc[common, score_columns].to_numpy(dtype=np.float64, copy=False)
    umi = frame.loc[common, 'bin_total_umi'].to_numpy(dtype=np.float64, copy=False)
    finite = np.isfinite(umi) & (umi >= 0) & np.isfinite(scores).all(axis=1)
    filtered = frame.loc[common].loc[finite].copy()
    filtered['common_layer'] = layers.loc[common].loc[finite].to_numpy()
    return (filtered, scores[finite], umi[finite])

def pair_rows_for_section(species: str, section: str, scores: np.ndarray, umi: np.ndarray, pair_index: pd.DataFrame, program_to_index: dict[str, int]) -> pd.DataFrame:
    if scores.shape[0] < MIN_SECTION_BINS:
        correlations = np.full(len(pair_index), np.nan)
        status = 'insufficient_bins'
    else:
        matrix = partial_spearman_matrix(scores, np.log1p(umi))
        ia = pair_index['program_a'].map(program_to_index).to_numpy()
        ib = pair_index['program_b'].map(program_to_index).to_numpy()
        correlations = matrix[ia, ib]
        status = 'available'
    out = pair_index.copy()
    out.insert(0, 'species', species)
    out.insert(1, 'section', section)
    out['n_bins_complete'] = scores.shape[0]
    out['partial_spearman'] = correlations
    out['availability'] = status
    return out

def aggregate_layer_stats(layer_section: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = layer_section.groupby(['species', 'new_P', 'common_layer'], observed=True)['section_layer_mean'].agg(layer_median='median', layer_mean='mean', n_sections='count').reset_index()
    grouped['layer_z_within_program_species'] = np.nan
    for ((_, _), index) in grouped.groupby(['species', 'new_P']).groups.items():
        rows = grouped.loc[index].set_index('common_layer').reindex(COMMON_LAYERS)
        values = rows['layer_median'].to_numpy(dtype=float)
        if np.isfinite(values).all() and np.nanstd(values) > 0:
            z = (values - np.mean(values)) / np.std(values, ddof=0)
            z_by_layer = dict(zip(COMMON_LAYERS, z))
            grouped.loc[index, 'layer_z_within_program_species'] = grouped.loc[index, 'common_layer'].map(z_by_layer)
    wide = grouped.pivot(index='new_P', columns=['species', 'common_layer'], values='layer_z_within_program_species')
    rows = []
    for program in sorted(wide.index, key=lambda value: int(value[1:])):
        prefs: dict[str, float] = {}
        rec: dict[str, object] = {'new_P': program, 'tau': 0.2}
        for species in ('human', 'macaque', 'mouse'):
            values = {layer: wide.loc[program].get((species, layer), np.nan) for layer in COMMON_LAYERS}
            pref = np.mean([values['L1'], values['L2']]) - np.mean([values['L5'], values['L6']])
            prefs[species] = float(pref)
            rec[f'{species}_preference'] = float(pref)
            rec[f'{species}_sign'] = int(np.sign(pref)) if np.isfinite(pref) else 0
        (rec['class_ABC'], rec['divergence_type']) = classify_laminar(prefs['human'], prefs['macaque'], prefs['mouse'], tau=0.2)
        rows.append(rec)
    return (grouped, pd.DataFrame(rows))

def aggregate_area_stats(area_section: pd.DataFrame, table_s4: pd.DataFrame) -> pd.DataFrame:
    area_summary = area_section.groupby(['new_P', 'area'], observed=True)['section_area_median'].agg(area_median='median', n_section_area='count').reset_index()
    rows = []
    for (program, group) in area_section.groupby('new_P', observed=True):
        area_values = area_summary.loc[area_summary['new_P'] == program, 'area_median'].to_numpy(float)
        anova_groups = [values['section_area_median'].to_numpy(float) for (_, values) in group.groupby('area', observed=True) if len(values) >= 2]
        if len(anova_groups) >= 2 and sum(map(len, anova_groups)) > len(anova_groups):
            test = f_oneway(*anova_groups)
            (f_stat, p_value) = (float(test.statistic), float(test.pvalue))
            df_between = len(anova_groups) - 1
            df_within = sum(map(len, anova_groups)) - len(anova_groups)
        else:
            f_stat = p_value = np.nan
            df_between = df_within = 0
        rows.append({'new_P': program, 'area_median_sd': float(np.std(area_values, ddof=1)) if len(area_values) > 1 else np.nan, 'n_areas': int(len(area_values)), 'n_section_area_observations': int(len(group)), 'anova_n_areas_with_replicates': int(len(anova_groups)), 'anova_df_between': int(df_between), 'anova_df_within': int(df_within), 'anova_F': f_stat, 'anova_p': p_value, 'observation_unit': 'section_area_median', 'min_bins_per_section_area': MIN_AREA_BINS})
    out = pd.DataFrame(rows)
    out['anova_q_bh_retained54'] = benjamini_hochberg(out['anova_p'].to_numpy())
    out = out.sort_values('area_median_sd', ascending=False).reset_index(drop=True)
    out['area_variability_rank'] = np.arange(1, len(out) + 1)
    out['top14_display'] = out['area_variability_rank'] <= 14
    flags = table_s4.dropna(subset=['new_P']).set_index('new_P')
    if not flags.index.is_unique:
        raise ValueError()
    robust = flags['validity_flag'].astype(str).str.startswith('cohort-robust')
    variable = flags['region_variability_class'].astype(str).str.lower().eq('true')
    out['human_region_variable'] = out['new_P'].map(variable).fillna(False)
    out['human_primary_flag'] = out['new_P'].map(variable & robust).fillna(False)
    return out

def aggregate_pair_universe(pair_section: pd.DataFrame, pair_index: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    aggregate = pair_section.groupby(['species', 'pair_id'], observed=True)['partial_spearman'].agg(section_median='median', section_mean='mean', n_sections='count').reset_index()
    values = pair_index.copy()
    for species in ('human', 'macaque', 'mouse'):
        sp = aggregate.loc[aggregate['species'] == species].set_index('pair_id')
        values[f'{species}_section_median'] = values['pair_id'].map(sp['section_median'])
        values[f'{species}_section_mean'] = values['pair_id'].map(sp['section_mean'])
        values[f'{species}_n_sections'] = values['pair_id'].map(sp['n_sections']).fillna(0).astype(int)
    for species in ('macaque', 'mouse'):
        h = values['human_section_median'].to_numpy(float)
        x = values[f'{species}_section_median'].to_numpy(float)
        available = np.isfinite(h) & np.isfinite(x) & (h != 0) & (x != 0)
        values[f'{species}_available'] = available
        values[f'{species}_same_sign'] = np.where(available, np.sign(h) == np.sign(x), pd.NA)
    summary_rows = []
    for species in ('macaque', 'mouse'):
        for (universe, mask) in (('inclusive_with_diagonal', np.ones(len(values), dtype=bool)), ('off_diagonal', ~values['is_diagonal'].to_numpy())):
            available = mask & values[f'{species}_available'].to_numpy()
            same = values.loc[available, f'{species}_same_sign'].astype(bool)
            (numerator, denominator) = (int(same.sum()), int(len(same)))
            (lo, hi) = wilson_interval(numerator, denominator)
            summary_rows.append({'comparison': f'human_vs_{species}', 'pair_universe': universe, 'pair_count': denominator, 'same_sign_numerator': numerator, 'same_sign_fraction': numerator / denominator if denominator else np.nan, 'wilson_95ci_low': lo, 'wilson_95ci_high': hi})
    return (values, pd.DataFrame(summary_rows))

def aggregate_category(pair_section: pd.DataFrame, category_pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = pair_section.merge(category_pairs, on=['pair_id', 'program_a', 'program_b'], how='inner')
    per_pair = detail.groupby(['species', 'contrast', 'pair_id'], observed=True)['partial_spearman'].agg(pair_section_median='median', pair_section_mean='mean', n_available_sections='count').reset_index()
    stats = per_pair.groupby(['species', 'contrast'], observed=True)['pair_section_median'].agg(category_pair_median='median', category_pair_mean='mean', n_pairs='count').reset_index()
    human = stats.loc[stats['species'] == 'human'].set_index('contrast')['category_pair_median']
    stats['human_reference_sign'] = stats['contrast'].map(np.sign(human)).astype(int)
    stats['sign_matches_human'] = np.sign(stats['category_pair_median']) == stats['human_reference_sign']
    return (detail, stats)

def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--table-s4', required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    names = pd.read_csv(PROJECT / 'results/crossregion_v1/program_names.tsv', sep='\t')
    program_map = retained_program_map(names)
    if len(program_map) != 54:
        raise RuntimeError()
    table_s4 = pd.read_csv(args.table_s4, sep='\t')
    (membership, category_pairs) = category_membership(program_map)
    membership.to_csv(args.outdir / 'fig7_category_membership.tsv', sep='\t', index=False)
    programs = program_map['new_P'].tolist()
    score_columns = program_map['score_column'].tolist()
    program_to_index = {program: index for (index, program) in enumerate(programs)}
    pair_index = build_pair_index(programs, include_diagonal=True)
    pair_section_parts = []
    layer_section_parts = []
    area_section_parts = []
    section_inventory = []
    for species in ('human', 'macaque', 'mouse'):
        for path in section_files(species):
            section = path.stem
            (metadata, scores, umi) = load_section(path, species, score_columns)
            section_inventory.append({'species': species, 'section': section, 'source_path': str(path), 'n_complete_common_layer_bins': len(scores), 'n_common_layers': metadata['common_layer'].nunique(), 'median_bin_total_umi': float(np.median(umi))})
            pair_section_parts.append(pair_rows_for_section(species, section, scores, umi, pair_index, program_to_index))
            scaled = section_scale(scores)
            for layer in COMMON_LAYERS:
                layer_mask = metadata['common_layer'].to_numpy() == layer
                if not layer_mask.any():
                    continue
                means = np.nanmean(scaled[layer_mask], axis=0)
                layer_section_parts.append(pd.DataFrame({'species': species, 'section': section, 'common_layer': layer, 'new_P': programs, 'section_layer_mean': means, 'n_bins': int(layer_mask.sum())}))
            if species == 'mouse':
                metadata = metadata.reset_index(drop=True)
                areas = metadata['area'].astype(str).to_numpy()
                for area in sorted(set(areas)):
                    area_mask = areas == area
                    if int(area_mask.sum()) < MIN_AREA_BINS:
                        continue
                    medians = np.nanmedian(scaled[area_mask], axis=0)
                    area_section_parts.append(pd.DataFrame({'section': section, 'area': area, 'new_P': programs, 'section_area_median': medians, 'n_bins': int(area_mask.sum())}))
    pair_section = pd.concat(pair_section_parts, ignore_index=True)
    layer_section = pd.concat(layer_section_parts, ignore_index=True)
    area_section = pd.concat(area_section_parts, ignore_index=True)
    pd.DataFrame(section_inventory).to_csv(args.outdir / 'fig7_section_inventory.tsv', sep='\t', index=False)
    pair_section.to_csv(args.outdir / 'fig7_pair_section.tsv', sep='\t', index=False)
    layer_section.to_csv(args.outdir / 'fig7_layer_section.tsv', sep='\t', index=False)
    area_section.to_csv(args.outdir / 'fig7_mouse_area_section.tsv', sep='\t', index=False)
    (layer_stats, classification) = aggregate_layer_stats(layer_section)
    layer_stats.to_csv(args.outdir / 'fig7_layer_stats.tsv', sep='\t', index=False)
    classification.to_csv(args.outdir / 'fig7_laminar_classification.tsv', sep='\t', index=False)
    (category_detail, category_stats) = aggregate_category(pair_section, category_pairs)
    category_detail.to_csv(args.outdir / 'fig7_category_pair_section.tsv', sep='\t', index=False)
    category_stats.to_csv(args.outdir / 'fig7_category_stats.tsv', sep='\t', index=False)
    area_stats = aggregate_area_stats(area_section, table_s4)
    area_stats.to_csv(args.outdir / 'fig7_area_stats.tsv', sep='\t', index=False)
    (pair_universe, pair_summary) = aggregate_pair_universe(pair_section, pair_index)
    pair_universe.to_csv(args.outdir / 'fig7_pair_universe.tsv', sep='\t', index=False)
    pair_summary.to_csv(args.outdir / 'fig7_pair_summary.tsv', sep='\t', index=False)
    summary = {'retained_programs': len(program_map), 'cnmf_replicates': 100, 'sections': pd.DataFrame(section_inventory).groupby('species')['section'].nunique().to_dict(), 'pair_universe_inclusive': len(pair_universe), 'pair_universe_off_diagonal': int((~pair_universe['is_diagonal']).sum()), 'human_primary_region_variable_n': int(area_stats['human_primary_flag'].sum())}
    (args.outdir / 'fig7_numeric_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
if __name__ == '__main__':
    main()
