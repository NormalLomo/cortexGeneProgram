#!/usr/bin/env python3
from __future__ import annotations

def _annotation_group(value):
    return 'class_a' if str(value).endswith('sig') else 'class_b'

def _annotation_table(table):
    if hasattr(table, 'columns') and 'confidence' in table.columns:
        table = table.copy()
        table['confidence'] = table['confidence'].map(_annotation_group)
    return table
import argparse
from pathlib import Path
import matplotlib
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
LAYERS = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']
COMMON_LAYERS = ['L1', 'L2', 'L4', 'L5', 'L6']
SPECIES = ['human', 'macaque', 'mouse']
SPECIES_PREFIX = {'human': 'human', 'macaque': 'monkey', 'mouse': 'mouse'}
SPECIES_COLORS = {'human': '#1f6f78', 'macaque': '#5ba3aa', 'mouse': '#e08a3c'}
DISPLAY_PROGRAMS = ['P19', 'P28', 'P41', 'P23']
PROGRAM_TITLES = {'P19': 'Axon guidance/neural crest', 'P28': 'Axon guid./adh. (L6b)', 'P41': 'Myelination', 'P23': 'Reg. myelination'}
(FS_TITLE, FS_LAB, FS_TICK, FS_SM) = (7.5, 6.5, 5.5, 5.2)
FS_CAV = 5.0
C_MOUSE = '#e08a3c'
C_GREY = '#9aa0a6'
plt.rcParams.update({'svg.fonttype': 'none', 'font.family': 'Liberation Sans', 'font.size': 6, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2, 'axes.edgecolor': '#444444', 'pdf.fonttype': 42, 'ps.fonttype': 42})

def _read_matrix(path: Path, required_layers: list[str]) -> pd.DataFrame:
    frame = _annotation_table(pd.read_csv(path, sep='\t', index_col=0))
    absent = [layer for layer in required_layers if layer not in frame.columns]
    if absent:
        raise ValueError()
    return frame.reindex(columns=LAYERS).apply(pd.to_numeric, errors='coerce')

def build_panel_b_source(project_root: Path) -> pd.DataFrame:
    aggregate = project_root / 'results/xspecies_humanmap_v1/spatial_xspecies/_aggregate'
    names_path = project_root / 'results/crossregion_v1/program_names.tsv'
    names = _annotation_table(pd.read_csv(names_path, sep='\t'))
    retained = names.loc[names['new_P'].astype(str).str.upper().ne('EXCLUDED')].copy()
    retained = retained.set_index('new_P')
    if not retained.index.is_unique:
        raise ValueError()
    matrices: dict[str, dict[str, pd.DataFrame]] = {}
    matrix_paths: dict[str, dict[str, Path]] = {}
    for species in SPECIES:
        prefix = SPECIES_PREFIX[species]
        paths = {'native_baseline_removed': aggregate / f'{prefix}_program_x_layer_baselinerm.tsv', 'raw_layer_summary': aggregate / f'{prefix}_program_x_layer.tsv', 'within_program_z_5layer': aggregate / f'{prefix}_program_x_layer_perprogz.tsv'}
        matrix_paths[species] = paths
        matrices[species] = {key: _read_matrix(path, COMMON_LAYERS if key == 'within_program_z_5layer' else LAYERS) for (key, path) in paths.items()}
    rows: list[dict[str, object]] = []
    for (display_order, program) in enumerate(DISPLAY_PROGRAMS, start=1):
        if program not in retained.index:
            raise ValueError()
        component = int(retained.loc[program, 'cnmf_component'])
        source_program = f'program_{component}'
        for (species_order, species) in enumerate(SPECIES, start=1):
            species_matrices = matrices[species]
            if any((source_program not in matrix.index for matrix in species_matrices.values())):
                raise ValueError()
            for (layer_order, layer) in enumerate(LAYERS, start=1):
                native_value = float(species_matrices['native_baseline_removed'].loc[source_program, layer])
                raw_value = float(species_matrices['raw_layer_summary'].loc[source_program, layer])
                z_value = float(species_matrices['within_program_z_5layer'].loc[source_program, layer])
                rows.append({'panel': 'b', 'display_order': display_order, 'display_program': program, 'cnmf_component': component, 'source_program_id': source_program, 'species': species, 'species_order': species_order, 'layer': layer, 'layer_order': layer_order, 'plot_metric': 'native_layer_median_baseline_removed_score', 'plotted_value': native_value, 'native_baseline_removed_value': native_value, 'raw_layer_summary_value': raw_value, 'within_program_z_5layer': z_value, 'plotted_source_path': str(matrix_paths[species]['native_baseline_removed'])})
    result = pd.DataFrame(rows)
    if len(result) != len(DISPLAY_PROGRAMS) * len(SPECIES) * len(LAYERS):
        raise RuntimeError()
    return result

def _coerce_bool(series: pd.Series, column: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    valid = normalized.isin({'true', 'false', '1', '0'})
    if not valid.all():
        bad = sorted(normalized.loc[~valid].unique())
        raise ValueError()
    return normalized.isin({'true', '1'})

def load_panel_d_source(area_stats: Path, names_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    area = _annotation_table(pd.read_csv(area_stats, sep='\t'))
    area = area.drop(columns=['name_short', 'confidence', 'display_label'], errors='ignore')
    required = {'new_P', 'area_median_sd', 'area_variability_rank', 'top14_display', 'human_region_variable', 'human_primary_flag'}
    absent = sorted(required - set(area.columns))
    if absent:
        raise ValueError()
    if len(area) != 54 or area['new_P'].nunique() != 54:
        raise ValueError()
    expected_programs = {f'P{i}' for i in range(1, 55)}
    if set(area['new_P']) != expected_programs:
        raise ValueError()
    for column in ('top14_display', 'human_region_variable', 'human_primary_flag'):
        area[column] = _coerce_bool(area[column], column)
    area['area_variability_rank'] = pd.to_numeric(area['area_variability_rank'], errors='raise').astype(int)
    if sorted(area['area_variability_rank'].tolist()) != list(range(1, 55)):
        raise ValueError()
    expected_top14 = area['area_variability_rank'].le(14)
    if not area['top14_display'].equals(expected_top14):
        raise ValueError()
    if not area.loc[area['human_primary_flag'], 'human_region_variable'].all():
        raise ValueError()
    if 'observation_unit' in area.columns:
        units = set(area['observation_unit'].dropna().astype(str))
        if units != {'section_area_median'}:
            raise ValueError()
    names = _annotation_table(pd.read_csv(names_path, sep='\t'))
    keep = [column for column in ('new_P', 'name_short', 'confidence') if column in names]
    names = names.loc[names['new_P'].astype(str).str.upper().ne('EXCLUDED'), keep]
    names = names.drop_duplicates('new_P')
    full = area.merge(names, on='new_P', how='left', validate='one_to_one')
    if full['name_short'].isna().any():
        absent = full.loc[full['name_short'].isna(), 'new_P'].tolist()
        raise ValueError()
    if 'confidence' not in full.columns:
        full['confidence'] = ''
    stars = np.where(full['confidence'].astype(str).str.strip().eq('class_b'), ' *', '')
    full['display_label'] = full['new_P'] + ' ' + full['name_short'].astype(str) + stars
    full['area_stats_source'] = str(area_stats)
    full = full.sort_values('area_variability_rank').reset_index(drop=True)
    display = full.loc[full['human_primary_flag']].copy()
    display = display.sort_values('area_variability_rank').reset_index(drop=True)
    return (full, display)

def _save_all(fig: plt.Figure, output_dir: Path, stem: str, *, close: bool=True) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    common = {'bbox_inches': 'tight', 'pad_inches': 0.01}
    fig.savefig(output_dir / f'{stem}.svg', format='svg', transparent=True, **common)
    fig.savefig(output_dir / f'{stem}.pdf', format='pdf', transparent=True, **common)
    fig.savefig(output_dir / f'{stem}.png', format='png', dpi=600, facecolor='white', transparent=False, **common)
    if close:
        plt.close(fig)

def _add_standalone_panel_letter(fig: plt.Figure, letter: str, *, y: float=0.992) -> None:
    fig.text(0.012, y, letter, ha='left', va='top', fontsize=13, fontweight='bold', fontfamily='Liberation Sans', color='black')

def render_panel_b(source: pd.DataFrame, output_dir: Path) -> None:
    (fig, axes) = plt.subplots(4, 1, figsize=(1.55, 4.05), sharey=True, sharex=True)
    x = np.arange(len(LAYERS))
    for (ax, program) in zip(axes.ravel(), DISPLAY_PROGRAMS):
        for species in SPECIES:
            values = source.loc[(source['display_program'] == program) & (source['species'] == species)].set_index('layer').reindex(LAYERS)['plotted_value'].to_numpy(float)
            linewidth = 1.0 if species == 'macaque' else 1.3
            alpha = 0.65 if species == 'macaque' else 0.95
            ax.plot(x, values, color=SPECIES_COLORS[species], lw=linewidth, alpha=alpha, marker='o', ms=1.8, mec='white', mew=0.3, zorder=2 if species == 'macaque' else 3)
        ax.axhline(0, color='#bbb', lw=0.5, ls='--', zorder=0)
        ax.axvspan(3.5, 5.5, color='#1f6f78', alpha=0.05, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels(LAYERS, fontsize=FS_CAV, rotation=0)
        ax.tick_params(axis='x', length=1.5, pad=1)
        ax.set_xlim(-0.4, 5.4)
        ax.set_title(f'{program} · {PROGRAM_TITLES[program]}', fontsize=FS_SM, fontweight='bold', pad=2, color='#222')
        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)
        ax.set_ylabel('laminar score\n(depth-corr.)', fontsize=FS_CAV)
        ax.tick_params(axis='y', length=1.5, labelsize=FS_CAV)
    handles = [Line2D([0], [0], color=SPECIES_COLORS[species], lw=1.5, label=species) for species in SPECIES]
    fig.legend(handles=handles, fontsize=FS_CAV, loc='upper center', bbox_to_anchor=(0.55, 1.005), ncol=3, frameon=False, handlelength=1.1, handletextpad=0.3, columnspacing=0.9)
    fig.suptitle('Laminar profiles across species', fontsize=5.6, fontweight='bold', y=1.052, x=0.55, linespacing=1.05)
    fig.subplots_adjust(left=0.205, right=0.97, top=0.915, bottom=0.045, hspace=0.55)
    _save_all(fig, output_dir, 'panel_b', close=False)
    _add_standalone_panel_letter(fig, 'b')
    _save_all(fig, output_dir, 'panel_b_labeled')

def panel_d_summary_text(n_top: int) -> str:
    return f'{n_top}/7 in mouse top-14\nall 7 ANOVA q < 0.05'

def render_panel_d(display: pd.DataFrame, output_dir: Path) -> None:
    (fig, ax) = plt.subplots(figsize=(3.0, 1.98))
    y = np.arange(len(display))[::-1]
    ranks = display['area_variability_rank'].to_numpy(int)
    top14 = ranks <= 14
    colors = [C_MOUSE if value else C_GREY for value in top14]
    ax.hlines(y, 0, ranks, color=colors, lw=1.3, alpha=0.8, zorder=1)
    ax.scatter(ranks, y, s=34, c=colors, edgecolor='white', linewidth=0.6, zorder=3)
    for (y_value, rank) in zip(y, ranks):
        ax.text(rank + 1.4, y_value, f'#{rank}', va='center', fontsize=FS_SM, color='#333')
    ax.axvline(14, color='#c0492f', lw=0.7, ls='--', zorder=0)
    ax.text(14, len(display) - 0.35, 'top-14\nmost area-variable', fontsize=FS_CAV, color='#c0492f', ha='center', va='bottom')
    ax.set_yticks(y)
    ax.set_yticklabels(display['display_label'], fontsize=FS_CAV)
    ax.set_xlabel('mouse areal-variability rank (1 = most variable, of 54)', fontsize=FS_SM)
    ax.set_xlim(0, 58)
    ax.set_ylim(-0.7, len(display) + 0.2)
    ax.set_title('Human region-variable programs in mouse areas', fontsize=FS_TITLE, fontweight='bold', pad=12)
    n_top = int(top14.sum())
    ax.text(0.985, 0.78, panel_d_summary_text(n_top), transform=ax.transAxes, fontsize=FS_CAV, ha='right', va='center', bbox={'boxstyle': 'round,pad=0.25', 'fc': 'white', 'ec': '#999', 'lw': 0.4})
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    fig.subplots_adjust(left=0.5, right=0.97, top=0.85, bottom=0.18)
    _save_all(fig, output_dir, 'panel_d', close=False)
    _add_standalone_panel_letter(fig, 'd', y=0.94)
    _save_all(fig, output_dir, 'panel_d_labeled')

def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--project-root', required=True, type=Path)
    parser.add_argument('--area-stats', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--source-data-dir', required=True, type=Path)
    args = parser.parse_args()
    args.source_data_dir.mkdir(parents=True, exist_ok=True)
    panel_b = build_panel_b_source(args.project_root)
    names_path = args.project_root / 'results/crossregion_v1/program_names.tsv'
    (panel_d_full, panel_d_display) = load_panel_d_source(args.area_stats, names_path)
    panel_b.to_csv(args.source_data_dir / 'fig7b_source_data.tsv', sep='\t', index=False)
    panel_d_full.to_csv(args.source_data_dir / 'fig7d_area_stats_retained54.tsv', sep='\t', index=False)
    panel_d_display.to_csv(args.source_data_dir / 'fig7d_display_primary7.tsv', sep='\t', index=False)
    render_panel_b(panel_b, args.output_dir)
    render_panel_d(panel_d_display, args.output_dir)
if __name__ == '__main__':
    main()
