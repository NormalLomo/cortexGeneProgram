#!/usr/bin/env python3
from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import os
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
SEED = 12345
N_PERM = 2000
SIG_FDR = 0.05
def bh_fdr(values: np.ndarray | pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or (not np.isfinite(p).all()):
        raise ValueError()
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty(len(p), dtype=float)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output

def cosine_rows(vector: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    numerator = matrix @ vector
    denominator = np.linalg.norm(matrix, axis=1) * np.linalg.norm(vector)
    with np.errstate(invalid='ignore', divide='ignore'):
        result = numerator / denominator
    result[~np.isfinite(result)] = np.nan
    return result

def program_number(value: Any) -> int:
    text = str(value)
    if text.startswith('P'):
        text = text[1:]
    return int(text)

def support_category(permutation_q: float, identity_gate: bool) -> str:
    permutation_supported = permutation_q < SIG_FDR
    if permutation_supported and identity_gate:
        return 'permutation_fdr_and_top_third'
    if permutation_supported:
        return 'permutation_fdr_only'
    if identity_gate:
        return 'top_third_only'
    return 'neither'

def load_universe(inputs: Path) -> tuple[pd.DataFrame, list[int], list[int]]:
    mapping = pd.read_csv(inputs / 'program_renumber_map.tsv', sep='\t')
    identity = pd.read_csv(inputs / 'TableS1_program_annotation.tsv', sep='\t')
    if len(mapping) != 60 or mapping['old_P'].nunique() != 60:
        raise ValueError()
    kept = mapping.loc[mapping['status'] == 'kept'].copy()
    excluded = set(mapping.loc[mapping['status'] != 'kept', 'old_P'].astype(int))
    if len(kept) != 54 or sorted(kept['new_P'].astype(int)) != list(range(1, 55)):
        raise ValueError()
    kept['old_P'] = kept['old_P'].astype(int)
    kept['new_P_int'] = kept['new_P'].astype(int)
    kept = kept.sort_values('new_P_int')
    identity = identity.copy()
    identity['new_P_int'] = identity['new_P'].map(program_number)
    identity['cnmf_component'] = identity['cnmf_component'].astype(int)
    identity = identity.sort_values('new_P_int')
    expected_release = list(range(1, 55))
    if len(identity) != 54 or identity['new_P_int'].tolist() != expected_release:
        raise ValueError()
    if kept['new_P_int'].tolist() != expected_release:
        raise ValueError()
    if kept['old_P'].tolist() != identity['cnmf_component'].tolist():
        raise ValueError()
    mapping_table = identity.merge(kept[['old_P', 'new_P_int']], left_on=['cnmf_component', 'new_P_int'], right_on=['old_P', 'new_P_int'], how='inner', validate='one_to_one')
    mapping_table = mapping_table.drop(columns=['old_P'])
    retained_sources = mapping_table['cnmf_component'].astype(int).tolist()
    return (mapping_table, retained_sources, sorted(excluded))

def load_npz(path: Path) -> tuple[np.ndarray, np.ndarray, list[int]]:
    archive = np.load(path, allow_pickle=True)
    human = archive['Hload'].astype(np.float64)
    refit = archive['refit'].astype(np.float64)
    programs = [program_number(value) for value in archive['programs']]
    if human.shape != refit.shape or human.shape[0] != len(programs):
        raise ValueError()
    if programs != list(range(1, 61)):
        raise ValueError()
    return (human, refit, programs)

def compute_species(human: np.ndarray, refit: np.ndarray, source_programs: list[int], universe: list[int], release_by_source: dict[int, str], species: str, rng: np.random.Generator | None=None, permutation_p_by_source: dict[int, float] | None=None) -> pd.DataFrame:
    index_by_source = {source: index for (index, source) in enumerate(source_programs)}
    selected_indices = [index_by_source[source] for source in universe]
    selected_human = human[selected_indices]
    rows: list[dict[str, Any]] = []
    for (diagonal_position, source) in enumerate(universe):
        index = index_by_source[source]
        comparisons = cosine_rows(refit[index], selected_human)
        observed = float(comparisons[diagonal_position])
        if not np.isfinite(observed):
            raise ValueError()
        off_diagonal = np.delete(comparisons, diagonal_position)
        if not np.isfinite(off_diagonal).all():
            raise ValueError()
        diag_rank = int(np.sum(comparisons >= observed))
        cross_program_p = (int(np.sum(off_diagonal >= observed)) + 1.0) / len(universe)
        percentile = 100.0 * float(np.mean(off_diagonal < observed))
        if permutation_p_by_source is None:
            if rng is None:
                raise ValueError()
            hi = human[index]
            ri = refit[index]
            denominator = np.linalg.norm(hi) * np.linalg.norm(ri) + 1e-300
            exceedances = 0
            for _ in range(N_PERM):
                permutation = rng.permutation(human.shape[1])
                null_cosine = float(ri[permutation] @ hi / denominator)
                if null_cosine >= observed:
                    exceedances += 1
            permutation_p = (exceedances + 1.0) / (N_PERM + 1.0)
        else:
            permutation_p = float(permutation_p_by_source[source])
            exceedances_float = permutation_p * (N_PERM + 1.0) - 1.0
            exceedances = int(round(exceedances_float))
            if not math.isclose(exceedances_float, exceedances, abs_tol=1e-09):
                raise ValueError()
        rows.append({'release_program_id': release_by_source.get(source, ''), 'source_program_id': f'P{source}', 'original_component': source, 'species': species, 'cosine': observed, 'permutation_exceedances': exceedances, 'permutation_p': permutation_p, 'cross_program_p': cross_program_p, 'cross_program_percentile': percentile, 'diag_rank': diag_rank, 'rank_family_n': len(universe), 'fdr_family_n': len(universe), 'n_permutations': N_PERM, 'seed': SEED})
    result = pd.DataFrame(rows)
    result['permutation_q'] = bh_fdr(result['permutation_p'])
    result['cross_program_q'] = bh_fdr(result['cross_program_p'])
    gate_top_n = max(1, len(result) // 3)
    result['identity_gate_top_n'] = gate_top_n
    result['identity_gate'] = result['diag_rank'] <= gate_top_n
    result['permutation_supported'] = result['permutation_q'] < SIG_FDR
    result['supported_by_permutation_and_rank'] = result['permutation_supported'] & result['identity_gate']
    result['display_category'] = [support_category(q_value, gate) for (q_value, gate) in zip(result['permutation_q'], result['identity_gate'])]
    result['display_rank1_circle'] = result['diag_rank'] == 1
    return result

def make_panel_source(new: pd.DataFrame, identity: pd.DataFrame, decay: pd.DataFrame) -> pd.DataFrame:
    identity_cols = ['new_P', 'cnmf_component', 'functional_name', 'dominant_class', 'dominant_subclass', 'confidence']
    base = identity[identity_cols].rename(columns={'new_P': 'release_program_id', 'cnmf_component': 'original_component'})
    decay = decay.copy()
    decay['original_component'] = decay['program'].map(program_number)
    base = base.merge(decay[['original_component', 'human_macaque', 'human_mouse', 'mouse_macaque']], on='original_component', how='left', validate='one_to_one')
    metric_cols = ['release_program_id', 'original_component', 'cosine', 'permutation_exceedances', 'permutation_p', 'permutation_q', 'cross_program_p', 'cross_program_q', 'cross_program_percentile', 'diag_rank', 'rank_family_n', 'fdr_family_n', 'identity_gate_top_n', 'identity_gate', 'permutation_supported', 'supported_by_permutation_and_rank', 'display_category', 'display_rank1_circle']
    for species in ['macaque', 'mouse']:
        block = new.loc[new['species'] == species, metric_cols].copy()
        block = block.rename(columns={column: f'{species}_{column}' for column in metric_cols[2:]})
        base = base.merge(block, on=['release_program_id', 'original_component'], how='left', validate='one_to_one')
    if not np.allclose(base['human_macaque'], base['macaque_cosine'], atol=1e-12):
        raise ValueError()
    if not np.allclose(base['human_mouse'], base['mouse_cosine'], atol=1e-12):
        raise ValueError()
    base = base.sort_values('macaque_cosine', ascending=False).reset_index(drop=True)
    base.insert(0, 'panel_order', np.arange(1, len(base) + 1))
    base['display_extreme_label'] = (base['panel_order'] <= 6) | (base['panel_order'] >= len(base) - 5)
    return base

def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    inputs = root / 'inputs'
    outputs = root / 'outputs'
    outputs.mkdir(parents=True, exist_ok=True)
    (identity, retained_sources, excluded_sources) = load_universe(inputs)
    release_by_source = dict(zip(identity['cnmf_component'], identity['new_P']))
    (macaque_h, macaque_r, macaque_programs) = load_npz(inputs / 'loadings_macaque.npz')
    (mouse_h, mouse_r, mouse_programs) = load_npz(inputs / 'loadings_mouse.npz')
    decay = pd.read_csv(inputs / 'decay_per_program_full.csv')
    if len(decay) != 60 or decay['program'].map(program_number).nunique() != 60:
        raise ValueError()
    rng = np.random.default_rng(SEED)
    new_macaque = compute_species(macaque_h, macaque_r, macaque_programs, retained_sources, release_by_source, 'macaque', rng=rng)
    new_mouse = compute_species(mouse_h, mouse_r, mouse_programs, retained_sources, release_by_source, 'mouse', rng=rng)
    new = pd.concat([new_macaque, new_mouse], ignore_index=True)
    panel_source = make_panel_source(new, identity, decay)
    new.to_csv(outputs / 'significance_direct54.tsv', sep='\t', index=False)
    panel_source.to_csv(outputs / 'panel_e_source_data.tsv', sep='\t', index=False)
    if set(new['original_component'].astype(int)) & set(excluded_sources):
        raise ValueError()
    if any((new.loc[new['species'] == species, 'fdr_family_n'].nunique() != 1 or new.loc[new['species'] == species, 'fdr_family_n'].iloc[0] != 54 for species in ['macaque', 'mouse'])):
        raise ValueError()
    if any((new.loc[new['species'] == species, 'rank_family_n'].nunique() != 1 or new.loc[new['species'] == species, 'rank_family_n'].iloc[0] != 54 for species in ['macaque', 'mouse'])):
        raise ValueError()
    if any((new.loc[new['species'] == species, 'identity_gate_top_n'].nunique() != 1 or new.loc[new['species'] == species, 'identity_gate_top_n'].iloc[0] != 18 for species in ['macaque', 'mouse'])):
        raise ValueError()
if __name__ == '__main__':
    main()
