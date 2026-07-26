from __future__ import annotations
from itertools import combinations, combinations_with_replacement
from statistics import NormalDist
from typing import Iterable
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import os
from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program'))
COMMON_LAYERS = ('L1', 'L2', 'L4', 'L5', 'L6')

def normalize_layer(species: str, value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    species = species.lower()
    layer = str(value).strip().lower()
    mapping = {'human': {'l1': 'L1', 'l2': 'L2', 'l4': 'L4', 'l5': 'L5', 'l6': 'L6'}, 'macaque': {'l1': 'L1', 'l2': 'L2', 'l4': 'L4', 'l5': 'L5', 'l6': 'L6'}, 'mouse': {'1': 'L1', '2': 'L2', '2/3': 'L2', '4': 'L4', '5': 'L5', '6': 'L6', '6a': 'L6'}}
    if species not in mapping:
        raise ValueError()
    return mapping[species].get(layer)

def classify_laminar(human_preference: float, macaque_preference: float, mouse_preference: float, tau: float=0.2) -> tuple[str, str]:
    values = np.asarray([human_preference, macaque_preference, mouse_preference], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError()
    signs = np.sign(values).astype(int)
    same_nonzero_sign = np.all(signs == signs[0]) and signs[0] != 0
    if same_nonzero_sign and np.all(np.abs(values) >= tau):
        return ('A', '')
    if same_nonzero_sign:
        return ('B', '')
    if signs[1] == signs[2] and signs[0] != signs[1]:
        divergence = 'human-specific'
    elif signs[0] == signs[2] and signs[1] != signs[0]:
        divergence = 'macaque-only'
    elif signs[0] == signs[1] and signs[2] != signs[0]:
        divergence = 'mouse-only'
    else:
        divergence = 'three-way'
    return ('C', divergence)

def retained_program_map(names: pd.DataFrame) -> pd.DataFrame:
    required = {'new_P', 'cnmf_component', 'name_short', 'confidence'}
    absent = required.difference(names.columns)
    if absent:
        raise ValueError()
    out = names.loc[names['new_P'].astype(str).str.upper() != 'EXCLUDED'].copy()
    out['program_number'] = out['new_P'].astype(str).str.removeprefix('P').astype(int)
    out['new_P'] = 'P' + out['program_number'].astype(str)
    out['cnmf_component'] = out['cnmf_component'].astype(int)
    out = out.sort_values('program_number').reset_index(drop=True)
    if out['new_P'].duplicated().any() or out['cnmf_component'].duplicated().any():
        raise ValueError()
    expected = list(range(1, len(out) + 1))
    if out['program_number'].tolist() != expected:
        raise ValueError()
    out['score_column'] = 'program_' + out['cnmf_component'].astype(str)
    out['display_name'] = out['new_P'].astype(str) + ' ' + out['name_short'].astype(str)
    return out

def build_pair_index(programs: Iterable[str], include_diagonal: bool) -> pd.DataFrame:
    values = list(programs)
    iterator = combinations_with_replacement(values, 2) if include_diagonal else combinations(values, 2)
    rows = []
    for (a, b) in iterator:
        rows.append({'pair_id': f'{a}__{b}', 'program_a': a, 'program_b': b, 'is_diagonal': a == b})
    return pd.DataFrame(rows)

def partial_spearman_matrix(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    covariate = np.asarray(covariate, dtype=float)
    if values.ndim != 2 or covariate.ndim != 1 or values.shape[0] != covariate.size:
        raise ValueError()
    finite = np.isfinite(covariate) & np.isfinite(values).all(axis=1)
    values = values[finite]
    covariate = covariate[finite]
    if values.shape[0] < 4:
        return np.full((values.shape[1], values.shape[1]), np.nan)
    ranked_values = rankdata(values, axis=0, method='average')
    ranked_covariate = rankdata(covariate, method='average')
    design = np.column_stack([np.ones(ranked_covariate.size), ranked_covariate])
    residuals = ranked_values - design @ np.linalg.lstsq(design, ranked_values, rcond=None)[0]
    corr = np.corrcoef(residuals, rowvar=False)
    if corr.ndim == 0:
        corr = np.array([[float(corr)]])
    corr = np.asarray(corr, dtype=float)
    np.fill_diagonal(corr, 1.0)
    return corr

def wilson_interval(successes: int, total: int, confidence: float=0.95) -> tuple[float, float]:
    if total <= 0 or not 0 <= successes <= total:
        return (float('nan'), float('nan'))
    z = NormalDist().inv_cdf(0.5 + confidence / 2.0)
    p = successes / total
    denom = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denom
    half = z * np.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denom
    return (max(0.0, center - half), min(1.0, center + half))

def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full(pvalues.shape, np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    p = pvalues[finite]
    if p.size == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * p.size / np.arange(1, p.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out[finite] = restored
    return out
