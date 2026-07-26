#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path
from typing import NamedTuple
import os
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
EXPECTED_SOURCE_PROGRAMS = 60
EXPECTED_RETAINED_PROGRAMS = 54

class ProgramContract(NamedTuple):
    old_to_new: dict[int, int]
    retained_old_ids: tuple[int, ...]
    excluded_old_ids: tuple[int, ...]

def _parse_id(value: str) -> int | None:
    normalized = value.strip().removeprefix('P')
    if not normalized or normalized.upper() == 'EXCLUDED':
        return None
    return int(normalized)

def load_program_contract(path: str | Path) -> ProgramContract:
    path = Path(path)
    with path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or not {'old_P', 'new_P'}.issubset(rows[0]):
        raise ValueError()
    old_to_new: dict[int, int] = {}
    excluded: list[int] = []
    source_ids: list[int] = []
    for row in rows:
        old_id = _parse_id(row['old_P'])
        if old_id is None:
            raise ValueError()
        source_ids.append(old_id)
        new_id = _parse_id(row['new_P'])
        if new_id is None:
            excluded.append(old_id)
        else:
            old_to_new[old_id] = new_id
    if sorted(source_ids) != list(range(1, EXPECTED_SOURCE_PROGRAMS + 1)):
        raise ValueError()
    if len(old_to_new) != EXPECTED_RETAINED_PROGRAMS:
        raise ValueError()
    if sorted(old_to_new.values()) != list(range(1, EXPECTED_RETAINED_PROGRAMS + 1)):
        raise ValueError()
    retained = tuple((old_id for old_id in source_ids if old_id in old_to_new))
    return ProgramContract(old_to_new, retained, tuple(sorted(excluded)))
