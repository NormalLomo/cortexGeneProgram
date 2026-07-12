#!/usr/bin/env python3
"""Validate the source-component to retained-program mapping used by Fig. 5."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import NamedTuple


EXCLUDED_OLD_IDS = (9, 18, 19, 35, 52, 57)
EXPECTED_SOURCE_PROGRAMS = 60
EXPECTED_RETAINED_PROGRAMS = 54


class ProgramContract(NamedTuple):
    old_to_new: dict[int, int]
    retained_old_ids: tuple[int, ...]
    excluded_old_ids: tuple[int, ...]


def _parse_id(value: str) -> int | None:
    normalized = value.strip().removeprefix("P")
    if not normalized or normalized.upper() == "EXCLUDED":
        return None
    return int(normalized)


def load_program_contract(path: str | Path) -> ProgramContract:
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"old_P", "new_P"}.issubset(rows[0]):
        raise ValueError("program map must contain old_P and new_P columns")

    old_to_new: dict[int, int] = {}
    excluded: list[int] = []
    source_ids: list[int] = []
    for row in rows:
        old_id = _parse_id(row["old_P"])
        if old_id is None:
            raise ValueError("old_P must identify every source component")
        source_ids.append(old_id)
        new_id = _parse_id(row["new_P"])
        if new_id is None:
            excluded.append(old_id)
        else:
            old_to_new[old_id] = new_id

    if sorted(source_ids) != list(range(1, EXPECTED_SOURCE_PROGRAMS + 1)):
        raise ValueError("program map must contain source IDs 1 through 60 exactly once")
    if tuple(sorted(excluded)) != EXCLUDED_OLD_IDS:
        raise ValueError(f"excluded old IDs must be {EXCLUDED_OLD_IDS}")
    if len(old_to_new) != EXPECTED_RETAINED_PROGRAMS:
        raise ValueError("program map must retain exactly 54 programs")
    if sorted(old_to_new.values()) != list(range(1, EXPECTED_RETAINED_PROGRAMS + 1)):
        raise ValueError("retained IDs must be unique and contiguous from 1 through 54")

    retained = tuple(old_id for old_id in source_ids if old_id in old_to_new)
    return ProgramContract(old_to_new, retained, tuple(sorted(excluded)))
