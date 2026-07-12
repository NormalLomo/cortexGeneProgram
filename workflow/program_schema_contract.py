"""Validate equivalence between the main and supplementary retained-program maps."""
from __future__ import annotations

import csv
from pathlib import Path


EXCLUDED_SOURCE_IDS = (9, 18, 19, 35, 52, 57)
SOURCE_PROGRAM_COUNT = 60
RETAINED_PROGRAM_COUNT = 54


def _parse_program_id(value: str) -> int | None:
    normalized = value.strip().removeprefix("P")
    if not normalized or normalized.upper() == "EXCLUDED":
        return None
    return int(normalized)


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_main_schema(path: str | Path) -> dict[int, int]:
    rows = _read_rows(path)
    if not rows or not {"old_P", "new_P"}.issubset(rows[0]):
        raise ValueError("main schema requires old_P and new_P columns")

    source_ids: list[int] = []
    mapping: dict[int, int] = {}
    excluded: list[int] = []
    for row in rows:
        old_id = _parse_program_id(row["old_P"])
        if old_id is None:
            raise ValueError("main schema has an invalid old_P value")
        source_ids.append(old_id)
        new_id = _parse_program_id(row["new_P"])
        if new_id is None:
            excluded.append(old_id)
        else:
            mapping[old_id] = new_id

    if sorted(source_ids) != list(range(1, SOURCE_PROGRAM_COUNT + 1)):
        raise ValueError("main schema must contain source IDs 1 through 60 exactly once")
    if tuple(sorted(excluded)) != EXCLUDED_SOURCE_IDS:
        raise ValueError(f"main schema exclusions must be {EXCLUDED_SOURCE_IDS}")
    return _validate_retained_mapping(mapping, "main schema")


def load_supplementary_schema(path: str | Path) -> dict[int, int]:
    rows = _read_rows(path)
    if not rows or not {"cnmf_component", "new_P"}.issubset(rows[0]):
        raise ValueError("supplementary schema requires cnmf_component and new_P columns")

    mapping: dict[int, int] = {}
    for row in rows:
        old_id = int(row["cnmf_component"])
        new_id = _parse_program_id(row["new_P"])
        if new_id is None:
            raise ValueError("supplementary schema cannot contain excluded rows")
        mapping[old_id] = new_id

    expected_old_ids = set(range(1, SOURCE_PROGRAM_COUNT + 1)) - set(EXCLUDED_SOURCE_IDS)
    if set(mapping) != expected_old_ids:
        raise ValueError("supplementary schema does not retain the expected 54 source IDs")
    return _validate_retained_mapping(mapping, "supplementary schema")


def _validate_retained_mapping(mapping: dict[int, int], name: str) -> dict[int, int]:
    if len(mapping) != RETAINED_PROGRAM_COUNT:
        raise ValueError(f"{name} must contain exactly 54 retained programs")
    if sorted(mapping.values()) != list(range(1, RETAINED_PROGRAM_COUNT + 1)):
        raise ValueError(f"{name} retained IDs must be unique and contiguous from 1 through 54")
    return mapping


def assert_equivalent_54_program_schemas(main_path: str | Path, supplementary_path: str | Path) -> int:
    main_mapping = load_main_schema(main_path)
    supplementary_mapping = load_supplementary_schema(supplementary_path)
    if main_mapping != supplementary_mapping:
        raise ValueError("main and supplementary retained-program schemas are not equivalent")
    return RETAINED_PROGRAM_COUNT
