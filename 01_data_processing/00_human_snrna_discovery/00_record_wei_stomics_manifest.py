#!/usr/bin/env python3
"""Record the exact public Wei STOMICS file listing without downloading source data."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_URL = "https://db.cngb.org/stomics/ajax/get_files/"
DATASET_ID = "STDS0000242"
EXPECTED_PROVIDER_FILES = 69
AGGREGATE_H5AD_ID = "STDF000000000022418"
AGGREGATE_H5AD_NAME = "snRNA.h5ad"
METADATA_ID = "STDF000000000022423"
METADATA_NAME = "SnRNA_Meta_with_EdLein_dataset_included.csv"
RDS_ID = "STDF000000000022422"
RDS_NAME = "SnRNA_seurat.RDS"


def api_page_url(page: int, per_page: int = 500) -> str:
    query = urlencode(
        {
            "filters": "[]",
            "pagination": json.dumps({"page": page, "per_page": per_page}, separators=(",", ":")),
        }
    )
    return f"{API_URL}?{query}"


def request_json(url: str) -> dict[str, object]:
    request = Request(url, headers={"User-Agent": "cortex-program-release-source-recorder/1.0"})
    with urlopen(request, timeout=60) as response:
        payload = json.load(response)
    if payload.get("code") != 0:
        raise RuntimeError(f"STOMICS API returned a nonzero status: {payload.get('code')}")
    return payload


def page_records(payload: dict[str, object]) -> tuple[list[dict[str, object]], int]:
    data = payload.get("data")
    if not isinstance(data, dict):
        raise ValueError("STOMICS API response is missing data")
    records = data.get("data")
    metadata = data.get("meta")
    if not isinstance(records, list) or not isinstance(metadata, dict):
        raise ValueError("STOMICS API response is missing file records or pagination metadata")
    pagination = metadata.get("pagination")
    if not isinstance(pagination, dict) or not isinstance(pagination.get("total_pages"), int):
        raise ValueError("STOMICS API response is missing total_pages")
    return [record for record in records if isinstance(record, dict)], pagination["total_pages"]


def belongs_to_dataset(record: dict[str, object]) -> bool:
    dataset_ids = record.get("dataset_ids", [])
    if isinstance(dataset_ids, str):
        dataset_ids = [dataset_ids]
    return isinstance(dataset_ids, list) and DATASET_ID in dataset_ids


def classify_provider_rows(records: list[dict[str, object]]) -> list[dict[str, str]]:
    if len(records) != EXPECTED_PROVIDER_FILES:
        raise ValueError(f"Expected {EXPECTED_PROVIDER_FILES} files for {DATASET_ID}, found {len(records)}")
    file_ids = [str(record.get("file_id", "")) for record in records]
    if len(file_ids) != len(set(file_ids)) or any(not file_id for file_id in file_ids):
        raise ValueError("STOMICS file listing must contain unique, nonempty file IDs")

    def matches(file_id: str, file_name: str) -> list[dict[str, object]]:
        return [
            record
            for record in records
            if str(record.get("file_id")) == file_id and str(record.get("file_name")) == file_name
        ]

    aggregate = matches(AGGREGATE_H5AD_ID, AGGREGATE_H5AD_NAME)
    metadata = matches(METADATA_ID, METADATA_NAME)
    rds = matches(RDS_ID, RDS_NAME)
    if len(aggregate) != 1 or len(metadata) != 1 or len(rds) != 1:
        raise ValueError("STOMICS listing no longer contains the required aggregate H5AD, metadata CSV, and RDS")

    spatial = [
        record
        for record in records
        if str(record.get("data_type", "")).lower() == "stomics"
        and str(record.get("file_type", "")).lower() == "h5ad"
    ]
    snrna_h5ad = [
        record
        for record in records
        if str(record.get("data_type", "")).lower() == "snrna"
        and str(record.get("file_type", "")).lower() == "h5ad"
    ]
    if len(spatial) != 44 or len(snrna_h5ad) != 23:
        raise ValueError("STOMICS listing must contain 44 spatial and 23 snRNA H5AD files")

    role_by_id = {str(record["file_id"]): "spatial_cellbin_h5ad" for record in spatial}
    role_by_id.update({str(record["file_id"]): "snrna_component_h5ad" for record in snrna_h5ad})
    role_by_id[AGGREGATE_H5AD_ID] = "human_cnmf_input"
    role_by_id[METADATA_ID] = "human_metadata"
    role_by_id[RDS_ID] = "snrna_reference_rds"
    if len(role_by_id) != EXPECTED_PROVIDER_FILES:
        raise ValueError("STOMICS listing contains files outside the declared public source contract")

    return [
        {
            "file_id": str(record["file_id"]),
            "file_name": str(record.get("file_name", "")),
            "file_path": str(record.get("file_path", "")),
            "file_size": str(record.get("file_size", "")),
            "data_type": str(record.get("data_type", "")),
            "provider_download": str(bool(record.get("download"))).lower(),
            "role": role_by_id[str(record["file_id"])],
        }
        for record in sorted(records, key=lambda item: str(item["file_id"]))
    ]


def fetch_dataset_rows() -> list[dict[str, object]]:
    first_records, total_pages = page_records(request_json(api_page_url(1)))
    records = first_records
    for page in range(2, total_pages + 1):
        page_rows, reported_pages = page_records(request_json(api_page_url(page)))
        if reported_pages != total_pages:
            raise RuntimeError("STOMICS pagination changed while recording the source listing")
        records.extend(page_rows)
    return [record for record in records if belongs_to_dataset(record)]


def write_rows(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = classify_provider_rows(fetch_dataset_rows())
    write_rows(rows, args.output_tsv)
    disabled = sum(row["provider_download"] == "false" for row in rows)
    print(f"Recorded {len(rows)} files for {DATASET_ID} in {args.output_tsv}")
    print(
        f"Provider download=false for {disabled}/{len(rows)} rows; this receipt does not claim local materialization.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
