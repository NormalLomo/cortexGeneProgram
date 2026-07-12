#!/usr/bin/env python3
"""Convert raw 60-component inputs from the released summary to retained 54-program outputs for FigS5a."""

import argparse
import zipfile
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-zip", required=True)
    parser.add_argument("--retained-map", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mapping = pd.read_csv(args.retained_map, sep="\t")
    assert len(mapping) == 54
    mapping["new_int"] = mapping["new_P"].astype(str).str.removeprefix("P").astype(int)
    mapping["old_int"] = mapping["cnmf_component"].astype(int)
    assert mapping["new_int"].tolist() == list(range(1, 55))
    assert set(range(1, 61)) - set(mapping["old_int"]) == {9, 18, 19, 35, 52, 57}

    with zipfile.ZipFile(args.summary_zip) as archive:
        with archive.open("mat_program_subclass.tsv") as handle:
            long = pd.read_csv(handle, sep="\t")
    assert long.shape[0] == 22 * 60
    retained = long.merge(
        mapping[["old_int", "new_int"]],
        left_on="program",
        right_on="old_int",
        how="inner",
        validate="many_to_one",
    )
    matrix = retained.pivot(index="subclass", columns="new_int", values="mean")
    matrix = matrix.reindex(columns=range(1, 55))
    matrix.columns = matrix.columns.astype(str)
    assert matrix.shape == (22, 54)
    assert matrix.notna().all().all()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output, sep="\t")
    print(f"WROTE {output} shape={matrix.shape}")


if __name__ == "__main__":
    main()
