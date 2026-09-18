#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--chunksize", type=int, default=100000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(args.mapping, sep="\t")
    old_col = "old_P" if "old_P" in mapping.columns else "cnmf_component"
    mapping = mapping.loc[mapping["status"].astype(str).str.lower().eq("kept")].copy()
    mapping[old_col] = pd.to_numeric(mapping[old_col]).astype(int)
    mapping["new_program"] = mapping["new_P"].astype(str).str.replace(r"^P", "", regex=True).astype(int)
    mapping = mapping.sort_values("new_program")
    mapping["new_program"] = "P" + mapping["new_program"].astype(str)
    if mapping["new_program"].tolist() != [f"P{i}" for i in range(1, 55)]:
        raise ValueError("Retained mapping is not P1-P54")

    old_ids = mapping[old_col].astype(str).tolist()
    rename = dict(zip(old_ids, mapping["new_program"]))
    programs = mapping["new_program"].tolist()
    meta_columns = ["donor", "batch", "region", "subclass"]

    obs = pd.read_csv(
        args.obs,
        usecols=lambda c: c == "Unnamed: 0" or c in set(meta_columns),
        index_col=0,
    )
    if not obs.index.is_unique:
        raise ValueError("Observation identifiers are not unique")

    counts = None
    sums = None
    sums_sq = None
    matched = 0
    reader = pd.read_csv(
        args.scores,
        sep="\t",
        usecols=lambda c: c == "Unnamed: 0" or c in set(old_ids),
        index_col=0,
        chunksize=args.chunksize,
    )
    group_columns = ["donor", "batch", "region", "subclass"]
    for chunk in reader:
        if not chunk.index.is_unique:
            raise ValueError("Score identifiers are not unique within a chunk")
        meta = obs.reindex(chunk.index)
        if meta[meta_columns].isna().any(axis=None):
            missing = int(meta[meta_columns].isna().any(axis=1).sum())
            raise ValueError(f"{missing} score rows lack required metadata")
        chunk = chunk.rename(columns=rename)[programs].astype(float)
        joined = meta.join(chunk)
        grouped = joined.groupby(group_columns, observed=True)
        c = grouped.size().astype(np.int64)
        s = grouped[programs].sum()
        squares = joined[programs].pow(2)
        squares[group_columns] = joined[group_columns]
        ss = squares.groupby(group_columns, observed=True)[programs].sum()
        counts = c if counts is None else counts.add(c, fill_value=0)
        sums = s if sums is None else sums.add(s, fill_value=0)
        sums_sq = ss if sums_sq is None else sums_sq.add(ss, fill_value=0)
        matched += len(chunk)

    counts = counts.astype(np.int64).sort_index()
    sums = sums.reindex(counts.index).sort_index()
    sums_sq = sums_sq.reindex(counts.index).sort_index()
    means = sums.div(counts, axis=0)
    variances = (sums_sq - sums.pow(2).div(counts, axis=0)).div((counts - 1).clip(lower=1), axis=0)
    variances = variances.clip(lower=0)
    sds = np.sqrt(variances)
    sems = sds.div(np.sqrt(counts), axis=0)

    wide = means.copy()
    wide["n_cells"] = counts
    wide.reset_index().to_csv(outdir / "donor_region_subclass_means.tsv", sep="\t", index=False)

    long = []
    for program in programs:
        part = pd.DataFrame(
            {
                "n_cells": counts,
                "mean": means[program],
                "sd": sds[program],
                "sem": sems[program],
            }
        ).reset_index()
        part.insert(4, "program", program)
        long.append(part)
    long = pd.concat(long, ignore_index=True)
    long.to_csv(outdir / "donor_region_subclass_program_summary.tsv", sep="\t", index=False)

    count_frame = counts.reset_index(name="n_cells")
    coverage = count_frame.groupby("subclass", observed=True).agg(
        donor_region_combinations=("n_cells", "size"),
        combinations_ge10=("n_cells", lambda x: int((x >= 10).sum())),
        donors=("donor", "nunique"),
        regions=("region", "nunique"),
        cells=("n_cells", "sum"),
    ).reset_index()
    coverage.to_csv(outdir / "subclass_coverage.tsv", sep="\t", index=False)



if __name__ == "__main__":
    main()
