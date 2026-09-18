#!/usr/bin/env python3
"""Build the audited K60 raw-component to retained-P1--P54 map."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EXCLUSION_COLUMNS = [
    "cohort_partial_eta2",
    "subclass_eta2",
    "depth_corr_log10nCount",
    "validity_flag",
    "redundancy_cluster",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-s2", type=Path, required=True)
    parser.add_argument("--table-s3", type=Path, required=True)
    parser.add_argument("--table-s4", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    s2 = pd.read_csv(args.table_s2, sep="\t")
    s3 = pd.read_csv(args.table_s3, sep="\t")
    s4 = pd.read_csv(args.table_s4, sep="\t")

    required_s4 = {"cnmf_component", "new_P", "status", "name_short", *EXCLUSION_COLUMNS}
    missing = sorted(required_s4.difference(s4.columns))
    if missing:
        raise ValueError(f"Table S4 is missing required columns: {missing}")
    if s4.shape[0] != 60 or sorted(s4["cnmf_component"].tolist()) != list(range(1, 61)):
        raise ValueError("Table S4 must define the complete raw K60 component universe 1..60.")

    retained = s4["status"].eq("kept")
    if retained.sum() != 54 or (~retained).sum() != 6:
        raise ValueError("Expected exactly 54 retained and 6 excluded raw K60 components.")
    retained_p = s4.loc[retained, "new_P"].tolist()
    if retained_p != [f"P{i}" for i in range(1, 55)]:
        raise ValueError("Retained components are not in the canonical P1--P54 order.")

    annotation = s3.rename(columns={"functional_name": "table_s3_functional_name"})[
        ["new_P", "dominant_class", "dominant_subclass", "confidence", "table_s3_functional_name"]
    ]
    # Excluded rows intentionally have new_P=NA, while the 54 retained IDs are unique.
    out = s4.merge(annotation, how="left", on="new_P", validate="many_to_one")
    out.insert(0, "raw_component", out.pop("cnmf_component"))
    out.insert(1, "final_analysis_set", out["status"].eq("kept").map({True: "retained_54", False: "excluded"}))
    out.insert(2, "retained_program_id", out["new_P"].where(out["status"].eq("kept"), pd.NA))
    out.insert(
        3,
        "exclusion_basis",
        out["status"].eq("kept").map(
            {
                True: "retained after cohort/technical validity review",
                False: "cohort-technical-excluded; operational numeric threshold unresolved",
            }
        ),
    )
    out = out.drop(columns=["new_P"])
    out.to_csv(args.output_dir / "canonical_program_map.tsv", sep="\t", index=False)

    # Preserve the distinct diagnostic grids instead of implying one shared K sweep.
    k_rows = [
        {
            "diagnostic_family": "native_cNMF_stability_and_error",
            "candidate_K": ";".join(map(str, [30, 40, 50, 60, 70, 80, 90])),
            "coverage_statement": "exact configured grid",
            "status": "certain",
            "source": "frozen-SCT release configuration and native cNMF scan",
            "K60_interpretation": "balanced working resolution within the broad 50-80 stable band; not a unique optimum",
        },
        {
            "diagnostic_family": "biological_sensitivity_TableS2",
            "candidate_K": ";".join(map(str, s2["K"].astype(str).str.extract(r"(\d+)", expand=False).astype(int).tolist())),
            "coverage_statement": "exact submitted Table S2 grid",
            "status": "certain",
            "source": "submitted Table S2 plus traced K80 extension output",
            "K60_interpretation": "regional-variability conclusions remain concordant at tested ranks",
        },
        {
            "diagnostic_family": "information_criterion_parsimony",
            "candidate_K": "20-200",
            "coverage_statement": "documented range only; exact vector unavailable",
            "status": "unresolved",
            "source": "submitted manuscript and Supplementary Data 2 source",
            "K60_interpretation": "criteria disagree and do not identify a unique optimum",
        },
    ]
    pd.DataFrame(k_rows).to_csv(args.output_dir / "k_grid_reconciled.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
