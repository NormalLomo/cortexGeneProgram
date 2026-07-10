#!/usr/bin/env python3
"""Calculate program enrichment for a public gene-set association table."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectra", type=Path, required=True, help="Programs by genes TSV.")
    parser.add_argument("--gene-sets", type=Path, required=True, help="TSV with set_id and gene columns.")
    parser.add_argument("--top-genes", type=int, default=150)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    spectra = pd.read_csv(args.spectra, sep="\t", index_col=0)
    gene_sets = pd.read_csv(args.gene_sets, sep="\t")
    if not {"set_id", "gene"}.issubset(gene_sets.columns):
        raise ValueError("Gene-set table requires set_id and gene columns")
    background = set(spectra.columns.astype(str))
    rows: list[dict[str, object]] = []
    for program, values in spectra.iterrows():
        selected = set(values.sort_values(ascending=False).head(args.top_genes).index.astype(str))
        for set_id, group in gene_sets.groupby("set_id"):
            members = set(group["gene"].astype(str)) & background
            overlap = len(selected & members)
            p_value = hypergeom.sf(overlap - 1, len(background), len(members), len(selected)) if overlap else 1.0
            rows.append({"program": program, "set_id": set_id, "overlap": overlap, "p_value": p_value})
    output = pd.DataFrame(rows).sort_values("p_value")
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote {args.output_tsv}")


if __name__ == "__main__":
    main()
