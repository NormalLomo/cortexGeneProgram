#!/usr/bin/env python3
"""Project human spectra onto a species expression matrix through one-to-one orthologs."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-spectra", type=Path, required=True, help="Programs by human genes TSV.")
    parser.add_argument("--species-expression", type=Path, required=True, help="Samples by species genes TSV.")
    parser.add_argument("--ortholog-tsv", type=Path, required=True)
    parser.add_argument("--species-gene-column", required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    spectra = pd.read_csv(args.human_spectra, sep="\t", index_col=0)
    expression = pd.read_csv(args.species_expression, sep="\t", index_col=0)
    ortholog = pd.read_csv(args.ortholog_tsv, sep="\t")
    mapping = ortholog.loc[ortholog["orthology_type"].eq("one_to_one")].drop_duplicates(["human_gene", args.species_gene_column])
    mapping = mapping.set_index(args.species_gene_column)["human_gene"]
    species_genes = expression.columns.intersection(mapping.index)
    human_genes = mapping.loc[species_genes].to_numpy()
    keep = pd.Index(human_genes).isin(spectra.columns)
    species_genes = species_genes[keep]
    human_genes = human_genes[keep]
    if len(species_genes) == 0:
        raise ValueError("No one-to-one orthologs overlap both matrices")
    scores = expression.loc[:, species_genes].to_numpy() @ spectra.loc[:, human_genes].to_numpy().T
    result = pd.DataFrame(scores, index=expression.index, columns=spectra.index)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_tsv, sep="\t")
    print(f"Wrote {args.output_tsv}")


if __name__ == "__main__":
    main()
