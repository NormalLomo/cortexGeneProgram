#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import os
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()

def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--macaque-h5ad', type=Path, required=True)
    parser.add_argument('--human-spectra', type=Path, required=True)
    parser.add_argument('--ortholog-tsv', type=Path, required=True)
    parser.add_argument('--region-column', required=True)
    parser.add_argument('--output-tsv', type=Path, required=True)
    args = parser.parse_args()
    adata = ad.read_h5ad(args.macaque_h5ad)
    if args.region_column not in adata.obs:
        raise KeyError()
    spectra = pd.read_csv(args.human_spectra, sep='\t', index_col=0)
    ortholog = pd.read_csv(args.ortholog_tsv, sep='\t')
    required = {'human_gene', 'macaque_gene', 'orthology_type'}
    if not required.issubset(ortholog.columns):
        raise ValueError()
    ortholog = ortholog.loc[ortholog['orthology_type'].eq('one_to_one')].drop_duplicates(['human_gene', 'macaque_gene'])
    mapping = ortholog.set_index('macaque_gene')['human_gene']
    macaque_genes = pd.Index(adata.var_names.astype(str))
    usable = macaque_genes[macaque_genes.isin(mapping.index)]
    human_genes = mapping.loc[usable].to_numpy()
    keep = pd.Index(human_genes).isin(spectra.columns)
    usable = usable[keep]
    human_genes = human_genes[keep]
    if len(usable) == 0:
        raise ValueError()
    matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    expression = matrix[:, macaque_genes.get_indexer(usable)]
    loading = spectra.loc[:, human_genes].to_numpy(dtype=float).T
    scores = expression @ loading
    scores = pd.DataFrame(np.asarray(scores), index=adata.obs[args.region_column].astype(str), columns=spectra.index.astype(str))
    summary = scores.groupby(level=0).mean().reset_index().rename(columns={'index': args.region_column})
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_tsv, sep='\t', index=False)
if __name__ == '__main__':
    main()
