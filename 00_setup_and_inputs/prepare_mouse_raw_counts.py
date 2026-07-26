#!/usr/bin/env python3
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
import os
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda : handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--input-h5ad', type=Path, required=True)
    parser.add_argument('--output-h5ad', type=Path, required=True)
    parser.add_argument('--expected-nuclei', type=int, default=761378)
    args = parser.parse_args()
    source = ad.read_h5ad(args.input_h5ad)
    if source.n_obs != args.expected_nuclei:
        raise ValueError()
    matrix = source.X.tocsr() if sparse.issparse(source.X) else sparse.csr_matrix(source.X)
    if np.any(matrix.data < 0) or not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ValueError()
    result = ad.AnnData(X=matrix.astype(np.float32), obs=source.obs.copy(), var=source.var.copy())
    result.uns['input_contract'] = 'Macosko Isocortex raw integer UMI subset'
    result.uns['source_sha256'] = sha256(args.input_h5ad)
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    result.write(args.output_h5ad)
if __name__ == '__main__':
    main()
