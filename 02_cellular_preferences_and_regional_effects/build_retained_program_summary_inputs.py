#!/usr/bin/env python3
import argparse
import zipfile
from pathlib import Path
import pandas as pd
import os
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--summary-zip', required=True)
    parser.add_argument('--retained-map', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    mapping = pd.read_csv(args.retained_map, sep='\t')
    mapping['new_int'] = mapping['new_P'].astype(str).str.removeprefix('P').astype(int)
    mapping['old_int'] = mapping['cnmf_component'].astype(int)
    with zipfile.ZipFile(args.summary_zip) as archive:
        with archive.open('mat_program_subclass.tsv') as handle:
            long = pd.read_csv(handle, sep='\t')
    retained = long.merge(mapping[['old_int', 'new_int']], left_on='program', right_on='old_int', how='inner', validate='many_to_one')
    matrix = retained.pivot(index='subclass', columns='new_int', values='mean')
    matrix = matrix.reindex(columns=range(1, 55))
    matrix.columns = matrix.columns.astype(str)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output, sep='\t')
if __name__ == '__main__':
    main()
