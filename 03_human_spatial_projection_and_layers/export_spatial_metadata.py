import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd
import os
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get('CORTEX_PROGRAM_ROOT', '/DATA/cortex_nmf_program')).resolve()
base = str(PROJECT_ROOT / 'results/crossregion_v1/')
for f in ['spatial_bin50_meta.parquet', 'spatial_bin50_rctd_weights.parquet', 'spatial_bin50_program_score_SCT.parquet']:
    p = pq.ParquetFile(base + f)
meta = pq.read_table(base + 'spatial_bin50_meta.parquet').to_pandas()
g = meta.groupby('chip').agg(n=('x', 'size'), nlayers=('majorDomain', lambda s: s.nunique()))
g = g.sort_values('n', ascending=False)
meta.to_parquet(str(PROJECT_ROOT / 'scripts/fig2/_meta_cache.parquet'))
