import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import pandas as pd

base = "CORTEX_PROGRAM_ROOT/results/crossregion_v1/"

for f in ["spatial_bin50_meta.parquet",
          "spatial_bin50_rctd_weights.parquet",
          "spatial_bin50_program_score.parquet"]:
    p = pq.ParquetFile(base + f)
    print("====", f)
    print("rows:", p.metadata.num_rows, "cols:", p.metadata.num_columns)
    print("schema:", [field.name for field in p.schema_arrow])
    print()

# meta: distribution of chips, majorDomain, region
meta = pq.read_table(base + "spatial_bin50_meta.parquet").to_pandas()
print("META head:\n", meta.head())
print("\nchip counts (top):\n", meta["chip"].value_counts().head(10))
print("\nn chips:", meta["chip"].nunique())
print("\nmajorDomain values:\n", meta["majorDomain"].value_counts())
print("\ndomain values:\n", meta["domain"].value_counts())
print("\nregion values (top):\n", meta["region"].value_counts().head(20))
print("\nx range:", meta["x"].min(), meta["x"].max(), "y range:", meta["y"].min(), meta["y"].max())
print("\nNaN majorDomain:", meta["majorDomain"].isna().sum())

# per-chip bin counts to pick representative chip with all layers present
g = meta.groupby("chip").agg(n=("x","size"), nlayers=("majorDomain", lambda s: s.nunique()))
g = g.sort_values("n", ascending=False)
print("\nTop chips by size & layer coverage:\n", g.head(15))
meta.to_parquet("CORTEX_PROGRAM_ROOT/scripts/fig2/_meta_cache.parquet")
