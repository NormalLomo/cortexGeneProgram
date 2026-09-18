"""Retain the metadata cache consumed by the original spatial panel producers.

This is the input/output portion of the archived 01_probe.py; exploratory
printing and its reads of unused pre-SCT scores are not part of this entry.
"""
import os
from pathlib import Path
import pandas as pd

root = Path(os.environ[(__import__("os").environ["CORTEX_PROGRAM_ROOT"] + "")])
destination = root / "scripts/fig2/_meta_cache.parquet"
destination.parent.mkdir(parents=True, exist_ok=True)
pd.read_parquet(root / "results/crossregion_v1/spatial_bin50_meta.parquet").to_parquet(destination)
