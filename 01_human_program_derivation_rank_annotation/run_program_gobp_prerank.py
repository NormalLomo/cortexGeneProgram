#!/usr/bin/env python3
import os
from pathlib import Path

import gseapy
import pandas as pd


PROJECT_ROOT = Path(os.environ.get("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program")).expanduser().resolve()
INPUT_ROOT = PROJECT_ROOT / "inputs" / "gobp_prerank_runs"
OUTPUT_ROOT = Path(
    os.environ.get(
        "CORTEX_PROGRAM_GOBP_OUTPUT_ROOT",
        PROJECT_ROOT / "results" / "crossregion_v1" / "gobp_prerank_runs",
    )
).expanduser().resolve()


for program in range(1, 61):
    input_dir = INPUT_ROOT / f"P{program}"
    ranks = pd.read_csv(
        input_dir / "prerank_data.rnk",
        sep="\t",
        header=None,
        names=["gene_name", "ranking"],
    )
    output_dir = OUTPUT_ROOT / f"P{program}"
    output_dir.mkdir(parents=True, exist_ok=True)
    gseapy.prerank(
        rnk=ranks,
        gene_sets=str(input_dir / "gene_sets.gmt"),
        outdir=str(output_dir),
        min_size=10,
        max_size=1000,
        permutation_num=1000,
        seed=42,
        no_plot=True,
        verbose=True,
    )
