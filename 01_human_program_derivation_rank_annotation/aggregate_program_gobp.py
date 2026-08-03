#!/usr/bin/env python3
import os
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(os.environ.get("CORTEX_PROGRAM_ROOT", "/DATA/cortex_nmf_program")).expanduser().resolve()
REPORT_ROOT = Path(
    os.environ.get(
        "CORTEX_PROGRAM_GOBP_REPORT_ROOT",
        PROJECT_ROOT / "inputs" / "gobp_prerank_runs",
    )
).expanduser().resolve()
SOURCE = Path(
    os.environ.get(
        "CORTEX_PROGRAM_GOBP_SOURCE",
        PROJECT_ROOT / "results" / "crossregion_v1" / "program_annotation_gobp.tsv",
    )
).expanduser().resolve()
OUTPUT = Path(
    os.environ.get(
        "CORTEX_PROGRAM_GOBP_OUTPUT",
        PROJECT_ROOT / "results" / "crossregion_v1" / "program_annotation_gobp_rebuilt.tsv",
    )
).expanduser()


source = pd.read_csv(SOURCE, sep="\t", dtype=str).fillna("")
source["program"] = source["program"].astype(str)
source = source.set_index("program")
rows = []
for program in range(1, 61):
    report = pd.read_csv(REPORT_ROOT / f"P{program}" / "gseapy.gene_set.prerank.report.csv")
    report["NES"] = pd.to_numeric(report["NES"], errors="coerce")
    report["NOM p-val"] = pd.to_numeric(report["NOM p-val"], errors="coerce")
    report["FDR q-val"] = pd.to_numeric(report["FDR q-val"], errors="coerce")
    top = report.loc[report["NES"] > 0].sort_values("NES", ascending=False).head(3)
    first = top.iloc[0]
    source_row = source.loc[str(program)]
    term = str(first["Term"])
    match = re.search(r"\((GO:\d+)\)\s*$", term)
    rows.append(
        {
            "program": program,
            "top_BP_term": term,
            "term_id": match.group(1) if match else "",
            "NES_or_oddsratio": f"{float(first['NES']):.4f}",
            "pval": f"{float(first['NOM p-val']):.15g}",
            "fdr": f"{float(first['FDR q-val']):.15g}",
            "top3_BP_terms": " | ".join(top["Term"].astype(str)),
            "proposed_name": source_row["proposed_name"],
            "top10_loading_genes": source_row["top10_loading_genes"],
        }
    )
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(OUTPUT, sep="\t", index=False)
