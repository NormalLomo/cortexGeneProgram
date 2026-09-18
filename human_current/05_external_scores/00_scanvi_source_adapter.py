#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import gc
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROJECT = Path(os.environ["NMF_SOURCE_ROOT"]) / "inputs/cortex_nmf_program"
ROOT = PROJECT / "R2_Benchmark_Staging/01_human_single_cell_benchmark"
METHOD_DIR = ROOT / "06_scanvi"
SOURCE_RDS = PROJECT / "original_data/human/single_cell/0_original/3_SnRNA_seurat_merged_1m_Cells.RDS"
RUNTIME_CARRIER = ROOT / "00_rds_raw_counts_carrier.h5ad"
COHORT = os.environ.get("SEAAD_COHORT", "normal")
if COHORT not in {"normal", "dementia"}:
    raise RuntimeError("SEAAD_COHORT must be normal or dementia")
PER_LIBRARY_ROOT = ROOT / "input/SEAAD" / COHORT
REFERENCE_N = 1_036_039
REFERENCE_GENE_N = 36_547
REGION_LIBRARY_N = {
    region: len(sorted(PER_LIBRARY_ROOT.joinpath(region).glob("library_*.rds")))
    for region in ("dlpfc", "mtg")
}
if any(value == 0 for value in REGION_LIBRARY_N.values()):
    raise RuntimeError(f"missing retained SEAAD {COHORT} library inputs")
REGION_CELL_N = ({
    "normal": {"dlpfc": 706_189, "mtg": 710_579},
    "dementia": {"dlpfc": 598_833, "mtg": 536_808},
})[COHORT]


def load_reference() -> ad.AnnData:
    reference = ad.read_h5ad(RUNTIME_CARRIER)
    if reference.n_obs != REFERENCE_N or reference.n_vars != REFERENCE_GENE_N:
        raise RuntimeError("reference dimensions differ from the approved source")
    if "gene" in reference.var.columns:
        reference.var_names = reference.var["gene"].astype(str).to_numpy()
    if not reference.var_names.is_unique:
        raise RuntimeError("reference gene axis is not unique")
    for key in ("donor", "region", "subclass"):
        if key not in reference.obs:
            raise RuntimeError(f"reference field missing: {key}")
    reference.obs["canonical22"] = reference.obs["subclass"].astype(str).to_numpy()
    reference.obs["donor"] = reference.obs["donor"].astype(str).to_numpy()
    reference.obs["region"] = reference.obs["region"].astype(str).to_numpy()
    return reference


def library_path(region: str, ordinal: int) -> Path:
    if region not in REGION_LIBRARY_N or not 1 <= ordinal <= REGION_LIBRARY_N[region]:
        raise ValueError("invalid region/library ordinal")
    return PER_LIBRARY_ROOT / region / f"library_{ordinal:05d}.rds"


def _r_vector(obj, name: str, dtype=str) -> np.ndarray:
    return np.asarray(obj.rx2(name), dtype=dtype)


def load_library(region: str, ordinal: int, reference_genes: np.ndarray) -> ad.AnnData:
    from rpy2 import rinterface as ri
    if not ri.embedded.isready():
        ri.initr()
    if "package:utils" not in tuple(str(x) for x in ri.baseenv["search"]()):
        ri.baseenv["library"](ri.StrSexpVector(["utils"]))
    from rpy2 import robjects as ro

    obj = ro.r["readRDS"](str(library_path(region, ordinal)))
    if str(obj.rx2("query_key")[0]) != region or int(obj.rx2("library_ordinal")[0]) != ordinal:
        raise RuntimeError("library identity differs from its approved path")
    if bool(obj.rx2("truth_read")[0]) or bool(obj.rx2("subchunked")[0]):
        raise RuntimeError("library input violates truth/subchunk contract")
    counts = obj.rx2("counts")
    dims = np.asarray(counts.slots["Dim"], dtype=np.int64)
    matrix = sp.csc_matrix(
        (
            np.asarray(counts.slots["x"]).copy(),
            np.asarray(counts.slots["i"], dtype=np.int32),
            np.asarray(counts.slots["p"], dtype=np.int64),
        ),
        shape=(int(dims[0]), int(dims[1])),
    ).T.tocsr()
    matrix.indices = matrix.indices.astype(np.int64, copy=False)
    matrix.indptr = matrix.indptr.astype(np.int64, copy=False)
    genes = np.asarray(list(counts.slots["Dimnames"])[0], dtype=str)
    reference_genes = np.asarray(reference_genes, dtype=str)
    if matrix.shape[1] != reference_genes.size or not np.array_equal(genes, reference_genes):
        raise RuntimeError("library model gene axis differs from reference")
    cell_ids = _r_vector(obj, "cell_ids")
    source_rows = _r_vector(obj, "source_rows", np.int64)
    cell_n = int(obj.rx2("cell_n")[0])
    specimen_id = str(obj.rx2("specimen_id")[0])
    donor_id = str(obj.rx2("donor_id")[0])
    region_name = str(obj.rx2("region")[0])
    if matrix.shape[0] != cell_n or cell_ids.size != cell_n or source_rows.size != cell_n:
        raise RuntimeError("library row count differs from approved metadata")
    obs = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "source_row": source_rows,
            "specimen_id": specimen_id,
            "donor": donor_id,
            "region": region_name,
            "cohort": f"seaad_{COHORT}_{region}",
            "canonical22": "UNKNOWN",
        },
        index=pd.Index(cell_ids, name="_index"),
    )
    del obj, counts
    gc.collect()
    return ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(reference_genes, name="_index")))



DONOR_LIBRARIES = {'dlpfc': {'H20.33.035': [1, 5, 80],
           'H21.33.038': [2, 73],
           'H21.33.014': [3, 37],
           'H21.33.037': [4, 89],
           'H21.33.030': [6, 82],
           'H21.33.015': [7, 54],
           'H20.33.013': [8, 63, 95],
           'H21.33.047': [9, 17],
           'H21.33.022': [10, 22],
           'H21.33.032': [11, 44, 61, 77],
           'H20.33.027': [12, 68, 71],
           'H20.33.025': [13, 64, 79],
           'H21.33.006': [14, 67],
           'H20.33.039': [15, 81, 84],
           'H20.33.002': [16, 62],
           'H21.33.019': [18, 26, 34],
           'H21.33.026': [19, 65, 99],
           'H19.33.004': [20, 24],
           'H21.33.028': [21, 42],
           'H20.33.005': [23, 90, 96],
           'H20.33.014': [25, 43],
           'H21.33.003': [27, 47, 86, 98],
           'H20.33.012': [28, 46, 75],
           'H21.33.036': [29, 76],
           'H21.33.035': [30, 57],
           'H20.33.034': [31, 32],
           'H20.33.044': [33, 45],
           'H20.33.019': [35, 53, 88],
           'H21.33.011': [36, 50],
           'H20.33.024': [38, 39],
           'H20.33.030': [40, 55],
           'H21.33.004': [41, 49, 51],
           'H21.33.025': [48, 83],
           'H21.33.023': [52, 78],
           'H21.33.040': [56, 72],
           'H20.33.036': [58, 93],
           'H20.33.032': [59, 87],
           'H20.33.001': [60, 69, 91],
           'H20.33.008': [66, 85, 97],
           'H21.33.041': [70, 94],
           'H21.33.033': [74, 92]},
 'mtg': {'H20.33.034': [1, 81],
         'H21.33.035': [2, 15],
         'H20.33.014': [3, 19, 49],
         'H21.33.015': [4, 16],
         'H21.33.036': [5, 54],
         'H20.33.025': [6, 18, 30],
         'H21.33.047': [7, 83],
         'H21.33.040': [8, 50],
         'H20.33.001': [9, 40, 102],
         'H21.33.026': [10, 27],
         'H20.33.002': [11, 55, 91],
         'H21.33.025': [12, 13],
         'H19.33.004': [14, 24, 63, 77],
         'H21.33.019': [17, 31, 72],
         'H21.33.022': [20, 76],
         'H21.33.004': [21, 34],
         'H21.33.011': [22, 42],
         'H20.33.044': [23, 32, 105],
         'H20.33.012': [25, 71, 80],
         'H21.33.041': [26, 70],
         'H21.33.006': [28, 39, 46],
         'H21.33.014': [29, 37],
         'H21.33.037': [33, 53],
         'H20.33.005': [35, 45, 60, 86, 90],
         'H200.1023': [36],
         'H20.33.039': [38, 43, 56],
         'H21.33.032': [41, 75],
         'H20.33.036': [44, 95, 104],
         'H20.33.030': [47, 64, 67],
         'H21.33.023': [48, 94],
         'H20.33.032': [51, 97],
         'H20.33.043': [52, 58, 79],
         'H20.33.035': [57, 68],
         'H21.33.028': [59, 85],
         'H20.33.013': [61, 73, 84],
         'H20.33.027': [62, 93, 103],
         'H21.33.003': [65, 88, 98],
         'H20.33.008': [66, 74],
         'H20.33.019': [69, 82, 96],
         'H20.33.024': [78, 89],
         'H21.33.038': [87],
         'H21.33.033': [92, 99],
         'H21.33.030': [100, 101]}}

def load_donor(region: str, donor_id: str, reference_genes: np.ndarray) -> ad.AnnData:
    if donor_id not in DONOR_LIBRARIES.get(region, {}):
        raise ValueError("donor is outside approved retained donor set")
    libraries = [load_library(region, i, reference_genes) for i in DONOR_LIBRARIES[region][donor_id]]
    if any(not x.var_names.equals(libraries[0].var_names) for x in libraries[1:]):
        raise RuntimeError("library feature axes are not identical")
    combined = ad.concat(libraries, axis=0, join="inner", merge="same", index_unique=None)
    if combined.obs_names.duplicated().any() or combined.obs["source_row"].duplicated().any():
        raise RuntimeError("donor cell/source_row identity is not unique")
    if combined.obs["donor"].nunique() != 1 or str(combined.obs["donor"].iloc[0]) != donor_id:
        raise RuntimeError("donor input identity differs from approved mapping")
    return combined


def load_region(region: str, reference_genes: np.ndarray) -> ad.AnnData:
    libraries = [load_library(region, i, reference_genes) for i in range(1, REGION_LIBRARY_N[region] + 1)]
    if any(not x.var_names.equals(libraries[0].var_names) for x in libraries[1:]):
        raise RuntimeError("library feature axes are not identical")
    combined = ad.concat(libraries, axis=0, join="inner", merge="same", index_unique=None)
    if combined.n_obs != REGION_CELL_N[region]:
        raise RuntimeError("regional retained-cell denominator differs from approval")
    if combined.obs_names.duplicated().any() or combined.obs["source_row"].duplicated().any():
        raise RuntimeError("regional cell/source_row identity is not unique")
    return combined
