#!/usr/bin/env python3
"""Donor-aware human regional and section-aware spatial validation.

Inference units are donor-region pseudobulks for snRNA-seq and spatial sections
(chips) for the spatial assay. Nuclei and bins are only used to construct those
aggregates. The script reads frozen sources and writes only to ``--outdir``.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrix
from statsmodels.genmod.cov_struct import Independence
from statsmodels.genmod.families import Gaussian
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests


PREFIX = "HUMAN_VALIDATION_"
BOOTSTRAPS = 1000
RNG_SEED = 20260802


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--chip-donor-lookup", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAPS)
    return parser.parse_args()


def outpath(outdir: Path, basename: str) -> Path:
    return outdir / f"{PREFIX}{basename}"




def pnum(program: str) -> int:
    return int(str(program).removeprefix("P"))


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    std = np.nanstd(values, ddof=0)
    if not np.isfinite(std) or std == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.nanmean(values)) / std


def build_locked_mapping(mapping: pd.DataFrame) -> pd.DataFrame:
    """Return the formal raw-component to display-P map, with strict checks."""
    required = {"old_P", "new_P", "status"}
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"formal mapping misses columns: {sorted(missing)}")
    locked = mapping.copy()
    locked["raw_component"] = pd.to_numeric(locked["old_P"], errors="raise").astype(int)
    display_values = locked["new_P"].astype(str).str.strip()
    locked["display_program"] = np.where(
        display_values.str.fullmatch(r"[0-9]+", na=False),
        "P" + display_values,
        display_values,
    )
    locked = locked.loc[
        locked["status"].astype(str).eq("kept")
        & locked["display_program"].str.fullmatch(r"P[0-9]+", na=False)
    ].copy()
    locked["display_number"] = locked["display_program"].map(pnum)
    locked = locked.sort_values("display_number").reset_index(drop=True)
    if len(locked) != 54 or locked["display_number"].tolist() != list(range(1, 55)):
        raise ValueError("formal map does not contain contiguous retained P1-P54")
    expected_raw = [i for i in range(1, 61) if i not in {9, 18, 19, 35, 52, 57}]
    if locked["raw_component"].tolist() != expected_raw:
        raise ValueError("formal map raw-component contract does not match retained 54 set")
    return locked


def design_record(name: str, formula: str, data: pd.DataFrame) -> dict[str, object]:
    matrix = dmatrix(formula, data, return_type="dataframe")
    rank = int(np.linalg.matrix_rank(matrix.to_numpy(dtype=float)))
    return {
        "design": name,
        "formula": formula,
        "n_rows": len(data),
        "n_columns": matrix.shape[1],
        "matrix_rank": rank,
        "full_rank": rank == matrix.shape[1],
    }


def partial_eta_from_anova(anova: pd.DataFrame, term: str) -> float:
    ss_term = float(anova.loc[term, "sum_sq"])
    ss_error = float(anova.loc["Residual", "sum_sq"])
    return ss_term / (ss_term + ss_error) if ss_term + ss_error > 0 else math.nan


def ols_region_test(data: pd.DataFrame, value: str) -> tuple[sm.regression.linear_model.RegressionResultsWrapper, dict[str, float]]:
    fit = smf.ols(f"{value} ~ C(region) + C(donor)", data=data).fit()
    table = anova_lm(fit, typ=2)
    term = "C(region)"
    return fit, {
        "F_region": float(table.loc[term, "F"]),
        "p_region": float(table.loc[term, "PR(>F)"]),
        "df_region": float(table.loc[term, "df"]),
        "df_resid": float(table.loc["Residual", "df"]),
        "partial_eta_sq": partial_eta_from_anova(table, term),
    }


def subclass_sensitivity_test(data: pd.DataFrame, value: str) -> dict[str, float]:
    formula = f"{value} ~ C(region) + C(donor) + C(subclass)"
    fit = smf.ols(formula, data=data).fit()
    term = "C(region)"
    term_columns = [name for name in fit.params.index if name.startswith("C(region)")]
    contrast = np.zeros((len(term_columns), len(fit.params)), dtype=float)
    for row, name in enumerate(term_columns):
        contrast[row, list(fit.params.index).index(name)] = 1.0
    clustered = fit.get_robustcov_results(cov_type="cluster", groups=data["donor_region"])
    wald = clustered.wald_test(contrast, scalar=True)
    anova = anova_lm(fit, typ=2)
    return {
        "wald_chi2_region": float(np.asarray(wald.statistic).squeeze()),
        "p_region_clustered": float(np.asarray(wald.pvalue).squeeze()),
        "df_region": float(len(term_columns)),
        "partial_eta_sq": partial_eta_from_anova(anova, term),
        "n_donor_region_clusters": int(data["donor_region"].nunique()),
        "n_rows_subclass_means": int(len(data)),
    }


def gee_term_test(
    data: pd.DataFrame, value: str, use_donor: bool
) -> tuple[sm.genmod.generalized_estimating_equations.GEEResultsWrapper, dict[str, float]]:
    formula = f"{value} ~ C(region) + C(majorDomain)"
    if use_donor:
        formula += " + C(donor)"
    fit = sm.GEE.from_formula(
        formula,
        groups="chip",
        data=data,
        family=Gaussian(),
        cov_struct=Independence(),
    ).fit(cov_type="robust")

    def wald(prefix: str) -> tuple[float, float, int]:
        names = [name for name in fit.params.index if name.startswith(prefix)]
        if not names:
            return math.nan, math.nan, 0
        contrast = np.zeros((len(names), len(fit.params)), dtype=float)
        for row, name in enumerate(names):
            contrast[row, list(fit.params.index).index(name)] = 1.0
        test = fit.wald_test(contrast, scalar=True)
        return (
            float(np.asarray(test.statistic).squeeze()),
            float(np.asarray(test.pvalue).squeeze()),
            len(names),
        )

    region_stat, region_p, region_df = wald("C(region)")
    layer_stat, layer_p, layer_df = wald("C(majorDomain)")
    return fit, {
        "region_wald_chi2": region_stat,
        "region_p": region_p,
        "region_df": region_df,
        "layer_wald_chi2": layer_stat,
        "layer_p": layer_p,
        "layer_df": layer_df,
    }


def balanced_predictions(
    fit: object, fixed_name: str, levels: list[str], template: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for level in levels:
        new = template.copy()
        new[fixed_name] = level
        rows.append((level, float(np.mean(fit.predict(new)))))
    output = pd.DataFrame(rows, columns=[fixed_name, "adjusted_mean"])
    output["standardized_effect"] = zscore(output["adjusted_mean"].to_numpy())
    return output


def fdr_column(frame: pd.DataFrame, p_column: str, q_column: str) -> pd.DataFrame:
    output = frame.copy()
    valid = output[p_column].notna() & np.isfinite(output[p_column])
    output[q_column] = math.nan
    if valid.any():
        output.loc[valid, q_column] = multipletests(output.loc[valid, p_column], method="fdr_bh")[1]
    return output


def read_snrna_obs(obs_path: Path) -> pd.DataFrame:
    metadata_columns = {"donor", "region", "subclass", "batch", "library_prep"}
    header = pd.read_csv(obs_path, nrows=0)
    index_column = header.columns[0]
    header_missing = metadata_columns.difference(header.columns)
    if header_missing:
        raise ValueError(f"snRNA obs metadata missing columns: {sorted(header_missing)}")
    obs = pd.read_csv(
        obs_path,
        index_col=index_column,
        usecols=[index_column, *sorted(metadata_columns)],
    )
    missing = metadata_columns.difference(obs.columns)
    if missing:
        raise ValueError(f"snRNA obs metadata missing columns: {sorted(missing)}")
    return obs


def load_snrna_pseudobulk(root: Path, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    crossregion = root / "results/crossregion_v1"
    cell_path = crossregion / "cell_program_region_subclass.parquet"
    obs_path = root / "inputs/snRNA_1M_obs.csv"
    all_programs = mapping["display_program"].tolist()
    raw_columns = [str(raw) for raw in mapping["raw_component"]]
    rename = dict(zip(raw_columns, all_programs))

    cell = pd.read_parquet(cell_path, columns=raw_columns + ["region", "subclass"])
    cell = cell.rename(columns=rename)
    obs = read_snrna_obs(obs_path)
    overlap = cell.index.intersection(obs.index)
    if len(overlap) != len(cell) or len(cell) != len(obs):
        raise ValueError(f"snRNA score/obs barcode mismatch: score={len(cell)}, obs={len(obs)}, overlap={len(overlap)}")
    joined = cell[all_programs].join(obs[["donor", "region", "subclass", "batch", "library_prep"]], how="inner")
    if joined["donor"].isna().any() or joined["region"].isna().any():
        raise ValueError("donor or region missing after exact snRNA join")
    donor_region = joined.groupby(["donor", "region", "batch"], as_index=False)[all_programs].mean()
    n_nuclei = joined.groupby(["donor", "region"], as_index=False).size().rename(columns={"size": "n_nuclei"})
    donor_region = donor_region.merge(n_nuclei, on=["donor", "region"], validate="one_to_one")
    subclass = joined.groupby(["donor", "region", "subclass"], as_index=False)[all_programs].mean()
    subclass["donor_region"] = subclass["donor"].astype(str) + "|" + subclass["region"].astype(str)
    return donor_region, subclass, joined[["donor", "region", "subclass", "batch", "library_prep"]]


def run_snrna_models(
    donor_region: pd.DataFrame, subclass: pd.DataFrame, mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    programs = mapping["display_program"].tolist()
    base_region = sorted(donor_region["region"].unique())[0]
    donors = sorted(donor_region["donor"].unique())
    template = pd.DataFrame({"donor": donors, "region": base_region})
    effect_rows: list[dict[str, object]] = []
    profile_rows: list[pd.DataFrame] = []
    sensitivity_rows: list[dict[str, object]] = []
    adjusted = donor_region[["donor", "region"]].copy()

    for program in programs:
        fit, stats = ols_region_test(donor_region, program)
        effect_rows.append({"program": program, "raw_component": int(mapping.loc[mapping.display_program.eq(program), "raw_component"].iat[0]), **stats})
        profile = balanced_predictions(fit, "region", sorted(donor_region.region.unique()), template)
        profile.insert(0, "program", program)
        support = donor_region.groupby("region").agg(
            n_donors_observed=("donor", "nunique"), n_donor_region=("donor", "size"), n_nuclei=("n_nuclei", "sum")
        ).reset_index()
        profile = profile.merge(support, on="region", validate="one_to_one")
        profile_rows.append(profile)
        sensitivity_rows.append({"program": program, **subclass_sensitivity_test(subclass, program)})

        nuisance = donor_region.copy()
        nuisance["region"] = base_region
        adjusted[program] = donor_region[program].to_numpy(dtype=float) - np.asarray(
            fit.predict(nuisance), dtype=float
        )

    effects = fdr_column(pd.DataFrame(effect_rows), "p_region", "p_region_bh_fdr")
    effects = effects.sort_values("p_region_bh_fdr", na_position="last").reset_index(drop=True)
    profiles = pd.concat(profile_rows, ignore_index=True)
    sensitivity = fdr_column(pd.DataFrame(sensitivity_rows), "p_region_clustered", "p_region_clustered_bh_fdr")
    return effects, profiles, sensitivity, adjusted


def aggregate_spatial_sections(root: Path, mapping: pd.DataFrame, lookup_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    crossregion = root / "results/crossregion_v1"
    meta_path = crossregion / "spatial_bin50_meta.parquet"
    score_path = crossregion / "spatial_bin50_program_score_SCT.parquet"
    programs = mapping["display_program"].tolist()
    raw_columns = [f"program_{raw}" for raw in mapping["raw_component"]]
    rename = dict(zip(raw_columns, programs))

    metadata = pd.read_parquet(meta_path, columns=["bin", "chip", "region", "majorDomain"]).set_index("bin")
    if metadata.index.has_duplicates:
        raise ValueError("spatial bin metadata has duplicate bin IDs")
    totals: dict[tuple[str, str, str], tuple[np.ndarray, int]] = {}
    parquet = pq.ParquetFile(score_path)
    processed = 0
    for batch in parquet.iter_batches(batch_size=200_000, columns=["bin"] + raw_columns):
        block = batch.to_pandas()
        block_meta = metadata.reindex(block["bin"])
        if block_meta["chip"].isna().any():
            raise ValueError("spatial score bin absent from metadata")
        work = pd.concat([block_meta.reset_index(drop=True), block[raw_columns].rename(columns=rename)], axis=1)
        grouped = work.groupby(["chip", "region", "majorDomain"], sort=False)[programs].agg(["sum", "count"])
        for key, values in grouped.iterrows():
            sums = values.xs("sum", level=1).to_numpy(dtype=float)
            count_values = values.xs("count", level=1).to_numpy(dtype=int)
            if len(set(count_values.tolist())) != 1:
                raise ValueError(f"incomplete score vector for spatial group {key}")
            count = int(count_values[0])
            if key in totals:
                previous_sum, previous_count = totals[key]
                totals[key] = (previous_sum + sums, previous_count + count)
            else:
                totals[key] = (sums, count)
        processed += len(block)

    if processed != parquet.metadata.num_rows:
        raise ValueError(f"spatial row count mismatch: {processed} vs {parquet.metadata.num_rows}")
    records = []
    for (chip, region, domain), (sums, count) in totals.items():
        records.append([chip, region, domain, count, *(sums / count)])
    section_layer = pd.DataFrame(records, columns=["chip", "region", "majorDomain", "n_bins", *programs])
    lookup = pd.read_csv(lookup_path, sep="\t", usecols=["chip", "region", "donor", "age", "sex"])
    if lookup.chip.duplicated().any():
        raise ValueError("chip-donor lookup has duplicate chip keys")
    section_layer = section_layer.merge(lookup, on=["chip", "region"], how="left", validate="many_to_one")
    if section_layer.donor.isna().any():
        missing = sorted(section_layer.loc[section_layer.donor.isna(), "chip"].unique())
        raise ValueError(f"chip donor lookup missing spatial sections: {missing}")
    input_sections = set(metadata.chip.unique())
    lookup_sections = set(lookup.chip)
    if input_sections != lookup_sections:
        raise ValueError(f"chip lookup coverage mismatch: input={len(input_sections)}, lookup={len(lookup_sections)}")
    section_counts = section_layer.groupby("region")["chip"].nunique().rename("n_sections").reset_index()
    section_counts["inference_status"] = np.where(
        section_counts.n_sections.ge(2), "inferential_section_support", "display_only_single_section"
    )
    section_layer = section_layer.merge(section_counts, on="region", validate="many_to_one")
    return section_layer, section_counts


def run_spatial_models(
    section_layer: pd.DataFrame, mapping: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, bool]:
    programs = mapping["display_program"].tolist()
    infer = section_layer.loc[section_layer.n_sections.ge(2)].copy()
    if infer.chip.nunique() < 20:
        raise ValueError("too few sections for spatial clustered model")
    full_formula = "1 + C(region) + C(majorDomain) + C(donor)"
    design = design_record("spatial_inferential", full_formula, infer)
    use_donor = bool(design["full_rank"])
    if not use_donor:
        reduced = design_record("spatial_inferential_without_donor", "1 + C(region) + C(majorDomain)", infer)
        if not reduced["full_rank"]:
            raise ValueError("spatial region/layer design is not identifiable")
    regions = sorted(infer.region.unique())
    domains = sorted(infer.majorDomain.unique())
    donors = sorted(infer.donor.unique())
    region_template = pd.DataFrame(list(itertools.product(donors, domains)), columns=["donor", "majorDomain"])
    region_template["region"] = regions[0]
    layer_template = pd.DataFrame(list(itertools.product(donors, regions)), columns=["donor", "region"])
    layer_template["majorDomain"] = domains[0]

    effect_rows: list[dict[str, object]] = []
    region_profiles: list[pd.DataFrame] = []
    layer_profiles: list[pd.DataFrame] = []
    adjusted = infer[["chip", "region", "donor", "majorDomain"]].copy()
    base_region = regions[0]
    for program in programs:
        fit, stats = gee_term_test(infer, program, use_donor=use_donor)
        effect_rows.append(
            {
                "program": program,
                "raw_component": int(mapping.loc[mapping.display_program.eq(program), "raw_component"].iat[0]),
                "n_sections": int(infer.chip.nunique()),
                "n_section_layer_observations": int(len(infer)),
                "n_regions_inferential": len(regions),
                "donor_adjustment_included": use_donor,
                **stats,
            }
        )
        region_profile = balanced_predictions(fit, "region", regions, region_template)
        region_profile.insert(0, "program", program)
        support = infer.groupby("region").agg(
            n_sections=("chip", "nunique"), n_donors=("donor", "nunique"), n_section_layer=("chip", "size")
        ).reset_index()
        region_profiles.append(region_profile.merge(support, on="region", validate="one_to_one"))
        layer_profile = balanced_predictions(fit, "majorDomain", domains, layer_template)
        layer_profile.insert(0, "program", program)
        layer_profiles.append(layer_profile)

        nuisance = infer.copy()
        nuisance["region"] = base_region
        adjusted[program] = infer[program].to_numpy(dtype=float) - np.asarray(
            fit.predict(nuisance), dtype=float
        )

    effects = pd.DataFrame(effect_rows)
    effects = fdr_column(effects, "region_p", "region_p_bh_fdr")
    effects = fdr_column(effects, "layer_p", "layer_p_bh_fdr")
    effects = effects.sort_values("region_p_bh_fdr", na_position="last").reset_index(drop=True)
    chip_adjusted = adjusted.groupby(["chip", "region", "donor"], as_index=False)[programs].mean()
    return effects, pd.concat(region_profiles, ignore_index=True), pd.concat(layer_profiles, ignore_index=True), chip_adjusted, use_donor


def write_methods(outdir: Path, spatial_donor_included: bool) -> None:
    text = f"""# Human donor-aware and spatial validation

## Evidence boundary

- snRNA-seq regional inference uses donor x region aggregates. Nuclei only form those aggregates and are not treated as biological replicates.
- Spatial regional and layer analyses use chip-by-majorDomain aggregates with chip as the repeated section unit. Bins only form those aggregates and are not treated as biological replicates.
- The formal component map is read from `results/crossregion_v1/program_renumber_map.tsv`; it is not inferred from display order. P33 is raw component 37.
- snRNA-seq primary model is `score ~ region + donor`. Batch/cohort is not added because it is exactly collinear with donor in this dataset. A subclass-conditioned, donor-region-clustered sensitivity analysis is reported separately.
- Spatial primary model is a Gaussian GEE with robust chip-cluster covariance: `score ~ region + majorDomain {'+ donor' if spatial_donor_included else ''}`. Regions represented by one section are excluded from this inferential model and preserved only in display-oriented tables.
"""
    outpath(outdir, "methods_and_evidence_boundary.md").write_text(text)


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    crossregion = root / "results/crossregion_v1"
    map_path = crossregion / "program_renumber_map.tsv"
    mapping_source = pd.read_csv(map_path, sep="\t")
    mapping = build_locked_mapping(mapping_source)
    mapping.to_csv(outpath(outdir, "component_mapping_locked.tsv"), sep="\t", index=False)

    donor_region, subclass, snrna_obs = load_snrna_pseudobulk(root, mapping)
    snrna_design = [
        design_record("snRNA_donor_adjusted", "1 + C(region) + C(donor)", donor_region),
        design_record("snRNA_with_batch", "1 + C(region) + C(donor) + C(batch)", donor_region),
    ]
    if not bool(snrna_design[0]["full_rank"]):
        raise ValueError("snRNA donor-adjusted region model is not identifiable")
    if int(snrna_design[1]["matrix_rank"]) != int(snrna_design[0]["matrix_rank"]):
        raise ValueError("unexpected batch contribution beyond donor; inspect cohort structure")
    snrna_effects, snrna_profiles, subclass_sensitivity, _ = run_snrna_models(donor_region, subclass, mapping)
    donor_region.to_csv(outpath(outdir, "snrna_donor_region_pseudobulk_all54.tsv"), sep="\t", index=False)
    snrna_effects.to_csv(outpath(outdir, "snrna_all54_donor_adjusted_region_effects.tsv"), sep="\t", index=False)
    snrna_profiles.to_csv(outpath(outdir, "snrna_all54_donor_adjusted_region_profiles.tsv"), sep="\t", index=False)
    subclass_sensitivity.to_csv(outpath(outdir, "snrna_all54_subclass_conditioned_sensitivity.tsv"), sep="\t", index=False)

    section_layer, spatial_counts = aggregate_spatial_sections(root, mapping, args.chip_donor_lookup.resolve())
    spatial_effects, spatial_profiles, layer_profiles, _, spatial_donor_included = run_spatial_models(section_layer, mapping)
    section_layer.to_csv(outpath(outdir, "spatial_section_by_layer_aggregates_all54.tsv"), sep="\t", index=False)
    spatial_counts.to_csv(outpath(outdir, "spatial_section_support_by_region.tsv"), sep="\t", index=False)
    spatial_effects.to_csv(outpath(outdir, "spatial_all54_region_and_layer_effects.tsv"), sep="\t", index=False)
    spatial_profiles.to_csv(outpath(outdir, "spatial_all54_region_profiles.tsv"), sep="\t", index=False)
    layer_profiles.to_csv(outpath(outdir, "spatial_all54_layer_profiles.tsv"), sep="\t", index=False)

    write_methods(outdir, spatial_donor_included)



if __name__ == "__main__":
    main()
