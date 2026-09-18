#!/usr/bin/env python3
"""Regional donor analysis for the retained 54 human programs.

Saved as the actual production logic; NOT executed in the 2026-09-08 revision.
The primary association is the within-donor permutation F test (exchangeability
of regional observations within donor is assumed). Ordinary blocked F/p and
partial eta-squared are model summaries. Wild bootstrap and donor LOO are
sensitivity analyses. BH across TARGET8 is explicitly post-selection, not an
independent confirmatory family. LOO ranges are influence summaries, not CIs.
No CR omnibus, spatial analysis, or additional diagnostic artifacts are created.
"""
from __future__ import annotations
import argparse
import itertools
import math
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import scipy.stats as st
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

PROGRAMS = [f"P{i}" for i in range(1, 55)]
TARGET8 = ["P1", "P3", "P4", "P6", "P8", "P9", "P13", "P33"]
RNG_SEED = 20260802
WEBB_REPLICATES = 99_999
PERMUTATION_REPLICATES = 99_999
ROOT = Path(__import__("os").environ["NMF_WORK_ROOT"])

def bh(values: pd.Series) -> np.ndarray:
    array = values.to_numpy(dtype=float)
    output = np.full(len(array), np.nan)
    valid = np.isfinite(array)
    if valid.any():
        output[valid] = multipletests(array[valid], method="fdr_bh")[1]
    return output


def matrices(data: pd.DataFrame) -> dict[str, object]:
    x0_frame = patsy.dmatrix("1 + C(donor)", data, return_type="dataframe")
    x1_frame = patsy.dmatrix("1 + C(donor) + C(region)", data, return_type="dataframe")
    x0 = x0_frame.to_numpy(dtype=float)
    x1 = x1_frame.to_numpy(dtype=float)
    q0, _ = np.linalg.qr(x0, mode="reduced")
    q1, _ = np.linalg.qr(x1, mode="reduced")
    m0 = np.eye(len(data)) - q0 @ q0.T
    m1 = np.eye(len(data)) - q1 @ q1.T
    rank0 = int(np.linalg.matrix_rank(x0))
    rank1 = int(np.linalg.matrix_rank(x1))
    donors = list(pd.unique(data["donor"]))
    donor_index = np.array([donors.index(value) for value in data["donor"]], dtype=int)
    return {
        "x0": x0,
        "x1": x1,
        "q0": q0,
        "q1": q1,
        "m0": m0,
        "m1": m1,
        "difference": m0 - m1,
        "rank0": rank0,
        "rank1": rank1,
        "df_num": rank1 - rank0,
        "df_denom": len(data) - rank1,
        "donors": donors,
        "donor_index": donor_index,
    }


def blocked_ols(y: np.ndarray, design: dict[str, object]) -> dict[str, float]:
    m0 = design["m0"]
    m1 = design["m1"]
    rss0 = float(y @ m0 @ y)
    rss1 = float(y @ m1 @ y)
    df_num = int(design["df_num"])
    df_denom = int(design["df_denom"])
    f_stat = ((rss0 - rss1) / df_num) / (rss1 / df_denom)
    p_value = float(st.f.sf(f_stat, df_num, df_denom))
    partial_eta = (rss0 - rss1) / rss0 if rss0 > 0 else math.nan
    return {
        "blocked_ols_F": f_stat,
        "blocked_ols_df_num": df_num,
        "blocked_ols_df_denom": df_denom,
        "blocked_ols_p_descriptive": p_value,
        "partial_eta_sq": partial_eta,
        "rss_restricted": rss0,
        "rss_full": rss1,
    }


def exact_rademacher_p(y: np.ndarray, design: dict[str, object], f_observed: float) -> dict[str, object]:
    residual = design["m0"] @ y
    donor_index = design["donor_index"]
    donor_count = len(design["donors"])
    residual_blocks = np.column_stack([
        residual * (donor_index == cluster) for cluster in range(donor_count)
    ])
    informative = np.flatnonzero(np.linalg.norm(residual_blocks, axis=0) > 1e-12)
    weights_small = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(informative))), dtype=float)
    weights = np.ones((len(weights_small), donor_count), dtype=float)
    weights[:, informative] = weights_small
    a = residual_blocks.T @ design["difference"] @ residual_blocks
    b = residual_blocks.T @ design["m1"] @ residual_blocks
    numerator = np.einsum("bi,ij,bj->b", weights, a, weights)
    denominator = np.einsum("bi,ij,bj->b", weights, b, weights)
    valid = denominator > 1e-18
    f_boot = (numerator[valid] / design["df_num"]) / (denominator[valid] / design["df_denom"])
    exceed = int(np.sum(f_boot >= f_observed - 1e-12))
    return {
        "wild_rademacher_F_p": exceed / len(f_boot),
        "wild_rademacher_exceedances": exceed,
        "wild_rademacher_draws": int(len(f_boot)),
        "wild_rademacher_informative_clusters": int(len(informative)),
        "wild_rademacher_effective_unique_statistics_max": int(2 ** max(len(informative) - 1, 0)),
    }


def webb_p(
    y: np.ndarray,
    design: dict[str, object],
    f_observed: float,
    weights: np.ndarray,
) -> dict[str, object]:
    residual = design["m0"] @ y
    donor_index = design["donor_index"]
    donor_count = len(design["donors"])
    residual_blocks = np.column_stack([
        residual * (donor_index == cluster) for cluster in range(donor_count)
    ])
    a = residual_blocks.T @ design["difference"] @ residual_blocks
    b = residual_blocks.T @ design["m1"] @ residual_blocks
    numerator = np.einsum("bi,ij,bj->b", weights, a, weights)
    denominator = np.einsum("bi,ij,bj->b", weights, b, weights)
    valid = denominator > 1e-18
    f_boot = (numerator[valid] / design["df_num"]) / (denominator[valid] / design["df_denom"])
    exceed = int(np.sum(f_boot >= f_observed))
    return {
        "wild_webb_F_p": (exceed + 1) / (len(f_boot) + 1),
        "wild_webb_exceedances": exceed,
        "wild_webb_draws": int(len(f_boot)),
    }


def donor_block_permutation(
    data: pd.DataFrame,
    design: dict[str, object],
    observed_f: np.ndarray,
    replicates: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    values = data[PROGRAMS].to_numpy(dtype=float)
    groups = [np.flatnonzero(data["donor"].to_numpy() == donor) for donor in design["donors"]]
    q1 = design["q1"]
    rss0 = np.sum((design["m0"] @ values) ** 2, axis=0)
    exceedances = np.zeros(len(PROGRAMS), dtype=np.int64)
    completed = 0
    chunk_size = 1000
    while completed < replicates:
        chunk = min(chunk_size, replicates - completed)
        permuted = np.empty((chunk, len(data), len(PROGRAMS)), dtype=float)
        for indices in groups:
            if len(indices) == 1:
                permuted[:, indices, :] = values[indices][None, :, :]
            else:
                order = np.argsort(rng.random((chunk, len(indices))), axis=1)
                permuted[:, indices, :] = values[indices][order]
        norm_sq = np.sum(permuted * permuted, axis=1)
        projection = np.einsum("bnk,nq->bqk", permuted, q1, optimize=True)
        rss1 = norm_sq - np.sum(projection * projection, axis=1)
        f_values = ((rss0[None, :] - rss1) / design["df_num"]) / (rss1 / design["df_denom"])
        exceedances += np.sum(f_values >= observed_f[None, :], axis=0)
        completed += chunk
    p_values = (exceedances + 1) / (replicates + 1)
    return p_values, exceedances


def adjusted_region_profile(data: pd.DataFrame, program: str) -> tuple[pd.Series, dict[str, float]]:
    fit = smf.ols(f"{program} ~ C(region) + C(donor)", data=data).fit()
    regions = sorted(data["region"].unique())
    donors = sorted(data["donor"].unique())
    grid = pd.DataFrame(list(itertools.product(regions, donors)), columns=["region", "donor"])
    grid["prediction"] = fit.predict(grid)
    profile = grid.groupby("region")["prediction"].mean().reindex(regions)
    profile = profile - profile.mean()
    design = matrices(data)
    stats = blocked_ols(data[program].to_numpy(dtype=float), design)
    return profile, stats


def leave_one_donor_out(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_profiles: dict[str, pd.Series] = {}
    full_stats: dict[str, dict[str, float]] = {}
    for program in PROGRAMS:
        full_profiles[program], full_stats[program] = adjusted_region_profile(data, program)
    rows: list[dict[str, object]] = []
    for omitted in sorted(data["donor"].unique()):
        subset = data.loc[~data["donor"].eq(omitted)].copy()
        if subset["region"].nunique() != data["region"].nunique():
            raise ValueError(f"leaving out {omitted} removes a region level")
        for program in PROGRAMS:
            profile, stats = adjusted_region_profile(subset, program)
            full = full_profiles[program].reindex(profile.index)
            pearson = float(st.pearsonr(full, profile).statistic)
            spearman = float(st.spearmanr(full, profile).statistic)
            tolerance = max(float(np.max(np.abs(full))) * 1e-10, 1e-15)
            comparable = (np.abs(full) > tolerance) | (np.abs(profile) > tolerance)
            direction_agreement = float(np.mean(np.sign(full[comparable]) == np.sign(profile[comparable])))
            full_eta = full_stats[program]["partial_eta_sq"]
            rows.append({
                "program": program,
                "omitted_donor": omitted,
                "n_donors": int(subset["donor"].nunique()),
                "n_donor_region_observations": int(len(subset)),
                "n_regions": int(subset["region"].nunique()),
                "partial_eta_sq": stats["partial_eta_sq"],
                "full_partial_eta_sq": full_eta,
                "partial_eta_ratio_to_full": stats["partial_eta_sq"] / full_eta if full_eta > 0 else math.nan,
                "profile_pearson_r_to_full": pearson,
                "profile_spearman_rho_to_full": spearman,
                "centered_region_direction_agreement": direction_agreement,
                "top_region": profile.idxmax(),
                "top_region_matches_full": bool(profile.idxmax() == full.idxmax()),
                "bottom_region": profile.idxmin(),
                "bottom_region_matches_full": bool(profile.idxmin() == full.idxmin()),
            })
    loo = pd.DataFrame(rows)
    summary = loo.groupby("program", as_index=False).agg(
        loo_partial_eta_min=("partial_eta_sq", "min"),
        loo_partial_eta_median=("partial_eta_sq", "median"),
        loo_partial_eta_max=("partial_eta_sq", "max"),
        loo_partial_eta_ratio_min=("partial_eta_ratio_to_full", "min"),
        loo_profile_pearson_min=("profile_pearson_r_to_full", "min"),
        loo_profile_spearman_min=("profile_spearman_rho_to_full", "min"),
        loo_direction_agreement_min=("centered_region_direction_agreement", "min"),
        loo_top_region_match_fraction=("top_region_matches_full", "mean"),
        loo_bottom_region_match_fraction=("bottom_region_matches_full", "mean"),
    )
    return loo, summary



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "analysis/01_existing_regional_evidence/HUMAN_VALIDATION_snrna_donor_region_pseudobulk_all54.tsv")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Use a separately approved destination; this entry point was not run in the current revision.")
    args = parser.parse_args()
    donor = pd.read_csv(args.input, sep="\t")
    design = matrices(donor)
    rng = np.random.default_rng(RNG_SEED)
    webb_values = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5),
                            math.sqrt(0.5), 1.0, math.sqrt(1.5)])
    weights = rng.choice(webb_values, size=(WEBB_REPLICATES, len(design["donors"])))
    rows = []
    observed_f = np.empty(len(PROGRAMS), dtype=float)
    for index, program in enumerate(PROGRAMS):
        y = donor[program].to_numpy(dtype=float)
        ols = blocked_ols(y, design)
        observed_f[index] = ols["blocked_ols_F"]
        row = {"program": program, **ols}
        row.update(exact_rademacher_p(y, design, ols["blocked_ols_F"]))
        row.update(webb_p(y, design, ols["blocked_ols_F"], weights))
        rows.append(row)
    robust = pd.DataFrame(rows)
    perm_p, exceed = donor_block_permutation(donor, design, observed_f,
                                            PERMUTATION_REPLICATES, RNG_SEED)
    robust["donor_block_permutation_F_p"] = perm_p
    robust["donor_block_permutation_exceedances"] = exceed
    robust["donor_block_permutation_draws"] = PERMUTATION_REPLICATES
    target = robust["program"].isin(TARGET8)
    for prefix in ["wild_rademacher_F", "wild_webb_F", "donor_block_permutation_F"]:
        robust[prefix + "_q_all54"] = bh(robust[prefix + "_p"])
        robust[prefix + "_q_target8"] = np.nan
        robust.loc[target, prefix + "_q_target8"] = bh(robust.loc[target, prefix + "_p"])
        robust[prefix + "_df_num"] = int(design["df_num"])
        robust[prefix + "_df_denom_reference"] = int(design["df_denom"])
    loo, summary = leave_one_donor_out(donor)
    robust = robust.merge(summary, on="program", validate="one_to_one")
    # No conjunction of sensitivity methods defines the primary regional set.
    args.outdir.mkdir(parents=True, exist_ok=True)
    robust.to_csv(args.outdir / "DONOR_ROBUST_ALL54.tsv", sep="\t", index=False)
    loo.to_csv(args.outdir / "DONOR_LOO_STABILITY.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
