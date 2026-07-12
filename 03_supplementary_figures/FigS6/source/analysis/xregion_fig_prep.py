#!/usr/bin/env python3
"""Prepare retained-hit network and spatial-field data for FigS6 e-j."""

import glob
import os

import pandas as pd


ROOT = os.environ.get("CORTEX_NMF_ROOT", "CORTEX_PROGRAM_ROOT")
B = f"{ROOT}/results/crossregion_v1"
OUT = os.environ.get("XREGION_OUTPUT_DIR", f"{B}/xregion_auroc")
RETAIN_MAP = os.environ.get("RETAIN_MAP", f"{ROOT}/tables/TableS3_program_annotation.tsv")
TOPN = 8


def main():
    mapping = pd.read_csv(RETAIN_MAP, sep="\t")
    assert len(mapping) == 54
    mapping["new_int"] = mapping["new_P"].astype(str).str.removeprefix("P").astype(int)
    mapping["old_int"] = mapping["cnmf_component"].astype(int)
    assert mapping["new_int"].tolist() == list(range(1, 55))
    new_to_old = dict(zip(mapping["new_int"], mapping["old_int"]))
    names = dict(zip(mapping["new_int"], mapping["functional_name"].astype(str)))

    hits = pd.read_csv(f"{OUT}/m4_rewiring_hits.tsv", sep="\t")
    assert len(hits) > 0, "strict funnel has no hits; do not synthesize network panels"
    hit_programs = hits["program"].astype(int).tolist()

    regions = sorted(
        os.path.basename(path).replace("m2b_coact_corr_", "").replace(".tsv", "")
        for path in glob.glob(f"{OUT}/m2b_coact_corr_*.tsv")
    )
    corr = {
        region: pd.read_csv(f"{OUT}/m2b_coact_corr_{region}.tsv", sep="\t", index_col=0)
        for region in regions
    }
    for matrix in corr.values():
        matrix.index = matrix.index.astype(int)
        matrix.columns = matrix.columns.astype(int)
        assert matrix.shape == (54, 54)

    edge_rows = []
    for hit in hit_programs:
        partner_union = set()
        for region in regions:
            partner_union.update(
                corr[region].loc[hit].drop(hit).nlargest(TOPN).index.astype(int)
            )
        for region in regions:
            for partner in sorted(partner_union):
                edge_rows.append(
                    (hit, partner, region, float(corr[region].loc[hit, partner]),
                     names[hit], names[partner])
                )
    edges = pd.DataFrame(
        edge_rows,
        columns=["hit", "partner", "region", "corr", "hit_name", "partner_name"],
    )
    edges.to_csv(f"{OUT}/fig_partner_network_edges.tsv", sep="\t", index=False)

    neigh = pd.read_csv(f"{OUT}/m2b_neigh_program_region_self_auroc.tsv", sep="\t")
    plan_rows = []
    slots = 1
    for hit in hit_programs:
        by_target = (
            neigh[neigh["program"].eq(hit)]
            .groupby("regionB")["self_auroc"]
            .mean()
            .sort_values()
        )
        reference = str(by_target.index[-1])
        for target in by_target.index[:3]:
            plan_rows.append(
                (slots, hit, new_to_old[hit], reference, str(target), str(target))
            )
            slots += 1
            if slots > 3:
                break
        if slots > 3:
            break
    assert len(plan_rows) == 3
    plan = pd.DataFrame(
        plan_rows,
        columns=["slot", "program", "old_component", "reference", "target", "field_region"],
    )
    plan.to_csv(f"{OUT}/fig_hit_panel_plan.tsv", sep="\t", index=False)

    old_components = sorted(plan["old_component"].unique())
    meta = pd.read_parquet(
        f"{B}/spatial_bin50_meta.parquet",
        columns=["bin", "chip", "x", "y", "region", "majorDomain"],
    )
    score_cols = ["bin"] + [f"program_{old}" for old in old_components]
    score = pd.read_parquet(f"{B}/spatial_bin_program_score.parquet", columns=score_cols)
    rename = {
        f"program_{old}": f"program_{new}"
        for new, old in new_to_old.items()
        if old in old_components
    }
    score = score.rename(columns=rename)
    fields = score.merge(meta, on="bin", how="inner")
    selected_regions = plan["field_region"].unique().tolist()
    fields = fields[fields["region"].isin(selected_regions)]
    chip_counts = (
        fields.groupby(["region", "chip"]).size().reset_index(name="n")
        .sort_values("n", ascending=False)
    )
    selected_chips = chip_counts.groupby("region").first().reset_index()[["region", "chip"]]
    fields = fields.merge(selected_chips, on=["region", "chip"], how="inner")
    assert set(fields["region"]) == set(selected_regions)
    fields.to_csv(f"{OUT}/fig_hit_spatial_fields.tsv", sep="\t", index=False)

    print(f"[prep] strict hits={hit_programs} network_rows={len(edges)}")
    print(plan.to_string(index=False))
    print(f"[prep] spatial rows={len(fields)} regions={selected_regions}")


if __name__ == "__main__":
    main()
