#!/usr/bin/env python3
"""Build Supplementary Table S1 from the fixed public resource-comparison record."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "existing_resource_axis",
    "representative_resources",
    "what_they_enable",
    "remaining_limitation",
    "what_this_study_adds",
    "source_note",
]

ROWS = [
    [
        "Cellular taxonomies and whole-brain cell references",
        "Adult human brain and cortical cell-reference studies",
        "Cell-type and cell-state taxonomies across brain regions.",
        "The primary coordinate is generally cell type, cluster, or dataset identity rather than a fixed adult human cortical program vocabulary.",
        "Defines 54 adult human cortical gene programs and relates each to cell preference, area, laminar context, tissue neighborhood, and cross-species recovery.",
        "Covered by the cited public source studies.",
    ],
    [
        "Spatial molecular architecture and cytoarchitectural maps",
        "Adult cortical spatial-transcriptomic and cytoarchitectural resources",
        "Placement of molecular states and cell organization in laminar and anatomical context.",
        "Spatial organization is not generally organized around one fixed program basis recovered from nuclei and carried through matched tissue.",
        "Projects one human program basis into matched spatial tissue and summarizes laminar domains and program neighborhoods.",
        "Covered by the cited public source studies.",
    ],
    [
        "Cross-species cellular references",
        "Human, macaque, and mouse cortical comparative resources",
        "Cross-species cell correspondence and species-specific cortical maps.",
        "Comparisons are usually cell-type indexed, species-specific, or de novo rather than fixed-basis projection of adult human programs.",
        "Tests macaque and mouse recovery, conservation, and divergence in a fixed human program basis.",
        "Covered by the cited public source studies.",
    ],
    [
        "Adult multi-area human cortical spatial resources",
        "Adult human cortical spatial and matched-tissue resources",
        "Adult human cortical spatial organization across regions and laminar contexts.",
        "They do not alone establish a reusable fixed program vocabulary with nuclei-to-tissue-to-species continuity.",
        "Combines multi-area nuclei, matched spatial tissue, program projection, and cross-species analyses in one vocabulary.",
        "Covered by the cited public source studies.",
    ],
    [
        "Disease and aging annotation layers",
        "Public disease- and aging-oriented gene-set resources",
        "External annotation of molecular results against disease and aging contexts.",
        "These resources are annotation contexts rather than adult cortical fixed-program coordinate systems.",
        "Maps public disease and aging annotations onto the same cortical program vocabulary after core cellular, spatial, and species analyses.",
        "Public source versions and checksums are recorded by the acquisition receipt workflow.",
    ],
    [
        "Broad multidimensional brain database integration",
        "scBrainScope",
        "Broad query across genes, gene sets, cell types, regions, development, disease, aging, and spatial datasets.",
        "Gene- and cell-level database overlays do not replace a controlled adult human cortical fixed-basis program resource.",
        "Provides a constrained adult human cortical program-coordinate layer with matched tissue and species projection.",
        "DOI 10.1093/nar/gkaf1092; https://www.brainscopes.org.",
    ],
    [
        "Neocortical compendium and projection resources",
        "NeMO Analytics",
        "Large neocortical data compendium and cross-modality projection tooling.",
        "The focus is developmental or compendium-level projection rather than an all-major-cell-class adult human cortical fixed-program resource with matched spatial recovery.",
        "Provides a complementary adult-cortex program-coordinate analysis across nuclei, matched tissue, macaque, mouse, spatial organization, aging, and disease.",
        "DOI 10.1038/s41593-026-02204-4; https://nemoanalytics.org/landing/neocortex/.",
    ],
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(FIELDS)
        writer.writerows(ROWS)
    print(f"Wrote {args.output_tsv}: {len(ROWS)} resource axes")


if __name__ == "__main__":
    main()
