# Cortex gene program analysis code

This repository contains the analysis scripts for the human cortical gene program study. The analyses use single nucleus and spatial transcriptomic data from human, macaque and mouse cortex to define human gene programs, examine their cellular and tissue organization, compare their molecular and spatial patterns across species, and relate selected programs to external disease, aging and functional data.

## Current human-only release

`human_current/` contains the current human-only analysis and figure source for
“Large-scale cellular and spatial census of human cortical gene programs”. Its
README defines the adopted P1–P54 program mapping, the single-task 100-start
cNMF public route, the fixed-H and spatial-score semantics, the environment
lanes and the boundaries of the recoverable historical source.

The numbered directories below are retained historical project stages. They
include earlier cross-species work and are not the entry point for the current
human-only manuscript.

## Methods stages

- `00_setup_and_inputs`: source data preparation.
- `01_human_program_derivation_rank_annotation`: human program derivation, rank assessment and annotation.
- `02_cellular_preferences_and_regional_effects`: cellular preferences and regional effects.
- `03_human_spatial_projection_and_layers`: human spatial projection and laminar organization.
- `04_local_subclass_program_relationships`: local relationships between cortical subclasses and programs.
- `05_program_program_spatial_associations`: local relationships among programs.
- `06_cross_species_nuclear_projection`: projection into macaque and mouse nuclei.
- `07_cross_species_spatial_laminar_areal`: cross species spatial, laminar and areal comparisons.
- `08_external_annotations`: disease, aging and functional annotations.

The scripts are organized by Methods stage. `Fig1A.png` summarizes the study design.

GO Biological Process enrichment uses each program's saved ranked genes and frozen GMT file. Report aggregation rebuilds the enrichment fields, while program naming and curated term choices remain in downstream naming steps.
