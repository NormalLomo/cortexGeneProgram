# Cortex gene program analysis code

This repository contains the analysis scripts for the human cortical gene program study. The current release defines human cortical gene programs from single-nucleus data, maps the fixed programs into human spatial transcriptomic sections, and relates their cellular, regional, spatial, disease, aging and functional contexts.

## Current human-only release

`human_current/` contains the current human-only analysis and figure source for
“Large-scale cellular and spatial census of human cortical gene programs”. Its
README defines the adopted P1–P54 program mapping, the single-task 100-start
cNMF public route, the fixed-H and spatial-score semantics, the environment
lanes and the boundaries of the recoverable historical source.

The current executable source is `human_current/`. It contains the adopted
P1–P54 mapping, the public 100-start cNMF route, fixed-H scoring,
donor-aware regional models, spatial relationships, external annotations and
current figure/table producers.

GO Biological Process enrichment uses each program's saved ranked genes and frozen GMT file. Report aggregation rebuilds the enrichment fields, while program naming and curated term choices remain in downstream naming steps.
