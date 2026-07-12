# Human Cortical Gene Programs

This repository contains the analysis code for human cortical gene programs. It accepts externally obtained source files, produces derived analysis outputs, and builds numbered figures, tables, and supplementary data. Scientific data payloads are not bundled with the source code.

## Research Background

Human cortical areas share many cellular programs while differing in how those programs are coordinated across regions, layers, spatial neighborhoods, species, and disease-related contexts. This code archive documents the analysis path used to identify and evaluate those coordinated cortical gene programs across single-nucleus and spatial datasets.

## Fig. 1A

Fig. 1A summarizes the analysis workflow, from externally acquired single-nucleus and spatial inputs through cNMF discovery, program annotation, spatial validation, and figure assembly. It is a guide to the executable stages below rather than a replacement for the input contracts in `DATA_SOURCES.tsv` and `workflow/PIPELINE_DAG.tsv`.

![Fig. 1A workflow overview](assets/Fig1A_workflow_overview.png)

## Setup

```bash
conda env create -f environment/environment.yml
conda activate cortex-gene-program
Rscript environment/install_spatial_r_packages.R
cp config/paths.example.env config/paths.env
source environment/activate_paths.sh config/paths.env
bash environment/bootstrap_legacy_path_aliases.sh
```

Set `CORTEX_PROGRAM_DATA_ROOT` to the directory holding externally obtained inputs and `CORTEX_PROGRAM_RESULTS_ROOT` to an empty results directory. Set `CORTEX_PROGRAM_CANONICAL_ROOT` only when running Fig. 2, Fig. 7, or Fig. 8 source entrypoints; it is the analysis root containing the required `results/`, `data/`, and `figures/` paths and can be overridden by `--canonical-root`. `environment/install_spatial_r_packages.R` installs SeuratDisk and spacexr from the exact commits stated in that script. `bootstrap_legacy_path_aliases.sh` creates the generic path aliases used by the analysis scripts. Official records, access routes, input contracts, and current acquisition status are listed in `DATA_SOURCES.tsv`.

## Execution Order

1. Run the Jorstad downloader to retrieve and verify all five collection assets. Run `00_record_wei_stomics_manifest.py --output-tsv WEI_PROVIDER_FILES.tsv` to record the exact 69-file Wei provider contract; it identifies the aggregate `snRNA.h5ad` rather than merging the 23 cell-type H5AD records. `download=false` means manual provider access is required, not that the acquisition code is missing. Do not substitute an arbitrary single H5AD for the Jorstad manifest.
2. Build the human raw-UMI input with `01_data_processing/00_human_snrna_discovery/01_prepare_jorstad_wei.py`. For Macosko, `00_fetch_macosko_nemo_bag.sh --output-dir "$CORTEX_PROGRAM_DATA_ROOT/macosko_raw" --resolve-fetch` materializes and validates the public 25 TB Raw_10x BDBag, but it does not provide the reconstructed H5AD and metadata required by `00_fix_subset_mouse.py`; that 761,378-nucleus edge remains explicitly blocked.
3. Run human cNMF with `python 01_data_processing/00_human_snrna_discovery/02_run_cnmf_discovery.py --config config/cnmf_discovery.yaml`. The workflow performs exactly 100 replicates for every configured K, then exports the K=60 package-native spectra and consensus usages to `results/cnmf_human_k60/`.
4. Run spatial preprocessing in this order: `00_bin_spatial_counts.py`, `01_convert_h5ad_to_h5seurat.R`, `01_sct_program_scoring.R`, `02_prepare_rctd_objects.R`, and `02_run_rctd.R`. The H5Seurat conversion and RCTD object preparation make the public bin50-to-RCTD boundary explicit. Han section identifiers remain provider-selected; validate a user-acquired nonempty list with `03_validate_han_sections.py`. Chen snRNA files are obtained through manual provider access and can be checked with `00_validate_chen_provider_files.py`.
5. Run cross-species projection, association analyses, and null models in their numbered directories under `01_data_processing/`.
6. Use `python run_release.py --list` or `python run_release.py --asset Fig1` to inspect each asset-specific chain. `run_release.py` is a catalog and validator; it does not synthesize figures or copy arbitrary tables.

The final asset map is `final_asset_map.tsv`; it is the only numbering authority for Fig1-Fig8, FigS1-FigS10, TableS1-TableS6, and SuppData1-SuppData6. In particular, FigS9 is disease sensitivity and FigS10 is cognition association. SuppData1, SuppData2, SuppData4, and SuppData6 include their renderer and finalization step. Rendered bytes can vary across software environments.

## Figure Input Boundaries

Fig. 5 uses a 60-row `old_P/new_P` map, excludes source component IDs 9, 18, 19, 35, 52, and 57, and asserts a retained 54 by 54 matrix before rendering panel a. Supplementary producers use a 54-row `cnmf_component/new_P` map. `workflow/program_schema_contract.py` and its fixture-based test require the two schemas to encode the same retained source-to-display mapping. Raw component IDs are exposed only through the generated `program_id_provenance.tsv` translation table.

Fig. 6a cannot be reproduced from this repository. The required derived input `results/xspecies_humanmap_v1/decay_per_program_full.csv` is absent, and no producer is available in the public DAG. The required input schema is one unique `program` identifier plus finite numeric `human_macaque`, `human_mouse`, and `mouse_macaque` columns. Integration also requires an immutable source URI or producer command, byte size, SHA256, row-count validation, and a declared relationship to `program_renumber_map.tsv`. A valid producer must generate `panelA_pairs_long.csv` with columns `program,pair,cosine` and `panelA_medians.csv` with columns `pair,median`. The existing panel renderer is retained, but its presence is not a reproducibility claim.

FigS3 reads the released `programs_master.tsv`, joins it to the retained map, and validates exactly 54 source components before rendering cohort and batch diagnostics. FigS5a uses a released subclass summary that is reduced from 22 by 60 to 22 by 54; FigS5b-e use canonical cell and spatial parquet inputs and the same retained map. FigS6 evaluates 54 retained programs across 14 regions and 22 subclasses, with 182 region pairs and 53 alternative programs in each comparison. FigS8 reads the canonical region tensor, subsets both program axes from 60 to 54, and validates the 1,431 off-diagonal retained-program pairs. The required external and canonical-derived inputs are specified in `DATA_SOURCES.tsv` and `workflow/PIPELINE_DAG.tsv`.

## Release Metadata

This source release is version `2026.07.12`, released on `2026-07-12` under the MIT License. The archival collection is identified by Zenodo concept DOI `10.5281/zenodo.21245200`; a version-specific DOI is recorded only after Zenodo mints it for an immutable archive.

## Verification

```bash
bash SMOKE_TESTS.sh
```

The check validates syntax, mapped source chains, public-boundary rules, stable BLAKE2b seeding, cNMF export naming, and `SHA256SUMS`.
