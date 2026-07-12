# Public Source Status

The release contains code and no raw data. `DATA_SOURCES.tsv` is the machine-readable route authority; the source-access READMEs next to the acquisition and validation scripts carry the same boundaries.

Jorstad/CELLxGENE is an executable public API route. The downloader queries the collection, selects exactly the five verified Supercluster objects, and records object IDs, provider URLs, bytes, and ETags at acquisition. The provider does not expose SHA256 for those current objects, so the receipt does not invent one.

Wei/STOmics is a verified 69-file API listing. The provider currently returns `download=false` for all 69 records. This is **manual provider access**, not a missing code chain: `00_record_wei_stomics_manifest.py` records the exact listing and materialized files must receive a local size and SHA256 receipt before preparation. The human aggregate is `STDF000000000022418` (`snRNA.h5ad`); metadata is `STDF000000000022423`. The 44 spatial H5AD records additionally require a verified domain TSV.

The NeMO BDBag route is executable for the Raw_10x archive and its 1,068 provider-manifested FASTQ payloads. It does not establish an exact public object route for the reconstructed H5AD and metadata that the 761,378-nucleus Isocortex branch requires. That derived edge remains unresolved.

Han/CNP0003837 is a public project with provider-controlled section selection. The project-level record is sufficient to direct users to the provider flow, but it is not evidence for an invented section filename or checksum. `03_validate_han_sections.py` verifies only a nonempty user-acquired receipt.

Chen/CNP0002035 has exact, provider-controlled snRNA and spatial route evidence. The release records verified object IDs, filenames, byte counts, and hashes where known; `00_validate_chen_provider_files.py` checks the three exact snRNA RDS files after manual provider access. `source_receipts/chen_cnp0002035_spatial/` contains a provider metadata receipt (`provider_manifest.json`) and an existing MD5 receipt (`chips.md5`), not spatial payloads. No recheck log, host-independent hash result, or timestamped command receipt is included, so this release makes no independent recheck claim. The redistribution basis for the provider metadata receipt is not evidenced in this release.

Fig. 6a requires `results/xspecies_humanmap_v1/decay_per_program_full.csv`. Neither that object nor a producer, immutable source location, byte size, SHA256, validated schema, or retained-program relationship is available in this source tree. The panel renderer therefore remains blocked and is not evidence of reproducibility.

FigS3, FigS5, FigS6, and FigS8 have explicit retained-program source contracts. FigS3 requires the released `programs_master.tsv`; FigS5 requires the released subclass summary plus canonical cell and spatial parquet inputs; FigS6 requires canonical cross-region derived matrices; and FigS8 requires the canonical region tensor and chip map. These scripts validate the 54-program map before analysis. Their external and derived-input requirements are recorded in `DATA_SOURCES.tsv`, and their executable ordering is recorded in `workflow/PIPELINE_DAG.tsv`.

All provider routes are subject to current provider terms, attribution requirements, and applicable human or nonhuman genomic-data conditions. No third-party scientific data payloads are included in this repository. The only provider metadata receipt retained here is the explicitly identified Chen object listing; its redistribution basis remains undocumented.
