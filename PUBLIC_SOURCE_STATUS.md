# Public Source Status

The release contains code and no raw data. `DATA_SOURCES.tsv` is the machine-readable route authority; the source-access READMEs next to the acquisition and validation scripts carry the same boundaries.

Jorstad/CELLxGENE is an executable public API route. The downloader queries the collection, selects exactly the five verified Supercluster objects, and records object IDs, provider URLs, bytes, and ETags at acquisition. The provider does not expose SHA256 for those current objects, so the receipt does not invent one.

Wei/STOmics is a verified 69-file API listing. The provider currently returns `download=false` for all 69 records. This is **manual provider access**, not a missing code chain: `00_record_wei_stomics_manifest.py` records the exact listing and materialized files must receive a local size and SHA256 receipt before preparation. The human aggregate is `STDF000000000022418` (`snRNA.h5ad`); metadata is `STDF000000000022423`. The 44 spatial H5AD records additionally require a verified domain TSV.

The NeMO BDBag route is executable for the Raw_10x archive and its 1,068 provider-manifested FASTQ payloads. It does not establish an exact public object route for the reconstructed H5AD and metadata that the 761,378-nucleus Isocortex branch requires. That derived edge remains unresolved.

Han/CNP0003837 is a public project with provider-controlled section selection. The project-level record is sufficient to direct users to the provider flow, but it is not evidence for an invented section filename or checksum. `03_validate_han_sections.py` verifies only a nonempty user-acquired receipt.

Chen/CNP0002035 has exact, provider-controlled snRNA and spatial route evidence. The release records verified object IDs, filenames, byte counts, and hashes where known; `00_validate_chen_provider_files.py` checks the three exact snRNA RDS files after manual provider access. It does not fetch, bundle, or redistribute provider data.

All provider routes are subject to current provider terms, attribution requirements, and applicable human or nonhuman genomic-data conditions. No third-party data or provider metadata payloads are included in this repository.
