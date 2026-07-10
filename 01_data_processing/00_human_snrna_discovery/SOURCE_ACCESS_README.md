# Human Source Access

## Jorstad/CELLxGENE

`00_download_jorstad_cellxgene.R` queries collection `d17249d2-0e6e-4500-abb8-e6c93fa1ac6f` and accepts only five Supercluster datasets: `5346f9c6-755e-4336-94cc-38706ec00c2f`, `9c63201d-bfd9-41a8-bbbc-18d947556f3d`, `c3aa4f95-7a18-4a7d-8dd8-ca324d714363`, `d01c9dff-abd1-4825-bf30-2eb2ba74597e`, and `e4ddac12-f48f-4455-8e8d-c2a48a683437`. It records dynamic provider object URLs, bytes, and ETags in `SOURCE_CHECKSUMS.tsv`; the provider does not publish SHA256 for those objects. Do not treat an older same-title copy as the current object and do not redistribute H5AD files.

## Wei/STOmics

`00_record_wei_stomics_manifest.py` records the exact `STDS0000242` 69-row listing: 44 spatial H5AD files, 23 snRNA H5AD files, `STDF000000000022418` / `snRNA.h5ad`, `STDF000000000022422` / `SnRNA_seurat.RDS`, and `STDF000000000022423` / `SnRNA_Meta_with_EdLein_dataset_included.csv`. The provider returns `download=false`; that is manual provider access, not a missing acquisition implementation. After provider-controlled retrieval, record local bytes and SHA256 before invoking `01_prepare_jorstad_wei.py`. No human payload belongs in this release.
