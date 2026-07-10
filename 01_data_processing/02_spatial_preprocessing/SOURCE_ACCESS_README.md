# Spatial Source Access

## Han/CNP0003837

The official project reports 255 samples, 1,168 experiments/runs, and approximately 125.58 TB. Use the provider's data, metadata, or FTP-link flow to select exact section files. No filename, object ID, size, or hash is supplied here until returned by that flow. Put each selected file in a user receipt with `section_id`, `relative_file`, and `sha256`, then run `03_validate_han_sections.py`. This is manual provider access, not missing release code.

## Chen/CNP0002035

The verified provider route identifies dataset `1663381185152036865`, `dataCode=7cae09e7b7494af988cadcc7999e20c6`, and 174 files totaling 605,767,495,166 bytes. Exact snRNA files are `snRNA.metadata.2monkeys.rds` (31,466,747 bytes), `snRNA.sparseMatrix_Monkey1.counts.rds` (4,129,681,166 bytes), and `snRNA.sparseMatrix_Monkey2.counts.rds` (3,544,696,257 bytes); use `00_validate_chen_provider_files.py` after manual provider access. Known spatial records include object `1839599677524664321` / `regions-macaque1.csv`, the corresponding macaque2 and macaque3 region files, `chips.md5`, and contour files. Preserve the provider-returned metadata and MD5 receipts locally; do not redistribute them through this repository.
