# Source Receipts

This directory contains metadata receipts only. It does not contain source matrices, spatial count files, chip payloads, or other third-party data.

`chen_cnp0002035_spatial/provider_manifest.json` is a provider metadata receipt for dataset `1663381185152036865` and data code `7cae09e7b7494af988cadcc7999e20c6`. It contains 174 file records, including the three region CSV objects identified in `DATA_SOURCES.tsv`. The file is retained because its bytes and SHA256 are included below; no provider redistribution authorization is included in this release.

`chen_cnp0002035_spatial/chips.md5` is the existing 161-row checksum receipt for macaque spatial text objects (119 macaque1, 19 macaque2, and 23 macaque3 entries). No recheck log, command receipt, or independent hash result is included. Full payload rehashing was not performed; inclusion of this receipt does not claim that every third-party payload was rematerialized or revalidated for this release.

Receipt SHA256 values:

```text
5135fd89419e53cf29a6b945b5808a0b4394ae4017a73e6c23837675e30e1996  chen_cnp0002035_spatial/provider_manifest.json
73817d9aa14fd1c30f6fb1f2d536e77e6859dfd48cd52ca80fbdb17a5dd305e1  chen_cnp0002035_spatial/chips.md5
```

No qualifying Macosko 761,378-nucleus or Han 49-section receipt was found during the read-only source-host search. Their unresolved statuses remain in `DATA_SOURCES.tsv`.
