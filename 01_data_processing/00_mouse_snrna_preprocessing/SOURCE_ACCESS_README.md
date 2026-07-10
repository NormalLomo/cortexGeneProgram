# Mouse Source Access

`00_fetch_macosko_nemo_bag.sh` retrieves the exact NeMO BDBag archive `Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz` for `nemo:dat-y5zxh0y`. The downloaded archive is a 44,972-byte bag descriptor; `bdbag --resolve-fetch all --validate full` resolves the 1,068 provider-manifested FASTQ payloads and requires approximately 25 TB.

The BDBag is not the reconstructed `Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad` plus `snRNA_cellMeta.csv` required by `00_fix_subset_mouse.py`. The 761,378-cell Isocortex edge therefore remains unresolved rather than being inferred from the raw bag. Do not add FASTQs, derived H5AD files, or metadata to the release.
