#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${CORTEX_PROGRAM_ROOT:-/DATA/cortex_nmf_program}"
BAG_URL='https://data.nemoarchive.org/publication_release/Langlieb_Macosko_WMB_Atlas_2023/Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz'
ARCHIVE_NAME='Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz'
BAG_NAME='Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x'
OUTPUT_DIR=''
RESOLVE_FETCH=0
while (($#)); do
  case "$1" in
    --output-dir)
      OUTPUT_DIR=${2:?value required for --output-dir}
      shift 2
      ;;
    --resolve-fetch)
      RESOLVE_FETCH=1
      shift
      ;;
	    *)
	      exit 2
	      ;;
	  esac
done
[[ -n "$OUTPUT_DIR" ]] || { exit 2; }
command -v curl >/dev/null || { exit 127; }
command -v tar >/dev/null || { exit 127; }
command -v bdbag >/dev/null || { exit 127; }
mkdir -p "$OUTPUT_DIR"
ARCHIVE="$OUTPUT_DIR/$ARCHIVE_NAME"
BAG_DIR="$OUTPUT_DIR/$BAG_NAME"
if [[ ! -f "$ARCHIVE" ]]; then
  curl --fail --location --retry 3 --output "$ARCHIVE" "$BAG_URL"
fi
if [[ ! -d "$BAG_DIR" ]]; then
  tar -xzf "$ARCHIVE" -C "$OUTPUT_DIR"
fi
bdbag --validate structure "$BAG_DIR"
if [[ "$RESOLVE_FETCH" -eq 1 ]]; then
  bdbag --resolve-fetch all --validate full "$BAG_DIR"
else
  :
fi
