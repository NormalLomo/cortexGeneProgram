#!/usr/bin/env bash
set -euo pipefail

# Download the public NeMO Raw_10x BDBag and validate its provider manifests.
# The bag is 25 TB when fully resolved and does not contain the derived H5AD or
# cell metadata required by 00_fix_subset_mouse.py.

BAG_URL='https://data.nemoarchive.org/publication_release/Langlieb_Macosko_WMB_Atlas_2023/Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz'
ARCHIVE_NAME='Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x.tgz'
BAG_NAME='Raw_data_bag_1_Langlieb_Macosko_WMB_Atlas_2023_Raw_10x'
OUTPUT_DIR=''
RESOLVE_FETCH=0

usage() {
  cat <<'EOF'
Usage: 00_fetch_macosko_nemo_bag.sh --output-dir DIR [--resolve-fetch]

Downloads the exact public NeMO Raw_10x BDBag, validates its structure, and
writes a local receipt. --resolve-fetch materializes all 1,068 raw FASTQ
payloads (about 25 TB) and performs a full provider-manifest validation.
EOF
}

while (($#)); do
  case "$1" in
    --output-dir)
      OUTPUT_DIR=${2:?missing value for --output-dir}
      shift 2
      ;;
    --resolve-fetch)
      RESOLVE_FETCH=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -n "$OUTPUT_DIR" ]] || { usage >&2; exit 2; }
command -v curl >/dev/null || { echo 'curl is required' >&2; exit 127; }
command -v tar >/dev/null || { echo 'tar is required' >&2; exit 127; }
command -v bdbag >/dev/null || { echo 'Install bdbag before running this script' >&2; exit 127; }

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
{
  printf 'collection\tnemo:dat-y5zxh0y\n'
  printf 'bag_url\t%s\n' "$BAG_URL"
  printf 'archive_sha256\t%s\n' "$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')"
  printf 'payload_contract\t1068 raw FASTQ payloads; provider MD5 manifest\n'
  printf 'resolution\t%s\n' "$([[ "$RESOLVE_FETCH" -eq 1 ]] && echo full || echo structure_only)"
} > "$OUTPUT_DIR/macosko_nemo_bdbag_receipt.tsv"

if [[ "$RESOLVE_FETCH" -eq 1 ]]; then
  bdbag --resolve-fetch all --validate full "$BAG_DIR"
  echo "Resolved and validated $BAG_DIR"
else
  echo "Validated BDBag structure. Re-run with --resolve-fetch to materialize the 25 TB payload."
fi
