#!/usr/bin/env bash
set -euo pipefail

: "${CORTEX_PROGRAM_DATA_ROOT:?Source environment/activate_paths.sh before this command}"
root="$(cd "$(dirname "$0")/.." && pwd)"

link_alias() {
  local target="$1"
  local alias="$2"
  if [ -e "$alias" ] && [ ! -L "$alias" ]; then
    printf 'Refusing to replace non-symlink path: %s\n' "$alias" >&2
    exit 1
  fi
  ln -sfn "$target" "$alias"
}

link_alias "$root" "$root/CORTEX_PROGRAM_ROOT"
link_alias "$CORTEX_PROGRAM_DATA_ROOT" "$root/CORTEX_PROGRAM_DATA_ROOT"
printf 'Configured source aliases under %s\n' "$root"
