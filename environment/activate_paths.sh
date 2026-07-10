#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  printf 'Usage: source environment/activate_paths.sh config/paths.env\n' >&2
  return 2 2>/dev/null || exit 2
fi

set -a
. "$1"
set +a
: "${CORTEX_PROGRAM_DATA_ROOT:?CORTEX_PROGRAM_DATA_ROOT is required}"
: "${CORTEX_PROGRAM_RESULTS_ROOT:?CORTEX_PROGRAM_RESULTS_ROOT is required}"
export CORTEX_PROGRAM_DATA_ROOT CORTEX_PROGRAM_RESULTS_ROOT
