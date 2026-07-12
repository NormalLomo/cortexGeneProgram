#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"

if [[ "${CORTEX_SMOKE_DISPOSABLE:-0}" != "1" ]]; then
  scratch="$(mktemp -d "${TMPDIR:-/tmp}/cortex-gene-program-smoke.XXXXXX")"
  trap 'rm -rf "$scratch"' EXIT
  mkdir -p "$scratch/tree"
  tar -C "$root" \
    --exclude='.git' \
    --exclude='__pycache__' \
    --exclude='.pytest_cache' \
    --exclude='*.pyc' \
    -cf - . | tar -C "$scratch/tree" -xf -
  CORTEX_SMOKE_DISPOSABLE=1 bash "$scratch/tree/SMOKE_TESTS.sh"
  exit 0
fi

cd "$root"
export PYTHONDONTWRITEBYTECODE=1

while IFS= read -r -d '' file; do
  python -c 'import ast, pathlib, sys; ast.parse(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"), filename=sys.argv[1])' "$file"
done < <(find . -type f -name '*.py' -print0)

if command -v Rscript >/dev/null; then
  while IFS= read -r -d '' file; do
    Rscript -e 'parse(file=commandArgs(trailingOnly=TRUE)[1])' "$file" >/dev/null
  done < <(find . -type f -name '*.R' -print0)
fi

while IFS= read -r -d '' file; do
  bash -n "$file"
done < <(find . -type f -name '*.sh' -print0)

python -m unittest discover -s tests
python workflow/release_text_scan.py
python run_release.py --list | test "$(wc -l)" -eq 30
python run_release.py --validate
python 01_data_processing/00_human_snrna_discovery/02_run_cnmf_discovery.py --config config/cnmf_discovery.yaml --dry-run >/dev/null
python - <<'PY'
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

path = Path("01_data_processing/06_null_models/stable_seed.py")
spec = spec_from_file_location("stable_seed", path)
module = module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)
assert module.stable_seed("a", 1) == module.stable_seed("a", 1)
assert module.stable_seed("a", 1) != module.stable_seed("b", 1)
PY

for name in BUILD_REPORT.md LICENSE_DECISION_REQUIRED.txt NUMERIC_20_REVIEW.tsv; do
  test ! -e "$name"
done
test ! -e AGENT_MANIFEST.md
test ! -d .git
test -z "$(find . -type d -name __pycache__ -print -quit)"
test -z "$(find . -type f \( -iname '*.h5ad' -o -iname '*.rds' -o -iname '*.parquet' -o -iname '*.csv' -o -iname '*.xlsx' -o -iname '*.sqlite' -o -iname '*.db' -o -iname '*.zip' -o -iname '*.tar' -o -iname '*.gz' \) -print -quit)"

scan() {
  local pattern="$1"
  local matches
  matches="$(rg -l -i --pcre2 --glob '*.{py,R,r,sh,md,txt,tsv,yaml,yml,json,cff}' "$pattern" . -g '!SMOKE_TESTS.sh' -g '!tests/**' || true)"
  test -z "$matches"
}

scan_sensitive() {
  local pattern="$1"
  local matches
  matches="$(rg -l --pcre2 --glob '*.{py,R,r,sh,md,txt,tsv,yaml,yml,json,cff}' "$pattern" . -g '!SMOKE_TESTS.sh' -g '!tests/**' || true)"
  test -z "$matches"
}

scan '\\bdraft\\b|manuscript|portal|website|\\breview\\b|reviewer|audit|internal'
scan_sensitive 'PROJECT_ROOT|(?<!CORTEX_PROGRAM_)DATA_ROOT|/Users/|/home/|/mnt/|192\.168\.|fsfy|gpuserver|luomeng@|[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}'
scan 'check-only|NUMERIC_20|20[[:space:]]+replicates|n_iter[[:space:]]*[:=][[:space:]]*20'
scan 'worker candidate|release candidate|submitted master|criticfix|\[restore\]|old[[:space:]_-]+version|previous[[:space:]_-]+version'
scan 'owner content revision|deprecated|restored legacy'
scan_sensitive '\bP0\b'

if command -v sha256sum >/dev/null; then
  sha256sum -c SHA256SUMS
else
  shasum -a 256 -c SHA256SUMS
fi
echo 'Public release smoke checks passed.'
