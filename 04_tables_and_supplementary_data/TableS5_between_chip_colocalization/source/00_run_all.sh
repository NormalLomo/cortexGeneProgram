#!/usr/bin/env bash
# 00_run_all.sh — Full between-chip Stouffer null pipeline orchestrator
# Usage: bash scripts/00_run_all.sh [--smoke]
# Produces outputs in: results/crossregion_v1/markcorr_betweenchip_v1/
#
# Steps:
#   01 cellprog torus-shift (44 chips × 1188 pairs × 1000 perms) ~3h GPU
#   01 progprog torus-shift (44 chips × 1431 pairs × 1000 perms) ~2.5h GPU
#   02 Stouffer combine (both modes)
#   03 BH-FDR (both modes)
#   04 Binomial floor (both modes)
#   05 Xspecies aggregate (placeholder if no xspecies data)

set -euo pipefail

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJROOT="CORTEX_PROGRAM_ROOT"
PYTHON="python"
LOGDIR="${PROJROOT}/analysis/null_betweenchip/logs"
mkdir -p "${LOGDIR}"

SMOKE=""
if [[ "${1:-}" == "--smoke" ]]; then
    SMOKE="--smoke"
    echo "[$(date +%H:%M:%S)] SMOKE MODE — 1 chip, n_perm=10"
fi

log() { echo "[$(date +%H:%M:%S)] $*"; }

# Determinism
export CUBLAS_WORKSPACE_CONFIG=":4096:8"

log "=== Phase 1: cellprog torus-shift null ==="
${PYTHON} "${SCRIPTDIR}/01_torus_shift_perchip.py" --mode cellprog ${SMOKE} \
    2>&1 | tee "${LOGDIR}/01_cellprog.log"
log "cellprog per-chip Z done"

log "=== Phase 1: progprog torus-shift null ==="
${PYTHON} "${SCRIPTDIR}/01_torus_shift_perchip.py" --mode progprog ${SMOKE} \
    2>&1 | tee "${LOGDIR}/01_progprog.log"
log "progprog per-chip Z done"

log "=== Phase 2: Stouffer combine (cellprog) ==="
${PYTHON} "${SCRIPTDIR}/02_stouffer_combine.py" --mode cellprog ${SMOKE} \
    2>&1 | tee "${LOGDIR}/02_cellprog.log"

log "=== Phase 2: Stouffer combine (progprog) ==="
${PYTHON} "${SCRIPTDIR}/02_stouffer_combine.py" --mode progprog ${SMOKE} \
    2>&1 | tee "${LOGDIR}/02_progprog.log"

log "=== Phase 3: BH-FDR ==="
${PYTHON} "${SCRIPTDIR}/03_bh_fdr.py" ${SMOKE} \
    2>&1 | tee "${LOGDIR}/03_bh_fdr.log"

log "=== Phase 4: Binomial floor ==="
${PYTHON} "${SCRIPTDIR}/04_binomial_floor.py" ${SMOKE} \
    2>&1 | tee "${LOGDIR}/04_binomial.log"

log "=== Phase 5: Xspecies aggregate ==="
${PYTHON} "${SCRIPTDIR}/05_aggregate_samesign_xspecies.py" ${SMOKE} \
    2>&1 | tee "${LOGDIR}/05_xspecies.log"

log "=== ALL DONE ==="
echo "Results in: ${PROJROOT}/results/crossregion_v1/markcorr_betweenchip_v1/"
ls -lh "${PROJROOT}/results/crossregion_v1/markcorr_betweenchip_v1/"
