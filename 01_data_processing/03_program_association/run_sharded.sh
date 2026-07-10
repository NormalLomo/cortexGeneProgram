#!/usr/bin/env bash
# Sharded parallel production run for markcorr (cellprog + progprog).
# K=16 shards/mode; each shard process handles chips where (index % 16 == shard)
# and writes ONLY per_chip/<mode>_<chip>.npz (sharding => no in-line aggregate).
# Resumable: runner skips any chip whose per_chip dump already exists.
# After both modes' shards complete, aggregate offline (all 4 methods) and mark
# _DONE_shardrun (or _FAILED_shardrun on any error).
set -u

ANALYSIS_DIR=CORTEX_PROGRAM_ROOT/scripts/analysis
OUTDIR=CORTEX_PROGRAM_ROOT/results/crossregion_v1/markcorr
SHARDLOGS=$OUTDIR/_shardlogs
PY=python
K=16

# per-process cupy backstop: single-chip peak ~1.7GB + cublas workspace ~0.5-1GB
# -> 4500MB/proc covers one chip + overhead (incl progprog 60x60 > cellprog 22x60)
# while bounding pool growth; 16 x 4500MB = 72GB < ~80GB free GPU (headroom).
export CUPY_GPU_MEMORY_LIMIT=4500000000   # 4500 MB (bytes)
export PYTHONNOUSERSITE=1

mkdir -p "$SHARDLOGS"
cd "$ANALYSIS_DIR" || { touch "$OUTDIR/_FAILED_shardrun"; exit 1; }

# clear any stale shard-run markers (NOT the per_chip dumps; those drive resume)
rm -f "$OUTDIR/_DONE_shardrun" "$OUTDIR/_FAILED_shardrun"

run_mode () {
    local mode=$1
    echo "[$(date +%H:%M:%S)] launching $K shards for mode=$mode"
    local pids=()
    for i in $(seq 0 $((K-1))); do
        "$PY" "${mode_script}" --shard "$i" --nshards "$K" \
            > "$SHARDLOGS/${mode}_shard${i}.log" 2>&1 &
        pids+=($!)
    done
    local rc=0
    for p in "${pids[@]}"; do
        wait "$p" || rc=1
    done
    echo "[$(date +%H:%M:%S)] mode=$mode shards done rc=$rc"
    return $rc
}

overall_rc=0

mode_script=07_markcorr_cellprog.py
run_mode cellprog || overall_rc=1

mode_script=08_markcorr_progprog.py
run_mode progprog || overall_rc=1

if [ "$overall_rc" -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] a shard process FAILED -> _FAILED_shardrun"
    touch "$OUTDIR/_FAILED_shardrun"
    exit 1
fi

echo "[$(date +%H:%M:%S)] all shards ok; aggregating (all methods) ..."
"$PY" markcorr_aggregate_from_dumps.py --mode cellprog --method all \
    > "$SHARDLOGS/aggregate_cellprog.log" 2>&1 \
  && "$PY" markcorr_aggregate_from_dumps.py --mode progprog --method all \
    > "$SHARDLOGS/aggregate_progprog.log" 2>&1
if [ $? -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] aggregate FAILED -> _FAILED_shardrun"
    touch "$OUTDIR/_FAILED_shardrun"
    exit 1
fi

touch "$OUTDIR/_DONE_shardrun"
echo "[$(date +%H:%M:%S)] ALL DONE -> _DONE_shardrun"
