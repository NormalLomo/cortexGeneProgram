#!/bin/bash
# Run per-chip cross-corr over all 44 row groups serially. Logs to crosscorr_run.log
cd CORTEX_PROGRAM_ROOT
PY=python
SCRIPT=scripts/extended/spatial_crosscorr_perchip.py
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
echo "START $(date)"
$PY - "$@" <<'PYEOF'
import json,subprocess,sys,os
m=json.load(open('CORTEX_PROGRAM_ROOT/results/crossregion_v1/spatial_crosscorr/_chipmap.json'))
# list of (chip, region, rg)
for chip,region,rg in m['order']:
    out=f"CORTEX_PROGRAM_ROOT/results/crossregion_v1/spatial_crosscorr/_perchip/{chip}.npz"
    if os.path.exists(out):
        print(f"skip {chip} (exists)",flush=True); continue
    r=subprocess.run(['python','CORTEX_PROGRAM_ROOT/scripts/extended/spatial_crosscorr_perchip.py',chip,region,str(rg)])
    if r.returncode!=0:
        print(f"FAILED {chip} rc={r.returncode}",flush=True); sys.exit(1)
print("ALL_CHIPS_DONE",flush=True)
PYEOF
echo "END $(date)"
echo "EXIT_MARKER_DONE"
