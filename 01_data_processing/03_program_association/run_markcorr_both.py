#!/usr/bin/env python
"""Chain 07 (cellprog) then 08 (progprog) sequentially; single _DONE on success,
_FAILED (with traceback) on error. Minimal stdout."""
import os, sys, time, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from markcorr_runner_lib import run, OUTDIR

if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    for m in ("_DONE", "_FAILED"):
        p = os.path.join(OUTDIR, m)
        if os.path.exists(p):
            os.remove(p)
    try:
        t0 = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] START 07 cellprog", flush=True)
        run(mode="cellprog", smoke=False, n_perm=1000)
        print(f"[{time.strftime('%H:%M:%S')}] 07 done {(time.time()-t0)/60:.1f}min; START 08 progprog", flush=True)
        run(mode="progprog", smoke=False, n_perm=1000)
        with open(os.path.join(OUTDIR, "_DONE"), "w") as f:
            f.write(f"cellprog+progprog done {time.strftime('%Y-%m-%d %H:%M:%S')} "
                    f"total {(time.time()-t0)/60:.1f}min\n")
        print(f"[{time.strftime('%H:%M:%S')}] ALL DONE {(time.time()-t0)/60:.1f}min", flush=True)
    except Exception:
        tb = traceback.format_exc()
        with open(os.path.join(OUTDIR, "_FAILED"), "w") as f:
            f.write(tb)
        print("FAILED:\n" + tb, flush=True)
        sys.exit(1)
