#!/usr/bin/env python3
"""
make_tableS3_betweenchip.py
生成新 TableS3_between_chip_colocalization.tsv + .xlsx
用新 between-chip Stouffer 結果取代舊 within-chip null
2026-06-25
"""

import pandas as pd
import numpy as np
from pathlib import Path
import datetime

# ── 路徑設定 ──────────────────────────────────────────────────────────────────
BASE = Path("CORTEX_PROGRAM_ROOT")
SRC_DIR = BASE / "results/crossregion_v1/markcorr_betweenchip_v1"
OUT_DIR = BASE / "figure_release/SUBMISSION_final/supplementary"
MANIFEST = BASE / "RUN_LOG_W-tableS3_2026-06-25.md"

cellprog_tsv = SRC_DIR / "betweenchip_cellprog_stouffer_q.tsv"
progprog_tsv = SRC_DIR / "betweenchip_progprog_stouffer_q.tsv"
out_tsv  = OUT_DIR / "TableS3_between_chip_colocalization.tsv"
out_xlsx = OUT_DIR / "TableS3_between_chip_colocalization.xlsx"

def log_manifest(msg: str):
    """Append timestamped entry to RUN_LOG."""
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    with open(MANIFEST, "a") as f:
        f.write(f"\n[{ts}] {msg}")

log_manifest("START make_tableS3_betweenchip.py")

# ── 讀入新 TSV ────────────────────────────────────────────────────────────────
log_manifest("讀 betweenchip_cellprog_stouffer_q.tsv")
cp = pd.read_csv(cellprog_tsv, sep="\t")
log_manifest(f"  cellprog rows={len(cp)}, cols={list(cp.columns)}")

log_manifest("讀 betweenchip_progprog_stouffer_q.tsv")
pp = pd.read_csv(progprog_tsv, sep="\t")
log_manifest(f"  progprog rows={len(pp)}, cols={list(pp.columns)}")

# ── 驗證欄位存在 ─────────────────────────────────────────────────────────────
required_cols = [
    "mode", "A_name", "B_name", "n_chips",
    "Z_combined", "p_stouffer", "Z_combined_unweighted", "p_stouffer_unweighted",
    "median_log2g", "iqr_log2g", "frac_same_sign", "n_same_sign",
    "q_bh", "is_headline", "p_binom", "at_perm_floor", "p_final", "p_final_note"
]
for col in required_cols:
    assert col in cp.columns, f"Missing col in cellprog: {col}"
    assert col in pp.columns, f"Missing col in progprog: {col}"
log_manifest("欄位驗證通過")

# ── 合併 ──────────────────────────────────────────────────────────────────────
df = pd.concat([cp, pp], ignore_index=True)
log_manifest(f"合併後 rows={len(df)}")
assert len(df) == 2619, f"Expected 2619 rows, got {len(df)}"

# ── 欄位重命名（schema §9）───────────────────────────────────────────────────
# pair_type: "cellprog" → "cell-type x program", "progprog" → "program x program"
pair_type_map = {
    "cellprog": "cell-type x program",
    "progprog": "program x program"
}
df["pair_type"] = df["mode"].map(pair_type_map)
df = df.rename(columns={
    "A_name": "A_label",
    "B_name": "B_label",
})

# ── 欄位排序（schema §9 順序）────────────────────────────────────────────────
col_order = [
    "pair_type",
    "A_label",
    "B_label",
    "n_chips",
    "Z_combined",
    "p_stouffer",
    "Z_combined_unweighted",
    "p_stouffer_unweighted",
    "median_log2g",
    "iqr_log2g",
    "frac_same_sign",
    "n_same_sign",
    "q_bh",
    "is_headline",
    "p_binom",
    "at_perm_floor",
    "p_final",
    "p_final_note",
    "mode",          # 保留原始 mode 欄（snake_case）在最後
]
# 確認所有欄都在 df 中
for c in col_order:
    assert c in df.columns, f"Missing col in merged df: {c}"

df = df[col_order]
log_manifest(f"欄位排序完成: {list(df.columns)}")

# ── 驗證 is_headline 數量 ─────────────────────────────────────────────────────
n_headline_cp = df[df["mode"] == "cellprog"]["is_headline"].sum()
n_headline_pp = df[df["mode"] == "progprog"]["is_headline"].sum()
log_manifest(f"is_headline: cellprog={n_headline_cp}, progprog={n_headline_pp}, total={n_headline_cp+n_headline_pp}")
print(f"is_headline: cellprog={n_headline_cp}, progprog={n_headline_pp}, total={n_headline_cp+n_headline_pp}")
# spec 預期: 318 + 1208 = 1526
if n_headline_cp != 318:
    print(f"WARNING: cellprog headline={n_headline_cp}, expected 318")
    log_manifest(f"WARNING: cellprog headline={n_headline_cp}, expected 318")
if n_headline_pp != 1208:
    print(f"WARNING: progprog headline={n_headline_pp}, expected 1208")
    log_manifest(f"WARNING: progprog headline={n_headline_pp}, expected 1208")

# ── 輸出 TSV ─────────────────────────────────────────────────────────────────
log_manifest(f"寫出 TSV → {out_tsv}")
df.to_csv(out_tsv, sep="\t", index=False)
log_manifest(f"TSV 寫出完成，size={out_tsv.stat().st_size} bytes")

# ── 輸出 XLSX（兩 sheet）──────────────────────────────────────────────────────
log_manifest(f"寫出 XLSX → {out_xlsx}")

df_cp = df[df["mode"] == "cellprog"].copy()
df_pp = df[df["mode"] == "progprog"].copy()

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    df_cp.to_excel(writer, sheet_name="cellprog_pairs", index=False)
    df_pp.to_excel(writer, sheet_name="progprog_pairs", index=False)

    # 加一個 combined sheet
    df.to_excel(writer, sheet_name="all_pairs", index=False)

log_manifest(f"XLSX 寫出完成（3 sheets: cellprog_pairs / progprog_pairs / all_pairs）")
print(f"XLSX 寫出完成，size={out_xlsx.stat().st_size} bytes")

# ── 最終驗證 ─────────────────────────────────────────────────────────────────
log_manifest("=== 最終驗證 ===")

# 讀回驗證
df_verify = pd.read_csv(out_tsv, sep="\t")
assert len(df_verify) == 2619, f"TSV row count mismatch: {len(df_verify)}"
required_new_cols = ["Z_combined", "p_stouffer", "p_binom", "p_final", "frac_same_sign", "q_bh", "is_headline"]
for c in required_new_cols:
    assert c in df_verify.columns, f"Missing col in output: {c}"
log_manifest(f"驗證通過：rows={len(df_verify)}, 新欄全在")

print("\n=== 完成 ===")
print(f"TSV: {out_tsv}")
print(f"XLSX: {out_xlsx}")
print(f"Rows: {len(df_verify)}")
print(f"Cols: {list(df_verify.columns)}")
log_manifest("DONE make_tableS3_betweenchip.py")
