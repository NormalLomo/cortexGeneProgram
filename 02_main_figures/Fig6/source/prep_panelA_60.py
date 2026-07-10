#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
WorkerB P0-4 fix: Fig 6 panel a 切回 60-母體 cosine median
（與正文/Fig 1b/Methods 一致：0.499 / 0.294 / 0.560）

Fig 6 其他 panel (B/E/F/G) 仍用 54 母體不動。
panel a 是「per-species-pair distribution」, 60 vs 54 數值差很小
（最大差在 mouse-macaque: 0.5431 vs 0.5596），但精確匹配正文 = 採 60。

讀 decay_per_program_full.csv（未 DROP6, 60 program）→ 重生 panelA_pairs_long.csv + panelA_medians.csv
"""
import os, pandas as pd

BASE = "CORTEX_PROGRAM_ROOT/results/xspecies_humanmap_v1"
OUT = os.path.join(BASE, "figures", "Fig10_v2", "data")

decay = pd.read_csv(os.path.join(BASE, "decay_per_program_full.csv"))
# 不剔 DROP6, 60 母體

a = decay.melt(id_vars="program",
               value_vars=["human_macaque", "human_mouse", "mouse_macaque"],
               var_name="pair", value_name="cosine")
pair_label = {"human_macaque": "human-macaque",
              "human_mouse":   "human-mouse",
              "mouse_macaque": "mouse-macaque"}
a["pair"] = a["pair"].map(pair_label)
a.to_csv(os.path.join(OUT, "panelA_pairs_long.csv"), index=False)
med = a.groupby("pair")["cosine"].median().reset_index().rename(columns={"cosine": "median"})
med.to_csv(os.path.join(OUT, "panelA_medians.csv"), index=False)

print("WorkerB P0-4 — Fig 6 panel a 已重生為 60 母體：")
print(med)
print(f"n per pair = {a.groupby('pair').size().iloc[0]} (應為 60)")
