#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Update supplementary-figure references and the S2 ANCOVA entry in public article text."""
import re
TEXT_ROOT="CORTEX_PROGRAM_ROOT/publication_text"

# ============ EN ============
fn=f"{TEXT_ROOT}/article_en.md"; t=open(fn,encoding="utf-8").read(); n=0

# 1) panel-ref: match with flexible hyphen class
pat = re.compile(r"reported jointly in the multi[-‐‑‒–]criterion K[-‐‑‒–]selection panel \(Supplementary Fig\. K[-‐‑‒–]selection multi[-‐‑‒–]criterion\)")
t2,c = pat.subn("reported jointly in the multi-criterion K-selection panel (Fig. S3)", t);
if c: t=t2; n+=c
print("EN panel-ref replaced:",c)

# 2) supp-material boxed sentence L303: remove native, use contiguous numbering wording
old_sm = ("supplementary figures (K-robustness, spatial cross-correlation, identity-backbone "
          "cluster-confusion, cohort QC, cNMF K-selection diagnostic; final numbering by figure assembly)")
new_sm = ("supplementary figures (Fig. S1 K-robustness; Fig. S2 cohort/donor-ANCOVA classification; "
          "Fig. S3 multi-criterion K-selection; Fig. S4 spatial cross-correlation; Fig. S5 "
          "identity-backbone cluster-confusion; Fig. S6 program–disease recall with N-sensitivity)")
if old_sm in t: t=t.replace(old_sm,new_sm); n+=1; print("EN supp-material sentence updated")
else: print("WARN EN supp-material sentence not found")

# 3) insert S2 (ANCOVA) bullet into in-text list between S1 and S3 bullets
s2_bullet_en = ("\n- **Fig. S2: Cohort/donor-ANCOVA classification of the region-variable set.** "
 "The 14 region-variable programs (per-cell ANOVA, R3) are re-assessed by ANCOVA controlling "
 "for source cohort and donor: seven remain cohort-robust (P1, P3, P4, P6, P8, P10, P14) and "
 "seven are cohort-sensitive (P9, P18, P19, P35, P37, P52, P57); the cohort-robust set is the "
 "one highlighted throughout Fig. 3a–c,g.")
# anchor: end of S1 bullet (right before the S3 bullet line)
anchor = "\n- **Fig. S3: K-selection (multi-criterion, K = 30–200).**"
if anchor in t:
    t=t.replace(anchor, s2_bullet_en + anchor, 1); n+=1; print("EN S2 bullet inserted")
else:
    print("WARN EN S3 anchor not found for S2 insert")

open(fn,"w",encoding="utf-8").write(t); print("EN edits:",n)

# ============ ZH ============
fn=f"{TEXT_ROOT}/article_zh.md"; t=open(fn,encoding="utf-8").read(); n=0

# supp-material ZH sentence L328
old_zh=("补充图（K 稳健性、空间互相关、身份骨架聚类混淆、队列质控、cNMF K 选择诊断；最终编号由图表组装环节指定）")
new_zh=("补充图（图 S1 K 稳健性；图 S2 队列/供体 ANCOVA 加固；图 S3 多准则 K 选择；图 S4 空间互相关；"
        "图 S5 身份骨架聚类混淆；图 S6 程序–疾病召回与 N 敏感性）")
if old_zh in t: t=t.replace(old_zh,new_zh); n+=1; print("ZH supp-material sentence updated")
else: print("WARN ZH supp-material sentence not found")

# insert ZH S2 (ANCOVA) bullet between 图 S1 and 图 S3 bullets
s2_bullet_zh=("\n- **图 S2 ｜ 区域可变程序集合的队列与供体 ANCOVA 加固。** 对来自逐细胞 ANOVA 初筛（R3）的 "
 "14 个区域可变程序，在控制来源队列与供体（ANCOVA）后重新评估：7 个保持队列稳健"
 "（P1、P3、P4、P6、P8、P10、P14），7 个为队列敏感（P9、P18、P19、P35、P37、P52、P57）；"
 "队列稳健集合即全图 3a–c,g 中高亮的集合。")
anchor_zh="\n- **图 S3 ｜ K 选择（多准则，K = 30–200）。**"
if anchor_zh in t:
    t=t.replace(anchor_zh, s2_bullet_zh+anchor_zh, 1); n+=1; print("ZH S2 bullet inserted")
else:
    print("WARN ZH S3 anchor not found for S2 insert")

open(fn,"w",encoding="utf-8").write(t); print("ZH edits:",n)
