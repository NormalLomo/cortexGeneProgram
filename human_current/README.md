# 人类皮层 NMF：现稿分析代码

本目录整理原始计数入口、SCT、cNMF、程序分数、区域模型、空间关系、外部数据和图表代码。它不是一次新的分析运行；本次没有执行这些生产脚本或任何测试，也没有复制表达矩阵、模型或大结果。原始目录未改动。

## 使用边界

- 当前科学对象为人类 K60 参考中的保留 54 程序；正式区域 21 程序是 LMM、limma duplicateCorrelation、dream **三方法**分别在 54 程序内 BH 后的交集。“六类输出”指三方法×总体/亚类两层，不是六方法共识。原组件 9、18、19、35、52、57 排除，原组件 ID 与 P1–P54 不可直接混用。
- 公开交付的 cNMF 主路线为 `public_n100`，在一次任务中对每个候选 K 分别进行 **100 次独立初始化**。100 是初始化次数，不是一次优化的迭代步数。论文既有共识来自分两批完成的 80+20 次运行，原 80 次结果保存在另一台机器，本次不迁移，也不作为公开入口的运行依赖。`historical_n20` 保留历史 20 次支线供来源追溯，不是公开主入口。现存 K60 权重/usage 仍是现稿下游分析的固定输入。
- 最早可追踪的实际发现入口是已合并 `3_SnRNA_seurat_merged_1m_Cells.RDS`。源论文原始文件到该 RDS 的 QC/合并 producer 尚未定位。单次任务 100 次初始化的代码路线与历史 80+20 的原运行记录分别说明，本次没有执行新的发现任务，也不声称已从公共原始文件完整复现论文。
- `07_figures_and_tables/fig3_layer_region` 已同步当前 Fig3 的局部统计源与无 b 黑框绘图入口；两页 FigS10 的绘图分支仍保留。固定九 panel PDF 是输入资产，**不是当前主图仍有九个 panel**。本次仅整理代码，未重新拟合统计或生产图件。
- 不包含跨物种、129 程序、NetRep、全文跨数据集注释重算、benchmark 模型训练。外部打分使用随包提供的 fixed-H NMF scorer 与 reference-scoring implementation，不依赖未随包提供的外部运行时。

## 路径与依赖

`shared/paths.example.env` 只列必要环境变量；没有自动总入口。正式执行应使用另行批准的工作位置，不能将可写根设为原始来源仓库。

主要来源包括历史 SCT producer、采用的 cNMF 与 formal21 结果树、旧发布代码、空间 source snapshots 及 standalone 工作树，均作为只读输入。SEA-AD reader源于既有 human single-cell benchmark adapter。旧文件名是来源谱系标识，不是当前图号；实际根目录通过下列环境变量提供。

|变量|接口|
|---|---|
|`NMF_CODE_ROOT`|此代码目录在执行机器上的位置，仅用于找到代码|
|`CORTEX_PROGRAM_ROOT`|可写的 human54 工作树，保留 `inputs/`、`results/`、`figures/`、`scripts/fig2/` 的既有结果布局|
|`NMF_WORK_ROOT`|可写的 standalone 工作树，保留 `analysis/`、`tables/`、`inputs/`、`figures/` 的既有布局|
|`NMF_SOURCE_ROOT`|只读 standalone 来源根；外部脚本从其 `inputs/cortex_nmf_program` 读参考、外部输入和标签映射|
|`NMF_ARCHIVE_SOURCE_ROOT`|Fig3 读取既有空间分数、元数据和 `scripts/fig2/` 表的只读 archived 根，与可写的 `CORTEX_PROGRAM_ROOT` 分开|
|`SNRNA_MERGED_RDS`|上述合并 RNA-count RDS，保留原始 metadata 与 library_prep|
|`SNRNA_WORK_DIR`|SCT 对象、donor MTX、obs/var CSV、snRNA_1M.h5ad 的独立输出位置|
|`SEURAT_OBJECT_SOURCE`、`SEURAT_SOURCE`|历史 producer 使用的本地 SeuratObject/Seurat 源目录；不自动安装或替换版本|
|`SPATIAL_SCT_DIR`|已有每片 `<chip>_sct.h5ad`，`.X` 为校正非负 counts|
|`CNMF_COUNTS_H5AD`|已过滤的 SCT 校正 counts，供历史 n20 入口显式读取|
|`SCSLAT_SITE_PACKAGES`|历史 SLAT dual-PCA 所需安装目录；仅描述性嵌入使用|
|`FORMAL21_R_LIBS`|可选的额外 R library 路径列表，使用系统路径分隔符连接；formal21 默认仍读取项目内 `Rlib` 与常规 `.libPaths()`|

工作树中被读取的现成表、固定 K60 谱与 usage、逐核/空间分数、markcorr tensor、元数据和 native panel 资产需按各脚本接口提供；本包没有复制这些大型数据和全部图形输入。同目录的 `Fig3_selected_nine_panel.pdf` 是已随代码保留的固定资产。环境变量不负责补齐输入，也不会自动重排文件。当前图表脚本的 `NMF_WORK_ROOT/inputs` 和 `tables` 仍是既有目录布局；不要将整个只读源目录误当可写工作根。TSV 与工作簿也不能仅凭同号表名互换。

Python 依赖按模块使用：numpy、pandas、scipy、anndata、scanpy、cnmf、pyarrow、statsmodels、patsy、scikit-learn、PyYAML、openpyxl；绘图与 PDF 处理用 matplotlib、seaborn、networkx、PyMuPDF、pypdf、Pillow、svgutils/相关 SVG 工具；外部 fixed-H 用 torch/CUDA，SEAAD 用 rpy2；历史空间内核可使用 cupy；认知分支用 nilearn/NiMARE；dual-PCA 用 scSLAT。R 依赖包括 Matrix、sctransform、Seurat/SeuratObject、BPCells、qs、flock、tidyverse、lme4、lmerTest、pbkrtest、limma、variancePartition、ComplexHeatmap、circlize、svglite、igraph/ggraph/tidygraph。Fig3 的局部区域对比实际直接调用 lme4 与 pbkrtest，不以 lmerTest 包装入口或其它自由度近似替代。

PDF 预览还依赖 Poppler；Fig4 当前入口调用 `pdftoppm`，候选布局分支另用 `pdftocairo`。绘图沿用系统字体，包含 DejaVu Sans 与 Liberation Sans/Narrow；Fig5 需要可用的 Liberation Sans Narrow 常规/粗体字体。Fig2 原生裁剪使用 PyMuPDF/MuPDF 低层接口与 pypdf 内容流接口，不是仅安装任意 PDF 阅读库即可替代的通用转换。`environment/` 保存分lane环境说明；没有一个虚构的统一锁文件覆盖全部脚本。

所有 R 入口保留明确的 `future.globals.maxSize = 50 * 1024^3`。这不启用并行；原科学计算/线程设置未调优。后续正式运行仍需单独批准输入、脚本及输出。

## 顺序与输入输出

下列顺序用于区分已采用结果的依赖关系，不是一个需要全部重跑的总命令。公开发现入口使用 `public_n100`。复用现稿则从既有 K60 谱/usage、保留 54 映射及各章结果进入相应下游入口，历史 n20 单列保留。新任务得到的组件不能自动继承历史原组件及 P1–P54 编号。三种活动量分别保留：发现集单核 cNMF usage、空间 SCT 校正 counts 对固定谱的投影、外部数据 fixed-H 拟合 usage。

### 1. 合并计数 → SCT → 单核矩阵

`01_input_and_sct/snrna/01_prepareData_snRNA.R` → `02_prep_full_1M.R` → `03_build_h5ad.py`。

第一步按 `library_prep` 拆分，`sctransform::vst` 后以各 library UMI 中位数的中位数作为共同深度调用 `correct(do_round=TRUE, do_pos=TRUE)`。**RNA counts 层是校正计数，不是原始 UMI，不是 residual 或 log1p。** 第二步读取该 counts 层，按 donor 输出 genes×cells MTX 及条码；第三步按 donor 文件顺序合并、转置成 cells×genes，按条码排列 obs，输出 `snRNA_1M.h5ad`。残差另存，不能作为 NMF 的非负输入替代。

### 2. SCT 矩阵 → 发现/固定 K60 → 54

`02_cnmf_and_retained54/historical_n20/01_build_counts.py --counts-only` 完整写出过滤后矩阵后结束；不是试算。实际沿原代码完成 ENSG→symbol、去未映射/重复、V3 blacklist、min_cells=20、min_counts=500。输出为 `results/cnmf_snrna_joint_full1M_v1/snrna_joint_full1M_v1_counts.h5ad`。

历史脚本不加该选项还会进入原 20 次发现分支，**不要把这一分支接成论文 100 次历史**。`02_run_cnmf_main.py` 也是历史 n20，保留当时 K、seed、HVG、density 设置；仅改路径/去除删除旧结果的操作。

补入的 `historical_n20/03_run_cnmf_k65.py` 读取 `CNMF_COUNTS_H5AD`，保留 K=65、20 次初始化的独立分支，工作目录为 `CORTEX_PROGRAM_ROOT/results/cnmf_snrna_joint_full1M_v1/cnmf_work_k65`，导出原命名的载荷与 cell_scores TSV。已有工作目录不会被递归删除，不再写无用的 done 标记。`consensus_k80_90.py` 则读取同一结果根下既有 `cnmf_work/snrna_joint_full1M_v1` 的分解结果，仅补 K=80、K=90 共识及对应载荷/usage 导出，不重新 factorize。**这里的 80 和 90 是秩，不是初始化次数；这两个脚本都没有补回原 80 次批次或 80+20 融合入口。**

**公开主入口**为 `public_n100/02_run_cnmf_discovery.py --config <config> --data-root <input-root> --results-root <output-root>`，使用同目录 `config/cnmf_discovery.yaml`，并由相邻 `03_export_cnmf_outputs.py` 导出结果。入口将 `n_iter=100` 传给 cNMF prepare，继而完成 factorize、combine 和 consensus，每个候选 K 在同一次任务中分别进行100次独立初始化。现有配置保留 K=30、40、50、60、70、80、90，选择 K60，seed=42、HVG=3000、density threshold=0.15，本次未修改这些设置。`counts_h5ad` 必须指向上述过滤后的 SCT 校正矩阵。该入口不需要另一台机器的历史80次文件，不等同于宣称新任务已运行或与原80+20共识逐组件相同。

复用论文结果时，输入已采用的 K60 gene_spectra_score、gene_spectra_tpm、consensus usages 与 `program_renumber_map.tsv`。`build_canonical_map.py` 的 table-s2/s3/s4 参数采用旧发布表 schema；其中 TableS4 必须是覆盖原组件 1–60、含 `cnmf_component`、`new_P`、`status`、`name_short` 及排除依据字段的 **60 行 TSV**，不是现稿仅列保留 54 程序的 Table S4 工作簿。该脚本也未补出原先未明确的技术排除数值阈值。不要互换同号表，或用重新发现的组件次序硬套此映射。

### 3. 单核活动 → 区域、亚类与 formal21

`03_regional_models/reference_profiles/01_region_program_matrix.py` 从 cell_scores+obs 形成 region/subclass 均值和逐核 joined parquet。`02_variability_classify.py`、`03_within_subclass_anova.py`、`annotate_programs.py`、`make_program_names.py` 是原有描述性/注释链，不能取代 donor 推断。

`formal21/input/aggregate_donor_region_scores.py` 与 `aggregate_donor_region_subclass_scores.py` 均接受 `--obs --scores --mapping --outdir`；分别输出 donor×region、donor×region×subclass 均值及对应长表。保留 nucleus counts 与原分组字段，不把核作为独立 donor。

正式总体入口为 `celltype_regional/overall_nature_lmm/run_overall_nature_lmm.R`、`limma_dupcor_validation/run_limma_dupcor_validation.R`、`dream_validation/run_dream_overall.R`。对应亚类入口位于 `celltype_regional/nature_lmm` 及 `regional_donor_validation_v2/methods`。目录名中的 validation 是历史命名，内含论文实际采用的统计模型，**不是本次运行的技术测试**。模型保留 region+age+sex 与 donor 随机/相关结构；亚类有效性条件和 multiple-testing family 不变。

三方法在总体/亚类两层先写各自原有结果文件；`scripts/analysis/integrate_regional_donor_models.R <root>` 直接读取这些原输出，按下游接口导出不改变内容的六表 `results/crossregion_v2/{overall,subclass}_{nature_lmm,limma_dupcor,dream}.tsv`，再写 `program_method_matrix.tsv` 和 `regional_consensus_21.tsv`；其交集不使用目标8校正。`release_candidate/scripts/build/build_table_s4.py` 读取 `FORMAL21_PROJECT_ROOT` 生成现有表格格式。亚类limma使用 `FORMAL21_MODEL_ROOT`（或 `LIMMA_DUCPOR_SOURCE_ROOT`）；默认写该工作根下原methods子目录，不再向代码目录写研究结果。

`HUMAN_VALIDATION_run_analysis.py --root --outdir --chip-donor-lookup` 保留原单核/空间附加分析；空间聚合是 chip×domain **均值**和原加性 GEE，不是新的 region×layer 交互。`01_regional_donor_analysis.py --input --outdir` 保留原 donor robustness：Rademacher 全信息 donor 符号枚举、permutation/Webb 各99,999；现稿使用完整 P1–P54 家族的 `q_all54`。

### 4. 空间 counts → SCT 分数 → 层与局部关系

`01_input_and_sct/spatial/04_sct.R <chip> <counts.mtx> <genes> <bcs> 1352 <corrected.mtx>` 读取 genes×bins raw counts，写同轴 corrected MTX 与基因/条码及另存 residual；正式程序评分只取 corrected counts。

历史 scale 来源为原空间流程的 `04_sct_scale_factor.py`：全片列表各片 `bin_total_umi` 中位数，再取中位数并 round。当前采用值固定 **1352**；不以新目录 glob 或新样本重算来替代，未复制或执行该prepass。

`01_input_and_sct/spatial/04_sct_io.py` 已纳入：`export --chip --rich --tmpdir` 优先读取raw `layers['counts']`、无此层才用 `.X`，按原代码写genes×bins MTX；`build --chip --rich --corrected --scale_factor 1352 --output` 按条码和基因轴回装校正矩阵，未建模基因补0，保留raw层/metadata/RCTD。新副本要求独立输出，拒绝覆盖rich源。将输出指定为评分输入目录的 `<chip>_sct.h5ad`；不把残差层用于程序分数。

`04_spatial_scoring_and_relations/04_spatial_score_sct.py <chip>`：SCT `.X` ×共享基因 TPM 谱，非负 dot product，无额外 CPM/log1p/每-bin归一化；`04_spatial_score_sct_aggregate.py` 按所有片的所有 bin 逐 program 标准化。原 60 列按固定映射取54。**不使用**旧 Pearson-residual 04b/04c 路线。

`07_markcorr_cellprog.py`、`08_markcorr_progprog.py` 调用同目录 `markcorr_runner_lib.py`、`markcorr_core.py`、`markcorr_aggregate.py`，保留 mask、非负截断 marks、10 distance rings、方向与物理坐标换算；`markcorr_median_iqr.py`、`markcorr_byarea_median_iqr.py` 汇总既有 per-chip tensor。`01_torus_shift_perchip.py --mode ...` → `02_stouffer_combine.py --mode ...` → `03_bh_fdr.py` 保留原跨片 null/inference。`02_true_donor_spatial_support.py` 从既有 per-chip tensors 和真实 donor lookup 生成 donor 与 whole-donor LOO，未新增 P/FDR。

50 DNB pixels 对应25µm；[0,25)是同一 measurement bin。主统计 mask 为成功RCTD、非arachnoid且UMI≥100；某些显示用UMI≥200。不能将显示 k25 中位数平滑应用到所有统计。

### 5. 外部数据 → 固定参考活动

`05_external_scores/05_external_single_cell_score_correlation.py` 调用同目录、已随包提供的 `gep_score.py` 中的 `score_counts_fixed_h`；该内核再使用同目录的 `contracts.py`。因此本入口使用随包固定-H NMF scorer/reference scoring implementation，不需要未随包提供的外部运行时。运行仍需 NumPy、SciPy、PyTorch/CUDA 及本节和“路径与依赖”中列明的其它依赖，并须提供既有固定参考、profile、表达数据和标签映射输入。随包代码采用本目录 `LICENSE` 所列 MIT 许可证；第三方原始数据仍遵循来源条款。原始输入的具体 `.X`/raw选择、疾病标签、library RDS及crosswalk路径仍由该脚本按研究定义，**不得一律转成发现集SCT**。

共享基因、重复符号求和，study级ddof=0 SD，零SD贡献置零；固定60谱，非负乘法更新，max_iter1000/tol1e-4、每10步判断改善；先60归一化再取54。nearest-reference再L2并搜索全参考，不按真实标签限制。SEAAD source adapter只读RDS，不进行scANVI训练。

原入口无动作选项时执行其已存在的研究打分顺序；后处理使用 `--study-celltype-program-source`、`--broad-program-consistency-source`、`--cell-nearest-reference` 等已存在选项。输出在 `NMF_WORK_ROOT/analysis/external_single_cell_program_scores`，包括逐细胞parquet及profile/correlation TSV；现稿热图/散点由07中的对应R入口读取。

### 6. 疾病、衰老、认知

`06_disease_aging_cognition/fig8_analysis.py` 生成疾病富集；`fig8_program_disease.py`、`ed_program_disease_supp.py` 消费富集并绘图；`program_aging_ora.py` 对原aging GMT/metadata进行ORA；`fig7_program_cognition.py` 是原认知补充链。保留原gene universe、top-N、BH family和panel显示选择，不能把旧文件名Fig8当作现稿图号。Neurosynth/NiMARE及相关基因集是外部参考依赖，本次未下载。

### 7. 图表

#### Fig3：既有局部统计 → 当前主图

当前绘图入口是 `07_figures_and_tables/fig3_layer_region/fig3_layer_region_revision.py --local-significance`。它读取已保留的 `NMF_WORK_ROOT/analysis/fig3_layer_region/layer_region_profiles.tsv` 与 `fig3_local_donor_contrasts.tsv`，结合下述固定输入，只写该分析目录下的 `Fig3.pdf`、`Fig3.png`，随后返回；**不调用 R、不重拟合、不改两张原始汇总表，也不回写 FigS10**。固定 native 输入仍是代码同目录的 `Fig3_selected_nine_panel.pdf`，不能以当前主图/补图替换。

- Fig3 现为 a–d：a 保留原 program/层域定位；b 是七域 L1–L6+WM×54 程序热图，按核参考亚类偏好分组，短功能名、P 编号和逐列 Z 定标保持。**b 不再画任何黑框，也不再采用列最大值/argmax 标记**；已算出的 b 统计留在配套 TSV，不因去框而删除。星号只表示较低功能注释置信，参考亚类偏好不等于空间细胞定位。
- c完整显示14区域×L1–L6的**84行×54列**，每个program列在全部84显示行上计算z（ddof=1），不分区域单独标准化。随后对行、列分别使用SciPy `linkage/leaves_list/dendrogram`：Euclidean距离、complete linkage、`optimal_ordering=False`，无权重。聚类输入为**未clip的z值**；84行可全打散，区域与层两条annotation仍保留。缺失值不填补，当前源要求用于聚类的84×54矩阵有限。
- c 的双向聚类热图、逐列 Z、布局及色限 −3..3 均保留；色限只是显示饱和，不改变聚类输入。绘图代码读取 c 的局部 BY 结果；**本次已保留结果中没有 c 对比达到 BY q<0.05，因此当前 c 无显著性框**，不是取消了 c 或把其统计删去。b 的对称色限仍沿用七域列 Z 值的原定标规则。d 仍为 P16 的三个真实切片及两个区域 profile，原分数、未平滑全 bin 地图和真实曲线不变，不恢复已排除六组件。

局部统计来源为同目录 `fig3_local_donor_contrasts.R`。它从 `NMF_WORK_ROOT/tables/TableS3_program_annotation.tsv`、工作根下 `inputs/cortex_nmf_program/archived/revisions/2026-08-02_v49_professor_review/human_validation/HUMAN_VALIDATION_spatial_section_by_layer_aggregates_all54.tsv`、上述 `layer_region_profiles.tsv` 以及 `NMF_ARCHIVE_SOURCE_ROOT/scripts/fig2/prog_x_layer_per_chip.tsv` 读取既有值，输出 `NMF_WORK_ROOT/analysis/fig3_layer_region/fig3_local_donor_contrasts.tsv`。R 入口与绘图入口读取注释/身份表的根不同：前者为 `NMF_WORK_ROOT`，后者为 `NMF_SOURCE_ROOT`；两处须提供同一已采用来源的对应文件，不能假定路径变量会自动复制它们。

- **b 的供体内差值**：源 `mean_z` 是切片×域内 bin 分数的中位数；同供体同脑区同域先等权平均切片，再在该供体该脑区计算“目标域−其余六域的等权平均”，随后在供体内等权平均实际观察到且七域完整的脑区。五名供体各贡献一个差值，以双侧单样本 t 对比零；54 程序×7 域的 378 个比较构成独立 BY 家族，不补齐未观察到的脑区。
- **c 的同层区域差值**：每个程序、每层分别以供体×脑区分数拟合 `score ~ region + (1 | donor)`，使用 `lme4::lmer` 的 REML；九个多供体脑区为 AG、DLPFC、FPPFC、M1、S1、SMG、SPL、V1、VLPFC。目标脑区与同层其余八区的等权调整均值比较，直接使用 `pbkrtest::vcovAdj`、`Lb_ddf` 的 Kenward–Roger 协方差与自由度做双侧对比；54×6×9 的 2916 个比较构成另一 BY 家族。不自动改用其它自由度近似；可估计且 KR 有效的边界方差拟合保留。单供体 ACC、ITG、PoCG、S1E、STG 仍显示在热图中，但不进入该推断家族，相应 NA **不是“不显著”结论**。

复用当前主图直接读上述已保留 TSV；只有单独获准重算时才执行统计 R 源。该统计与绘图链不把显示用 Z 值作为新统计响应。本次同步没有重跑上述统计。

#### FigS10 与 Fig3 的早期汇总分支

- 原完整54程序×34片八域profile的r及已存median/切片身份移至**FigS10h**，读取既有值而不重算。FigS10仍是同一PDF两页，**每页170×225mm**；原a–g保留（首页从固定九panel原c/d/e/i重编a–d；四原分布e、P40 f、P48 g），两页含caption。对应 `FigS10.png` 是完整两页纵向预览，依赖Pillow，不是第一页缩图。
- 主图源设置仍170×200mm、图字至少7pt；这是源设置，不是整刊排版合规认证。S10从固定源重建，不读取当前S10再append，因此不会重复追加；无新固定源、S21或S11–S20改动。旧OLIGO距离仅在固定九panel输入保留，不声称在S10。
- 外部只读输入：`NMF_SOURCE_ROOT/tables/TableS3_program_annotation.tsv`、该source根下既有 `HUMAN_VALIDATION_spatial_section_by_layer_aggregates_all54.tsv`；`NMF_ARCHIVE_SOURCE_ROOT/results/crossregion_v1` 下 `spatial_bin50_meta.parquet`、`spatial_bin50_program_score_SCT.parquet`；同archive根 `scripts/fig2` 下 `prog_x_layer_per_chip.tsv`、`prog_x_layer_global.tsv`、`panelg_summary.tsv`、`panelg_reproducibility.tsv`、`repchip_meta.tsv`。数据未复制进代码包。`NMF_ARCHIVE_SOURCE_ROOT` 必须是只读来源根，区别于可写的 `CORTEX_PROGRAM_ROOT`。
- 两张L1–L6表位于 `NMF_WORK_ROOT/analysis/fig3_layer_region/layer_bin_distributions.tsv`、`layer_region_profiles.tsv`；原值与早期 `display_example/illustrated_layer` 字段不改、不扩WM。S10e四例仍对应早期P5/P13/P28/P8快照，P16/P40/P48按当前producer选择。缺失槽保持 `observed=False`、score空，不补0。原生产任务表意为44片/5donor/14区；distribution 14,580行、profile 27,216行，本整理任务未重算。
- `--render-only` 是保留的主图/补图联合绘图分支：读两张原值表，在上述分析目录写 `Fig3.pdf`、`Fig3.png`、`FigS10.pdf`、`FigS10.png`，不写两表，但**会覆盖 Fig3，且不读取局部统计 TSV**，因此不能代替当前主图的 `--local-significance` 入口。无选项默认入口还会生成并写入两张汇总表，不能误作仅重绘命令。

#### Fig2、Fig4、Fig5 的现行及历史入口

以下文件均在 `07_figures_and_tables/current`；“已有源”不代表本次重新执行或重画。

- **Fig2 正式区域面板**：`render_fig2_donor_panels.R` 读取已有 formal21 模型结果和注释；其中 `--within-subclass-only` 分支只生产当前亚类支持面板。新增的 `render_fig2_donor_evidence.py` 保存较早的 donor-robust 图源，不能把其历史八程序显示或 permutation 结果替换成 formal21 三方法定义。其默认入口读取 `NMF_WORK_ROOT/analysis/01_existing_regional_evidence/` 中的 `DONOR_ROBUST_ALL54.tsv`、`DONOR_LOO_STABILITY.tsv`、`HUMAN_VALIDATION_snrna_all54_donor_adjusted_region_profiles.tsv`、`HUMAN_VALIDATION_snrna_donor_region_pseudobulk_all54.tsv`，并读取工作根下 `inputs/current_six_figures/Supplementary_Tables_S1-S6.xlsx` 的 Table S3；相应 PDF/PNG 写入 `figures/human_revision/panels/fig2_donor_evidence/`，用于该历史链及 S18/S19 面板。`--main-panels-only` 排除完整八程序供体点图和方向敏感性面板；`--heatmap-only` 则仅读已存区域 profile、robust 表及 `tables/TableS3_program_annotation.tsv`，输出早期 `Fig2_donor_adjusted_regional_heatmap.pdf/.png`。这些选项均不重新拟合推断或 LOO，也不是自动重建当前 Fig2 的总入口。
- **Fig2 原生内容整理**：`reduce_main_figure_layouts.py --fig2-clean-native-only` 读取并原位更新 `NMF_WORK_ROOT/figures/human_revision/main_figures_pdf/Fig2.pdf`，只去除实际放置区域外内容及被后续不透明白矩形完整遮盖的旧字，保留可见 vector/text、Form、字体、clip 和绘制顺序。它不重画统计面板、不生成 PNG，也**不读取或回写稿件阅读 PDF**；不能把 PDF 内容裁剪当作上游分析。
- **Fig4 当前构图**：`render_fig4_cross_cell_preferences.py --reorganized` 从标准输入读取对应固定 native 底图的 PDF 字节，结合 `NMF_WORK_ROOT/tables/TableS3_program_annotation.tsv`、`TableS6_between_chip_colocalization.tsv` 及现成空间/距离关系输入，写入 `NMF_WORK_ROOT/analysis/fig4_cross_cell_preferences/Fig4_reorganized.pdf/.png`。`--compact-candidate` 是另存的候选布局，不是当前正式构图；默认入口生产跨参考偏好–目标亚类概览，也不等同于 `--reorganized`。
- **Fig5 当前构图**：`render_fig5_revision.py` 同样从标准输入读取对应固定 native 底图，使用工作根的 TableS3、真实供体/整供体留出结果，以及 `CORTEX_PROGRAM_ROOT/results/crossregion_v1/` 下既有 SCT 分数、空间元数据、RCTD、program similarity 和 program–program markcorr 汇总，写 `NMF_WORK_ROOT/analysis/fig5_revision/Fig5_revised.pdf/.png`。Fig4/5 底图须符合原 panel 布局，不能用任意同名或已重排 PDF 顶替。两份核心 producer 本来已与最新生产源同步，本次没有重画。
- **Fig5/S20 的真实供体面板源**：新增 `render_fig5_true_donor_evidence.py` 无动作选项。它读取 `NMF_WORK_ROOT/analysis/02_true_donor_spatial_support/` 的 `same_bin_per_donor_effects.tsv`、`same_bin_true_donor_loo.tsv`、`P16_distance_0_500_per_donor_effects.tsv`、`P16_distance_0_500_true_donor_loo.tsv`、`P16_section_ring_effects.tsv`，并从上述旧工作簿 Table S3 读取注释，写入 `figures/human_revision/panels/fig5_true_donor_evidence/`。输出 stem 为 `Fig5_P16_same_bin`、`Fig5_P16_distance_support`、`Fig5_P16_section_distance_profiles`、`Fig5_P16_section_points`，各含 PDF/PNG；保留同 bin 与 0–500 µm 的区别，不重新计算效应或 LOO。整供体留出范围是影响敏感性，不是置信区间，也不改变当前 Fig5 的例证选择。

#### 其它保留图源

- `reference_panels`：发现集usage的中间表、原60维描述嵌入与Fig1源panel。固定150k barcode选择依赖已有 `umap_embedding.csv`，不能在缺失时临时换采样。SLAT嵌入不替代正式54程序矩阵。
- `within_subclass_S7`：旧FigS4谱系，现稿S7相关材料。`fig4_prep_c.py` 用固定全局亚类比例/全局亚类均值的两类counterfactual，两个加权方差归一化再乘原η²；`fig4_prep.py` 的20次bootstrap是**核水平**，不是donor、更不是cNMF重复。
- `rank_S1/rank_S2`：既有rank诊断源。其已保存的高K/信息准则表是读取输入；不声称缺失历史重复已由这两个脚本补齐。
- `spatial_legacy_panels`：旧空间图源；`01_prepare_metadata.py` 只保留原01_probe的生产性cache导出，取消探查输出。`02_aggregate.py` 中 `mean_z` 实为bin中位数；`05_panel_fgh_data.py` 的抽样/域汇总保留，不与均值GEE混淆。
- `current`：standalone现有Fig1/2/4/5修订入口、外部图及布局。`reduce_main_figure_layouts.py` 依赖原 `inputs/current_six_figures`、其他原panel PDF与tables；它包含多个历史布局动作，**不是“全部重画最新版”的总入口**。旧Fig3分支不能覆盖本轮新的Fig3。

本包有意保留实际分析与原脚本名称，不把旧图号自动批量重命名成当前图号；最终版构图可继续读取既有vector panel，代码没有伪造缺失原图层或将PDF裁剪冒充上游统计。

## 本次对副本的调整

既有副本整理保留路径环境变量/参数、同目录代码引用和 SCT 表示说明，取消危险的旧结果递归删除及不必要记录；R 内存上限为 50 GiB。原独立 dream 尺度比较拟合已移除，实际采用的 response×1000 模型与结果生成保留。最新同步仅补入已存在的历史 rank/供体图源、更新 Fig3 已批准的生产源，并将 Fig2 native 内容处理隔离为源图入口；没有新设科学阈值、模型、分母、映射或计算并行参数。本包不存在“运行所有入口即可重建最新版”的总控脚本。

GO:BP 的来源边界已进一步明确。历史 GO 表合并、54 程序筛选及绘表脚本已纳入 `07_figures_and_tables/gobp_program_table/`。上游富集计算 producer 与当前 TableS3 的精确最终装配脚本仍未定位。已找到的129程序 prerank 脚本属于另一分析对象，未混入。`make_program_names.py` 自带警告，自动全名会覆盖后来人工修正，故不是现稿名称的再生成入口。现稿采用 TableS3 及固定 program_names，不盲跑历史命名脚本。

原始文件→合并 RDS 的 QC/merge 以及上述 GO:BP 来源链仍有未定位环节。历史80次结果位于另一台机器，公开单次任务100次初始化路线不以迁移这批结果为前提。各下游入口仍需外置的原始数据、固定参考、既有分析结果和 native 图形资产；这些输入依赖与源码缺口分别保留。本包交付了当前可追踪的代码及使用边界，**不声称已完成新一轮100次任务、端到端复现或运行测试**。
