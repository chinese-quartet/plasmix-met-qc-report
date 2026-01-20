```markdown
# PlasmixMetaboQC: Metabolomics Quality Control Toolkit
**PlasmixMetaboQC** 是专为 Plasmix 计划开发的代谢组学数据质量控制（QC）R 包。该工具基于比例（Ratio-based）定量策略，利用 Plasmix 血浆参考物质（M, Y, P, X, F）评估代谢组数据的信噪比（SNR）和与参考数据集的一致性（PCC）。

该包旨在简化非靶向和靶向代谢组学数据的质控流程，支持一键生成中文 Word 报告。

## ✨ 核心功能

* **平台无关性 (Platform Agnostic)**：适用于各类代谢组学检测平台，无需指定平台参数，统一进行 Scale 归一化处理。
* **智能数据清洗**：
    * **自动去重**：自动识别并重命名重复的 Feature ID（如 `HMDB001` -> `HMDB001.1`），防止分析报错。
    * **旧格式兼容**：自动识别并移除包含 `Type` 的冗余列。
* **关键指标计算**：
    * **SNR (Signal-to-Noise Ratio)**：评估样本组间生物学差异（信号）与技术重复变异（噪声）的比值。
    * **PCC (Pearson Correlation Coeffiency)**：评估实测数据的比值指纹与标准参考数据集的 Pearson 相关性。
* **自动化报告**：直接生成包含散点图、PCA 图和详细结论的 `.docx` 格式质控报告。

## 📦 安装说明

请确保您已安装 `devtools` 包。您可以通过以下方式安装 `PlasmixMetaboQC`：

```r
# 方法 1: 如果项目已上传至 GitHub (请替换您的用户名)
# devtools::install_github("https://github.com/markx945/Plasmix_Metabolite_QC_Report/tree/main/PlasmixMetQC")

# 方法 2: 如果您在本地项目根目录下
devtools::install_local(".")

```

## 📂 数据准备

运行 QC 流程需要准备两个输入文件：

### 1. 表达矩阵 (Expression Matrix)

* **格式**：`.txt`, `.csv`, 或 `.tsv`
* **结构**：
* 第一列必须是 **Feature ID**（列名推荐为 `Feature`，也支持 `Metabolite` 或 `Compound`）。
* 其余列为样本的定量数据（建议为 Log2 转换后的值）。

* **重复处理**：如果存在重复的 Feature ID，程序会自动添加后缀（如 `.1`, `.2`）进行区分。

| Feature | Sample_01 | Sample_02 | Sample_03 | ... |
| --- | --- | --- | --- | --- |
| HMDB00001 | 14.2 | 13.8 | 14.1 | ... |
| HMDB00002 | 10.5 | 10.1 | 10.9 | ... |

### 2. 元数据 (Metadata)

* **格式**：`.txt`, `.csv`, 或 `.tsv`
* **必需列**（列名不区分大小写）：
* `library`: 必须与表达矩阵的列名（样本ID）完全一致。
* `sample`: 样本组别名称。必须包含 **P** 样（作为参考分母），以及 **M, Y, X, F** 中的任意组合。

* **可选列**：
* `platform`: 不再强制要求。


| library | sample |
| --- | --- |
| Sample_01 | P |
| Sample_02 | M |
| Sample_03 | Y |

## 🚀 快速上手

以下代码展示了如何从读取数据到生成报告的完整流程：

```r
library(PlasmixMetaboQC)

# 1. 设置文件路径
expr_file  <- "./data/expression_matrix.txt"
meta_file  <- "./data/metadata.txt"
output_dir <- "./output"

# 2. 运行核心分析
# 此步骤会计算 SNR 和 RC，并生成中间图表
qc_res <- qc_conclusion(
  exp_path = expr_file, 
  meta_path = meta_file, 
  output_dir = output_dir,
  plot = TRUE  # 生成 PCA 和 相关性图
)

# 3. 查看简要结果
print(qc_res$qc_metrics_table)

# 4. 生成中文 Word 报告
generate_metabo_report(
  qc_result = qc_res,            # 上一步的分析结果对象
  report_dir = output_dir,       # 报告保存目录
  report_name = "Metabo_QC_Report.docx",
  batch_name = "Batch_2026_01"   # 指定批次名称（显示在报告中）
)

```

## 📊 结果解读

自动生成的报告将依据以下标准对数据质量进行判定：

| 指标 (Metric) | 阈值 (Threshold) | 说明 |
| --- | --- | --- |
| **信噪比 (SNR)** | **≥ 5** | 数值越高，表示组间差异信号越显著，技术噪音越低。 |
| **相对相关性 (PCC)** | **≥ 0.80** | 数值越高，表示实测数据的相对定量趋势与标准参考值越一致。 |
| **整体质量** | **Pass** | 当 **SNR** 达标时，判定为通过。 PCC仅作为参考|

> **注意**：若报告中的数值旁出现红色箭头 (↓)，表示该指标未达到推荐标准。

## 🛠️ 常见问题 (Troubleshooting)

**Q: 报错 `Could not load 'reference_dataset' from PlasmixMetQC**`

* **原因**：代码中引用的包名与实际安装的包名不一致。
* **解决**：请检查 `DESCRIPTION` 文件中的 `Package:` 字段，确保其为 `PlasmixMetaboQC`，并重新安装。

**Q: 警告 `Duplicates detected in Feature IDs...**`

* **原因**：表达矩阵中存在完全相同的代谢物名称。
* **解决**：无需手动处理。本包内置了自动清洗功能，会自动将其重命名为唯一 ID（例如 `Name.1`）并继续运行。

## 📚 参考文献

1. Zheng Y, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. *Nature Biotechnology*, 2024.
2. 上海临床队列组学检测工作指引（征求意见稿）, 2025.

---
