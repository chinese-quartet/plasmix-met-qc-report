library(PlasmixMetQC)

# 2. 定义文件路径
# 请确保这两个文件在你的工作目录，或者填写绝对路径
expr_file <- system.file("extdata", "plasmix_met_test_expr.csv", package = "PlasmixMetQC")
meta_file <- system.file("extdata", "plasmix_met_test_meta.csv", package = "PlasmixMetQC")

template <- system.file("extdata", "Plasmix_template.docx", package = "PlasmixMetQC")

qc_result <- qc_conclusion(
  exp_path = expr_file,
  meta_path = meta_file,
)

generate_metabo_report(
  qc_result = qc_result,
  report_template = template,
  # report_dir = output_dir
  # report_name = "我的测试报告.docx"
  # batch_name = "20260108测试批次" # 可选，指定表格里的批次名
)
