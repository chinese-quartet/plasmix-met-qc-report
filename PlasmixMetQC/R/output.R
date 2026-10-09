#' Generating a table of conclusion
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @importFrom utils write.table
#' @export
qc_conclusion <- function(exp_path, meta_path, output_dir = NULL, plot = TRUE) {
  # 1. Load the input data
    # 🔧 添加输出目录检查
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        message(sprintf("📁 输出目录不存在，正在创建：%s", output_dir))
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        if (!dir.exists(output_dir)) {
          stop(sprintf(
            "❌ 无法创建输出目录：%s\n请检查：\n  1. 路径是否有写入权限\n  2. 路径格式是否正确",
            output_dir
          ))
        }
      }
    }
    
  data_list <- input_data(exp_path, meta_path)

  expr_dt <- data_list$expr_dt
  meta <- data_list$metadata

  # 2. Run the QC pipelines
  allmetrics_results <- qc_allmetrics(expr_dt, meta, output_dir, plot)

  # 提取结果表
  output_table <- allmetrics_results$output_table

  # 3. Save & Output
  if (!is.null(output_dir)) {
    output_dir_final <- file.path(output_dir, "conclusion_table.tsv")
    write.table(output_table, output_dir_final, sep = "\t", row.names = F, quote = F)
  }

  final_list <- list(
    qc_metrics_table = output_table,
    snr_plot = allmetrics_results$snr_results$snr_plot,
    cor_plot = allmetrics_results$cor_results$cor_plot,
    detail_results = allmetrics_results,
    input_meta = meta
  )

  return(final_list)
}