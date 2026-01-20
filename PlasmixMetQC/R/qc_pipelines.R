#' Calculating all QC metrics
#' @param expr_dt A generic expression table (Feature level)
#' @param meta_dt A metadata file
#' @param output_dir A directory for results
#' @param plot if True, a plot will be output.
#' @export
qc_allmetrics <- function(expr_dt, meta_dt,
                          output_dir = NULL, plot = TRUE) {
  # 1. 基础信息
  pro_info <- qc_info(expr_dt, meta_dt)

  # 2. 计算 SNR
  snr_results <- qc_snr(expr_dt, meta_dt, output_dir, plot)
  snr_value <- snr_results$SNR

  # 3. 计算 RC
  cor_results <- qc_cor(expr_dt, meta_dt, output_dir, plot)
  cor_value <- cor_results$COR

  # 4. 组装结果 (Only SNR & RC)
  metrics <- c(
    "Signal-to-Noise Ratio (SNR)",
    "Relative Correlation with Reference Datasets (RC)"
  )

  qc_values <- c(snr_value, cor_value)

  output_table <- data.frame(
    "Quality_Metrics" = metrics,
    "Value" = qc_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  all_results <- list(
    snr_results = snr_results,
    cor_results = cor_results,
    output_table = output_table
  )

  return(all_results)
}
