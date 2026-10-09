#' @title Generate Metabolomics Quality Report (Chinese Version)
#' @description Generates a report following the ShenKang Plasmix Metabolomics template structure.
#' @param qc_result list, the output from qc_conclusion()
#' @param report_template character, path to a .docx file to use as style template (optional)
#' @param report_dir character, directory to save the report
#' @param report_name character, filename of the report
#' @param batch_name character, optional, manual batch name. If NULL, tries to infer from metadata.
#' @importFrom dplyr %>%
#' @importFrom flextable flextable theme_box align width bold bg color body_add_flextable set_header_labels
#' @importFrom officer read_docx body_add_par body_add_gg body_add_break body_add_fpar fpar fp_text prop_section
#' @export
generate_metabo_report <- function(qc_result,
                                   report_template = NULL,
                                   report_dir = NULL,
                                   report_name = NULL,
                                   batch_name = NULL) {
  # --- 1. 参数检查与初始化 ---
  if (is.null(qc_result)) stop("qc_result is required.")

  if (is.null(report_dir)) {
    report_dir <- file.path(getwd(), "output")
    message(sprintf("📁 未指定报告输出目录，使用默认路径：%s", report_dir))
  }
  if (!dir.exists(report_dir)) {
    message(sprintf("📁 输出目录不存在，正在创建：%s", report_dir))
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(report_dir)) {
      stop(sprintf(
        "❌ 无法创建报告输出目录：%s\n请检查：\n  1. 路径是否有写入权限\n  2. 磁盘空间是否充足",
        report_dir
      ))
    }
  }
  if (is.null(report_name)) {
    report_name <- paste0("Plasmix_Metabo_QC_Report_", format(Sys.Date(), "%Y%m%d"), ".docx")
  }
  output_file <- file.path(report_dir, report_name)

  # --- 2. 提取并处理数据 ---
  raw_table <- qc_result$qc_metrics_table
  meta <- qc_result$input_meta

  if (is.null(batch_name)) {
    if ("study_id" %in% colnames(meta) && length(unique(meta$study_id)) == 1) {
      batch_name <- unique(meta$study_id)[1]
    } else {
      batch_name <- paste0(format(Sys.Date(), "%Y%m%d"), "_data")
    }
  }

  # 提取数值 (Only SNR and RC)
  snr_val <- raw_table$Value[grep("Signal-to-Noise Ratio", raw_table$Quality_Metrics)]
  rc_val <- raw_table$Value[grep("Relative Correlation", raw_table$Quality_Metrics)]

  val_snr <- as.numeric(snr_val)
  val_rc <- as.numeric(rc_val)

  # --- 3. 质量判定与格式化 ---

  # 3.1 SNR (Target >= 5)
  if (length(val_snr) == 0 || is.na(val_snr)) {
    txt_snr <- "-"
    pass_snr <- FALSE
  } else {
    txt_snr <- sprintf("%.2f", val_snr)
    if (val_snr < 5) {
      txt_snr <- paste0(txt_snr, " ↓")
      pass_snr <- FALSE
    } else {
      pass_snr <- TRUE
    }
  }

  # 3.2 RC (Target >= 0.80)
  if (length(val_rc) == 0 || is.na(val_rc)) {
    txt_rc <- "-"
    # pass_rc <- FALSE
  } else {
    txt_rc <- sprintf("%.3f", val_rc)
  }
  #   if (val_rc < 0.80) {
  #     txt_rc <- paste0(txt_rc, " ↓")
  #     pass_rc <- FALSE
  #   } else {
  #     pass_rc <- TRUE
  #   }
  # }

  # 3.3 整体判定 (SNR 通过即可)
  # is_all_pass <- pass_snr && pass_rc
  is_all_pass <- pass_snr
  overall_quality <- ifelse(is_all_pass, "Yes", "No")

  # --- 4. 构建表格数据 ---
  cols <- c("样本组", "信噪比", "Pearson相关系数", "是否通过")

  # 第一行：推荐标准
  row_std <- c("推荐质量标准", "≥5", "≥0.80", "-")

  # 第二行：实际数据
  row_dat <- c(batch_name, txt_snr, txt_rc, overall_quality)

  df_rep <- data.frame(rbind(row_std, row_dat), stringsAsFactors = FALSE)
  colnames(df_rep) <- cols

  # 创建 Flextable
  ft <- flextable(df_rep) %>%
    theme_box() %>%
    flextable::font(part = "all", fontname = "Times New Roman") %>%
    align(align = "center", part = "all") %>%
    width(width = 1.2) %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#EFEFEF") %>%
    color(i = 1, color = "gray40") %>%
    bold(i = 2) %>%
    # 列索引: 2=SNR, 3=RC, 4=Quality
    color(i = 2, j = 2, color = ifelse(grepl("↓", txt_snr), "red", "black")) %>%
    color(i = 2, j = 3, color = ifelse(grepl("↓", txt_rc), "red", "black")) %>%
    color(i = 2, j = 4, color = ifelse(overall_quality == "Yes", "black", "red"))

  # --- 5. 文档文本内容 ---
  txt_summary <- "本报告基于多项组学关键质量控制指标，总结了 Plasmix血浆参考物质所生成代谢组数据的质量情况。质量控制流程从用户输入血浆代谢组表达矩阵开始，分别计算外部质控品的信噪比 (Signal-to-Noise Ratio, SNR)、与参考数据集的Pearson相关系数（Pearson correlation coefficient, PCC）及整体质量判断。"

  txt_def_title <- "质量控制指标"
  txt_def_1 <- "信噪比（Signal-to-Noise Ratio, SNR）：用于刻画某一检测平台、实验室或批次区分不同生物样本组之间内在代谢谱差异（“信号”）与同一样本组技术重复变异（“噪声”）的能力。"
  txt_def_2 <- "与参考数据集的Pearson相关系数（Pearson correlation coefficient, PCC）：定义为在给定样本对之间，测试数据集中比值型代谢表达水平与对应比值型参考数据集之间的 Pearson 相关系数，用于表征比值表达谱在数值层面的整体一致性趋势。"

  txt_refs <- c(
    "1. Zheng Y, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024.",
    "2. Liu, Y. et al. Harmonizing plasma proteomics data with the sample-to-reference ratio approach. BioRxiv (2026)",
    "3. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26"
  )

  txt_disclaimer <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。"

  # --- 6. 生成 Word 文档 ---

  if (!is.null(report_template) && file.exists(report_template)) {
    doc <- read_docx(path = report_template)
  } else {
    doc <- read_docx()
  }

  doc <- doc %>%
    body_add_par("Plasmix血浆代谢组质量报告", style = "heading 1") %>%
    body_add_par("", style = "Normal")

  doc <- doc %>%
    body_add_par("摘要", style = "heading 2") %>%
    body_add_par(txt_summary, style = "Normal") %>%
    body_add_par("", style = "Normal")

  doc <- doc %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")

  doc <- doc %>%
    body_add_par(txt_def_title, style = "heading 2") %>%
    body_add_par(txt_def_1, style = "Normal") %>%
    body_add_par(txt_def_2, style = "Normal") %>%
    body_add_par("", style = "Normal")

  doc <- doc %>%
    body_add_par("参考文献", style = "heading 2")
  for (ref in txt_refs) {
    doc <- doc %>% body_add_par(ref, style = "Normal")
  }
  doc <- doc %>% body_add_par("", style = "Normal")

  doc <- doc %>%
    body_add_par("免责声明", style = "heading 2") %>%
    body_add_par(txt_disclaimer, style = "Normal") %>%
    body_add_break()

  # # Plots
  # doc <- doc %>%
  #   body_add_par("Signal-to-Noise Ratio", style = "heading 2") %>%
  #   body_add_gg(value = qc_result$snr_plot, style = "centered", width = 6, height = 5) %>%
  #   body_add_par("", style = "Normal")
  #
  # doc <- doc %>%
  #   body_add_par("Pearson Correlation Coefficient", style = "heading 2") %>%
  #   body_add_gg(value = qc_result$cor_plot, style = "centered", width = 6, height = 5) %>%
  #   body_add_par("", style = "Normal")
  # Plots

  doc <- doc %>%
    body_add_par("Signal-to-Noise Ratio", style = "heading 2")

  if (inherits(qc_result$snr_plot, "gg")) {
    doc <- doc %>%
      body_add_gg(value = qc_result$snr_plot, style = "centered", width = 6, height = 5)
  } else {
    doc <- doc %>%
      body_add_par("SNR plot unavailable.", style = "Normal")
  }

  doc <- doc %>%
    body_add_par("", style = "Normal") %>%
    body_add_par("Pearson Correlation Coefficient", style = "heading 2")

  if (inherits(qc_result$cor_plot, "gg")) {
    doc <- doc %>%
      body_add_gg(value = qc_result$cor_plot, style = "centered", width = 6, height = 5)
  } else {
    cor_text <- NA_character_
    if (!is.null(raw_table) &&
      "Quality_Metrics" %in% colnames(raw_table) &&
      "Value" %in% colnames(raw_table)) {
      cor_text <- raw_table$Value[grep("Relative Correlation", raw_table$Quality_Metrics)]
      cor_text <- if (length(cor_text) > 0) as.character(cor_text[1]) else NA_character_
    }

    msg_rc <- if (is.na(cor_text) || cor_text == "NA") {
      "PCC unavailable (COR = NA). Possible reason: insufficient matched features with reference dataset."
    } else {
      paste0("PCC plot unavailable (COR = ", cor_text, ").")
    }

    doc <- doc %>%
      body_add_par(msg_rc, style = "Normal")
  }

  doc <- doc %>%
    body_add_par("", style = "Normal")

  print(doc, target = output_file)
  message(paste("Metabolomics Report generated successfully at:", output_file))
}
