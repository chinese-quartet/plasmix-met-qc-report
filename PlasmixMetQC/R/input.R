#' Input files
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @import stats
#' @import utils
#' @importFrom data.table fread
#' @importFrom dplyr %>%
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @export
input_data <- function(exp_path, meta_path) {
  # Load data ------------------------------------------------
  # expr <- fread(exp_path)
  # meta <- fread(meta_path)
  expr <- read.table(exp_path, sep = ",", header = TRUE, na.strings = c("#num!", "#NUM!" , "N/A", "MISSING", "NA"), comment.char = "", check.names = FALSE)
  meta <- read.table(meta_path, sep = ",", header = TRUE, check.names = FALSE)

  expr <- as.data.frame(expr)
  meta <- meta %>%
    rename_with(tolower)

  # Check metadata format ------------------------------------
  # [修改] platform 不再是必须列
  req_cols <- c("library", "sample")

  if ("name" %in% colnames(meta)) {
    meta <- meta %>% dplyr::rename(library = name)
  }

  if (!all(req_cols %in% colnames(meta))) {
    stop('The columns named "library" and "sample" are required in metadata.')
  }

  # [修改] 如果没有 platform 列，自动补充一个默认值，保证后续流程兼容性
  if (!"platform" %in% colnames(meta)) {
    meta$platform <- "Metabolomics"
  }

  # 选择列
  meta_final <- meta %>% select(library, sample, platform)

  if (any(duplicated(meta_final$library))) {
    stop("Duplicated library IDs in metadata.")
  }

  # Check expression data format -----------------------------

  # 兼容旧格式 (Type + Feature/Metabolite)
  if (ncol(expr) >= 2) {
    col1 <- colnames(expr)[1]
    col2 <- colnames(expr)[2]

    if (col1 == "Type" && (col2 == "Feature" || col2 == "Metabolite" || col2 == "Compound")) {
      message("Detected legacy format with 'Type' column. Removing 'Type'.")
      expr <- expr[, -1]
    }
  }

  if (length(which(duplicated(colnames(expr))))) {
    stop("Duplicated column names in data.")
  }

  # 确定 Feature 列
  feature_col_name <- colnames(expr)[1]

  if (feature_col_name != "Feature") {
    message(sprintf("Renaming the first column '%s' to 'Feature' (Metabolite ID).", feature_col_name))
    colnames(expr)[1] <- "Feature"
  }

  # 检查 Library 列匹配
  expr_libs <- colnames(expr)[2:ncol(expr)]

  missing_libs <- setdiff(meta_final$library, expr_libs)
  if (length(missing_libs) > 0) {
    warning(sprintf("%d libraries in metadata are missing from expression data.", length(missing_libs)))
  }

  valid_libs <- intersect(expr_libs, meta_final$library)

  if (length(valid_libs) == 0) {
    stop("No common library IDs found between expression data and metadata.")
  }

  expr_final <- expr[, c("Feature", valid_libs)]

  data_list <- list(
    "expr_dt" = expr_final,
    "metadata" = meta_final
  )

  return(data_list)
}
