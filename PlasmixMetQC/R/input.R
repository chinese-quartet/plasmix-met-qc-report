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
  expr <- fread(exp_path)
  meta <- fread(meta_path)
  
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
  
  # Handle duplicated Features -------------------------------
  if (any(duplicated(expr$Feature))) {
    dup_features <- unique(expr$Feature[duplicated(expr$Feature)])
    message(sprintf("Found %d duplicated Feature(s). Keeping rows with larger mean values.", length(dup_features)))
    
    # Calculate row means for numeric columns
    numeric_cols <- 2:ncol(expr)
    expr$row_mean <- rowMeans(expr[, numeric_cols], na.rm = TRUE)
    
    # For each duplicated feature, keep only the row with max mean
    keep_rows <- sapply(1:nrow(expr), function(i) {
      feat <- expr$Feature[i]
      if (feat %in% dup_features) {
        # Find all rows with this feature
        dup_indices <- which(expr$Feature == feat)
        # Keep only if this row has the maximum mean among duplicates
        return(i == dup_indices[which.max(expr$row_mean[dup_indices])])
      } else {
        return(TRUE)
      }
    })
    
    expr <- expr[keep_rows, ]
    expr$row_mean <- NULL  # Remove temporary column
    
    message(sprintf("Retained %d unique features after removing duplicates.", nrow(expr)))
  }
  
  # 检查 Library 列匹配
  expr_libs <- colnames(expr)[2:ncol(expr)]
  missing_libs <- setdiff(meta_final$library, expr_libs)
  if (length(missing_libs) > 0) {
    stop(sprintf(
      "❌ 样本信息不匹配！\n\n以下 %d 个样本在样本信息表（metadata）中存在，但在表达数据（expression）中缺失：\n%s\n\n请检查：\n  1. 样本信息表中的 library 列名是否正确\n  2. 表达数据中的列名是否与 library 一致",
      length(missing_libs),
      paste("  -", head(missing_libs, 10), collapse = "\n")
    ))
  }
  # 同时添加反向检查（表达数据中有但元数据中没有的样本）：
  extra_libs <- setdiff(expr_libs, meta_final$library)
  if (length(extra_libs) > 0) {
    warning(sprintf(
      "⚠️  以下 %d 个样本在表达数据中存在，但在样本信息表中缺失，将被忽略：\n%s",
      length(extra_libs),
      paste("  -", head(extra_libs, 10), collapse = "\n")
    ))
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
