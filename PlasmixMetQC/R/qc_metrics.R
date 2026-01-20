#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @importFrom ggplot2 geom_abline
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample

  d[d == 0] <- NA

  # 1. Feature Count (Internal check only)
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)

  # 2. Missing %
  d_vals <- d[, 2:ncol(d)]
  d_all_num <- nrow(d_vals) * ncol(d_vals)
  d_missing_num <- length(which(is.na(d_vals)))
  prop_missing <- 0
  if (d_all_num > 0) prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)

  # 3. Correlation & CV
  samples <- table(s)
  rep_samples <- samples[samples > 1]

  stat_acor <- NA
  stat_cv <- NA

  if (length(rep_samples) > 0) {
    d_mtx <- d[, 2:ncol(d)]
    common_libs <- intersect(colnames(d_mtx), m$library)
    if (length(common_libs) >= 2) {
      d_mtx <- d_mtx[, common_libs]
      cor_res <- cor(d_mtx, use = "pairwise.complete.obs", method = "pearson")
      stat_acor <- round(median(cor_res[upper.tri(cor_res)], na.rm = TRUE), 3)
    }

    # CV Calculation
    d_long <- melt(d, id.vars = colnames(d)[1])
    d_long <- na.omit(d_long)
    d_long <- merge(d_long, m, by.x = "variable", by.y = "library")

    # 假设输入是 Log2，转回 Linear 算 CV
    d_long$value_linear <- 2^(d_long$value)

    d_cv <- aggregate(
      value_linear ~ Feature + sample,
      data = d_long,
      FUN = function(x) sd(x) / mean(x)
    )
    stat_cv <- round(median(d_cv$value_linear, na.rm = T) * 100, 3)
  }

  return(c(stat_num, stat_missing, stat_acor, stat_cv))
}

#' Get color mapping
#' @export
get_sample_colors <- function(samples) {
  cols <- c("M" = "#3171b8", "Y" = "#53a949", "P" = "#6d3390", "X" = "#ffc65d", "F" = "#ae231c")
  res <- cols[samples]
  res[is.na(res)] <- "gray50"
  return(res)
}

#' Calculating SNR
#' @import stats
#' @import utils
#' @importFrom data.table data.table setkey
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme labs guides guide_legend ggsave scale_x_continuous scale_y_continuous
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # 数据准备
  # [修复] 处理重复行名问题
  feat_names <- expr_dt[[1]]
  if (any(duplicated(feat_names))) {
    warning("Duplicates detected in Feature IDs during SNR calculation. Making them unique (e.g., ID.1).")
    feat_names <- make.unique(as.character(feat_names))
  }

  expr_df <- data.frame(expr_dt[, 2:ncol(expr_dt)], row.names = feat_names)
  expr_df[is.na(expr_df)] <- 0

  libs <- colnames(expr_df)
  meta_sub <- meta_dt[match(libs, meta_dt$library), ]
  group <- meta_sub$sample

  # 代谢组默认 Scale
  do_scale <- TRUE

  # PCA
  expr_t <- t(expr_df)
  vars <- apply(expr_t, 2, var)
  expr_t <- expr_t[, vars > 0 & !is.na(vars), drop = FALSE]

  if (ncol(expr_t) < 2) stop("Too few features for PCA.")

  pca <- prcomp(expr_t, scale. = do_scale)

  pcs <- as.data.frame(pca$x)
  pcs$sample <- group

  # 计算 SNR
  prop <- summary(pca)$importance[2, 1:2]

  n <- nrow(pcs)
  dist_val <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      d2 <- prop[1] * (pcs[i, 1] - pcs[j, 1])^2 + prop[2] * (pcs[i, 2] - pcs[j, 2])^2
      dist_val[i, j] <- d2
    }
  }

  intra_dists <- c()
  inter_dists <- c()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (group[i] == group[j]) {
        intra_dists <- c(intra_dists, dist_val[i, j])
      } else {
        inter_dists <- c(inter_dists, dist_val[i, j])
      }
    }
  }

  snr <- mean(inter_dists) / mean(intra_dists)
  snr_db <- round(10 * log10(snr), 2)

  # Plot
  p <- NULL
  if (plot) {
    cols <- get_sample_colors(unique(group))
    p <- ggplot(pcs, aes(x = PC1, y = PC2, color = sample)) +
      geom_point(size = 6) +
      scale_color_manual(values = cols) +
      theme_few() +
      labs(
        title = paste0("SNR = ", snr_db),
        x = sprintf("PC1 (%.1f%%)", prop[1] * 100),
        y = sprintf("PC2 (%.1f%%)", prop[2] * 100)
      )
  }

  if (!is.null(output_dir)) {
    if (plot) ggsave(file.path(output_dir, "pca_plot.png"), p, width = 6, height = 5)
    write.table(pcs, file.path(output_dir, "pca_table.tsv"), sep = "\t", quote = F)
  }

  return(list(SNR = snr_db, snr_plot = p))
}

#' Calculating RC
#' @importFrom utils data
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic labs coord_fixed
#' @export
qc_cor <- function(expr_dt, meta_dt, output_dir = NULL, plot = FALSE, show_sample_pairs = TRUE) {
  tryCatch(
    {
      utils::data("reference_dataset", package = "PlasmixMetQC", envir = environment())
    },
    error = function(e) {
      stop("Could not load 'reference_dataset' from PlasmixMetQC.")
    }
  )

  if (!exists("reference_dataset")) stop("reference_dataset not found.")

  ref_dt <- reference_dataset
  if (nrow(ref_dt) == 0) warning("Reference dataset is empty.")

  # [修复] 处理重复行名问题
  mat <- as.matrix(expr_dt[, -1])

  feat_names <- expr_dt[[1]]
  if (any(duplicated(feat_names))) {
    warning("Duplicates detected in Feature IDs during RC calculation. Making them unique.")
    feat_names <- make.unique(as.character(feat_names))
  }

  rownames(mat) <- feat_names
  mat <- mat[, meta_dt$library]

  mat_linear <- 2^mat

  samps <- unique(meta_dt$sample)
  ref_samp <- "P"
  if (!"P" %in% samps) {
    ref_samp <- samps[1]
    warning("No 'P' sample. Using ", ref_samp)
  }
  test_samps <- setdiff(samps, ref_samp)

  res_all <- data.frame()

  for (ts in test_samps) {
    cols_test <- which(meta_dt$sample == ts)
    cols_ref <- which(meta_dt$sample == ref_samp)

    dat_ref <- mat_linear[, cols_ref, drop = FALSE]
    dat_test <- mat_linear[, cols_test, drop = FALSE]

    count_ref <- rowSums(!is.na(dat_ref))
    count_test <- rowSums(!is.na(dat_test))
    keep_rows <- (count_ref > 0) & (count_test > 0)

    current_features <- rownames(dat_ref)

    pair_name <- paste(ts, ref_samp, sep = "/")
    feats_ref <- ref_dt$feature[ref_dt$sample_pair == pair_name]
    common_feats <- intersect(current_features, feats_ref)

    if (length(common_feats) < 3) next

    dat_ref <- dat_ref[common_feats, , drop = FALSE]
    dat_test <- dat_test[common_feats, , drop = FALSE]

    mean_ref <- rowMeans(dat_ref, na.rm = TRUE)
    mean_test <- rowMeans(dat_test, na.rm = TRUE)

    ratio_val <- mean_test / mean_ref

    df_tmp <- data.frame(
      Feature = common_feats,
      Sample.Pair = pair_name,
      logFC.Test = ratio_val
    )

    res_all <- rbind(res_all, df_tmp)
  }

  if (nrow(res_all) == 0) {
    return(list(COR = NA, cor_plot = NULL))
  }

  res_all$key <- paste(res_all$Feature, res_all$Sample.Pair)
  ref_dt$key <- paste(ref_dt$feature, ref_dt$sample_pair)

  merged <- merge(res_all, ref_dt, by = "key")

  cor_val <- cor(merged$logFC.Test, merged$assigned_value, use = "complete.obs")
  cor_val <- round(cor_val, 3)

  p <- NULL
  if (plot) {
    unique_pairs <- unique(merged$Sample.Pair)
    # numerators <- unique(sapply(strsplit(as.character(unique_pairs), "/"), `[`, 1))
    # sample_cols_base <- get_sample_colors(numerators)
    # 1. 提取所有涉及的分子样本 (Test Samples)
    numerators <- unique(sapply(strsplit(as.character(unique_pairs), "/"), `[`, 1))

    # 2. 获取这些样本的标准颜色
    sample_cols_base <- get_sample_colors(numerators) # named vector

    # 3. 构建 Pair -> Color 的映射向量
    pair_colors_vec <- c()
    for (pair in unique_pairs) {
      num <- strsplit(as.character(pair), "/")[[1]][1]
      # 安全查找：如果 num 存在于 sample_cols_base 的名字中
      if (num %in% names(sample_cols_base)) {
        pair_colors_vec[pair] <- sample_cols_base[num]
      } else {
        pair_colors_vec[pair] <- "gray"
      }
    }

    p <- ggplot(merged, aes(x = assigned_value, y = logFC.Test)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(aes(color = Sample.Pair), alpha = 0.6, size = 2) +
      scale_color_manual(values = pair_colors_vec) +
      theme_few() +
      coord_fixed(xlim = c(0, 5), ylim = c(0, 5)) +
      scale_x_continuous(breaks = 0:5) +
      scale_y_continuous(breaks = 0:5) +
      labs(
        title = paste0("PCC = ", cor_val),
        x = "Reference Value",
        y = "Test Value"
      )
  }

  return(list(COR = cor_val, cor_plot = p, logfc = merged))
}
