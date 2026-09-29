# scripts/plot/heatmap.R — 表达热图（P3）
## 规范：行 Z-score + RdBu 发散色尺（红=上调语义与火山图同源）；
## 行名显隐按标签预算算术判定（画布高度能容纳几行就显示几行），
## 决策与算式落 label_decisions.csv —— 不硬编码 show_rownames。
plot_heatmap <- function(mat, path, title = 'Heatmap',
                         fontsize_row = 5, w = W_ONE_HALF, h = mm(110),
                         label_decisions_path = NULL) {
  m <- as.matrix(mat)
  z <- t(scale(t(m)))
  z[!is.finite(z)] <- 0
  z <- pmin(pmax(z, -3), 3)   # 夹到 ±3 防单点拉爆色尺

  # ---- 标签预算（老框架算术判据）：面板可容纳行数 = 高(in)*72*panel_frac / (字号+min_gap)
  panel_frac <- 0.82          # pheatmap 无副标题、右侧细色条
  min_gap    <- 2.5           # pt；约 10px @300dpi，"一眼能分开两行"的下限
  budget     <- floor(h * 72 * panel_frac / (fontsize_row + min_gap))
  show_rows  <- nrow(z) <= budget

  # 决策留痕（不落盘则只留日志；"为什么没行名"要可核对）
  message(sprintf('行名标签: %d 行 / 画布可容纳 %d 行（高 %.2fin，字号 %gpt，余量 %gpt）-> %s',
                  nrow(z), budget, h, fontsize_row, min_gap,
                  if (show_rows) '显示行名' else '不显示行名（会糊成一片）'))
  if (!is.null(label_decisions_path)) {
    df <- data.frame(figure = basename(path), label = 'gene',
                     n_labels = nrow(z), capacity = budget,
                     height_in = h, fontsize = fontsize_row,
                     panel_frac = panel_frac, min_gap = min_gap,
                     shown = show_rows)
    write.csv(df, label_decisions_path, row.names = FALSE)
  }

  ord <- hclust(stats::dist(z))$order          # 行按层次聚类排序
  z <- z[ord, , drop = FALSE]
  df <- data.frame(gene = factor(rownames(z), levels = rownames(z)),
                   as.data.frame(z), check.names = FALSE)
  long <- stats::reshape(df, direction = 'long', varying = 2:ncol(df),
                         v.names = 'z', idvar = 'gene',
                         timevar = 'sample', times = colnames(z))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = sample, y = gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = pal_diverging(100),
                                  limits = c(-3, 3), name = 'Z-score') +
    ggplot2::scale_y_discrete(breaks = if (show_rows) rownames(z) else NULL) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_paper() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7))
  save_png(p, path, w = w, h = h)
  invisible(path)
}
