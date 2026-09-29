# scripts/plot/enrichment.R — ORA 点图（P3）
## 规范：文献三点式（点大小=Count、色=-log10(p)、x=GeneRatio），
## 红=显著（-log10p 渐变），图例框外右侧纵向。
## 输入 ego 为含 ID / Description / GeneRatio / p.adjust / Count 列的数据框。
plot_enrichment <- function(ego, path, title = 'Enrichment', top_n = 15,
                            w = W_ONE_HALF, h = mm(100)) {
  df <- as.data.frame(ego)
  need <- c('Description', 'GeneRatio', 'p.adjust', 'Count')
  missing <- setdiff(need, colnames(df))
  if (length(missing) > 0) stop('plot_enrichment: 缺少必需列: ', paste(missing, collapse = ', '))

  # GeneRatio 形如 "5/200" -> 数值 0.025
  gr <- do.call(rbind, strsplit(as.character(df$GeneRatio), '/'))
  df$gene_ratio <- as.numeric(gr[, 1]) / as.numeric(gr[, 2])
  df$neglog10p  <- -log10(pmax(df$p.adjust, 1e-300))

  df <- df[order(df$p.adjust), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$Description <- factor(df$Description, levels = rev(df$Description))  # 最显著在顶部

  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene_ratio, y = Description,
                                        size = Count, colour = neglog10p)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradientn(colours = pal_sequential(100), name = '-log10(p.adj)') +
    ggplot2::scale_size_continuous(range = c(2, 7), name = 'Count') +
    ggplot2::labs(title = title, x = 'GeneRatio', y = NULL) +
    theme_paper()
  save_png(p, path, w = w, h = h)
  invisible(path)
}
