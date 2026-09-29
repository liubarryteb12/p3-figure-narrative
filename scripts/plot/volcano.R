# scripts/plot/volcano.R — 火山图（P3）
## 规范：红=Up / 蓝=Down / 灰=ns（语义色板全篇一致）；
## 阈值线虚线标出；图例框外右侧纵向；mm 图幅 183mm 封顶；PNG 300dpi + PDF 矢量。
plot_volcano <- function(deg, path, title = 'Volcano', pcol = 'adj.P.Val', lfc_col = 'logFC',
                         p_cut = 0.05, lfc_cut = 1.0,
                         w = W_DOUBLE, h = mm(90)) {
  deg <- as.data.frame(deg)
  deg$sig <- ifelse(!is.na(deg[[pcol]]) & deg[[pcol]] < p_cut & abs(deg[[lfc_col]]) > lfc_cut,
                    ifelse(deg[[lfc_col]] > 0, 'Up', 'Down'), 'NS')
  deg$neglog10p <- -log10(pmax(deg[[pcol]], 1e-300))
  n_up <- sum(deg$sig == 'Up'); n_dn <- sum(deg$sig == 'Down'); n_ns <- sum(deg$sig == 'NS')
  cols <- c(Up = PAL$up, Down = PAL$down, NS = PAL$ns)
  subtitle <- sprintf('Up %d / Down %d / NS %d  (adj.P < %g, |log2FC| > %g)',
                      n_up, n_dn, n_ns, p_cut, lfc_cut)
  p <- ggplot2::ggplot(deg, ggplot2::aes(x = .data[[lfc_col]], y = neglog10p, color = sig)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.2) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed',
                        colour = PAL$muted, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = 'dashed',
                        colour = PAL$muted, linewidth = 0.3) +
    ggplot2::scale_color_manual(values = cols,
                                breaks = c('Up', 'Down', 'NS')) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = 'log2 FC', y = '-log10 adj.P', color = '') +
    theme_paper()
  save_png(p, path, w = w, h = h)
  invisible(path)
}
