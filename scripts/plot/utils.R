# scripts/plot/utils.R — P3 出图层共享工具
## 语义色板（红=高/实验、蓝=低/对照，全篇一致）+ 毫米制画幅 + 出版主题
## 依据：老框架 geo-normal-pipeline-skill/scripts/lib/common.R（已验证实现）
##       + governance/08_SCI_FIGURE_FEATURES.md（画法特征库）

## ---- 语义色板 --------------------------------------------------------------
## up=红(#B2182B ColorBrewer RdBu 红端) / down=蓝(#2166AC RdBu 蓝端)，
## 与发散色尺同源：热图上的红就是火山图上"上调"的红，读者不用重学。
## ns 用中性灰（不承载色相语义）；primary 与 Python 侧同名同值。
PAL <- list(
  up       = "#B2182B",
  down     = "#2166AC",
  ns       = "#BDBDBD",
  primary  = "#0072B2",
  orange   = "#E69F00",
  mid      = "#F7F7F7",
  ink      = "#1A1A1A",
  muted    = "#666666",
  grid     = "#E8E8E8"
)

## 发散色板：RdBu（蓝→近白→红），热图 Z-score 用
pal_diverging <- function(n = 100) {
  grDevices::colorRampPalette(c(PAL$down, PAL$mid, PAL$up))(n)
}

## 序列色板：viridis 截去暗端（#365C8D 距 down 蓝仅 3.1° 色相，会读成"下调"）
pal_sequential <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#277F8E", "#1FA187", "#4AC16D", "#A0DA39", "#FDE725"))(n)
}

## 分类色板：色盲安全四色（穷举验证：最小两两 OKLab 距离最大）
## 超过 4 类时插值降级（插值色未经验证，图注需说明）
pal_categorical <- function(k) {
  base <- c("#2E6B3E", "#B5527A", "#56B4E9", "#E69F00")  # 绿/玫红/天蓝/橙
  if (k <= length(base)) return(base[seq_len(k)])
  grDevices::colorRampPalette(base)(k)
}

## 组间对比配色：两端红/蓝固定（tumor/上调=红、normal/对照=蓝），由对比语义
## 决定而非分组出现顺序（修"PCA 蓝、火山红"不一致）
pal_condition <- function(arms) {
  arms <- as.character(arms)
  if (length(arms) == 2L) return(stats::setNames(c(PAL$up, PAL$down), arms))
  stats::setNames(c(PAL$up, PAL$down,
                    pal_categorical(max(0L, length(arms) - 2L)))[seq_along(arms)], arms)
}

## ---- 画幅（毫米制）----------------------------------------------------------
## 期刊栏宽按毫米规定；宽度一律夹到标准三档，不随类别数无限增长
W_SINGLE   <- 89   / 25.4   # 单栏
W_ONE_HALF <- 136  / 25.4   # 一栏半
W_DOUBLE   <- 183  / 25.4   # 双栏（通栏）—— 上限
mm <- function(...) c(...) / 25.4

## ---- 主题：图例一律框外右侧、纵向 -------------------------------------------
## legend.direction 管"单个图例内部键竖排"；legend.box 管"多个图例之间竖摞"，
## 两个都要设纵向——并排后总宽超出 183mm 会被 ggplot 静默裁掉且不报警。
theme_paper <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(colour = PAL$grid, linewidth = 0.25),
      panel.border      = ggplot2::element_rect(colour = PAL$grid, linewidth = 0.4),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(colour = PAL$muted, size = base_size - 1.5),
      plot.margin       = ggplot2::margin(5, 6, 4, 5),
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.box        = "vertical",
      legend.title      = ggplot2::element_text(size = base_size - 1),
      legend.text       = ggplot2::element_text(size = base_size - 1.5),
      legend.key.size   = ggplot2::unit(0.7, "lines"),
      legend.margin     = ggplot2::margin(1, 1, 1, 1),
      legend.box.spacing = ggplot2::unit(3, "pt"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

## ---- 出图：PNG（300dpi 起步，投稿下限）+ 同名 PDF 矢量 ------------------------
## 宽度默认双栏 183mm；超过 W_DOUBLE 夹到 W_DOUBLE（不画装不进一页的图）。
## PNG 走 cairo 抗锯齿；失败不中断（PDF 已拿到）。
save_png <- function(p, path, w = W_DOUBLE, h = mm(90), dpi = 300) {
  w <- min(w, W_DOUBLE)
  pdf_path <- sub("\\.png$", ".pdf", path)
  tryCatch({
    grDevices::pdf(pdf_path, width = w, height = h)
    print(p)
    grDevices::dev.off()
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    message('WARN: PDF 输出失败（继续出 PNG）: ', conditionMessage(e))
  })
  tryCatch({
    grDevices::png(path, width = w, height = h, units = "in", res = dpi, type = "cairo")
    print(p)
    grDevices::dev.off()
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    stop('PNG 输出失败: ', conditionMessage(e))
  })
  invisible(path)
}

## ---- 图名五段坐标 <阶段>-<模块>-<图>-unit<N>-<slug> ---------------------------
## 阶段号本模板统一 01（geo 约定）；fig_names.csv 是 P3/P4 交接的图清单。
fig_filename <- function(phase, module, kind, unit, slug, ext = "png") {
  sprintf("%s-%s-%s-%s-%s.%s", phase, module, kind, unit, slug, ext)
}

## ---- 统计三件套进图（HR(95%CI)+p / AUC / 显著性星）的 p 显示规则 -------------
## p < 0.001 显示 "< 0.001"，不显示 "0.000"（文献硬惯例）
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  if (p < 0.01)  return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

## 显著性星号：*<0.05 **<0.01 ***<0.001，不显著显式标 "ns"（诚实惯例）
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}
