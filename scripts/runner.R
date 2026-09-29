#!/usr/bin/env Rscript
## ============================================================
## P3 出图/叙事层：真实绘图实现（volcano / heatmap / enrichment）。
## 分析数值仍为确定性模拟数据（附录 D#1：真实 limma/GEOquery 属后续迭代），
## 但图件全部为真实 ggplot 渲染：语义色板 + mm 图幅 + 图例外置 + 图名五段坐标。
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) { i <- match(paste0('--', name), args); if (is.na(i)) NULL else args[i + 1] }
job_file <- get_arg('job'); outdir <- get_arg('outdir'); logdir <- get_arg('logdir')
if (is.null(job_file) || is.null(outdir)) stop("Usage: runner.R --job <f> --outdir <d> --logdir <d>")

`%||%` <- function(a, b) if (is.null(a)) b else a

suppressPackageStartupMessages({ library(jsonlite) })
job <- fromJSON(job_file, simplifyVector = FALSE)
wo <- job$work_order
seed <- as.integer(wo$params$random_seed %||% 42)

dir.create(file.path(outdir, 'figures'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, 'tables'),   recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

if (file.exists('scripts/plot/utils.R'))      source('scripts/plot/utils.R')
if (file.exists('scripts/plot/volcano.R'))    source('scripts/plot/volcano.R')
if (file.exists('scripts/plot/heatmap.R'))    source('scripts/plot/heatmap.R')
if (file.exists('scripts/plot/enrichment.R')) source('scripts/plot/enrichment.R')

deg_df <- NULL
status <- list()

for (unit in wo$units) {
  uid <- unit$id
  if (isTRUE(unit$enabled == FALSE)) { status[[uid]] <- 'skipped'; next }
  message('Running unit: ', uid)
  status[[uid]] <- 'pending'
  tryCatch({
    if (uid %in% c('CHIP-01', 'RNA-01')) {
      # 数据下载（占位：真实实现接 GEOquery）
      status[[uid]] <- 'ok'
    } else if (uid %in% c('CHIP-04', 'RNA-05')) {
      # 归一化（占位）
      status[[uid]] <- 'ok'
    } else if (uid %in% c('CHIP-06', 'RNA-07')) {
      # 差异分析（占位：确定性模拟数据，供全链路验证绘图与交接）
      set.seed(seed)
      n <- 500
      lfc <- c(rnorm(40, 2.2, 0.4), rnorm(40, -2.1, 0.4), rnorm(n - 80, 0, 0.3))
      pval <- c(rbeta(80, 1, 30), rbeta(n - 80, 1, 1))
      deg_df <- data.frame(
        gene      = paste0('GENE', sprintf('%04d', seq_len(n))),
        logFC     = round(lfc, 4),
        P.Value   = round(pval, 6),
        adj.P.Val = round(p.adjust(pval, 'BH'), 6),
        stringsAsFactors = FALSE
      )
      write.csv(deg_df, file.path(outdir, 'tables', 'degs.csv'), row.names = FALSE)

      # P3 附带产物：top 热图矩阵 + mock 富集表（供 heatmap / enrichment 出图）
      top_genes <- deg_df[order(deg_df$P.Value), ][seq_len(min(40, nrow(deg_df))), ]
      expr_mat <- matrix(NA_real_, nrow = nrow(top_genes), ncol = 6)
      colnames(expr_mat) <- c(paste0('T', 1:3), paste0('N', 1:3))
      rownames(expr_mat) <- top_genes$gene
      set.seed(seed + 1)
      eff <- ifelse(top_genes$logFC > 0, 1, -1)
      for (g in seq_len(nrow(expr_mat))) {
        expr_mat[g, 1:3] <- rnorm(3, 2 * eff[g], 0.4)
        expr_mat[g, 4:6] <- rnorm(3, -2 * eff[g], 0.4)
      }
      assign('.top_mat', expr_mat, envir = globalenv())

      set.seed(seed + 2)
      go_terms <- c('mitotic cell cycle', 'cell adhesion', 'signal transduction',
                    'apoptotic process', 'immune response', 'metabolic process',
                    'DNA repair', 'angiogenesis', 'cell proliferation',
                    'inflammatory response', 'chromatin remodeling', 'EMT')
      enrich_df <- data.frame(
        ID         = paste0('GO:', sprintf('%07d', seq_along(go_terms))),
        Description = go_terms,
        GeneRatio  = sprintf('%d/%d', sample(3:15, length(go_terms)), 500),
        p.adjust   = round(p.adjust(rbeta(length(go_terms), 1, 12), 'BH'), 5),
        Count      = sample(3:15, length(go_terms)),
        stringsAsFactors = FALSE)
      write.csv(enrich_df, file.path(outdir, 'tables', 'enrichment.csv'), row.names = FALSE)
      assign('.enrich_df', enrich_df, envir = globalenv())
      status[[uid]] <- 'ok'
    } else if (uid %in% c('CHIP-09', 'RNA-10')) {
      # 筛选报告（占位）
      status[[uid]] <- 'ok'
    } else {
      message('  (placeholder, no-op)')
      status[[uid]] <- 'ok'
    }
  }, error = function(e) {
    status[[uid]] <<- paste0('failed: ', conditionMessage(e))
  })
}

# ---- P3 出图（真实渲染）：图名五段 <阶段>-<模块>-<图>-unit<N>-<slug> ----------
if (!is.null(deg_df) && exists('plot_volcano')) {
  figdir <- file.path(outdir, 'figures')
  fig_names <- character(0)
  tryCatch({
    p1 <- fig_filename('01', 'deg', 'volcano', 'unit1', 'plot')
    plot_volcano(deg_df, file.path(figdir, p1),
                 title = 'Volcano (simulated data)')
    fig_names <- c(fig_names, p1)
    status[['plot_volcano']] <- 'ok'
  }, error = function(e) status[['plot_volcano']] <<- paste0('failed: ', conditionMessage(e)))
  if (exists('.top_mat') && exists('plot_heatmap')) {
    tryCatch({
      p2 <- fig_filename('01', 'deg', 'heatmap', 'unit1', 'top40')
      plot_heatmap(.top_mat, file.path(figdir, p2), title = 'Top40 DEG Heatmap',
                   label_decisions_path = file.path(outdir, 'label_decisions.csv'))
      fig_names <- c(fig_names, p2)
      status[['plot_heatmap']] <- 'ok'
    }, error = function(e) status[['plot_heatmap']] <<- paste0('failed: ', conditionMessage(e)))
  }
  if (exists('.enrich_df') && exists('plot_enrichment')) {
    tryCatch({
      p3 <- fig_filename('01', 'enrich', 'dotplot', 'unit1', 'go')
      plot_enrichment(.enrich_df, file.path(figdir, p3), title = 'GO Enrichment')
      fig_names <- c(fig_names, p3)
      status[['plot_enrichment']] <- 'ok'
    }, error = function(e) status[['plot_enrichment']] <<- paste0('failed: ', conditionMessage(e)))
  }
  # 图清单落盘（P3/P4 交接用）
  write.csv(data.frame(figure_id = sprintf('fig-%03d', seq_along(fig_names)),
                       filename = fig_names), file.path(figdir, 'fig_names.csv'),
            row.names = FALSE)
}

export_meta <- list(
  format = 'export_matrix_v1',
  matrix_spec = list(
    scale_type = 'placeholder', is_log_transformed = FALSE,
    batch_corrected = FALSE, batch_correction_scope = 'none',
    gene_id_type = 'SYMBOL_GENCODE_v44', reference_genome = 'GRCh38',
    n_genes = if (is.null(deg_df)) 0 else nrow(deg_df), n_samples = 1
  ),
  metadata_spec = list(
    group_variable = 'group', reference_level = 'normal',
    sample_id_column = 'sample_id',
    sample_filtering = list(removed_samples = list(), reason = NULL, removed_at_unit = NULL)
  ),
  provenance = list(
    source_domain = wo$domain_id, source_unit = 'CHIP-06',
    source_run_id = job$run_id %||% 'unknown',
    unit_version = 'v0.1-placeholder',
    container_digest = Sys.getenv('CONTAINER_DIGEST', 'unknown'),   # P-09：由 workflow 注入真实值
    created_at = format(Sys.time(), tz = 'UTC', format = '%Y-%m-%dT%H:%M:%SZ')
  )
)
write(toJSON(export_meta, auto_unbox = TRUE, pretty = TRUE),
      file.path(outdir, 'export_metadata.json'))

# ==== 补丁 12（B5）：runner_status.json 主落 outdir（进 artifact，P4 可读），logdir 镜像 ====
status_payload <- toJSON(list(run_id = job$run_id %||% 'unknown', status = status),
                         auto_unbox = TRUE, pretty = TRUE)
write(status_payload, file.path(outdir, 'runner_status.json'))
write(status_payload, file.path(logdir,  'runner_status.json'))

# ==== 补丁 13（B1 次生）：显式退出语义 ====
failed <- Filter(function(s) grepl('^failed', s), status)
if (length(failed) > 0) {
  for (n in names(failed)) message('FAILED unit ', n, ': ', failed[[n]])
  quit(save = 'no', status = 1)
}
message('runner.R finished OK')
