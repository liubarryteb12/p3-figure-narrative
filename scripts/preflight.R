#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: preflight.R <job.json>")
suppressPackageStartupMessages({ library(jsonlite) })
job <- fromJSON(args[1], simplifyVector = FALSE)
wo <- job$work_order
domain <- wo$domain_id
if (!domain %in% c("geo_chip", "geo_rna", "ml")) stop("unknown domain: ", domain)
if (domain %in% c("geo_chip", "geo_rna")) {
  has_entry <- any(sapply(wo$units, function(u) grepl("^CHIP-01$|^RNA-01$", u$id)))
  if (!has_entry) stop("geo domain requires CHIP-01 or RNA-01")
}
cat("Pre-flight passed for domain:", domain, "\n")
