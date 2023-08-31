#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument(
    "--subdir",
    type = "character"
)
parser$add_argument(
    "--outfile",
    type = "character"
)
args <- parser$parse_args()
stopifnot(dir.exists(args$subdir))

annotate_left_right_counts <- function(dt) { # nolint start
    lefts <- dt[, unique(unlist(strsplit(left, '|', fixed = TRUE)))]
    x <- rowSums(dt[, lapply(.SD, function(r) !is.na(r)), .SDcols = lefts])
    dt[, nLeft := x]

    y <- rowSums(dt[, lapply(tstrsplit(right, "|", fixed = TRUE),
        function(r) !is.na(r))])
    dt[, nRight := y]

    return(dt)
} # nolint end

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(data.table)
})
result_files <- list.files(args$subdir, pattern = "*.RDS", full.names = TRUE)
r <- lapply(result_files, readRDS)
suppressWarnings({
  for (i in seq_along(r)) {
    adm_rank <- r[[i]]$rank_estimate
    wave_rank <- r[[i]]$wave_rank_estimate
    wave_ranks_mt <- paste(r[[i]]$multiple_testing_wave_rank_estimates,
        collapse = "|")
    pattern <- paste("p", r[[i]]$pat, sep = "=")
    left <- paste(r[[i]]$left, collapse = "|")
    right <- paste(r[[i]]$right, collapse = "|")
    target <- r[[i]]$target
    r[[i]]$popdrop %<>% mutate(i = i, adm_rank = adm_rank,
        wave_rank = wave_rank, wave_ranks_mt = wave_ranks_mt,
        pattern = pattern, left = left, right = right, target = target)
    r[[i]]$popdrop <- as.data.table(r[[i]]$popdrop)
  }
})

t <- rbindlist(lapply(r, `[[`, "popdrop"), fill = TRUE)
setDT(t)

setorder(t, -feasible, -p)

desired_colorder <- c(
    setdiff(sort(colnames(t)), c("left", "right", "target")),
    c("left", "right", "target")
)

setcolorder(t, desired_colorder)
setcolorder(t, c(
    "i", "adm_rank", "wave_rank", "wave_ranks_mt", "pat", "pattern", "wt",
    "dof", "chisq", "p", "f4rank", "feasible", "best", "dofdiff", "chisqdiff",
    "p_nested"
))

t <- annotate_left_right_counts(t)

fwrite(t, args$outfile,
    sep = "\t", na = "NA", quote = FALSE,
    row.names = FALSE
)
