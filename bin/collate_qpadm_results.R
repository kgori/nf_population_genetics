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

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(data.table)
})
result_files <- list.files(args$subdir, pattern = "*.RDS", full.names = TRUE)
r <- lapply(result_files, readRDS)
suppressWarnings({
  for (i in seq_along(r)) {
    waverank <- r[[i]]$wave.rankdrop
    rank <- r[[i]]$rankdrop
    max_wr <- max(waverank[waverank$p_nested < 0.05, ]$f4rank, na.rm = TRUE)
    max_r <- max(rank[rank$p_nested < 0.05, ]$f4rank, na.rm = TRUE)

    left <- paste(r[[i]]$left, collapse = "|")
    right <- paste(r[[i]]$right, collapse = "|")
    target <- r[[i]]$target
    r[[i]]$popdrop %<>% mutate(i = i, max_wr = max_wr, max_r = max_r,
        left = left, right = right, target = target)
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
    "i", "max_wr", "max_r", "pat", "wt", "dof", "chisq", "p", "f4rank",
    "feasible", "best", "dofdiff", "chisqdiff", "p_nested"
))

fwrite(t, args$outfile,
    sep = "\t", na = "NA", quote = FALSE,
    row.names = FALSE
)
