#!/usr/bin/env Rscript

library(argparse)
parser <- argparse::ArgumentParser()
parser$add_argument("--input", help = "Input file")
parser$add_argument("--output", help = "Prefix for output PDF files")
parser$add_argument("--groupings",
    help = "CSV file with population groupings for plot")
args <- parser$parse_args()

if (!file.exists(args$input)) {
    stop(sprintf("Input file does not exist: %s", args$input))
}

if (!file.exists(args$groupings)) {
    stop(sprintf("Groupings file does not exist: %s", args$groupings))
}

library(data.table)

load_f4_data <- function(datafilename, groupingfilename) { # nolint start
    require(data.table)
    stopifnot(all(file.exists(datafilename, groupingfilename)))
    f4s <- fread(datafilename)
    groupings <- fread(groupingfilename)
    dt <- f4s[groupings, , on = "pop4"]
    setorder(dt, pop2, Group, -est)
    dt[, I := seq_len(.N)]
    dt[, ci_l := est - 3 * se]
    dt[, ci_u := est + 3 * se]
    dt
} # nolint end

# Reorder sample data to be in descending order of f4 within population groups,
# as measured for HT.
reorder_table <- function(dt) { # nolint start
    setorder(dt, pop2, Group, -est)
    dt[, I := seq_len(.N)]
    truth <- dt[pop2 == "HT", .(pop4, I)]
    dt[, I := NULL]
    dt[truth, I := i.I, on = "pop4"]
    dt
} # nolint end

f4 <- reorder_table(load_f4_data(args$input, args$groupings))
fwrite(f4[, .(pop4, I)], file.path(args$output))
