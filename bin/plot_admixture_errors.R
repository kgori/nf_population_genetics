#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description = "Plot admixture errors")
parser$add_argument("-i", "--input",
    help = "Input file", required = TRUE)
parser$add_argument("-o", "--output",
    help = "Output file", required = TRUE)
args <- parser$parse_args()

if (!file.exists(args$input)) {
    stop("Input file does not exist")
}

library(data.table)

dt <- fread(args$input)
setorder(dt, rep, k)

colours <- colorspace::lighten(palette.colors(10, pal = "Set3"))

pdf(args$output, width = 10, height = 8)
dt[, plot(k, error, type = "n")]
for (rep_id in 1:10) {
    dt[rep == rep_id, lines(k, error, type = "l", col = colours[rep])]
}
dt[, .(mean_err = mean(error)), by = k][,
  lines(k, mean_err, type = "o", pch = 20,
    lwd = 2, cex = 1.2, col = "grey10")]
dev.off()
