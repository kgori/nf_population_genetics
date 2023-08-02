#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description = "Run BIONJ on a distance matrix")
parser$add_argument("-d", "--distances", required = TRUE, help = "Distance matrix")
parser$add_argument("-i", "--ids", required = TRUE, help = "IDs for the distance matrix")
parser$add_argument("-p", "--prefix", required = TRUE, help = "Prefix for output files")
args <- parser$parse_args()

stopifnot(file.exists(args$distances))
stopifnot(file.exists(args$ids))

library(ape)

m <- as.matrix(
    read.delim(args$distances, header = FALSE)
)

n <- read.delim(args$ids, header = FALSE)[, 2]

colnames(m) <- rownames(m) <- n

d <- as.dist(m)

t <- bionj(d)

if ("AndeanFox01" %in% t$tip.label) t <- root(t, "AndeanFox01")

t <- ladderize(t)

pdf(paste0(args$prefix, "_bionj.pdf"), width = 18, height = 128)
plot(t)
dev.off()

write.tree(t, paste0(args$prefix, "_bionj.nwk"))
