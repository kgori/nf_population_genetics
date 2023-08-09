#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description = "Annotate a tree with bootstrap values")
parser$add_argument("tree", help = "Newick tree file")
parser$add_argument("bootstrap", help = "Directory containing bootstrap trees")

args <- parser$parse_args()

stopifnot(file.exists(args$tree))
stopifnot(dir.exists(args$bootstrap))

library(ape)

files <- list.files(args$bootstrap, pattern = "nwk$", full.names = TRUE)
stopifnot(length(files) > 0)
bootstrap_trees <- lapply(files, read.tree)

tree <- read.tree(args$tree)

part <- prop.part(bootstrap_trees)
partition_counts <- prop.clades(tree, part = part)
partition_counts[is.na(partition_counts)] <- 0

bootstrap_pct <- 100 * partition_counts / length(bootstrap_trees)

tree$node.labels <- bootstrap_pct

write.tree(tree, file = paste(tools::file_path_sans_ext(args$tree),
                              "annotated.nwk", sep = "."))
