#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description = "Repopulate plink fam file with new population labels")
parser$add_argument("--input-prefix",
    help = "Input prefix of the .bed, .bim and .fam files")
parser$add_argument("--pooling",
    help = "CSV file matching samples to population pools")
parser$add_argument("--output-prefix",
    help = "Output prefix of the .bed, .bim and .fam files")
args <- parser$parse_args()

for (suffix in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(args$input_prefix, suffix))) {
        stop(paste0("File ", args$input_prefix, suffix, " does not exist"))
    }
}
if (!file.exists(args$pooling)) {
    stop(paste0("File ", args$pooling, " does not exist"))
}

if (args$input_prefix == args$output_prefix) {
    stop("Input and output prefixes must be different")
}

library(data.table)
in_set <- fread(args$pooling)

fam <- fread(paste0(args$input_prefix, ".fam"))
setnames(fam, old = c("V1", "V2"), new = c("pop", "sample"))
fam[, pop := as.character(pop)]
fam[, lineorder := .I]
setkey(fam, sample)
fam[in_set, pop := pop1]
fam[, V6 := 1]
fwrite(fam[order(lineorder), .(pop, sample, V3, V4, V5, V6)],
    paste0(args$output_prefix, ".fam"),
    col.names = FALSE, sep = "\t"
)
file.copy(paste0(args$input_prefix, ".bed"),
    paste0(args$output_prefix, ".bed"))
file.copy(paste0(args$input_prefix, ".bim"),
    paste0(args$output_prefix, ".bim"))
