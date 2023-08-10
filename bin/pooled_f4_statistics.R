#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description = "Calculate pooled F4 statistics")
parser$add_argument("--input",
    help = "Prefix of .bed, .bim and .fam input files",
    required = TRUE)
parser$add_argument("--output", help = "Output CSV file",
    required = TRUE)
parser$add_argument("--subset",
    help = "Subset of data to use: all, chr1, chr7 or chr21",
    default = "all",
    choices = c("all", "chr1", "chr7", "chr21"))
args <- parser$parse_args()

# Check input file exists
file_checks <- c(
    file.exists(paste0(args$input, ".bed")),
    file.exists(paste0(args$input, ".bim")),
    file.exists(paste0(args$input, ".fam"))
)
if (!(all(file_checks))) {
    stop("Input files do not exist with this prefix:", args$input)
}

library(admixtools)
library(data.table)


compute_f4s <- function(data_prefix, chrom = NULL, blocksize = 100000) { # nolint start
    snps <- if (is.null(chrom)) {
        bim[, V2]
    } else {
        bim[V1 == chrom, V2]
    }
    f4s <- f4(data_prefix,
        pop1 = "AndeanFox",
        pop2 = c("HT", "CTVT"),
        pop3 = modern,
        pop4 = pops,
        auto_only = FALSE,
        blgsize = blocksize,
        adjust_pseudohaploid = TRUE,
        keepsnps = snps,
        allsnps = TRUE)

    setDT(f4s)
    setorder(f4s, pop2, est)
    f4s
} # nolint end


bim <- fread(paste0(args$input, ".bim"))
fam <- fread(paste0(args$input, ".fam"))
pops <- fam[, 1][[1]]
outgroup <- "AndeanFox"
targets <- c("HT", "CTVT")
modern <- "GermanShepherdDog"
stopifnot(all(c(outgroup, targets, modern) %in% pops))
pops <- setdiff(pops, c(outgroup, modern))
f4s <- if (args$subset == "all") {
    compute_f4s(args$input, blocksize = 50000)
} else {
    compute_f4s(
        args$input,
        chrom = sub("^chr", "", args$subset),
        blocksize = 50000)
}

fwrite(f4s, args$output)
