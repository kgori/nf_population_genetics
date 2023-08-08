#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(
    description = "Run admixture analysis on a single admixture combination"
)

# qpAdm populations - left, right and target
parser$add_argument("--target",
    type = "character", help = "Target population (single value)"
)

parser$add_argument("--left",
    type = "character", help = "Left populations (comma-separated list)"
)

parser$add_argument("--right",
    type = "character", help = "Right populations (comma-separated list)"
)

parser$add_argument("--sampleset",
    type = "character",
    help = "Prefix of the sampleset file",
    default = "plink.pop1"
)

parser$add_argument("--ancient",
    action = "store_true",
    help = "Is the target an ancient sample?"
)

args <- parser$parse_args()

if (!(file.exists(paste0(args$sampleset, ".bed")))) {
    stop("Sampleset file %s does not exist", paste0(args$sampleset, ".bed"))
}

library(data.table)
library(admixtools)

left <- strsplit(args$left, ",")[[1]]
right <- strsplit(args$right, ",")[[1]]

if (args$ancient) {
    right <- setdiff(right, "GermanShepherdDog")
    left <- setdiff(left, "GermanShepherdDog")
}

if (args$target %in% c(left, right)) {
    left <- setdiff(left, args$target)
    right <- setdiff(right, args$target)
    warning("Target population is in the admixture combination")
}

adm_result <- qpadm(args$sampleset,
    left = left,
    right = right,
    target = args$target,
    allsnps = TRUE,
    blgsize = 50000,
    auto_only = FALSE
)

if (length(left) > 1) {
    wave_result <- qpadm(args$sampleset,
        left = left,
        right = right,
        target = NULL,
        allsnps = TRUE,
        blgsize = 50000,
        auto_only = FALSE
    )
    adm_result$wave.rankdrop <- wave_result$rankdrop
}

adm_result$left <- left
adm_result$right <- right
adm_result$target <- args$target

saveRDS(adm_result,
    file = paste0(
        args$target,
        ".",
        paste0(left, collapse = "_"),
        ".",
        paste0(right, collapse = "_"),
        ".RDS")
)
