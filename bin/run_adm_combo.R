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

print(sprintf("Target = %s", args$target))
print(sprintf("Left = %s", args$left))
print(sprintf("Right = %s", args$right))

if (!(file.exists(paste0(args$sampleset, ".bed")))) {
    stop("Sampleset file %s does not exist", paste0(args$sampleset, ".bed"))
}

estimate_rank <- function(rank_table, alpha = 0.05) { # nolint start
    if (is.null(rank_table)) {
        return (0)
    }
    dt <- as.data.table(copy(rank_table))
    rank <- dt[, max(f4rank) + 1]
    for (i in seq_len(nrow(dt))) {
        if (dt[i, p] < alpha) {
            break
        }
        rank <- dt[i, f4rank]
        if (dt[i, is.na(p_nested)] | dt[i, p_nested] < alpha) {
            break
        }
    }
    return (rank)
} # nolint end


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

adm_result$rank_estimate <- estimate_rank(adm_result$rankdrop, alpha = 0.05)

# qpadm does not return the f4 matrix when using genotype data files,
# so we need to calculate it ourselves
adm_f4blockdat <- f4blockdat_from_geno(
    args$sampleset,
    left = c(args$target, setdiff(left, args$target)),
    right = right,
    allsnps = TRUE,
    blgsize = 50000,
    auto_only = FALSE
)

adm_result$f4 <- admixtools:::f4blockdat_to_f4out(adm_f4blockdat, boot = FALSE)

if (length(left) > 1) {
    # There are enough left pops to do a qpwave analysis
    wave_result <- qpadm(args$sampleset,
        left = left,
        right = right,
        target = NULL,
        allsnps = TRUE,
        blgsize = 50000,
        auto_only = FALSE,
        verbose = FALSE
    )
    adm_result$wave_rankdrop <- wave_result$rankdrop

    # Store the f4 table for the wave analysis
    wave_f4blockdat <- f4blockdat_from_geno(
        args$sampleset,
        left = left,
        right = right,
        allsnps = TRUE,
        blgsize = 50000,
        auto_only = FALSE,
        verbose = FALSE
    )
    adm_result$wavef4 <- admixtools:::f4blockdat_to_f4out(wave_f4blockdat,
        boot = FALSE)

    # Estimate the rank
    adm_result$wave_rank_estimate <- estimate_rank(wave_result$rankdrop,
        alpha = 0.05)

    # The wave result is affected by the order the left pops are considered.
    # This seems arbitrary, so we will run the analysis for all possible
    # permutations of the left pops and work out the rank estimate for each.
    # Change the significance level to account for multiple tests.
    res <- sapply(seq_along(left),
        function(i) {
        wave_result_i <- qpadm(
            args$sampleset,
            left = c(left[i], setdiff(left, left[i])),
            right = right,
            target = NULL,
            allsnps = TRUE,
            blgsize = 50000,
            auto_only = FALSE,
            verbose = FALSE
        )
        # Estimate rank with and without multiple testing correction
        list(nomt = estimate_rank(wave_result_i$rankdrop,
                alpha = 0.05),
             mt = estimate_rank(wave_result_i$rankdrop,
                    alpha = 0.05 / length(left)))
    })

    adm_result$permutation_wave_rank_estimates <- res["nomt", ]
    adm_result$permutation_wave_rank_estimates_mt <- res["mt", ]
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
