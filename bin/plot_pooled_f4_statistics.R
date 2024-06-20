#!/usr/bin/env Rscript

library(argparse)
parser <- argparse::ArgumentParser()
parser$add_argument("--input", help = "Input file", required = TRUE)
parser$add_argument("--output", help = "Prefix for output PDF files",
    required = TRUE)
parser$add_argument("--groupings",
    help = "CSV file with population groupings for plot", required = TRUE)
parser$add_argument("--order",
    help = "CSV file with population order for plot", required = TRUE)
parser$add_argument("--plot_height", help = "Height of plot in inches",
    type = "numeric", default = 5)
parser$add_argument("--plot_width", help = "Width of plot in inches",
    type = "numeric", default = 10)
args <- parser$parse_args()

if (!file.exists(args$input)) {
    stop(sprintf("Input file does not exist: %s", args$input))
}

if (!file.exists(args$groupings)) {
    stop(sprintf("Groupings file does not exist: %s", args$groupings))
}

if (!file.exists(args$order)) {
    stop(sprintf("Order file does not exist: %s", args$order))
}

library(ggplot2)
library(data.table)

load_f4_data <- function(datafilename, groupingfilename, orderfilename) { # nolint start
    require(data.table)
    stopifnot(all(file.exists(datafilename, groupingfilename, orderfilename)))
    f4s <- fread(datafilename)
    groupings <- fread(groupingfilename)
    ordering <- fread(orderfilename)
    dt <- f4s[groupings, , on = "pop4"]
    dt[ordering, I := i.I, on = "pop4"]
    dt[, ci_l := est - 3 * se]
    dt[, ci_u := est + 3 * se]
    dt
} # nolint end

plot_f4_data <- function(data, target = c("HT", "CTVT")) { # nolint start
    require(ggplot2)
    targetpop <- match.arg(target)
    plot <- ggplot(data[pop2 == targetpop],
                   aes(x = reorder(pop4, I), y = est,
                       # colour = as.factor(Group),
                       ymin = ci_l, ymax = ci_u)) +
        geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey20") +
        geom_point() +
        geom_errorbar(width = 0, linewidth = 1) +
        xlab("Test population") +
        ylab("f4") +
        theme_classic() +
            theme(axis.text.x = element_text(
                angle = 90, hjust = 0.98, vjust = 0.5, size = 6),
              legend.position = "none",
              axis.title.y = element_text(angle = 0, vjust = 0.5),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank())
    plot
} # nolint end

f4 <- load_f4_data(args$input, args$groupings, args$order)
p_ht <- plot_f4_data(f4, target = "HT")
p_ctvt <- plot_f4_data(f4, target = "CTVT")

ggsave(paste(tools::file_path_sans_ext(args$output), "ht.pdf", sep = "_"),
    plot = p_ht, width = args$plot_width, height = args$plot_height)
ggsave(paste(tools::file_path_sans_ext(args$output), "ctvt.pdf", sep = "_"),
    plot = p_ctvt, width = args$plot_width, height = args$plot_height)
