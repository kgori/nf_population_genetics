#!/usr/bin/env Rscript

library(argparse)
parser <- argparse::ArgumentParser()
parser$add_argument("--input_file", type = "character",
    help = "path to admixture results Q file")
parser$add_argument("--pops", type = "character",
    help = "path to population ids")
parser$add_argument("--output_file", type = "character",
    help = "path to output pdf")
args <- parser$parse_args()

if (!file.exists(args$input_file)) {
    stop(sprintf("Input file does not exist: %s", args$input_file))
}

if (!file.exists(args$pops)) {
    stop(sprintf("Population file does not exist: %s", args$pops))
}

library(dplyr)
library(forcats)
library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)
library(vroom)
library(RColorBrewer)
library(data.table)

# specify color palette (replace number in parentheses with max k value)
cols <- colorRampPalette(brewer.pal(8, "Dark2"))(15)

# # population ids for each sample
# # text file with two columns (sample and population)
groups <- read.table(args$pops) %>%
    rename(sample = V1,
           population = V2) %>% as.data.table()
groups <- as.data.frame(groups[sample != "Sloughi02"])

# load admixture results --------------------------------------------------


# read in admixture results to a single df
q_df <- vroom(args$input_file,
    col_names = FALSE,
    id = "k",
    delim = " ",
    show_col_types = FALSE) %>%
    # may need to adjust this - wanted to extract k value from file name
    mutate(k = as.integer(word(basename(k), 2, sep = "[.]"))) %>%
    # combine results with population ids
    cbind(groups, .) %>%
    # convert to long format
    gather(key = cluster,
           value = prop,
           -k, -sample, -population) %>%
    # remove missing values for lower values of k
    drop_na()


# sort samples by proportion ----------------------------------------------
# https://stackoverflow.com/questions/41679888/how-to-group-order-data-in-r-for-a-barplot

q_ordered <- q_df %>%
    group_by(sample) %>%
    # determine most likely assignment
    mutate(likely_assignment = population[which.max(prop)],
           assignment_prob = max(prop)) %>%
    # sort samples by assignment
    arrange(likely_assignment, desc(assignment_prob)) %>%
    ungroup() %>%
    # make sure samples are ordered
    mutate(sample = fct_inorder(sample))

setDT(q_ordered)
q_ordered[, population := factor(population,
    levels = c("Canids", "EurasianWolves", "Sahul",
        "EastAsia", "European", "SiberianArctic",
        "African", "MESA", "Ancient", "TVT"))]
setDT(q_df)
q_df[, population := factor(population,
    levels = c("Canids", "EurasianWolves", "Sahul",
        "EastAsia", "European", "SiberianArctic",
        "African", "MESA", "Ancient", "TVT"))]

# plot the results! -------------------------------------------------------
for (k_ in q_df[, unique(k)]) {
    palette_choice <- if (k_ > 9) { "Polychrome 36" } else { "Okabe-Ito" }
    palette_ <- palette.colors(k_, pal = palette_choice)
    plot <- q_df %>%
        filter(k == k_) %>%
        ggplot(., aes(sample, prop, fill = cluster)) +
    # stacked barplot
    geom_bar(stat = "identity") +
    # use custom colors
    scale_fill_manual(values = unname(rev(palette_))) +
    # rotate x-axis labels
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    # facet by k value and population id
    facet_grid(k ~ population,
        scales = "free_x",
        space = "free") +
    theme_bw() +
    # rotate strip text for easier reading
    theme(strip.text.x = element_text(angle = 90,
        hjust = 0),
        axis.text.x = element_text(size = 6),
        # reduce spacing between panels
        panel.spacing.x = unit(0.001, "lines"),
        legend.position = "none")

    ggsave(args$output_file,
        plot = plot,
        height = 4,
        width = 18)
}
