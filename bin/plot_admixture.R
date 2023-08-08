#!/usr/bin/env Rscript

library(argparse)
parser <- argparse::ArgumentParser()
parser$add_argument("--input_dir", type = "character",
    help = "path to admixture results")
parser$add_argument("--pops", type = "character",
    help = "path to population ids")
parser$add_argument("--output_dir", type = "character",
    help = "path to output directory")
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

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

# load sample ids ---------------------------------------------------------
# Want to make sure samples are in correct order & combine with population ids
# I have population ids in a separate file, prob should just replace the fam file

# Show that the best fit number of clusters averaged over all runs is 5
#-------------------------------------------------------------------------
# Read all data
# cv <- data.table(k = 1:15)
# cv$k2 <- cv$k^2
# cv$cv5 <- fread("./admixture_results/cv5/CV_errors.txt")$V2
# cv$cv10 <- fread("./admixture_results/cv10/CV_errors.txt")$V2
# cv$bcv5 <- fread("./admixture_results/boot200_cv5/CV_errors.txt")$V2
# cv$bcv10 <- fread("./admixture_results/boot200_cv10/CV_errors.txt")$V2
# # Melt the data tall and narrow
# mydata = melt(cv, measure.vars = c("cv5", "cv10", "bcv5", "bcv10"))
# # Fit a quadratic model of how value varies with k^2, k and intercept
# mymodel <- lm(value ~ k + k2 + 1, data = mydata)
# mydata[, plot(value ~ k, col = variable)]
# v <- seq(0, 16, length.out = 200)
# lines(v, predict(mymodel, list(k = v, k2 = v^2)), col = "dodgerblue")
# # Solve for derivative of quadratic being 0: min_x ~= 5
# min_x <- -coef(mymodel)["k"] / (2 * coef(mymodel)["k2"])
# min_y <- predict(mymodel, list(k = min_x, k2 = min_x^2))
# abline(h = y, v = x, col = "grey80", lty = 2)

# # fam file with sample order
# # this is a .fam formatted file from PLINK
# fam <- read.table(file.path(result, "plink.fam")) %>%
#     select(V2) %>%
#     rename(sample = V2)

# # population ids for each sample
# # text file with two columns (sample and population)
groups <- read.table(args$pops) %>%
    rename(sample = V1,
           population = V2) %>% as.data.table()
groups <- as.data.frame(groups[sample != "Sloughi02"])

# load admixture results --------------------------------------------------

# list admixture output files
# assumes files are located in current working directory
afiles <- list.files(path = args$input_dir,
                     pattern = "*.Q$",
                     full.names = TRUE)

afiles <- afiles[!grepl("\\.1\\.Q$", afiles)]

# read in admixture results to a single df
q_df <- map_df(afiles, ~vroom(.x,
                              col_names = FALSE,
                              id = "k",
                              delim = " ")) %>%
    # may need to adjust this - wanted to extract k value from file name
    mutate(k = as.integer(word(k, 2, sep = "[.]"))) %>%
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
q_ordered[, population := factor(population, levels = c("Canids", "EurasianWolves", "Sahul",
                                                        "EastAsia", "European", "SiberianArctic",
                                                        "African", "MESA", "Ancient", "TVT"))]
setDT(q_df)
q_df[, population := factor(population, levels = c("Canids", "EurasianWolves", "Sahul",
                                                        "EastAsia", "European", "SiberianArctic",
                                                        "African", "MESA", "Ancient", "TVT"))]

# plot the results! -------------------------------------------------------
# palette_ <- rev(dutchmasters::dutchmasters$pearl_earring)[1:5]
# palette_[1] <- adjustcolor(palette_[1], blue.f = 0.4, green.f = 0.6, red.f = 2)
# palette_[3] <- adjustcolor(palette_[3], blue.f = 1.1, red.f = 0.9)
for (k_ in 2:15) {
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

    ggsave(file.path(args$output_dir, sprintf("admixture_plot_k%d.pdf", k_)),
        plot = plot,
        height = 4,
        width = 18)
}

palette_ <- palette.colors(15, palette = "Polychrome 36")
# palette_ <- rev(dutchmasters::dutchmasters$pearl_earring)[1:5]
# palette_[1] <- adjustcolor(palette_[1], blue.f = 0.4, green.f = 0.6, red.f = 2)
# palette_[3] <- adjustcolor(palette_[3], blue.f = 1.1, red.f = 0.9)
plot <- q_df %>%
    ggplot(., aes(sample, prop, fill = cluster)) +
    # stacked barplot
    geom_bar(stat = "identity") +
    # use custom colors
    scale_fill_manual(values = unname(palette_)) +
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
plot

# save the output
ggsave(file.path(args$output_dir, "admixture_plot.pdf"),
       plot = plot,
       height = 12,
       width = 18)
