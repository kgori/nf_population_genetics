#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--input", help = "Input file")
parser$add_argument("--filter",
    help = "Filter file, containing names of samples to filter out")
parser$add_argument("--second-filter",
    help = "Filter file, containing names of samples to omit from smaller plots")
parser$add_argument("--populations-file",
    help = "File containing population information")
parser$add_argument("--output", help = "Output prefix")
parser$add_argument("--plot-width", type = "integer",
    help = "Width of plot in inches")
default_args <- list(input = "../final_merged_biallelic_sites_filtered_annotated_transversions",
    filter = "./excluded_samples_from_f4_plot.txt",
    second_filter = "./further_excluded_samples_from_f4_plot.txt",
    output = "./Aug4",
    plot_width = 18)
args <- parser$parse_args()
if (all(sapply(args, is.null))) {
    args <- default_args
}
print(args)

data_prefix <- tools::file_path_sans_ext(args$input)

library(admixtools)
library(data.table)
library(ggplot2)
library(patchwork)

fam <- fread(paste0(data_prefix, ".fam"))
bim <- fread(paste0(data_prefix, ".bim"))
pops <- fam[, 1][[1]]
outgroup <- "AndeanFox01"
targets <- c("HT", "CTVT")
modern <- "GermanShepherd01"
stopifnot(all(c(outgroup, targets, modern) %in% pops))
pops <- setdiff(pops, c(outgroup, modern))
if (!is.null(args$filter)) {
    if (file.exists(args$filter)) {
        filter <- fread(args$filter, header = FALSE)
        pops <- setdiff(pops, filter[, V1])
    }
}


compute_f4s <- function(chrom = NULL, blocksize = 100000) { # nolint start
    snps <- if (is.null(chrom)) {
        bim[, V2]
    } else {
        bim[V1 == chrom, V2]
    }
    f4s <- f4(data_prefix,
        pop1 = outgroup,
        pop2 = targets,
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

make_plot <- function(table, pointsize = 2.0, linewidth = 1.0, title = "f4(AndeanFox, N-HT1, German Shepherd, 'sample')") { # nolint start   table <- as.data.table(table)
    require(ggplot2)
    require(data.table)
    l <- tables()
    table[
        l$metadata,
        population := i.population,
        on = c("pop4" = "sample")
    ]

    table[, sample := pop4]


    table[, ci_l := est - 3 * se]
    table[, ci_u := est + 3 * se]
    table[, zcat := abs(z) < 3]
    setorder(table, est)
    table[, I := .I]


    p <- ggplot(
        data = table,
        mapping = aes(
            x = reorder(sample, -I),
            y = est,
            ymin = ci_l,
            ymax = ci_u,
            colour = population
        )
    )
    p +
        geom_hline(yintercept = 0, colour = "grey80") +
        geom_errorbar(width = 0, linewidth = linewidth) +
        geom_point(size = pointsize) +
        theme_bw() +
        theme(axis.text.y = element_text(
            size = 5, angle = 0,
            hjust = 0.98, vjust = 0.5
        ),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(
            size = 5, angle = 90,
            hjust = 0.98, vjust = 0.5
        )
        ) +
        scale_colour_manual(values = l$colours) +
        xlab("Sample") + ylab("f4 statistic") + ggtitle(title)
} # nolint stop

tables <- function() {
    require(data.table)
    metadata <- fread(args$populations_file)
    metadata[ancient == TRUE, population := "Ancient"]
    metadata[startsWith(pop2, "Dog") & pop2 == "DogEurope", population := "Europe"]
    metadata[startsWith(pop2, "Dog") & pop2 %in% c("DogMix", "DogUnknown", "DogIntermedWE"), population := "Mix/Unknown"]
    metadata[startsWith(pop2, "Dog") & pop2 == "DogOceania", population := "Oceania"]
    metadata[startsWith(pop2, "Dog") & pop2 == "DogSiberiaAmerica", population := "Siberia/America"]
    metadata[startsWith(pop2, "Dog") & pop2 %like% "Africa", population := "Africa"]
    metadata[startsWith(pop2, "Dog") & pop2 %like% "NorthAmerica", population := "North America"]
    metadata[startsWith(pop2, "Dog") & pop2 %like% "SouthAmerica", population := "South America"]
    metadata[startsWith(pop2, "Dog") & pop2 %like% "M.EastS.Asia", population := "Middle East/South Asia"]
    metadata[startsWith(pop2, "Dog") & pop2 %like% "EastEurasia", population := "East Asia"]
    metadata[startsWith(pop2, "Coyote"), population := "Outgroup"]
    metadata[startsWith(pop2, "Wolf"), population := "Outgroup"]
    metadata[endsWith(pop2, "Wolf"), population := "Outgroup"]
    metadata[startsWith(pop2, "AndeanFox"), population := "Outgroup"]
    metadata[startsWith(pop2, "Dhole"), population := "Outgroup"]
    metadata[startsWith(pop2, "AfricanHuntingDog"), population := "Outgroup"]
    metadata[startsWith(pop2, "GoldenJackal"), population := "Outgroup"]
    metadata[sample == "HT", population := "N-HT1"]
    metadata[sample == "CTVT", population := "CTVT"]
    metadata[sample %in% c("ASHQ01", "ASHQ06", "ASHQ08", "THRZ02", "TGEZ06", "UZAA01"), population := "AncientIsrael"]
    metadata[sample %in% c("OL4223", "AL3194", "AL3223", "CGG6", "C26", "C27"), population := "AncientArctic"]

    colors = c("Siberia/America" = rgb(150, 210, 235, maxColorValue = 255),
        "Ancient" = rgb(152, 152, 152, maxColorValue = 255),
        "Europe" =rgb(0, 121, 177, maxColorValue = 255),
        "Oceania" = rgb(178, 132, 185, maxColorValue = 255),
        "Africa" = rgb(255, 166, 100, maxColorValue = 255),
        "Middle East/South Asia" = rgb(249, 122, 127, maxColorValue = 255),
        "East Asia" = rgb(178, 216, 143, maxColorValue = 255),
        "N-HT1" = "gold",
        "CTVT" = "red",
        "Outgroup" = "grey20",
        "Mix/Unknown" = "seagreen",
        "AncientIsrael" = rgb(40, 200, 40, maxColorValue = 255),
        "AncientArctic" = rgb(40, 40, 200, maxColorValue = 255))
    for (pop in names(colors)) {
        metadata[population == pop, colour := colors[pop]]
    }

    groups <- lapply(sort(metadata[, unique(population)]), function(popname) {
        metadata[population == popname, sample]
    })
    names(groups) <- sort(metadata[, unique(population)])
    colours <- unique(metadata[, .(colour, population)])[, setNames(colour, population)]
    groups["N-HT1"] <- "HT"
    colours["N-HT1"] <- "gold"
    groups[["Coyote"]] <- c("AndeanFox01", "Coyote01")
    names(groups)[names(groups) == "Coyote"] <- "Outgroup"
    names(colours)[names(colours) == "Coyote"] <- "Outgroup"
    list(metadata = metadata, colours = colours)
}

f4s_all <- compute_f4s(blocksize = 50000)
f4s_chr1 <- compute_f4s(1, blocksize = 50000)
f4s_chr7 <- compute_f4s(7, blocksize = 50000)
f4s_chr21 <- compute_f4s(21, blocksize = 50000)

top_plot_theme <- theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 24),
          panel.border = element_rect(linewidth = 2),
          title = element_text(size = 28),
          legend.position = c(0.6, 0.92),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          legend.background = element_blank(),
          legend.direction = "horizontal")

bottom_left_plot_theme <- theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 28),
          axis.text.y = element_text(size = 24),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(linewidth = 2),
          title = element_text(size = 24),
          legend.position = "none")

bottom_plot_theme <- theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          #axis.text.y = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(linewidth = 2),
          title = element_text(size = 24),
          legend.position = "none")

plots <- list()
pt_top <- 4.0
lw_top <- 2.4
pt_btm <- 3.6
lw_btm <- 2.0

hide_from_plot <- NULL
if (!is.null(args$second_filter)) {
    if (file.exists(args$second_filter)) {
        hide_from_plot <- fread(args$second_filter, header = FALSE)[, V1]
    }
}

for (target in c("CTVT", "HT")) {
    for (i in 1:4) {
        f4s_table <- list(f4s_all, f4s_chr1, f4s_chr7, f4s_chr21)[[i]]
        if (i == 1) {
            title_text <- paste0("f4(Andean Fox, ", names(f4s_table)[[1]], ", GermanShepherd, *) All chromosomes")
            p <- make_plot(f4s_table[pop2 == target & pop4 != target], title = title_text, pointsize = pt_top, linewidth = lw_top) + top_plot_theme
            plots[["All"]] <- p
        } else if (i == 2) {
            title_text <- paste0("Chr 1")
            p <- make_plot(f4s_table[pop2 == target & pop4 != target & !(pop4 %in% hide_from_plot)], title = title_text, pointsize = pt_btm, linewidth = lw_btm) + bottom_plot_theme
            plots[["Chr1"]] <- p
        } else if (i == 3) {
            title_text <- paste0("Chr 7")
            p <- make_plot(f4s_table[pop2 == target & pop4 != target & !(pop4 %in% hide_from_plot)], title = title_text, pointsize = pt_btm, linewidth = lw_btm) + bottom_plot_theme
            plots[["Chr7"]] <- p
        } else if (i == 4) {
            title_text <- paste0("Chr 21")
            p <- make_plot(f4s_table[pop2 == target & pop4 != target & !(pop4 %in% hide_from_plot)], title = title_text, pointsize = pt_btm, linewidth = lw_btm) + bottom_left_plot_theme
            plots[["Chr21"]] <- p
        }
    }
    ggsave(paste0(args$output, "_", target, "_combined.pdf"),
    plots[["All"]] / (plots[["Chr21"]] + plots[["Chr1"]] + plots[["Chr7"]]), width = args$plot_width, height = 10)
}
