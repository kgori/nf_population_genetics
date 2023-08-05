#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description = "Plot tree with populations")
parser$add_argument("-p", "--populations", help = "Populations file", required = TRUE)
parser$add_argument("-d", "--pcadata", help = "PCA file", required = TRUE)
parser$add_argument("-i", "--ingroup", help = "Ingroup", required = TRUE)
parser$add_argument("-o", "--output", help = "Output file", required = TRUE)
args <- parser$parse_args()

stopifnot(file.exists(args$populations))
stopifnot(file.exists(args$ingroup))
if (!file.exists(args$pcadata)) {
    stop("PCA data file does not exist - ", args$pcadata)
}

suppressPackageStartupMessages({
    library(plotly)
    library(data.table)
    library(ggplot2)
})

metadata <- fread(args$populations)
metadata[ancient == TRUE, population := "Ancient"]
metadata[startsWith(pop2, "Dog") & pop2 == "DogEurope",
    population := "Europe"]
metadata[startsWith(pop2, "Dog") & pop2 %in% c("DogMix", "DogUnknown",
                                               "DogIntermedWE"),
    population := "Mix/Unknown"]
metadata[startsWith(pop2, "Dog") & pop2 == "DogOceania",
    population := "Oceania"]
metadata[startsWith(pop2, "Dog") & pop2 == "DogSiberiaAmerica",
    population := "Siberia/America"]
metadata[startsWith(pop2, "Dog") & pop2 %like% "Africa",
    population := "Africa"]
metadata[startsWith(pop2, "Dog") & pop2 %like% "NorthAmerica",
    population := "North America"]
metadata[startsWith(pop2, "Dog") & pop2 %like% "SouthAmerica",
    population := "South America"]
metadata[startsWith(pop2, "Dog") & pop2 %like% "M.EastS.Asia",
    population := "Middle East/South Asia"]
metadata[startsWith(pop2, "Dog") & pop2 %like% "EastEurasia",
    population := "East Asia"]
metadata[startsWith(pop2, "Coyote"), population := "Outgroup"]
metadata[startsWith(pop2, "Wolf"), population := "Outgroup"]
metadata[endsWith(pop2, "Wolf"), population := "Outgroup"]
metadata[startsWith(pop2, "AndeanFox"), population := "Outgroup"]
metadata[startsWith(pop2, "Dhole"), population := "Outgroup"]
metadata[startsWith(pop2, "AfricanHuntingDog"), population := "Outgroup"]
metadata[startsWith(pop2, "GoldenJackal"), population := "Outgroup"]
metadata[sample == "HT", population := "N-HT1"]
metadata[sample == "CTVT", population := "CTVT"]
metadata[sample %in% c("ASHQ01", "ASHQ06", "ASHQ08",
                       "THRZ02", "TGEZ06", "UZAA01"),
    population := "AncientIsrael"]
metadata[sample %in% c("OL4223", "AL3194", "AL3223", "CGG6", "C26", "C27"),
    population := "AncientArctic"]

colors <- c("Siberia/America" = rgb(150, 210, 235, maxColorValue = 255),
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
colours <- unique(metadata[, .(colour, population)])[,
    setNames(colour, population)]
groups["N-HT1"] <- "HT"
colours["N-HT1"] <- "gold"
groups[["Coyote"]] <- c("AndeanFox01", "Coyote01")
names(groups)[names(groups) == "Coyote"] <- "Outgroup"
names(colours)[names(colours) == "Coyote"] <- "Outgroup"

pca <- fread(cmd = paste("grep -v '#'",  args$pcadata))[, 1:11]

setnames(pca, c("sample", paste0("PC", seq_len(ncol(pca) - 1))))
ingroup <- fread(args$ingroup)
pca <- pca[ingroup, , on = "sample"]
pca[metadata, colour := i.colour, on = "sample"]
pca[metadata, population := i.population, on = "sample"]
pca <- pca[population != "Unknown"]

library(patchwork)

pc12 <- ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = population)) +
    scale_colour_manual(values = colours) +
    theme_bw()

pc13 <- ggplot(pca, aes(x = PC1, y = PC3)) +
    geom_point(aes(colour = population)) +
    scale_colour_manual(values = colours) +
    theme_bw()

pc23 <- ggplot(pca, aes(x = (-PC2), y = PC3)) +
    geom_point(aes(colour = population)) +
    scale_colour_manual(values = colours) +
    theme_bw()

pc32 <- ggplot(pca, aes(x = PC3, y = PC2)) +
    geom_point(aes(colour = population)) +
    scale_colour_manual(values = colours) +
    theme_bw()

plot <- pc12 | ((pc12 + theme(legend.position = "none") | pc32 + theme(legend.position = "none")) / (pc13  + theme(legend.position = "none")| pc23 + theme(legend.position = "none")))
output_filename <- if (endsWith(args$output, ".pdf")) {
    args$output
} else {
    paste0(args$output, ".pdf")
}
ggsave(output_filename, plot, width = 18, height = 8)
