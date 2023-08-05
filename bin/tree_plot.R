#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description = "Plot tree with populations")
parser$add_argument("-p", "--populations", help = "Populations file", required = TRUE)
parser$add_argument("-t", "--tree", help = "Tree file", required = TRUE)
parser$add_argument("-o", "--output", help = "Output file", required = TRUE)
args <- parser$parse_args()

stopifnot(file.exists(args$populations))
stopifnot(file.exists(args$tree))

suppressPackageStartupMessages({
    library(ape)
    library(data.table)
    library(treeio)
    library(ggplot2)
    library(ggtree)
})

metadata <- fread(args$populations)
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
           "Mix/Unknown" = "seagreen")
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

# Loading the tree
tree <- treeio::read.newick(args$tree)
tree$edge.length[tree$edge.length < 0] <- 0

outgroup <- getMRCA(tree, c("AndeanFox01", "Coyote01"))
tree <- ladderize(root(tree, node = outgroup))

tree_table <- as.data.table(as_tibble(tree))
ctvt_node <- tree_table[label == "CTVT", node]
ht_node <- tree_table[label == "HT", node]

tree2 <- groupOTU(tree,
                  groups[rev(order(names(groups)))],
                  group_name = "population")
tree2$tip.label[tree2$tip.label == "HT"] <- "N-HT1"
p <- ggtree(tree2, aes(colour = population), ladderize = TRUE, right = TRUE) +
    scale_colour_manual(values = colours) + geom_tiplab(size = 1.5) +
    geom_nodelab(aes(label = round(as.integer(label))), size  = 1., colour = "black", nudge_x = 0.0001)
ggsave(args$output, plot = p, width = 12, height = 48)
