#!/usr/bin/env Rscript

# Prepare all combinations of populations
comprehensive <- c(
    "America_pool",
    "Ashkelon_pool",
    "Baikal_pool",
    "Croatia_Eneolithic.ALPO01",
    "Croatia_Eneolithic.SOTN01",
    "GermanShepherdDog",
    "Germany_CordedWare.CTC",
    "Germany_Neolithic.HXH",
    "Greece_Neolithic.OL4222",
    "Iran_Chalcolithic.AL2571",
    "Ireland_Neolithic.Newgrange",
    "Israel_Neolithic.THRZ02",
    "Italy_BronzeAge.AL2397",
    "Karelia_Mesolithic.OL4061",
    "NewGuineaSingingDog",
    "Samara_BronzeAge.C5",
    "Serbia_Neolithic.AL2946",
    "Spain_Neolithic.OL4029",
    "SwedenPWC_pool",
    "Sweden_4k.C94",
    "Sweden_BronzeAge.C62",
    "Sweden_Neolithic.C88"
)

core <- c(
    "America_pool",
    "Baikal_pool",
    "GermanShepherdDog",
    "Iran_Chalcolithic.AL2571",
    "Israel_Neolithic.THRZ02",
    "Karelia_Mesolithic.OL4061",
    "NewGuineaSingingDog",
    "Samara_BronzeAge.C5",
    "Sweden_Neolithic.C88"
)

chosen_subset <- c(
    "America_pool",
    "Baikal_pool",
    "Israel_Neolithic.THRZ02",
    "Iran_Chalcolithic.AL2571",
    "Samara_BronzeAge.C5",
    "NewGuineaSingingDog",
    "GermanShepherdDog"
)

lefts <- lapply(1:4, function(n) t(combn(chosen_subset, n)))


lines <- c("LEFT\tRIGHT")
for (sublist in lefts) {
    for (i in seq_len(nrow(sublist))) {
        l <- sublist[i, ]
        r <- c("CoyoteCalifornia", setdiff(comprehensive, l))
        line <- paste(paste(l, collapse = ","), paste(r, collapse = ","), sep = "\t")
        lines <- append(lines, line)
    }
}

conn <- file("adm_combinations_comprehensive.txt")
writeLines(lines, conn)
close(conn)


lines <- c("LEFT\tRIGHT")
for (sublist in lefts) {
    for (i in seq_len(nrow(sublist))) {
        l <- sublist[i, ]
        r <- c("CoyoteCalifornia", setdiff(core, l))
        line <- paste(paste(l, collapse = ","), paste(r, collapse = ","), sep = "\t")
        lines <- append(lines, line)
    }
}

conn <- file("adm_combinations_core.txt")
writeLines(lines, conn)
close(conn)
