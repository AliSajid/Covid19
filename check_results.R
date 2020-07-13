library(tidyverse)

diseases <- c("influenza", "mers", "sars", "covidc", "covidm", "dNHBE_1",
              "dA549_2", "dA549_3", "dA549_p", "dACE2_4", "dCalu3_5")

thresholds <- c(0, 0.26, 0.5, 0.85, 1)

prefix <- "dataset/disease/HA1E"

names <- paste("l1000", "HA1E", thresholds, "consensus", "connected", sep = "-")
names <- paste(rep(diseases, each = length(names)), names, sep = "-")

names <- paste(names, "tsv", sep = ".")

names <- paste(prefix, names, sep = "/")

file.exists(names)

data <- list()

for (i in 1:length(diseases)) {
  for (j in 1:length(thresholds)) {
    n <- (i*5) -(5-j)
    data[[diseases[i]]][[as.character(thresholds[j])]] <- read_tsv(names[n])
  }
}

lapply(data, dim)
