# This script applies the threshold values for each consensus signature and
# counts how many drugs are pulled up after the concordant drug search.

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

col_spec <- cols(
  signatureid = col_character(),
  compound = col_character(),
  similarity = col_double()
)

for (i in 1:length(diseases)) {
  for (j in 1:length(thresholds)) {
    n <- (i*5) -(5-j)
    data[[diseases[i]]][[as.character(thresholds[j])]] <- read_tsv(names[n], col_types = col_spec)
  }
}

uldata <- unlist(data, recursive = F)

counts <- as_tibble(lapply(uldata, dim))%>% 
  t %>% 
  as_tibble(rownames = NA) %>% 
  rename(rows = V1, cols = V2) %>% 
  rownames_to_column("rowname") %>% 
  mutate(newname = sub("\\.", "/", rowname)) %>% 
  separate(col = newname, into = c("disease", "threshold"), sep = "/") %>% 
  select(disease, threshold, rows) %>% 
  pivot_wider(names_from = threshold, values_from = rows)

counts %>% write_csv(file.path("results", "disease_threshold_counts.csv"))
