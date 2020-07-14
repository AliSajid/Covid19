# This script applies the threshold values for each consensus signature for drug and
# counts how many drugs are pulled up after the concordant drug search.

library(tidyverse)
source("pipeline.R")

groups <- c("Azt", "L", "Los", "Nib", "Quine", "Vir")

thresholds <- c(0, 0.26, 0.5, 0.85, 1)

ha1e_prefix <- file.path("dataset", "groups", "HA1E", "consensus")

ha1e_names <- paste("Group", groups, rep(thresholds, times=length(groups)), "l1000", "consensus", "connected", sep = "-")

ha1e_names <- paste(ha1e_names, "tsv", sep = ".")

ha1e_names <- sort(paste(ha1e_prefix, ha1e_names, sep = "/"))

file.exists(ha1e_names)

ha1e_data <- list()

col_spec <- cols(
  signatureid = col_character(),
  compound = col_character(),
  similarity = col_double()
)

for (i in 1:length(groups)) {
  for (j in 1:length(thresholds)) {
    n <- (i*5) -(5-j)
    ha1e_data[[groups[i]]][[as.character(thresholds[j])]] <- read_tsv(ha1e_names[n], col_types = col_spec)
  }
}

ha1e_uldata <- unlist(ha1e_data, recursive = F)

ha1e_counts <- as_tibble(lapply(ha1e_uldata, dim))%>% 
  t %>% 
  as_tibble(rownames = NA) %>% 
  rename(rows = V1, cols = V2) %>% 
  rownames_to_column("rowname") %>% 
  mutate(newname = sub("\\.", "/", rowname)) %>% 
  separate(col = newname, into = c("group", "threshold"), sep = "/") %>% 
  select(group, threshold, rows) %>% 
  pivot_wider(names_from = threshold, values_from = rows)

ha1e_counts %>% write_csv(file.path("results", "group_threshold_counts_ha1e.csv"))

mcf7_prefix <- file.path("dataset", "groups", "MCF7", "consensus")

mcf7_names <- paste("Group", groups, rep(thresholds, times=length(groups)), "l1000", "consensus", "connected", sep = "-")

mcf7_names <- paste(mcf7_names, "tsv", sep = ".")

mcf7_names <- sort(paste(mcf7_prefix, mcf7_names, sep = "/"))

file.exists(mcf7_names)

mcf7_data <- list()

col_spec <- cols(
  signatureid = col_character(),
  compound = col_character(),
  similarity = col_double()
)

for (i in 1:length(groups)) {
  for (j in 1:length(thresholds)) {
    n <- (i*5) -(5-j)
    mcf7_data[[groups[i]]][[as.character(thresholds[j])]] <- read_tsv(mcf7_names[n], col_types = col_spec)
  }
}

mcf7_uldata <- unlist(mcf7_data, recursive = F)

mcf7_counts <- as_tibble(lapply(mcf7_uldata, dim))%>% 
  t %>% 
  as_tibble(rownames = NA) %>% 
  rename(rows = V1, cols = V2) %>% 
  rownames_to_column("rowname") %>% 
  mutate(newname = sub("\\.", "/", rowname)) %>% 
  separate(col = newname, into = c("group", "threshold"), sep = "/") %>% 
  select(group, threshold, rows) %>% 
  pivot_wider(names_from = threshold, values_from = rows)

mcf7_counts %>% write_csv(file.path("results", "group_threshold_counts_mcf7.csv"))
