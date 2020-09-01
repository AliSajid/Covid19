# This script tries to find the best possible combinations for drugs

library(tidyverse)
source("pipeline.R")

files <- list.files("raw/drug_signature_list/", full.names = T)

df <- lapply(files, read_tsv)

df[[7]]$Perturbagen <- "Losartan"
df[[8]]$Perturbagen <- "Ritonavir"
df[[3]] <- df[[3]] %>%
  filter(Perturbagen == "Chloroquine")

drug_data <- bind_rows(df) %>%
  select(SignatureId, Perturbagen, Concentration, CellLine, Time) %>%
  group_by(Concentration, CellLine, Time) %>%
  write_csv("data/signatures/sig_lists/drug-signature-attribute-mapping.csv")

counts <- drug_data %>%
  summarise(count = n(),
            unique = n_distinct(Perturbagen)) %>%
  mutate(complete = if_else(unique > 5, "Yes", "No")) %>%
  filter(complete == "Yes") %>%
  write_csv("data/signatures/sig_lists/maximum_representation_combinations.csv")


filter_drug_data <- function(variables) {
  file <- paste(variables["CellLine"], variables["Concentration"], variables["Time"], "Signature_Map", sep = "-")
  file_name <- file.path("maps", paste(file, "tsv", sep = "."))

  drug_data %>%
    filter(Concentration == variables["Concentration"],
           Time == variables["Time"],
           CellLine == variables["CellLine"]) %>%
    write_tsv(file_name)
}

filtered_lists <- apply(counts, 1, filter_drug_data)


