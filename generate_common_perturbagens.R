# Final Result Generation Script
# This script takes in the connected perturbagens from the grouped drug signatures and compares them with the 
# ones in the SARS connected perturbagen list at 0.5 threshold.

library(tidyverse)

data_col_spec <- cols(
  signatureid = col_character(),
  compound = col_character(),
  similarity = col_double()
)

sars_concordants_HA1E <- read_tsv(
  file.path(
    "dataset", "disease", "HA1E", "sars-l1000-HA1E-0.5-consensus-connected.tsv"
    ),
  col_types = data_col_spec
  )

sars_concordants_MCF7 <- read_tsv(
  file.path(
    "dataset", "disease", "MCF7", "sars-l1000-MCF7-0.5-consensus-connected.tsv"
  ),
  col_types = data_col_spec
)

HA1E_prefix <- file.path("dataset", "groups", "HA1E", "consensus")
HA1E_files <- list.files(file.path("dataset", "groups", "HA1E", "consensus"), pattern = "0.85")
MCF7_prefix <- file.path("dataset", "groups", "MCF7", "consensus")
MCF7_files <- list.files(file.path("dataset", "groups", "MCF7", "consensus"), pattern = "0.85")

ids <- c("Azt", "L", "Los", "Nib", "Quine", "Vir")

remove_duplicates <- function(filename) {
  group_name <- str_match(filename, "Group-(.*)-0.85")[,2]
  dataframe <- read_tsv(filename, col_types = data_col_spec)
  output <- dataframe %>% 
    group_by(compound) %>% 
    filter(similarity == max(abs(similarity))) %>%
    slice_head(n = 1) %>% 
    select(compound, similarity) %>% 
    ungroup %>% 
    rename(!!group_name := similarity)
  
  return(output)
}

HA1E_data <- lapply(file.path(HA1E_prefix, HA1E_files), remove_duplicates)
names(HA1E_data) <- ids
MCF7_data <- lapply(file.path(MCF7_prefix, MCF7_files), remove_duplicates)
names(MCF7_data) <- ids

HA1E_cleaned <- reduce(HA1E_data, full_join, by = "compound")
HA1E_cleaned <- inner_join(HA1E_cleaned, sars_concordants_HA1E, by = "compound")
MCF7_cleaned <- reduce(MCF7_data, full_join, by = "compound")
MCF7_cleaned <- inner_join(MCF7_cleaned, sars_concordants_MCF7, by = "compound")

HA1E_passed <- HA1E_cleaned %>% 
  select(-L, -signatureid) %>% 
  rename(sars = similarity) %>% 
  mutate(nas = rowSums(is.na(.))) %>% 
  filter(nas < 4) %>% 
  write_csv(file.path("results", "HA1E-filtered-candidates.csv")) %>% 
  rowwise() %>% 
  mutate(meansim = mean(c(Azt, Los, Nib, Quine, Vir), na.rm = T))

MCF7_passed <-  MCF7_cleaned %>% 
  select(-L, -signatureid) %>% 
  rename(sars = similarity) %>% 
  mutate(nas = rowSums(is.na(.))) %>% 
  filter(nas < 4) %>% 
  write_csv(file.path("results", "MCF7-filtered-candidaates.csv")) %>% 
  rowwise() %>% 
  mutate(meansim = mean(c(Azt, Los, Nib, Quine, Vir), na.rm = T))

common <- intersect(MCF7_passed$compound, HA1E_passed$compound)
