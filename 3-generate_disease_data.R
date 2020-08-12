# This script takes the disease signatures and processses them one by one.
# It filters the signatures, gets concordants and the generates a consensus list of drugs.
#

library(tidyverse)
library(glue)
source("pipeline.R")

diseases <- c("covidc", "covidm", "dNHBE_1", "dA549_2", "dA549_3", "dA549_p",
              "dACE2_4", "dCalu3_5", "influenza", "mers", "sars")

process_disease <- function(disease, threshold = 0.85, library = "LIB_5", cell_line) {
  basepath <- file.path("data", "disease", disease)
  sig_path <- file.path("data", "signatures", "disease")
  filtered_up_path <- file.path(basepath, "filtered", "up")
  filtered_down_path <- file.path(basepath, "filtered", "down")
  connected_up_path <- file.path(basepath, "connected", "up")
  connected_down_path <- file.path(basepath, "connected", "down")
  consensus_path <- file.path(basepath, "consensus", cell_line)

  name_sig <- paste(disease, "signature", sep = "-")
  name_filtered <- paste(disease, threshold, "filtered", sep = "-")
  name_connected <- paste(disease, threshold, "connected", sep = "-")
  name_consensus <- paste(disease, threshold, "consensus", sep = "-")

  file_sig <- generate_name(sig_path, name_sig, "tsv")
  file_filtered_up <- generate_name(filtered_up_path, name_filtered, "tsv")
  file_filtered_down <- generate_name(filtered_down_path, name_filtered, "tsv")
  file_connected_up <- generate_name(connected_up_path, name_connected, "tsv")
  file_connected_down <- generate_name(connected_down_path, name_connected, "tsv")
  file_consensus <- generate_name(consensus_path, name_consensus, "tsv")

  src_path <- file.path("raw")
  src_file <- generate_name(src_path, name_sig, "tsv")

  print(glue("Processing {disease} for threshold {threshold} and cell line {cell_line}"))
  sig <- read_tsv(src_file) %>%
    write_tsv(file_sig)

  print(glue("Generating filtered signatures for {disease}"))
  filtered_up <- generate_filtered_signature(sig, direction = "up", threshold = threshold)
  write_tsv(filtered_up, file_filtered_up)
  filtered_down <- generate_filtered_signature(sig, direction = "down", threshold = threshold)
  write_tsv(filtered_down, file_filtered_down)

  print(glue("Getting connected signatures for {disease}"))
  connected_up <- get_concordant_signatures(filtered_up, library = library)
  write_tsv(connected_up, file_connected_up)
  connected_down <- get_concordant_signatures(filtered_down, library = library)
  write_tsv(connected_down, file_connected_down)

  print(glue("Generating consensus perturbagen list for {disease}"))
  consensus <- generate_consensus_signature(connected_up, connected_down, cell_line = cell_line, discordant = TRUE)
  write_tsv(consensus, file_consensus)

}

process_diseases <- function(diseases, threshold = 0.85, library = "LIB_5") {
  for (dis in diseases) {
    for (cell_line in c("HA1E", "MCF7")) {
      process_disease(dis, threshold = threshold, library = library, cell_line = cell_line)
    }
  }
}

cutoffs <- c(0, 0.26, 0.5, 0.85, 1)
results <- sapply(cutoffs, process_diseases, diseases = diseases)
