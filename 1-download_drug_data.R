# This script takes a list of signatures and downloads their data.
# It Then processes their data, generating the up and down signatures with our given threshold and finds concordant perturbagens to it.


library(tidyverse)
library(glue)
source("pipeline.R")


# This function takes a pertid and optionally thresholds and a library and then processes the signature from start to finish.
# This processing involves the following:
# 1. It downloads the designated L1000 Signature
# 2. It generates separate up and down filtered signatures for the given threshold
# 3. It then gets concordant drugs for those signatures
# 4. It then combines the results in one file
#
# At each step, it writes the intermediate results to a file

process_perturbagen <- function(pertid, threshold = 0.85, library = "LIB_5", cell_line) {
  basepath <- file.path("data", cell_line)
  sig_path <- file.path("data", "signatures", "drug")
  filtered_up_path <- file.path(basepath, "filtered", "drug", "up")
  filtered_down_path <- file.path(basepath, "filtered", "drug", "down")
  connected_up_path <- file.path(basepath, "connected", "drug", "up")
  connected_down_path <- file.path(basepath, "connected", "drug", "down")
  consensus_path <- file.path(basepath, "consensus", "drug")

  name_sig <- paste(pertid, "signature", sep = "-")
  name_filtered <- paste(pertid, threshold, "filtered", sep = "-")
  name_connected <- paste(pertid, threshold, "connected", sep = "-")
  name_consensus <- paste(pertid, threshold, "consensus", sep = "-")

  file_sig <- generate_name(sig_path, name_sig, "tsv")
  file_filtered_up <- generate_name(filtered_up_path, name_filtered, "tsv")
  file_filtered_down <- generate_name(filtered_down_path, name_filtered, "tsv")
  file_connected_up <- generate_name(connected_up_path, name_connected, "tsv")
  file_connected_down <- generate_name(connected_down_path, name_connected, "tsv")
  file_consensus <- generate_name(consensus_path, name_consensus, "tsv")

  print(glue("Processing {pertid} for threshold {threshold}"))
  if (!file.exists(file_sig)) {
    sig <- get_l1000_signature(pertid)
    write_tsv(sig, file_sig)
  } else {
    print(glue("{file_sig} already exists"))
    sig <- read_tsv(file_sig)
  }

  print(glue("Generating filtered signatures for {pertid}"))
  if (!file.exists(file_filtered_up)) {
    filtered_up <- generate_filtered_signature(sig, direction = "up", threshold = threshold)
    write_tsv(filtered_up, file_filtered_up)
  } else {
    print(glue("{file_filtered_up} already exists"))
    filtered_up <- read_tsv(file_filtered_up)
  }
  if (!file.exists(file_filtered_down)) {
    filtered_down <- generate_filtered_signature(sig, direction = "down", threshold = threshold)
    write_tsv(filtered_down, file_filtered_down)
  } else {
    print(glue("{file_filtered_down} already exists"))
    filtered_down <- read_tsv(file_filtered_down)
  }

  print(glue("Getting connected signatures for {pertid}"))
  if (!file.exists(file_connected_up)) {
    connected_up <- get_concordant_signatures(filtered_up, library = library)
    write_tsv(connected_up, file_connected_up)
  } else {
    print(glue("{file_connected_up} already exists"))
    connected_up <- read_tsv(file_connected_up)
  }
  if (!file.exists(file_connected_down)) {
    connected_down <- get_concordant_signatures(filtered_down, library = library)
    write_tsv(connected_down, file_connected_down)
  } else {
    print(glue("{file_connected_down} already exists"))
    connected_down <- read_tsv(file_connected_down)
  }

  print(glue("Generating consensus perturbagen list for {pertid}"))
  if (!file.exists(file_consensus)) {
    consensus <- generate_consensus_signature(connected_up, connected_down, cell_line = cell_line)
    write_tsv(consensus, file_consensus)
  } else {
    print(glue("{file_consensus} already exists"))
  }

}

# This function take a list of signature-to-drug maps and then processes each signature id.
process_seed_drugs <- function(datafile, threshold = 0.85, library = "LIB_5", cell_line) {
  pertids <- read_tsv(datafile) %>%
    pull(SignatureId)

  for (pert in pertids) {
    process_perturbagen(pert, threshold = threshold, library = library, cell_line = cell_line)
  }
}

cutoffs <- c(0, 0.26, 0.5, 0.85, 1)

ha1e <- sapply(cutoffs, process_seed_drugs, datafile = "maps/HA1E-Drug-Signature_Map.tsv", cell_line = "HA1E", library = "LIB_5")

mcf7 <- sapply(cutoffs, process_seed_drugs, datafile = "maps/MCF7-Drug-Signature_Map.tsv", cell_line = "MCF7", library = "LIB_5")
