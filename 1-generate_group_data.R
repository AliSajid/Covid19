# This script takes the list of individual drugs and then sorts them into groups.
# Once sorted into groups, it generates an averaged signature, filters it and then
# gets the connected signatures

library(tidyverse)
library(glue)
source("pipeline.R")

names <- c("quine", "nib", "vir", "azt", "los")

quine_group_mcf7 <- c("Chloroquine", "Hydroxychloroquine")
nib_group_mcf7 <- c("Baricitinib", "Fedratinib", "Ruxolitinib")
vir_group_mcf7 <- c("Lopinavir", "Ritonavir")
azt_group_mcf7 <- c("Azithromycin")
los_group_mcf7 <- c("Losartan")

mcf7_groups <- list(quine_group_mcf7, nib_group_mcf7, vir_group_mcf7, azt_group_mcf7, los_group_mcf7)

mcf7_sig_map <- read_tsv("maps/MCF7-Drug-Signature_Map.tsv")

quine_group_ha1e <- c("Chloroquine")
nib_group_ha1e <- c("Baricitinib", "Fedratinib", "Ruxolitinib")
vir_group_ha1e <- c("Lopinavir", "Ritonavir")
azt_group_ha1e <- c("Azithromycin")
los_group_ha1e <- c("Losartan")

ha1e_groups <- list(quine_group_ha1e, nib_group_ha1e, vir_group_ha1e, azt_group_ha1e, los_group_ha1e)

ha1e_sig_map <- read_tsv("maps/HA1E-Drug-Signature_Map.tsv")

strip_name <- function(call) {
  name <- as.character(call)[2]
  out <- str_extract(name, "[A-Za-z]+")
  return(out)
}

process_group <- function(members, name, signature_map, cell_line, threshold = 0.85, library = "LIB_5") {

  basepath <- file.path("data", cell_line)
  sig_path <- file.path("data", "signatures", "group")
  filtered_up_path <- file.path(basepath, "filtered", "group", "up")
  filtered_down_path <- file.path(basepath, "filtered", "group", "down")
  connected_up_path <- file.path(basepath, "connected", "group", "up")
  connected_down_path <- file.path(basepath, "connected", "group", "down")
  consensus_path <- file.path(basepath, "consensus", "group")

  name_sig <- paste(name, cell_line, "signature", sep = "-")
  name_filtered <- paste(name, threshold, "filtered", sep = "-")
  name_connected <- paste(name, threshold, "connected", sep = "-")
  name_consensus <- paste(name, threshold, "consensus", sep = "-")

  file_sig <- generate_name(sig_path, name_sig, "tsv")
  file_filtered_up <- generate_name(filtered_up_path, name_filtered, "tsv")
  file_filtered_down <- generate_name(filtered_down_path, name_filtered, "tsv")
  file_connected_up <- generate_name(connected_up_path, name_connected, "tsv")
  file_connected_down <- generate_name(connected_down_path, name_connected, "tsv")
  file_consensus <- generate_name(consensus_path, name_consensus, "tsv")

  sigs <- signature_map %>%
    filter(Drug %in% members) %>%
    pull(SignatureId)

  print(glue("Processing {name} group for threshold {threshold}"))
  if (!file.exists(file_sig)) {
    sigs_data <- lapply(sigs, function(x) get_l1000_signature(x) %>% select(Name_GeneSymbol, Value_LogDiffExp))
    group_sig <- reduce(sigs_data, inner_join, by = "Name_GeneSymbol") %>%
      pivot_longer(cols = starts_with("Value")) %>%
      group_by(Name_GeneSymbol) %>%
      summarise(Value_LogDiffExp = mean(value)) %>%
      write_tsv(file_sig)
  } else {
    print(glue("{file_sig} already exists"))
    group_sig <- read_tsv(file_sig)
  }

  print(glue("Generating filtered signatures for {name} group"))
  if (!file.exists(file_filtered_up)) {
    filtered_up <- generate_filtered_signature(group_sig, direction = "up", threshold = threshold)
    write_tsv(filtered_up, file_filtered_up)
  } else {
    print(glue("{file_filtered_up} already exists"))
    filtered_up <- read_tsv(file_filtered_up)
  }
  if (!file.exists(file_filtered_down)) {
    filtered_down <- generate_filtered_signature(group_sig, direction = "down", threshold = threshold)
    write_tsv(filtered_down, file_filtered_down)
  } else {
    print(glue("{file_filtered_down} already exists"))
    filtered_down <- read_tsv(file_filtered_down)
  }

  print(glue("Getting connected signatures for {name} group"))
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

  print(glue("Generating consensus perturbagen list for {name} group"))
  if (!file.exists(file_consensus)) {
    consensus <- generate_consensus_signature(connected_up, connected_down, cell_line = cell_line)
    write_tsv(consensus, file_consensus)
  } else {
    print(glue("{file_consensus} already exists"))
  }
}

process_seed_groups <- function(groups, names, signature_map, cell_line, threshold = 0.85, library = "LIB_5") {

  for (ind in 1:length(names)) {
    print(glue("Processing {names[ind]} group of {cell_line}"))
    process_group(members = groups[[ind]], name = names[[ind]], signature_map = signature_map, cell_line = cell_line, threshold = threshold, library = library)
  }
}

cutoffs <- c(0, 0.26, 0.5, 0.85, 1)

ha1e <- sapply(cutoffs, process_seed_groups, groups = ha1e_groups, names = names, signature_map = ha1e_sig_map, cell_line = "HA1E")

mcf7 <- sapply(cutoffs, process_seed_groups, groups = mcf7_groups, names = names, signature_map = mcf7_sig_map, cell_line = "MCF7")
