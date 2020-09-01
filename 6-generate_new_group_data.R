# This script takes the list of individual drugs and then sorts them into groups.
# Once sorted into groups, it generates an averaged signature, filters it and then
# gets the connected signatures

library(tidyverse)
library(glue)
source("pipeline.R")

names <- c("quine", "nib", "vir", "azt", "los")

quine_group <- c("Chloroquine", "Hydroxychloroquine")
nib_group <- c("Baricitinib", "Fedratinib", "Ruxolitinib")
vir_group <- c("Lopinavir", "Ritonavir")
azt_group <- c("Azithromycin")
los_group <- c("Losartan")

all_groups <- list(quine_group, nib_group, vir_group, azt_group, los_group)

sig_map_lists <- setdiff(list.files("maps"), list.files("maps", "Drug"))
sig_map_names <- str_match(sig_map_lists, "(.*)\\-Signature*")[,2]
map_prefix <- "maps"

maps <- lapply(file.path(map_prefix, sig_map_lists), read_tsv)
names(maps) <- sig_map_names

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

a5491 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`A549-10uM-24h`, cell_line = sig_map_names[1])

a5492 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`A549-10uM-6h`, cell_line = sig_map_names[2])

ha1e2 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`HA1E-10uM-24h`, cell_line = sig_map_names[3])

ht29 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`HT29-10uM-24h`, cell_line = sig_map_names[4])

mcf72 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`MCF7-10uM-24h`, cell_line = sig_map_names[5])

pc3 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`PC3-10uM-24h`, cell_line = sig_map_names[6])

vcap1 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`VCAP-10uM-24h`, cell_line = sig_map_names[7])

vcap2 <- sapply(cutoffs, process_seed_groups, groups = all_groups, names = names, signature_map = maps$`VCAP-10uM-6h`, cell_line = sig_map_names[8])
