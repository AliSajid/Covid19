library(tidyverse)
source("pipeline.R")

quine_group_mcf7 <- c("Chloroquine", "Hydroxychloroquine")
nib_group_mcf7 <- c("Baricitinib", "Fedratinib", "Ruxolitinib")
vir_group_mcf7 <- c("Lopinavir", "Ritonavir")
l01xe_group_mcf7 <- c("Fedratinib", "Ruxolitinib")
azt_group_mcf7 <- c("Azithromycin")
los_group_mcf7 <- c("Losartan")

quine_group_ha1e <- c("Chloroquine")
nib_group_ha1e <- c("Baricitinib", "Fedratinib", "Ruxolitinib")
vir_group_ha1e <- c("Lopinavir", "Ritonavir")
l01xe_group_ha1e <- c("Fedratinib", "Ruxolitinib")
azt_group_ha1e <- c("Azithromycin")
los_group_ha1e <- c("Losartan")

prefix <- "dataset/groups"
raw_prefix <- paste(prefix, c("HA1E", "MCF7"), "raw", sep = "/")
up_prefix <- paste(prefix, c("HA1E", "MCF7"), "up", sep = "/")
down_prefix <- paste(prefix, c("HA1E", "MCF7"), "down", sep = "/")
concordant_prefix <- paste(prefix, c("HA1E", "MCF7"), "concordant", sep = "/")
consensus_prefix <- paste(prefix, c("HA1E", "MCF7"), "consensus", sep = "/")

HA1E_signatures <- read_tsv("HA1E-Drug-Signature_Map.tsv") %>% select(SignatureId, Drug)
MCF7_signatures <- read_tsv("MCF7-Drug-Signature_Map.tsv") %>% select(SignatureId, Drug)

HA1E_groups <- c("quine_group_ha1e", "nib_group_ha1e", "vir_group_ha1e",
                 "l01xe_group_ha1e", "azt_group_ha1e", "los_group_ha1e")

MCF7_groups <- c("quine_group_mcf7", "nib_group_mcf7", "vir_group_mcf7",
                 "l01xe_group_mcf7", "azt_group_mcf7", "los_group_mcf7")

cutoffs <- c(0, 0.26, 0.5, 0.85, 1)

transform_name_to_signatureid <- function(drug_name, cell_line) {
  name_map_file <- paste(cell_line, "Drug-Signature_Map.tsv", sep = "-")
  result <- read_tsv(name_map_file) %>%
    filter(Drug == drug_name) %>%
    select(SignatureId) %>%
    unlist
  
  return(result)
}

get_group_data <- function(members, cell_line) {
  prefix <- paste("dataset", 
                  "drugs",
                  cell_line,
                  "raw",
                  sep = "/")
  number <- length(members)
  dfs <- list()
  for (index in 1:number) {
    signature_name <- transform_name_to_signatureid(members[index], cell_line)
    filename <- generate_name(prefix, 
                              paste(signature_name, "l1000", "signature", sep = "-"),
                              "tsv")
    dfs[[index]] <- read_tsv(filename)
  }
  return(dfs)
}

average_data <- function(dfs) {
  reduced <- reduce(dfs, inner_join, by = "Name_GeneSymbol")
  avg_logexp <- reduced %>% 
    select(Name_GeneSymbol, starts_with("Value")) %>% 
    column_to_rownames("Name_GeneSymbol") %>% 
    rowMeans(na.rm = T)
  avged <- tibble(Name_GeneSymbol = names(avg_logexp),
                  Value_LogDiffExp = avg_logexp)
  return(avged)
}



# Process Signatures for HA1E Cell Lines
if (!dir.exists(raw_prefix[1])) {
  dir.create(raw_prefix[1], recursive = T)
}

if (!dir.exists(up_prefix[1])) {
  dir.create(up_prefix[1], recursive = T)
}

if (!dir.exists(down_prefix[1])) {
  dir.create(down_prefix[1], recursive = T)
}

if (!dir.exists(consensus_prefix[1])) {
  dir.create(consensus_prefix[1], recursive = T)
}

if (!dir.exists(concordant_prefix[1])) {
  dir.create(concordant_prefix[1], recursive = T)
}

for (index in 1:length(HA1E_groups)) {
  
  print(paste("Processing", HA1E_groups[index]))
  
  for (cutoff in cutoffs) {
    
  print(paste("Processing", cutoff))
  group <- HA1E_groups[index]
  name <- str_extract(group, "[A-Za-z]*")
  members <- get(group)
  
  raw_output_file_name <- generate_name(raw_prefix[1], paste("Group", str_to_title(name),
                                                         "l1000", "signature", sep = "-"), "tsv")
  up_output_file_name <- generate_name(up_prefix[1], paste("Group", str_to_title(name), cutoff,
                                                             "l1000", "up", "signature", sep = "-"), "tsv")
  down_output_file_name <- generate_name(down_prefix[1], paste("Group", str_to_title(name), cutoff,
                                                             "l1000", "down", "signature", sep = "-"), "tsv")
  concordant_up_output_file_name <- generate_name(concordant_prefix[1], paste("Group", str_to_title(name), cutoff,
                                                                              "l1000", "concordant", "up", "connected", sep = "-"), "tsv")
  concordant_down_output_file_name <- generate_name(concordant_prefix[1], paste("Group", str_to_title(name), cutoff,
                                                                              "l1000", "concordant", "down", "connected", sep = "-"), "tsv")
  consensus_output_file_name <- generate_name(consensus_prefix[1], paste("Group", str_to_title(name), cutoff,
                                                             "l1000", "consensus", "connected", sep = "-"), "tsv")
  group_data <- get_group_data(members, "HA1E")
  avg_data <- average_data(group_data)
  write_tsv(avg_data, raw_output_file_name)
  
  up_data <- generate_filtered_signature(avg_data, direction = "up", threshold = cutoff)
  write_tsv(up_data, up_output_file_name)
  
  down_data <- generate_filtered_signature(avg_data, direction = "down", threshold = cutoff)
  write_tsv(down_data, down_output_file_name)
  
  
  up_concordant <- get_concordant_signatures(signature_df = up_data)
  down_concordant <- get_concordant_signatures(signature_df = down_data)
  write_tsv(up_concordant, concordant_up_output_file_name)
  write_tsv(down_concordant, concordant_down_output_file_name)
  
  consensus <- generate_consensus_signature(up_concordant, down_concordant, cell_line = "HA1E")
  write_tsv(consensus, consensus_output_file_name)
  }
}

# Process Signatures for MCF7 Cell Lines
if (!dir.exists(raw_prefix[2])) {
  dir.create(raw_prefix[2], recursive = T)
}

if (!dir.exists(up_prefix[2])) {
  dir.create(up_prefix[2], recursive = T)
}

if (!dir.exists(down_prefix[2])) {
  dir.create(down_prefix[2], recursive = T)
}

if (!dir.exists(consensus_prefix[2])) {
  dir.create(consensus_prefix[2], recursive = T)
}

if (!dir.exists(concordant_prefix[2])) {
  dir.create(concordant_prefix[2], recursive = T)
}

for (index in 1:length(MCF7_groups)) {
  print(paste("Processing", MCF7_groups[index]))
  for (cutoff in cutoffs) {
    print(paste("Processing", cutoff))
  group <- MCF7_groups[index]
  name <- str_extract(group, "[A-Za-z]*")
  members <- get(group)
  
  raw_output_file_name <- generate_name(raw_prefix[2], paste("Group", str_to_title(name),
                                                             "l1000", "signature", sep = "-"), "tsv")
  up_output_file_name <- generate_name(up_prefix[2], paste("Group", str_to_title(name), cutoff,
                                                           "l1000", "up", "signature", sep = "-"), "tsv")
  down_output_file_name <- generate_name(down_prefix[2], paste("Group", str_to_title(name), cutoff,
                                                               "l1000", "down", "signature", sep = "-"), "tsv")
  concordant_up_output_file_name <- generate_name(concordant_prefix[2], paste("Group", str_to_title(name), cutoff,
                                                                              "l1000", "concordant", "up", "connected", sep = "-"), "tsv")
  concordant_down_output_file_name <- generate_name(concordant_prefix[2], paste("Group", str_to_title(name), cutoff,
                                                                                "l1000", "concordant", "down", "connected", sep = "-"), "tsv")
  consensus_output_file_name <- generate_name(consensus_prefix[2], paste("Group", str_to_title(name), cutoff,
                                                                         "l1000", "consensus", "connected", sep = "-"), "tsv")
  group_data <- get_group_data(members, "MCF7")
  avg_data <- average_data(group_data)
  write_tsv(avg_data, raw_output_file_name)
  
  up_data <- generate_filtered_signature(avg_data, direction = "up", threshold = cutoff)
  write_tsv(up_data, up_output_file_name)
  
  down_data <- generate_filtered_signature(avg_data, direction = "down", threshold = cutoff)
  write_tsv(down_data, down_output_file_name)
  
  
  up_concordant <- get_concordant_signatures(signature_df = up_data)
  down_concordant <- get_concordant_signatures(signature_df = down_data)
  write_tsv(up_concordant, concordant_up_output_file_name)
  write_tsv(down_concordant, concordant_down_output_file_name)
  
  consensus <- generate_consensus_signature(up_concordant, down_concordant, cell_line = "MCF7")
  write_tsv(consensus, consensus_output_file_name)
  }
}

