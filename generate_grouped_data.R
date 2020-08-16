library(tidyverse)

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

get_drug_data <- function(members, cell_line, direction) {
  direction_dir = paste(str_to_lower(direction), "processed", sep = "-")
  prefix <- paste("data", 
                  "L1000-Signature",
                  direction_dir,
                  sep = "/")
  number <- length(members)
  dfs <- list()
    for (index in 1:number) {
      filename <- paste(
        paste(members[index], cell_line, direction, "Signature", sep = "_")
        , "csv", sep = ".")
      dfs[[index]] <- read_csv(paste(prefix, filename, sep = "/"))
    }
    return(dfs)
}

average_data <- function(dfs) {
  reduced <- reduce(dfs, inner_join, by = "Name_GeneSymbol")
  symbols <- reduced$Name_GeneSymbol
  avg_logexp <- reduced %>% 
    select(starts_with("Value")) %>% 
    rowMeans()
  avged <- tibble(Name_GeneSymbol = symbols,
                  Value_LogDiffExp = avg_logexp)
  return(avged)
}

process_grouping <- function(members, cell_line, name) {
  
  direction_dir = paste(str_to_lower(c("up", "down")), "processed", sep = "-")
  prefix <- paste("data", 
                  "L1000-Signature",
                  direction_dir,
                  sep = "/")
  
  up_filename_csv <- filename <- paste(
    paste("Group", name, cell_line, "Up", "Signature", sep = "_")
    , "csv", sep = ".")
  
  up_filename_tsv <- filename <- paste(
    paste("Group", name, cell_line, "Up", "Signature", sep = "_")
    , "tsv", sep = ".")
  
  down_filename_csv <- filename <- paste(
    paste("Group", name, cell_line, "Down", "Signature", sep = "_")
    , "csv", sep = ".")
  
  down_filename_tsv <- filename <- paste(
    paste("Group", name, cell_line, "Down", "Signature", sep = "_")
    , "tsv", sep = ".")
  
  
  data_up <- get_drug_data(members, cell_line, "Up")
  data_down <- get_drug_data(members, cell_line, "Down")
  
  avg_data_up <- average_data(data_up)
  avg_data_down <- average_data(data_down)
  
  write_csv(avg_data_up, paste(prefix[1], up_filename_csv, sep = "/"))
  write_tsv(avg_data_up, paste(prefix[1], up_filename_tsv, sep = "/"))
  
  write_csv(avg_data_down, paste(prefix[2], down_filename_csv, sep = "/"))
  write_tsv(avg_data_down, paste(prefix[2], down_filename_tsv, sep = "/"))
  
}

process_grouping_all <- function(members, cell_line, name) {
  
  # direction_dir = paste(str_to_lower(c("up", "down")), "processed", sep = "-")
  prefix <- paste("data", 
                  "L1000-Signature",
                  "raw",
                  sep = "/")
  
  out_filename_csv <- filename <- paste(
    paste("Group", name, cell_line, "Signature", sep = "_")
    , "csv", sep = ".")
  
  out_filename_tsv <- filename <- paste(
    paste("Group", name, cell_line, "Signature", sep = "_")
    , "tsv", sep = ".")
  
  data <- get_drug_data_all(members, cell_line)
  
  avg_data <- average_data(data)
  
  write_csv(avg_data, paste(prefix, out_filename_csv, sep = "/"))
  write_tsv(avg_data, paste(prefix, out_filename_tsv, sep = "/"))
  
}

get_drug_data_all <- function(members, cell_line) {
  prefix <- paste("data", 
                  "L1000-Signature",
                  "raw",
                  sep = "/")
  number <- length(members)
  dfs <- list()
  for (index in 1:number) {
    filename <- paste(
      paste(members[index], "L1000", cell_line, "Signature", sep = "-")
      , "csv", sep = ".")
    dfs[[index]] <- read_csv(paste(prefix, filename, sep = "/"))
  }
  return(dfs)
}

groups_ha1e <- c(quine_group_ha1e, nib_group_ha1e, vir_group_ha1e,
                 l01xe_group_ha1e, azt_group_ha1e, los_group_ha1e)

groups_mcf7 <- c(quine_group_mcf7, nib_group_mcf7, vir_group_mcf7,
                 l01xe_group_mcf7, azt_group_mcf7, los_group_mcf7)

process_grouping(quine_group_ha1e, "HA1E", "Quine")
process_grouping(nib_group_ha1e, "HA1E", "Nib")
process_grouping(vir_group_ha1e, "HA1E", "Vir")
process_grouping(l01xe_group_ha1e, "HA1E", "L01XE")
process_grouping(azt_group_ha1e, "HA1E", "Azt")
process_grouping(los_group_ha1e, "HA1E", "Los")


process_grouping(quine_group_mcf7, "MCF7", "Quine")
process_grouping(nib_group_mcf7, "MCF7", "Nib")
process_grouping(vir_group_mcf7, "MCF7", "Vir")
process_grouping(l01xe_group_mcf7, "MCF7", "L01XE")
process_grouping(azt_group_mcf7, "MCF7", "Azt")
process_grouping(los_group_mcf7, "MCF7", "Los")


# Complete
process_grouping_all(quine_group_ha1e, "HA1E", "Quine")
process_grouping_all(nib_group_ha1e, "HA1E", "Nib")
process_grouping_all(vir_group_ha1e, "HA1E", "Vir")
process_grouping_all(l01xe_group_ha1e, "HA1E", "L01XE")
process_grouping_all(azt_group_ha1e, "HA1E", "Azt")
process_grouping_all(los_group_ha1e, "HA1E", "Los")


process_grouping_all(quine_group_mcf7, "MCF7", "Quine")
process_grouping_all(nib_group_mcf7, "MCF7", "Nib")
process_grouping_all(vir_group_mcf7, "MCF7", "Vir")
process_grouping_all(l01xe_group_mcf7, "MCF7", "L01XE")
process_grouping_all(azt_group_mcf7, "MCF7", "Azt")
process_grouping_all(los_group_mcf7, "MCF7", "Los")
