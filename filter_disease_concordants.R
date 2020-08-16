# Obsoleted

library(tidyverse)

filter_data <-
  function(dataframe_up,
           dataframe_down,
           cell_line,
           cutoff) {
    dataframe <- bind_rows(dataframe_up, dataframe_down)
    output <- dataframe %>%
      filter(CellLine == cell_line, Concordance < cutoff) %>%
      group_by(Perturbagen) %>%
      filter(Concordance == max(Concordance)) %>%
      ungroup() %>%
      select(SignatureId, Perturbagen, Concordance)
    return(output)
  }


generate_input_names <- function(disease, cutoff) {
  
  directions <- c("Up", "Down")
  prefix <- paste("data", "Corona-Signature", "connected_signatures", sep = "/")
  
  filenames <- expand.grid(disease, cutoff, directions, "Connected", "tsv") %>% 
    unite(name, Var1:Var4, sep = "_") %>% 
    unite(name, sep = ".")
  
  filenames <- filenames$name
  
  files <- paste(prefix, filenames, sep = "/")
  return(files)
}

generate_output_names <- function(disease, cutoff, cell_line, extension) {
  prefix <- paste("data", "Corona-Signature", "connected_signatures", sep = "/")
  
  filenames <- expand.grid(disease, cutoff, cell_line, "Consensus", "Connected", extension) %>% 
    unite(name, Var1:Var5, sep = "_") %>% 
    unite(name, sep = ".")
  
  filenames <- filenames$name
  
  files <- paste(prefix, filenames, sep = "/")
  return(files)
}

process_cutoff <- function(cutoff, disease, threshold) {
  input_file_names <- generate_input_names(disease, cutoff)
  data_up <- read_tsv(input_file_names[1])
  data_down <- read_tsv(input_file_names[2])
  
  out_mcf7_tsv_name <- generate_output_names(disease, cutoff, "MCF7", "tsv")
  out_mcf7_csv_name <- generate_output_names(disease, cutoff, "MCF7", "csv")
  
  out_ha1e_tsv_name <- generate_output_names(disease, cutoff, "HA1E", "tsv")
  out_ha1e_csv_name <- generate_output_names(disease, cutoff, "HA1E", "csv")
  
  consensus_mcf7 <- filter_data(data_up, data_down, "MCF7", threshold)
  consensus_ha1e <- filter_data(data_up, data_down, "HA1E", threshold)
  
  write_csv(consensus_mcf7, out_mcf7_csv_name)
  write_tsv(consensus_mcf7, out_mcf7_tsv_name)

  write_csv(consensus_ha1e, out_ha1e_csv_name)
  write_tsv(consensus_ha1e, out_ha1e_tsv_name)
  
}

process_disease <- function(disease, threshold = -0.321) {
  cutoffs <- c("0", "0.5", "0.26")
  sapply(cutoffs, process_cutoff, disease = disease, threshold = threshold)
}

diseases <- c("Influenza", "MERS", "SARS")

sapply(diseases, process_disease)
