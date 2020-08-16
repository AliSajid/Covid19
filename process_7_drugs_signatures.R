library(tidyverse)

col_spec <- cols(
  signatureID = col_skip(),
  PROBE = col_skip(),
  ID_geneid = col_skip(),
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double(),
  Significance_pvalue = col_double()
)

generate_up_down_signatures <- function(file, cutoff_up, cutoff_down) {
  name_components <- str_split(str_split(file, pattern = "/")[[1]][5], "-")[[1]]
  up_outfilename <- paste(name_components[1], name_components[3], "Up", "Signature", sep = "_")
  up_file <- paste("data", "L1000-Signature", "up-processed", up_outfilename, sep = "/")
  up_file_csv <- paste(up_file, "csv", sep = ".")
  up_file_tsv <- paste(up_file, "tsv", sep = ".")
  down_outfilename <- paste(name_components[1], name_components[3], "Down", "Signature", sep = "_")
  down_file <- paste("data", "L1000-Signature", "down-processed", down_outfilename, sep = "/")
  down_file_csv <- paste(down_file, "csv", sep = ".")
  down_file_tsv <- paste(down_file, "tsv", sep = ".")
  
  df <- read_csv(file, col_types = col_spec)
  up_data <- df %>% 
    filter(Value_LogDiffExp > cutoff_up)
  
  write_csv(up_data, up_file_csv)
  write_tsv(up_data, up_file_tsv)
  
  down_data <- read_csv(file, col_types = col_spec) %>% 
    filter(Value_LogDiffExp < cutoff_down)
  
  write_csv(down_data, down_file_csv)
  write_tsv(down_data, down_file_tsv)
  
}

files <- list.files("data/L1000-Signature/raw/", full.names = TRUE, pattern = "*.csv")
lapply(files, generate_up_down_signatures, cutoff_up = 0.85, cutoff_down = -0.85)

#TODO Add iLINCS API Support to automatically download the connected signatures