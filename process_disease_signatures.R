library(tidyverse)

col_spec <- cols(
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double(),
  Significance_pValue = col_double()
)

generate_outfile_name <- function(name, direction, cutoff, extension) {
  components <- str_split(str_split(name, pattern = "/")[[1]][5], "_")[[1]][c(1,3)]
  dataset <- components[1]
  disease <- str_split(components[2], pattern = "\\.")[[1]][1]
  outname <- paste(dataset, "L1000", disease, abs(cutoff), direction, "Signature", sep = "_")
  outname <- paste(outname, extension, sep = ".")
  return(outname)
}

generate_up_down_signatures <- function(file, cutoff_up, cutoff_down) {
  prefix_up <- paste("data", "Corona-Signature", "up-processed", sep = "/")
  prefix_down <- paste("data", "Corona-Signature", "down-processed", sep = "/")

  up_outfilename_tsv <- generate_outfile_name(file, "Up", cutoff_up, "tsv")
  up_outfilename_csv <- generate_outfile_name(file, "Up", cutoff_up, "csv")

  down_outfilename_tsv <- generate_outfile_name(file, "Down", cutoff_down, "tsv")
  down_outfilename_csv <- generate_outfile_name(file, "Down", cutoff_down, "csv")
  
  df <- read_csv(file, col_types = col_spec)
  up_data <- df %>% 
    filter(Value_LogDiffExp > cutoff_up)
  
  write_csv(up_data, paste(prefix_up, up_outfilename_csv, sep = "/"))
  write_tsv(up_data, paste(prefix_up, up_outfilename_tsv, sep = "/"))
  
  down_data <- read_csv(file, col_types = col_spec) %>% 
    filter(Value_LogDiffExp < cutoff_down)
  
  write_csv(down_data, paste(prefix_down, down_outfilename_csv, sep = "/"))
  write_tsv(down_data, paste(prefix_down, down_outfilename_tsv, sep = "/"))
  
}

files <- list.files("data/Corona-Signature/raw/", full.names = TRUE, pattern = "*.csv")
lapply(files, generate_up_down_signatures, cutoff_up = 0.26, cutoff_down = -0.26)
lapply(files, generate_up_down_signatures, cutoff_up = 0.5, cutoff_down = -0.5)
lapply(files, generate_up_down_signatures, cutoff_up = 0, cutoff_down = 0)

#TODO Add iLINCS API Support to automatically download the connected signatures