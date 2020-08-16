source("pipeline.R")

raw_prefix <- paste("data", "Corona-Signature", "raw", sep = "/")

diseases <- paste(c("Influenza", "MERS", "SARS"), "csv", sep = ".")

files <- paste(raw_prefix, diseases, sep = "/")

new_col_names <- c("ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp", "Significance_pvalue")

process_disease <- function(disease_data, cutoff, cell_line) {
  up <- generate_filtered_signature(disease_data, direction = "up", threshold = cutoff)
  down <- generate_filtered_signature(disease_data, direction = "down", threshold = cutoff)
  
  up_concordant <- get_concordant_signatures(signature_df = up)
  down_concordant <- get_concordant_signatures(signature_df = down)
  
  concordant <- generate_consensus_signature(up_result = up_concordant,
                                             down_result = down_concordant,
                                            cell_line = cell_line,
                                            discordant = TRUE)
  
  return(concordant)
}

influenza <- read_csv(files[1], col_names = new_col_names)
mers <- read_csv(files[2], col_names = new_col_names)
sars <- read_csv(files[3], col_names = new_col_names)

dfs <- c("influenza", "mers", "sars")

cutoffs <- c(0, 0.26, 0.5)

cell_lines <- c("MCF7", "HA1E")

for (df in dfs) {
  for (cutoff in cutoffs) {
    for (cell_line in cell_lines) {
      file_name <- paste(paste(df, cutoff, cell_line, "Perturbagens", sep = "-"), "csv", sep = ".")
      data <- get(df)
      out <- process_disease(data, cutoff = cutoff, cell_line = cell_line)
      write.csv(out, paste(raw_prefix, file_name, sep = "/"))
    }
  }
}
