# Summarize the number of drugs for each threshold for diseases

library(tidyverse)

data_dir <- file.path("data", "disease")

diseases <- c("covidc", "covidm", "dA549_2")

consensus_dir <- "consensus"

cell_line_dirs <- rep(c("A549-10uM-24h", "A549-10uM-6h", "HA1E-10uM-24h", "HT29-10uM-24h", "MCF7-10uM-24h", "PC3-10uM-24h", "VCAP-10uM-24h", "VCAP-10uM-6h"), each = length(diseases))

dirs <- file.path(data_dir, diseases, consensus_dir, cell_line_dirs)

extract_files <- function(dir) {
  files <- list.files(dir, full.names = T)
  pat <- "data/disease/(\\w+)/consensus/([A-Za-z0-9-]+)/\\w+-([0-9.]+)-consensus.tsv"
  matched <- str_match(files, pat)
  disease <- matched[,2]
  cell_line <- matched[,3]
  threshold <- matched[,4]

  col_spec <- cols(
    signatureid = col_character(),
    compound = col_character(),
    similarity = col_double()
  )

  num_drugs <- files %>%
    map(~ read_tsv(.x, col_types = col_spec)) %>%
    map_dbl(~ nrow(.x))

  out_df <- tibble(cell_line = cell_line,
                   disease = disease,
                   threshold = threshold,
                   count = num_drugs)

  return(out_df)
}

complete <- dirs %>%
  map(~ extract_files(.x)) %>%
  bind_rows %>%
  pivot_wider(names_from = threshold, values_from = count) %>%
  write_csv("results/disease_at_threshold_map.csv")
