# Summarize the number of drugs for each threshold

library(tidyverse)

dirs <- list.dirs("data", recursive = F)
dirs <- dirs[str_detect(dirs, "\\d+")]

extract_files <- function(cell_line) {
  files <- list.files(file.path(cell_line, "consensus", "group"), full.names = T)
  pat <- "data/([A-Za-z0-9-]+)/consensus/group/([a-z]+)-(.+)-consensus.tsv"
  matched <- str_match(files, pat)
  cell_line <- matched[,2]
  group <- matched[,3]
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
                   drug_group = group,
                   threshold = threshold,
                   count = num_drugs)

  return(out_df)
}

complete <- dirs %>%
  map(~ extract_files(.x)) %>%
  bind_rows %>%
  pivot_wider(names_from = threshold, values_from = count) %>%
  write_csv("results/drugs_at_threshold_map.csv")
