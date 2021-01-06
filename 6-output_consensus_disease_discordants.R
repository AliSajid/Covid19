# This script takes all the discordant drugs for each disease signature and creates
# consensus lists

library(tidyverse)
library(glue)
library(UpSetR)

mark_nonfda <- function(df) {
  out <- df %>%
    mutate(
      fda_status = case_when(
        str_detect(compound, "CHEMBL") ~ "non_fda",
        str_detect(compound, "SCHEMBL") ~ "non_fda",
        str_detect(compound, "^\\d+") ~ "non_fda",
        str_detect(compound, "^[A-Z]\\d*\\w*\\-?\\s?\\d+") ~ "non_fda",
        str_detect(compound, "[Ii]nhibitor") ~ "non_fda",
        str_detect(compound, "^Broad") ~ "non_fda",
        str_detect(compound, "^BRD*") ~ "non_fda",
        str_detect(compound, "^UNII") ~ "non_fda",
        str_detect(compound, "omer") ~ "non_fda",
        str_detect(compound, "^Tyrphostin") ~ "non_fda",
        TRUE ~ "fda"
      )
    )
  out
}

diseases <- c("covidc", "covidm", "dA549_2", "dACE2_4")

cell_line <- "A549-10uM-24h"

thresholds <- c(0, 0.26, 0.5, 0.85, 1)

col_spec <- cols(
  signatureid = col_skip(),
  compound = col_character(),
  similarity = col_double()
)

paths <- expand_grid(diseases, cell_line, thresholds) %>%
  pmap_chr(~ file.path("data", "disease", ..1, "consensus", ..2, str_glue(..1, ..3, "consensus.tsv", .sep = "-")))

simple <- paths %>%
  map(~ read_tsv(.x, col_types = col_spec)) %>%
  map(~ mark_nonfda(.x)) %>%
  map2(rep(thresholds, 4), ~ mutate(.x, threshold = .y)) %>%
  map2(rep(diseases, each = 5), ~ mutate(.x, disease = .y)) %>%
  map2_dfr(paths, ~ mutate(.x, file = .y)) %>%
  write_csv(file.path("results", "combined_disease_discordant_dataset.csv"))


common <- simple %>%
  group_by(threshold, fda_status, compound) %>%
  summarise(num = n()) %>%
  ungroup() %>%
  group_by(threshold, fda_status) %>%
  arrange(desc(num)) %>%
  write_csv(file.path("results", "summarized_disease_discordant_dataset.csv"))

results <- read_csv("results/sars2-summarized-dataset.csv") %>%
  filter(avg > 0.5, sdev < 0.06) %>%
  mark_nonfda %>%
  filter(fda_status == "fda") %>%
  select(-fda_status) %>%
  pull(compound)

comparison <- common %>%
  filter(compound %in% results)

#names(simple) <- diseases

#upset(fromList(simple), nsets = length(diseases))
