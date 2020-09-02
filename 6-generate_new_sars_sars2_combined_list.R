library(tidyverse)

sars <- read_csv("results/sars-new_groups-combined-drugs.csv")
sars2 <- read_csv("results/sars2-new_groups-combined-drugs.csv") %>% select(compound)

combined <- sars %>%
  inner_join(sars2, by = "compound") %>%
  arrange(across(where(is.numeric), desc)) %>%

  mutate(ha1e_rank = rank(desc(ha1e)),
         mcf7_rank = rank(desc(mcf7)),
         avg_rank = (ha1e_rank + mcf7_rank) / 2,
         dist = abs(ha1e_rank - mcf7_rank)) %>%
  arrange(avg_rank, dist) %>%
  mutate(rank = seq_along(compound)) %>%
  write_csv("results/sars-sars2-combined-drugs.csv")
