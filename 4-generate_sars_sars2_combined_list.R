library(tidyverse)

sars <- read_csv("results/sars-ha1e-mcf7-combined-drugs.csv")
sars2 <- read_csv("results/dA549_2-ha1e-mcf7-combined-drugs.csv") %>% select(-ha1e, -mcf7)

combined <- sars %>%
  inner_join(sars2, by = "compound") %>%
  arrange(desc(mcf7), desc(ha1e)) %>%
  mutate(ha1e_rank = rank(desc(ha1e)),
         mcf7_rank = rank(desc(mcf7)),
         avg_rank = (ha1e_rank + mcf7_rank) / 2,
         dist = abs(ha1e_rank - mcf7_rank)) %>%
  arrange(avg_rank, dist) %>%
  mutate(rank = seq_along(compound)) %>%
  write_csv("results/sars-sars2-combined-drugs.csv")
