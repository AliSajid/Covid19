library(tidyverse)

sars <- read_csv("results/sars-groups-combined-drugs.csv")
sars2 <- read_csv("results/sars2-groups-combined-drugs.csv") %>% select(compound)

combined <- sars %>%
  inner_join(sars2, by = "compound") %>%
  arrange(desc(across(where(is.numeric)))) %>%
  write_csv("results/sars-sars2-common-drugs.csv")
