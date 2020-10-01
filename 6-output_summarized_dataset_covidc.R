library(tidyverse)

covidc <- read_csv("results/covidc-groups-combined-drugs.csv")

transformed <- covidc %>%
  pivot_longer(cols = where(is.numeric)) %>%
  group_by(compound) %>%
  summarize(avg = mean(value), sdev = sd(value)) %>%
  mutate(slog10 = -log10(sdev),
         slog2 = -log2(sdev)) %>%
  write_csv("results/covidc-summarized-dataset.csv")

fda_approved <- transformed %>%
  filter(str_detect(compound, "CHEMBL", negate = T),
         str_detect(compound, "SCHEMBL", negate = T),
         str_detect(compound, "^\\d+", negate = T),
         str_detect(compound, "^[A-Z]\\d*\\w*\\-?\\s?\\d+", negate = T),
         str_detect(compound, "[Ii]nhibitor", negate = T),
         str_detect(compound, "^Broad", negate = T),
         str_detect(compound, "^BRD*", negate = T),
         str_detect(compound, "^UNII", negate = T),
         str_detect(compound, "omer", negate = T),
         str_detect(compound, "^Tyrphostin", negate = T)
  )
