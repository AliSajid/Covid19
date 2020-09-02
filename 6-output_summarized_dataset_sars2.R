library(tidyverse)

sars2 <- read_csv("results/sars2-groups-combined-drugs.csv")

transformed <- sars2 %>%
  pivot_longer(cols = where(is.numeric)) %>%
  group_by(compound) %>%
  summarize(avg = mean(value), sdev = sd(value)) %>%
  mutate(slog10 = -log10(sdev),
         slog2 = -log2(sdev)) %>%
  write_csv("results/sars2-summarized-dataset.csv")
