library(UpSetR)
library(tidyverse)

col_spec <- cols(
  compound = col_character(),
  avg = col_double(),
  sdev = col_double(),
  slog10 = col_double(),
  slog2 = col_double()
)

sars2 <- read_csv("results/ace2-summarized-dataset.csv", col_types = col_spec) %>%
  pull(compound)

sars2_fda <-   read_csv("results/ace2-summarized-dataset.csv", col_types = col_spec) %>%
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
) %>%
  pull(compound)

covidc <- read_csv("results/covidc-summarized-dataset.csv", col_types = col_spec) %>%
  pull(compound)

covidc_fda <- read_csv("results/covidc-summarized-dataset.csv", col_types = col_spec)  %>%
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
  ) %>%
  pull(compound)

covidm <- read_csv("results/covidm-summarized-dataset.csv", col_types = col_spec) %>%
  pull(compound)

covidm_fda <- read_csv("results/covidm-summarized-dataset.csv", col_types = col_spec)  %>%
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
  ) %>%
  pull(compound)

dataset_all <- list(sars2 = sars2,
                sars2_fda = sars2_fda,
                covidc = covidc,
                covidc_fda = covidc_fda,
                covidm = covidm,
                covidm_fda = covidm_fda)

dataset_unfiltered <- list(sars2 = sars2,
                           covidc = covidc,
                           covidm = covidm)

dataset_fda <- list(`GSE147507 CL` = sars2_fda,
                `GSE145926` = covidc_fda,
                `GSE147507 PS` = covidm_fda)

#upset(fromList(dataset_all), nsets = length(dataset_all))
#upset(fromList(dataset_unfiltered), nsets = length(dataset_unfiltered))

png("figures/sars2-china-mtsinai_comparison_upset.png", width = 10, height = 6.2, units = "in", res = 300)
upset(fromList(dataset_fda), nsets = length(dataset_fda),
      mainbar.y.label = "# of Common Drugs",
      sets.x.label = "# of Identified Drugs",
      scale.intersections = "identity")
dev.off()

jpeg("figures/sars2-china-mtsinai_comparison_upset.jpeg", width = 10, height = 6.2, units = "in", res = 300)
upset(fromList(dataset_fda), nsets = length(dataset_fda),
      mainbar.y.label = "# of Common Drugs",
      sets.x.label = "# of Identified Drugs",
      scale.intersections = "identity")
dev.off()

pdf("figures/sars2-china-mtsinai_comparison_upset.pdf", width = 10, height = 6.2)
upset(fromList(dataset_fda), nsets = length(dataset_fda),
      mainbar.y.label = "# of Common Drugs",
      sets.x.label = "# of Identified Drugs",
      scale.intersections = "identity")
dev.off()
