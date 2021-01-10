# This scripts creates a list of all pairwise comparisons between covidc, covidm,
# the a549-sars2 and the ace2 consensus drug lists.


library(tidyverse)


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

top_drugs <-
  c(
    "Gemcitabine",
    "Trametinib",
    "Withaferin A",
    "Saracatinib",
    "Erlotinib",
    "Alvocidib",
    "Itrazole",
    "Elesclomol",
    "Dasatinib",
    "Panobinostat",
    "Gallocatechin Gallate",
    "Genistein",
    "Imatinib",
    "Dexamethasone Acetate",
    "Simvastatin",
    "Sirolimus",
    "Tamoxifen"
  )

covidc <- read_csv("results/covidc-summarized-dataset.csv") %>%
  select(compound, avg, sdev) %>%
  mark_nonfda() %>%
  mutate(covidc = 2 ^ 0)
covidm <- read_csv("results/covidm-summarized-dataset.csv") %>%
  select(compound, avg, sdev) %>%
  mark_nonfda() %>%
  mutate(covidm = 2 ^ 1)
a549 <- read_csv("results/sars2-summarized-dataset.csv") %>%
  select(compound, avg, sdev) %>%
  mark_nonfda() %>%
  mutate(a549 = 2 ^ 2)
ace2 <- read_csv("results/ace2-summarized-dataset.csv") %>%
  select(compound, avg, sdev) %>%
  mark_nonfda() %>%
  mutate(ace2 = 2 ^ 3)
selected_list <- read_csv("results/sars2-summarized-dataset.csv") %>%
  select(compound, avg, sdev) %>%
  mark_nonfda() %>%
  filter(avg > 0.5, sdev < 0.06, fda_status == "fda") %>%
  mutate(selected = 2 ^ 4,
         top = if_else(compound %in% top_drugs, 2 ^ 5, 0))




all_drugs <- covidc %>%
  full_join(covidm, by = c("compound", "avg", "sdev", "fda_status")) %>%
  full_join(a549, by = c("compound", "avg", "sdev", "fda_status")) %>%
  full_join(ace2, by = c("compound", "avg", "sdev", "fda_status")) %>%
  full_join(selected_list, by = c("compound", "avg", "sdev", "fda_status")) %>%
  mutate(covidc = if_else(is.na(covidc), 0, covidc),
         covidm = if_else(is.na(covidm), 0, covidm),
         a549 = if_else(is.na(a549), 0, a549),
         ace2 = if_else(is.na(ace2), 0, ace2),
         selected = if_else(is.na(selected), 0, selected),
         top = if_else(is.na(top), 0, top),
         ident = covidc + covidm + a549 + ace2 + selected + top) %>%
  write_csv("results/all_drugs_combinations_comparison.csv")

