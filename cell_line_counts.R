library(tidyverse)

groups_HA1E <- list.files("data/L1000-Signature/connected_signatures/HA1E/", full.names = T)
groups_HA1E <- groups_HA1E[str_detect(groups_HA1E, "Group")]
groups_HA1E_consensus <- groups_HA1E[str_detect(groups_HA1E, "Consensus", negate = F)]
groups_HA1E_consensus <- groups_HA1E_consensus[str_detect(groups_HA1E_consensus, "tsv", negate = F)]
groups_HA1E <- groups_HA1E[str_detect(groups_HA1E, "Consensus", negate = T)]


groups_MCF7 <- list.files("data/L1000-Signature/connected_signatures/MCF7/", full.names = T)
groups_MCF7 <- groups_MCF7[str_detect(groups_MCF7, "Group")]
groups_MCF7_consensus <- groups_MCF7[str_detect(groups_MCF7, "Consensus", negate = F)]
groups_MCF7_consensus <- groups_MCF7_consensus[str_detect(groups_MCF7_consensus, "tsv")]
groups_MCF7 <- groups_MCF7[str_detect(groups_MCF7, "Consensus", negate = T)]

col_spec <- cols(
  .default = col_skip(),
  SignatureId = col_character(),
  Perturbagen = col_character(),
  CellLine = col_character()
  )

ha1e <- list()

for (i in 1:length(groups_HA1E)) {
  df <- read_tsv(groups_HA1E[i], col_types = col_spec)
  ha1e[[i]] <- df
}

ha1e_counts <- reduce(ha1e, union) %>% 
  filter(CellLine == "HA1E")

mcf7 <- list()

for (i in 1:length(groups_MCF7)) {
  df <- read_tsv(groups_MCF7[i], col_types = col_spec)
  mcf7[[i]] <- df
}

mcf7_counts <- reduce(mcf7, union) %>% 
  filter(CellLine == "MCF7")

col_spec_cons <- cols(
  .default = col_skip(),
  SignatureId = col_character(),
  Perturbagen = col_character()
)

ha1e_cons <- list()

for (i in 1:length(groups_HA1E_consensus)) {
  df <- read_tsv(groups_HA1E_consensus[i], col_types = col_spec_cons)
  ha1e_cons[[i]] <- df
}

ha1e_cons_counts <- reduce(ha1e_cons, union)

mcf7_cons <- list()

for (i in 1:length(groups_MCF7_consensus)) {
  df <- read_tsv(groups_MCF7_consensus[i], col_types = col_spec_cons)
  mcf7_cons[[i]] <- df
}

mcf7_cons_counts <- reduce(mcf7_cons, union)

