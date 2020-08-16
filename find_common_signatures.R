library(tidyverse)

mcf7_files <- list.files("data/L1000-Signature/connected_signatures/MCF7", pattern = "*Consensus_Connected.csv", full.names = T)
ha1e_files <- list.files("data/L1000-Signature/connected_signatures/HA1E", pattern = "*Consensus_Connected.csv", full.names = T)

