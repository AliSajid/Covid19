# This script runs first and creates the necessary directory structure
#

source("pipeline.R")

data_dir <- file.path("data")

ha1e_dir <- file.path(data_dir, "HA1E")
mcf7_dir <- file.path(data_dir, "MCF7")

ha1e_filtered_dir <- file.path(ha1e_dir, "filtered")
ha1e_filtered_drug_dir <- file.path(ha1e_filtered_dir, "drug")
ha1e_filtered_group_dir <- file.path(ha1e_filtered_dir, "group")
mcf7_filtered_dir <- file.path(mcf7_dir, "filtered")
mcf7_filtered_drug_dir <- file.path(mcf7_filtered_dir, "drug")
mcf7_filtered_group_dir <- file.path(mcf7_filtered_dir, "group")

ha1e_filtered_drug_up_dir <- file.path(ha1e_filtered_drug_dir, "up")
ha1e_filtered_drug_down_dir <- file.path(ha1e_filtered_drug_dir, "down")
mcf7_filtered_drug_up_dir <- file.path(mcf7_filtered_drug_dir, "up")
mcf7_filtered_drug_down_dir <- file.path(mcf7_filtered_drug_dir, "down")

ha1e_filtered_group_up_dir <- file.path(ha1e_filtered_group_dir, "up")
ha1e_filtered_group_down_dir <- file.path(ha1e_filtered_group_dir, "down")
mcf7_filtered_group_up_dir <- file.path(mcf7_filtered_group_dir, "up")
mcf7_filtered_group_down_dir <- file.path(mcf7_filtered_group_dir, "down")

ha1e_connected_dir <- file.path(ha1e_dir, "connected")
ha1e_connected_drug_dir <- file.path(ha1e_connected_dir, "drug")
ha1e_connected_group_dir <- file.path(ha1e_connected_dir, "group")
mcf7_connected_dir <- file.path(mcf7_dir, "connected")
mcf7_connected_drug_dir <- file.path(mcf7_connected_dir, "drug")
mcf7_connected_group_dir <- file.path(mcf7_connected_dir, "group")


ha1e_connected_drug_up_dir <- file.path(ha1e_connected_drug_dir, "up")
ha1e_connected_drug_down_dir <- file.path(ha1e_connected_drug_dir, "down")
mcf7_connected_drug_up_dir <- file.path(mcf7_connected_drug_dir, "up")
mcf7_connected_drug_down_dir <- file.path(mcf7_connected_drug_dir, "down")

ha1e_connected_group_up_dir <- file.path(ha1e_connected_group_dir, "up")
ha1e_connected_group_down_dir <- file.path(ha1e_connected_group_dir, "down")
mcf7_connected_group_up_dir <- file.path(mcf7_connected_group_dir, "up")
mcf7_connected_group_down_dir <- file.path(mcf7_connected_group_dir, "down")

ha1e_consensus_dir <- file.path(ha1e_dir, "consensus")
ha1e_consensus_drug_dir <- file.path(ha1e_consensus_dir, "drug")
ha1e_consensus_group_dir <- file.path(ha1e_consensus_dir, "group")
mcf7_consensus_dir <- file.path(mcf7_dir, "consensus")
mcf7_consensus_drug_dir <- file.path(mcf7_consensus_dir, "drug")
mcf7_consensus_group_dir <- file.path(mcf7_consensus_dir, "group")

disease_dir <- file.path(data_dir, "disease")
diseases <- c("influenza", "mers", "sars", "covidc", "covidm", "dNHBE_1",
              "dA549_2", "dA549_3", "dA549_p", "dACE2_4", "dCalu3_5")

diseases_dirs <- file.path(disease_dir, diseases)

diseases_filtered <- file.path(diseases_dirs, "filtered")
diseases_filtered_up <- file.path(diseases_filtered, "up")
diseases_filtered_down <- file.path(diseases_filtered, "down")

diseases_connected <- file.path(diseases_dirs, "connected")
diseases_connected_up <- file.path(diseases_connected, "up")
diseases_connected_down <- file.path(diseases_connected, "down")

diseases_consensus <- file.path(diseases_dirs, "consensus")
diseases_consensus_ha1e <- file.path(diseases_consensus, "HA1E")
diseases_consensus_mcf7 <- file.path(diseases_consensus, "MCF7")

figures_dir <- file.path("figures")

results_dir <- file.path("results")

signatures_dir <- file.path(data_dir, "signatures")
signatures_subdirs <- file.path(signatures_dir, c("drug", "group", "disease"))


all_dirs <- c(signatures_subdirs, figures_dir, results_dir,
              diseases_filtered_up, diseases_filtered_down,
              diseases_connected_up, diseases_connected_down, diseases_consensus,
              ha1e_filtered_drug_up_dir, ha1e_filtered_drug_down_dir,
              ha1e_filtered_group_up_dir, ha1e_filtered_group_down_dir,
              ha1e_connected_drug_up_dir, ha1e_connected_drug_down_dir,
              ha1e_connected_group_up_dir, ha1e_connected_group_down_dir,
              ha1e_consensus_group_dir, ha1e_consensus_drug_dir,
              mcf7_filtered_drug_up_dir, mcf7_filtered_drug_down_dir,
              mcf7_filtered_group_up_dir, mcf7_filtered_group_down_dir,
              mcf7_connected_drug_up_dir, mcf7_connected_drug_down_dir,
              mcf7_connected_group_up_dir, mcf7_connected_group_down_dir,
              mcf7_consensus_group_dir, mcf7_consensus_drug_dir,
              diseases_consensus_ha1e, diseases_consensus_mcf7)

sapply(all_dirs, create_dir)
