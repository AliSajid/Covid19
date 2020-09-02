# This script runs first and creates the necessary directory structure
#

source("pipeline.R")

data_dir <- file.path("data")
map_dir <- file.path("maps")

ha1e_dir <- file.path(data_dir, "HA1E")
mcf7_dir <- file.path(data_dir, "MCF7")

a5491_dir <- file.path(data_dir, "A549-10uM-24h")
a5492_dir <- file.path(data_dir, "A549-10uM-6h")
ha1e2_dir <- file.path(data_dir, "HA1E-10uM-24h")
ht29_dir <- file.path(data_dir, "HT29-10uM-24h")
mcf72_dir <- file.path(data_dir, "MCF7-10uM-24h")
pc3_dir <- file.path(data_dir, "PC3-10uM-24h")
vcap1_dir <- file.path(data_dir, "VCAP-10uM-24h")
vcap2_dir <- file.path(data_dir, "VCAP-10uM-6h")

ha1e_filtered_dir <- file.path(ha1e_dir, "filtered")
ha1e_filtered_drug_dir <- file.path(ha1e_filtered_dir, "drug")
ha1e_filtered_group_dir <- file.path(ha1e_filtered_dir, "group")
mcf7_filtered_dir <- file.path(mcf7_dir, "filtered")
mcf7_filtered_drug_dir <- file.path(mcf7_filtered_dir, "drug")
mcf7_filtered_group_dir <- file.path(mcf7_filtered_dir, "group")

a5491_filtered_dir <- file.path(a5491_dir, "filtered")
a5491_filtered_drug_dir <- file.path(a5491_filtered_dir, "drug")
a5491_filtered_group_dir <- file.path(a5491_filtered_dir, "group")
a5492_filtered_dir <- file.path(a5492_dir, "filtered")
a5492_filtered_drug_dir <- file.path(a5492_filtered_dir, "drug")
a5492_filtered_group_dir <- file.path(a5492_filtered_dir, "group")
ha1e2_filtered_dir <- file.path(ha1e2_dir, "filtered")
ha1e2_filtered_drug_dir <- file.path(ha1e2_filtered_dir, "drug")
ha1e2_filtered_group_dir <- file.path(ha1e2_filtered_dir, "group")
ht29_filtered_dir <- file.path(ht29_dir, "filtered")
ht29_filtered_drug_dir <- file.path(ht29_filtered_dir, "drug")
ht29_filtered_group_dir <- file.path(ht29_filtered_dir, "group")
mcf72_filtered_dir <- file.path(mcf72_dir, "filtered")
mcf72_filtered_drug_dir <- file.path(mcf72_filtered_dir, "drug")
mcf72_filtered_group_dir <- file.path(mcf72_filtered_dir, "group")
pc3_filtered_dir <- file.path(pc3_dir, "filtered")
pc3_filtered_drug_dir <- file.path(pc3_filtered_dir, "drug")
pc3_filtered_group_dir <- file.path(pc3_filtered_dir, "group")
vcap1_filtered_dir <- file.path(vcap1_dir, "filtered")
vcap1_filtered_drug_dir <- file.path(vcap1_filtered_dir, "drug")
vcap1_filtered_group_dir <- file.path(vcap1_filtered_dir, "group")
vcap2_filtered_dir <- file.path(vcap2_dir, "filtered")
vcap2_filtered_drug_dir <- file.path(vcap2_filtered_dir, "drug")
vcap2_filtered_group_dir <- file.path(vcap2_filtered_dir, "group")

ha1e_filtered_drug_up_dir <- file.path(ha1e_filtered_drug_dir, "up")
ha1e_filtered_drug_down_dir <- file.path(ha1e_filtered_drug_dir, "down")
mcf7_filtered_drug_up_dir <- file.path(mcf7_filtered_drug_dir, "up")
mcf7_filtered_drug_down_dir <- file.path(mcf7_filtered_drug_dir, "down")

a5491_filtered_drug_up_dir <- file.path(a5491_filtered_drug_dir, "up")
a5491_filtered_drug_down_dir <- file.path(a5491_filtered_drug_dir, "down")
a5492_filtered_drug_up_dir <- file.path(a5492_filtered_drug_dir, "up")
a5492_filtered_drug_down_dir <- file.path(a5492_filtered_drug_dir, "down")
ha1e2_filtered_drug_up_dir <- file.path(ha1e2_filtered_drug_dir, "up")
ha1e2_filtered_drug_down_dir <- file.path(ha1e2_filtered_drug_dir, "down")
ht29_filtered_drug_up_dir <- file.path(ht29_filtered_drug_dir, "up")
ht29_filtered_drug_down_dir <- file.path(ht29_filtered_drug_dir, "down")
mcf72_filtered_drug_up_dir <- file.path(mcf72_filtered_drug_dir, "up")
mcf72_filtered_drug_down_dir <- file.path(mcf72_filtered_drug_dir, "down")
pc3_filtered_drug_up_dir <- file.path(pc3_filtered_drug_dir, "up")
pc3_filtered_drug_down_dir <- file.path(pc3_filtered_drug_dir, "down")
vcap1_filtered_drug_up_dir <- file.path(vcap1_filtered_drug_dir, "up")
vcap1_filtered_drug_down_dir <- file.path(vcap1_filtered_drug_dir, "down")
vcap2_filtered_drug_up_dir <- file.path(vcap2_filtered_drug_dir, "up")
vcap2_filtered_drug_down_dir <- file.path(vcap2_filtered_drug_dir, "down")


ha1e_filtered_group_up_dir <- file.path(ha1e_filtered_group_dir, "up")
ha1e_filtered_group_down_dir <- file.path(ha1e_filtered_group_dir, "down")
mcf7_filtered_group_up_dir <- file.path(mcf7_filtered_group_dir, "up")
mcf7_filtered_group_down_dir <- file.path(mcf7_filtered_group_dir, "down")

a5491_filtered_group_up_dir <- file.path(a5491_filtered_group_dir, "up")
a5491_filtered_group_down_dir <- file.path(a5491_filtered_group_dir, "down")
a5492_filtered_group_up_dir <- file.path(a5492_filtered_group_dir, "up")
a5492_filtered_group_down_dir <- file.path(a5492_filtered_group_dir, "down")
ha1e2_filtered_group_up_dir <- file.path(ha1e2_filtered_group_dir, "up")
ha1e2_filtered_group_down_dir <- file.path(ha1e2_filtered_group_dir, "down")
ht29_filtered_group_up_dir <- file.path(ht29_filtered_group_dir, "up")
ht29_filtered_group_down_dir <- file.path(ht29_filtered_group_dir, "down")
mcf72_filtered_group_up_dir <- file.path(mcf72_filtered_group_dir, "up")
mcf72_filtered_group_down_dir <- file.path(mcf72_filtered_group_dir, "down")
pc3_filtered_group_up_dir <- file.path(pc3_filtered_group_dir, "up")
pc3_filtered_group_down_dir <- file.path(pc3_filtered_group_dir, "down")
vcap1_filtered_group_up_dir <- file.path(vcap1_filtered_group_dir, "up")
vcap1_filtered_group_down_dir <- file.path(vcap1_filtered_group_dir, "down")
vcap2_filtered_group_up_dir <- file.path(vcap2_filtered_group_dir, "up")
vcap2_filtered_group_down_dir <- file.path(vcap2_filtered_group_dir, "down")


ha1e_connected_dir <- file.path(ha1e_dir, "connected")
ha1e_connected_drug_dir <- file.path(ha1e_connected_dir, "drug")
ha1e_connected_group_dir <- file.path(ha1e_connected_dir, "group")
mcf7_connected_dir <- file.path(mcf7_dir, "connected")
mcf7_connected_drug_dir <- file.path(mcf7_connected_dir, "drug")
mcf7_connected_group_dir <- file.path(mcf7_connected_dir, "group")

a5491_connected_dir <- file.path(a5491_dir, "connected")
a5491_connected_drug_dir <- file.path(a5491_connected_dir, "drug")
a5491_connected_group_dir <- file.path(a5491_connected_dir, "group")
a5492_connected_dir <- file.path(a5492_dir, "connected")
a5492_connected_drug_dir <- file.path(a5492_connected_dir, "drug")
a5492_connected_group_dir <- file.path(a5492_connected_dir, "group")
ha1e2_connected_dir <- file.path(ha1e2_dir, "connected")
ha1e2_connected_drug_dir <- file.path(ha1e2_connected_dir, "drug")
ha1e2_connected_group_dir <- file.path(ha1e2_connected_dir, "group")
ht29_connected_dir <- file.path(ht29_dir, "connected")
ht29_connected_drug_dir <- file.path(ht29_connected_dir, "drug")
ht29_connected_group_dir <- file.path(ht29_connected_dir, "group")
mcf72_connected_dir <- file.path(mcf72_dir, "connected")
mcf72_connected_drug_dir <- file.path(mcf72_connected_dir, "drug")
mcf72_connected_group_dir <- file.path(mcf72_connected_dir, "group")
pc3_connected_dir <- file.path(pc3_dir, "connected")
pc3_connected_drug_dir <- file.path(pc3_connected_dir, "drug")
pc3_connected_group_dir <- file.path(pc3_connected_dir, "group")
vcap1_connected_dir <- file.path(vcap1_dir, "connected")
vcap1_connected_drug_dir <- file.path(vcap1_connected_dir, "drug")
vcap1_connected_group_dir <- file.path(vcap1_connected_dir, "group")
vcap2_connected_dir <- file.path(vcap2_dir, "connected")
vcap2_connected_drug_dir <- file.path(vcap2_connected_dir, "drug")
vcap2_connected_group_dir <- file.path(vcap2_connected_dir, "group")

ha1e_connected_drug_up_dir <- file.path(ha1e_connected_drug_dir, "up")
ha1e_connected_drug_down_dir <- file.path(ha1e_connected_drug_dir, "down")
mcf7_connected_drug_up_dir <- file.path(mcf7_connected_drug_dir, "up")
mcf7_connected_drug_down_dir <- file.path(mcf7_connected_drug_dir, "down")

a5491_connected_drug_up_dir <- file.path(a5491_connected_drug_dir, "up")
a5491_connected_drug_down_dir <- file.path(a5491_connected_drug_dir, "down")
a5492_connected_drug_up_dir <- file.path(a5492_connected_drug_dir, "up")
a5492_connected_drug_down_dir <- file.path(a5492_connected_drug_dir, "down")
ha1e2_connected_drug_up_dir <- file.path(ha1e2_connected_drug_dir, "up")
ha1e2_connected_drug_down_dir <- file.path(ha1e2_connected_drug_dir, "down")
ht29_connected_drug_up_dir <- file.path(ht29_connected_drug_dir, "up")
ht29_connected_drug_down_dir <- file.path(ht29_connected_drug_dir, "down")
mcf72_connected_drug_up_dir <- file.path(mcf72_connected_drug_dir, "up")
mcf72_connected_drug_down_dir <- file.path(mcf72_connected_drug_dir, "down")
pc3_connected_drug_up_dir <- file.path(pc3_connected_drug_dir, "up")
pc3_connected_drug_down_dir <- file.path(pc3_connected_drug_dir, "down")
vcap1_connected_drug_up_dir <- file.path(vcap1_connected_drug_dir, "up")
vcap1_connected_drug_down_dir <- file.path(vcap1_connected_drug_dir, "down")
vcap2_connected_drug_up_dir <- file.path(vcap2_connected_drug_dir, "up")
vcap2_connected_drug_down_dir <- file.path(vcap2_connected_drug_dir, "down")

ha1e_connected_group_up_dir <- file.path(ha1e_connected_group_dir, "up")
ha1e_connected_group_down_dir <- file.path(ha1e_connected_group_dir, "down")
mcf7_connected_group_up_dir <- file.path(mcf7_connected_group_dir, "up")
mcf7_connected_group_down_dir <- file.path(mcf7_connected_group_dir, "down")

a5491_connected_group_up_dir <- file.path(a5491_connected_group_dir, "up")
a5491_connected_group_down_dir <- file.path(a5491_connected_group_dir, "down")
a5492_connected_group_up_dir <- file.path(a5492_connected_group_dir, "up")
a5492_connected_group_down_dir <- file.path(a5492_connected_group_dir, "down")
ha1e2_connected_group_up_dir <- file.path(ha1e2_connected_group_dir, "up")
ha1e2_connected_group_down_dir <- file.path(ha1e2_connected_group_dir, "down")
ht29_connected_group_up_dir <- file.path(ht29_connected_group_dir, "up")
ht29_connected_group_down_dir <- file.path(ht29_connected_group_dir, "down")
mcf72_connected_group_up_dir <- file.path(mcf72_connected_group_dir, "up")
mcf72_connected_group_down_dir <- file.path(mcf72_connected_group_dir, "down")
pc3_connected_group_up_dir <- file.path(pc3_connected_group_dir, "up")
pc3_connected_group_down_dir <- file.path(pc3_connected_group_dir, "down")
vcap1_connected_group_up_dir <- file.path(vcap1_connected_group_dir, "up")
vcap1_connected_group_down_dir <- file.path(vcap1_connected_group_dir, "down")
vcap2_connected_group_up_dir <- file.path(vcap2_connected_group_dir, "up")
vcap2_connected_group_down_dir <- file.path(vcap2_connected_group_dir, "down")

ha1e_consensus_dir <- file.path(ha1e_dir, "consensus")
ha1e_consensus_drug_dir <- file.path(ha1e_consensus_dir, "drug")
ha1e_consensus_group_dir <- file.path(ha1e_consensus_dir, "group")
mcf7_consensus_dir <- file.path(mcf7_dir, "consensus")
mcf7_consensus_drug_dir <- file.path(mcf7_consensus_dir, "drug")
mcf7_consensus_group_dir <- file.path(mcf7_consensus_dir, "group")

a5491_consensus_dir <- file.path(a5491_dir, "consensus")
a5491_consensus_drug_dir <- file.path(a5491_consensus_dir, "drug")
a5491_consensus_group_dir <- file.path(a5491_consensus_dir, "group")
a5492_consensus_dir <- file.path(a5492_dir, "consensus")
a5492_consensus_drug_dir <- file.path(a5492_consensus_dir, "drug")
a5492_consensus_group_dir <- file.path(a5492_consensus_dir, "group")
ha1e2_consensus_dir <- file.path(ha1e2_dir, "consensus")
ha1e2_consensus_drug_dir <- file.path(ha1e2_consensus_dir, "drug")
ha1e2_consensus_group_dir <- file.path(ha1e2_consensus_dir, "group")
ht29_consensus_dir <- file.path(ht29_dir, "consensus")
ht29_consensus_drug_dir <- file.path(ht29_consensus_dir, "drug")
ht29_consensus_group_dir <- file.path(ht29_consensus_dir, "group")
mcf72_consensus_dir <- file.path(mcf72_dir, "consensus")
mcf72_consensus_drug_dir <- file.path(mcf72_consensus_dir, "drug")
mcf72_consensus_group_dir <- file.path(mcf72_consensus_dir, "group")
pc3_consensus_dir <- file.path(pc3_dir, "consensus")
pc3_consensus_drug_dir <- file.path(pc3_consensus_dir, "drug")
pc3_consensus_group_dir <- file.path(pc3_consensus_dir, "group")
vcap1_consensus_dir <- file.path(vcap1_dir, "consensus")
vcap1_consensus_drug_dir <- file.path(vcap1_consensus_dir, "drug")
vcap1_consensus_group_dir <- file.path(vcap1_consensus_dir, "group")
vcap2_consensus_dir <- file.path(vcap2_dir, "consensus")
vcap2_consensus_drug_dir <- file.path(vcap2_consensus_dir, "drug")
vcap2_consensus_group_dir <- file.path(vcap2_consensus_dir, "group")


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

diseases_consensus_a5491 <- file.path(diseases_consensus, "A549-10uM-24h")
diseases_consensus_a5492 <- file.path(diseases_consensus, "A549-10uM-6h")
diseases_consensus_ha1e2 <- file.path(diseases_consensus, "HA1E-10uM-24h")
diseases_consensus_ht29 <- file.path(diseases_consensus, "HT29-10uM-24h")
diseases_consensus_mcf72 <- file.path(diseases_consensus, "MCF7-10uM-24h")
diseases_consensus_pc3 <- file.path(diseases_consensus, "PC3-10uM-24h")
diseases_consensus_vcap1 <- file.path(diseases_consensus, "VCAP-10uM-24h")
diseases_consensus_vcap2 <- file.path(diseases_consensus, "VCAP-10uM-6h")
diseases_consensus_all <- file.path(diseases_consensus, "all")

figures_dir <- file.path("figures")

results_dir <- file.path("results")

signatures_dir <- file.path(data_dir, "signatures")
signatures_subdirs <- file.path(signatures_dir, c("drug", "group", "disease", "sig_lists"))


all_dirs <- c(map_dir, signatures_subdirs, figures_dir, results_dir,
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
              diseases_consensus_ha1e, diseases_consensus_mcf7,
              diseases_consensus_a5491, diseases_consensus_a5492,
              diseases_consensus_ha1e2, diseases_consensus_ht29,
              diseases_consensus_mcf72, diseases_consensus_pc3,
              diseases_consensus_vcap1, diseases_consensus_vcap2,
              a5491_filtered_drug_up_dir, a5491_filtered_drug_down_dir,
              a5491_filtered_group_up_dir, a5491_filtered_group_down_dir,
              a5491_connected_drug_up_dir, a5491_connected_drug_down_dir,
              a5491_connected_group_up_dir, a5491_connected_group_down_dir,
              a5491_consensus_group_dir, a5491_consensus_drug_dir,
              a5492_filtered_drug_up_dir, a5492_filtered_drug_down_dir,
              a5492_filtered_group_up_dir, a5492_filtered_group_down_dir,
              a5492_connected_drug_up_dir, a5492_connected_drug_down_dir,
              a5492_connected_group_up_dir, a5492_connected_group_down_dir,
              a5492_consensus_group_dir, a5492_consensus_drug_dir,
              ha1e2_filtered_drug_up_dir, ha1e2_filtered_drug_down_dir,
              ha1e2_filtered_group_up_dir, ha1e2_filtered_group_down_dir,
              ha1e2_connected_drug_up_dir, ha1e2_connected_drug_down_dir,
              ha1e2_connected_group_up_dir, ha1e2_connected_group_down_dir,
              ha1e2_consensus_group_dir, ha1e2_consensus_drug_dir,
              ht29_filtered_drug_up_dir, ht29_filtered_drug_down_dir,
              ht29_filtered_group_up_dir, ht29_filtered_group_down_dir,
              ht29_connected_drug_up_dir, ht29_connected_drug_down_dir,
              ht29_connected_group_up_dir, ht29_connected_group_down_dir,
              ht29_consensus_group_dir, ht29_consensus_drug_dir,
              mcf72_filtered_drug_up_dir, mcf72_filtered_drug_down_dir,
              mcf72_filtered_group_up_dir, mcf72_filtered_group_down_dir,
              mcf72_connected_drug_up_dir, mcf72_connected_drug_down_dir,
              mcf72_connected_group_up_dir, mcf72_connected_group_down_dir,
              mcf72_consensus_group_dir, mcf72_consensus_drug_dir,
              pc3_filtered_drug_up_dir, pc3_filtered_drug_down_dir,
              pc3_filtered_group_up_dir, pc3_filtered_group_down_dir,
              pc3_connected_drug_up_dir, pc3_connected_drug_down_dir,
              pc3_connected_group_up_dir, pc3_connected_group_down_dir,
              pc3_consensus_group_dir, pc3_consensus_drug_dir,
              vcap1_filtered_drug_up_dir, vcap1_filtered_drug_down_dir,
              vcap1_filtered_group_up_dir, vcap1_filtered_group_down_dir,
              vcap1_connected_drug_up_dir, vcap1_connected_drug_down_dir,
              vcap1_connected_group_up_dir, vcap1_connected_group_down_dir,
              vcap1_consensus_group_dir, vcap1_consensus_drug_dir,
              vcap2_filtered_drug_up_dir, vcap2_filtered_drug_down_dir,
              vcap2_filtered_group_up_dir, vcap2_filtered_group_down_dir,
              vcap2_connected_drug_up_dir, vcap2_connected_drug_down_dir,
              vcap2_connected_group_up_dir, vcap2_connected_group_down_dir,
              vcap2_consensus_group_dir, vcap2_consensus_drug_dir, diseases_consensus_all
              )

sapply(all_dirs, create_dir)
