library(tidyverse)
source("pipeline.R")

generate_common_perturbagens <-
  function(disease, disease_thresh, groups_thresh, cell_line) {

    #    group_names <- c("azt", "los", "nib", "quine", "vir")

    disease_prefix <-
      file.path("data", "disease", disease, "consensus", "all")
    disease_filename <-
      paste(disease, disease_thresh, "consensus", sep = "-")
    disease_file <- generate_name(disease_prefix, disease_filename, "tsv")

    data_prefix <- file.path("data", cell_line, "consensus", "group")
    data_files <-
      list.files(data_prefix, pattern = as.character(groups_thresh), full.names = T)

    outfile_prefix <- "results"
    outfile_name <- paste(cell_line, disease_thresh, groups_thresh, "common", "filtered", sep = "-")
    outfile <- generate_name(outfile_prefix, outfile_name, "csv")

    col_spec <- cols(
      .default = col_double(),
      signatureid = col_skip(),
      compound = col_character()
    )

    group_names <- str_extract(basename(data_files), "(\\w+)")

    disease_data <- read_tsv(disease_file, col_types = col_spec)

    group_data <-
      map2(
        data_files, group_names, function(x, y)
          read_tsv(
            x,
            col_types = col_spec,
            col_names =  c("signatureid", "compound", y),
            skip = 1
          ))

    mean_data <- reduce(group_data, full_join, by = "compound") %>%
      inner_join(disease_data, by = "compound") %>%
      pivot_longer(cols = any_of(group_names)) %>%
      group_by(compound) %>%
      summarise(meansim = mean(value, na.rm = T))

    complete_data <- reduce(group_data, full_join, by = "compound") %>%
      inner_join(disease_data, by = "compound") %>%
      inner_join(mean_data, by = "compound") %>%
      mutate(nas = rowSums(is.na(.))) %>%
      filter(nas < (length(group_names) - 1)) %>%
      write_csv(outfile)

    return(complete_data)
  }

a5491_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "A549-10um-24h")
a5492_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "A549-10um-6h")
ha1e2_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "HA1E-10um-24h")
ht29_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "HT29-10um-24h")
mcf72_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "MCF7-10um-24h")
pc3_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "PC3-10um-24h")
vcap1_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "VCAP-10um-24h")
vcap2_all <- generate_common_perturbagens("covidm", 0.5, 0.85, "VCAP-10um-6h")

a5491 <- a5491_all %>% select(compound, meansim) %>% rename(a5491 = meansim)
a5492 <- a5492_all %>% select(compound, meansim) %>% rename(a5492 = meansim)
ha1e2 <- ha1e2_all %>% select(compound, meansim) %>% rename(ha1e2 = meansim)
ht29 <- ht29_all %>% select(compound, meansim) %>% rename(ht29 = meansim)
mcf72 <- mcf72_all %>% select(compound, meansim) %>% rename(mcf72 = meansim)
pc3 <- pc3_all %>% select(compound, meansim) %>% rename(pc3 = meansim)
vcap1 <- vcap1_all %>% select(compound, meansim) %>% rename(vcap1 = meansim)
vcap2 <- vcap2_all %>% select(compound, meansim) %>% rename(vcap2 = meansim)

all_data <- list(a5491, a5492, ha1e2, ht29, mcf72, pc3, vcap1, vcap2)

combined <- reduce(all_data, inner_join, by = "compound") %>%
  unique %>%
  write_csv("results/covidm-new_groups-combined-drugs.csv")
