library(tidyverse)
source("pipeline.R")

generate_common_perturbagens <-
  function(sars_thresh, groups_thresh, cell_line) {

    group_names <- c("azt", "los", "nib", "quine", "vir")

    sars_prefix <-
      file.path("data", "disease", "sars", "consensus", cell_line)
    sars_filename <-
      paste("sars", sars_thresh, "consensus", sep = "-")
    sars_file <- generate_name(sars_prefix, sars_filename, "tsv")

    data_prefix <- file.path("data", cell_line, "consensus", "group")
    data_files <-
      list.files(data_prefix, pattern = as.character(groups_thresh), full.names = T)

    outfile_prefix <- "results"
    outfile_name <- paste(cell_line, sars_thresh, groups_thresh, "common", "filtered", sep = "-")
    outfile <- generate_name(outfile_prefix, outfile_name, "csv")

    col_spec <- cols(
      .default = col_double(),
      signatureid = col_skip(),
      compound = col_character()
    )

    sars_data <- read_tsv(sars_file, col_types = col_spec)

    group_data <-
      map2(
        data_files, group_names, function(x, y)
        read_tsv(
          x,
          col_types = col_spec,
          col_names =  c("signatureid", "compound", y),
          skip = 1
        ))

    complete_data <- reduce(group_data, full_join, by = "compound") %>%
      inner_join(sars_data, by = "compound") %>%
      mutate(nas = rowSums(is.na(.)), presence = 5 - nas) %>%
      filter(nas < 4) %>%
      rowwise %>%
      mutate(meansim = mean(c(azt, los, nib, quine, vir), na.rm = T)) %>%
      write_csv(outfile)

    return(complete_data)
  }

mcf7 <- generate_common_perturbagens(0.5, 0.85, "MCF7") %>% select(compound, meansim) %>% rename(mcf7 = meansim)
ha1e <- generate_common_perturbagens(0.5, 0.85, "HA1E") %>% select(compound, meansim) %>% rename(ha1e = meansim)

combined <- inner_join(mcf7, ha1e, by = "compound") %>% write_csv("results/ha1e-mcf7-combined-drugs.csv")
