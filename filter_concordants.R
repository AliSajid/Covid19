# Obsoleted

library(tidyverse)

filter_data <-
  function(dataframe_up,
           dataframe_down,
           cell_line,
           cutoff) {
    dataframe <- bind_rows(dataframe_up, dataframe_down)
    output <- dataframe %>%
      filter(CellLine == cell_line, Concordance > cutoff) %>%
      group_by(Perturbagen) %>%
      filter(Concordance == max(Concordance)) %>%
      ungroup() %>%
      select(SignatureId, Perturbagen, Concordance)
    return(output)
  }

drugs_mcf7 <-
  c(
    "Azithromycin",
    "Baricitinib",
    "Chloroquine",
    "Fedratinib",
    "Lopinavir",
    "Losartan",
    "Ritonavir",
    "Ruxolitinib",
    "Hydroxychloroquine",
    "Group-Azt",
    "Group-L01XE",
    "Group-Los",
    "Group-Nib",
    "Group-Quine",
    "Group-Vir"
  )

drugs_ha1e <-
  c(
    "Azithromycin",
    "Baricitinib",
    "Chloroquine",
    "Fedratinib",
    "Lopinavir",
    "Losartan",
    "Ritonavir",
    "Ruxolitinib",
    "Group-Azt",
    "Group-L01XE",
    "Group-Los",
    "Group-Nib",
    "Group-Quine",
    "Group-Vir"
  )

process_cell_line <-
  function(cell_line, drugs_list, cutoff = 0.321) {
    prefix <-
      paste("data",
            "L1000-Signature",
            "connected_signatures",
            cell_line,
            sep = "/")
    for (drug in drugs_list) {
      upfile <-
        paste(paste(drug, cell_line, "Up", "Connected", sep = "_"),
              "tsv",
              sep = ".")
      downfile <-
        paste(paste(drug, cell_line, "Down", "Connected", sep = "_"),
              "tsv",
              sep = ".")
      up <- read_tsv(paste(prefix, upfile, sep = "/"))
      down <- read_tsv(paste(prefix, downfile, sep = "/"))
      consensus <- filter_data(up, down, cell_line, cutoff)
      write_tsv(consensus, path = paste(prefix,
                                        paste(
                                          paste(drug,
                                                cell_line,
                                                "Consensus",
                                                "Connected",
                                                sep = "_"),
                                          "tsv",
                                          sep = "."
                                        ),
                                        sep = "/"))
      write_csv(consensus, path = paste(prefix, paste(
        paste(drug, cell_line, "Consensus", "Connected", sep = "_"),
        "csv",
        sep = "."
      ), sep = "/"))
    }
  }

process_cell_line("MCF7", drugs_mcf7)
process_cell_line("HA1E", drugs_ha1e)