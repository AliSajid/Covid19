library(tidyverse)

diseases <- c("Influenza", "SARS", "MERS")

cutoffs <- c("0", "0.5", "0.26")

cell_lines <- c("HA1E", "MCF7")

suffix_cons <- "Consensus"
suffix_conn <- "Connected"

extensions <- "csv"

prefix <-
  paste("data",
        "Corona-Signature",
        "connected_signatures",
        sep = "/")

f <-
  expand.grid(prefix, diseases, cutoffs, cell_lines, suffix_cons, suffix_conn, extensions)

f <- f %>% 
  unite(name, Var2:Var6) %>% 
  unite(name, name, Var7, sep = ".") %>% 
  unite(name, sep = "/")

files <- f$name

return_perturbagens <- function(files, cell_line, cutoff) {
  selected_files <- files[str_detect(files, cell_line)]
  selected_files <- selected_files[str_detect(selected_files, paste0(cutoff, "_"))]
  
  pattern = "connected_signatures/([A-Za-z0-9]+)"
  group_names <- str_match(selected_files, pattern)
  
  for (index in 1:length(selected_files)) {
    data <- read_csv(files[index])$Perturbagen
    assign(group_names[,2][index], data)
  }
  
  groups <- sapply(group_names, function(x) {return(get(x, envir = parent.frame(n=3)))})
  return(groups)
}

data_mcf7_0 <- return_perturbagens(files, "MCF7", "0")
data_ha1e_0 <- return_perturbagens(files, "HA1E", "0")

data_mcf7_0.5 <- return_perturbagens(files, "MCF7", "0.5")
data_ha1e_0.5 <- return_perturbagens(files, "HA1E", "0.5")

data_mcf7_0.26 <- return_perturbagens(files, "MCF7", "0.26")
data_ha1e_0.26 <- return_perturbagens(files, "HA1E", "0.26")
