library(tidyverse)

groups_l01 <- c("Group-Azt",
            "Group-L01XE",
            "Group-Los",
            "Group-Nib",
            "Group-Quine",
            "Group-Vir")

groups <- c("Group-Azt",
                "Group-Los",
                "Group-Nib",
                "Group-Quine",
                "Group-Vir")

cell_lines <- c("HA1E", "MCF7")

suffix_cons <- "Consensus"
suffix_conn <- "Connected"

extensions <- "csv"

prefix <-
  paste("data",
        "L1000-Signature",
        "connected_signatures",
        sep = "/")

f <-
  expand.grid(prefix, cell_lines, groups, cell_lines, suffix_cons, suffix_conn, extensions)

f <- f %>% 
  filter(Var2 == Var4) %>% 
  unite(name, Var3:Var6) %>% 
  unite(name, name, Var7, sep = ".") %>% 
  unite(name, sep = "/")

files <- f$name

f_l01 <-
  expand.grid(prefix, cell_lines, groups_l01, cell_lines, suffix_cons, suffix_conn, extensions)

f_l01 <- f_l01 %>% 
  filter(Var2 == Var4) %>% 
  unite(name, Var3:Var6) %>% 
  unite(name, name, Var7, sep = ".") %>% 
  unite(name, sep = "/")

files_l01 <- f_l01$name

return_perturbagens <- function(files, cell_line) {
  selected_files <- files[str_detect(files, cell_line)]
  pattern = "(Group-[A-Za-z0-9]+)"
  group_names <- str_extract(selected_files, pattern)
  
  for (index in 1:length(selected_files)) {
    data <- read_csv(files[index])$Perturbagen
    assign(group_names[index], data)
  }
  
  groups <- sapply(group_names, function(x) {return(get(x, envir = parent.frame(n=3)))})
  return(groups)
}

data_mcf7 <- return_perturbagens(files, "MCF7")
data_ha1e <- return_perturbagens(files, "HA1E")

data_l01_mcf7 <- return_perturbagens(files_l01, "MCF7")
data_l01_ha1e <- return_perturbagens(files_l01, "HA1E")