library(tidyverse)

col_names <- c("Perturbagen", "Synonym", "CellLine", "signatureID",
               "PROBE", "ID_geneid", "Name_GeneSymbol", "Value_LogDiffExp",
               "Significance_pvalue")
col_spec <- cols(
  .default = col_character(),
  PROBE = col_skip(),
  Value_LogDiffExp = col_double(),
  Significance_pvalue = col_double()
)

df <- read_csv("data/Los-Azt-All-Signature.csv", col_names = col_names, col_types = col_spec, skip = 1)

azithromycin_HA1E_up <- df %>% 
  filter(Value_LogDiffExp > 0.85 , Synonym == "Azithromycin", CellLine == "HA1E") %>% 
  select(6:8)

azithromycin_MCF7_up <- df %>% 
  filter(Value_LogDiffExp > 0.85 , Synonym == "Azithromycin", CellLine == "MCF7") %>% 
  select(6:8)

losartan_HA1E_up <- df %>% 
  filter(Value_LogDiffExp > 0.85 , Synonym == "Losartan", CellLine == "HA1E") %>% 
  select(6:8)

losartan_MCF7_up <- df %>% 
  filter(Value_LogDiffExp > 0.85 , Synonym == "Losartan", CellLine == "MCF7") %>% 
  select(6:8)

azithromycin_HA1E_down <- df %>% 
  filter(Value_LogDiffExp < -0.85, Synonym == "Azithromycin", CellLine == "HA1E") %>% 
  select(6:8)

azithromycin_MCF7_down <- df %>% 
  filter(Value_LogDiffExp < -0.85, Synonym == "Azithromycin", CellLine == "MCF7") %>% 
  select(6:8)

losartan_HA1E_down <- df %>% 
  filter(Value_LogDiffExp < -0.85, Synonym == "Losartan", CellLine == "HA1E") %>% 
  select(6:8)

losartan_MCF7_down <- df %>% 
  filter(Value_LogDiffExp < -0.85, Synonym == "Losartan", CellLine == "MCF7") %>% 
  select(6:8)


write_csv(path = "data/Azithromycin-HA1E-Signature-Up.csv", x = azithromycin_HA1E_up)
write_csv(path = "data/Azithromycin-MCF7-Signature-Up.csv", x = azithromycin_MCF7_up)
write_csv(path = "data/Azithromycin-HA1E-Signature-Down.csv", x = azithromycin_HA1E_down)
write_csv(path = "data/Azithromycin-MCF7-Signature-Down.csv", x = azithromycin_MCF7_down)
write_csv(path = "data/Losartan-HA1E-Signature-Up.csv", x = losartan_HA1E_up)
write_csv(path = "data/Losartan-MCF7-Signature-Up.csv", x = losartan_MCF7_up)
write_csv(path = "data/Losartan-HA1E-Signature-Down.csv", x = losartan_HA1E_down)
write_csv(path = "data/Losartan-MCF7-Signature-Down.csv", x = losartan_MCF7_down)

write_tsv(path = "data/Azithromycin-HA1E-Signature-Up.tsv", x = azithromycin_HA1E_up)
write_tsv(path = "data/Azithromycin-MCF7-Signature-Up.tsv", x = azithromycin_MCF7_up)
write_tsv(path = "data/Azithromycin-HA1E-Signature-Down.tsv", x = azithromycin_HA1E_down)
write_tsv(path = "data/Azithromycin-MCF7-Signature-Down.tsv", x = azithromycin_MCF7_down)
write_tsv(path = "data/Losartan-HA1E-Signature-Up.tsv", x = losartan_HA1E_up)
write_tsv(path = "data/Losartan-MCF7-Signature-Up.tsv", x = losartan_MCF7_up)
write_tsv(path = "data/Losartan-HA1E-Signature-Down.tsv", x = losartan_HA1E_down)
write_tsv(path = "data/Losartan-MCF7-Signature-Down.tsv", x = losartan_MCF7_down)
