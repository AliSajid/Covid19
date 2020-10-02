# This script creates a file containing the skew, and other statistics for normality for each drug signatures


library(tidyverse)
library(broom)
library(e1071)

files <- list.files("data/signatures/group/", full.names = T)

col_spec <- cols(
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double()
)

output_statistics <- function(df, name) {
  pat <- "([a-z]+)-(\\w+-\\w+-\\w+)-signature.tsv"
  matches <- str_match(name, pat)

  group <- matches[,2]
  cell_line <- matches[,3]
  normal <- shapiro.test(df$Value_LogDiffExp)
  skew <- skewness(df$Value_LogDiffExp)
  out <- list(w = round(normal$statistic, 3),
              pval = round(normal$p.value, 3),
              line = cell_line,
              skewness = round(skew, 3))
}

output <- files %>%
  map(~ read_tsv(.x, col_types = col_spec)) %>%
  map2_dfr(basename(files), ~ output_statistics(.x, .y)) %>%
  select(line, skewness, w, pval) %>%
  rename(CellLine = line,
         Skewness = skewness,
         W = w,
         `P-Value` = pval) %>%
  write_csv("results/cell-line_normality_statistics.csv")
