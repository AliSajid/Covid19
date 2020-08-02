library(tidyverse)
source("pipeline.R")

prefix <- "dataset/drugs"
raw_prefix <- paste(prefix, c("HA1E", "MCF7"), "raw", sep = "/")
up_prefix <- paste(prefix, c("HA1E", "MCF7"), "up", sep = "/")
down_prefix <- paste(prefix, c("HA1E", "MCF7"), "down", sep = "/")
concordant_prefix <- paste(prefix, c("HA1E", "MCF7"), "concordant", sep = "/")
consensus_prefix <- paste(prefix, c("HA1E", "MCF7"), "consensus", sep = "/")

HA1E_signatures <- read_tsv("HA1E-Drug-Signature_Map.tsv") %>% select(SignatureId)
MCF7_signatures <- read_tsv("MCF7-Drug-Signature_Map.tsv") %>% select(SignatureId)


# Process Signatures for HA1E Cell Lines
if (!dir.exists(raw_prefix[1])) {
  dir.create(raw_prefix[1], recursive = T)
}

for (index in 1:dim(HA1E_signatures)[1]) {
  sig <- HA1E_signatures[index,]
  prefix <- paste(prefix, "HA1E", sep = "/")
  output_file_name <- generate_name(raw_prefix[1], paste(sig, "l1000", "signature", sep = "-"), "tsv")
  d <- get_l1000_signature(sig)
  write_tsv(d, output_file_name)
}

# Process Signatures for MCF7 Cell Lines
if (!dir.exists(raw_prefix[2])) {
  dir.create(raw_prefix[2], recursive = T)
}

for (index in 1:dim(MCF7_signatures)[1]) {
  sig <- MCF7_signatures[index,]
  output_file_name <- generate_name(raw_prefix[2], paste(sig, "l1000", "signature", sep = "-"), "tsv")
  d <- get_l1000_signature(sig)
  write_tsv(d, output_file_name)
}

# Process Signatures for HA1E Cell Lines
if (!dir.exists(up_prefix[1])) {
  dir.create(up_prefix[1])
}

if (!dir.exists(down_prefix[1])) {
  dir.create(down_prefix[1])
}

if (!dir.exists(concordant_prefix[1])) {
  dir.create(concordant_prefix[1])
}

if (!dir.exists(consensus_prefix[1])) {
  dir.create(consensus_prefix[1])
}

for (index in 1:dim(HA1E_signatures)[1]) {
  sig <- HA1E_signatures[index,]
  infile <- generate_name(raw_prefix[1], paste(sig, "l1000", "signature", sep = "-"), "tsv")
  outfile_up <- generate_name(up_prefix[1], paste(sig, "l1000", "up", "signature", sep = "-"), "tsv")
  outfile_down <- generate_name(down_prefix[1], paste(sig, "l1000", "down", "signature", sep = "-"), "tsv")
  outfile_concordant_up <- generate_name(concordant_prefix[1], 
                                      paste(sig, "l1000", "concordant", "up", "connected", sep = "-"),
                                      "tsv")
  outfile_concordant_down <- generate_name(concordant_prefix[1], 
                                         paste(sig, "l1000", "concordant", "down", "connected", sep = "-"),
                                         "tsv")
  outfile_consensus <- generate_name(consensus_prefix[1], 
                                      paste(sig, "l1000", "consensus","connected", sep = "-"),
                                      "tsv")
  
  df <- read_tsv(infile)
  df_up <- generate_filtered_signature(df, direction = "up")
  write_tsv(df_up, outfile_up)
  df_down <- generate_filtered_signature(df, direction = "down")
  write_tsv(df_down, outfile_down)
  
  up_concordant <- get_concordant_signatures(signature_df = df_up)
  down_concordant <- get_concordant_signatures(signature_df = df_down)
  write_tsv(up_concordant, outfile_concordant_up)
  write_tsv(down_concordant, outfile_concordant_down)
  
  consensus <- generate_consensus_signature(up_concordant, down_concordant, cell_line = "HA1E")
  write_tsv(consensus, outfile_consensus)
}

# Process Signatures for MCF7 Cell Lines
if (!dir.exists(up_prefix[2])) {
  dir.create(up_prefix[2])
}

if (!dir.exists(down_prefix[2])) {
  dir.create(down_prefix[2])
}

if (!dir.exists(concordant_prefix[2])) {
  dir.create(concordant_prefix[2])
}

if (!dir.exists(consensus_prefix[2])) {
  dir.create(consensus_prefix[2])
}

for (index in 1:dim(MCF7_signatures)[1]) {
  sig <- MCF7_signatures[index,]
  infile <- generate_name(raw_prefix[2], paste(sig, "l1000", "signature", sep = "-"), "tsv")
  outfile_up <- generate_name(up_prefix[2], paste(sig, "l1000", "up", "signature", sep = "-"), "tsv")
  outfile_down <- generate_name(down_prefix[2], paste(sig, "l1000", "down", "signature", sep = "-"), "tsv")
  outfile_concordant_up <- generate_name(concordant_prefix[2], 
                                         paste(sig, "l1000", "concordant", "up", "connected", sep = "-"),
                                         "tsv")
  outfile_concordant_down <- generate_name(concordant_prefix[2], 
                                           paste(sig, "l1000", "concordant", "down", "connected", sep = "-"),
                                           "tsv")
  outfile_consensus <- generate_name(consensus_prefix[2], 
                                     paste(sig, "l1000", "consensus","connected", sep = "-"),
                                     "tsv")
  
  df <- read_tsv(infile)
  df_up <- generate_filtered_signature(df, direction = "up")
  write_tsv(df_up, outfile_up)
  df_down <- generate_filtered_signature(df, direction = "down")
  write_tsv(df_down, outfile_down)
  
  up_concordant <- get_concordant_signatures(signature_df = df_up)
  down_concordant <- get_concordant_signatures(signature_df = df_down)
  write_tsv(up_concordant, outfile_concordant_up)
  write_tsv(down_concordant, outfile_concordant_down)
  
  consensus <- generate_consensus_signature(up_concordant, down_concordant, cell_line = "MCF7")
  write_tsv(consensus, outfile_consensus)
}
