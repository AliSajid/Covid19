library(tidyverse)
source("pipeline.R")

cutoffs <- c(0, 0.26, 0.5)

diseases <- c("influenza", "mers", "sars")

cell_lines <- c("HA1E", "MCF7")

prefix <- "dataset/disease"
HA1E_prefix <- paste(prefix, "HA1E", sep = "/")
MCF7_prefix <- paste(prefix, "MCF7", sep = "/")
up_prefix <- paste(prefix, diseases, "up", sep = "/")
down_prefix <- paste(prefix, diseases, "down", sep = "/")

if (!dir.exists(HA1E_prefix)) {
  dir.create(HA1E_prefix)
}

if (!dir.exists(MCF7_prefix)) {
  dir.create(MCF7_prefix)
}

for (pre in c(down_prefix, up_prefix)) {
  if (!dir.exists(pre)) {
    dir.create(pre)
  }
}

# Process HA1E Cell Line
for (index in 1:length(diseases)) {
  for (cutoff in cutoffs) {
    dis <- diseases[index]
    infile <- generate_name(paste(prefix, dis, sep = "/"),
                            paste(dis, "l1000", "signature", sep = "-"),
                            "tsv")
    outfile_up <- generate_name(up_prefix[index],
                                paste(dis, "l1000", "up", cutoff, "signature", sep = "-"),
                                "tsv")
    outfile_down <- generate_name(down_prefix[index],
                                paste(dis, "l1000", "down", cutoff, "signature", sep = "-"),
                                "tsv")
    outfile_concordant_up <- generate_name(up_prefix[index],
                                           paste(dis, "l1000", "up", cutoff, "connected", sep = "-"),
                                           "tsv")
    outfile_concordant_down <- generate_name(down_prefix[index],
                                           paste(dis, "l1000", "down", cutoff, "connected", sep = "-"),
                                           "tsv")
    outfile_consensus <- generate_name(HA1E_prefix,
                                       paste(dis, "l1000", "HA1E", cutoff, "consensus", "connected", sep = "-"),
                                       "tsv")
    
    df <- read_tsv(infile)
    df_up <- generate_filtered_signature(df, direction = "up", threshold = cutoff)
    write_tsv(df_up, outfile_up)
    df_down <- generate_filtered_signature(df, direction = "down", threshold = cutoff)
    write_tsv(df_down, outfile_down)
    
    up_concordant <- get_concordant_signatures(signature_df = df_up)
    write_tsv(up_concordant, outfile_concordant_up)
    down_concordant <- get_concordant_signatures(signature_df = df_down)
    write_tsv(down_concordant, outfile_concordant_down)
    
    consensus <- generate_consensus_signature(up_concordant, down_concordant,
                                              cell_line = "HA1E", discordant = TRUE)
    write_tsv(consensus, outfile_consensus)
  }
}

# Process MCF7 Cell Line
for (index in 1:length(diseases)) {
  for (cutoff in cutoffs) {
    dis <- diseases[index]
    infile <- generate_name(paste(prefix, dis, sep = "/"),
                            paste(dis, "l1000", "signature", sep = "-"),
                            "tsv")
    outfile_up <- generate_name(up_prefix[index],
                                paste(dis, "l1000", "up", cutoff, "signature", sep = "-"),
                                "tsv")
    outfile_down <- generate_name(down_prefix[index],
                                  paste(dis, "l1000", "down", cutoff, "signature", sep = "-"),
                                  "tsv")
    outfile_concordant_up <- generate_name(up_prefix[index],
                                           paste(dis, "l1000", "up", cutoff, "connected", sep = "-"),
                                           "tsv")
    outfile_concordant_down <- generate_name(down_prefix[index],
                                             paste(dis, "l1000", "up", cutoff, "connected", sep = "-"),
                                             "tsv")
    outfile_consensus <- generate_name(MCF7_prefix,
                                       paste(dis, "l1000", "MCF7", cutoff, "consensus", "connected", sep = "-"),
                                       "tsv")
    
    df <- read_tsv(infile)
    df_up <- generate_filtered_signature(df, direction = "up", threshold = cutoff)
    write_tsv(df_up, outfile_up)
    df_down <- generate_filtered_signature(df, direction = "down", threshold = cutoff)
    write_tsv(df_down, outfile_down)
    
    up_concordant <- get_concordant_signatures(signature_df = df_up)
    write_tsv(up_concordant, outfile_concordant_up)
    down_concordant <- get_concordant_signatures(signature_df = df_down)
    write_tsv(down_concordant, outfile_concordant_down)
    
    consensus <- generate_consensus_signature(up_concordant, down_concordant,
                                              cell_line = "MCF7", discordant = TRUE)
    write_tsv(consensus, outfile_consensus)
  }
}