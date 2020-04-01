library(httr)
library(jsonlite)
library(tidyverse)

get_l1000_signature <- function(signature_id) {
  url <- "http://www.ilincs.org/api/ilincsR/downloadSignature"
  query = list(sigID = signature_id, noOfTopGenes = 978)
  
  request <- POST(url, query = query)
  
  if (status_code(request) == 200) {
    signature <- fromJSON(
      toJSON(
        content(request)$data$signature
        )
      )
    signature <- signature %>% 
      mutate_all(unlist)
    return(signature)
  }
}

get_concordant_signatures <- function(signature_file, signature_df) {
  if (missing(signature_file)) {
    if (missing(signature_df)) {
      stop("Either signature_file or a data-frame with the signature should be supplied")
    } else {
      signature_file <- tempfile()
      write_tsv(signature_df, signature_file)
    }
  } else {
    if (missing(signature_df) == FALSE) {
      stop("Only one of signature_file or signature_df should be supplied")
    }
  }
  
  url <- "http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze"
  query = list(lib = "LIB_5")
  body = list(file = upload_file(signature_file))
  
  request <- POST(url, query = query, body = body)
  
  if (status_code(request) == 200) {
    signature <- fromJSON(
      toJSON(
        content(request)$status$concordanceTable
      )
    )
    signature <- signature %>% 
      mutate_all(unlist)
    return(signature)
  }
}

generate_filtered_signature <- function(signature_df, signature_file, direction, threshold = 0.85) {
  
  if (missing(signature_df)) {
    if (missing(signature_file)) {
      stop("Either signature_file or a data-frame with the signature should be supplied")
    } else {
      if (str_detect(signature_file, ".csv")) {
        signature_df <- read_csv(signature_file)
      } else if (str_detect(signature_file, ".tsv")) {
        signature_df <- read_tsv(signature_file)
      } else {
        stop("Input file should either be tab or comma separated values")
      }
    }
  } else {
    if (missing(signature_file) == FALSE) {
      stop("Only one of signature_file or signature_df should be supplied")
    }
  }
  
  direction <- str_to_upper(direction)
  
  if (direction == "DOWN") {
    result <- signature_df %>% 
      filter(Value_LogDiffExp < -threshold)
  } else if (direction == "UP") {
    result <- signature_df %>% 
      filter(Value_LogDiffExp > threshold)
  } else {
    stop("Only valid directions are Up or Down")
  }
  
  result <- result %>% 
    select(ID_geneid, Name_GeneSymbol, Value_LogDiffExp, Significance_pvalue)
  
  return(result)
}

generate_consensus_signature <- function(up_result, down_result,
                                         cutoff = 0.321, cell_line = "all",
                                         discordant = FALSE) {
  
  if (missing(up_result)) {
    stop("Please specify the up-regulated genes connected signatures")
  } else if (missing(down_result)) {
    stop("Please specify the down-regulated genes connected signatures")
  }
  
  data <- bind_rows(up_result, down_result)
  
  if (cell_line != "all") {
    data <- data %>% 
      filter(cellline %in% cell_line)
  } else {
    
  }
  
  if (discordant) {
    result <- data %>% 
      filter(similarity < -cutoff) %>% 
      group_by(compound) %>%
      filter(similarity == max(similarity)) %>%
      ungroup() %>%
      select(signatureid, compound, similarity)
  } else {
    result <- data %>% 
      filter(similarity > cutoff) %>% 
      group_by(compound) %>%
      filter(similarity == max(similarity)) %>%
      ungroup() %>%
      select(signatureid, compound, similarity)
  }
  
  return(result)
  
}

process_groups <- function(...) {
  signatures <- unlist(list(...))
  dfs <- list()
  for (index in 1:length(signatures)) {
    dfs[[index]] <- get_l1000_signature(signatures[index])
  }
  reduced <- reduce(dfs, inner_join, by = "Name_GeneSymbol")
  symbols <- reduced$Name_GeneSymbol
  avg_logexp <- reduced %>% 
    select(starts_with("Value")) %>% 
    rowMeans()
  avged <- tibble(Name_GeneSymbol = symbols,
                  Value_LogDiffExp = avg_logexp)
  return(avged)
}


filter_grouped_signature <- function(df, direction,  threshold = 0.85) {
  
  direction <- str_to_upper(direction)
  
  if (direction == "DOWN") {
    result <- signature_df %>% 
      filter(Value_LogDiffExp < -threshold)
  } else if (direction == "UP") {
    result <- signature_df %>% 
      filter(Value_LogDiffExp > threshold)
  } else {
    stop("Only valid directions are Up or Down")
  }
  
  return(result)
}

combine_results <- function(x, y, prefix = NULL) {
  xd <- read_tsv(x)
  yd <- read_tsv(y)
  
  file_name <- str_extract(x, pattern = regex("(Group_.*[HA1E|MCF7])"))
  file_name <- paste(file_name, "tsv", sep = ".")
  path <- paste(prefix, file_name, sep = "/")
  
  d <- bind_rows(xd, yd)
  
  write_csv(d, path)
  invisible(d)
}
