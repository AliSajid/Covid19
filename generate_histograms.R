library(tidyverse)
library(ggpubr)
library(Cairo)

HA1E_prefix <- file.path("dataset", "groups", "HA1E", "raw")
MCF7_prefix <- file.path("dataset", "groups", "MCF7", "raw")
HA1E_files <- list.files(HA1E_prefix)
MCF7_files <- list.files(MCF7_prefix)

generate_graph <- function(filename, prefix) {
  name <- str_match(filename, "Group\\-(.*)\\-signature.*")[2]
  data <- read_tsv(file.path(prefix, filename))
  cell_line <- str_split(prefix, "/")[[1]][3]
  
  normal <- shapiro.test(data$Value_LogDiffExp)
  skew <- e1071::skewness(data$Value_LogDiffExp)
  stat <- paste("W =", round(normal$statistic, 3))
  pval <- paste("p =", round(normal$p.value, 3))
  line <- paste("Cell Line:", cell_line)
  skewness <- paste("Skewness:", skew)
  
  p <- ggplot(data, aes(x = Value_LogDiffExp))
  pl <- p + geom_histogram(binwidth = 0.5,
                     fill = "lightblue",
                     color = "black") + 
    scale_x_continuous(limits = c(-10, 10)) +
    xlab("Differential Expression") + 
    ylab("Frequency") + 
    geom_vline(xintercept = c(-0.85, 0.85), 
               color = "red", 
               lwd = 2) +
    ggtitle(label = paste("Distribution of LogFC for", name), 
            subtitle = paste(line, stat, pval, skewness, sep = "; "))
  
  d <- ggdensity(data$Value_LogDiffExp) + 
       geom_vline(xintercept = c(-0.85, 0.85), color = "red", lwd = 2) + 
       theme_bw() + 
       xlab("Log Fold Change") +
       ylab("Density") +
    ggtitle(label = paste("Density plot of LogFC for", name), 
            subtitle = paste(line, stat, pval, skewness, sep = "; "))
  
  q <- ggqqplot(data$Value_LogDiffExp) + 
    theme_bw() +
    ggtitle(label = paste("Density plot of LogFC for", name), 
            subtitle = paste(line, stat, pval, sep = "; "))
  
  out <- list(histogram = pl,
              density = d,
              qqplot = q,
              w = round(normal$statistic, 3),
              pval = round(normal$p.value, 3),
              line = cell_line,
              skewness = skew)
  return(out)
}

ha1e_res <- c()

for (h in HA1E_files) {
  g <- generate_graph(h, HA1E_prefix)
  ha1e_res <- c(ha1e_res, g$skewness)
  out_name_png <- gsub("tsv", "png", h)
  out_name_pdf <- gsub("tsv", "pdf", h)
  out_names_png <- str_replace(out_name_png, "Group\\-", 
                           paste("HA1E", "Group", c("histogram", "density", "qq"), "", sep = "-"))
  out_names_pdf <- str_replace(out_name_pdf, "Group\\-", 
                               paste("HA1E", "Group", c("histogram", "density", "qq"), "", sep = "-"))
  ggsave(file.path("figures", out_names_png[1]), plot = g$histogram)
  ggsave(file.path("figures", out_names_png[2]), plot = g$density)
  ggsave(file.path("figures", out_names_png[3]), plot = g$qqplot)
  ggsave(file.path("figures", out_names_pdf[1]), plot = g$histogram)
  ggsave(file.path("figures", out_names_pdf[2]), plot = g$density)
  ggsave(file.path("figures", out_names_pdf[3]), plot = g$qqplot)
}

mcf7_res <- c()

for (m in MCF7_files) {
  g <- generate_graph(m, MCF7_prefix)
  mcf7_res <- c(mcf7_res, g$skewness)
  out_name_png <- gsub("tsv", "png", m)
  out_name_pdf <- gsub("tsv", "pdf", m)
  out_names_png <- str_replace(out_name_png, "Group\\-", 
                               paste("MCF7", "Group", c("histogram", "density", "qq"), "", sep = "-"))
  out_names_pdf <- str_replace(out_name_pdf, "Group\\-", 
                               paste("MCF7", "Group", c("histogram", "density", "qq"), "", sep = "-"))
  ggsave(file.path("figures", out_names_png[1]), plot = g$histogram)
  ggsave(file.path("figures", out_names_png[2]), plot = g$density)
  ggsave(file.path("figures", out_names_png[3]), plot = g$qqplot)
  ggsave(file.path("figures", out_names_pdf[1]), plot = g$histogram)
  ggsave(file.path("figures", out_names_pdf[2]), plot = g$density)
  ggsave(file.path("figures", out_names_pdf[3]), plot = g$qqplot)
  }
