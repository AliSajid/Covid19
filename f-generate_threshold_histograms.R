library(tidyverse)
library(ggpubr)
library(e1071)
source("pipeline.R")

generate_graph <- function(filename, prefix) {
  parts <- str_match(filename, "(\\w+)\\-(\\w+)-signature.*")
  group <- parts[2]
  cell_line <- parts[3]
  data <- read_tsv(file.path(prefix, filename))

  normal <- shapiro.test(data$Value_LogDiffExp)
  skew <- skewness(data$Value_LogDiffExp)
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
    ggtitle(label = paste("Distribution of LogFC for", group),
            subtitle = paste(line, stat, pval, skewness, sep = "; ")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  d <- ggdensity(data$Value_LogDiffExp) +
    geom_vline(xintercept = c(-0.85, 0.85), color = "red", lwd = 2) +
    theme_bw() +
    xlab("Log Fold Change") +
    ylab("Density") +
    ggtitle(label = paste("Density plot of LogFC for", group),
            subtitle = paste(line, stat, pval, skewness, sep = "; ")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  q <- ggqqplot(data$Value_LogDiffExp) +
    theme_bw() +
    ggtitle(label = paste("Density plot of LogFC for", group),
            subtitle = paste(line, stat, pval, sep = "; ")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  out <- list(histogram = pl,
              density = d,
              qqplot = q,
              w = round(normal$statistic, 3),
              pval = round(normal$p.value, 3),
              line = cell_line,
              skewness = skew)
  return(out)
}

prefix <- file.path("data", "signatures", "group")
filenames <- list.files("data/signatures/group/")
out_prefix <- "figures"
create_dir(file.path(out_prefix, "histograms"))
create_dir(file.path(out_prefix, "densityplots"))
create_dir(file.path(out_prefix, "qqplots"))

data_skewness <- list()

for (file in filenames) {
  g <- generate_graph(file, prefix)
  data_skewness[[file]] <- g$skewness
  out_name_png <- gsub("tsv", "png", file)
  out_filenames <- str_replace(out_name_png, pattern = "signature", c("histogram", "densityplot", "qqplot"))

  ggsave(file.path(out_prefix, "histograms", out_filenames[1]), plot = g$histogram)
  ggsave(file.path(out_prefix, "densityplots", out_filenames[2]), plot = g$density)
  ggsave(file.path(out_prefix, "qqplots", out_filenames[3]), plot = g$qqplot)

}
