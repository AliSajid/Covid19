library(tidyverse)
library(RColorBrewer)
library(extrafont)

loadfonts()

col_spec <- cols(
  .default = col_double(),
  compound = col_character()
)

antiviral <- c("Gemcitabine", "Trametinib", "Withaferin A", "Saracatinib",
               "Erlotinib", "Alvocidib", "Itrazole", "Elesclomol",
               "Dasatinib", "Panobinostat")

trial <- c("Gallocatechin Gallate", "Genistein", "Imatinib",
           "Dexamethasone Acetate", "Simvastatin", "Sirolimus",
           "Tamoxifen")

data <-
  read_csv("results/sars2-summarized-dataset.csv", col_types = col_spec) %>%
  mutate(
    class = if_else((avg >= 0.5 & sdev <= 0.06), "novel", "filtered"),
    class = if_else(compound %in% antiviral, "antiviral", class),
    class = if_else(compound %in% trial, "trial", class),
  ) %>%
  select(compound, avg, sdev, class) %>%
  filter(
    str_detect(compound, "CHEMBL", negate = T),
    str_detect(compound, "SCHEMBL", negate = T),
    str_detect(compound, "^\\d+", negate = T),
    str_detect(compound, "^[A-Z]\\d*\\w*\\-?\\s?\\d+", negate = T),
    str_detect(compound, "[Ii]nhibitor", negate = T),
    str_detect(compound, "^Broad", negate = T),
    str_detect(compound, "^BRD*", negate = T),
    str_detect(compound, "^UNII", negate = T),
    str_detect(compound, "omer", negate = T),
    str_detect(compound, "^Tyrphostin", negate = T)
  )

g <- ggplot(data, aes(x = avg, y = sdev, group = class, color = class, shape = class))

p <- g + geom_point(size = 5) +
  scale_x_continuous(breaks = seq(0.32, 0.7, 0.02), limits = c(0.4, 0.7)) +
  scale_y_reverse(breaks = seq(0, 0.15, 0.01), limits = c(0.10, 0)) +
  geom_vline(xintercept = 0.5, color = "darkred", lwd = 1.5) +
  geom_hline(yintercept = 0.06, color = "darkred", lwd = 1.5) +
  scale_color_manual(breaks = c("antiviral", "trial", "novel", "filtered"),
                     labels = c("Known Anti-viral Properties",
                                "Currently in trial",
                                "Novel identified",
                                "Did not pass Threshold"),
                     name = "Drug Status",
                     values = c("darkgreen", "maroon",
                                "darkblue", "lightgrey")) +
  scale_shape_manual(breaks = c("antiviral", "trial", "novel", "filtered"),
                     labels = c("Known Anti-viral Properties",
                                "Currently in trial",
                                "Novel identified",
                                "Did not pass Threshold"),
                     name = "Drug Status",
                     values = c(17, 18, 15, 16)) +
  xlab("Mean Concordance") +
  ylab("Standard Deviation of Concordance") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30, family = "ArialMT", hjust = 0.5),
    axis.title = element_text(size = 30, family = "ArialMT"),
    axis.text = element_text(size = 24, family = "ArialMT"),
    legend.text = element_text(size = 24, family = "ArialMT"),
    legend.title = element_text(size = 24, family = "ArialMT"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(1.5, "cm")
  )

p

full <- p + ggtitle("Concordance Plot of Identified Drugs")

ggsave("figures/SARS2-Mean-SD-Plot.png", plot = full, device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")

ggsave("figures/SARS2-Mean-SD-Plot.jpg", plot = full, device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")
