library(tidyverse)
library(RColorBrewer)
library(extrafont)

loadfonts()

col_spec <- cols(
  Candidates = col_character(),
  MCF7 = col_double(),
  HA1E = col_double()
)

specials <- c("Ruxolitinib", "Sirolimus")

top <- c("Ivermectin", "Genistein", "Alvocidib")

df <- read_csv("results/sars-ha1e-mcf7-combined-drugs.csv",
               col_names = c("Candidates", "MCF7", "HA1E"),
               col_types = col_spec,
               skip = 1) %>%
  mutate(Selected = if_else(HA1E >= 0.5 & MCF7 >= 0.5, "Yes", "No"),
         Selected = if_else(Candidates %in% specials, "Special", Selected),
         Selected = if_else(Candidates %in% top, "Top", Selected))

p <- ggplot(data=df,
            mapping = aes(x = MCF7,
                          y = HA1E,
                          size = 4,
                          shape = Selected,
                          color = Selected)
)


bp <- p +
  geom_hline(yintercept = 0.5, size = 1) +
  geom_hline(yintercept = 0.8, size = 1) +
  geom_vline(xintercept = 0.5, size = 1) +
  geom_vline(xintercept = 0.8, size = 1) +
  geom_point(stroke = 2) +
  xlab("\n\nAverage Reported Concordance Score for MCF7") +
  ylab("Average Reported Concordance Score for HA1E\n\n")

fp <- bp +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30, family = "Arial", hjust = 0.5),
    axis.title = element_text(size = 30, family = "Arial"),
    axis.text = element_text(size = 24, family = "Arial"),
    legend.text = element_text(size = 24, family = "Arial"),
    legend.title = element_text(size = 24, family = "Arial"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(1.5, "cm")
  ) +
  scale_shape_manual(name = element_blank(),
                     values = c(16, 1, 17, 15),
                     breaks = c("Yes", "No", "Top", "Special"),
                     labels = c("Passed Threshold", "Did Not Pass Threshold",
                                "Top Candidates",  substitute('Previously Identified '~italic(x), list(x="in silico")))) +
  scale_color_brewer(name = element_blank(),
                     palette = "Dark2",
                     labels = c("Passed Threshold", "Did Not Pass Threshold",
                                "Top Candidates", substitute('Previously Identified '~italic(x), list(x="in silico"))),
                     breaks = c("Yes", "No", "Top", "Special")) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0.3,1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    limits = c(0.3,1)
  ) +
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size=10))
  )

fp + ggtitle("Concordance Scatter plot for SARS")

ggsave("figures/SARS-Concordance-Scatterplot.png", device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")
