library(tidyverse)
library(RColorBrewer)
library(extrafont)

loadfonts()

col_spec <- cols(
  .default = col_double(),
  Candidates = col_character()
)

specials <- c("Ruxolitinib", "Sirolimus")

top <- c("Ivermectin", "Genistein", "Alvocidib")

df <- read_csv("results/sars2-summarized-dataset.csv",
               col_names = c("Candidates", "Avg", "SDev", "SDLog10", "SDLog2"),
               col_types = col_spec,
               skip = 1) %>%
  mutate(Selected = if_else(Avg >= 0.5 & SDev <= 0.06, "Yes", "No"),
         Selected = if_else(Candidates %in% specials, "Special", Selected),
         Selected = if_else(Candidates %in% top, "Top", Selected))

p <- ggplot(data=df,
            mapping = aes(x = Avg,
                          y = SDev,
                          size = 4,
                          shape = Selected,
                          color = Selected)
)


bp <- p +
  geom_hline(yintercept = 0.3, size = 1) +
  geom_hline(yintercept = 0.5, size = 1) +
  geom_vline(xintercept = 0.05, size = 1) +
  geom_vline(xintercept = 0.10, size = 1) +
  geom_point(stroke = 2) +
  xlab("\n\nAverage Reported Concordance Scores across Cell Line Datasets") +
  ylab("Standard Deviation of Reportedd Concordance Scores across Cell Line Datasets\n\n")

fp <- bp +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30, family = "ArialMT", hjust = 0.5),
    axis.title = element_text(size = 30, family = "ArialMT"),
    axis.text = element_text(size = 24, family = "ArialMT"),
    legend.text = element_text(size = 24, family = "ArialMT"),
    legend.title = element_text(size = 24, family = "ArialMT"),
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
    limits = c(0.3,0.8)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.15, 0.01),
    limits = c(0, 0.15)
  ) +
  guides(size=FALSE,
         color = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size=10))
  )

full <- fp + ggtitle("Concordance Scatter plot for SARS-CoV-2")

ggsave("figures/SARS2-Concordance-Scatterplot.png", plot = full, device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")
