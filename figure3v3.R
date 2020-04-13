library(ggplot2)
library(tidyverse)
library(RColorBrewer)

col_spec <- cols(
  Candidate = col_character(),
  Concordance = col_double(),
  Discordance = col_double(),
  CellLine = col_factor(levels = NULL),
  Selected = col_character()
)

specials <- c("Ruxolitinib", "Sirolimus")

df <- read_csv("data/combined_drug_lists.csv",
               col_types = col_spec) %>% 
  mutate(Selected = if_else(Candidate %in% specials, "Special", Selected)) %>% 
  filter(Selected == "Yes")

p <- ggplot(data=df, 
            mapping = aes(x = Concordance,
                          y = Discordance,
                          size = 4,
                          shape = Selected,
                          color = Selected,
                          group = CellLine,
                          label = Candidate)
) + facet_grid(rows = vars(CellLine))


bp <- p + 
  geom_hline(yintercept = -0.3, size = 1) +
  #geom_hline(yintercept = -0.4, size = 1) +
  geom_hline(yintercept = -0.5, size = 1) +
  #geom_hline(yintercept = -0.6, size = 1) +
  geom_hline(yintercept = -0.7, size = 1) +
  #geom_hline(yintercept = -0.8, size = 1) +
  geom_vline(xintercept = 0.8, size = 1) +
  geom_vline(xintercept = 0.9, size = 1) +
  geom_vline(xintercept = 1, size = 1) +
  geom_point(stroke = 2) + 
  xlab("Average Reported Concordance Score for Each Cell Line") + 
  ylab("Average Reported Discordance Score for SARS")

fp <- bp + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 30, family = "Arial"),
    axis.text = element_text(size = 24, family = "Arial"),
    legend.text = element_text(size = 24, family = "Arial"),
    legend.title = element_text(size = 24, family = "Arial"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(1.5, "cm")
  ) +
  scale_shape_manual(name = element_blank(),
                     values = c(16, 1, 17),
                     breaks = c("Yes", "No", "Special"),
                     labels = c("Passed Threshold", "Did Not Pass Threshold", "Special Consideration")) +
  scale_color_brewer(name = element_blank(),
                     palette = "Dark2",
                     labels = c("Passed Threshold", "Did Not Pass Threshold", "Special Consideration"),
                     breaks = c("Yes", "No", "Special")) +
  scale_y_continuous(
    breaks = seq(-1, 0, 0.1),
    limits = c(-1, -0.3)
  ) + 
  scale_x_continuous(
    breaks = seq(0.70, 1, 0.05),
    limits = c(0.7,1)
  ) +
  guides(size=FALSE,
         color = FALSE,
         shape = FALSE
  )

fp + geom_text()

ggsave("figures/Concordance-Discordance-Scatterplot.png", device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")
