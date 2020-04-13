library(ggplot2)
library(tidyverse)
library(RColorBrewer)

df <- read_csv("data/combined_results.csv")

df <- df %>% 
  mutate(Presence = if_else(Common == T, "Common", Cell))

p <- ggplot(data=df, 
            mapping = aes( x = Concordance,
                           y = Tanimoto,
                           color = Drug,
                           size = 4,
                           shape = Presence)
            )


bp <- p + 
  geom_hline(yintercept = 0.5, size = 2) + 
  geom_hline(yintercept = 0.75, size = 2) + 
  geom_vline(xintercept = 0.321, size = 2) + 
  geom_point(stroke = 2) + 
  xlab("Maximum Reported Concordance Score") + 
  ylab("Tanimoto Score")

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
  guides(size=FALSE,
         colour = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size=10, color="black", fill = NULL))
         ) +
  scale_shape_manual(name = "Cell Lines",
                       values = c(19, 0, 2),
                       breaks = c("Common", "HA1E", "MCF7"),
                       labels = c("Both Cell Lines", "HA1E Only", "MCF7 Only")) +
  scale_color_brewer(name = "Reference Drugs", palette = "Dark2")

fp

ggsave("figures/Tanimoto-Concordance-Drugs.png", device = "png", width = 11.69 * 2, height = 8.27 * 2, units = "in")
