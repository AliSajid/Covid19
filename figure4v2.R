library(tidyverse)
library(RColorBrewer)

special <-
  c("Niclosamide",
    "Ruxolitinib",
    "Sirolimus",
    "Valsartan",
    "Everolimus",
    "Ivermectin",
    "Sitagliptin")

col_spec <- cols(
  .default = col_double(),
  perturbagen = col_factor(),
  threshold = col_factor()
)

data <- read_csv("results/covid_cor_analysis_min.csv", col_types = col_spec) %>% 
  select(perturbagen, threshold, min_n1, min_a2) %>% 
  mutate(threshold = as_factor(threshold),
         arank = as.integer(rank(min_a2)),
         nrank = as.integer(rank(min_n1))) %>% 
  filter(threshold == 0.85) %>% 
  select(perturbagen, min_a2, arank) %>%
  mutate(passed = if_else(min_a2 < -0.321, "yes", "no"),
         passed = if_else(perturbagen %in% special, "special", passed),
         arank = as.integer(rank(arank)),
         perturbagen = reorder(perturbagen, arank)) %>% 
  arrange(min_a2)


p <- ggplot(data = data, mapping = aes(x = perturbagen, y = min_a2, color = passed, shape = passed))

mp <- p + geom_point() +
  coord_flip() +
  geom_hline(yintercept = -0.321, color = "blue", lwd = 0.5) +
  xlab("Drugs") +
  ylab("Concordance") +
  theme_bw() + 
#  scale_x_discrete(breaks = data$perturbagen[scale]) +
  scale_y_continuous(breaks = seq(-1, 0, by = 0.05)) +
  scale_shape_manual(name = element_blank(),
                     values = c(16, 1, 17),
                     breaks = c("yes", "no", "special"),
                     labels = c("Passed Threshold", "Did Not Pass Threshold", "Special Consideration")) +
  scale_color_brewer(name = element_blank(),
                     palette = "Dark2",
                     labels = c("Passed Threshold", "Did Not Pass Threshold", "Special Consideration"),
                     breaks = c("yes", "no", "special")) +
  ggtitle("Plot of Concordance Scores for A549") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

mp

ggsave(file.path("figures", "covid_drug_A549_0.85_concordance.png"), plot = mp, height = 11.8, width = 11.8)
ggsave(file.path("figures", "covid_drug_A549_0.85_concordance.pdf"), plot = mp, height = 11.8, width = 11.8)
