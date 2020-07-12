library(tidyverse)

data <- read_csv("results/covid_cor_analysis_min.csv") %>% 
  select(perturbagen, threshold, min_n1, min_a2) %>% 
  mutate(threshold = as_factor(threshold))

p <- ggplot(data = data, mapping = aes(x= min_n1,y = min_a2, color = threshold))

p <- p + geom_point() +
  geom_vline(xintercept = -0.321, color = "blue", lwd = 0.5) +
  geom_hline(yintercept = -0.321, color = "blue", lwd = 0.5) +
  xlab("Concordant Drugs within the NHBE Cell Line") +
  ylab("Concordant Drugs within the A549 Cell Line") +
  theme_bw() +
  scale_x_continuous(breaks = seq(-1, 0.2, by = 0.1)) +
  scale_y_continuous(breaks = seq(-1, 0.2, by = 0.1)) +
  scale_color_brewer(name = "Threshold", palette = "Set1") +
  ggtitle("Scatter Plot of Concordance Scores") +
  theme(
     plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path("figures", "covid_drug_all_threshold_concordance.png"), plot = p)
ggsave(file.path("figures", "covid_drug_all_threshold_concordance.pdf"), plot = p)

p2 <- data %>% 
  filter(threshold == 0.85) %>% 
  ggplot(data = ., mapping = aes(x= min_n1, y = min_a2, color = threshold))

p2 <- p2 + geom_point() +
  geom_vline(xintercept = -0.321, color = "blue", lwd = 0.5) +
  geom_hline(yintercept = -0.321, color = "blue", lwd = 0.5) +
  xlab("Concordant Drugs within the NHBE Cell Line") +
  ylab("Concordant Drugs within the A549 Cell Line") +
  theme_bw() +
  scale_x_continuous(breaks = seq(-1, 0.2, by = 0.1)) +
  scale_y_continuous(breaks = seq(-1, 0.2, by = 0.1)) +
  scale_color_brewer(name = "Threshold", palette = "Set1") +
  ggtitle("Scatter Plot of Concordance Scores") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
  
ggsave(file.path("figures", "covid_drug_0.85_concordance.png"), plot = p2)
ggsave(file.path("figures", "covid_drug_0.85_concordance.pdf"), plot = p2)

p3 <- data %>% 
  filter(threshold == 0.85) %>% 
  select(perturbagen, min_a2) %>%
  mutate(passed = min_a2 < -0.321) %>% 
  ggplot(data = ., mapping = aes(x = perturbagen, y = min_a2, color = passed))

p3 <- p3 + geom_point() +
  coord_flip() +
  geom_hline(yintercept = -0.321, color = "blue", lwd = 0.5) +
  xlab("Drugs") +
  ylab("Concordance") +
  theme_bw() + 
  scale_y_continuous(breaks = seq(-1, 0, by = 0.05)) +
  scale_color_brewer(palette = "Set2", name = "Concordance < -0.321", labels = c("No", "Yes")) +
  ggtitle("Plot of Concordance Scores for A549") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path("figures", "covid_drug_A549_0.85_concordance.png"), plot = p3, height = 11.8, width = 11.8)
ggsave(file.path("figures", "covid_drug_A549_0.85_concordance.pdf"), plot = p3, height = 11.8, width = 11.8)

p4 <- data %>% 
  filter(threshold == 0.85) %>% 
  select(perturbagen, min_n1) %>%
  mutate(passed = min_n1 < -0.321) %>% 
  ggplot(data = ., mapping = aes(x = perturbagen, y = min_n1, color = passed))

p4 <- p4 + geom_point() +
  coord_flip() +
  geom_hline(yintercept = -0.321, color = "blue", lwd = 0.5) +
  xlab("Drugs") +
  ylab("Concordance") +
  theme_bw() + 
  scale_y_continuous(breaks = seq(-1, 0, by = 0.05)) +
  scale_color_brewer(palette = "Set2", name = "Concordance < -0.321", labels = c("No", "Yes")) +
  ggtitle("Plot of Concordance Scores for NHBE") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path("figures", "covid_drug_NHBE_0.85_concordance.png"), plot = p4, height = 11.8, width = 11.8)
ggsave(file.path("figures", "covid_drug_NHBE_0.85_concordance.pdf"), plot = p4, height = 11.8, width = 11.8)
