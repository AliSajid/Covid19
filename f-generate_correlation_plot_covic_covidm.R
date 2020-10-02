# This script creates a scatterplot of concordance between covidc and covidm

library(tidyverse)

col_spec <- cols(
  Name_GeneSymbol = col_character(),
  Value_LogDiffExp = col_double(),
  Significance_pvalue = col_skip()
)

covidc <- read_tsv("data/signatures/disease/covidc-signature.tsv", col_types = col_spec) %>%
  rename(covidc = Value_LogDiffExp)
covidm <- read_tsv("data/signatures/disease/covidm-signature.tsv", col_types = col_spec) %>%
  rename(covidm = Value_LogDiffExp)

combined <- inner_join(covidc, covidm, by = "Name_GeneSymbol")

correlation <- cor.test(combined$covidc, combined$covidm)

g <- ggplot(combined, aes(x = covidc, y = covidm))

op1 <- g +
  geom_abline(slope = 1, intercept = 0, color = "grey80", lwd = 2) +
  geom_point() +
  xlab("LFC of DEG from China Human Dataset") +
  ylab("LFC of DEG from MtSinai Human Dataset") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(-15, 5, 2.5), limits = c(-15, 5)) +
  scale_y_continuous(breaks = seq(-15, 5, 2.5), limits = c(-15, 5)) +
  ggtitle("Correlation between DEG for China and Mt. Sinai Human Dataset",
          subtitle = str_glue("Pearson's r: {round(correlation$estimate, 3)}, p-value: {round(correlation$p.val, 3)}")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  coord_fixed()

ggsave("figures/covidc-covidm_correlation_plot.png", plot = op1)
ggsave("figures/covidc-covidm_correlation_plot.pdf", plot = op1)
