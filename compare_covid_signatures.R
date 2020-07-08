library(tidyverse)
source("pipeline.R")

files <- list.files("raw/sig_lists/")
in_prefix <- file.path("raw", "sig_lists")
out_prefix <- file.path("dataset", "drug_signatures")
dir.create(out_prefix, recursive = T)

for (drug in files) {
  infile <- file.path(in_prefix, drug)
  name <- gsub(".tsv", "", drug)
  print(paste("Processing drug:", drug))
  if (!dir.exists(file.path(out_prefix, name))) {
    dir.create(file.path(out_prefix, name))
  }
  sig_ids <- read_tsv(infile) %>% filter(CellLine %in% c("HA1E", "MCF7")) %>% pull(SignatureId)
  for (id in sig_ids) {
    if (!file.exists(file.path(out_prefix, name, paste(id, "tsv", sep=".")))) {
      print(paste("Downloading signature for", id))
      sig <- get_l1000_signature(id)
      sig %>% write_tsv(file.path(out_prefix, name, paste(id, "tsv", sep=".")))
    } else {
      print(paste(id, "signature already exists"))
    }
  }
}

covidc <- read_tsv("dataset/disease/covidc/covidc-l1000-signature.tsv") %>% arrange(Name_GeneSymbol)
covidm <- read_tsv("dataset/disease/covidm/covidm-l1000-signature.tsv") %>% arrange(Name_GeneSymbol)
dNHBE_1 <- read_tsv("dataset/disease/dNHBE_1/dNHBE_1_l1000_de.tsv") %>% arrange(Name_GeneSymbol)
dA549_2 <- read_tsv("dataset/disease/dA549_2/dA549_2_l1000_de.tsv") %>% arrange(Name_GeneSymbol)
dA549_3 <- read_tsv("dataset/disease/dA549_3/dA549_3_l1000_de.tsv") %>% arrange(Name_GeneSymbol)
dA549_p <- read_tsv("dataset/disease/dA549_p/dA549_p_l1000_de.tsv") %>% arrange(Name_GeneSymbol)
dACE2_4 <- read_tsv("dataset/disease/dACE2_4/dACE2_4_l1000_de.tsv") %>% arrange(Name_GeneSymbol)
dCalu3_5 <- read_tsv("dataset/disease/dCalu3_5/dCalu3_5_l1000_de.tsv") %>% arrange(Name_GeneSymbol)

get_correlation <- function(covid, drug, threshold = 0) {
  c <- covid %>% 
       group_by(Name_GeneSymbol) %>% 
       summarise(Value_LogDiffExp = max(Value_LogDiffExp), .groups = "drop_last") %>% 
       filter(abs(Value_LogDiffExp) > threshold) %>% 
       select(Name_GeneSymbol, Value_LogDiffExp) %>% 
       arrange(Name_GeneSymbol)
  
  common <- intersect(c$Name_GeneSymbol, drug$Name_GeneSymbol)
  
  d <- drug %>% 
       filter(Name_GeneSymbol %in% common) %>% 
       group_by(Name_GeneSymbol) %>% 
       summarise(Value_LogDiffExp = max(Value_LogDiffExp), .groups = "drop_last") %>% 
       select(Name_GeneSymbol, Value_LogDiffExp) %>% 
       arrange(Name_GeneSymbol)
  
  c <- c %>% 
    filter(Name_GeneSymbol %in% common) %>% arrange(Name_GeneSymbol)
  
  res <- cor(c$Value_LogDiffExp, d$Value_LogDiffExp, method="pearson")
  return(res)
}

thresholds <- c(0, 0.26, 0.5, 0.85, 1)

sig <- c()
pert <- c()
thresh <- c()
covidc_cor <- c()
covidm_cor <- c()
dNHBE_1_cor <- c()
dA549_2_cor <- c()
dA549_3_cor <- c()
dA549_p_cor <- c()
dACE2_4_cor <- c()
dCalu3_5_cor <- c()

for (th in thresholds) {
  print(paste("Processing threshold:", th))
  for (drug in files) {
    name <- gsub(".tsv", "", drug)
    print(paste("Processing drug:", drug))
    signatures <- list.files(file.path(out_prefix, name))
    for (sid in signatures) {
      spec <- cols(
        signatureID = col_skip(),
        PROBE = col_skip(),
        ID_geneid = col_skip(),
        Name_GeneSymbol = col_character(),
        Value_LogDiffExp = col_double(),
        Significance_pvalue = col_skip()
      )
      print(paste("Processing signature:", sig))
      d <- read_tsv(file.path(out_prefix, name, sid), col_types = spec)
      cc <- get_correlation(covidc, d, th)
      cm <- get_correlation(covidm, d, th)
      cn1 <- get_correlation(dNHBE_1, d, th)
      ca2 <- get_correlation(dA549_2, d, th)
      ca3 <- get_correlation(dA549_3, d, th)
      cap <- get_correlation(dA549_p, d, th)
      ca4 <- get_correlation(dACE2_4, d, th)
      cc5 <- get_correlation(dCalu3_5, d, th)
      pert <- c(pert, name)
      sig <- c(sig, sid)
      thresh <- c(thresh, th)
      covidc_cor <- c(covidc_cor, cc)
      covidm_cor <- c(covidm_cor, cm)
      dNHBE_1_cor <- c(dNHBE_1_cor, cn1)
      dA549_2_cor <- c(dA549_2_cor, ca2)
      dA549_3_cor <- c(dA549_3_cor, ca3)
      dA549_p_cor <- c(dA549_p_cor, cap)
      dACE2_4_cor <- c(dACE2_4_cor, ca4)
      dCalu3_5_cor <- c(dCalu3_5_cor, cc5)
    }
  }
}

cor_data <- tibble(signature = sig,
                   perturbagen = pert,
                   threshold = thresh,
                   mtsinai_cor = covidm_cor,
                   china_cor = covidc_cor,
                   NHBE_1_cor = dNHBE_1_cor,
                   A549_2_cor = dA549_2_cor,
                   A549_3_cor = dA549_3_cor,
                   A549_p_cor = dA549_p_cor,
                   ACE2_4_cor = dACE2_4_cor,
                   Calu3_5_cor = dCalu3_5_cor)

cor_analysis <- cor_data %>% 
  group_by(perturbagen, threshold) %>% 
  summarize(pos_m = sum(mtsinai_cor > 0), 
            neg_m = sum(mtsinai_cor <= 0),
            min_m = min(mtsinai_cor),
            max_m = max(mtsinai_cor),
            mean_m = mean(mtsinai_cor),
            sd_m = sd(mtsinai_cor),
            pos_c = sum(china_cor > 0), 
            neg_c = sum(china_cor <= 0),
            min_c = min(china_cor),
            max_c = max(china_cor),
            mean_c = mean(china_cor),
            sd_c = sd(china_cor),
            pos_n1 = sum(NHBE_1_cor > 0), 
            neg_n1 = sum(NHBE_1_cor <= 0),
            min_n1 = min(NHBE_1_cor),
            max_n1 = max(NHBE_1_cor),
            mean_n1 = mean(NHBE_1_cor),
            sd_n1 = sd(NHBE_1_cor),
            pos_a2 = sum(A549_2_cor > 0), 
            neg_a2 = sum(A549_2_cor <= 0),
            min_a2 = min(A549_2_cor),
            max_a2 = max(A549_2_cor),
            mean_a2 = mean(A549_2_cor),
            sd_a2 = sd(A549_2_cor),
            pos_a3 = sum(A549_3_cor > 0), 
            neg_a3 = sum(A549_3_cor <= 0),
            min_a3 = min(A549_3_cor),
            max_a3 = max(A549_3_cor),
            mean_a3 = mean(A549_3_cor),
            sd_a3 = sd(A549_3_cor),
            pos_ap = sum(A549_p_cor > 0), 
            neg_ap = sum(A549_p_cor <= 0),
            min_ap = min(A549_p_cor),
            max_ap = max(A549_p_cor),
            mean_ap = mean(A549_p_cor),
            sd_ap = sd(A549_p_cor),
            pos_a4 = sum(ACE2_4_cor > 0), 
            neg_a4 = sum(ACE2_4_cor <= 0),
            min_a4 = min(ACE2_4_cor),
            max_a4 = max(ACE2_4_cor),
            mean_a4 = mean(ACE2_4_cor),
            sd_a4 = sd(ACE2_4_cor),
            pos_c5 = sum(Calu3_5_cor > 0), 
            neg_c5 = sum(Calu3_5_cor <= 0),
            min_c5 = min(Calu3_5_cor),
            max_c5 = max(Calu3_5_cor),
            mean_c5 = mean(Calu3_5_cor),
            sd_c5 = sd(Calu3_5_cor)
            ) %>% 
  mutate(m_ratio = pos_m/neg_m, c_ratio = pos_c/neg_c,
         n1_ratio = pos_n1/neg_n1, a2_ratio = pos_a2/neg_a2,
         a3_ratio = pos_a3/neg_a3, ap_ratio = pos_ap/neg_ap,
         a4_ratio = pos_a4/neg_a4, c5_ratio = pos_c5/neg_c5
         )

cor_analysis %>% 
  write_csv("results/covid_cor_analysis.csv")

cor_min <- cor_analysis %>% 
  select(perturbagen, threshold, starts_with("min"))

cor_min %>% write_csv("results/covid_cor_min_analysis.csv")

cor_max <- cor_analysis %>% 
  select(perturbagen, threshold,  starts_with("max"))

cor_max %>% write_csv("results/covid_cor_max_analysis.csv")

cor_mean <- cor_analysis %>% 
  select(perturbagen, threshold,  starts_with("mean"))

cor_mean %>% write_csv("results/covid_cor_mean_analysis.csv")

cor_sd <- cor_analysis %>% 
  select(perturbagen, threshold,  starts_with("sd"))

cor_sd %>% write_csv("results/covid_cor_sd_analysis.csv")

cor_ratio <- cor_analysis %>% 
  select(perturbagen, threshold,  ends_with("ratio"))

cor_ratio %>% write_csv("results/covid_cor_ratio_analysis.csv")

cor_min %>% 
  pivot_longer(starts_with("min"), names_to = "sample", values_to = "concordance") %>% 
  mutate(threshold = as_factor(threshold),
         pass = as_factor(if_else(concordance < -0.321, "Passed", "Failed"))) %>% 
  ggplot(data = ., aes(x = sample, y = concordance, color = pass)) -> g

gp <- g +  geom_point() + 
  facet_grid(~ threshold) + 
  coord_flip() + 
  geom_hline(yintercept = -0.321, color = "darkgreen", linetype = 2) +
  xlab("Concordance Scores") +
  ylab("Samples") +
  scale_y_continuous(breaks = seq(-0.8, 0.2, 0.2)) + 
  ggtitle("Distribution of Min Concordance scores") +
  theme_bw() +
  scale_color_brewer(name = "Threshold (-0.321)", palette = "Set1") +
  scale_x_discrete(breaks = c("min_m", "min_c", "min_n1", "min_a2", "min_a3", "min_a4", "min_c5", "min_ap"),
                   limits = c("min_m", "min_c", "min_n1", "min_a2", "min_a3", "min_a4", "min_c5", "min_ap"),
                   labels = c("MtSinai", "China", "NHBE", "A549-1", "A549-2", "ACE2-ko", "Calu3", "A549-pooled"))

ggsave("figures/min-concordance-plot-covid.png", plot = gp, units = "in", width = 11, height = 7)
