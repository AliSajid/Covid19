# This file uses the feature counts file for GSE56192 available and creates the
# differential expression with MERS

library(edgeR)
library(org.Hs.eg.db)
library(tidyverse)

metadata <- read_tsv("raw/annotation/GSE56192-RunMapping.txt") %>%
  filter(str_detect(sample_alias, "\\_MERS.*High.*24hr.*") | str_detect(sample_alias, "\\_MOCK.*")) %>%
  extract(sample_alias, into = c("condition"), regex = "VMERS\\_(\\w+)\\-.*")

raw <- read_tsv("raw/annotation/GSE56192-featurecounts.tsv")
symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = raw$`Gene ID`, columns = c("SYMBOL", "ENTREZID"), keytype = "ENSEMBL")

counts <- raw %>%
  column_to_rownames("Gene ID") %>%
  rownames_to_column("geneID") %>%
  select(geneID, all_of(metadata$run_accession)) %>%
  inner_join(symbols, by = c("geneID" = "ENSEMBL")) %>%
  select(-geneID) %>%
  group_by(SYMBOL, ENTREZID) %>%
  summarise(across(where(is.numeric), sum)) %>%
  relocate(all_of(metadata$run_accession), .after = c("SYMBOL", "ENTREZID"))


dge <- DGEList(counts = counts[, 3:14], genes = counts[, 1:2], group = metadata$condition)

dge_all <- sumTechReps(dge, ID = metadata$sample_accession)

keep <- filterByExpr(dge_all)

dge_filtered <- dge_all[keep, , keep.lib.sizes = FALSE]

dge_filtered <- calcNormFactors(dge_filtered)

treatment <- as_factor(dge_filtered$samples$group)

design <- model.matrix(~ 0 + treatment)

dge_filtered <- estimateDisp(dge_filtered, design = design)

fit <- glmQLFit(dge_filtered, design)
qlf <- glmQLFTest(fit, contrast = c(1, -1))

l1000 <- read_tsv("raw/annotation/l1000-list.tsv") %>% pull(names)

top <- topTags(qlf, n = Inf)

table <- top$table %>%
  select(ENTREZID, SYMBOL, logFC, PValue) %>%
  rename(ID_geneid = ENTREZID,	Name_GeneSymbol = SYMBOL,
         Value_LogDiffExp = logFC,	Significance_pvalue = PValue) %>%
  filter(Name_GeneSymbol %in% l1000) %>%
  write_tsv("raw/mers-signature.tsv")
