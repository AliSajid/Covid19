# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Mar 24 11:59:03 EDT 2020

################################################################
#   Differential expression analysis with limma
# library(Biobase)
# library(GEOquery)
# library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE56677", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17077", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

sml <- c("cov0", "mock7", "cov3", "mock3",
         "cov18", "cov7", "cov0", "cov7",
         "cov12", "mock3", "mock7", "mock12",
         "mock0", "cov18", "cov3", "cov18",
         "cov12", "cov24", "cov0", "mock18",
         "mock12", "mock7", "mock18", "cov24",
         "mock3", "cov3", "cov24", "mock12",
         "mock18", "cov12", "mock0", "cov7",
         "mock0")

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
#sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(g2g0=G2-G0, g1g0=G1-G0, g2g1=G2-G1, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GENE_SYMBOL","ENSEMBL_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
