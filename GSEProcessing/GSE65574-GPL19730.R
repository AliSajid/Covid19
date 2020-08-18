# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Mar 24 11:59:03 EDT 2020

################################################################
#   Differential expression analysis with limma
# library(Biobase)
# library(GEOquery)
# library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE65574", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL19730", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

sml <- c(
  paste0(rep("mock", 5), c(rep("7", 3), rep("12", 2))),
  paste0(rep("ic", 5), c("12", rep("7", 2), rep("12", 2))),
  paste0(rep("fp", 6), c(rep("7", 3), rep("12", 3))),
  paste0(rep("dnsp16", 5), c(rep("7", 3), rep("12", 2))),
  paste0(rep("d4b", 6), c(rep("7", 3), rep("12", 3))),
  paste0(rep("d35", 4), c(rep("7", 2), rep("12", 2)))
)

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
cont.matrix <- makeContrasts(mock7_12=mock7-mock12, mock7_ic12=mock7-ic12, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250, coef="mock7_ic12")

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GENE_NAME"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
