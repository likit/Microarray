# go to a directory containing CEL files
setwd('~/projects/mdv/affy_data/24hrs/')

# load required packages
library('affy')
library('limma')

# read all CEL files in the directory
data <- ReadAffy()

# obtain expression data
eset <- rma(data)

# read phenotype data from a file
targets <- read.table('phenoData.txt', header=T)

# select a particular set of data, according to time
target <- targets[targets$Hour=='24',]

# assign a phenotypic data to an expression data set
pData(eset) <- target

# create a design matrix
TS <- paste(target$Meq, target$Line, sep='.')
TS
TS <- factor(TS, levels=c("P.6", "N.6", "P.7", "N.7"))
design <- model.matrix(~0+TS)
design
colnames(design) <- levels(TS)

# fit the model
fit <- lmFit(eset, design)

# create a contrast matrix
cont.matrix <- makeContrasts(P.6vs7 = P.6 - P.7, N.6vs7 = N.6 - N.7, levels=design)

# fit a model to a contrast matrix
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#select genes with two-fold change and p.value >= 0.001
results <- decideTests(fit2, p.value=0.001, lfc=1)

# draw a venn diagram to a plot
pdf('6vs7_plot.pdf')
vennDiagram(results)
title(main='DGE at 24 hours, p.value=0.001, lfc=1')
dev.off()

# display top genes

P.table <- topTable(fit2, coef='P.6vs7', sort.by="logFC", number=nrow(eset))
N.table <- topTable(fit2, coef='N.6vs7', sort.by="logFC", number=nrow(eset))

write.table(P.table, file='Meq_6vs7.txt', sep='\t', row.names=F)
write.table(N.table, file='noMeq_6vs7.txt', sep='\t', row.names=F)

# create a contrast matrix
cont.matrix <- makeContrasts(PvsN.6 = P.6 - N.6, PvsN.7 = P.7 - N.7, levels=design)

# fit a model to a contrast matrix
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#select genes with two-fold change and p.value >= 0.001
results <- decideTests(fit2, p.value=0.001, lfc=1)

# draw a venn diagram to a plot
pdf('MeqvsNoMeq_plot.pdf')
vennDiagram(results)
title(main='DGE at 24 hours, p.value=0.001, lfc=1')
dev.off()

# display top genes

P.table <- topTable(fit2, coef='PvsN.6', sort.by="logFC", number=nrow(eset))
N.table <- topTable(fit2, coef='PvsN.7', sort.by="logFC", number=nrow(eset))

write.table(P.table, file='MeqvsNoMeq_6.txt', sep='\t', row.names=F)
write.table(P.table, file='MeqvsNoMeq_7.txt', sep='\t', row.names=F)
