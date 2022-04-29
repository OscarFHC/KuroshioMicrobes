library(vegan)
data(mite)
data(mite.env)
data(mite.pcnm)
mod <- varpart(mite, mite.env, mite.pcnm, data=mite.env, transfo="hel")

mite.all <- cbind(mite, mite.env, mite.pcnm)

data(dune)
data(dune.env)
mod <- ?adonis(mite ~ SubsDens + WatrCont, data = mite.env)
mod
summary(mod)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasilla")
#BiocManager::install("DESeq")
BiocManager::install("DESeq2")
BiocManager::install("vsn")

library(pasilla)
#library(DESeq)

datafile <- system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
datafile
pasillaCountTable <- read.table( datafile, header=TRUE, row.names=1 )
head( pasillaCountTable )

pasillaDesign = data.frame(
  row.names = colnames( pasillaCountTable ),
  condition = c("untreated", "untreated", "untreated",
                "untreated", "treated", "treated", "treated" ),
  libType = c("single-end", "single-end", "paired-end",
              "paired-end", "single-end", "paired-end", "paired-end" ) )
pasillaDesign

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]


cds <- newCountDataSet( countTable, condition )
cds <- estimateSizeFactors(cds) #%>% estimateDispersions
sizeFactors(cds)
cds_N <- counts( cds, normalized=TRUE )
cdsBlind <- estimateDispersions(cds, method = "blind")
varianceStabilizingTransformation(cds)
abund <- getVarianceStabilizedData(cds)
head(abund)

cdsFull = newCountDataSet( pasillaCountTable, pasillaDesign)
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )
fit1 = fitNbinomGLMs( cdsFull, count ~ libType + condition )
fit0 = fitNbinomGLMs( cdsFull, count ~ libType )
summary(fit1)
pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )

varianceStabilizingTransformation( cdsBlind )


Bac <- read.table(file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_rareSeqXSt.csv", 
                     sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)
Bac_Surf <- Bac[, substr(colnames(Bac), nchar(colnames(Bac)), nchar(colnames(Bac))) == "S"]
Bac_Chla <- Bac[, substr(colnames(Bac), nchar(colnames(Bac)), nchar(colnames(Bac))) == "C"]


Bac_cds <- newCountDataSet(Bac, BacDesign$condition)
Bac_cds <- estimateSizeFactors(Bac_cds)
sizeFactors(Bac_cds)
head(counts(Bac_cds, normalize = TRUE))

library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
str(coldata)

dds <- makeExampleDESeqDataSet(m=1000)
vsd <- varianceStabilizingTransformation(dds)

vsd2 <- vst(dds, blind = FALSE)

