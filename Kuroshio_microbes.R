library(DESeq2)
library(tidyverse)
library(iNEXT)

Bac_un <- read.table(file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_unrareSeqXSt.csv", 
                     sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)
condition = factor(substr(colnames(Bac_un), nchar(colnames(Bac_un)), nchar(colnames(Bac_un))))

coldata = data.frame(
  row.names = colnames(Bac_un),
  condition = as.factor(substr(colnames(Bac_un), nchar(colnames(Bac_un)), nchar(colnames(Bac_un)))),
  type = as.factor(rep("paired-end", ncol(Bac_un)))
)
Bac_un_dds <- DESeqDataSetFromMatrix(as.matrix(Bac_un), 
                                     colData = coldata, 
                                     design = ~ condition)

ntd <- normTransform(Bac_un_dds)
library("vsn")
meanSdPlot(assay(ntd))

Bac_un_dds <- estimateSizeFactors(Bac_un_dds)
sizeFactors(Bac_un_dds)
Bac_un_dds <- estimateDispersions(Bac_un_dds)
Bac_vsd <- varianceStabilizingTransformation(Bac_un_dds, blind = FALSE)
write.table(assay(Bac_vsd), file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_unrareVSSeqXSt_2.csv",
            sep = ",", col.names = TRUE, row.names = TRUE)

# Bac_vsd <- vst(Bac_un_dds, blind = FALSE, nsub = 500)
# meanSdPlot(assay(Bac_vsd))
# write.table(assay(Bac_vsd), file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_unrareVSSeqXSt_2.csv",
#             sep = ",", col.names = TRUE, row.names = TRUE)

Bac_un_cds <- estimateSizeFactors(Bac_un_cds)
sizeFactors(Bac_un_cds)
Bac_Normcd <- counts(Bac_un_cds, normalized = TRUE)
write.table(Bac_Normcd, file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_unrareNormSeqXSt.csv",
            sep = ",", col.names = TRUE, row.names = TRUE)
# Bac_un_cds <- estimateDispersions(Bac_un_cds)
# plotDispEsts(Bac_un_cds)
# Bac_VScd <- ?getVarianceStabilizedData(Bac_un_cds)
# write.table(Bac_VScd, file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_unrareVSSeqXSt.csv",
#             sep = ",", col.names = TRUE, row.names = TRUE)
# 
# head(Bac_VScd)

library(DESeq2)
library(tidyverse)

Bac_un <- read.table(file = "https://raw.githubusercontent.com/OscarFHC/KuroshioBac/main/data/Bac_unrareSeqXSt.csv?token=AFL4RI5SZERQJ33SGIUDJ5S7UJTLO", 
                     sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)

Bac <- Bac_un
Mincomm_size <- min(colSums(Bac_un))
Meancomm_size <- mean(colSums(Bac_un))
Medcomm_size <- median(colSums(Bac_un))
Maxcomm_size <- max(colSums(Bac_un))

Bac_MinRare <- Bac_un
Bac_MeanRare <- Bac_un
Bac_MedRare <- Bac_un
Bac_MaxRare <- Bac_un
t <- Sys.time()
for (i in 1:ncol(Bac)){
  samp_temp <- matrix(0, nrow(Bac), 1000)
  for (k in 2:999){
    sp_list <- which(Bac[, i] != 0)
    fake_comm <- sample(sp_list, size = Medcomm_size, replace = TRUE, prob = Bac[Bac[, i] != 0, i]/sum(Bac[Bac[, i] != 0, i]))
    
    samp_temp[, 1] <- Bac[, i]
    for (j in 1:length(sp_list)){
      samp_temp[sp_list[j], k] <- length(which(fake_comm == sp_list[j]))
    }
    Bac_meanRare[, i] <- rowMeans(samp_temp)
  }
} 
Sys.time() - t

write.table(Bac_meanRare, file = "D:/Dropbox/Research/KuroshioBac_Assemb/data/Bac_MeanRareSeqXSt.csv",
            sep = ",", col.names = TRUE, row.names = TRUE)

  
  









