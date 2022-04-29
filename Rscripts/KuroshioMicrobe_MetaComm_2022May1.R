###############################################################################################
##### Loading packages ########################################################################
###############################################################################################
if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(vegan)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

if (!require(mgcv)) {
  install.packages("mgcv", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(mgcv)
}else{library(mgcv)}

if (!require(ape)) {
  install.packages("ape", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ape)
}else{library(ape)}

if (!require(parallel)) {
  install.packages("vegan", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(parallel)
}else{library(parallel)}

if (!require(SpadeR)) {
  install.packages("SpadeR", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(SpadeR)
}else{library(SpadeR)}

if (!require(iNEXT)) {
  install.packages("iNEXT", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(iNEXT)
}else{library(iNEXT)}

if (!require(picante)) {
  install.packages("picante", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(picante)
}else{library(picante)}

if (!require(geiger)) {
  install.packages("geiger", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(geiger)
}else{library(geiger)}

if (!require(GUniFrac)) {
  install.packages("GUniFrac", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GUniFrac)
}else{library(GUniFrac)}

if (!require(phytools)) {
  install.packages("phytools", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(phytools)
}else{library(phytools)}

if (!require(ecodist)) {
  install.packages("ecodist", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ecodist)
}else{library(ecodist)}

if (!require(ade4)) {
  install.packages("ade4", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ade4)
}else{library(ade4)}

if (!require(nlme)) {
  install.packages("nlme", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(nlme)
}else{library(nlme)}

if (!require(geosphere)) {
  install.packages("geosphere", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(geosphere)
}else{library(geosphere)}

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(viridis)) {
  install.packages("viridis", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(viridis)
}else{library(viridis)}

if (!require(GGally)) {
  install.packages("GGally", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(GGally)
}else{library(GGally)}

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}
###############################################################################################
##### Loading packages ########################################################################
###############################################################################################

###############################################################################################
##### Loading functions #######################################################################
###############################################################################################
cor_fun <- function(data, mapping, method="pearson", ndp=2, sz=5, stars=TRUE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  x <- Kuro.ADiv$Temperature
  y <- Kuro.ADiv$Salinity
  corr <- cor.test(x, y, method = "pearson")
  est <- corr$estimate

  if(stars){
    stars <- c("*", "")[findInterval(corr$p.value, c(0, 0.1, 1))]
    lbl <- paste0(round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data = data, mapping = mapping) + 
    annotate("text", x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE), label = lbl, size = lb.size, ...)+
    theme(panel.grid = element_blank())
}

fit_fun <- function(data, mapping, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point()
  
  if (cor.test(x, y, method = "pearson")$p.value < 0.1){
    p <- p + 
      geom_smooth(method = loess, fill = "red", color = "red", ...) +
      geom_smooth(method = lm, fill = "blue", color = "blue", ...)
    p
  }
  
  p
  
}
###############################################################################################
##### Loading functions #######################################################################
###############################################################################################

###############################################################################################
##### Loading and preping data ################################################################
###############################################################################################
HNF.raw <- as.data.frame(t(read.table(
  file = "D:/Dropbox/Research/KuroshioMicrobes/data/KuroHNF_SeqXSt.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
HNF.ra <- HNF.raw / rowSums(HNF.raw)

Bac.raw <- as.data.frame(t(read.table(
  file = "D:/Dropbox/Research/KuroshioMicrobes/data/KuroBac_rareSeqXSt.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac.ra <- Bac.raw / rowSums(Bac.raw)

MetaData <- read.table(file = "D:/Dropbox/Research/KuroshioMicrobes/data/Kuroshio_metadata.csv", sep = ",", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
###############################################################################################
##### Loading and preping data ################################################################
###############################################################################################


#################################
##### Surface layer #############
#################################
##### prepare distance for calculation #####
which(substr(rownames(Bac.raw), nchar(rownames(Bac.raw)), nchar(rownames(Bac.raw))) == "S")

BacSurf <- Bac.raw[which(substr(rownames(Bac.raw), nchar(rownames(Bac.raw)), nchar(rownames(Bac.raw))) == "S"), ]
HNFSurf <- HNF.raw[which(substr(rownames(HNF.raw), nchar(rownames(HNF.raw)), nchar(rownames(HNF.raw))) == "S"), ]
BacSurf.dist <- matrix(as.matrix(vegdist(BacSurf, method = "bray")), ncol = 1)
HNFSurf.dist <- matrix(as.matrix(vegdist(HNFSurf, method = "bray")), ncol = 1)

EnviS <- MetaData %>%
  filter(Layer == "Surface") %>%
  select(Cruise, Station, Line, Longitude, Latitude, Temperature, Salinity, DIN, P, Chla) %>%
  mutate(Temp.s = as.numeric(scale(Temperature)),
         Sal.s = as.numeric(scale(Salinity)),
         TN.s = as.numeric(scale(DIN)),
         TP.s = as.numeric(scale(P)),
         Chla.s = as.numeric(scale(Chla))
         )
EnviS.dist <- matrix(as.matrix(vegdist(EnviS[, c("Temp.s", "Sal.s", "TN.s", "TP.s", "Chla.s")], method = "euclidean")), ncol = 1)
GDist <-     distm(EnviS[, c("Longitude", "Latitude")], fun = distGeo)
GeoS.dist <- matrix(distm(EnviS[, c("Longitude", "Latitude")], fun = distGeo), ncol = 1)

Surf.dist <- data.frame(BacDist = BacSurf.dist, 
                        HNFDist = HNFSurf.dist, 
                        EnviDist = EnviS.dist, 
                        GeoDist = GeoS.dist) %>%
  mutate(focal = rep(rownames(BacSurf), each = length(rownames(BacSurf))),
         to = rep(rownames(BacSurf), length(rownames(BacSurf))),
         transect = rep(EnviS$Line, each = length(EnviS$Line)),
         Cr = rep(EnviS$Line, each = length(EnviS$Line))) %>%
  filter(BacDist > 0) %>%
  mutate(BacDist.s = as.numeric(scale(BacDist)),
         HNFDist.s = as.numeric(scale(HNFDist)),
         EnviDist.s = as.numeric(scale(EnviDist)),
         GeoDist.s = as.numeric(scale(GeoDist))) %>%
  mutate(BacDist.Res = residuals(lm(BacDist.s ~ EnviDist.s + GeoDist.s)),
         HNFDist.Res = residuals(lm(HNFDist.s ~ EnviDist.s + GeoDist.s)),
         EnviDist.BacRes = residuals(lm(EnviDist.s ~ BacDist.s + GeoDist.s)),
         EnviDist.HNFRes = residuals(lm(EnviDist.s ~ HNFDist.s + GeoDist.s)),
         GeoDist.BacRes = residuals(lm(GeoDist.s ~ BacDist.s + EnviDist.s)),
         GeoDist.HNFRes = residuals(lm(GeoDist.s ~ HNFDist.s + EnviDist.s))
  )
##### prepare distance for calculation #####

########## for bacteria ##########
##### linear model #####
# BacSurf.lme <- lme(BacDist ~ HNFDist + EnviDist + GeoDist, random = ~ 1 | transect, method = "ML", 
#                    data = Surf.dist, na.action = na.omit)
# summary(BacSurf.lme)
# performance::r2(BacSurf.lme) 
BacSurf.lmeAll <- lm(BacDist.s ~ HNFDist.s + EnviDist.s + GeoDist.s,# random = ~ 1 | transect, method = "ML", 
                     data = Surf.dist, na.action = na.omit)
summary(BacSurf.lmeAll) 
performance::r2(BacSurf.lmeAll)

BacSurf.lmeHNF <- lm(BacDist.s ~ HNFDist.Res, #random = ~ 1 | transect, method = "ML", 
                    data = Surf.dist, na.action = na.omit)
summary(BacSurf.lmeHNF) 
performance::r2(BacSurf.lmeHNF)

BacSurf.lmeEnvi <- lm(BacDist.s ~ EnviDist.HNFRes, #random = ~ 1 | transect, method = "ML", 
                      data = Surf.dist, na.action = na.omit)
summary(BacSurf.lmeEnvi) 
performance::r2(BacSurf.lmeEnvi)

BacSurf.lmeGeo <- lm(BacDist.s ~ GeoDist.HNFRes, #random = ~ 1 | transect, method = "ML", 
                       data = Surf.dist, na.action = na.omit)
summary(BacSurf.lmeGeo) 
performance::r2(BacSurf.lmeGeo)
##### linear model #####

##### PERMANOVA #####
BacSurf.perMANOVA <- adonis(BacDist.s ~ HNFDist.s + EnviDist.s + GeoDist.s, data = Surf.dist, method = "euclidean")
#
BacSurf.perMANOVA
summary(BacSurf.perMANOVA)
##### PERMANOVA #####

##### Varpart #####
BacSurf.varpart <- varpart(Surf.dist$BacDist.s, Surf.dist$HNFDist.s, Surf.dist$EnviDist.s, Surf.dist$GeoDist.s)
BacSurf.varpart
rda(Surf.dist$BacDist.s ~ Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)
rda(Surf.dist$BacDist.s ~ Surf.dist$HNFDist.s + Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)
anova.cca(rda(Surf.dist$BacDist.s ~ Surf.dist$HNFDist.s + Condition(Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)))
anova.cca(rda(Surf.dist$BacDist.s ~ Surf.dist$EnviDist.s + Condition(Surf.dist$HNFDist.s + Surf.dist$GeoDist.s)))
anova.cca(rda(Surf.dist$BacDist.s ~ Surf.dist$GeoDist.s + Condition(Surf.dist$HNFDist.s + Surf.dist$EnviDist.s)))
##### Varpart #####

##### Plotting #####
BacSurfDist.p <- Surf.dist %>%
  select(BacDist, HNFDist, EnviDist, GeoDist) %>%
  mutate(GeoDist = GeoDist/1000) %>%
  gather(key = "ExpVar", value = "Distance", -BacDist) %>%
  mutate(ExpVar = factor(ExpVar, 
                         levels = c("HNFDist", "EnviDist", "GeoDist"),
                         labels = c(expression(paste("HNF community dissimilarity")),
                                    expression(paste("Environmental dissimilarity")),
                                    expression(paste("Geographical distance (km)")) )
  )
  ) %>%
  ggplot(aes(x = Distance, y = BacDist)) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  facet_grid(cols = vars(ExpVar), scales = "free", labeller = "label_parsed") +
  labs(x = expression(""),
       y = expression("Bacterial community dissimilarity")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 24)
  )
BacSurfDist.p
ggsave(BacSurfDist.p, file = "D:/Dropbox/Manuscripts/KuroshioMicrobes/MetaComm/Figs/BacSurfDist.png",
       dpi = 600, width = 60, height = 32, units = "cm")
##### Plotting #####
########## for bacteria ##########

########## for HNF ##########
##### linear model #####
# BacSurf.lme <- lme(BacDist ~ HNFDist + EnviDist + GeoDist, random = ~ 1 | transect, method = "ML", 
#                    data = Surf.dist, na.action = na.omit)
# summary(BacSurf.lme)
# performance::r2(BacSurf.lme) 
HNFSurf.lmeAll <- lm(HNFDist.s ~ BacDist.s + EnviDist.s + GeoDist.s,# random = ~ 1 | transect, method = "ML", 
                     data = Surf.dist, na.action = na.omit)
summary(HNFSurf.lmeAll) 
performance::r2(HNFSurf.lmeAll)

HNFSurf.lmeBac <- lm(HNFDist.s ~ BacDist.Res, #random = ~ 1 | transect, method = "ML", 
                     data = Surf.dist, na.action = na.omit)
summary(HNFSurf.lmeBac) 
performance::r2(HNFSurf.lmeBac)

HNFSurf.lmeEnvi <- lm(HNFDist.s ~ EnviDist.BacRes, #random = ~ 1 | transect, method = "ML", 
                      data = Surf.dist, na.action = na.omit)
summary(HNFSurf.lmeEnvi) 
performance::r2(HNFSurf.lmeEnvi)

HNFSurf.lmeGeo <- lm(HNFDist.s ~ GeoDist.BacRes, #random = ~ 1 | transect, method = "ML", 
                     data = Surf.dist, na.action = na.omit)
summary(HNFSurf.lmeGeo) 
performance::r2(HNFSurf.lmeGeo)
##### linear model #####

##### PERMANOVA #####
HNFSurf.perMANOVA <- adonis(HNFDist.s ~ BacDist.s + EnviDist.s + GeoDist.s, data = Surf.dist, method = "euclidean")
#
HNFSurf.perMANOVA
summary(HNFSurf.perMANOVA)
##### PERMANOVA #####

##### Varpart #####
HNFSurf.varpart <- varpart(Surf.dist$HNFDist.s, Surf.dist$BacDist.s, Surf.dist$EnviDist.s, Surf.dist$GeoDist.s)
HNFSurf.varpart
rda(Surf.dist$HNFDist.s ~ Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)
rda(Surf.dist$HNFDist.s ~ Surf.dist$BacDist.s + Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)
anova.cca(rda(Surf.dist$HNFDist.s ~ Surf.dist$BacDist.s + Condition(Surf.dist$EnviDist.s + Surf.dist$GeoDist.s)))
anova.cca(rda(Surf.dist$HNFDist.s ~ Surf.dist$EnviDist.s + Condition(Surf.dist$BacDist.s + Surf.dist$GeoDist.s)))
anova.cca(rda(Surf.dist$HNFDist.s ~ Surf.dist$GeoDist.s + Condition(Surf.dist$BacDist.s + Surf.dist$EnviDist.s)))
##### Varpart #####

##### Plotting #####
HNFSurfDist.p <- Surf.dist %>%
  select(BacDist, HNFDist, EnviDist, GeoDist) %>%
  mutate(GeoDist = GeoDist/1000) %>%
  gather(key = "ExpVar", value = "Distance", -HNFDist) %>%
  mutate(ExpVar = factor(ExpVar, 
                         levels = c("BacDist", "EnviDist", "GeoDist"),
                         labels = c(expression(paste("Bacteria community dissimilarity")),
                                    expression(paste("Environmental dissimilarity")),
                                    expression(paste("Geographical distance (km)")) )
  )
  ) %>%
  ggplot(aes(x = Distance, y = HNFDist)) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  facet_grid(cols = vars(ExpVar), scales = "free", labeller = "label_parsed") +
  labs(x = expression(""),
       y = expression("HNF community dissimilarity")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 24)
  )
HNFSurfDist.p
ggsave(HNFSurfDist.p, file = "D:/Dropbox/Manuscripts/KuroshioMicrobes/MetaComm/Figs/HNFSurfDist.png",
       dpi = 600, width = 60, height = 32, units = "cm")

##### Plotting #####

#################################
##### Surface layer #############
#################################

#################################
##### DCM layer #################
#################################

##### prepare distance for calculation #####
DCM.list <- rownames(HNF.raw)[which(substr(rownames(HNF.raw), nchar(rownames(HNF.raw)), nchar(rownames(HNF.raw))) == "C")]

BacDCM <- Bac.raw[which(substr(rownames(Bac.raw), nchar(rownames(Bac.raw)), nchar(rownames(Bac.raw))) == "C"), ]
HNFDCM <- HNF.raw[which(substr(rownames(HNF.raw), nchar(rownames(HNF.raw)), nchar(rownames(HNF.raw))) == "C"), ]
BacDCM.dist <- matrix(as.matrix(vegdist(BacDCM, method = "bray")), ncol = 1)
HNFDCM.dist <- matrix(as.matrix(vegdist(HNFDCM, method = "bray")), ncol = 1)

EnviDCM <- MetaData %>%
  filter(SampleID %in% DCM.list) %>%
  select(Cruise, Station, Line, Longitude, Latitude, Temperature, Salinity, DIN, P, Chla) %>%
  mutate(Temp.s = as.numeric(scale(Temperature)),
         Sal.s = as.numeric(scale(Salinity)),
         TN.s = as.numeric(scale(DIN)),
         TP.s = as.numeric(scale(P)),
         Chla.s = as.numeric(scale(Chla))
  )
EnviDCM.dist <- matrix(as.matrix(vegdist(EnviDCM[, c("Temp.s", "Sal.s", "TN.s", "TP.s", "Chla.s")], method = "euclidean")), ncol = 1)
GDist <-     distm(EnviDCM[, c("Longitude", "Latitude")], fun = distGeo)
GeoDCM.dist <- matrix(distm(EnviDCM[, c("Longitude", "Latitude")], fun = distGeo), ncol = 1)

DCM.dist <- data.frame(BacDist = BacDCM.dist, 
                       HNFDist = HNFDCM.dist, 
                       EnviDist = EnviDCM.dist, 
                       GeoDist = GeoDCM.dist) %>%
  mutate(focal = rep(rownames(BacDCM), each = length(rownames(BacDCM))),
         to = rep(rownames(BacDCM), length(rownames(BacDCM))),
         transect = rep(EnviDCM$Line, each = length(EnviDCM$Line)),
         Cr = rep(EnviDCM$Line, each = length(EnviDCM$Line))) %>%
  filter(BacDist > 0) %>%
  mutate(BacDist.s = as.numeric(scale(BacDist)),
         HNFDist.s = as.numeric(scale(HNFDist)),
         EnviDist.s = as.numeric(scale(EnviDist)),
         GeoDist.s = as.numeric(scale(GeoDist))) %>%
  mutate(BacDist.Res = residuals(lm(BacDist.s ~ EnviDist.s + GeoDist.s)),
         HNFDist.Res = residuals(lm(HNFDist.s ~ EnviDist.s + GeoDist.s)),
         EnviDist.BacRes = residuals(lm(EnviDist.s ~ BacDist.s + GeoDist.s)),
         EnviDist.HNFRes = residuals(lm(EnviDist.s ~ HNFDist.s + GeoDist.s)),
         GeoDist.BacRes = residuals(lm(GeoDist.s ~ BacDist.s + EnviDist.s)),
         GeoDist.HNFRes = residuals(lm(GeoDist.s ~ HNFDist.s + EnviDist.s))
  )
##### prepare distance for calculation #####

########## for bacteria ##########
##### linear model #####
# BacDCM.lme <- lme(BacDist ~ HNFDist + EnviDist + GeoDist, random = ~ 1 | transect, method = "ML", 
#                    data = DCM.dist, na.action = na.omit)
# summary(BacDCM.lme)
# performance::r2(BacDCM.lme) 
BacDCM.lmeAll <- lm(BacDist.s ~ HNFDist.s + EnviDist.s + GeoDist.s,# random = ~ 1 | transect, method = "ML", 
                     data = DCM.dist, na.action = na.omit)
summary(BacDCM.lmeAll) 
performance::r2(BacDCM.lmeAll)

BacDCM.lmeHNF <- lm(BacDist.s ~ HNFDist.Res, #random = ~ 1 | transect, method = "ML", 
                     data = DCM.dist, na.action = na.omit)
summary(BacDCM.lmeHNF) 
performance::r2(BacDCM.lmeHNF)

BacDCM.lmeEnvi <- lm(BacDist.s ~ EnviDist.HNFRes, #random = ~ 1 | transect, method = "ML", 
                      data = DCM.dist, na.action = na.omit)
summary(BacDCM.lmeEnvi) 
performance::r2(BacDCM.lmeEnvi)

BacDCM.lmeGeo <- lm(BacDist.s ~ GeoDist.HNFRes, #random = ~ 1 | transect, method = "ML", 
                     data = DCM.dist, na.action = na.omit)
summary(BacDCM.lmeGeo) 
performance::r2(BacDCM.lmeGeo)
##### linear model #####

##### PERMANOVA #####
BacDCM.perMANOVA <- adonis(BacDist.s ~ HNFDist.s + EnviDist.s + GeoDist.s, data = DCM.dist, method = "euclidean")
#
BacDCM.perMANOVA
summary(BacDCM.perMANOVA)
##### PERMANOVA #####

##### Varpart #####
BacDCM.varpart <- varpart(DCM.dist$BacDist.s, DCM.dist$HNFDist.s, DCM.dist$EnviDist.s, DCM.dist$GeoDist.s)
BacDCM.varpart
rda(DCM.dist$BacDist.s ~ DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)
rda(DCM.dist$BacDist.s ~ DCM.dist$HNFDist.s + DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)
anova.cca(rda(DCM.dist$BacDist.s ~ DCM.dist$HNFDist.s + Condition(DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)))
anova.cca(rda(DCM.dist$BacDist.s ~ DCM.dist$EnviDist.s + Condition(DCM.dist$HNFDist.s + DCM.dist$GeoDist.s)))
anova.cca(rda(DCM.dist$BacDist.s ~ DCM.dist$GeoDist.s + Condition(DCM.dist$HNFDist.s + DCM.dist$EnviDist.s)))
##### Varpart #####

##### Plotting #####
BacDCMDist.p <- DCM.dist %>%
  select(BacDist, HNFDist, EnviDist, GeoDist) %>%
  mutate(GeoDist = GeoDist/1000) %>%
  gather(key = "ExpVar", value = "Distance", -BacDist) %>%
  mutate(ExpVar = factor(ExpVar, 
                         levels = c("HNFDist", "EnviDist", "GeoDist"),
                         labels = c(expression(paste("HNF community dissimilarity")),
                                    expression(paste("Environmental dissimilarity")),
                                    expression(paste("Geographical distance (km)")) )
  )
  ) %>%
  ggplot(aes(x = Distance, y = BacDist)) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  facet_grid(cols = vars(ExpVar), scales = "free", labeller = "label_parsed") +
  labs(x = expression(""),
       y = expression("Bacterial community dissimilarity")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 24)
  )
BacDCMDist.p
ggsave(BacDCMDist.p, file = "D:/Dropbox/Manuscripts/KuroshioMicrobes/MetaComm/Figs/BacDCMDist.png",
       dpi = 600, width = 60, height = 32, units = "cm")
##### Plotting #####
########## for bacteria ##########

########## for HNF ##########
##### linear model #####
# BacDCM.lme <- lme(BacDist ~ HNFDist + EnviDist + GeoDist, random = ~ 1 | transect, method = "ML", 
#                    data = DCM.dist, na.action = na.omit)
# summary(BacDCM.lme)
# performance::r2(BacDCM.lme) 
HNFDCM.lmeAll <- lm(HNFDist.s ~ BacDist.s + EnviDist.s + GeoDist.s,# random = ~ 1 | transect, method = "ML", 
                    data = DCM.dist, na.action = na.omit)
summary(HNFDCM.lmeAll) 
performance::r2(HNFDCM.lmeAll)

HNFDCM.lmeHNF <- lm(HNFDist.s ~ BacDist.Res, #random = ~ 1 | transect, method = "ML", 
                    data = DCM.dist, na.action = na.omit)
summary(HNFDCM.lmeHNF) 
performance::r2(HNFDCM.lmeHNF)

HNFDCM.lmeEnvi <- lm(HNFDist.s ~ EnviDist.BacRes, #random = ~ 1 | transect, method = "ML", 
                     data = DCM.dist, na.action = na.omit)
summary(HNFDCM.lmeEnvi) 
performance::r2(HNFDCM.lmeEnvi)

HNFDCM.lmeGeo <- lm(HNFDist.s ~ GeoDist.BacRes, #random = ~ 1 | transect, method = "ML", 
                    data = DCM.dist, na.action = na.omit)
summary(HNFDCM.lmeGeo) 
performance::r2(HNFDCM.lmeGeo)
##### linear model #####

##### PERMANOVA #####
HNFDCM.perMANOVA <- adonis(HNFDist.s ~ BacDist.s + EnviDist.s + GeoDist.s, data = DCM.dist, method = "euclidean")
#
HNFDCM.perMANOVA
summary(HNFDCM.perMANOVA)
##### PERMANOVA #####

##### Varpart #####
HNFDCM.varpart <- varpart(DCM.dist$HNFDist.s, DCM.dist$BacDist.s, DCM.dist$EnviDist.s, DCM.dist$GeoDist.s)
HNFDCM.varpart
rda(DCM.dist$HNFDist.s ~ DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)
rda(DCM.dist$HNFDist.s ~ DCM.dist$BacDist.s + DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)
anova.cca(rda(DCM.dist$HNFDist.s ~ DCM.dist$BacDist.s + Condition(DCM.dist$EnviDist.s + DCM.dist$GeoDist.s)))
anova.cca(rda(DCM.dist$HNFDist.s ~ DCM.dist$EnviDist.s + Condition(DCM.dist$HNFDist.s + DCM.dist$GeoDist.s)))
anova.cca(rda(DCM.dist$HNFDist.s ~ DCM.dist$GeoDist.s + Condition(DCM.dist$HNFDist.s + DCM.dist$EnviDist.s)))
##### Varpart #####

##### Plotting #####
HNFDCMDist.p <- DCM.dist %>%
  select(BacDist, HNFDist, EnviDist, GeoDist) %>%
  mutate(GeoDist = GeoDist/1000) %>%
  gather(key = "ExpVar", value = "Distance", -HNFDist) %>%
  mutate(ExpVar = factor(ExpVar, 
                         levels = c("BacDist", "EnviDist", "GeoDist"),
                         labels = c(expression(paste("HNF community dissimilarity")),
                                    expression(paste("Environmental dissimilarity")),
                                    expression(paste("Geographical distance (km)")) )
  )
  ) %>%
  ggplot(aes(x = Distance, y = HNFDist)) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  facet_grid(cols = vars(ExpVar), scales = "free", labeller = "label_parsed") +
  labs(x = expression(""),
       y = expression("Bacterial community dissimilarity")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 32),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 24)
  )
HNFDCMDist.p
ggsave(HNFDCMDist.p, file = "D:/Dropbox/Manuscripts/KuroshioMicrobes/MetaComm/Figs/HNFDCMDist.png",
       dpi = 600, width = 60, height = 32, units = "cm")
##### Plotting #####
########## for HNF ##########

#################################
##### DCM layer #################
#################################



###############################################################################################
##### Exploratory analyses ####################################################################
###############################################################################################
Corrs <- Kuro.ADiv %>%
  ggpairs(columns = c("Temperature", "Salinity", "DIN", "P", "Chla", "HNF_q0", "HNF_q1", "HNF_q2", "Bac_q0", "Bac_q1", "Bac_q2"),
          columnLabels = c("Temperature", "Salinity", "Total nitrogen", "Total phosphate", "Chlorophyll a",
                           "HNF richness", "HNF q1", "HNF q2", "Bac richness", "Bac q1", "Bac q2"),
          upper = list(continuous = cor_fun),
          lower = list(continuous = fit_fun)) +
  theme(strip.text.x = element_text(color = "black", size = 14),
        strip.text.y = element_text(angle = 45, color = "black", size = 14))
Corrs

lab.names <- as_labeller(c(ln.Temp = "Temperature~(degree*C)", 
                           Salinity = "Salinity~(PSU)",
                           DIN = "Total~nitrogen~(mu*M)", 
                           P = "Total~phosphate~(mu*M)", 
                           Chla = "Chlorophyll~a~concentration~(mg/m^3)",
                           HNF_q0 = "HNF~richness", 
                           HNF_q1 = "HNF~Shannon~diversity",
                           HNF_q2 = "HNF~Simpson~diversity",
                           Bac_q0 = "Bacteria~richness", 
                           Bac_q1 = "Bacteria~Shannon~diversity",
                           Bac_q2 = "Bacteria~Simpson~diversity"),
                         default = label_parsed)  
Kuro.ADiv %>%
  mutate(ln.Temp = log(Temperature)) %>%
  select(Cruise, Station, Line, Latitude, Temperature, Salinity, DIN, P, Chla, HNF_q0, HNF_q1, HNF_q2, Bac_q0, Bac_q1, Bac_q2) %>%
  arrange(Latitude) %>%
  gather(key = "Variable", value = "value", -c(Cruise, Station, Line, Latitude)) %>%
  mutate(Variable = factor(Variable, levels = c("Temperature", "Salinity", "DIN", "P", "Chla", "", 
                                                "HNF_q0", "HNF_q1", "HNF_q2", "Bac_q0", "Bac_q1", "Bac_q2")),
         Line = factor(Line, levels = c("TE", "OKTV", "OK", "TK", "TS", "BS"))) %>%
  ggplot(aes(x = Latitude)) + 
  geom_point(aes(y = value, color = Line), size = 3) +
  scale_color_viridis(alpha = 0.8, discrete = TRUE, begin = 1, end = 0) +
  facet_wrap(~ Variable, ncol = 3, scales = "free_y", drop = FALSE,
             labeller = lab.names) + 
  labs(y = expression("")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)
  )

Kuro.ADiv %>%
  mutate(ln.Temp = log(Temperature)) %>%
  select(Cruise, Station, Line, Latitude, ln.Temp, Salinity, DIN, P, Chla, HNF_q0, HNF_q1, HNF_q2, Bac_q0, Bac_q1, Bac_q2) %>%
  gather(key = "Variable", value = "value", -c(Cruise, Station, Line, Latitude)) %>%
  mutate(Variable = factor(Variable, levels = c("ln.Temp", "Salinity", "DIN", "P", "Chla", "", 
                                                "HNF_q0", "HNF_q1", "HNF_q2", "Bac_q0", "Bac_q1", "Bac_q2")),
         Line = factor(Line, levels = c("TE", "OKTV", "OK", "TK", "TS", "BS"))) %>%
  ggplot(aes(x = Line)) + 
  geom_boxplot(aes(y = value, fill = Line)) +
  scale_fill_viridis(alpha = 0.8, discrete = TRUE, begin = 1, end = 0) +
  facet_wrap(~ Variable, ncol = 3, scales = "free_y", drop = FALSE,
             labeller = lab.names) + 
  labs(y = expression("")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)
  )
###############################################################################################
##### Exploratory analyses ####################################################################
###############################################################################################