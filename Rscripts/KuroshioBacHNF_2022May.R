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

HNF.phylo<- read.tree(file = "D:/Dropbox/Research/KuroshioMicrobes/data/KuroHNF_rare_treeNJ.tree")

Bac.raw <- as.data.frame(t(read.table(
  file = "D:/Dropbox/Research/KuroshioMicrobes/data/KuroBac_rareSeqXSt.csv", 
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac.ra <- Bac.raw / rowSums(Bac.raw)

Bac.phylo<- read.tree(file = "D:/Dropbox/Research/KuroshioMicrobes/data/KuroBac_Cyano_rare_treeNJ.tree")

MetaData <- read.table(file = "D:/Dropbox/Research/KuroshioMicrobes/data/Kuroshio_metadata.csv", sep = ",", 
                       header = TRUE, stringsAsFactors = FALSE, fill = TRUE)


HNF.ADiv <- iNEXT(t(HNF.raw), q = 0, datatype = "abundance", size = max(colSums(HNF.raw)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(HNF.raw)) %>%
  inner_join(MetaData, by = c("Site" = "SampleID"))

Bac.ADiv <- iNEXT(t(Bac.raw), q = 0, datatype = "abundance", size = max(colSums(Bac.raw)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(Bac.raw))

Kuro.ADiv <- HNF.ADiv %>%
  inner_join(Bac.ADiv, by = c("Site" = "Site"))
  
###############################################################################################
##### Loading and preping data ################################################################
###############################################################################################

which(substr(rownames(Bac.raw), nchar(rownames(Bac.raw)), nchar(rownames(Bac.raw))) == "S")

BacSurf <- Bac.raw[which(substr(rownames(Bac.raw), nchar(rownames(Bac.raw)), nchar(rownames(Bac.raw))) == "S"), ]
HNFSurf <- HNF.raw[which(substr(rownames(HNF.raw), nchar(rownames(HNF.raw)), nchar(rownames(HNF.raw))) == "S"), ]
BacSurf.dist <- matrix(as.matrix(vegdist(BacSurf, method = "bray")), ncol = 1)
HNFSurf.dist <- matrix(as.matrix(vegdist(HNFSurf, method = "bray")), ncol = 1)

EnviS <- MetaData %>%
  filter(Layer == "Surface") %>%
  select(Cruise, Station, Line, Longitude, Latitude, Temperature, Salinity, DIN, P, Chla) %>%
  mutate(Temp.s = scale(Temperature),
         Sal.s = scale(Salinity),
         TN.s = scale(DIN),
         TP.s = scale(P),
         Chla.s = scale(Chla)
         )
EnviS.dist <- matrix(as.matrix(vegdist(EnviS[, c("Temp.s", "Sal.s", "TN.s", "TP.s", "Chla.s")], method = "euclidean")), ncol = 1)
#GDist <- matrix(distm(EnviS[, c("Longitude", "Latitude")], fun = distGeo), ncol = 1)
GeoS.dist <- matrix(distm(EnviS[, c("Longitude", "Latitude")], fun = distGeo), ncol = 1)

Surf.dist <- data.frame(BacDist = BacSurf.dist,
                        HNFDist = HNFSurf.dist,
                        EnviDist = EnviS.dist, 
                        GeoDist = GeoS.dist) %>%
  filter(BacDist > 0) %>%
  mutate(GeoDist = as.numeric(scale(GeoDist)))

BacDist.p <- Surf.dist %>%
  gather(key = "ExpVar", value = "Distance", -BacDist) %>%
  mutate(ExpVar = factor(ExpVar, 
                         levels = c("HNFDist", "EnviDist", "GeoDist"),
                         labels = c(expression(paste("HNF community dissimilarity")),
                                    expression(paste("Environmental dissimilarity")),
                                    expression(paste("Geographical distance")) )
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
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    strip.text.x = element_text(size = 16)
  )
BacDist.p
ggsave(BacDist.p, file = "D:/Dropbox/Research/KuroshioMicrobes/Figs/BacDist_HNF_Envi_Dist.png",
       dpi = 600, width = 32, height = 16, units = "cm")

summary(lm(BacDist ~ HNFDist + EnviDist + GeoDist, data = Surf.dist))
varpart(Surf.dist$BacDist, Surf.dist$HNFDist, Surf.dist$EnviDist, Surf.dist$GeoDist)

varpart(BacSurf.dist, HNFSurf.dist, EnviS.dist, GeoS.dist)
plot(varpart(BacSurf.dist, HNFSurf.dist, EnviS.dist, GeoS.dist))

HNFrda <- rda(BacSurf.dist ~ HNFSurf.dist + Condition(EnviS.dist) + Condition(GeoS.dist))
anova(HNFrda, step=200, perm.max=200)


MetaData %>%
  ggplot() +
  geom_point(aes(x = Latitude, y = Temperature, color = Line))

cor.test(Kuro.ADiv$DIN, Kuro.ADiv$P)

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


summary(lm(HNF_q0 ~ Latitude, data = Kuro.ADiv))
summary(lm(HNF_q1 ~ Latitude, data = Kuro.ADiv))
summary(lm(HNF_q2 ~ Latitude, data = Kuro.ADiv))

summary(lme(HNF_q0 ~ Temperature + Salinity + DIN + P + Chla + Bac_q0, random = ~ 1 | Line, data = Kuro.ADiv, method = "ML"))
summary(lme(HNF_q0 ~ Temperature + Salinity + DIN + Chla + Bac_q0, random = ~ 1 | Line, data = Kuro.ADiv, method = "ML"))
summary(lme(HNF_q0 ~ Temperature + DIN + Chla + Bac_q0, random = ~ 1 | Line, data = Kuro.ADiv, method = "ML"))
summary(lme(HNF_q0 ~ DIN + Chla + Bac_q0, random = ~ 1 | Line, data = Kuro.ADiv, method = "ML"))
summary(lme(HNF_q0 ~ DIN + Bac_q0, random = ~ 1 | Line, data = Kuro.ADiv, method = "ML"))


Surf.p <- Kuro.ADiv %>%
  filter(Layer == "Surface") %>%
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
  labs(y = expression(""),
       title = "Surface") + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 20)
  )
Surf.p













Bac_A <- iNEXT(t(Bac_comm), q = 0, datatype = "abundance", size = max(colSums(Bac_comm)) + 100000)$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(Bac_q0 = "Species richness", Bac_q1 = "Shannon diversity", Bac_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(Bac_comm))

HNF_A <- iNEXT(t(HNF_comm), q = 0, datatype = "abundance", size = max(colSums(HNF_comm)))$AsyEst %>% 
  select(Site, Diversity, Estimator) %>% 
  spread(Diversity, Estimator) %>%
  rename(HNF_q0 = "Species richness", HNF_q1 = "Shannon diversity", HNF_q2 = "Simpson diversity") %>%
  mutate(Site = rownames(HNF_comm))

HNF_Bac_A <- Bac_A %>%
  inner_join(Bac_Ampd, by = c("Site" = "Site")) %>%
  inner_join(HNF_A, by = c("Site" = "Site")) %>%
  inner_join(HNF_Ampd, by = c("Site" = "Site")) %>%
  inner_join(Vars, by = c("Site" = "SampleID")) %>%
  filter(!is.na(NF_Biom)) %>%
  mutate(ln.Bac_q0 = log(Bac_q0),
         ln.HNF_q0 = log(HNF_q0),
         ln.Bac_q1 = log(Bac_q1),
         ln.HNF_q1 = log(HNF_q1),
         ln.Bac_q2 = log(Bac_q2),
         ln.HNF_q2 = log(HNF_q2),
         ln.Bac_Biom = log(Bac_Biom),
         ln.HNF_Biom = log(HNF_Biom),
         ln.Temp = log(Temp),
         ln.Sal = log(Sal),
         ln.PAR = log(PAR),
         ln.NO2 = log(NO2 + 0.0001),
         ln.NO3 = log(NO3 + 0.0001),
         ln.DIN = log(DIN + 0.0001),
         ln.PO3 = log(PO3 + 0.0001), 
         ln.Chla = log(Chla + 0.00001),
         Month = as.character(Month))
HNF_Bac_A <- as.data.frame(HNF_Bac_A)
head(HNF_Bac_A)

### Prepare spatial autocorrelation among stations

### Build build Moran's Eigenvector Maps (MEM, Dray, Legendre, and Peres-Neto (2006)) 
### that are orthogonal vectors maximizing the spatial autocorrelation (measured by Moran's coefficient)
### Ref: https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#selection-of-swm-and-mem
library(geodist)
Bac.raw <- as.data.frame(t(read.table(
  file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/sECS_4/sECS_Bac_seqXst_PR2_4.csv",
  sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE, fill = TRUE)))
Bac.comm <- Bac.raw %>%
  mutate(Cruise = substr(row.names(Bac.raw), 4, 9),
         Station = substr(row.names(Bac.raw), 10, 13))
St_sECS = read.table(file = "https://raw.githubusercontent.com/OscarFHC/PdPy_Div/master/data/raw/St_Location.csv",
                     header = T, fill = T, sep = ",")[1:6, ] %>% select(-Station)
row.names(St_sECS) <- unique(Vars$Station)
geo <- as.matrix(geodist(St_sECS, measure = "geodesic")) / 1000 # geographical distance in kilo-meter
diag(geo) <- NA

Stxy <- as.matrix(St_sECS)

Bac.avg <- matrix(0, length(unique(Bac.comm$Station)), ncol(Bac.raw))
for (i in 1:length(unique(Bac.comm$Station))){
  Bac.temp <- Bac.comm[Bac.comm$Station == unique(Bac.comm$Station)[i], -which(names(Bac.comm) %in% c("Cruise", "Station"))]
  Bac.avg[i,] <- colMeans(Bac.temp)
}
# 
# listw.cand <- listw.candidates(coord = coordinates(Stxy), style = "W", nb = "dnear", d1 = 0, d2 = max(geo[which(geo != "NA")]), 
#                                weight = c("flin", "fup", "fdown"), y_fup = seq(0.1, 1, by = 0.01), y_fdown = seq(1, 10, by = 0.1))
# listw.sel <- listw.select(Bac.avg, listw.cand, MEM.autocor = "all", method = "FWD", 
#                           p.adjust = TRUE, nperm = 999)
# listw.cand[which(listw.sel$candidates$R2Adj == max(listw.sel$candidates$R2Adj))]

#listw.explore()

library(adespatial);library(sp);library(spdep)
nb <- chooseCN(coordinates(Stxy), type = 5, d1 = 0, d2 = max(geo[which(geo != "NA")]), plot.nb = FALSE)
distnb <- nbdists(nb, coordinates(Stxy), longlat = TRUE)
fdist <- lapply(distnb, function(x) 1 / x^0.22)
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)
mem.select(Bac.avg, lw)
mem.fup <- mem(lw)

MEM <- data.frame(Spatial.corr = mem.fup$MEM1) %>%
  mutate(Station = unique(Vars$Station))

HNF_Bac_A <- HNF_Bac_A %>% 
  left_join(MEM, by = c("Station" = "Station"))

###############################################################################################
##### Preping data ############################################################################
###############################################################################################