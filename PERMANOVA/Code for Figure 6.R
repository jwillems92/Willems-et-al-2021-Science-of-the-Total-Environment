library(tidyverse)
library(vegan)
library(vegan3d)
library(multcomp)
library(lme4)
library(performance)
library(MuMIn)
library(blmeco)
library(matrixStats)
library(ecodist)
library(MASS)
library(rgl)
library(ggeffects)
library(patchwork)
options(na.action = "na.fail")



d <- read.csv("Data/CamSppBin2.csv", header=TRUE)
head(d)



# how many zeros per treatment
d_zero <- subset(d, d$rich < 1)     # Very few zeroes now, only 4 cameras
d_zero

# cameras with zero detections cannot be included
d2 <- subset(d, d$rich > 0)

# how many non zero cameras are left per treatment?
sapply(split(d2$Treatment, d2$Treatment), length) # slightly fewer cameras at combined and light treatments


rankindex(d2$Lux, d2[2:24],indices = c("euclidean", "manhattan", "gower", "altGower", "canberra", "clark", "kulczynski", 
                                       "horn", "binomial", "jaccard", "bray"))

# for lux: euclidean, manhattan, gower, binomial, altGower, canberra
# for laeq: kulczynski, altGower, canberra, clark, horn, jaccard, bray
# for treatment: euclidean, manhattan, gower, binomial, jaccard


betaD1 <-vegdist(d2[2:24], method ="gower", binary=TRUE) # for Treatment model

betaD2 <-vegdist(d2[2:24], method ="altGower", binary=TRUE) # for Lux/LAeq model






# Separate models

adonis2(betaD2 ~ Lux*Laeq +
          Lux + 
          Laeq +
          mIndex + 
          Landcover + 
          JDay, 
        strata = d2$Cluster, 
        by = "margin", data = d2) 


adonis2(betaD2 ~ Lux + 
          Laeq +
          mIndex + 
          Landcover + 
          JDay, 
        strata = d2$Cluster, 
        by = "margin", data = d2) 


adonis2(betaD1 ~ Treatment +
          mIndex + 
          Landcover + 
          JDay, 
        strata = d2$Cluster, 
        by = "margin", data = d2) 



# FIGURES


leqD <- vegdist(d2$Laeq, method = "altGower")
luxD <- vegdist(d2$Lux, method = "altGower")
mIndexD <- vegdist(d2$mIndex, method = "altGower")
JDayD <- dist(d2$JDay)
betaDv <- c(betaD2) 

envD <- data.frame(c(leqD),c(luxD),c(JDayD), c(mIndexD), c(betaDv))
names(envD) <- c("LeqD", "LuxD", "JDayD", "MoonIndexD", "betaDdf")

betaDMod <- lm(betaDdf ~ LeqD* LuxD + JDayD + MoonIndexD, data = envD)

summary(betaDMod)



# Noise
DFsdB <- ggpredict(betaDMod, terms = c("LeqD", "LuxD"), back.transform =F, type ="fe")

C1 <- plot(DFsdB, colors="eight", add.data=F)  + 
        labs(x="Sound difference (dB)", y ="Dissimilarity", title ="") + 
        theme_classic(base_size = 18) + 
        theme(panel.grid.minor = element_blank()) + 
        theme(axis.line = element_line(colour = "black")) +
        theme(legend.position = "none") +
        labs(subtitle = "(a)") +
        theme(plot.subtitle = element_text(size = 20, face = "bold")) +
        theme(plot.subtitle = element_text(hjust = -0.13))
C1


# Light
DFsLux <- ggpredict(betaDMod, terms = c("LuxD", "LeqD"), back.transform =F, type ="fe")

D1 <- plot(DFsLux, colors="eight", add.data=F)  + 
        labs(x="Light difference (Lux)", y ="Dissimilarity", title ="") + 
        theme_classic(base_size = 18) + 
        theme(panel.grid.minor = element_blank()) + 
        theme(axis.line = element_line(colour = "black")) +
        theme(legend.position = "none") +
        labs(subtitle = "(b)") +
        theme(plot.subtitle = element_text(size = 20, face = "bold")) +
        theme(plot.subtitle = element_text(hjust = -0.13))

D1


# Moon
DFsMoon <- ggpredict(betaDMod, terms = c("MoonIndexD", "LeqD"), back.transform =F, type ="fe")

M1 <- plot(DFsMoon, colors="eight", add.data=F)  + 
        labs(x="Moonlight difference", y ="Dissimilarity", title ="") + 
        theme_classic(base_size = 18) + 
        theme(panel.grid.minor = element_blank()) + 
        theme(axis.line = element_line(colour = "black")) +
        theme(legend.position = "none") +
        labs(subtitle = "(c)") +
        theme(plot.subtitle = element_text(size = 20, face = "bold")) +
        theme(plot.subtitle = element_text(hjust = -0.13))
M1


# All combined
C1 / D1 / M1
