# R code for the analyses in Mena et al. in prep. Vertical stratification and horizontal diversity

#### Libraries ####
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(forcats)
library(iNEXT)
library(ggplot2)
library(vegan)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(SpadeR)
library(geosphere)
library(reshape2)
library(scales)
library(bestNormalize)

#### Diversity analyses ####

setwd("~/Overall diversity")
getwd()
datalalo<-read.delim("LaloDiversity.txt", header = TRUE)
lalodiv<- iNEXT(datalalo, q=c(0,1,2,3), datatype="abundance", size=NULL, 
                endpoint=1000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datalalo, datatype="abundance", base="coverage", level=0.99, conf=0.95)
lalodiv$DataInfo
lalodiv$iNextEst
lalodiv$AsyEst

datacanande<-read.delim("CanandeDiv.txt", header = TRUE)
canandediv<- iNEXT(datacanande, q=c(0,1,2,3), datatype="abundance", size=NULL, 
                   endpoint=1000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datacanande, datatype="abundance", base="coverage", level=0.99, conf=0.95)
canandediv$DataInfo
canandediv$iNextEst
canandediv$AsyEst

datadurango<-read.delim("DurangoDiv.txt", header = TRUE)
durangodiv<- iNEXT(datadurango, q=2, datatype="abundance", size=NULL, 
                   endpoint=1000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datadurango, datatype="abundance", base="coverage", level=0.99, conf=0.95)
durangodiv$DataInfo
durangodiv$iNextEst
durangodiv$AsyEst

datajorupe<-read.delim("JorupeDiv.txt", header = TRUE)
jorupediv<- iNEXT(datajorupe, q=c(0,1,2,3), datatype="abundance", size=NULL, 
                  endpoint=2000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datajorupe, datatype="abundance", base="coverage", level=0.99, conf=0.95)
jorupediv$DataInfo
jorupediv$iNextEst
jorupediv$AsyEst

datamashpi<-read.delim("MashpiDiv.txt", header = TRUE)
mashpidiv<- iNEXT(datamashpi, q=c(0,1,2,3), datatype="abundance", size=NULL, 
                  endpoint=2000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datamashpi, datatype="abundance", base="coverage", level=0.99, conf=0.95)
mashpidiv$DataInfo
mashpidiv$iNextEst
mashpidiv$AsyEst

datasumaco<-read.delim("SumacoDiv.txt", header = TRUE)
sumacodiv<- iNEXT(datasumaco, q=c(0,1,2,3), datatype="abundance", size=NULL, 
                  endpoint=5000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(datasumaco, datatype="abundance", base="coverage", level=0.99, conf=0.95)
sumacodiv$DataInfo
sumacodiv$iNextEst
sumacodiv$AsyEst

data<-read.delim("Yasuni_50hplot.txt", header = TRUE)
output<- iNEXT(data, q=c(0,1,2,3), datatype="abundance", size=NULL, 
               endpoint=15000, knots=40, se=TRUE, conf=0.95, nboot=1000)
estimateD(data, datatype="abundance", base="coverage", level=0.99, conf=0.95)
output$DataInfo
output$iNextEst
output$AsyEst

####### NMDS-ANOSIM #######

#NMDS 
#carrion
setwd("~/NMDS")
getwd()
pccc = read.csv("CanandeNMDSdd.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comcc = pccc[,3:ncol(pccc)]
#turn abundance data frame into a matrix
m_comcc = as.matrix(comcc)
nmdscc = metaMDS(m_comcc, distance = "bray")
nmdscc
#extract NMDS scores (x and y coordinates)
data.scorescc = as.data.frame(scores(nmdscc))
#add columns to data frame 
data.scorescc$sample = pccc$trap
data.scorescc$strata = pccc$strata
#to see first 6 data
head(data.scorescc)
#to plot

nmdsccplot<-ggplot()+
  geom_point(data = data.scorescc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",label="S=0.14", x=1, y=-0.8,  hjust = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anocc = anosim(m_comcc, pccc$strata, distance = "bray", permutations = 9999)
anocc

#banana
pccf = read.csv("CanandeNMDSddf.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comcf = pccf[,3:ncol(pccf)]
#turn abundance data frame into a matrix
m_comcf = as.matrix(comcf)
nmdscf = metaMDS(m_comcf, distance = "bray")
nmdscf
#extract NMDS scores (x and y coordinates)
data.scorescf = as.data.frame(scores(nmdscf))
#add columns to data frame 
data.scorescf$sample = pccf$trap
data.scorescf$strata = pccf$strata
#to see first 6 data
head(data.scorescf)
#to plot

nmdscfplot<-ggplot()+
  geom_point(data = data.scorescf, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",label="S=0.11", x=0.8, y=-0.5,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anocf = anosim(m_comcf, pccf$strata, distance = "bray", permutations = 9999)
anocf

#YASUNI 
#carrion
pcyc = read.csv("Yasuni_50hplot-carrion-NMDS.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comyc = pcyc[,3:ncol(pcyc)]
#turn abundance data frame into a matrix
m_comyc = as.matrix(comyc)
nmdsyc = metaMDS(m_comyc, distance = "bray")
nmdsyc
###used bray here, found identical results for jaccard and horn
#extract NMDS scores (x and y coordinates)
data.scoresyc = as.data.frame(scores(nmdsyc))
#add columns to data frame 
data.scoresyc$sample = pcyc$trap
data.scoresyc$strata = pcyc$strata
#to see first 6 data
head(data.scoresyc)
#to plot

nmdsycplot<-ggplot()+
  geom_point(data = data.scoresyc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+ 
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",label="S=0.12", x=0.8, y=-0.4,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anoyc = anosim(m_comyc, pcyc$strata, distance = "bray", permutations = 9999)
anoyc

#fruit

pcyf = read.csv("Yasuni_50hplot-fruit-NMDS.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comyf = pcyf[,3:ncol(pcyf)]
#turn abundance data frame into a matrix
m_comyf = as.matrix(comyf)
nmdsyf = metaMDS(m_comyf, distance = "bray")
nmdsyf
#extract NMDS scores (x and y coordinates)
data.scoresyf = as.data.frame(scores(nmdsyf))
#add columns to data frame 
data.scoresyf$sample = pcyf$trap
data.scoresyf$strata = pcyf$strata
#to see first 6 data
head(data.scoresyf)
#to plot

nmdsyfplot<-ggplot()+
  geom_point(data = data.scoresyf, aes(x=NMDS1, y=NMDS2, color=strata), size=5)+
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",label="S=0.15", x=1, y=-0.8,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anoyf = anosim(m_comyf, pcyf$strata, distance = "bray", permutations = 9999)
ano


#LALO
#carrion
pclc = read.csv("NMDSlc.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comlc = pclc[,3:ncol(pclc)]
#turn abundance data frame into a matrix
m_comlc = as.matrix(comlc)
nmdslc = metaMDS(m_comlc, distance = "bray")
nmdslc
#extract NMDS scores (x and y coordinates)
data.scoreslc = as.data.frame(scores(nmdslc))
#add columns to data frame 
data.scoreslc$sample = pclc$trap
data.scoreslc$strata = pclc$strata
#to see first 6 data
head(data.scoreslc)
#to plot
nmdslcplot<-ggplot()+
  geom_point(data = data.scoreslc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",label="S=0.14", x=1, y=-0.8,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anolc = anosim(m_comlc, pclc$strata, distance = "bray", permutations = 9999)
anolc

#NMDS 
#banana
pclf = read.csv("NMDSlf.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comlf = pclf[,3:ncol(pclf)]
#turn abundance data frame into a matrix
m_comlf = as.matrix(comlf)
nmdslf = metaMDS(m_comlf, distance = "bray")
nmdslf
#extract NMDS scores (x and y coordinates)
data.scoreslf = as.data.frame(scores(nmdslf))
#add columns to data frame 
data.scoreslf$sample = pclf$trap
data.scoreslf$strata = pclf$strata
#to see first 6 data
head(data.scoreslf)
#to plot

nmdslfplot<-ggplot()+
  geom_point(data = data.scoreslf, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  annotate("text",label="S=0.12", x=1, y=-0.6,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anolf = anosim(m_comlf, pclf$strata, distance = "bray", permutations = 9999)
anolf

#MASPI 
#carrion
pcmc = read.csv("MashpiNMDSc.csv", header= TRUE)
#make community matrix - extract columns with abundance information
commc = pcmc[,3:ncol(pcmc)]
#turn abundance data frame into a matrix
m_commc = as.matrix(commc)
nmdsmc = metaMDS(m_commc, distance = "bray")
nmdsmc
#extract NMDS scores (x and y coordinates)
data.scoresmc = as.data.frame(scores(nmdsmc))
#add columns to data frame 
data.scoresmc$sample = pcmc$trap
data.scoresmc$strata = pcmc$strata
#to see first 6 data
head(data.scoresmc)
#to plo

nmdsmcplot<-ggplot()+
  geom_point(data = data.scoresmc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  annotate("text",label="S=0.18", x=1.5, y=-1,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anomc = anosim(m_commc, pcmc$strata, distance = "bray", permutations = 9999)
anomc

#NMDS 
#banana-----update: working when drop canopy trap 8

pcmf1 <- read.csv("MashpiNMDSb2.csv", header= TRUE)
pcmf <- pcmf1[-8,]
#make community matrix - extract columns with abundance information
commf <- pcmf[,3:ncol(pcmf)]
#turn abundance data frame into a matrix
m_commf <- as.matrix(commf)
nmdsmf <- metaMDS(m_commf, distance = "bray", trymax=1000, k=2, weakties = TRUE)
nmdsmf
#extract NMDS scores (x and y coordinates)
data.scoresmf = as.data.frame(scores(nmdsmf))
#add columns to data frame 
data.scoresmf$sample = pcmf$trap
data.scoresmf$strata = pcmf$strata
#to see first 6 data
head(data.scoresmf)
#to plot

nmdsmfplot<-ggplot()+
  geom_point(data = data.scoresmf, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  annotate("text",label="S=0.12", x=0.8, y=-0.6,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anomf = anosim(m_commf, pcmf$strata, distance = "bray", permutations = 9999)
anomf


#SUMACO
#carrion
pcsc = read.csv("SumacoNMDSc.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comsc = pcsc[,3:ncol(pcsc)]
#turn abundance data frame into a matrix
m_comsc = as.matrix(comsc)
nmdssc = metaMDS(m_comsc, distance = "bray", trymax=1000)
nmdssc
#extract NMDS scores (x and y coordinates)
data.scoressc = as.data.frame(scores(nmdssc))
#add columns to data frame 
data.scoressc$sample = pcsc$trap
data.scoressc$strata = pcsc$strata
#to see first 6 data
head(data.scoressc)
#to plot

nmdsscplot<-ggplot()+
  geom_point(data = data.scoressc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank(),
        axis.ticks = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  annotate("text",label="S=0.17", x=1, y=-0.8,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anosc = anosim(m_comsc, pcsc$strata, distance = "bray", permutations = 9999)
anosc

#NMDS 
#banana
pcsf = read.csv("SumacoNMDSb.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comsf = pcsf[,3:ncol(pcsf)]
#turn abundance data frame into a matrix
m_comsf = as.matrix(comsf)
nmdssf = metaMDS(m_comsf, distance = "bray", trymax=1000)
nmdssf
#extract NMDS scores (x and y coordinates)
data.scoressf = as.data.frame(scores(nmdssf))
#add columns to data frame 
data.scoressf$sample = pcsf$trap
data.scoressf$strata = pcsf$strata
#to see first 6 data
head(data.scoressf)
#to plot

nmdssfplot<-ggplot()+
  geom_point(data = data.scoressf, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank(),axis.ticks = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  annotate("text",label="S=0.21", x=1.5, y=-1,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anosf = anosim(m_comsf, pcsf$strata, distance = "bray", permutations = 9999)
anosf

#DURANGO
#carrion
pcdc = read.csv("NMDSdurangoC.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comdc = pcdc[,3:ncol(pcdc)]
#turn abundance data frame into a matrix
m_comdc = as.matrix(comdc)
nmdsdc = metaMDS(m_comdc, distance = "bray", trymax=1000)
nmdsdc
#extract NMDS scores (x and y coordinates)
data.scoresdc = as.data.frame(scores(nmdsdc))
#add columns to data frame 
data.scoresdc$sample = pcdc$trap
data.scoresdc$strata = pcdc$strata
#to see first 6 data
head(data.scoresdc)
#to plot

nmdsdcplot<-ggplot()+
  geom_point(data = data.scoresdc, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank(),
        axis.ticks = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  annotate("text",label="S=0.21", x=1, y=-0.6,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anodc = anosim(m_comdc, pcdc$strata, distance = "bray", permutations = 9999)
anodc

#banana
pcdf = read.csv("NMDSdurangoF.csv", header= TRUE)
#make community matrix - extract columns with abundance information
comdf = pcdf[,3:ncol(pcdf)]
#turn abundance data frame into a matrix
m_comdf = as.matrix(comdf)
nmdsdf = metaMDS(m_comdf, distance = "bray", trymax=1000)
nmdsdf
#extract NMDS scores (x and y coordinates)
data.scoresdf = as.data.frame(scores(nmdsdf))
#add columns to data frame 
data.scoresdf$sample = pcdf$trap
data.scoresdf$strata = pcdf$strata
#to see first 6 data
head(data.scoresdf)
#to plot

nmdsdfplot <- ggplot()+
  geom_point(data = data.scoresdf, aes(x=NMDS1, y=NMDS2, color=strata),size=5)+ 
  scale_color_manual(values=c("royalblue4", "orange"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_blank(),
        axis.ticks = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  annotate("text",label="S=0.16", x=1.5, y=-1.2,  hjust = 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

anodf <- anosim(m_comdf, pcdf$strata, distance = "bray", permutations = 9999)
anodf

# NMDS all plots

nmdsycjld_cf<-ggarrange(nmdsycplot, 
                        nmdsccplot, 
                        nmdsscplot, 
                        nmdsjcplot, 
                        nmdslcplot, 
                        nmdsdcplot, 
                        nmdsmcplot,
                        nmdsyfplot, 
                        nmdscfplot, 
                        nmdssfplot, 
                        nmdsjfplot, 
                        nmdslfplot, 
                        nmdsdfplot,
                        nmdsmfplot,
                        common.legend = TRUE, legend = "bottom", 
                        #labels = c("Carrion", "Fruit"), 
                        ncol = 7, nrow = 2)
annotate_figure(nmdsycjld_cf, top = text_grob("", color = "black", 
                                              face = "bold", size = 14),
                left = text_grob("      Fruit                                        Carrion", color = "black", face="bold", rot = 90))+
  draw_plot_label(label = c("Yasuni", "Canande", "Sumaco", 
                            "Jorupe", "Lalo", "Durango", "Mashpi"), 
                  size = 14,
                  x = c( 0.04, 0.17, 0.32, 0.47, 0.62, 0.73, 0.88), 
                  y = c(0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99))


#### Horizontal diversity ####
setwd("~/Horizontal diversity")
#YASUNI
#q2
##canopy
yasunicc2<-read.csv("Yasuni_50hplot-canopy-carrion.csv", header = TRUE)
simmatrixyascc2<-SimilarityMult(yasunicc2, datatype = "abundance", q=2, nboot=200, "relative")
write.table(simmatrixyascc2$similarity.matrix,file="simmatixyasccq2.txt")
#q2 understorey
yasuniuc2<-read.csv("Yasuni_50hplot-understorey-carrion.csv", header = TRUE)
simmatrixyasuc2<-SimilarityMult(yasuniuc2, datatype = "abundance", q=2, nboot=200, "relative")
write.table(simmatrixyasuc2$similarity.matrix,file="simmatixyasucq2.txt")

#pairwises for fruit traps only
#q2canopy
yasunicf2<-read.csv("Yasuni_50hplot-canopy-fruit.csv", header = TRUE)
simmatrixyascf2<-SimilarityMult(yasunicf2, datatype = "abundance", q=2, nboot=200, "relative")
write.table(simmatrixyascf2$similarity.matrix,file="simmatixyascfq2.txt")

#q2 understorey
yasuniuf2<-read.csv("Yasuni_50hplot-understorey-fruit.csv", header = TRUE)
simmatrixyasuf2<-SimilarityMult(yasuniuf2, datatype = "abundance", q=2, nboot=200, "relative")
write.table(simmatrixyasuf2$similarity.matrix,file="simmatixyasufq2.txt")

#Canande
#Canande-canopy-carrion
canandecc0<-read.csv("CanandeCC.csv", header = TRUE)
simmatrixccc2<-SimilarityMult(canandecc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixccc2$pairwise
write.table(simmatrixccc2$pairwise,file="simmatixycccq2.txt")

#Canande-understorey-carrion
canandeuc0<-read.csv("CanandeUC.csv", header = TRUE)
simmatrixcuc2<-SimilarityMult(canandeuc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixcuc2$pairwise
write.table(simmatrixcuc2$pairwise,file="simmatixycucq2.txt")

#Canande-canopy-fruit
canandecf0<-read.csv("CanandeCF.csv", header = TRUE)
simmatrixccf2<-SimilarityMult(canandecf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixccf2$pairwise
write.table(simmatrixccf2$pairwise,file="simmatixyccfq2.txt")

#Canande-understorey-fruit
canandeuf0<-read.csv("CanandeUF.csv", header = TRUE)
simmatrixcuf2<-SimilarityMult(canandeuf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixcuf2$pairwise
write.table(simmatrixcuf2$pairwise,file="simmatixycufq2.txt")

#LALO
###lalo canopy carrion
lalocc0<-read.csv("laloCC.csv", header = TRUE)
simmatrixlcc2<-SimilarityMult(lalocc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixlcc2$pairwise
write.table(simmatrixlcc2$pairwise,file="simmatixlccq2.txt")

##lalo understorey carrion
lalouc0<-read.csv("LaloUC.csv", header = TRUE)
simmatrixluc2<-SimilarityMult(lalouc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixluc2$pairwise
write.table(simmatrixluc2$pairwise,file="simmatixlucq2.txt")

#lalo-canopy-fruit
lalocf0<-read.csv("LaloCF.csv", header = TRUE)
simmatrixlcf2<-SimilarityMult(lalocf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixlcf2$pairwise
write.table(simmatrixlcf2$pairwise,file="simmatixlcfq2.txt")

#lalo-understorey-fruit
lalouf0<-read.csv("LaloUF.csv", header = TRUE)
simmatrixluf2<-SimilarityMult(lalouf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixluf2$pairwise
write.table(simmatrixluf2$pairwise,file="simmatixlufq2.txt")

#JORUPE
###jorpue canopy carrion
jorupecc0<-read.csv("JorupeCC.csv", header = TRUE)
simmatrixjcc2<-SimilarityMult(jorupecc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixjcc2$pairwise
write.table(simmatrixjcc2$pairwise,file="simmatixjccq2.txt")

##jorupe understorey carrion
jorupeuc0<-read.csv("JorupeUC.csv", header = TRUE)
simmatrixjuc2<-SimilarityMult(jorupeuc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixjuc2$pairwise
write.table(simmatrixjuc2$pairwise,file="simmatixjucq2.txt")

#jorupe-canopy-fruit
jorupecf0<-read.csv("JorupeCF.csv", header = TRUE)
simmatrixjcf2<-SimilarityMult(jorupecf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixjcf2$pairwise
write.table(simmatrixjcf2$pairwise,file="simmatixjcfq2.txt")

#jorupe-understorey-fruit
jorupeuf0<-read.csv("JorupeUF.csv", header = TRUE)
simmatrixjuf2<-SimilarityMult(jorupeuf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixjuf2$pairwise
write.table(simmatrixjuf2$pairwise,file="simmatixjufq2.txt")

#SUMACO
#sumaco canopy carrion
sumacocc0<-read.csv("SumacoCC.csv", header = TRUE)
simmatrixscc2<-SimilarityMult(sumacocc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixscc2$pairwise
write.table(simmatrixscc2$pairwise,file="simmatixsccq2.txt")

##sumaco understorey carrion
sumacouc0<-read.csv("SumacoUC.csv", header = TRUE)
simmatrixsuc2<-SimilarityMult(sumacouc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixsuc2$pairwise
write.table(simmatrixsuc2$pairwise,file="simmatixsucq2.txt")

#sumaco-canopy-fruit
sumacocf0 <- read.csv("SumacoCF.csv", header = TRUE)
simmatrixscf2<-SimilarityMult(sumacocf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixscf2$pairwise
write.table(simmatrixscf2$pairwise,file="simmatixscfq2.txt")

#sumaco-understorey-fruit-
sumacouf0<-read.csv("SumacoUF.csv", header = TRUE)
simmatrixsuf2<-SimilarityMult(sumacouf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixsuf2$pairwise
write.table(simmatrixsuf2$pairwise,file="simmatixsufq2.txt")

#DURANGO
###durango canopy carrion
durangocc0<-read.csv("durangoCC.csv", header = TRUE)
simmatrixdcc2<-SimilarityMult(durangocc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixdcc2$pairwise
write.table(simmatrixdcc2$pairwise,file="simmatixdccq2.txt")

## understorey carrion
durangouc0<-read.csv("durangoUC.csv", header = TRUE)
simmatrixduc2<-SimilarityMult(durangouc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixduc2$pairwise
write.table(simmatrixduc2$pairwise,file="simmatixducq2.txt")

#durango-canopy-fruit--not working-> error in rmultinom negative probability-
#solucion 16 enero 2025- aumente un individuo para todas las filas y columnas de la matriz y funciono
#por eso el archivo utilizado fue CF2, no solo CF.
durangocf0<-read.csv("durangoCF2.csv", header = TRUE)
simmatrixdcf2<-SimilarityMult(durangocf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixdcf2$pairwise
write.table(simmatrixdcf2$pairwise,file="simmatixdcfq2.txt")

#durango-understorey-fruit
durangouf0<-read.csv("durangoUF.csv", header = TRUE)
simmatrixduf2<-SimilarityMult(durangouf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixduf2$pairwise
write.table(simmatrixduf2$pairwise,file="simmatixdufq2.txt")

#MASHPI
###mashpi canopy carrion
mashpicc0<-read.csv("mashpiCC.csv", header = TRUE)
simmatrixmcc2<-SimilarityMult(mashpicc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixmcc2$pairwise
write.table(simmatrixmcc2$pairwise,file="simmatixmccq2.txt")

## understorey carrion
mashpiuc0<-read.csv("mashpiCU.csv", header = TRUE)
simmatrixmuc2<-SimilarityMult(mashpiuc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixmuc2$pairwise
write.table(simmatrixmuc2$pairwise,file="simmatixmucq2.txt")

#mashpi-canopy-fruit
mashpicf0<-read.csv("mashpiCF2.csv", header = TRUE)
simmatrixmcf2<-SimilarityMult(mashpicf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixmcf2$pairwise
write.table(simmatrixmcf2$pairwise,file="simmatixmcfq2.txt")

#mashpi-understorey-fruit
mashpiuf0<-read.csv("mashpiUF2.csv", header = TRUE)
simmatrixmuf2<-SimilarityMult(mashpiuf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixmuf2$pairwise
write.table(simmatrixmuf2$pairwise,file="simmatixmufq2.txt")

#for OVERALL
#Overall-canopy-carrion
overallcc0<-read.csv("OverallCC.csv", header = TRUE)
simmatrixocc1<-SimilarityMult(overallcc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixocc1$estimated_richness
simmatrixocc1$estimated_relative

#Overall-understorey-carrion
overalluc0<-read.csv("OverallUC.csv", header = TRUE)
simmatrixouc1<-SimilarityMult(overalluc0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixouc1$estimated_richness
simmatrixouc1$estimated_relative

#Overall-canopy-fruit
overallcf0<-read.csv("OverallCF.csv", header = TRUE)
simmatrixocf1<-SimilarityMult(overallcf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixocf1$estimated_richness
simmatrixocf1$estimated_relative

#Overall-understorey-fruit
overalluf0<-read.csv("OverallUF.csv", header = TRUE)
simmatrixouf1<-SimilarityMult(overalluf0, datatype = "abundance", q=2, nboot=200, "relative")
simmatrixouf1$estimated_richness
simmatrixouf1$estimated_relative

#### Lollipops ####
# data
lolicarrion<-read.csv("lolipop.csv", header = TRUE)
summary(lolicarrion)

# Plot q2
lolq2<-ggplot(lolicarrion) +
  geom_point( aes(x=ï..x, y=q2c), color="royalblue4", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=c2LL, yend=c2UL), color="royalblue4", alpha=0.8) +
  geom_point( aes(x=ï..x, y=q2u), color="orange", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=u2LL, yend=u2UL), color="orange", alpha=0.8) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none",) +
  xlab("") +
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  #theme(axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,1))

# data
lolifruit<-read.csv("lolipopF.csv", header = TRUE)

# Plot q2
Flolq2<-ggplot(lolifruit) +
  geom_point( aes(x=ï..x, y=q2c), color="royalblue4", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=c2LL, yend=c2UL), color="royalblue4", 
                alpha=0.8) +
  geom_point( aes(x=ï..x, y=q2u), color="orange", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=u2LL, yend=u2UL), color="orange", 
                alpha=0.8) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none",) +
  xlab("") +
  #ylab("Simliarity (q2)")+
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  #theme(axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,1))

#allgraphs

loliqc2<-annotate_figure(Flolq2, top = text_grob(""))

#all plots (except overall)
#ggarrange(Lolic2,Lolif2,
#          labels=c("Carrion", "  Fruit"),
#          nrow = 2)
loliq2<-annotate_figure(loliq2, top = text_grob(""))

loliqc2<-annotate_figure(loliqc2, top = text_grob(""))


ggarrange(loliq2,loliqc2,
          labels=c("Carrion", "  Fruit"),
          nrow = 2)


# OVERALL
loliOcarrion <- read.csv("lolipopOverall.csv", header = TRUE)

# Plot q2
lolOq2<-ggplot(loliOcarrion) +
  geom_point( aes(x=ï..x, y=q2c), color="royalblue4", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=c2LL, yend=c2UL), color="royalblue4", alpha=0.8) +
  geom_point( aes(x=ï..x, y=q2u), color="orange", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=u2LL, yend=u2UL), color="orange", alpha=0.8) +
  coord_flip()+
  theme_bw() +
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab("Simliarity (q2)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  # theme(axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,0.5))

#allgraphs
library(ggplot2)
library(ggpubr)
library(cowplot)

LoliOc<-ggarrange(lolq2, lolOq2,
                  ncol = 1, nrow = 2)
LoliOc2<-annotate_figure(LoliOc)

# fruit
loliOfruit<-read.csv("lolipopOverallF.csv", header = TRUE)
# Plot q2
lolOq2f<-ggplot(loliOfruit) +
  geom_point( aes(x=ï..x, y=q2c), color="royalblue4", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=c2LL, yend=c2UL), color="royalblue4", 
                alpha=0.8) +
  geom_point( aes(x=ï..x, y=q2u), color="orange", size=3, alpha=1 ) +
  geom_segment( aes(x=ï..x, xend=ï..x, y=u2LL, yend=u2UL), color="orange", 
                alpha=0.8) +
  coord_flip()+
  theme_bw() +
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab("Simliarity (q2)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  # theme(axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,0.5))


# Add top margin to each plot
lolq2 <- lolq2 + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
Flolq2 <- Flolq2 + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
lolOq2 <- lolOq2 + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"))
lolOq2f <- lolOq2f + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"))

# Arrange the plots
lolipopsgood<-ggarrange(lolq2, Flolq2, lolOq2, lolOq2f,
                        heights = c(2, 0.7, 2, 0.7),
                        labels = c("Carrion", "Fruit", "", ""),
                        nrow = 2, ncol = 2)

# Create a blank plot with legend
legend_plot <- ggplot() +
  theme_void() +  # Blank theme
  annotate("point", x = 1, y = 1, color = "royalblue4", size = 3) +
  annotate("text", x = 1, y = 1, label = "Canopy", hjust = -.5) +
  annotate("point", x = 2, y = 1, color = "orange", size = 3) +
  annotate("text", x = 2, y = 1, label = "Understorey", hjust = -0.3)+
  xlim(-0.8, 4) + ylim(0.2, 1.5)

##put it together

ggarrange(lolipopsgood, legend_plot, nrow = 2, heights = c(15, 1))


#### Distance decay ####

setwd("~/Distance decay")

#Pairwise distances between traps
#YASUNI

#between carrion traps only
matrixy3c<-read.csv("matrix_y3-carrion.csv", header = TRUE)
distancesyasc<-distm(matrixy3c)
lower.tri(distancesyasc, diag = T)
dist <- melt(distancesyasc)
dist2 <- dist$value[-(which(dist$value == 0))]

write.table(distancesyasc,file="distancesyasc.txt")

#fruit traps only
matrixy3f<-read.csv("matrix_y3-fruit.csv", header = TRUE)
matrixy3f
distancesyasf<-distm(matrixy3f)
write.table(distancesyasf,file="distancesyasf.txt")

## CANANDE
#carrion
matrixcc<-read.csv("Puntos CanandeC.csv", header = TRUE)
matrixcc
distancescan<-distm(matrixcc)
write.table(distancescan,file="distancescc.txt")

#fruit
matrixcf<-read.csv("Puntos CanandeF.csv", header = TRUE)
matrixcf
distancescanf<-distm(matrixcf)
write.table(distancescanf,file="distancescf.txt")

##LALO
#carrion
matrixlc<-read.csv("Puntos LaloC.csv", header = TRUE)
matrixlc
distanceslal<-distm(matrixlc)
write.table(distanceslal,file="distanceslc.txt")

#fruit
matrixlf<-read.csv("Puntos LaloF.csv", header = TRUE)
matrixlf
distanceslalf<-distm(matrixlf)
write.table(distanceslalf,file="distanceslf.txt")

###JORUPE
#carrion
matrixjc<-read.csv("Points JorupeC.csv", header = TRUE)
matrixjc
distancesjorc<-distm(matrixjc)
write.table(distancesjorc,file="distancesjc.txt")
#fruit
matrixjf<-read.csv("Points JorupeF.csv", header = TRUE)
matrixjf
distancesjorf<-distm(matrixjf)
write.table(distancesjorf,file="distancesjf.txt")


####SUMACO
#carrion
matrixsc<-read.csv("Points SumacoC.csv", header = TRUE)
matrixsc
distancessumc<-distm(matrixsc)
write.table(distancessumc,file="distancessc.txt")
#fruit
matrixsf<-read.csv("Points SumacoF.csv", header = TRUE)
matrixsf
distancessumf<-distm(matrixsf)
write.table(distancessumf,file="distancessf.txt")

###DURANGO
#carrion
matrixdc<-read.csv("PointsDurangoC.csv", header = TRUE)
matrixdc
distancesdurc<-distm(matrixdc)
write.table(distancesdurc,file="distancesdc.txt")
#fruit
matrixdf<-read.csv("PointsDurangoF.csv", header = TRUE)
matrixdf
distancesdurf<-distm(matrixdf)
write.table(distancesdurf,file="distancesdf.txt")

#MASHPI
#all, because here carrion and fruit were shifted everyday so same coords for both
coord_mashpi<-read.csv("mashpi_traps.csv", header = TRUE)
coord_mashpi
distances_mashpi<-distm(subset(coord_mashpi, select = -trap))
write.table(distances_mashpi, file="distm1.csv")
#not working 4 some reason
dist_1 <- lower.tri(distances_mashpi, diag = T)
dist_m <- melt(dist_1)
dist_m <- dist_m$value[-(which(dist_m$value <1))]
write.table(dist_m,file="distances_mashpi.txt")

### PLOTS and MODELS
#YASUNI
ddyasuniq2 <- read.csv("preliminary analysisq2.csv", header = TRUE)

q2c <- ggplot(data=ddyasuniq2, aes(x=ï..distance, y=carrion_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", fill="royalblue4",
              alpha=0.2)+
  geom_smooth(data = ddyasuniq2, aes(x=ï..distance, y=carrion_u), 
              method="lm", se=TRUE, col="orange", fill="orange", alpha=0.2)+
  theme_bw()+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,600))

modelYq2c <- lm(carrion_c ~ ï..distance, data = ddyasuniq2)
summary(modelYq2c)$coefficients
modelYq2cu <- lm(carrion_u ~ ï..distance, data = ddyasuniq2)
summary(modelYq2cu)$coefficients

q2f <- ggplot(data=ddyasuniq2, aes(x=distancef, y=fruit_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4",  fill="royalblue4", 
              alpha=0.2)+
  geom_smooth(data = ddyasuniq2, aes(x=distancef, y=fruit_u), method="lm", 
              se=TRUE, col="orange",  fill="orange", alpha=0.2)+
  theme_bw()+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,600))

modelYq2f <- lm(fruit_c ~ ï..distance, data = ddyasuniq2)
summary(modelYq2f)$coefficients
modelYq2fu <- lm(fruit_u ~ ï..distance, data = ddyasuniq2)
summary(modelYq2fu)$coefficients

# LALO LOOR
ddlaloq2<-read.csv("lalodd2.csv", header = TRUE)

#plots
lq2c <- ggplot(data=ddlaloq2, aes(x=ï..distance, y=carrion_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", 
              fill="royalblue4", alpha=0.2, linetype = "dashed")+
  geom_smooth(data = ddlaloq2, aes(x=ï..distance, y=carrion_u), method="lm", 
              se=TRUE, col="orange", fill="orange", alpha=0.2, 
              linetype = "dashed")+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,500))

modelLLq2c <- lm(carrion_c ~ ï..distance, data = ddlaloq2)
summary(modelLLq2c)$coefficients
modelLLq2cu <- lm(carrion_u ~ ï..distance, data = ddlaloq2)
summary(modelLLq2cu)$coefficients

lq2f <- ggplot(data = ddlaloq2, aes(x = distancef, y = fruit_c))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4",  fill="royalblue4", 
              alpha=0.2, linetype = "dashed")+
  geom_smooth(data = ddlaloq2, aes(x=distancef, y=fruit_u), method = "lm", 
              se = TRUE, col = "orange",  fill = "orange", alpha = 0.2,
              fullrange=T)+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,650))

modelLLq2f <- lm(fruit_c ~ ï..distance, data = ddlaloq2)
summary(modelLLq2f)$coefficients
modelLLq2fu <- lm(fruit_u ~ ï..distance, data = ddlaloq2)
summary(modelLLq2fu)$coefficients

#JORUPE
ddjorupeq2<-read.csv("jorupedd2.csv", header = TRUE)

jq2c<-ggplot(data=ddjorupeq2, aes(x=ï..distance, y=carrion_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", 
              fill="royalblue4", alpha=0.2)+
  geom_smooth(data = ddjorupeq2, aes(x=ï..distance, y=carrion_u),
              method="lm", se=TRUE, col="orange", fill="orange", alpha=0.2)+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,700))

modelJq2f <- lm(carrion_c ~ ï..distance, data = ddjorupeq2)
summary(modelJq2f)$coefficients
modelJq2fu <- lm(carrion_u ~ ï..distance, data = ddjorupeq2)
summary(modelJq2fu)$coefficients

jq2f<-ggplot(data=ddjorupeq2, aes(x=distancef, y=fruit_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4",  fill="royalblue4", alpha=0.2)+
  geom_smooth(data = ddjorupeq2, aes(x=distancef, y=fruit_u), 
              method="lm", se=TRUE, col="orange",  fill="orange", alpha=0.2)+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,800))

modelJq2f <- lm(fruit_c ~ ï..distance, data = ddjorupeq2)
summary(modelJq2f)$coefficients
modelJq2fu <- lm(fruit_u ~ ï..distance, data = ddjorupeq2)
summary(modelJq2fu)$coefficients

#SUMACO
ddsumacoq2<-read.csv("sumacodd2.csv", header = TRUE)

sq2c<-ggplot(data=ddsumacoq2, aes(x= ï..distance, y=carrion_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", fill="royalblue4",
              alpha=0.2, linetype = "dashed")+
  geom_smooth(data = ddsumacoq2, aes(x= ï..distance, y=carrion_u), 
              method="lm", se=TRUE, col="orange", fill="orange", alpha=0.2)+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,1300))

modelSq2c <- lm(carrion_c ~ ï..distance, data = ddsumacoq2)
summary(modelSq2c)$coefficients
modelSq2cu <- lm(carrion_u ~ ï..distance, data = ddsumacoq2)
summary(modelSq2cu)$coefficients


sq2f<-ggplot(data=ddsumacoq2, aes(x=distancef, y=fruit_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", fill="royalblue4", 
              alpha=0.2, linetype = "dashed")+
  geom_smooth(data = ddsumacoq2, aes(x=distancef, y=fruit_u), method="lm", 
              se=TRUE, col="orange", fill="orange", alpha=0.2, 
              linetype = "dashed")+
  theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,1300))

sumaco<-ggarrange(sq0c, sq1c, sq2c,
                  ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
annotate_figure(sumaco, top = text_grob("Sumaco", color = "black", 
                                        face = "bold", size = 14),
                left = text_grob("Similarity", color = "black", 
                                 face="bold", rot = 90))+
  draw_plot_label(label = "Carrion", size = 10,
                  x = 0.1, y = 0.99)

#DURANGO
dddurangoq2<-read.csv("durangodd2.csv", header = TRUE)

dq2c <- ggplot(data = dddurangoq2, aes(x = ï..distance, y=carrion_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", fill="royalblue4", 
              alpha = 0.2)+
  geom_smooth(data = dddurangoq2, aes(x = ï..distance, y = carrion_u), method = "lm", 
              se = TRUE, col = "orange", fill = "orange", alpha = 0.2)+
  theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,1000))

modelDq2c <- lm(carrion_c ~ ï..distance, data = dddurangoq2)
summary(modelDq2c)$coefficients
modelDq2cu <- lm(carrion_u ~ ï..distance, data = dddurangoq2)
summary(modelDq2cu)$coefficients

dq2f <- ggplot(data = dddurangoq2, aes(x = ï..distance, y=fruit_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", fill="royalblue4", 
              alpha = 0.2)+
  geom_smooth(data = dddurangoq2, aes(x = ï..distance, y = fruit_u), method = "lm", 
              se = TRUE, col = "orange", fill = "orange", alpha = 0.2)+
  theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1.02))+
  scale_x_continuous(expand = c(0,0), breaks = c(100,1000))

modelDq2f <- lm(fruit_c ~ ï..distance, data = dddurangoq2)
summary(modelDq2f)$coefficients
modelDq2fu <- lm(fruit_u ~ ï..distance, data = dddurangoq2)
summary(modelDq2fu)$coefficients

durango<-ggarrange(dq0c, dq1c, dq2c, #labels = c("A", "B", "C", "D")
                   ncol = 1, nrow = 3, common.legend = TRUE, 
                   legend = "bottom")
annotate_figure(durango, top = text_grob("Durango", color = "black", 
                                         face = "bold", size = 14),
                left = text_grob("Similarity", color = "black", face="bold", 
                                 rot = 90))

#CANANDE
###q2
ddcanandeq2<-read.csv("canandedd2.csv", header = TRUE)

modelCq2c <- lm(carrion_c ~ ï..distance, data = ddcanandeq2)
summary(modelCq2c)$coefficients
modelCq2cu <- lm(carrion_u ~ ï..distance, data = ddcanandeq2)
summary(modelCq2cu)$coefficients

transf2<-orderNorm(ddcanandeq2$distancef)
ddcanandeq2$distancef<-transf2$x.t

norm3<-orderNorm(ddcanandeq2$ï..distance)
ddcanandeq2$ï..distance<-norm3$x.t

cq2c<-ggplot(data=ddcanandeq2, aes(x=ï..distance, y=carrion_u ))+ 
  geom_smooth(method = "lm", se = TRUE, col="royalblue4", 
              fill="royalblue4", alpha=0.2, fullrange = TRUE)+
  geom_smooth(data = ddcanandeq2, aes(x=ï..distance, y=carrion_c), 
              method="lm", se=TRUE, col="orange", fill="orange", 
              alpha=0.2, linetype = "dashed", fullrange = TRUE)+
  theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(-2,2), 
                     labels = c("100", "1000"))

cq2f<-ggplot(data = ddcanandeq2, aes(x = distancef, y = fruit_c ))+ 
  geom_smooth(method = "lm", se = TRUE, col = "royalblue4", 
              fill = "royalblue4", alpha = 0.2, linetype = "dashed")+
  geom_smooth(data = ddcanandeq2, aes(x = distancef, y = fruit_u), 
              method = "lm", se = TRUE, col = "orange", fill = "orange", 
              alpha = 0.2, linetype = "dashed")+
  theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))+
  scale_y_continuous(labels = label_number(accuracy = 1), breaks = c(0,1), 
                     expand = c(0,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0,0), breaks = c(-2,2), 
                     labels = c("100", "1000"))

modelCq2f <- lm(fruit_c ~ ï..distance, data = ddcanandeq2)
summary(modelCq2f)$coefficients
modelCq2fu <- lm(fruit_u ~ ï..distance, data = ddcanandeq2)
summary(modelCq2fu)$coefficients

###ALL PLOTS
#########carrion all plots

allcarrion <- ggarrange(lq2c, q2c, dq2c, jq2c, sq2c, cq2c,
                        ncol=6, nrow=1)
carrionfig <- annotate_figure(allcarrion,left= text_grob("Similarity q2",
                                                         color = "black", 
                                                         rot = 90), top = text_grob(" "), 
                              bottom = text_grob("Distance (m)"))+
  draw_plot_label(label = c("Lalo loor", "Yasuni", "Durango", 
                            "Jorupe", "Sumaco", "Canande"), size = 12,
                  x = c( 0.05, 0.22, 0.38, 0.55, 0.71, 0.86), 
                  y = c(0.99, 0.99, 0.99, 0.99, 0.99,0.99))

carrion3<-annotate_figure(carrionfig, fig.lab = "Carrion", 
                          fig.lab.face = "bold", 
                          fig.lab.size = 14,
                          fig.lab.pos = "top.left")

######## fruit all plots

empplot<-ggplot()

allfruit <- ggarrange(lq2f, q2f, 
                      dq2f, 
                      jq2f, 
                      sq2f,
                      cq2f,
                      ncol=6, 
                      nrow=1)

fruitfig<-annotate_figure(allfruit,
                          left= text_grob("Similarity q2",
                                          color = "black", rot = 90), 
                          top = text_grob(" "), 
                          bottom = text_grob("Distance (m)"))+
  draw_plot_label(label = c("Lalo loor", "Yasuni", "Durango", 
                            "Jorupe", "Sumaco", "Canande"), size = 12,
                  x = c( 0.05, 0.22, 0.38, 0.55, 0.71, 0.86), 
                  y = c(0.99, 0.99, 0.99, 0.99, 0.99,0.99))

fruit3<-annotate_figure(fruitfig, fig.lab = "Fruit", 
                        fig.lab.face = "bold", 
                        fig.lab.size = 14,
                        fig.lab.pos = "top.left")

#both carrion and fruit
dd_car_fru<-ggarrange(carrion3, fruit3, ncol = 1, nrow = 2)

legend_plot <- ggplot() +
  theme_void() +  # Blank theme
  annotate("segment", x = 1, xend = 1.1, y = 1, yend = 1, 
           color = "royalblue4", size = 1.2) +
  annotate("text", x = 1, y = 1, label = "Canopy", hjust = -.6) +
  annotate("segment", x = 2, xend = 2.1, y = 1, yend = 1, 
           color = "orange", size = 1.2) +
  annotate("text", x = 2, y = 1, label = "Understorey", hjust = -0.4)+
  xlim(-0.8, 4) + ylim(0.2, 1.5)

##put it together
ggarrange(dd_car_fru, legend_plot, nrow = 2, heights = c(15, 1))


