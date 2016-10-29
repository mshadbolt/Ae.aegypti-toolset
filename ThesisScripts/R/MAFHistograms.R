source("scripts/Functions.R")

# MAF histograms
mycols=c("#117733", "#88CCEE", "#332288", "#882255", "#40AA99")
# Load MAF outputs from vcftools
VietColMaf<-read.table("inputs/vcftoolsoutputs/VietCol27.frq", row.names = NULL)
VietWiMaf<-read.table("inputs/vcftoolsoutputs/VietWild27.frq", row.names=NULL)
AusWiViMaf<-read.table("inputs/vcftoolsoutputs/AusWildVi27.frq", row.names=NULL)

# Merge into a single dataframe
AusVietMaf<-data.frame(AusWiViMaf[,c(1,2,6)], VietColMaf[,6], VietWiMaf[,6])
colnames(AusVietMaf)<-c("CHROM", "POS", "AusMAF", "VietColMAF", "VietWiMAF")

# Foldover MAFvalues
AusVietMaf$AusWiFOLD<-sapply(AusVietMaf$AusMAF, minmaj)
AusVietMaf$VietWiFOLD<-sapply(AusVietMaf$VietWiMAF, minmaj)
AusVietMaf$VietColFOLD<-sapply(AusVietMaf$VietColMAF, minmaj)

AusViHists<-reshape(AusVietMaf,
                    varying = c("AusWiFOLD", "VietColFOLD", "VietWiFOLD"),
                    v.names = "MAF", 
                    timevar = "pop",
                    times = c("AusWiFOLD", "VietColFOLD", "VietWiFOLD"),
                    direction="long")


library(ggplot2)
labels = c(AusWiFOLD="AUwild", VietColFOLD="VNcol", VietWiFOLD="VNwild")
jpeg("AusViMAFHists.jpeg", 6,6, res=300, units="in")
AusViHist<-ggplot(data=AusViHists, aes(x=MAF, y=..density.., color=NULL, fill = as.factor(pop))) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(~pop, ncol = 1, labeller = labeller(pop=labels)) + 
  theme_bw() + 
  scale_fill_manual(values = c("#332288", "#40AA99",  "#882255")) + 
  scale_color_manual(values = c( "#332288", "#40AA99",  "#882255")) + 
  scale_x_continuous("Minor Allele Frequency") + 
  scale_y_continuous("Density") + 
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(fill="white"))
AusViHist
dev.off()

jpeg("AusViMAFSFS.jpeg", 9,6, res=300, units="in")
AusViSFS<-ggplot(data=AusViHists, aes(x=MAF, y=..density.., color=NULL, fill = as.factor(pop))) +
  geom_histogram(binwidth = 0.02, position="dodge") +
  theme_bw() + 
  scale_fill_manual(values = c("#332288", "#40AA99",  "#882255")) + 
  scale_color_manual(values = c( "#332288", "#40AA99",  "#882255")) + 
  scale_x_continuous("Minor Allele Frequency") + 
  scale_y_continuous("Density") + 
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(fill="white"))
AusViSFS
dev.off()


## Aus Braz
mycols=c("#117733", "#88CCEE", "#332288", "#882255", "#40AA99") 
AusBrMaf<-read.table("inputs/vcftoolsoutputs/AusWildBrz-20.frq", row.names = NULL)
BrColKD<-read.table("inputs/vcftoolsoutputs/BrazColKD-20.frq", row.names=NULL)
BrColW<-read.table("inputs/vcftoolsoutputs/BrazColW-20.frq", row.names=NULL)
BrWild<-read.table("inputs/vcftoolsoutputs/BrazWild-20.frq", row.names=NULL)

AusBrazMaf<-data.frame(AusBrMaf[,c(1,2,6)], BrColKD[,6], BrColW[,6], BrWild[,6])
colnames(AusBrazMaf)<-c("CHROM", "POS", "AusMAF", "KDColMAF", "WColMAF", "BrWildMAF")

AusBrazMaf$AusWiFOLD<-sapply(AusBrazMaf$AusMAF, minmaj)
AusBrazMaf$BrazWiFOLD<-sapply(AusBrazMaf$BrWildMAF, minmaj)
AusBrazMaf$KDColFOLD<-sapply(AusBrazMaf$KDColMAF, minmaj)
AusBrazMaf$WColFOLD<-sapply(AusBrazMaf$WColMAF, minmaj)


AusBrHists<-reshape(AusBrazMaf,
                    varying = c("AusWiFOLD", "BrazWiFOLD", "KDColFOLD", "WColFOLD"),
                    v.names = "MAF", 
                    timevar = "pop",
                    times = c("AusWiFOLD", "BrazWiFOLD", "KDColFOLD", "WColFOLD"),
                    direction="long")
AusBrHists$pop<-factor(AusBrHists$pop, levels = c("AusWiFOLD", "WColFOLD", "KDColFOLD", "BrazWiFOLD"))

library(ggplot2)
labels = c(AusWiFOLD="AUwild", KDColFOLD="BRcol2", WColFOLD="BRcol1", BrazWiFOLD="BRwild")
jpeg("AusBrMAFHists.jpeg", 6,8, res=300, units="in")
AusBrHist<-ggplot(data=AusBrHists, aes(x=MAF, y=..density.., color=NULL, fill = as.factor(pop))) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(~pop, ncol = 1, labeller = labeller(pop=labels)) + 
  theme_bw() + 
  scale_fill_manual(values = c("#332288",   "#117733", "#40AA99","#882255")) + 
  scale_x_continuous("Minor Allele Frequency") + 
  scale_y_continuous("Density") + 
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(fill="white"))
AusBrHist
dev.off()

jpeg("AusBrMAFSFS.jpeg", 9,6, res=300, units="in")
AusBrSFS<-ggplot(data=AusBrHists, aes(x=MAF, y=..density.., color=NULL, fill = as.factor(pop))) +
  geom_histogram(binwidth = 0.02, position="dodge") +
  theme_bw() + 
  scale_fill_manual(values = c("#332288", "#40AA99",  "#117733", "#882255")) + 
  scale_x_continuous("Minor Allele Frequency") + 
  scale_y_continuous("Density") + 
  theme(legend.position="none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(fill="white"))
AusBrSFS
dev.off()