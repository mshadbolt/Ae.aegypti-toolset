## Aus Vi Manhattan Plots

AusWiViWiFST<-read.delim("inputs/vcftoolsoutputs/AusVietWild27.weir.fst")
AusWiViColFST<-read.delim("inputs/vcftoolsoutputs/AusWildVietCol27.weir.fst")
ViWiViColFST<-read.delim("inputs/vcftoolsoutputs/VietWildVietCol27.weir.fst")

AusViet.df<-data.frame(AusWiViWiFST, AusWiViColFST$WEIR_AND_COCKERHAM_FST, ViWiViColFST$WEIR_AND_COCKERHAM_FST)
colnames(AusViet.df)<-c("CHROM", "POS", "AusWiViWiFST", "AusWiViColFST", "ViWiViColFST")

source("Scripts/MapFromDF_FSTplot.R")
MappedAusVi<-MapDF_FST(AusViet.df)
SNP<-1:length(MappedAusVi[,1])
MappedAusVi<-data.frame(SNP, MappedAusVi)

MappedAusBr<-MapDF_FST(AusBrz.df)
SNP<-1:length(MappedAusBr[,1])
MappedAusBr<-data.frame(SNP, MappedAusBr)



library(qqman)
library(scales)
library(ggExtra)
library(manhattanly)
library(ggsci)
colours<-pal_npg()(4)

jpeg("AusWiViColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusVi, col = alpha( c(colours[1], colours[4], colours[3], colours[1], colours[4], colours[3], colours[2]), 0.6), chr="JunChr", 
          bp="JunBP", p="AusWiViColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F')[ST]), ylim=c(-0.05, 1),
          chrlabs = c(1:3, "U1", "U2", "U3","Unassigned"))
dev.off()

makeManhattan<-function(df, colname, outliers) {
  library(ggsci)
  colours<-pal_npg()(5)
  filename<-paste0(colname, "FSTmanhattan.jpeg")
  source("Scripts/ModManhattan.R")
  if (missing(outliers)){
    jpeg(filename, 16, 4, units="in", res=300)
    modmanhattan(df, col = alpha( c(colours[1], colours[4], colours[3], colours[1], colours[4], colours[3], colours[2]), 0.6), 
              chr="JunChr", bp="JunBP", p=colname, snp="SNP", logp=FALSE, ylab=expression(italic('F')[ST]), ylim=c(-0.1, 1.1),
              chrlabs = c(1:3, "u1", "u2", "u3","Unassigned"), cex.label=0.1, 
              suggestiveline = mean(df[,colname], na.rm=T))
    dev.off()
  } else {
    outlierSNPs<-merge(outliers, df)
    View
    jpeg(filename, 16, 4, units="in", res=300)
    modmanhattan(df, col = alpha( c(colours[1], colours[4], colours[3], colours[1], colours[4], colours[3], colours[2]),0.6), 
              chr="JunChr", bp="JunBP", p=colname, snp="SNP", logp=FALSE, ylab=expression(italic('F')[ST]), ylim=c(-0.1, 1.1),
              highlight = outlierSNPs$SNP, suggestiveline = mean(df[,colname], na.rm=T))
    dev.off()
  }
}
makeManhattan(MappedAusVi, "AusWiViColFST", positions)
makeManhattan(MappedAusVi, "ViWiViColFST", positions)
makeManhattan(MappedAusVi, "AusWiViWiFST")

makeManhattan(MappedAusBr, "AusWiBrzWiFST")
makeManhattan(MappedAusBr, "AusWiKDColFST", KDColpositions)
makeManhattan(MappedAusBr, "KDColWColFST", KDColpositions)
makeManhattan(MappedAusBr, "BrzWiKDColFST", KDColpositions)

makeManhattan(MappedAusBr, "AusWiWColFST", WColpositions)
makeManhattan(MappedAusBr, "BrzWiWColFST", WColpositions)
makeManhattan(MappedAusBr, "AusWiBrzWiFST")

colours<-brewer.pal(4, "Dark2")
m<-manhattanr(x=MappedAusVi, chr="JunChr", bp="JunBP", p="AusWiViColFST", snp="SNP", logp=FALSE, 
              annotation1 = "JunSC", annotation2 = "JunBP")

plot<-manhattanly(m, genomewideline = F, 
            col = c(colours[1], colours[4], colours[3], colours[1], colours[4], colours[3], colours[2]),
            ylab='FST', labelChr = c(1:3, "U1", "U2", "U3","Unassigned"), ylim=c(-0.05, 1),
            suggestiveline = mean(MappedAusVi$AusWiViColFST, na.rm=T), 
            suggestiveline_color = "3084FE", title="",
            highlight = merge(positions, MappedAusVi)$SNP,
            highlight_color = "yellow") %>% 
  layout(yaxis = list(range = c(-0.05, 1)), width = 1500, height = 300)
plot

makeFSTmanhattanly<-function(df, FSTcolumn) {
  library(manhattanly)
  library(RColorBrewer)
  print(df$FSTColumn)
  colours<-brewer.pal(4, "Dark2")
  m<-manhattanr(x=df, chr="JunChr", bp="JunBP", p=FSTcolumn, snp="SNP", logp=FALSE, 
                annotation1 = "JunSC", annotation2 = "JunBP")
  plot<-manhattanly(m, genomewideline = F, 
                    col = c(colours[1], colours[4], colours[3], colours[1], colours[4], colours[3], colours[2]),
                    ylab='FST', labelChr = c(1:3, "U1", "U2", "U3","Unassigned"), ylim=c(-0.05, 1),
                    suggestiveline = mean(df$FSTcolumn, na.rm=T), 
                    suggestiveline_color = "3084FE", title="") %>% 
    layout(yaxis = list(range = c(-0.05, 1)), width = 1500, height = 300)
  plot
}
makeFSTmanhattanly(MappedAusVi, "AusWiViColFST")

MappedAusVi$"AusWiViColFST"  

jpeg("AusWiViWiFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusVi, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="AusWiViWiFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

jpeg("ViWiViColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusVi, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="ViWiViColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

## AusBraz Manhattan Plots

AusWiBrzWi<-read.delim("inputs/vcftoolsoutputs/AusBrazWild-20.weir.fst")
AusWiKDColFST<-read.delim("inputs/vcftoolsoutputs/AusColKD-20.weir.fst")
AusWiWColFST<-read.delim("inputs/vcftoolsoutputs/AusColW-20.weir.fst")
BrzWiKDColFST<-read.delim("inputs/vcftoolsoutputs/BrazWildColKD-20.weir.fst")
BrzWiWColFST<-read.delim("inputs/vcftoolsoutputs/BrazWildColW-20.weir.fst")
KDColWColFST<-read.delim("inputs/vcftoolsoutputs/ColKDColW-20.weir.fst")

AusBrz.df<-data.frame(AusWiBrzWi, AusWiKDColFST[,3], AusWiWColFST[,3], BrzWiWColFST[,3], BrzWiKDColFST[,3], KDColWColFST[,3])
colnames(AusBrz.df)<-c("CHROM", "POS", "AusWiBrzWiFST", "AusWiKDColFST", "AusWiWColFST","BrzWiWColFST", "BrzWiKDColFST", "KDColWColFST")

source("Scripts/MapFromDF.R")
MappedAusBr<-MapDF(AusBrz.df)
MappedAusBr<-MappedAusBr[-c(which(is.na(MappedAusBr$JunBP))),]
SNP<-1:length(MappedAusBr[,1])
MappedAusBr<-data.frame(SNP, MappedAusBr)

library(qqman)
library(scales)
library(ggExtra)
library(manhattanly)

jpeg("AusWiBrzWiFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="AusWiBrzWiFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

jpeg("AusWiColWFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="AusWiWColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

jpeg("AusWiKDColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="AusWiKDColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))

m<-manhattanr(x=MappedAusBr, chr="JunChr", bp="JunBP", p="KDColWColFST", snp="SNP", logp=FALSE, 
              annotation1 = "JunSC", annotation2 = "JunBP")
manhattanly(m, suggestiveline = F, genomewideline = F, col = c("#332288", "#88CCEE"),
            ylab='FST')
dev.off()

jpeg("BrzWiKDColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="BrzWiKDColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

jpeg("BrzWiWColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="BrzWiWColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

jpeg("WColKDColFSTManhattan.jpeg", 9,4, units="in", res=300)
manhattan(MappedAusBr, col = alpha(c("#332288", "#88CCEE"), 0.6), chr="JunChr", 
          bp="JunBP", p="KDColWColFST", snp="SNP", logp=FALSE, ylab=expression(italic('F'[ST])), 
          ylim=c(-0.05, 1))
dev.off()

