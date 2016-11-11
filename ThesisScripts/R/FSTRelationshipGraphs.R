# See Manhattan script to import FST data

## Aus Vi FST Relationship graph
library(ggExtra)
library(ggplot2)

p1 <- ggplot(AusViet.df, aes(ViWiViColFST, AusWiViColFST)) + 
  geom_point(colour="darkorchid4", alpha=0.5) + 
  theme_bw() + 
  geom_hline(aes(yintercept=mean(AusViet.df$AusWiViColFST, na.rm=T)), colour="Blue", linetype="dashed") +
  geom_vline(aes(xintercept=mean(AusViet.df$ViWiViColFST, na.rm=T)), colour="Green", linetype="dashed") +
  geom_abline(slope = 1, intercept = 0, linetype="dotted", colour="Black") +
  #geom_point(data=VietCbelow, aes(x=VietWi_VietC_FST, AusWi_VietC_FST), colour="magenta", alpha = 0.5) +
  scale_y_continuous(name=expression("AUwild-VNcol "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) +
  scale_x_continuous(name=expression("VNwild-VNcol "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
#p1
p1<-ggMarginal(p1, type = "histogram", colour = "darkorchid4", 
               binwidth = 0.02,
               xparams = list(colour = "white", fill="green", alpha=0.7),
               yparams = list(colour = "white", fill="blue", alpha=0.7))
p1
jpeg(file = "GvViWiColFSTRelGW.jpeg", width=6, height=6, units="in", res=300)
p1
dev.off()


## Aus Br FST Relationship Graphs

p2 <- ggplot(AusBrz.df, aes(BrzWiWColFST, AusWiWColFST)) + 
  geom_point(colour="darkorchid4", alpha=0.5) + 
  theme_bw() + 
  geom_hline(aes(yintercept=mean(AusBrz.df$AusWiWColFST, na.rm=T)), colour="Blue", linetype="dashed") +
  geom_vline(aes(xintercept=mean(AusBrz.df$BrzWiWColFST, na.rm=T)), colour="Green", linetype="dashed") +
  geom_abline(slope = 1, intercept = 0, linetype="dotted", colour="Black") +
  #geom_point(data=VietCbelow, aes(x=VietWi_VietC_FST, AusWi_VietC_FST), colour="magenta", alpha = 0.5) +
  scale_y_continuous(name=expression("AUwild-BRcol1 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) +
  scale_x_continuous(name=expression("BRwild-BRcol1 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p2
p2<-ggMarginal(p2, type = "histogram", colour = "darkorchid4", 
               binwidth = 0.02,
               xparams = list(colour = "white", fill="green", alpha=0.7),
               yparams = list(colour = "white", fill="blue", alpha=0.7))
p2
jpeg(file = "GvBrzWColFSTRelGW.jpeg", width=6, height=6, units="in", res=300)
p2
dev.off()

p3 <- ggplot(AusBrz.df, aes(BrzWiKDColFST, AusWiKDColFST)) + 
  geom_point(colour="darkorchid4", alpha=0.5) + 
  theme_bw() + 
  geom_hline(aes(yintercept=mean(AusBrz.df$AusWiKDColFST, na.rm=T)), colour="Blue", linetype="dashed") +
  geom_vline(aes(xintercept=mean(AusBrz.df$BrzWiKDColFST, na.rm=T)), colour="Green", linetype="dashed") +
  geom_abline(slope = 1, intercept = 0, linetype="dotted", colour="Black") +
  #geom_point(data=VietCbelow, aes(x=VietWi_VietC_FST, AusWi_VietC_FST), colour="magenta", alpha = 0.5) +
  scale_y_continuous(name=expression("AUwild-BRcol2 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) +
  scale_x_continuous(name=expression("BRwild-BRcol2 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p3<-ggMarginal(p3, type = "histogram", colour = "darkorchid4", 
               binwidth = 0.02,
               xparams = list(colour = "white", fill="green", alpha=0.7),
               yparams = list(colour = "white", fill="blue", alpha=0.7))
p3
jpeg(file = "GvBrzKDColFSTRelGW.jpeg", width=6, height=6, units="in", res=300)
p3
dev.off()

p4 <- ggplot(AusBrz.df, aes(BrzWiKDColFST, KDColWColFST)) + 
  geom_point(colour="darkorchid4", alpha=0.5) + 
  theme_bw() + 
  geom_hline(aes(yintercept=mean(AusBrz.df$KDColWColFST, na.rm=T)), colour="Blue", linetype="dashed") +
  geom_vline(aes(xintercept=mean(AusBrz.df$BrzWiKDColFST, na.rm=T)), colour="Green", linetype="dashed") +
  geom_abline(slope = 1, intercept = 0, linetype="dotted", colour="Black") +
  #geom_point(data=VietCbelow, aes(x=VietWi_VietC_FST, AusWi_VietC_FST), colour="magenta", alpha = 0.5) +
  scale_y_continuous(name=expression("BRcol1-BRcol2 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) +
  scale_x_continuous(name=expression("BRwild-BRcol2 "~italic(F[ST])), limits=c(-0.07,1), 
                     breaks=seq(from=-0.1, to=1, by=0.1), expand=c(0.01,0.01)) + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
p4<-ggMarginal(p4, type = "histogram", colour = "darkorchid4", 
               binwidth = 0.02,
               xparams = list(colour = "white", fill="green", alpha=0.7),
               yparams = list(colour = "white", fill="blue", alpha=0.7))
p4
jpeg(file = "BrzKD-WColFSTRelGW.jpeg", width=6, height=6, units="in", res=300)
p4
dev.off()

## Outliers below diagonal
# AusViet

belowAuVi<-AusViet.df[which(AusViet.df$ViWiViColFST > (0.05 + AusViet.df$AusWiViColFST)),c(1,2,5)]
ViWiViColQ95<-quantile(AusViet.df$ViWiViColFST, 0.95, na.rm=T)
AuViOutliers<-belowAuVi[belowAuVi$ViWiViColFST >= ViWiViColQ95,]

positions<-data.frame(AuViOutliers$CHROM, AuViOutliers$POS, stringsAsFactors = FALSE)
colnames(positions)<-c("CHROM", "POS")
write.table(positions, "ViWiColOutlierPositions.txt", sep="\t", quote = FALSE, row.names=FALSE, col.names = FALSE)

# AusBraz
#KDCol-BrzWild
belowAuBrzKD<-AusBrz.df[which(AusBrz.df$BrzWiKDColFST > (0.05 + AusBrz.df$AusWiKDColFST)), c(1,2,7)]
BrzWiKDColQ95<-quantile(AusBrz.df$BrzWiKDColFST, 0.95, na.rm=T)
AuBrzKDOutliers<-belowAuBrzKD[belowAuBrzKD$BrzWiKDColFST >= BrzWiKDColQ95,]
positions<-data.frame(AuBrzKDOutliers$CHROM, AuBrzKDOutliers$POS, stringsAsFactors = FALSE)
colnames(positions)<-c("CHROM", "POS")
write.table(positions, "KDBrWildOutlierPositions.txt", sep="\t", quote = FALSE, row.names=FALSE, col.names = FALSE)

#WCol-BrzWild
belowAuBrzW<-AusBrz.df[which(AusBrz.df$BrzWiWColFST > (0.05 + AusBrz.df$AusWiWColFST)), c(1,2,6)]
BrzWiWColQ95<-quantile(AusBrz.df$BrzWiWColFST, 0.95, na.rm=T)
AuBrzWOutliers<-belowAuBrzW[belowAuBrzW$BrzWiWColFST >= BrzWiWColQ95,]
WColpositions<-data.frame(AuBrzWOutliers$CHROM, AuBrzWOutliers$POS, stringsAsFactors = FALSE)
colnames(WColpositions)<-c("CHROM", "POS")
write.table(positions, "WColBrWildOutlierPositions.txt", sep = "\t", quote = FALSE, row.names=FALSE, col.names = FALSE)

#KDCol-WCol
belowWKD<-AusBrz.df[which(AusBrz.df$BrzWiKDColFST > (0.05 + AusBrz.df$KDColWColFST)), c(1,2,7)]
BrzWColBrzKDColOutliers<-belowWKD[belowWKD$BrzWiKDColFST >= BrzWiKDColQ95,]
KDColpositions<-data.frame(BrzWColBrzKDColOutliers$CHROM, BrzWColBrzKDColOutliers$POS, stringsAsFactors = FALSE)
colnames(KDColpositions)<-c("CHROM", "POS")
write.table(positions, "KDColWColOutlierPositions.txt", sep="\t", quote = FALSE, row.names=FALSE, col.names = FALSE)
