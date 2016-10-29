# Load required libraries
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

# Load required functions
source("scripts/Functions.R")
source("scripts/MappingFunctions.R") # requires JunejaGenetic assembly csv.

##################
## Vietnam Data ##
##################

AusWild<-as.matrix(read.table("inputs/subsampledvcfs/AusWildLD-27ALLRsquMatrix.txt", check.names=FALSE))
VietWild<-as.matrix(read.table("inputs/subsampledvcfs/VietWildLD-27ALLRsquMatrix.txt", check.names=FALSE))
VietCol<-as.matrix(read.table("inputs/subsampledvcfs/VietColLD-27ALLRsquMatrix.txt", check.names = FALSE))

AusWild<-melt(AusWild, na.rm = TRUE)
VietWild<-melt(VietWild, na.rm=TRUE)
VietCol<-melt(VietCol, na.rm=TRUE)

colnames(AusWild)<-c("Loc1","Loc2","LD")
colnames(VietWild)<-c("Loc1", "Loc2", "LD")
colnames(VietCol)<-c("Loc1", "Loc2", "LD")

MappedLoci<-MapVCFLoci("inputs/subsampledvcfs/VietColLD-27.recode.vcf")
Chrom1IDs<-MappedLoci[which(MappedLoci[,5]=='1'),3]
Chrom2IDs<-MappedLoci[which(MappedLoci[,5]=='2'),3]
Chrom3IDs<-MappedLoci[which(MappedLoci[,5]=='3'),3]

VietCol$Chrom1<-VietWild$Chrom1<-AusWild$Chrom1<-sapply(AusWild$Loc1, assignsChroms)
VietCol$Chrom2<-VietWild$Chrom2<-AusWild$Chrom2<-sapply(AusWild$Loc2, assignsChroms)
VietCol<-VietCol[complete.cases(VietCol),]
VietWild<-VietWild[complete.cases(VietWild),]
AusWild<-AusWild[complete.cases(AusWild),]

allchromsVi<-data.frame(AusWild, VietCol$LD, VietWild$LD)
colnames(allchromsVi)<-c("Loc1", "Loc2", "AusLD", "Chrom1", "Chrom2", "VietCol", "VietWild")
allchromsViSynt<-allchromsVi[which(allchromsVi$Chrom1 == allchromsVi$Chrom2),]
allChromsViSynt<-CalcDist(allchromsViSynt,MappedLoci)
allChromsViSynt5MB<-allChromsViSynt[allChromsViSynt$distances<5000000,]

AusViallChromsfit<-LDdecayFit(allChromsViSynt5MB$distances, allChromsViSynt5MB$AusLD, 54)
AusViallChromsLDPlot<-LDdecayPredict(AusViallChromsfit, allChromsViSynt5MB$distances, 54, "AUwild")

VietColallChromsfit<-LDdecayFit(allChromsViSynt5MB$distances, allChromsViSynt5MB$VietCol, 54)
VietColallChromsLDPlot<-LDdecayPredict(VietColallChromsfit, allChromsViSynt5MB$distances, 54, "VNcol")

VietWildallChromsfit<-LDdecayFit(allChromsViSynt5MB$distances, allChromsViSynt5MB$VietWild, 54)
VietWildallChromsLDPlot<-LDdecayPredict(VietWildallChromsfit, allChromsViSynt5MB$distances, 54, "VNwild")

allChromsLDPlot<-rbind(AusViallChromsLDPlot, VietColallChromsLDPlot, VietWildallChromsLDPlot)

p4all<-
  ggplot(data=allChromsViSynt5MB, aes(x=distances/1000000, y=AusLD)) +
  #ggtitle("LD r^2 versus Physical BrWColonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=AusViallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#332288") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("AUwild LD (" ~ r^{2} ~ ")"), limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p4all

p3all<-
  ggplot(data=allChromsViSynt5MB, aes(x=distances/1000000, y=VietCol)) +
  #ggtitle("LD r^2 versus Physical BrWColonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=VietColallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#44AA99") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("VNcol LD (" ~ r^{2} ~ ")"), limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p3all

p1all<-
  ggplot(data=allChromsViSynt5MB, aes(x=distances/1000000, y=VietWild)) +
  #ggtitle("LD r^2 versus Physical BrWildonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=VietWildallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#882255") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("VNwild LD (" ~ r^{2} ~ ")"), limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p1all

jpeg(file="GvViSubSamp27GWLDvsPhysical.jpeg", 10.5,3,units="in", res=300)
grid.arrange(p4all,p3all,p1all, nrow=1)
dev.off()

ALLLDDecay<-ggplot(data=allChromsLDPlot, aes(x=distances/1000000, y=fpoints,
                                             ymin=upperpoints, ymax=lowerpoints,
                                             group=pop, colour = pop, fill=pop)) +
  geom_line() +
  geom_ribbon(colour=NA, alpha=0.2) +
  scale_fill_manual(values=c("#332288", "#44AA99", "#882255")) +
  scale_color_manual(values=c("#332288", "#44AA99","#882255" )) +
  scale_x_continuous(name="Intermarker distance (mb)", limits = c(0,5)) +
  scale_y_continuous(name=expression("LD (" ~ r^{2} ~ ")"), limits = c(0,0.5)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(size=8),
        legend.title = element_blank(),
        legend.position = c(0.8,0.75),
        legend.key.size = unit(1,"cm"))
ALLLDDecay

jpeg(file="GvViSubSamp27GenomeWideLDdecaycurves.jpeg", 9,5, units="in", res=300)
ALLLDDecay
dev.off()

# Background LD

mean(AusWild$LD[which(AusWild$Chrom1 != AusWild$Chrom2)])
mean(VietWild$LD[which(VietWild$Chrom1 != VietWild$Chrom2)])
mean(VietCol$LD[which(VietCol$Chrom1 != VietCol$Chrom2)])
sd(AusWild$LD[which(AusWild$Chrom1 != AusWild$Chrom2)])
sd(VietWild$LD[which(VietWild$Chrom1 != VietWild$Chrom2)])
sd(VietCol$LD[which(VietCol$Chrom1 != VietCol$Chrom2)])

# Syntenic LD

mean(allChromsViSynt5MB$AusLD)
mean(allChromsViSynt5MB$VietCol)
mean(allChromsViSynt5MB$VietWild)

sd(allChromsViSynt5MB$AusLD)
sd(allChromsViSynt5MB$VietCol)
sd(allChromsViSynt5MB$VietWild)

####################
## Brazil dataset ##
####################

AusWildBr<-as.matrix(read.table("inputs/subsampledvcfs/AusWiBrz20ALLRsquMatrix.txt", check.names=FALSE))
BrWild<-as.matrix(read.table("inputs/subsampledvcfs/BrzWild20ALLRsquMatrix.txt", check.names=FALSE))
BrKDCol<-as.matrix(read.table("inputs/subsampledvcfs/BrzKDCol20ALLRsquMatrix.txt", check.names = FALSE))
BrWCol<-as.matrix(read.table("inputs/subsampledvcfs/BrzColW20ALLRsquMatrix.txt", check.names = FALSE))

AusWildBr<-melt(AusWildBr, na.rm=T)
BrWild<-melt(BrWild, na.rm=T)
BrKDCol<-melt(BrKDCol, na.rm=T)
BrWCol<-melt(BrWCol, na.rm=T)

colnames(BrWCol)<-colnames(BrKDCol)<-colnames(BrWild)<-colnames(AusWildBr)<-c("Loc1","Loc2","LD")

MappedLoci<-MapVCFLoci("inputs/subsampledvcfs/BrazWildLD-20.recode.vcf")
Chrom1IDs<-MappedLoci[which(MappedLoci[,5]=='1'),3]
Chrom2IDs<-MappedLoci[which(MappedLoci[,5]=='2'),3]
Chrom3IDs<-MappedLoci[which(MappedLoci[,5]=='3'),3]

BrKDCol$Chrom1<-BrWCol$Chrom1<-BrWild$Chrom1<-AusWildBr$Chrom1<-sapply(AusWildBr$Loc1, assignsChroms)
BrKDCol$Chrom2<-BrWCol$Chrom2<-BrWild$Chrom2<-AusWildBr$Chrom2<-sapply(AusWildBr$Loc2, assignsChroms)

BrKDCol<-BrKDCol[complete.cases(BrKDCol),]
BrWCol<-BrWCol[complete.cases(BrWCol),]
BrWild<-BrWild[complete.cases(BrWild),]
AusWildBr<-AusWildBr[complete.cases(AusWildBr),]

allchroms<-data.frame(AusWildBr, BrKDCol$LD, BrWild$LD, BrWCol$LD)
colnames(allchroms)<-c("Loc1", "Loc2", "AusLD", "Chrom1", "Chrom2", "BrKDCol", "BrWild", "BrWCol")

allchroms<-allchroms[which(allchroms$Chrom1 == allchroms$Chrom2),]
allChromsdist<-CalcDist(allchroms,MappedLoci)
allChromsdist<-allChromsdist[allChromsdist$distances<5000000,]

AusallChromsfit<-LDdecayFit(allChromsdist$distances, allChromsdist$AusLD, 40)
AusallChromsLDPlot<-LDdecayPredict(AusallChromsfit, allChromsdist$distances, 40, "AUwild")

KDColallChromsfit<-LDdecayFit(allChromsdist$distances, allChromsdist$BrKDCol, 40)
KDColallChromsLDPlot<-LDdecayPredict(KDColallChromsfit, allChromsdist$distances, 40, "BRcol2")

BrWColAllChromsFit<-LDdecayFit(allChromsdist$distances, allChromsdist$BrWCol, 40)
BrWColallChromsLDPlot<-LDdecayPredict(BrWColAllChromsFit, allChromsdist$distances, 40, "BRcol1")

BrWildAllChromsFit<-LDdecayFit(allChromsdist$distances, allChromsdist$BrWild, 40)
BrWildallChromsLDPlot<-LDdecayPredict(BrWildAllChromsFit, allChromsdist$distances, 40, "BRwild")

allchromsLDPlotBr<-rbind(AusallChromsLDPlot, KDColallChromsLDPlot, BrWColallChromsLDPlot, BrWildallChromsLDPlot)

p4all<-
  ggplot(data=allChromsdist, aes(x=distances/1000000, y=AusLD)) +
  #ggtitle("LD r^2 versus Physical BrWColonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=AusallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#332288") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("AUwild LD (" ~ r^{2} ~ ")"),limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p4all

p3all<-
  ggplot(data=allChromsdist, aes(x=distances/1000000, y=BrKDCol)) +
  #ggtitle("LD r^2 versus Physical BrWColonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=KDColallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#40AA99") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("BRcol2 LD (" ~ r^{2} ~ ")"),limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p3all

p2all<-
  ggplot(data=allChromsdist, aes(x=distances/1000000, y=BrWCol)) +
  #ggtitle("LD r^2 versus Physical BrWColonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=BrWColallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#117733") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("BRcol1 LD (" ~ r^{2} ~ ")"), limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p2all

p1all<-
  ggplot(data=allChromsdist, aes(x=distances/1000000, y=BrWild)) +
  #ggtitle("LD r^2 versus Physical BrWildonyChrom3") +
  geom_point(col="grey", alpha=4/10) +
  geom_line(data=BrWildallChromsLDPlot, aes(x=distances/1000000, y=fpoints), col="#882255") +
  #geom_line(data=BrWColChrom3MeanDecay, aes(plottingbins/1000000, means), col="#117733", linetype="dashed") +
  scale_y_continuous(expression("BRwild LD (" ~ r^{2} ~ ")"), limits=c(0,1)) +
  scale_x_continuous("Intermarker distance (mb)",breaks=c(seq(from=0, to = 5, by = 1)), limits=c(0,5), expand = c(0.01,0.01)) +
  theme_bw(base_size = 16)
p1all

jpeg(file="GvBrSubSamp20GenomeWideLDvsPhysical.jpeg", 14,3,units="in", res=300)
grid.arrange(p4all,p2all,p3all,p1all, nrow=1)
dev.off()

ALLLDDecay<-ggplot(data=allchromsLDPlotBr, aes(x=distances/1000000, y=fpoints,
                                             ymin=upperpoints, ymax=lowerpoints,
                                             group=pop, colour = pop, fill=pop)) +
  geom_line() +
  geom_ribbon(colour=NA, alpha=0.2) +
  scale_fill_manual(values=c("#332288", "#44AA99", "#117733","#882255"), 
                    breaks=c("AUwild", "BRcol1", "BRcol2", "BRwild")) +
  scale_color_manual(values=c("#332288", "#44AA99","#117733", "#882255"),
                     breaks=c("AUwild", "BRcol1", "BRcol2", "BRwild")) +
  scale_x_continuous(name="Intermarker distance (mb)", limits = c(0,5)) +
  scale_y_continuous(name=expression("LD (" ~ r^{2} ~ ")"), limits = c(0,0.5)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(size=8),
        legend.title = element_blank(),
        legend.position = c(0.8,0.72),
        legend.key.size = unit(1,"cm"))
ALLLDDecay

jpeg(file="GvBrSubSamp20GenomeWideLDdecaycurves.jpeg", 9,5, units="in", res=300)
ALLLDDecay
dev.off()
