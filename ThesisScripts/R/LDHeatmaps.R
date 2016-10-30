###################
## AU-VN Dataset ##
###################

# Load R2 matrices created using vcf2Rsquared()
AusWild<-as.matrix(read.table("inputs/subsampledvcfs/AusWildLD-27ALLRsquMatrix.txt", check.names=FALSE))
VietWild<-as.matrix(read.table("inputs/subsampledvcfs/VietWildLD-27ALLRsquMatrix.txt", check.names=FALSE))
VietCol<-as.matrix(read.table("inputs/subsampledvcfs/VietColLD-27ALLRsquMatrix.txt", check.names = FALSE))

library(pheatmap)
library(heatmaply)
library(plyr)
library(gplots)
library(LDheatmap)

# Transpose to make a symmetrical matrix
AusWild[lower.tri(AusWild)]<-t(AusWild)[lower.tri(AusWild)]
VietWild[lower.tri(VietWild)]<-t(VietWild)[lower.tri(VietWild)]
VietCol[lower.tri(VietCol)]<-t(VietCol)[lower.tri(VietCol)]

# Map markers to the Juneja assembly
source("scripts/MappingFunctions.R")
MappedLoci<-MapVCF("inputs/subsampledvcfs/VietWildLD-27.recode.vcf")

# Order markers by Chromosome, then basepair position
MappedLoci<-MappedLoci[order(MappedLoci$JunChr, MappedLoci$JunBP),]
# Remove markers mapped to chromosome without position
MappedLoci<-MappedLoci[complete.cases(MappedLoci),]

# Store IDs relating to each chromosome
Chrom1IDs<-MappedLoci[which(MappedLoci$JunChr=='1'),3]
Chrom2IDs<-MappedLoci[which(MappedLoci[,5]=='2'),3]
Chrom3IDs<-MappedLoci[which(MappedLoci[,5]=='3'),3]

# Create separate matrices for each chromosome
AusWildChrom1Mat<-AusWild[as.character(Chrom1IDs), as.character(Chrom1IDs)]
VietWildChrom1Mat<-VietWild[as.character(Chrom1IDs), as.character(Chrom1IDs)]
VietColChrom1Mat<-VietCol[as.character(Chrom1IDs), as.character(Chrom1IDs)]

# Save matrices to file (optional)
write.table(AusWildChrom1Mat, "app/AusWildChrom1Gv.txt")
write.table(VietWildChrom1Mat, "app/VietWildChrom1.txt")
write.table(VietColChrom1Mat, "app/VietColChrom1.txt")

### Chromosome 1 Heatmaps

# vector of superconts mapped to chrom 1
supercontsChrom1<-MappedLoci$JunSC[which(MappedLoci$JunChr==1)]

# Triangular heatmaps as used in thesis
jpeg("haploheatmapSubSamp27Chrom1VietCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietColChrom1Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapSubSamp27Chrom1AusWildViet.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(AusWildChrom1Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp27Chrom1VietWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietWildChrom1Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))
dev.off()

# Square heatmaps
pheatmap(VietColChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         #filename = "VietColChrom1.tiff",
         height = 6, width = 6)

pheatmap(AusWildChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         #filename = "AusWildViChrom1.tiff",
         height = 6, width = 6)

pheatmap(VietWildChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 12, fontsize_row=5, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "VietWildChrom1.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom1)

# interactive heatmaps
heatmaply(VietColChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1),
          title="VietCol")
heatmaply(VietWildChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))
heatmaply(AusWildChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))

### Chromosome 2 Heatmaps

# Create chromosome 2 matrices
AusWildChrom2Mat<-AusWild[as.character(Chrom2IDs), as.character(Chrom2IDs)]
VietWildChrom2Mat<-VietWild[as.character(Chrom2IDs), as.character(Chrom2IDs)]
VietColChrom2Mat<-VietCol[as.character(Chrom2IDs), as.character(Chrom2IDs)]
write.table(AusWildChrom2Mat, "app/AusWildChrom2Gv.txt")
write.table(VietWildChrom2Mat, "app/VietWildChrom2.txt")
write.table(VietColChrom2Mat, "app/VietColChrom2.txt")

supercontsChrom2<-MappedLoci$JunSC[which(MappedLoci$JunChr==2)]

# Triangular heatmaps as used in thesis

jpeg("haploheatmapsusbsamp27Chrom2VietCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietColChrom2Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp27Chrom2AusWildViet.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(AusWildChrom2Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp27Chrom2VietWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietWildChrom2Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))
dev.off()


# Square heatmaps
pheatmap(VietColChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "VietColChrom2.tiff",
         height = 6, width = 6)

pheatmap(AusWildChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "AusWildViChrom2.tiff",
         height = 6, width = 6)

pheatmap(VietWildChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 12, fontsize_row=5, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "VietWildChrom2.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom2,
         labels_col = supercontsChrom2)

# Interactive heatmaps

heatmaply(VietColChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(VietWildChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))
heatmaply(AusWildChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1), hide_colorbar = T)

## Chromosome 3 Heatmaps

AusWildChrom3Mat<-AusWild[as.character(Chrom3IDs), as.character(Chrom3IDs)]
VietWildChrom3Mat<-VietWild[as.character(Chrom3IDs), as.character(Chrom3IDs)]
VietColChrom3Mat<-VietCol[as.character(Chrom3IDs), as.character(Chrom3IDs)]

write.table(AusWildChrom3Mat, "app/AusWildChrom3Gv.txt")
write.table(VietWildChrom3Mat, "app/VietWildChrom3.txt")
write.table(VietColChrom3Mat, "app/VietColChrom3.txt")

supercontsChrom3<-MappedLoci$JunSC[which(MappedLoci$JunChr==3)]

# Triangular heatmaps as used in thesis

jpeg("haploheatmapsubsamp27Chrom3VietCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietColChrom3Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp27Chrom3AusWildViet.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(AusWildChrom3Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp27Chrom3VietWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(VietWildChrom3Mat, genetic.distances = MappedLoci$JunBP[which(MappedLoci$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))
dev.off()

# Square heatmaps

pheatmap(VietColChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "VietColChrom3.tiff",
         height = 6, width = 6)



pheatmap(AusWildChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "AusWildViChrom3.tiff",
         height = 6, width = 6)

pheatmap(VietWildChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 12, fontsize_row=5, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 1 wMelRio",
         filename = "VietWildChrom3.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom3,
         labels_col = supercontsChrom3)

# Interactive heatmaps

heatmaply(VietColChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(VietWildChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))


heatmaply(AusWildChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))

####################
### AU-BR Dataset ##
####################

# Load R2 matrices created using vcf2Rsquared()
AusWildBr<-as.matrix(read.table("inputs/subsampledvcfs/AusWiBrz20ALLRsquMatrix.txt", check.names=FALSE))
BrWild<-as.matrix(read.table("inputs/subsampledvcfs/BrzWild20ALLRsquMatrix.txt", check.names=FALSE))
BrKDCol<-as.matrix(read.table("inputs/subsampledvcfs/BrzKDCol20ALLRsquMatrix.txt", check.names = FALSE))
BrWCol<-as.matrix(read.table("inputs/subsampledvcfs/BrzColW20ALLRsquMatrix.txt", check.names = FALSE))

# Transpose to a symmetrical matrix
AusWildBr[lower.tri(AusWildBr)]<-t(AusWildBr)[lower.tri(AusWildBr)]
BrWild[lower.tri(BrWild)]<-t(BrWild)[lower.tri(BrWild)]
BrKDCol[lower.tri(BrKDCol)]<-t(BrKDCol)[lower.tri(BrKDCol)]
BrWCol[lower.tri(BrWCol)]<-t(BrWCol)[lower.tri(BrWCol)]

# Map loci to Juneja positions
MappedLociBr<-MapVCF("inputs/subsampledvcfs/BrazWildLD-20.recode.vcf")
MappedLociBr<-MappedLociBr[order(MappedLociBr$JunChr, MappedLociBr$JunBP),]
MappedLociBr<-MappedLociBr[complete.cases(MappedLociBr),]

# Store IDs for each chromosome
Chrom1IDs<-as.character(MappedLociBr[which(MappedLociBr$JunChr==1),3])
Chrom2IDs<-as.character(MappedLociBr[which(MappedLociBr[,5]=='2'),3])
Chrom3IDs<-as.character(MappedLociBr[which(MappedLociBr[,5]=='3'),3])

# Chromosome 1
AusWildBrChrom1Mat<-AusWildBr[as.character(Chrom1IDs), as.character(Chrom1IDs)]
BrWildChrom1Mat<-BrWild[as.character(Chrom1IDs), as.character(Chrom1IDs)]
BrKDColChrom1Mat<-BrKDCol[as.character(Chrom1IDs), as.character(Chrom1IDs)]
BrWColChrom1Mat<-BrWCol[as.character(Chrom1IDs), as.character(Chrom1IDs)]

# Write matrices to file (optional)
write.table(AusWildBrChrom1Mat, "app/AusWildBrChrom1.txt")
write.table(BrWildChrom1Mat, "app/BrWildChrom1.txt")
write.table(BrKDColChrom1Mat, "app/BrKDColChrom1.txt")
write.table(BrWColChrom1Mat, "app/BrWColChrom1.txt")

supercontsChrom1<-MappedLociBr$JunSC[which(MappedLociBr$JunChr==1)]

# Triangular heatmaps

jpeg("haploheatmapsubsamp20Chrom1BrKDCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrKDColChrom1Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsampChrom1BrWCol.jpeg", height=7, width = 7, units="in", res=300)
chrom1WCol<-LDheatmap(BrWColChrom1Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==1)],
                      color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
                      LDmeasure = "r", distances = "physical", title = NULL,
                      add.key=F, SNP.name=c("8431", "82246"))
LDheatmap.highlight(chrom1WCol, 16,67, col = "red", lwd = 1)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
grid.edit(gPath("ldheatmap", "geneMap", "SNPnames"), gp=gpar(cex=0.8, col="red"))
grid.edit(gPath("ldheatmap", "geneMap", "symbols"), gp=gpar(cex=1.5, col="red"))
dev.off()

jpeg("haploheatmapsubsamp20Chrom1AusWildBrBr.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(AusWildBrChrom1Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL,
          add.key=F)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp20Chrom1BrWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrWildChrom1Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==1)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))
dev.off()

# Square heatmaps

pheatmap(BrKDColChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         #filename = "wMelRioChrom1.tiff",
         height = 6, width = 6)

pheatmap(BrWColChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 1 wMelRio",
         #filename = "wMelRioChrom1.tiff",
         height = 6, width = 6)

pheatmap(AusWildBrChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F, main="Chromosome 1 wMelRio",
         #filename = "AusWildBrBRChrom1.tiff",
         height = 6, width = 6)

pheatmap(BrWildChrom1Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 14, fontsize_row=5, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 1 wMelRio",
         #filename = "BrWildChrom1.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom1)


# Interactive heatmaps

heatmaply(BrKDColChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWColChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWildChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))
heatmaply(AusWildBrChrom1Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))

# High LD section close up heatmap

# IDs of Loci in high LD section
highLDIDs<-MappedLociBr$ID[16:67]
# Create matrix with these loci
highLDMat<-BrWColChrom1Mat[c(16:67), c(16:67)]
# Create mapped table with only highLD loci
HighLDMappedLociBr<-MappedLociBr[16:67,]
head(HighLDMappedLociBr)

# Load sex-linked markers from Fontaine et al. 2016 (http://biorxiv.org/content/early/2016/06/21/060061.figures-only)
sexlinked<-read.csv("inputs/sexlinkedMarkers.csv")
# removed sexlinked marker outside region
sexlinked<-sexlinked[c(-1,-2),]
colnames(sexlinked)<-c("CHROM", "POS", "chrom_pos")
sexlinked$CHROM<-gsub("sc", "supercont", sexlinked$CHROM)
sexlinked<-MapDF(sexlinked)
sexlinked<-sexlinked[complete.cases(sexlinked),]

# Create dataframe for mapping sexlinked markers
sexlinkedMap<-data.frame(unique(sexlinked$JunSC))
colnames(sexlinkedMap)<-"JunSC"
right=left=1:length(sexlinkedMap[,1])

# Calculate left to right span of markers on a supercontig
for (i in 1:length(sexlinkedMap[,1])){
  positions<-sexlinked$JunBP[which(sexlinked$JunSC == sexlinkedMap$JunSC[i])]
  left[i]<-min(positions)
  right[i]<-max(positions)
}
sexlinkedMap<-data.frame(sexlinkedMap, left, right)

# Do the same for BRcol1 markers
supercontMap<-data.frame(unique(HighLDMappedLociBr$JunSC))
colnames(supercontMap)<-"JunSC"
left=c()
right=c()
for (i in 1:length(supercontMap[,1])){
  positions<-HighLDMappedLociBr$JunBP[which(HighLDMappedLociBr$JunSC == supercontMap$JunSC[i])]
  left<-c(left, min(positions))
  right<-c(right,max(positions))
}
supercontMap<-data.frame(supercontMap, left, right)

## create SNP locations plot
supercontsmap<-ggplot(data=supercontMap,aes(x=left, y=1, label=JunSC)) +
  geom_segment(aes(x=left, y=1, xend=right, yend=1), color="gray") +
  geom_segment(data=sexlinkedMap, aes(x=left, y=0.99, xend=right, yend=0.99), color="steelblue4", alpha=0.7) +
  geom_point(shape="|", size=3.5) +
  geom_point(data=HighLDMappedLociBr, aes(x=JunBP, y=1), shape="|", size=3.5) +
  geom_point(data=sexlinked, aes(x=JunBP, y=0.99),shape="|", size=3.5, colour="steelblue4") +
  geom_text(data=supercontMap, aes(x=left, y=1, label=JunSC), angle=90, nudge_y = 0.01, hjust=0, size=2.5, check_overlap = F) +
  geom_text(data=sexlinkedMap, aes(x=left, y=1, label=JunSC), angle=90, nudge_x=100000, nudge_y = -0.03,hjust=0, size=2.5, color="steelblue4") +
  scale_x_continuous(expand = c(0,1000000)) +
  scale_y_continuous(limits = c(0.97, 1.03), breaks=c(0.99,1,1.03), expand = c(0,0)) +
  theme_classic() +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        legend.position="none",
        axis.ticks=element_blank())
supercontsmap
supercontsmap<-ggplotGrob(supercontsmap)

## create heatmap with SNP locations, 
jpeg("supercontmap.jpeg", height=7, width = 7, units="in", res=300)
highLDheatmap<-LDheatmap(highLDMat, genetic.distances = MappedLociBr$JunBP[16:67],
                         color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
                         LDmeasure = "r", distances = "physical", title = NULL,
                         add.key=T, flip=T)
LDheatmap.addGrob(highLDheatmap, supercontsmap, height=0.4)
#pushViewport(viewport(x=1, y=1, width=.99, height=0.2))
dev.off()

# Final figure adjusted in Illustrator to tidy up labels and adjust spacing to match heatmap
postscript("HighLDheatmap.eps")
highLDheatmap<-LDheatmap(highLDMat, genetic.distances = MappedLociBr$JunBP[16:67],
                         color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
                         LDmeasure = "r", distances = "physical", title = NULL,
                         add.key=T, flip=T)

LDheatmap.addGrob(highLDheatmap, supercontsmap, height=0.4)
dev.off()


# Chromosome 2
# Load data
AusWildBrChrom2Mat<-AusWildBr[as.character(Chrom2IDs), as.character(Chrom2IDs)]
BrWildChrom2Mat<-BrWild[as.character(Chrom2IDs), as.character(Chrom2IDs)]
BrKDColChrom2Mat<-BrKDCol[as.character(Chrom2IDs), as.character(Chrom2IDs)]
BrWColChrom2Mat<-BrWCol[as.character(Chrom2IDs), as.character(Chrom2IDs)]
write.table(AusWildBrChrom2Mat, "app/AusWildBrChrom2.txt")
write.table(BrWildChrom2Mat, "app/BrWildChrom2.txt")
write.table(BrKDColChrom2Mat, "app/BrKDColChrom2.txt")
write.table(BrWColChrom2Mat, "app/BrWColChrom2.txt")

supercontsChrom2<-MappedLociBr$JunSC[which(MappedLociBr$JunChr==2)]

# Triangle heatmaps

jpeg("haploheatmapsubsamp20Chrom2BrKDCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrKDColChrom2Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp20Chrom2BrWCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrWColChrom2Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

jpeg("haploheatmapsubsamp20Chrom2AusWildBr.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(AusWildBrChrom2Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()


jpeg("haploheatmapsubsamp20Chrom2BrWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrWildChrom2Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==2)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=T)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))
dev.off()

# Symmetrical heatmaps

pheatmap(BrKDColChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 2 wMelRio",
         filename = "wMelRioChrom2.tiff",
         height = 6, width = 6)

pheatmap(BrWColChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 2 wMelBr",
         filename = "wMelBrChrom2.tiff",
         height = 6, width = 6)

pheatmap(AusWildBrChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 2 wMelRio",
         filename = "AusWildBrBRChrom2.tiff",
         height = 6, width = 6)


pheatmap(BrWildChrom2Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 14, fontsize_row=4, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 2 wMelRio",
         filename = "BrWildChrom2.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom2)

# Interactive heatmaps

heatmaply(BrKDColChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWColChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWildChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))
heatmaply(AusWildBrChrom2Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1), hide_colorbar = T)

# Chromosome 3

AusWildBrChrom3Mat<-AusWildBr[as.character(Chrom3IDs), as.character(Chrom3IDs)]
BrWildChrom3Mat<-BrWild[as.character(Chrom3IDs), as.character(Chrom3IDs)]
BrKDColChrom3Mat<-BrKDCol[as.character(Chrom3IDs), as.character(Chrom3IDs)]
BrWColChrom3Mat<-BrWCol[as.character(Chrom3IDs), as.character(Chrom3IDs)]

write.table(AusWildBrChrom3Mat, "app/AusWildBrChrom3.txt")
write.table(BrWildChrom3Mat, "app/BrWildChrom3.txt")
write.table(BrKDColChrom3Mat, "app/BrKDColChrom3.txt")
write.table(BrWColChrom3Mat, "app/BrWColChrom3.txt")

supercontsChrom3<-MappedLociBr$JunSC[which(MappedLociBr$JunChr==3)]
pheatmap(BrKDColChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend=F, drop_levels = F,
         #main="Chromosome 2 wMelRio",
         filename = "wMelRioChrom3.tiff",
         height = 6, width = 6)

jpeg("haploheatmapsubsamp20Chrom3BKDCol.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrKDColChrom3Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))

dev.off()


pheatmap(BrWColChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 2 wMelBr",
         filename = "wMelBrChrom3.tiff",
         height = 6, width = 6)

jpeg("haploheatmapsubsamp20Chrom3BWCol.jpeg", height=7, width = 7, units="in", res=300)
Chrom3Wheatmap<-LDheatmap(BrWColChrom3Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
LDheatmap.highlight(Chrom3Wheatmap, 111,125, col = "red", lwd = 1)
LDheatmap.highlight(Chrom3Wheatmap, 35,42, col = "blue", lwd = 1)
dev.off()

chrom3mapped<-MappedLociBr[which(MappedLociBr$JunChr==3),]
chrom3mapped[111:125,]
chrom3mapped[35:42,]

rownames(BrWColChrom3Mat[111,125])

pheatmap(AusWildBrChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 8, show_rownames = F, show_colnames = F,
         legend = F, drop_levels = F,
         #main="Chromosome 2 wMelRio",
         #filename = "AusWildBrBRChrom3.tiff",
         height = 6, width = 6)

jpeg("haploheatmapsubsamp20Chrom3AusWildBr.jpeg", height=7, width = 7, units="in", res=300)
wColChrom3<-LDheatmap(AusWildBrChrom3Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=F, flip = T)
LDheatmap.highlight(wColChrom3, i = )
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(col = "white"))
dev.off()

pheatmap(BrWildChrom3Mat, color=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
         cluster_rows = F, cluster_cols = F, border_color = NA,
         na.value="white", fontsize = 14, fontsize_row=4, show_colnames = F,
         legend_breaks = seq(0,1,0.25), drop_levels = F,
         #main="Chromosome 2 wMelRio",
         filename = "BrWildChrom3.tiff",
         height = 6, width = 7,
         labels_row=supercontsChrom3)

jpeg("haploheatmapsubsamp20Chrom3BrWild.jpeg", height=7, width = 7, units="in", res=300)
LDheatmap(BrWildChrom3Mat, genetic.distances = MappedLociBr$JunBP[which(MappedLociBr$JunChr==3)],
          color = rev(plasma(n = 256, alpha = 1, begin = 0, end = 1)),
          LDmeasure = "r", distances = "physical", title = NULL, add.key=T)#, flip = T)
grid.edit(gPath("ldheatmap", "geneMap", "title"), gp = gpar(cex = 1.2))
grid.edit(gPath("ldheatmap", "Key", "title"), gp=gpar(cex=1.2))
grid.edit(gPath("ldheatmap", "Key", "labels"), gp=gpar(cex=1.1))

dev.off()

heatmaply(BrKDColChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWColChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="4px",
          yaxis_font_size="4px",
          limits = c(0,1))
heatmaply(BrWildChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))
heatmaply(AusWildBrChrom3Mat, Rowv=FALSE, Colv=F, colors=plasma(n = 256, alpha = 1, begin = 0, end = 1), 
          xaxis_height=300, yaxis_height=300,
          xaxis_font_size="6px",
          yaxis_font_size="6px",
          limits = c(0,1))