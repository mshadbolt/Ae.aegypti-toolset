# Creates SNP density graph from VCFtools --SNPdensity output performed on Juneja converted VCF

library(ggplot2)
library(gridExtra)

snpdens<-read.table("inputs/vcftoolsoutputs/out.snpden", head=T)

jpeg("snpdensityGvBR1mb.jpg", units="in", 9,6, res=250)
ggplot(data=snpdens[which(snpdens$CHROM %in% c("chr1", "chr2", "chr3")),], aes(x=BIN_START/1000000, y=SNP_COUNT)) + 
  geom_bar(stat="identity") +
  facet_wrap(~CHROM, ncol=1) +
  theme_bw(base_size = 20) +
  scale_y_continuous("SNP count", breaks = seq(0,25,10), limits = c(0, 25)) +
  scale_x_continuous("Position along chromosome (mb)", limits = c(0, 375),breaks = seq(0,375,50), expand = c(0.01, 0.01))
dev.off()

