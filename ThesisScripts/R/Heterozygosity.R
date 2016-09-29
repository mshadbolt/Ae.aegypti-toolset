mycols=c("#117733", "#88CCEE", "#332288", "#882255", "#40AA99") 
GvVi012<-read.table("inputs/vcftoolsoutputs/Aus27Viet.012")
GvVi012<-GvVi012[,-1]
GvVi012[GvVi012==-1]<-NA

calcHetero<-function(x) {
  x<-x[!is.na(x)]
  return(sum(x == 1, na.rm=T)/length(x))
}

GvViHeteros<-apply(GvVi012, MARGIN = 1, calcHetero)

indivs<-read.table("inputs/vcftoolsoutputs/Aus27Viet.012.indv")
pops<-c(rep("VietWi", 27), rep("VietCol", 27), rep("AusWild", 27))
Hetero.df<-data.frame(pops, GvViHeteros)
colnames(Hetero.df)<-c("Pop", "Hetero")

GvBr012<-read.table("inputs/vcftoolsoutputs/Aus20Braz.012")
GvBr012<-GvBr012[,-1]
GvBr012[GvBr012==-1]<-NA

GvBrHeteros<-apply(GvBr012, MARGIN = 1, calcHetero)
pops<-c(rep("AusWild", 20), rep("WCol", 20), rep("KDCol", 20), rep("BrWild", 20))
GvBrHeteros<-data.frame(pops, GvBrHeteros)
colnames(GvBrHeteros)<-c("Pop", "Hetero")
Hetero.df<-rbind(Hetero.df, GvBrHeteros)
Hetero.df$Pop<-factor(Hetero.df$Pop, levels=c("VietWi", "BrWild", "VietCol", "KDCol",  "WCol", "AusWild"))
Hetero.df$Dataset<-c(rep("AusVi", 81), rep("AusBraz", 80))

labels = c(AusWild="AUwild", KDCol="BRcol2", WCol="BRcol1", BrWi="BRwild")
jpeg("AusBrazHetero.jpeg", 6,8, units="in", res=300)
ggplot(data=Hetero.df[which(Hetero.df$Dataset=="AusBraz"),], aes(y=Hetero, x = Pop, color=as.factor(Pop))) + 
  geom_boxplot(width = 0.6, outlier.shape = 8, outlier.colour = "white") + 
  geom_jitter(height = 0, width = 0.4) +
  coord_flip() +
  scale_color_manual(values = c( "#882255", "#40AA99","#117733", "#332288")) +
  theme_bw() +
  scale_x_discrete("Population", labels = c("BRwild", "BRcol2", "BRcol1", "AUwild")) +
  scale_y_continuous("Heterozygosity", limits = c(0.05, 0.2)) + 
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
dev.off()

jpeg("AusVietHetero.jpeg", 6,7, units="in", res=300)
ggplot(data=Hetero.df[which(Hetero.df$Dataset=="AusVi"),], aes(y=Hetero, x = Pop, color=as.factor(Pop))) + 
  geom_boxplot(width = 0.6,outlier.shape = 8, outlier.colour = "white") + 
  geom_jitter(height = 0, width = 0.4) +
  coord_flip() +
  scale_color_manual(values = c("#882255", "#40AA99", "#332288")) +
  theme_bw() +
  scale_x_discrete("Population", labels=c("VNwild","VNcol","AUwild")) +
  scale_y_continuous("Heterozygosity", limits = c(0.05, 0.2)) + 
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_blank())
dev.off()


