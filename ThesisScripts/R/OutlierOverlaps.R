library(VariantAnnotation)

VietColOutliers<-readVcf("inputs/vcftoolsoutputs/OutlierVcfs/juneja.VietColOutliers.recode.vcf", "AaegL1")
WColOutliers<-readVcf("inputs/vcftoolsoutputs/OutlierVcfs/juneja.WColBrWildOutliers.recode.vcf", "AaegL1")
#KDBrzWildOutliers<-readVcf("inputs/vcftoolsoutputs/OutlierVcfs/juneja.KDColBrWildOutliers.recode.vcf", "AaegL1")
KDWColOutliers<-readVcf("inputs/vcftoolsoutputs/OutlierVcfs/juneja.KDColWColOutliers.recode.vcf", "AaegL1")

# Find shared SNPs between VNcol and BRcol1, one hit
ViCWColIntersect<-findOverlaps(query = VietColOutliers, subject = WColOutliers, maxgap = 5000)
VietColOutliers@rowRanges[ViCWColIntersect@from]

# Find shared SNPs between VNcol and BRcol2, one hit
ViCKDIntersect<-findOverlaps(query = VietColOutliers, subject = KDWColOutliers, maxgap = 5000)
VietColOutliers@rowRanges[ViCKDIntersect@from]

# Find shared SNPs between BRcol1 and BRcol2, 17 hits
WColKDIntersect<-findOverlaps(query = WColOutliers, subject = KDWColOutliers, maxgap = 5000)
WColKDIntersectVCF<-subsetByOverlaps(query = WColOutliers, subject = KDWColOutliers, maxgap = 5000)
WColKDIntersectVCF@rowRanges

WColLociID<-names(WColOutliers@rowRanges[WColKDBrzWiIntersect@from])
as.character(WColOutliers@rowRanges@seqnames[WColKDBrzWiIntersect@from])
as.character(WColOutliers@rowRanges@ranges[WColKDBrzWiIntersect@from])


WColChrom<-as.character(WColOutliers@rowRanges@seqnames[WColKDBrzWiIntersect@from])
WKDColIntersect<-data.frame(WColChrom)
WKDColIntersect$WColPos<-as.character(WColOutliers@rowRanges@ranges[WColKDBrzWiIntersect@from])
WKDColIntersect$WColID<-WColLociID

WKDColIntersect$KDColPos<-as.character(KDBrzWildOutliers@rowRanges@ranges[WColKDBrzWiIntersect@to])
WKDColIntersect$KDColChrom<-as.character(KDBrzWildOutliers@rowRanges@seqnames[WColKDBrzWiIntersect@to])
WKDColIntersect$KDColID<-KDBrzWildLociID

KDBrzWildLociID<-names(KDBrzWildOutliers@rowRanges[WColKDBrzWiIntersect@to])

KDBrzWildLociID
WColLociID

library(UpSetR)
inputList<-list(BRcol1=names(WColOutliers@rowRanges), BRcol2 = names(KDWColOutliers@rowRanges), VNcol = names(VietColOutliers@rowRanges))
jpeg("FSTOutlierSNPintersections.jpg", units="in", width=8, height=6, res=300)
upset(fromList(inputList), main.bar.color = "steelblue", point.size = 8, text.scale = 2.5, mainbar.y.label = "SNP count")
dev.off()


require(ggplot2); require(plyr); require(gridExtra); require(grid);

between <- function(row, min, max){
  newData <- (row["ReleaseDate"] < max) & (row["ReleaseDate"] > min)
}

plot1 <- function(mydata, x){
  myplot <- (ggplot(mydata, aes_string(x= x, fill = "color"))
             + geom_histogram() + scale_fill_identity()
             + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

plot2 <- function(mydata, x, y){
  myplot <- (ggplot(data = mydata, aes_string(x=x, y=y, colour = "color"), alpha = 0.5)
             + geom_point() + scale_color_identity()
             + theme_bw() + theme(plot.margin = unit(c(0,0,0,0), "cm")))
}

attributeplots <- list(gridrows = 55,
                       plots = list(list(plot = plot1, x= "ReleaseDate",  queries = FALSE),
                                    list(plot = plot1, x= "ReleaseDate", queries = TRUE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = FALSE),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = TRUE)),
                       ncols = 3)

upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))

upset(movies, sets = c("Drama", "Comedy", "Action", "Thriller", "Western", "Documentary"),
      queries = list(list(query = intersects, params = list("Drama", "Action")),
                     list(query = between, params = list(1970, 1980), color = "red", active = TRUE)))

upset(movies, attribute.plots = attributeplots,
      queries = list(list(query = between, params = list(1920, 1940)),
                     list(query = intersects, params = list("Drama"), color= "red"),
                     list(query = elements, params = list("ReleaseDate", 1990, 1991, 1992))),
      main.bar.color = "yellow")
