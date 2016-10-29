# Calculate rsquared values from a vcf. Returns object LD object from genetics package
# saves Rsquared and Pvalue matrix to file
# Requires the packages pegas and genetics
vcf2Rsquared<-function(vcffilepath, pop, chrom)
{
    startTime<-Sys.time()
    library(pegas)
    vcf<-read.vcf(vcffilepath)
    vcf[vcf=="./."]<-NA
    vcf[vcf=="."]<-NA
    detach("package:pegas", unload=TRUE)
    library(genetics)
    genotypes<-makeGenotypes(vcf)
    LDgeno<-LD(genotypes)
    filename<-paste0(pop,chrom,"RsquMatrix.txt")
    write.table(LDgeno$`R^2`, file = filename)
    filename<-paste0(pop,chrom,"PvalMatrix.txt")
    write.table(LDgeno$`P-value`, filename)
    return(LDgeno)
}

#Calculate distance between markers
CalcDist<-function(meltedLDmatrix, MappedTable)
{
  distances<-1:length(meltedLDmatrix[,1])
  cM<-1:length(meltedLDmatrix[,1])
  for (i in 1:length(meltedLDmatrix[,1]))
  {
    pos1<-(MappedTable$JunBP[which(MappedTable$ID==toString(meltedLDmatrix[i,1]))])
    pos2<-(MappedTable$JunBP[which(MappedTable$ID==toString(meltedLDmatrix[i,2]))])
    distances[i]<-abs(pos1 - pos2)
    cM1<-(MappedTable$cM[which(MappedTable$ID==toString(meltedLDmatrix[i,1]))])
    cM2<-(MappedTable$cM[which(MappedTable$ID==toString(meltedLDmatrix[i,2]))])
    cM[i]<-abs(cM1-cM2)
  }
  meltedLDmatrix<-cbind(meltedLDmatrix, distances, cM)
}

assignsChroms<-function(x){
  if (x %in% Chrom1IDs){
    return("1")
  } else if (x %in% Chrom2IDs) {
    return("2")
  } else if (x %in% Chrom3IDs) {
    return("3")
  } else {
    return(NA)
  }
}

MeanDecayLine<-function(ChromCalcTable, binsize){
  ChromCalcTable<-ChromCalcTable[with(ChromCalcTable, order(distances)),]
  maxdist<-max(ChromCalcTable$distances)
  numbins<-ceiling(maxdist/binsize)
  means<-1:numbins
  i=1
  for (i in 1:numbins){
    binnedLDs<-ChromCalcTable$LD[which((ChromCalcTable$distances > (i-1)*binsize) & (ChromCalcTable$distances < i*binsize))]
    means[i]<-mean(binnedLDs)
  }
  plottingbins<-seq(from=binsize/2, to=maxdist+(binsize/2), by=binsize)
  plotting<-data.frame(plottingbins, means)
}


LDR2matrix<-function(LDgeno)
{
  R2matrix<-LDgeno$`R^2`
  R2matrix[lower.tri(R2matrix)]<-t(R2matrix)[lower.tri(R2matrix)]
  return(R2matrix)
}

LDPvalmatrix<-function(LDgeno)
{
  PValmatrix<-LDgeno$`P-value`
  PValmatrix[lower.tri(PValmatrix)]<-t(PValmatrix)[lower.tri(PValmatrix)]
  return(PValmatrix)
}

# Order loci in a matrix
OrderMatrix<-function(matrix, mapTable)
{
  library(plyr)
  mapTable<-as.data.frame(mapTable)
  orderedLoci<-arrange(mapTable, JunBP)[,3]
  locusnames<-colnames(matrix)
  rownames(matrix)<-locusnames
  colnames(matrix)<-locusnames
  levels(orderedLoci)<-orderedLoci
  orderedMatrix<-matrix[match(orderedLoci, rownames(matrix)), match(orderedLoci, colnames(matrix))]
  return(orderedMatrix)
}

#arguments are two numeric vectors, distance can be physical or cM
#output is the fitted model object
#n is number of haplotypes, i.e. 2 * number of samples
LDdecayFit<-function(distance, LD.data, n)
{
  HW.st<-c(C=0.0001)
  HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))), 
                    start=HW.st, control=nls.control(maxiter=100))
  return(HW.nonlinear)
}

#input is a fitted nls model and the distances to predict
#output is a dataframe with predicted fitted points and upper and lower 95% conf interval

LDDecayFitSved<-function(distance, LD.data, n)
{
  beta.start<-c(beta=0.03)
  nonlinear<-nls(LD.data~((1/(1+(4*beta*distance)))) + 1/n, start=beta.start,control=nls.control(maxiter=100))
  return(nonlinear)
}

# Make predicted decay line
LDdecayPredict<-function(model, distances, n, popName)
{
  fpoints<-predict(model, distances)
  lowerC<-confint(model, 'C', level=0.95)[1]
  lowerpoints<-((10+lowerC*distances)/((2+lowerC*distances)*(11+lowerC*distances)))*(1+((3+lowerC*distances)*(12+12*lowerC*distances+(lowerC*distances)^2))/(n*(2+lowerC*distances)*(11+lowerC*distances)))
  upperC<-confint(model,'C', level=0.95)[2]
  upperpoints<-((10+upperC*distances)/((2+upperC*distances)*(11+upperC*distances)))*(1+((3+upperC*distances)*(12+12*upperC*distances+(upperC*distances)^2))/(n*(2+upperC*distances)*(11+upperC*distances)))
  pop<-rep(popName, length(distances))
  ld.df<-data.frame(distances, fpoints, lowerpoints, upperpoints, pop)
  ld.df<-ld.df[order(ld.df$distances),]
}

# Calculate Mean decay line
MeanDecayLine<-function(ChromCalcTable, binsize){
  ChromCalcTable<-ChromCalcTable[with(ChromCalcTable, order(distances)),]
  maxdist<-max(ChromCalcTable$distances)
  numbins<-ceiling(maxdist/binsize)
  means<-1:numbins
  i=1
  for (i in 1:numbins){
    binnedLDs<-ChromCalcTable$LD[which((ChromCalcTable$distances > (i-1)*binsize) & (ChromCalcTable$distances < i*binsize))]
    means[i]<-mean(binnedLDs)
  }
  plottingbins<-seq(from=binsize/2, to=maxdist+(binsize/2), by=binsize)
  plotting<-data.frame(plottingbins, means)
}

# Calculate folded MAF
minmaj<-function(x){
    if(x > 0.5)
    {
        return(1 - x)
    }
    else{
        return(x)
    }
}