#create separate chroms vcf files
MapVCFLoci<-function(VCFfilepath)
{
  if(!require(VariantAnnotation)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("VariantAnnotation")
  }
  library(VariantAnnotation)
  origVcf<-readVcf(VCFfilepath, genome="AaegL1")
  assembly<-read.csv("JunejaGeneticAssemblyR.csv")
  supercontigJun<-JunBP<-JunChr<-cM<-numLoci<-length(origVcf@rowRanges@seqnames)
  supercontig<-as.vector(origVcf@rowRanges@seqnames)
  Pos<-origVcf@rowRanges@ranges@start
  ID<-names(origVcf@rowRanges)
  #assignment statistics
  numToMisassemb=0
  startTime<-Sys.time()
  for(i in 1:numLoci)
  { 
    if(supercontig[i] %in% assembly$old.scaffold.name)
    {
      AssemblyRow<-assembly[which(assembly$old.scaffold.name==supercontig[i]),]
      AssemblyRow<-AssemblyRow[order(AssemblyRow$new.scaffold.name),]
      #If scaffold is misassembled, assign new scaffold locations
      if(AssemblyRow$misassembled[1]=="yes") 
      { 
        if(Pos[i] <= AssemblyRow$end[1]) {
          JunChr[i]<-toString(AssemblyRow$chr[1])
          JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + Pos[i] - AssemblyRow[1,7])
          supercontigJun[i]<-supercontig[i]
          cM[i]<-AssemblyRow[1,2]
          numToMisassemb = numToMisassemb + 1
      } else if((Pos[i] >= AssemblyRow$start[2]) && (Pos[i] <= AssemblyRow$end[2])) {
          #assign locus 2nd row of assembly row
          JunChr[i]<-toString(AssemblyRow$chr[2])
          JunBP[i]<-((AssemblyRow$cum.bp.chr[2] - AssemblyRow$length[2]) + Pos[i] - AssemblyRow$start[2])
          supercontigJun[i]<-toString(AssemblyRow[2,11])
          cM[i]<-AssemblyRow[2,2]
          numToMisassemb = numToMisassemb + 1
      } else if((length(AssemblyRow[,1]) == 3) && (Pos[i] >= AssemblyRow$start[3]) && (Pos[i] <= AssemblyRow$end[3])) {
          JunChr[i]<-toString(AssemblyRow$chr[3])
          JunBP[i]<-((AssemblyRow$cum.bp.chr[3] - AssemblyRow$length[3]) + Pos[i] - AssemblyRow$start[3])
          supercontigJun[i]<-toString(AssemblyRow[3,11])
          cM[i]<-AssemblyRow[3,2]
          numToMisassemb = numToMisassemb + 1
      } else if((length(AssemblyRow[,1]) == 4) && (Pos[i] >= AssemblyRow$start[4]) && (Pos[i] <= AssemblyRow$end[4])) {
          JunChr[i]<-toString(AssemblyRow$chr[4])
          JunBP[i]<-((AssemblyRow$cum.bp.chr[4] - AssemblyRow$length[4]) + Pos[i] - AssemblyRow$start[4])
          supercontigJun[i]<-toString(AssemblyRow[4,11])
          cM[i]<-AssemblyRow[4,2]
          numToMisassemb = numToMisassemb + 1
      } else {
          JunChr[i]<-"unknown"
          JunBP[i]<-NA
          supercontigJun[i]<-NA
          cM[i]<-NA
      }
    }
      #If not misassembled
    else 
    {
      if (AssemblyRow$misassembled[1] == "unknown")
        {
          JunChr[i]<-toString(AssemblyRow$chr[1])
          JunBP[i]<-NA
          supercontigJun[i]<-toString(AssemblyRow[1,11])
          cM[i]<-AssemblyRow[1,2]
        }
      else
        {
          JunChr[i]<-toString(AssemblyRow$chr[1])
          JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + Pos[i])
          supercontigJun[i]<-toString(AssemblyRow[1,11])
          cM[i]<-AssemblyRow[1,2]
        }
      }
    }
    #if contig not found in assembly, assign unknown chrom and location
    else
    {
      JunChr[i]<-NA
      JunBP[i]<-NA
      supercontigJun[i]<-NA
      cM[i]<-NA
    }
  }
  #Store in matrix
  mappedLoci<-data.frame(supercontig, Pos, ID, supercontigJun, JunChr, JunBP, cM)
  #Create Loci IDs per chrom textfiles
  write.table(mappedLoci[which(mappedLoci[,5]=="1"), 3], "Chrom1Loci.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(mappedLoci[which(mappedLoci[,5]=="2"), 3], "Chrom2Loci.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(mappedLoci[which(mappedLoci[,5]=="3"), 3], "Chrom3Loci.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
  #Print assignment statistics
  cat("Total assigned:", length(which(mappedLoci[,4] != is.na(NA))),"\n")
  cat("Total assigned to misassembled contigs:", numToMisassemb,"\n")
  cat("Total unknown is", length(which(mappedLoci[,5]=="unknown")), "\n")
  cat("Total assigned to chrom 1:", length(which(mappedLoci[,5]==1)), "\n")
  cat("Total assigned to chrom 2:", length(which(mappedLoci[,5]==2)), "\n")
  cat("Total assigned to chrom 3:", length(which(mappedLoci[,5]==3)), "\n")
  cat("Time taken:", Sys.time()-startTime, "seconds.")
  return(mappedLoci)
}