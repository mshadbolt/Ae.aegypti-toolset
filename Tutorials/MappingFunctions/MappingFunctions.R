# Functions to perform mapping of markers to Juneja et al. 2014 assembly of Aedes aegypti genome
MapVCFLoci<-function(VCFfilepath)
{
  if(!file.exists("JunejaGeneticAssemblyR.csv")){
    stop("Can't find JunejaGeneticAssemblyR.csv, ensure it is in your working directory.")
  }
  if(!file.exists(VCFfilepath)){
    stop("Can't find VCF file. Ensure the filepath is correct.")
  }
  if(!require(VariantAnnotation)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("VariantAnnotation")
  }
  library(VariantAnnotation)
  origVcf<-readVcf(VCFfilepath, genome="AaegL1")
  assembly<-read.csv("JunejaGeneticAssemblyR.csv", stringsAsFactors = F)
  JunSC<-JunBP<-JunChr<-JuncM<-numLoci<-length(origVcf@rowRanges@seqnames)
  cat("Total Loci to be mapped: ", numLoci, "\n")
  CHROM<-as.vector(origVcf@rowRanges@seqnames)
  POS<-origVcf@rowRanges@ranges@start
  ID<-names(origVcf@rowRanges)
  #assignment statistics
  numToMisassemb=0
  startTime<-Sys.time()
  for(i in 1:numLoci)
  { 
    if(CHROM[i] %in% assembly$old.scaffold.name)
    {
      AssemblyRow<-assembly[which(assembly$old.scaffold.name==CHROM[i]),]
      AssemblyRow<-AssemblyRow[order(AssemblyRow$new.scaffold.name),]
      #If scaffold is misassembled, assign new scaffold locations
      if(AssemblyRow$misassembled[1]=="yes") 
      { 
        if(POS[i] <= AssemblyRow$end[1] && POS[i] >= AssemblyRow$start[1]) {
          JunChr[i]<-AssemblyRow$chr[1]
          JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + POS[i] - AssemblyRow[1,7])
          JunSC[i]<-toString(AssemblyRow[2,11])
          JuncM[i]<-AssemblyRow[1,2]
          numToMisassemb = numToMisassemb + 1
      } else if((POS[i] >= AssemblyRow$start[2]) && (POS[i] <= AssemblyRow$end[2])) {
          #assign locus 2nd row of assembly row
          JunChr[i]<-AssemblyRow$chr[2]
          JunBP[i]<-((AssemblyRow$cum.bp.chr[2] - AssemblyRow$length[2]) + POS[i] - AssemblyRow$start[2])
          JunSC[i]<-toString(AssemblyRow[2,11])
          JuncM[i]<-AssemblyRow[2,2]
          numToMisassemb = numToMisassemb + 1
      } else if((length(AssemblyRow[,1]) == 3) && (POS[i] >= AssemblyRow$start[3]) && (POS[i] <= AssemblyRow$end[3])) {
          JunChr[i]<-AssemblyRow$chr[3]
          JunBP[i]<-((AssemblyRow$cum.bp.chr[3] - AssemblyRow$length[3]) + POS[i] - AssemblyRow$start[3])
          JunSC[i]<-toString(AssemblyRow[3,11])
          JuncM[i]<-AssemblyRow[3,2]
          numToMisassemb = numToMisassemb + 1
      } else if((length(AssemblyRow[,1]) == 4) && (POS[i] >= AssemblyRow$start[4]) && (POS[i] <= AssemblyRow$end[4])) {
          JunChr[i]<-AssemblyRow$chr[4]
          JunBP[i]<-((AssemblyRow$cum.bp.chr[4] - AssemblyRow$length[4]) + POS[i] - AssemblyRow$start[4])
          JunSC[i]<-toString(AssemblyRow[4,11])
          JuncM[i]<-AssemblyRow[4,2]
          numToMisassemb = numToMisassemb + 1
      } else {
          JunChr[i]<-NA
          JunBP[i]<-NA
          JunSC[i]<-NA
          JuncM[i]<-NA
      }
    }
      #If not misassembled
    else 
    {
      if (AssemblyRow$misassembled[1]=="unknown")
        {
          JunChr[i]<-AssemblyRow$chr[1]
          JunBP[i]<-NA
          JunSC[i]<-toString(AssemblyRow[1,11])
          JuncM[i]<-AssemblyRow[1,2]
        }
      else
        {
          JunChr[i]<-AssemblyRow$chr[1]
          JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + POS[i])
          JunSC[i]<-toString(AssemblyRow[1,11])
          JuncM[i]<-AssemblyRow[1,2]
        }
      }
    }
    #if contig not found in assembly, assign unknown chrom and location
    else
    {
      JunChr[i]<-NA
      JunBP[i]<-NA
      JunSC[i]<-NA
      JuncM[i]<-NA
    }
  }
  #Store in matrix
  mappedLoci<-data.frame(CHROM, POS, ID, JunSC, JunChr, JunBP, JuncM)
  #Print assignment statistics
  cat("Total assigned to chromosomes:", length(which(!is.na(mappedLoci$JunChr))),"\n")
  cat("Total mapped to physical positions:", length(which(!is.na(mappedLoci$JunBP))), "\n")
  cat("Total assigned to misassembled contigs:", numToMisassemb,"\n")
  cat("Total unknown is", (length(which(mappedLoci[,5]=="unknown")) + 
                             length(which(is.na(mappedLoci$JunChr)))), "\n")
  cat("Total assigned to chrom 1:", length(which(mappedLoci[,5]==1)), "| chromosome only:", sum(is.na(mappedLoci$JunBP[which(mappedLoci$JunChr==1)])), "\n")
  cat("Total assigned to chrom 2:", length(which(mappedLoci[,5]==2)), "| chromosome only:", sum(is.na(mappedLoci$JunBP[which(mappedLoci$JunChr==2)])), "\n")
  cat("Total assigned to chrom 3:", length(which(mappedLoci[,5]==3)), "| chromosome only:", sum(is.na(mappedLoci$JunBP[which(mappedLoci$JunChr==3)])), "\n")
  cat("Time taken:", Sys.time()-startTime, "seconds.")
  return(mappedLoci)
}

### Function to map supercont locations to Juneja et al. 2014 chromosomes from a dataframe like that one outputted
# from vcftools --weir-fst option with columns CHROM, POS, and WEIR_AND_COCKERHAM_FST. 
# Appends additional columns JunChr, JunBP, JuncM, JunSC
# Can be used with any dataframe that contains the columns CHROM and POS with mapped AaegL1/2/3 superconts

MapDF<-function(df)
{
  startTime<-Sys.time()
  if(!file.exists("JunejaGeneticAssemblyR.csv")){
    stop("Can't find JunejaGeneticAssemblyR.csv, ensure it is in your working directory.")
  }
  assembly<-read.csv("JunejaGeneticAssemblyR.csv")
  numLoci<-nrow(df)
  cat("Total Loci to be mapped: ", numLoci, "\n")
  JuncM<-JunSC<-JunBP<-JunChr<-rep(NA, numLoci)
  df$CHROM<-as.character(df$CHROM)
  df<-data.frame(df, JunChr, JunBP, JunSC, JuncM)
  #assignment statistics
  numToMisassemb=0
  for(i in 1:numLoci)
  { 
    supercontmaxBP<-max(df$POS[df$CHROM==df$CHROM[i]])
    if(df$CHROM[i] %in% assembly$old.scaffold.name) {
      AssemblyRow<-assembly[which(assembly$old.scaffold.name==toString(df$CHROM[i])),]
      AssemblyRow<-AssemblyRow[order(AssemblyRow$new.scaffold.name),]
      #If scaffold is misassembled, assign new scaffold locations
      if(AssemblyRow$misassembled[1]=="yes") { 
        if(df$POS[i] <= AssemblyRow$end[1]) {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + df$POS[i] - AssemblyRow[1,7])
          df$JuncM[i]<-AssemblyRow[1,2]
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          numToMisassemb = numToMisassemb + 1
        } else if((df$POS[i] >= AssemblyRow$start[2]) && (df$POS[i] <= AssemblyRow$end[2])) {
          #assign locus 2nd row of assembly row
          df$JunChr[i]<-AssemblyRow$chr[2]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[2] - AssemblyRow$length[2]) + df$POS[i] - AssemblyRow$start[2])
          df$JuncM[i]<-AssemblyRow[2,2]
          df$JunSC[i]<-toString(AssemblyRow[2,11])
          numToMisassemb = numToMisassemb + 1
        } else if((length(AssemblyRow[,1]) == 3) && (df$POS[i] >= AssemblyRow$start[3]) && (df$POS[i] <= AssemblyRow$end[3])) {
          df$JunChr[i]<-AssemblyRow$chr[3]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[3] - AssemblyRow$length[3]) + df$POS[i] - AssemblyRow$start[3])
          df$JuncM[i]<-AssemblyRow[3,2]
          df$JunSC[i]<-toString(AssemblyRow[3,11])
          numToMisassemb = numToMisassemb + 1
        } else if((length(AssemblyRow[,1]) == 4) && (df$POS[i] >= AssemblyRow$start[4]) && (df$POS[i] <= AssemblyRow$end[4])) {
          df$JunChr[i]<-AssemblyRow$chr[4]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[4] - AssemblyRow$length[4]) + df$POS[i] - AssemblyRow$start[4])
          df$JuncM[i]<-AssemblyRow[4,2]
          df$JunSC[i]<-toString(AssemblyRow[4,11])
          numToMisassemb = numToMisassemb + 1
        } else {
          df$JunChr[i]<- NA
          df$JunBP[i]<- NA
          df$JunSC[i]<-NA
          df$JuncM[i]<-NA
        }
      } 
      
      else 
      {
        if (AssemblyRow$misassembled[1]=="unknown")
        {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-NA
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-AssemblyRow[1,2]
        }
        else
        {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + df$POS[i])
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-AssemblyRow[1,2]
        }
      }
    }
    #if contig not found in assembly, assign unknown chrom and location
    else
    {
      df$JunChr[i]<-NA
      df$JunBP[i]<-NA
      df$JunSC[i]<-NA
      df$JuncM[i]<-NA
    }
  }
  #Print assignment statistics
  cat("Total assigned to chromosomes:", length(which(!is.na(df$JunChr))),"\n")
  cat("Total mapped to physical positions:", length(which(!is.na(df$JunBP))), "\n")
  cat("Total assigned to misassembled contigs:", numToMisassemb,"\n")
  cat("Total unknown is", (length(which(df[,5]=="unknown")) + 
                             length(which(is.na(df$JunChr)))), "\n")
  cat("Total physically assigned to chrom 1:", length(which(df$JunChr==1)), "| chromosome only:", sum(is.na(df$JunBP[which(df$JunChr==1)])), "\n")
  cat("Total physically assigned to chrom 2:", length(which(df$JunChr==2)), "| chromosome only:", sum(is.na(df$JunBP[which(df$JunChr==2)])), "\n")
  cat("Total physically assigned to chrom 3:", length(which(df$JunChr==3)), "| chromosome only:", sum(is.na(df$JunBP[which(df$JunChr==3)])), "\n")
  cat("Time taken:", Sys.time()-startTime, "seconds.")
  return(df)
}

### Function to map supercont locations to Juneja et al. 2014 chromosomes from a dataframe like that one outputted
# from vcftools --weir-fst option with columns CHROM, POS, and WEIR_AND_COCKERHAM_FST. 
# Appends additional columns JunChr, JunBP, JuncM, JunSC
# Can be used with any dataframe that contains the columns CHROM and POS with mapped AaegL1/2/3 superconts
# Specifically designed for visualisation rather than analysis. Markers mapped to chromosomes 1, 2 and 3 without a physical position
# are assigned to false chromosomes 4, 5 and 6 respectively. All unassigned ccontigs are assigned to chromosome 7.
# Unassigned contigs are spaced arbitrarily and spaced by the integer value 'spacing', default set at 5kb.

MapDF_FST<-function(df, spacing=5000)
{
  startTime<-Sys.time()
  if(!file.exists("JunejaGeneticAssemblyR.csv")){
    stop("Can't find JunejaGeneticAssemblyR.csv, ensure it is in your working directory.")
  }
  assembly<-read.csv("JunejaGeneticAssemblyR.csv")
  numLoci<-nrow(df)
  cat("Total Loci to be mapped: ", numLoci, "\n")
  JuncM<-JunSC<-JunBP<-JunChr<-rep(NA, numLoci)
  df$CHROM<-as.character(df$CHROM)
  df<-data.frame(df, JunChr, JunBP, JunSC, JuncM)
  cumulchrom<-c(1000000,1000000,1000000,1000000)
  supercontmaxBP<-0
  #assignment statistics
  numToMisassemb=0
  for(i in 1:numLoci)
  { 
    supercontmaxBP<-max(df$POS[df$CHROM==df$CHROM[i]])
    if(df$CHROM[i] %in% assembly$old.scaffold.name) {
      AssemblyRow<-assembly[which(assembly$old.scaffold.name==toString(df$CHROM[i])),]
      AssemblyRow<-AssemblyRow[order(AssemblyRow$new.scaffold.name),]
      #If scaffold is misassembled, assign new scaffold locations
      if(AssemblyRow$misassembled[1]=="yes") { 
        if(df$POS[i] <= AssemblyRow$end[1]) {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + df$POS[i] - AssemblyRow[1,7])
          df$JuncM[i]<-AssemblyRow[1,2]
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          numToMisassemb = numToMisassemb + 1
        } else if((df$POS[i] >= AssemblyRow$start[2]) && (df$POS[i] <= AssemblyRow$end[2])) {
          #assign locus 2nd row of assembly row
          df$JunChr[i]<-AssemblyRow$chr[2]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[2] - AssemblyRow$length[2]) + df$POS[i] - AssemblyRow$start[2])
          df$JuncM[i]<-AssemblyRow[2,2]
          df$JunSC[i]<-toString(AssemblyRow[2,11])
          numToMisassemb = numToMisassemb + 1
        } else if((length(AssemblyRow[,1]) == 3) && (df$POS[i] >= AssemblyRow$start[3]) && (df$POS[i] <= AssemblyRow$end[3])) {
          df$JunChr[i]<-AssemblyRow$chr[3]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[3] - AssemblyRow$length[3]) + df$POS[i] - AssemblyRow$start[3])
          df$JuncM[i]<-AssemblyRow[3,2]
          df$JunSC[i]<-toString(AssemblyRow[3,11])
          numToMisassemb = numToMisassemb + 1
        } else if((length(AssemblyRow[,1]) == 4) && (df$POS[i] >= AssemblyRow$start[4]) && (df$POS[i] <= AssemblyRow$end[4])) {
          df$JunChr[i]<-AssemblyRow$chr[4]
          df$JunBP[i]<-((AssemblyRow$cum.bp.chr[4] - AssemblyRow$length[4]) + df$POS[i] - AssemblyRow$start[4])
          df$JuncM[i]<-AssemblyRow[4,2]
          df$JunSC[i]<-toString(AssemblyRow[4,11])
          numToMisassemb = numToMisassemb + 1
        } else {
          if(df$CHROM[i] == df$CHROM[i - 1]) {
            cumulchrom[4] <- cumulchrom[4] + df$POS[i]
            df$JunBP[i]<-cumulchrom[4]
          } else {
            cumulchrom[4] <- cumulchrom[4] + spacing
            df$JunBP[i]<-cumulchrom[4]
          }
          df$JunChr[i]<- 7
          df$JunSC[i]<-df$CHROM[i]
          df$JuncM[i]<-NA
        }
      } else {
        if (AssemblyRow$misassembled[1] == "unknown") {
          df$JunChr[i]<-AssemblyRow$chr[1] + 3
          if(df$CHROM[i] == df$CHROM[i - 1]) {
            cumulchrom[AssemblyRow$chr[1]] <- cumulchrom[AssemblyRow$chr[1]] + df$POS[i]
            df$JunBP[i]<-cumulchrom[AssemblyRow$chr[1]]
          } else {
            cumulchrom[AssemblyRow$chr[1]] <- cumulchrom[AssemblyRow$chr[1]] + spacing
            df$JunBP[i]<-cumulchrom[AssemblyRow$chr[1]]
          }
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-NA
        } else {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + df$POS[i])
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-AssemblyRow[1,2]
        }
      }
    } else {
      if(df$CHROM[i] == df$CHROM[i - 1]) {
        cumulchrom[4] <- cumulchrom[4] + df$POS[i]
        df$JunBP[i]<-cumulchrom[4]
      } else {
        cumulchrom[4] <- cumulchrom[4] + spacing
        df$JunBP[i]<-cumulchrom[4]
      }
      df$JunChr[i]<-7
      df$JunSC[i]<-df$CHROM[i]
      df$JuncM[i]<-NA
    }
  }
  # #Print assignment statistics
  cat("Total assigned:", length(which(df[,4] != is.na(NA))),"\n")
  cat("Total assigned to misassembled contigs:", numToMisassemb,"\n")
  cat("Total unknown is", length(which(is.na(df$JunBP))), "\n")
  cat("Total assigned to chrom 1:", length(which(df$JunChr==1)), "\n")
  cat("Total assigned to chrom 2:", length(which(df$JunChr==2)), "\n")
  cat("Total assigned to chrom 3:", length(which(df$JunChr==3)), "\n")
  cat("Time taken:", Sys.time()-startTime, "seconds.")
  return(df)
}