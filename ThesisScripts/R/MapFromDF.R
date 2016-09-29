### Function to map supercont locations to Juneja et al. 2014 chromosomes from a dataframe like that one outputted
# from vcftools --weir-fst option with columns CHROM, POS, and WEIR_AND_COCKERHAM_FST. 
# Appends additional columns JunChr, JunBP, JuncM, JunSC
# Can be used with any dataframe that contains the columns CHROM and POS with mapped AaegL1/2/3 superconts

MapDF<-function(df)
{
  startTime<-Sys.time()
  assembly<-read.csv("JunejaGeneticAssemblyR.csv")
  numLoci<-nrow(df)
  cat("Total Loci to be mapped: ", numLoci, "\n")
  JuncM<-JunSC<-JunBP<-JunChr<-rep(NA, numLoci)
  df<-data.frame(df, JunChr, JunBP, JunSC, JuncM)
  #assignment statistics
  numToMisassemb=0
  for(i in 1:numLoci)
  { 
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
          df$JunChr[i]<-NA
          df$JunBP[i]<-NA
          df$JunSC[i]<-NA
          df$JuncM[i]<-NA
        }
      } else {
        if (AssemblyRow$misassembled[1] == "unknown") {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-NA
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-NA
        } else {
          df$JunChr[i]<-AssemblyRow$chr[1]
          df$JunBP[i]<-((AssemblyRow[1, 9] - AssemblyRow[1,8]) + df$POS[i])
          df$JunSC[i]<-toString(AssemblyRow[1,11])
          df$JuncM[i]<-AssemblyRow[1,2]
        }
      }
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