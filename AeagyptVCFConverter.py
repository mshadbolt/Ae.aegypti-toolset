#!/usr/bin/env python

import sys
import csv

'''
Script to convert a vcf file aligned to AaegL1-3 i.e. supercont1.### into chromosomes according to
Juneja et al. 2014 assembly.
Expected input: Standard VCF file where header rows commence with # and the first two columns are:
    CHROM: i.e. supercont1.*
    POS: basepair position of SNP
Output is a copied vcf file prefixed with juneja, i.e. juneja.<inputfilename>
SNP loci are sorted in order of chromosome name, then bp position
Usage: python <inputfilename>
Dependencies: Needs to be run in the same directory as the genetic map assembly file in .csv format.
Expected filename: "JunejaGeneticAssembly.csv"
This file is from Juneja et al. 2014 Table S2 http://dx.doi.org/10.1371/journal.pntd.0002652

Author: Marion Shadbolt
Last updated: 23rd September 2016
email: marion.shadbolt@gmail.com
'''

class Converter:
    '''
    Class to build and store contig information from the Juneja 2014 et al. Ae. aegypti assembly
    '''
    def __init__(self):
        self.superContDict = {}

    def readContigs(self):
        print("Reading Juneja Assembly input file and constructing contig dictionary...")
        with open("JunejaGeneticAssembly.csv", "r") as mapfile:
            csvreader = csv.reader(mapfile, delimiter=',')
            cols = [0, 4, 5, 6, 7, 8]
            for row in csvreader:
                if not row[0].startswith("#"):
                    if row[3] not in self.superContDict.keys():
                        self.superContDict[row[3]] = [[row[i] for i in cols]]
                    else:
                        temp = self.superContDict[row[3]]
                        add = [row[i] for i in cols]
                        temp.insert(len(temp) - 1, add)
                        self.superContDict[row[3]] = temp

        print("done.")

    
    def convert(self, sc, bp):
        '''
        Method to convert supercont and pos information to Chromosome and pos information
        If location not found in the assembly, original information is returned
        '''
        try:
            assemblyInfo = self.superContDict[sc]
            if len(self.superContDict[sc]) > 1:
                for i in range(len(self.superContDict[sc])):
                    if bp >= int(assemblyInfo[i][2]) and bp <= int(assemblyInfo[i][3]):
                        chrombp = bp - int(assemblyInfo[i][2]) + int(assemblyInfo[i][5]) - int(assemblyInfo[i][4])
                        chrom = assemblyInfo[i][0]
                        return ((chrom, chrombp))
            elif int(assemblyInfo[0][3]) >= bp >= int(assemblyInfo[0][2]):
                chrombp = bp + int(assemblyInfo[0][5]) - int(assemblyInfo[0][4])
                chrom = assemblyInfo[0][0]
                return ((chrom, chrombp))
            else:
                return ((sc, bp))
            return ((sc, bp))
        except KeyError:
            return ((sc, bp))


if len(sys.argv) < 2:
    print("Check commandline input, usage: <vcf file>")
    sys.exit(1)

try:
    # Create the converter object and read assembly file
    PosConverter = Converter()
    PosConverter.readContigs()
    # Check appropriate input provided
    inputfilename = sys.argv[1]
    filenamecheck = inputfilename.split(".")
    if filenamecheck[len(filenamecheck) - 1].lower() != "vcf":
        print("Wrong input file, please use a .vcf file. Exiting.")
        sys.exit(0)
    convertedfilename = "juneja." + inputfilename
    print("Converting positions and writing to " + convertedfilename)
    # Open inputfile, read each line and write to output file with converted positions
    inputVcf = open(inputfilename, "r")
    outputVcf = open(convertedfilename, "w")
    linesToWrite = []
    for line in inputVcf:
        if line[0] is "#":
            outputVcf.write(line)
        else:
            linelist = line.split("\t")
            converted = PosConverter.convert(linelist[0], int(linelist[1]))
            linelist[0] = converted[0]
            linelist[1] = converted[1]
            linesToWrite.append(linelist)
    inputVcf.close()
    linesToWrite.sort()
    for item in linesToWrite:
        item[1] = str(item[1])
        outputVcf.write("\t".join(item))
    outputVcf.close()
    print("done.")
except IOError as e:
    print("Problem with reading input files, exiting see below: ")
    print(e)
    sys.exit(1)
sys.exit(0)
