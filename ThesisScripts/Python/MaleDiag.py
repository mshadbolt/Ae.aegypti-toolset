#!/usr/bin/env python

import glob
import sys
from subprocess import call
from collections import Counter
import multiprocessing
import time

start_time = time.time()

'''
Author: Marion Shadbolt
Email: marion.shadbolt@gmail.com

Last updated: 28/10/2016

Program Description:
Script to align pre-processed reads to male specific sequences to determine sex in Aedes aegypti NGS reads.
This script works with output files from the trim_align.sh bash script published by Rasic et al 2014.

Command line input:
1. Input directory of aligned reads. Should have sub-directories 'matched' and 'orphaned'
2. Path to population map file: tab delimited text file (same as provided to Stacks) saved in the input directory
Example use: python MaleDiag.py /path/to/input/files PopMap.txt or
python MaleDiag.py /path/to/input/files --single <PopName>
Number of processors can optionally be supplied as last commandline argument as an int

Output:
Script creates a subdirectory "MaleDiag" which it uses as a working directory
A summary tab-delimited txt file, 'MaleDiag.tsv' contains results of the diagnosis.

Dependencies:
Assumes bowtie aligner is installed and is in your PATH variable
Assumes same folder structure and file naming conventions as those outputted by trim_align.sh script
'''

class Sample:
    """
    Class to store all the relevant info about each sample processed.
    """
    def __init__(self, prefix, samplepop):
        self.name = prefix
        self.population = samplepop
        self.sex = "unknown"
        self.totalAligned = 0
        self.totalseqCount = {}
        self.seqCountPaired = {3526685: 0, 3076198: 0, 738613: 0}
        self.seqCountSingles = {3526685: 0, 3076198: 0, 738613: 0}

    def __repr__(self):
        return "\t".join([self.population, self.name, self.sex,
                          str(self.totalAligned),
                          str(self.seqCountPaired[3526685]), str(self.seqCountPaired[3076198]), str(self.seqCountPaired[738613]),
                          str(self.seqCountSingles[3526685]), str(self.seqCountSingles[3076198]), str(self.seqCountSingles[738613]),
                          ])

    def countReads(self):
        """Counts how many reads were aligned to each sequence."""
        self.totalAligned = sum(self.seqCountPaired.values()) + sum(self.seqCountSingles.values())
        self.totalseqCount = dict(Counter(self.seqCountPaired) + Counter(self.seqCountSingles))

    def diagnose(self):
        """Diagnoses whether the sample is Male or Female based on count of aligned sequences."""
        # incorporate read number threshold per sequence, default 4?
        myosex = 3526685
        if len(self.totalseqCount) == 0:
            self.sex = "Female"
        elif myosex in self.totalseqCount and len(self.totalseqCount) > 1:
            self.sex = "Male"
        else:
            self.sex = "unknown"

if len(sys.argv) < 3:
    print("Usage: <dir with subdirs with matched & orphaned trimmed reads> <pop map file txt> or --single <popname> ")
    sys.exit(0)

print("Commencing Male/Female diagnosis")

# Read the population map file

try:
    inputFilesPath = sys.argv[1]
    if not inputFilesPath.endswith("/"):
        inputFilesPath = inputFilesPath + "/"
    inputFiles = glob.glob(inputFilesPath + "matched/*")
    samples = []
    if sys.argv[2] != "--single":
        popmap = open(sys.argv[2], 'r')
        populations = {}
        for line in popmap:
            pre = line.split()[0]
            population = line.split()[1]
            populations[pre] = population
            samples.append(Sample(pre, population))
        populationsList = list(set(populations.values()))
        print(str(len(populationsList)) + " populations found in the map file. " )
        for pop in populationsList:
            popcount = 0
            for item in populations.items():
                if item[1] == pop:
                    popcount += 1
            print(str(pop) + ": " + str(popcount) + " individuals.")
        print
    elif sys.argv[2] == "--single":
        popname = sys.argv[3]
        samplenames = []
        for path in inputFiles:
            filename = path.split("/")[-1]
            sample = filename.split("_")[0]
            samplenames.append(sample)
        samplenames = set(samplenames)
        samples = [Sample(x, popname) for x in samplenames]
    else:
        print("You didnt specify a population, proceeding with anon")
        popname = "anon"
        for path in inputFiles:
            filename = path.split("/")[-1]
            sample = filename.split("_")[0]
            samples.append(Sample(sample, popname))

except IOError as e:
    print("Could not read pop map file (see below)!")
    print(e)
    sys.exit(1)

if len(inputFiles) == 0:
    print("No appropriate input files. Exiting.")
    sys.exit(1)

outputPath = inputFilesPath + "MaleDiag/"
call(['mkdir', outputPath])
bowtiedPath = outputPath + 'aligned'
call(['mkdir', bowtiedPath])
unalignedPath = outputPath+'unaligned'
call(['mkdir', unalignedPath])

# Arguments to supply to Bowtie
try:
    numberOfProcessors = int(sys.argv[-1])
    print("Using " + str(numberOfProcessors) + " processors.")
except ValueError as e:
    print("Number of processors not specified, using all available (" + str(multiprocessing.cpu_count()) + ").")
    numberOfProcessors = multiprocessing.cpu_count()
bowtie_path = "/home/working/bowtie"
bowtie_index_path = bowtie_path + "/indexes"
bowtieIndex = "Aaegypt-male"
min_size = 90
max_size = 400
maximum_number_of_mismatches_allowed = 3
maximum_number_of_reportable_alignments = 1

outputFile = open(outputPath + 'MaleFemaleDiagnosis.tsv', 'w')
outputFile.write("\t".join(["#population", "sample", "sex", "totalAligned",
                            "paired_myo", "paired_seq2", "paired_seq3",
                            "singles_myo", "singles_seq2", "singles_seq3"]) + "\n")
outputFile.close()

# Align pre-matched reads
print("Aligning reads ...")
counter = 0
totalSamples = len(samples)

for sample in samples:
    counter += 1
    try:
        # Align paired reads
        print("Aligning pairs from sample " + sample.name + " (sample " + str(counter) + " of " + str(totalSamples) + ")")
        
        # Unzip paired input files
        call(["gunzip", "-p", str(numberOfProcessors), inputFilesPath + "matched/" + sample.name + "_1.trimmed.matched.fq.gz"])
        call(["gunzip", "-p", str(numberOfProcessors), inputFilesPath + "matched/" + sample.name + "_2.trimmed.matched.fq.gz"])

        # Construct bowtie command
        bowtieCommandPaired = "-n " + str(maximum_number_of_mismatches_allowed) + \
            " -m " + str(maximum_number_of_reportable_alignments) + \
            " -p " + str(numberOfProcessors) + " -t --tryhard --chunkmbs 2048 --phred33 " + \
            "--minins " + str(min_size) + " --maxins " + str(max_size) + \
            " --al " + bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".paired.fq" \
            " --un " + unalignedPath + "/" + sample.name + "." + bowtieIndex + ".paired.fq" + \
            " -1 " + inputFilesPath + "matched/" + sample.name + "_1."+ "trimmed.matched.fq " \
            " -2 " + inputFilesPath + "matched/" + sample.name + "_2."+ "trimmed.matched.fq " \
            + bowtie_index_path + "/" + bowtieIndex + \
            " " + bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".paired.bowtie"
        # Align paired reads
        call("bowtie " + bowtieCommandPaired, shell=True)
        
        # Unzip orphan input files
        call(["gunzip", "-p", str(numberOfProcessors), inputFilesPath + "orphaned/" + sample.name + "_1.trimmed.orphan.fq.gz"])
        call(["gunzip", "-p", str(numberOfProcessors), inputFilesPath + "orphaned/" + sample.name + "_2.trimmed.orphan.fq.gz"])

        # Combine orphans and singles
        filenames = [unalignedPath + "/" + sample.name + "." + bowtieIndex + ".paired_1.fq",
                     unalignedPath + "/" + sample.name + "." + bowtieIndex + ".paired_2.fq",
                     inputFilesPath + "orphaned/" + sample.name + "_2.trimmed.orphan.fq",
                     inputFilesPath + "orphaned/" + sample.name + "_1.trimmed.orphan.fq"]
        singlesFile = unalignedPath + "/" + sample.name + "." + bowtieIndex + ".singles.fq"
        with open(singlesFile, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        # Align orphans and singles
        print("Aligning singles & orphans from sample " + sample.name + " (sample " + str(counter) + " of " +
              str(totalSamples) + ")")
        
        bowtieCommandSingles = " -n " + str(maximum_number_of_mismatches_allowed) + \
            " -m " + str(maximum_number_of_reportable_alignments) + \
            " -p " + str(numberOfProcessors) + " -t --tryhard --chunkmbs 2048 --phred33 " + \
            "--minins " + str(min_size) + " --maxins " + str(max_size) + \
            " --al " + bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".singles.fq " \
            " --un " + unalignedPath + "/" + sample.name + "." + bowtieIndex + ".Un.singles.fq " + \
            bowtie_index_path + "/" + bowtieIndex + " " + \
            singlesFile + " " + \
            bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".singles.bowtie"
        call("bowtie" + bowtieCommandSingles, shell=True)

        # Compress input files
        call("gzip -f -9 -p " + str(numberOfProcessors) + " " + inputFilesPath + "matched/*.fq", shell=True)
        call("gzip -f -9 -p " + str(numberOfProcessors) + " " + inputFilesPath + "orphaned/*.fq", shell=True)
        call("gzip -f -9 -p " + str(numberOfProcessors) + " " + bowtiedPath + "/*.fq", shell=True)

        # Remove unaligned files
        call("rm " + unalignedPath + "/*.fq", shell=True)

        print("Diagnosing " + sample.name + " (sample " + str(counter) + " of " + str(totalSamples) + ")")
        # count paired aligned reads
        bowtieFile = open(bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".paired.bowtie", "r")
        for line in bowtieFile:
            seq = line.split()[4]
            sample.seqCountPaired[int(seq)] += 1
        bowtieFile.close()

        # count single aligned reads
        bowtieFile = open(bowtiedPath + "/" + sample.name + "." + bowtieIndex + ".singles.bowtie", "r")
        for line in bowtieFile:
            seq = line.split()[4]
            sample.seqCountSingles[int(seq)] += 1
        bowtieFile.close()

        # call count and diagnose functions on the sample
        sample.countReads()
        sample.diagnose()
        # write to output file
        outputFile = open(outputPath + 'MaleFemaleDiagnosis.tsv', 'a')
        outputFile.write(str(sample) + "\n")
        outputFile.close()

        # Compress bowtie files
        print("Compressing bowties...")
        call("gzip -f -9 -p " + str(numberOfProcessors) + " " + bowtiedPath + "/*.bowtie", shell=True)

    except OSError as e:
        print(e)
        print("Could not find correct file or program, skipping.")
    except IOError as e:
        print(e)
        print("File not found, skipping.")
print("Diagnosis took: {0} minutes.".format(str(round((time.time() - start_time) / 60, 2))))
sys.exit(0)

