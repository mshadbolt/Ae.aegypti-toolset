#!/usr/bin/env python

from __future__ import with_statement
from collections import defaultdict
import wget
import requests
import sys
import gzip
import itertools

"""
DESCRIPTION:
Program to download and create the 'geneset' file for input into Gowinda for the species Aedes aegypti from the uniprot
and BioMart (Vectorbase) databases. The program first downloads the GO terms for all Aedes aegypti genes directly from
Uniprot and BioMart. It then parses this file to create a dictionary of GO terms with associated genes. Finally it
downloads the obo file from the gene ontology consortium website, outputting a tab-delimited textfile in the format
required by Gowinda i.e.:
<GO Accession>\t<underscore separated GO name>\t<space separated string of associated genes>.

If you already have a obo file and a file containing a tab delimited file of geneids and GO terms, they can optionally
be supplied on the commandline as follows:
--gene_file uniprot-organism.tab.gz
--obo_file go.obo

Notes:
* User may need to install the wget module if they don't already have it installed: pip install wget (https://pypi.python.org/pypi/wget)
* Could be easily adapted for any other organism in the uniprot database by modifying the uniprot query url, see:
http://www.uniprot.org/help/programmatic_access
* To ensure you are accessing the latest VectorBase version visit: https://www.vectorbase.org/releases and update vb_version below.
    should be a four digit number YYMM, e.g. if latest release is October, 2016, version is 1610

Author: Marion Shadbolt, apart from the first two methods for obo file parsing (citation information below).
Last updated: 29/10/2016
Report any issues on github https://github.com/mshadbolt/Ae_aegypti-toolset
"""

"""
Description for following two methods:
A constant-space parser for the GeneOntology OBO v1.2 format

Version 1.0

From: https://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/

__author__    = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__   = "Apache v2.0"

"""
vb_version = 1610

def processGOTerm(goTerm):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(goTerm)  # Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.items():
        if len(value) == 1:
            ret[key] = value[0]
    return ret


def parseGOOBO(filename):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        filename: The filename to read
    """
    with open(filename, "r") as infile:
        currentGOTerm = None
        for line in infile:
            line = line.strip()
            if not line: continue  # Skip empty
            if line == "[Term]":
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = defaultdict(list)
            elif line == "[Typedef]":
                # Skip [Typedef sections]
                currentGOTerm = None
            else:  # Not [Term]
                # Only process if we're inside a [Term] environment
                if currentGOTerm is None: continue
                key, sep, val = line.partition(":")
                currentGOTerm[key].append(val.strip())
        # Add last term
        if currentGOTerm is not None:
            yield processGOTerm(currentGOTerm)


''' Start Program '''

clinput = sys.argv


if len(clinput) > 1 and "--gene_file" in clinput:
    print("Parsing given geneset file...")
    flagindex = clinput.index("--gene_file")
    genesetfile = clinput[flagindex + 1]
    try:
        uniprotData = []
        extension = str(genesetfile.split(".")[-1])
        if extension == "gz":
            with gzip.open(genesetfile, "rb") as f:
                for line in f.readlines():
                    uniprotData.append(line.decode('UTF-8'))
                    #print(lines)
        else:  # assumes regular plain text file
            genesetfile = open(genesetfile, "r")
            for line in genesetfile.readlines():
                uniprotData.append(line)
    except IOError as e:
        print("Could not read geneset file, check path or remove --gene_file to download from UniProt")
        print(e)
        sys.exit(0)
else:
    print("Downloading UniProt genes and associated GO terms...")
    try:
        urlquery = "http://www.uniprot.org/uniprot/?query=organism:7159&format=tab&compress=yes&columns=genes,go-id"
        uniprotfile = wget.download(urlquery)
        print(str(uniprotfile) + " downloaded successfully.")
    except IOError as e:
        print("Could not access file, check internet connection and/or url")
        print(e)
        sys.exit(0)
    uniprotData = []
    with gzip.open(uniprotfile, "rb") as f:
        for line in f.readlines():
            uniprotData.append(line.decode('UTF-8'))
try:
    print("Downloading BioMart gene and GO ids from VectorBase...")
    biomartquery = 'http://biomart.vectorbase.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "vb_mart_' + str(vb_version) +'" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.7" ><Dataset name = "aaegypti_eg_gene" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "go_accession" /></Dataset></Query>'
    biomartfile = requests.get(biomartquery)
    # biomartinput = wget.download(biomartquery)
    print(str("BioMart file downloaded successfully."))
except IOError as e:
    print("Could not access file, check internet connection and/or url")
    print(e)
    sys.exit(0)



GenetoGOdict = {}
print("Parsing UniProt file and storing in dictionary...")
for line in uniprotData:
    try:
        line = line.replace(';', ' ')
        line = line.replace('\n', '')
        splitline = line.split()
        genes = set()
        GOs = []
        for word in splitline:
            if word.startswith("AAEL"):
                genes.add(word)
            if word.startswith("GO:"):
                GOs.append(word)
        if len(GOs) > 0:
            for gene in genes:
                if gene in GenetoGOdict:
                    for goAcc in GOs:
                        if goAcc not in GenetoGOdict[gene]:
                            GenetoGOdict[gene].append(goAcc)
                else:
                    GenetoGOdict[str(gene)] = GOs
    except IndexError as e:
        print(e)
        continue

print("Number of Genes from UniProt " + str(len(GenetoGOdict)))
print("Parsing BioMart file and storing in dictionary...")
biomartlines = biomartfile.text.split("\n")
for line in biomartlines:
    line.replace("\\t", "\t")
    splitline = line.split()
    if len(splitline) > 1 and splitline[0].startswith("AAEL"):
        if splitline[0] in GenetoGOdict:
            GOs = splitline[1:]
            for goAcc in GOs:
                if goAcc not in GenetoGOdict[splitline[0]]:
                    GenetoGOdict[splitline[0]].append(goAcc)
        else:
            GenetoGOdict[splitline[0]] = splitline[1:]

print("Total Genes in dictionary " + str(len(GenetoGOdict)))
# get unique list of all GO IDs
GOids = list(GenetoGOdict.values())
GOids = set(itertools.chain(*GOids))

print("Transforming dictionary, one GO term to many genes...")
GOtoGenedict = {}
for item in GOids:
    GOtoGenedict[item] = set()
    for gene in GenetoGOdict.keys():
        if item in GenetoGOdict[gene]:
            GOtoGenedict[item].add(gene)

print("Inserting GO names into dictionary...")

if len(clinput) > 1 and "--obo_file" in clinput:
    try:
        print("Reading given obo file...")
        flagindex = clinput.index("--obo_file")
        obofile = clinput[flagindex + 1]
    except IOError as e:
        print("Could not read obo file, check file path.")
        print(e)
        sys.exit(1)
else:
    try:
        print("Downloading GO term names")
        obourl = "http://purl.obolibrary.org/obo/go.obo"
        obofile = wget.download(obourl)
        print("Obo file downloaded successfully.")
    except IOError as e:
        print("Problem downloading file, check internet connection")
        print(e)
        sys.exit(0)

GODict = {}
try:
    print("Parsing GO definitions file...")
    obos = open(obofile, "rb")
    for goTerm in parseGOOBO(obofile):
        GODict[goTerm["id"]] = goTerm["name"]
        if "alt_id" in goTerm.keys():
            if type(goTerm["alt_id"]) is list:
                for item in goTerm["alt_id"]:
                    GODict[item] = goTerm["name"]
            else:
                GODict[goTerm["alt_id"]] = goTerm["name"]
        if "is_a" in goTerm.keys() and goTerm["id"] in GOtoGenedict:
            if type(goTerm["is_a"]) is list:
                for item in goTerm["is_a"]:
                    try:
                        GOtoGenedict[item.split()[0]].union(GOtoGenedict[goTerm["id"]])
                    except KeyError:
                        GOtoGenedict[item.split()[0]] = GOtoGenedict[goTerm["id"]]
            else:
                try:
                    GOtoGenedict[goTerm["is_a"].split()[0]].union(GOtoGenedict[goTerm["id"]])
                except KeyError:
                    GOtoGenedict[goTerm["is_a"].split()[0]] = GOtoGenedict[goTerm["id"]]
        if "part_of" in goTerm.keys() and goTerm["id"] in GOtoGenedict:
            if type(goTerm["part_of"]) is list:
                for item in goTerm["part_of"]:
                    if goTerm["part_of"] in GOtoGenedict:
                        try:
                            GOtoGenedict[item.split()[0]].union(GOtoGenedict[goTerm["id"]])
                        except KeyError:
                            GOtoGenedict[item.split()[0]] = GOtoGenedict[goTerm["id"]]
            else:
                try:
                    GOtoGenedict[goTerm["part_a"].split()[0]].union(GOtoGenedict[goTerm["id"]])
                except KeyError:
                    GOtoGenedict[goTerm["part_a"].split()[0]] = GOtoGenedict[goTerm["id"]]
except IOError as e:
    print("Could not access file, check file input, internet connection and/or url")
    print(e)
    sys.exit(0)

print("Writing go terms to file...")
GowindaOutput = open("UniprotBioMartGOassocGowinda.txt", "w")
for k,v in GOtoGenedict.items():
    GOname = GODict[k]
    linetoprint = k + "\t" + GOname + "\t" + " ".join(v) + "\n"
    # GOname = GOname.replace(' ', '_')
    GowindaOutput.write(linetoprint)

print("Total GO terms: " + str(len(GOtoGenedict)))
print("done.")

sys.exit(0)
