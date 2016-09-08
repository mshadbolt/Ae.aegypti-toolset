#!/usr/bin/env python

from __future__ import with_statement
from collections import defaultdict
import wget
import sys
import gzip
import itertools

"""
DESCRIPTION:
Program to download and create the 'geneset' file for input into Gowinda for the species Aedes aegypti from the uniprot
database. The program first downloads the GO terms for all Aedes aegypti genes directly from Uniprot. It then parses
this file to create a dictionary of GO terms with associated genes. Finally it downloads the obo file from the gene
ontology consortium website, outputting a tab-delimited textfile in the format required by Gowinda i.e.:
<GO Accession>\t<underscore separated GO name>\t<space separated string of associated genes>.

Could be easily adapted for any other organism in the uniprot database

Author: Marion Shadbolt, apart from the first two methods for obo file parsing (citation information below).
Last updated: 08/09/2016
"""

"""
Description for following two methods:
A constant-space parser for the GeneOntology OBO v1.2 format

Version 1.0

From: https://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/

"""
__author__    = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__   = "Apache v2.0"

def processGOTerm(goTerm):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
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
            if not line: continue #Skip empty
            if line == "[Term]":
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = defaultdict(list)
            elif line == "[Typedef]":
                #Skip [Typedef sections]
                currentGOTerm = None
            else: #Not [Term]
                #Only process if we're inside a [Term] environment
                if currentGOTerm is None: continue
                key, sep, val = line.partition(":")
                currentGOTerm[key].append(val.strip())
        #Add last term
        if currentGOTerm is not None:
            yield processGOTerm(currentGOTerm)

print("Downloading uniprot geneids to GO terms...")
try:
    urlquery = "http://www.uniprot.org/uniprot/?query=organism:7159&format=tab&compress=yes&columns=entry%20name,database(vectorbase),go-id"
    filename = wget.download(urlquery)
    print(str(filename) + " downloaded successfully.")
except IOError as e:
    print("Could not access file, check internet connection and/or url")
    print(e)
    sys.exit(0)

lines = []
with gzip.open(filename, "rb") as f:
    for line in f.readlines():
        lines.append(line.decode('UTF-8'))

GenetoGOdict = {}
print("Parsing file and storing in dictionary...")
for line in lines[1:]:
    try:
        line = line.replace(';', '')
        line = line.replace('\n', '')
        splitline = line.split()
        if splitline[1].startswith("AAEL") and splitline[2] is not '':
            GenetoGOdict[splitline[1]] = splitline[2:]
    except IndexError:
        continue

GOids = list(GenetoGOdict.values())
GOids = set(itertools.chain(*GOids))

print("Tranforming dictionary, one GO term to many genes...")
GOtoGenedict = {}
for item in GOids:
    GOtoGenedict[item] = []
    for gene in GenetoGOdict.keys():
        if item in GenetoGOdict[gene]:
            GOtoGenedict[item].append(gene)

print("Downloading gene ontology definitions...")
try:
    obourl = "http://purl.obolibrary.org/obo/go.obo"
    obofile = wget.download(obourl)
    print(str(obofile) + " downloaded successfully.")
except IOError as e:
    print("Could not access file, check internet connection and/or url")
    print(e)
    sys.exit(0)

GODict = {}

print("Parsing definitions file...")
for goTerm in parseGOOBO(obofile):
    GODict[goTerm["id"]] = goTerm["name"]
    if "alt_id" in goTerm.keys():
        if type(goTerm["alt_id"]) is list:
            for item in goTerm["alt_id"]:
                GODict[item] = goTerm["name"]
        else:
            GODict[goTerm["alt_id"]] = goTerm["name"]

print("Inserting GO names into dictionary...")

for item in GOtoGenedict.keys():
    GOname = GODict[item]
    GOname = GOname.replace(' ', '_')
    GOtoGenedict[item].insert(0, GOname)

print("Writing go terms to file...")
GowindaOutput = open("UniprotGOassocGowinda.txt", "w")
for k, v in GOtoGenedict.items():
    linetoprint = k + "\t" + v[0] + "\t" + " ".join(v[1:]) + "\n"
    GowindaOutput.write(linetoprint)
print("done.")
sys.exit(0)
