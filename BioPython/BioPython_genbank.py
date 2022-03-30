#!/usr/bin/env python3
# BioPython_genbank.py
"""Creates a list with given Seq objects and prints out sequences, type,
location and strand of each feature"""

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "weiss.rebe@northeastern.edu"

SEQ_LIST = []
# Fetch the seq object and write them to sequence list
with Entrez.efetch(db="nucleotide", retmode="text", rettype="gb", id="515056") as handle:
    SEQ_LIST.append(SeqIO.read(handle, "gb"))

with Entrez.efetch(db="nucleotide", retmode="text", rettype="gb", id="J01673.1") as handle:
    SEQ_LIST.append(SeqIO.read(handle, "gb"))

# collect feature and print for each sequence
for item in SEQ_LIST:
    print("Sequence: {}".format(item.name))
    for feature in item.features:
        print("Feature: {}, Location: {}, Strand: {}".format(feature.type, feature.location, feature.strand))


