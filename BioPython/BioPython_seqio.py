#!/usr/bin/env python3
# BioPython_seqio.py
"""Reads fasta file and outputs reverse complements of contents to new fasta file"""

import sys
from Bio import SeqIO

#Establish that command line must take in two arguments
if __name__ == "__main__":
    # Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count != 2:
        raise Exception("This script requires TWO arguments")
    FASTA_OLD = sys.argv[1]
    FASTA_NEW = sys.argv[2]

    #Iterate over records in old fasta file, generate reverse complement, and save as new output name
    for rec in SeqIO.parse(FASTA_OLD, "fasta"):
        records = [rec.reverse_complement(id="rc_" + rec.id, description="reverse complement")]
    SeqIO.write(records, FASTA_NEW, "fasta")



