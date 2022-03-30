#!/usr/bin/env python3
# sliding_window_fasta.py
"""Same as in module 1, but with SeqIO"""

import sys
from Bio import SeqIO


def sliding_window_fasta(kmer, fasta):
    '''Returns a list of kmers in the provided string size'''
    # Initialize an empty list
    dna = []

    # Create variable to determine last item in window
    last = len(fasta) - int(kmer) + 1
    # Iterate over string and return dna
    for first in range(0, last):
        dna.append(fasta[first:first + int(kmer)])
    return dna


def gc_content(dna):
    ''' Returns [0, 1], the fraction of GCs in the given string'''
    # Make the sequence lowercase
    dna = dna.lower()

    # Count the number of g's and c's
    gc = 0
    for nucleotide in dna:
        if nucleotide in ['g', 'c']:
            gc += 1
    return gc / len(dna)


if __name__ == "__main__":
    # Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires at least 2 arguments: kmer and fasta file")

    kmer = int(sys.argv[1])
    filehandle = sys.argv[2]
    # seq = ''
    # with open(filehandle, 'r') as fasta:
    #     for line in fasta:
    #         line = line.rstrip()
    #         if line.startswith('>'):
    #             print(line)
    #         else:
    #             seq += line

    #Open file with SeqIO instead
    for seq_record in SeqIO.parse(filehandle, "fasta"):
        print(">{}".format(seq_record.description))
        for string in sliding_window_fasta(kmer, seq_record.seq):
            print("{}\t{:.2f}".format(string, gc_content(string)))


