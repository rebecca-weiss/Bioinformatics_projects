#!/usr/bin/env python3
# BioPython_seq.py
"""Creates a SeqRecord object and writes file to BioPython_seq.gb"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO


#Converts string to upper case
SEQ_STR = "aaaatgggggggggggccccgtt".upper()

#Creates Seq object with DNAAlphabet from generic_dna
SEQ_OBJ = Seq(SEQ_STR, generic_dna)

# Creates SeqRecord
SEQ_EX1 = SeqRecord(SEQ_OBJ, id="#12345", description="example 1")

# Use SeqIO to write to gb file
SeqIO.write(SEQ_EX1, "BioPython_seq.gb", "genbank")
