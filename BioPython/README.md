# BINF6309 Module 2: BioPython

## Author: Rebecca Weiss 

## Overview
./basic_classes.py

Creates a class for "Circle" and "GraduateStudent" with its own attributes and methods. Each class has two examples
listed and will print output of those examples requested.


./BioPython_seq.py 

Creates a SeqRecord object based on provided parameters and then writes the output to a genbank file: "BioPython_seq.gb"


./BioPython_genbank.py

Creates a list with the Seq objects for 2 ids given, uses Entrez efetch, and prints out the sequence along with the type, 
location, and strand features for each sequence


./BioPython_seqio.py [fasta_orig] [fasta_reverse_complement] 

Reads a given (yeast.fasta as fasta_orig, for example) file, then generates a new fasta file saved as 
[fasta_reverse_complement] as the reverse complement of the sequence from the original fasta, using SeqIO
Documentation for this script with SeqIO: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec53


./sliding_window_fasta.py [kmer] [fasta]

Given [kmer] and [fasta], will print out header and print tab separated kmer and GC content 

