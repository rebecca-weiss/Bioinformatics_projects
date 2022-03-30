BINF 6200: Assignment 4

Author: Rebecca Weiss

chr21_gene_names.py
Overview:
This program takes one argument as infile (default: chr21_genes.txt), and then the user inputs a gene symbol --
program will then print the description for that gene based on data from the infile.
Example: python3 chr21_gene_names.py -i chr21_genes.txt


categories.py
Overview:
This program takes two arguments/infiles, and counts genes occurrence in each category (1.1, 1.2, 2.1 etc.)
from the chr21_genes.txt file to match with the description file (chr21_genes_categories.txt).
The output file is then stored in OUTPUT/categories.txt
Example: python3 categories.py -i1 chr21_genes.txt -i2 chr21_genes_categories.txt


intersection.py
Overview:
This program takes two arguments/infiles and looks for all gene symbols that appear in both of them.
The gene symbols are printed to a file in alphabetical order at OUTPUT/intersection_output.txt, and the number of unique
gene names in each file, along with the number of common symbols, are printed to STDOUT.
Example: python3 intersection.py -i1 chr21_genes.txt -i2 HUGO_genes.txt

assignment4/my_io.py
Overview: module that takes in filehandle that serves for input and output (writing to output, and reading in)
file handles for all programs.
Example of implementation in three programs above: from assignment 4 import my_io.
