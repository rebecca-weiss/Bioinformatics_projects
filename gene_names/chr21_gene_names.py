#!/usr/bin/env python3
# chr21_gene_names.py
"""Exercise 1: User inputs a gene symbol and then the description for that gene is printed
based on data from the chr21_genes.txt file, which is stored in a dictionary"""


import argparse
import sys
from assignment4 import my_io


def main():
    """Run the primary program"""
    args = get_cli_args()
    infile1 = args.infile

    # Filehandle to open
    fh_in = my_io.get_fh(infile1, "r")

    # Convert filehandle to dictionary
    gene_dict = file_2_dict(fh_in)

    # Collect user input, get description (if available), and print to STDOUT
    gene_input = input('\nEnter gene name of interest. Type quit to exit: ')
    while not (gene_input.upper() == "QUIT" or gene_input.upper() == "EXIT"):
        if gene_input.upper() in gene_dict:
            sys.stdout.write("\n{} found! Here is the description: \n{}\n".format(
                gene_input, gene_dict[gene_input.upper()]))
            gene_input = input('\nEnter gene name of interest. Type quit to exit: ')
        elif gene_input.upper() not in [symbol.upper() == gene_input.upper() for symbol in gene_dict.keys()]:
            sys.stdout.write('Not a valid gene name\n')
            gene_input = input('\nEnter gene name of interest. Type quit to exit: ')
    if gene_input.upper() == "QUIT" or gene_input.upper() == "EXIT":
        sys.stdout.write('Thanks for querying the data.\n')


def file_2_dict(fh_in):
    """Converts gene file into a dictionary where gene_symbol = key, gene_desc = value"""

    # Create empty dictionary for gene symbol/description
    gene_dict = {}

    with fh_in as file:
        next(file) #skip first line
        for line in file:
            gene_desc, gene_symbol = line.split('\t')[0:2]
            gene_dict.update({gene_desc: gene_symbol})

    return gene_dict


def get_cli_args():
    """Get the command line arguments (infile); takes no arguments"""
    parser = argparse.ArgumentParser(description='Open chr21_genes.txt, and ask user for a gene name')

    # Add the arguments
    parser.add_argument('-i', '--infile',
                        dest='infile', help='Path to the file to open',
                        type=str, default='chr21_genes.txt')

    # Execute the parse_args()
    return parser.parse_args()


if __name__ == "__main__":
    main()
