#!/usr/bin/env python3
# gene_information_query.py
"""Given a gene and species, query for a specific tissue expression"""

import argparse
import sys
import re
from assignment5 import my_io
from assignment5 import config


def main():
    """Run the primary program"""
    args = get_cli_args()
    temp_host = args.HOST
    host = modify_host_name(temp_host)
    temp_gene = args.GENE
    if temp_gene != None:
        gene = temp_gene
    else:
        gene = "TGM1"

    file = "/".join((config.get_unigene_directory(), host, gene + "." + config.get_unigene_extension()))

    # check for the existence of file
    if my_io.is_valid_gene_file_name(file):
        # using f-strings
        print(f"\nFound Gene {gene} for {host}")
    else:
        print("Not found")
        print(f"Gene {gene} does not exist for {host}. exiting now...", file=sys.stderr)
        sys.exit()

    tissue_list = get_gene_data(file)
    print_output(host, gene, tissue_list)


def modify_host_name(temp_host):
    """Takes the host name and checks for the available conversions from common
    to scientific names and return the scientific name."""
    keyword_dict = config.get_host_keywords()
    if temp_host is None:
        host = "homo sapiens"
    else:
        host = temp_host.lower()

    # return keyword_dict.keys(), host
    if host in keyword_dict.keys():
        return keyword_dict.get(host)
    if host not in keyword_dict.keys():
        _print_host_directories()
        sys.exit()


def _print_host_directories():
    """If hostname does not exist in dict, will notify user what directories exist"""
    sys.stdout.write('\nEither the Host Name you are searching for is not in the database\n\n'
                     'or If you are trying to use the scientific name please put the name in double quotes:\n\n'
                     '"Scientific name"')
    __print_scientific_name()
    __print_common_name()


def __print_common_name():
    """Will print the common name when the file name is unable to open"""
    keyword_dict = config.get_host_keywords()
    common_name = keyword_dict.keys()
    sys.stdout.write("\n\nHere is a (non-case sensitive) list of available Hosts by common name\n\n")
    list_common = tuple(enumerate(sorted(set(common_name)), 1))
    for index, name in list_common:
        sys.stdout.write("{:>2}. {}\n".format(index, name.capitalize()))


def __print_scientific_name():
    """Will print the scientific name when the file name is unable to open"""
    keyword_dict = config.get_host_keywords()
    scientific_name = keyword_dict.values()
    sys.stdout.write("\n\nHere is a (non-case sensitive) list of available Hosts by scientific name\n\n")
    list_scientific = tuple(enumerate(set(scientific_name), 1))
    for index, name in list_scientific:
        sys.stdout.write("{:>2}. {}\n".format(index, name.capitalize()))


def get_gene_data(gene_file):
    """opens the file for the host and gene, extracts the list of tissues in which this gene
    is expressed and returns a sorted list of the tissues"""
    fh_in = my_io.get_fh(gene_file, "r")
    for line in fh_in:
        match = re.search(r'^EXPRESS\s+(\D+)', line)
        if match:
            tissue_strig = match.group(1)
            temp_tissue_list = list(tissue_strig.split(sep='|'))
            tissue_list = sorted([tissue.strip() for tissue in temp_tissue_list])
            return tissue_list


def print_output(host, gene, tissue_list):
    """Prints tissue expression + index"""
    tissues_number = len(tissue_list)
    print(f"In {host}, There are {tissues_number} tissues that {gene} is expressed in:\n")
    for index, tissue in sorted(enumerate(tissue_list, 1)):
        print(f"{index:>2}. {tissue.capitalize()} ")
    print('\n')


def get_cli_args():
    """Get the command line arguments; takes no arguments"""
    parser = argparse.ArgumentParser(description='Give the Host and Gene name')

    # Add the arguments
    parser.add_argument('-host', dest='HOST', help='Name of Host', required=False)
    parser.add_argument('-gene', dest='GENE', help='Name of Gene', required=False)

    # Execute the parse_args()
    return parser.parse_args()


if __name__ == "__main__":
    main()
