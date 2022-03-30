#!/usr/bin/env python3
# categories.py
"""Exercise 2: counts how many genes are in each category (1.1, 1.2, 2.1 etc.) based on data from chr21_genes.txt.
The results will print so that categories are arranged in ascending order to an output file:
OUTPUT/categories.txt"""

import argparse
import collections
# import re
import os.path
from assignment4 import my_io


def main():
    """Run the primary program"""
    args = get_cli_args()
    infile1 = args.INFILE1
    infile2 = args.INFILE2
    combine_counts_cats_output(count_gene_cat(infile1), category_2_dict(infile2))


def count_gene_cat(infile1):
    """Determine gene category from a file and its number of occurrences"""
    list_cat = []
    with my_io.get_fh(infile1, "r") as fh_in1:
        for line in fh_in1:
            line = line.rstrip()
            category = line.split('\t')[2] if len(line.split('\t')) == 3 else "Other"
            if re.match(r'^\d', str(category)):
                list_cat.append(category)
        count_dict = collections.Counter(list_cat)
    return count_dict


def category_2_dict(infile2):
    """Converts gene categories text with file to dictionary {category; description}"""
    category_dict = {}
    with my_io.get_fh(infile2, "r") as fh_in2:
        for line in fh_in2:
            line = line.rstrip()
            cat_id, cat_desc = line.split('\t')[0:2]
            category_dict.update({cat_id: cat_desc})
    return category_dict


def combine_counts_cats_output(count_dict, category_dict):
    """Create the output file, sort the number of category occurrences dict, and join/output to file"""
    # Define output path (already made, per directions of this exercise)
    out_path = os.path.join("./OUTPUT", "categories.txt")

    with my_io.get_fh(out_path, "w") as fh_out:
        fh_out.write("\t".join(("Category", "Occurrence", "Definition", "\n")))
        for k in sorted(count_dict.keys()):
            if k in category_dict.keys():
                fh_out.write("\t".join((str(k), str(count_dict[k]), str(category_dict[k]), "\n")))


def get_cli_args():
    """Get the command line arguments (infile); takes no arguments"""
    parser = argparse.ArgumentParser(description='Combine on gene name and count the category occurrence')
    parser.add_argument('-i1', '--infile1', dest='INFILE1',
                        help='Path to the gene description file to open', default='chr21_genes.txt')
    parser.add_argument('-i2', '--infile2', dest='INFILE2',
                        help='Path to the gene category to open', default='chr21_genes_categories.txt')
    return parser.parse_args()


if __name__ == '__main__':
    main()
