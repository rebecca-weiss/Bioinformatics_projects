#!/usr/bin/env python3
# descriptive_statistics.py
"""User inputs a file name and column number and the script returns descriptive statistics"""

import sys
import math

# notes about the arguments/inputs needed to run this script; exit if two variables are not provided
if (len(sys.argv) - 1) < 2:
    sys.exit("\nTwo arguments are required: 1. an input file and 2. a column number to parse\n")
FILENAME = sys.argv[1]
COLUMN_TO_PARSE = int(sys.argv[2])

# establish empty list of NUMBERS to run statistics on, and counter of ALL lines including nan
NUMBERS = []
COUNT = 0

# opening the file and iterating over the lines; count lines and append valid numbers to list:
with open(FILENAME, 'r') as infile:
    for line in infile:
        COUNT += 1
        try:
            num = float(line.split("\t")[COLUMN_TO_PARSE])
            if not math.isnan(num):
                NUMBERS.append(float(num))
        except ValueError:
            print("\nSkipping line number {} : could not convert string to float: '{}'".format(
                COUNT, line.split("\t")[COLUMN_TO_PARSE]))
            continue
        except IndexError:
            print("\nExiting: There is no valid 'list index' in column {} in line {} "
                  "in file: {}\n".format(COLUMN_TO_PARSE, COUNT, FILENAME))
            sys.exit(1)


# Count the number of valid numbers, and if there are not any, exit the program with message
if len(NUMBERS) < 1:
    print("\nError: There were no valid number(s) in column {} in file: {}\n".format(
        COLUMN_TO_PARSE, FILENAME))
    sys.exit(1)
else:
    VALIDNUM = len(NUMBERS)


# Calculate the average of the valid numbers in the NUMBERS list
AVG = sum(NUMBERS)/len(NUMBERS)

def calc_variance(vals):
    """Calculates the difference of each item of NUMBERS from AVG, then calculates the variance"""
    stdevs = [(i - AVG) ** 2 for i in vals]

    # Calculate the variance, correcting for variables with only one valid number:
    if VALIDNUM == 1:
        variance = 0.0
    else:
        variance = sum(stdevs) / (VALIDNUM - 1)
    return variance

def calc_stdev(vals):
    """Calculate standard deviation using the output of calc_variance()"""
    stdev = calc_variance(vals) ** 0.5
    return stdev

def calc_median(vals):
    """Calculates the median"""
    value_sorted = sorted(vals)
    if VALIDNUM % 2 == 0:
        lower = value_sorted[VALIDNUM // 2]
        upper = value_sorted[(VALIDNUM // 2) - 1]
        median = (lower + upper) / 2
    else:
        median = value_sorted[VALIDNUM // 2]
    return median

# Output to console:
print("\n" + " "*4 + "{:<8} {} \n\n".format("Column:", COLUMN_TO_PARSE))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Count", "=", COUNT))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("ValidNum", "=", VALIDNUM))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Average", "=", AVG))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Maximum", "=", max(NUMBERS)))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Minimum", "=", min(NUMBERS)))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Variance", "=", calc_variance(NUMBERS)))
print(" "*8 + "{:<8} {:<4} {:>7.3f}".format("Std Dev", "=", calc_stdev(NUMBERS)))
print(" "*8 + "{:<8} {:<4} {:>7.3f}\n".format("Median", "=", calc_median(NUMBERS)))
