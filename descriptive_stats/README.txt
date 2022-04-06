
Author: Rebecca Weiss

Directory contains the python script descriptive_statistics.py, which when run with the commands below will display
descriptive statistics from a given column of the file:

python3 descriptive_statistics.py [data file] [column number]

The script will read in the lines of the input file entered and store the variables in the given column number
into a list. The output will notify the user if there are lines that are a string and skip them, and then display the
following output that is calculated in the script:
the column number that the user specified, count (number of items in that column number),
ValidNum (number of items in column that are numbers (excluding strings/NaN) and will be used for calculating the
remaining statistics)), Average, Maximum, Minimum, Variance, Std Dev (standard deviation), and Median.

