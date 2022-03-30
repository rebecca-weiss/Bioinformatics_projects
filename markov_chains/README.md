# BINF6309 Module 3: Markov Chains and HMMs in R

## Authors: Dr. Kaluziak (AnalyzeOutSpattfrom2Genomes.R) and Benjamin Tovar and Avril Coghlan, and Chesley Leslin (markov_chains.R)

## Overview
./AnalyzeOutSpattFrom2Genomes.R Pyrococcus_horikoshii.fasta.8.m1.allwords_Pyrococcus_abyssi.fasta.8.m1.allwords.txt
The output of two sspatt runs  from combining all kmers of length 8 using a first order Markov model from Pyrococcus_horikoshii.fasta and Pyrococcus_abyssi.fasta genomes; are pasted into one text file, that is the observed and expected values of the of sspatt runs. This script takes that text file and compares whether the occurence of a finding in one genome is what we would expect based on a binomial distribution, and the likelihood of the occurences based on the combined genomes.
The output of this file is stored in output.txt.



./markov_chains.R
"This program show you how to generate a DNA sequence using a multinomial model of DNA sequence evolution, and then moves onto a Markov model of DNA sequence evolution, plotting data as we go along" (markov_chains.R file)