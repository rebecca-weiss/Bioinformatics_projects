#!/usr/bin/env Rscript
# AnalyzeOutSpattFrom2Genomes.R

#limits final non-exponent values up to e-30 in our output
options("scipen" = 30) 

#usage statement to display script run with input file
usage <- "\nUsage: AnalyzeOutSpattFrom2Genomes.R <input file.txt>\n\n"

#tell R to get arguments
args <- commandArgs(trailingOnly = TRUE)

#ensure argument was passed to open the file
if(length(args) == 0) {
    cat(prompt=usage)
    q(save="no")
}

#create a data frame -- useful because columns can have different modes + names
spattMetrics <- read.table(args[1], sep="\t", header=FALSE)

#name columns in the data frame
names(spattMetrics) <- c("kmer1","occurrence1","expected1", "pvalue1",
                         "kmer2","occurrence2","expected2", "pvalue2") # variable names


#Add new column for probability occurence to data frame, and give all 0
spattMetrics["probablity_occurrence"] <- 0

#fill the column
spattMetrics$probablity_occurrence <- (spattMetrics$expected1 / (spattMetrics$expected1 + spattMetrics$expected2))

#Add new column for trials to data frame, and give all 0
spattMetrics["trials"] <- 0

#fill the column
spattMetrics$trials <-  spattMetrics$occurrence1 + spattMetrics$occurrence2

#create binomial upper tail column
spattMetrics["binomial_upper_tail"] <- 0

#Calculate the upper tale and fill the column
spattMetrics$binomial_upper_tail <- pbinom(spattMetrics$occurrence1-1,
                                spattMetrics$trials,
                                spattMetrics$probablity_occurrence,
                                lower.tail=FALSE)


#write output to a file I can open in Excel-- here a tab-separated text file
write.table(spattMetrics, file="out.txt", sep="\t", row.names = FALSE)