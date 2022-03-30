#!/usr/bin/env bash
# runSpades.sh

#absolute path for spades on Defiance
spades=/usr/local/programs/SPAdes-3.14.1-Linux/bin/spades.py
outPath="Rhodo/"

#make Rhodo dir
mkdir -p $outPath

#Commands for running spades with 4 threads, output Rhodo directory, to assemble the Rhodobacter genome using just the quality-trimmed reads in Paired.  
$spades \
--threads 4 \
--pe1-1 Paired/SRR522244_1.fastq --pe1-2 Paired/SRR522244_2.fastq \
-o $outPath \

1>runSpades.log 2>runSpades.err &