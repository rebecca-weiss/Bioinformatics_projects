#!/usr/bin/env bash
# alignAll.sh

outDir='quant/'
sample='/scratch/SampleDataFiles/Paired'

function alignAll {
    for fastq in $sample/Aip??.R1.paired.fastq; do
      #get the correct sample name from by iterating over the file list
      sampleName=`echo $fastq | rev | cut -d. -f4 | cut -d/ -f1 | rev`

    #use salmon to perform abundance estimation on each one
      salmon quant -l IU is \
        -1 $sample/$sampleName.R1.paired.fastq \
        -2 $sample/$sampleName.R2.paired.fastq \
        -i AipIndex \
        --validateMappings \
        -o $outDir$sampleName
    done
}

alignAll 1>align.log 2>align.err &

