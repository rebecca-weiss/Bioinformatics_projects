#!/usr/bin/env bash
# runQuast.sh

#make output dir
outPath="quast_output/"
mkdir -p $outPath

#reference genome file
refgen=/data/METHODS/Fall/Module10/GCF_000012905.2_ASM1290v2_genomic.fna

##commands to run quast
quast.py \
--threads 4 \
--gene-finding \
-R /data/METHODS/Fall/Module10/GCF_000012905.2_ASM1290v2_genomic.fna \
-s Rhodo/scaffolds.fasta \
-o $outPath 





