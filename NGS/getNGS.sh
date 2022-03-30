#!/usr/bin/env bash
# getNGS.sh

# Retrieve the Rhodobacter spheroides NGS reads.
fastq-dump --split-files SRR522244 1>getNGS.log 2>getNGS.err