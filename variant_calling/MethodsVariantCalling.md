## Methods: Variant Calling

Variant calling is used to identify genetic differences between
individuals of the same species (germline), or of tissues from the same
individual (somatic). The overall workflow for conducting variant
calling is to 1. conduct alignment with a reference genome, and 2.
analyze the alignment using a variant calling program such as GATK
(McKenna et. al) or DeepVariant (Poplin et. al). In this module, we run
DeepVariant protocol outlined in the Poplin et. al 2018 paper.

After the initial steps of getting the genome (GRCh38), indexing it,
trimming using Trimmomatic (Bolger et. al), and indexing the reads into
an alignment file, we are ready to begin the alignment. The next step of
this variant calling pipeline, the alignment, was done using Burrows
Wheeler aligner (BWA) (Li et. al). BWA consists of 3 algorithmic
methods: backtrack, SW, and MEM. In this module, we used BWA-MEM to
align the reads which we can then sort and create the bam file.
Following the alignment, the reads can then be indexed and DeepVariant
(Poplin et. al) is run to create the vcf file.

## References

<div id="refs" class="references">

<div id="ref-Bolger2014">

Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics* 30 (15):
2114–20. <https://doi.org/10.1093/bioinformatics/btu170>.

</div>

<div id="ref-Li2009">

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read
Alignment with Burrows-Wheeler Transform.” *Bioinformatics (Oxford,
England)* 25 (14): 1754–60.
<https://doi.org/10.1093/bioinformatics/btp324>.

</div>

<div id="ref-McKenna2010">

McKenna, Aaron, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian
Cibulskis, Andrew Kernytsky, Kiran Garimella, et al. 2010. “The Genome
Analysis Toolkit: A MapReduce Framework for Analyzing Next-Generation
DNA Sequencing Data.” *Genome Research* 20 (9): 1297–1303.
<https://doi.org/10.1101/gr.107524.110>.

</div>

<div id="ref-Poplin2018">

Poplin, Ryan, Pi-Chuan Chang, David Alexander, Scott Schwartz, Thomas
Colthurst, Alexander Ku, Dan Newburger, et al. 2018. “A Universal SNP
and Small-Indel Variant Caller Using Deep Neural Networks.” *Nature
Biotechnology* 36 (10): 983–87. <https://doi.org/10.1038/nbt.4235>.

</div>

</div>
