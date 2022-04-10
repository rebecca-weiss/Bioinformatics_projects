# Genome Assembly

## Rebecca Weiss

## Methods
The genomic reads from Rhodobacter spheroides were pulled from the Sequence Read Archive (SRA) toolkit using fastq-dump [1]. The fastq files were then quality-trimmed using PE trimmomatic [2]. The quality-trimmed paired reads from the Rhodobacter genome were assembled using SPAdes [3], with output in the Rhodo/directory containing contigs, scaffolds, etc. QUAST [4] was then used to run summary statistics on the scaffold assembly from the previous step, using Rhodobacter genome as a reference [5]. 

## Conclusions from Analysis
All output from the QUAST analysis [4] can be found in module10-rebecca-weiss/GenomeAssembly/quast_output/. The output in report.txt shows that the N50 is higher in the scaffolds than it is in the scaffolds_broken version. This makes sense because we expect to have the same number or fewer misassemblies than the scaffolded version [4], and can see this pattern visually in /quast_output/basic_stats/Nx_plot.pdf. The report.txt shows that both scaffolds and scaffolds_broken have the same GC% (68.80), which is nearly identical, but only slight lower than, the reference genome GC% (68.79); we can also see this visually in GC_content_plot.pdf, where the reference and scaffold/scaffold_broken are all nearly on top of one another. The total length of the scaffolds assembly, 453175, is smaller than the reference length of 460297, which can be gathered from report.txt. The QUAST output in report.txt shows that the number of predicted genes found by GeneMarkS is 4471, which is slightly higher than the 4464 that are listed on its NCBI page [5].


## Future Directions
Looking at the GFF file we calculated with QUAST (GenomeAssembly/quast_output/predicted_genes/scaffolds_genemark_genes.gff), future directions for work arise:
-Associating functions to these proteins, using GO terms. From the NCBI page [5], we know that Rhodobacter sphaeroides is used in metabolic processes, specifically in photosynthesis and chemotaxis, and can be used in aerobic and anaerobic conditions; thus, terms related to the these processes, as seen in module 9 can be used, such as ko00195 Photosynthesis, or ko02030	Bacterial chemotaxis, as examples 
-Parsing of the GFF file to obtain the FASTA file, via Python parsing



## References
1. Toolkit Documentation: Software: Sequence Read Archive: NCBI/NLM/NIH. (n.d.). Retrieved November 23, 2020, from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

2. Bolger, Anthony M, Lohse, Marc, and Usadel, Bjoern. "Trimmomatic: A Flexible Trimmer for Illumina Sequence Data." Bioinformatics (Oxford, England) 30.15 (2014): 2114-120. Web.

3. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

4. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

5. ASM1290v2—Genome—Assembly—NCBI. (n.d.). Retrieved November 23, 2020, from https://www.ncbi.nlm.nih.gov/assembly/GCF_000012905.2

