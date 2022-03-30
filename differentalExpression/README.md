# Differential Expression Analysis

## Author: Rebecca Weiss

## Overview
See methodsResults files

## Module 4: Differential Expression Analysis

## Methods

Differential expression analysis is done in two steps: estimate the
relative abundance of transcripts, and apply statistical models to
identify differental expression. The relative abundance of transcripts
are estimated by evaluating the number of matches there are between NGS
reads and a genome. To do this, Salmon \[Patro et al\], DESeq2 \[Love et
al\], and tximport \[Soneson et al\] packages were used.

The first step is to estimate relative abundance is to use Salmon
\[Patro et al\] on RNA-seq reads to calculate transcript abundance
overall. A de-novo transcriptome was used to create an index with kmer
length of 25, and then Salmon was used to align all Aip samples to
AipIndex. Before moving on to the next step, we created a table in R
that mapped annotated transcripts to genes. Once this was done, tximport
\[Soneson et al\] was used to import Salmon’s output of abundance
estimates. Lastly, salmon alignments were then imported into DESeq2
\[Love et al\], which was used to run statistics to identify the whether
there were differentially expressed genes. Results are below.

## Results

``` r
library(knitr)
deAnnotated <- read.csv("deAnnotated.csv", header=T)
kable(deAnnotated)  
```

|  X | pathway      | ko        | log2FoldChange |      padj | Factor                        | pathname                                  |
| -: | :----------- | :-------- | -------------: | --------: | :---------------------------- | :---------------------------------------- |
|  1 | path:ko03015 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | mRNA surveillance pathway                 |
|  2 | path:ko04071 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Sphingolipid signaling pathway            |
|  3 | path:ko04111 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Cell cycle - yeast                        |
|  4 | path:ko04144 | ko:K17920 |       1.118793 | 0.0000011 | Menthol\_Menthol\_vs\_Control | Endocytosis                               |
|  5 | path:ko04151 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | PI3K-Akt signaling pathway                |
|  6 | path:ko04152 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | AMPK signaling pathway                    |
|  7 | path:ko04261 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Adrenergic signaling in cardiomyocytes    |
|  8 | path:ko04390 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway                   |
|  9 | path:ko04391 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway - fly             |
| 10 | path:ko04530 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Tight junction                            |
| 11 | path:ko04728 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Dopaminergic synapse                      |
| 12 | path:ko05142 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Chagas disease (American trypanosomiasis) |
| 13 | path:ko05160 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Hepatitis C                               |
| 14 | path:ko05165 | ko:K04354 |     \-2.068973 | 0.0273098 | Menthol\_Menthol\_vs\_Control | Human papillomavirus infection            |

## References

<div id="refs" class="references">

<div id="ref-Love2014">

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550.
<https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Patro2017">

Patro, Rob, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl
Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of
Transcript Expression.” *Nature Methods* 14 (4): 417–19.
<https://doi.org/10.1038/nmeth.4197>.

</div>

<div id="ref-Soneson2015">

Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. 2015.
“Differential Analyses for RNA-Seq: Transcript-Level Estimates Improve
Gene-Level Inferences.” *F1000Research* 4: 1521.
<https://f1000research.com/articles/4-1521/v2/pdf>.

</div>

</div>
