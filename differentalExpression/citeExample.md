## Overview

The two main steps in performing differential expression analysis are to
estimate the relative abundance of transcripts, and to apply statistical
models to test for differential expression between treatment groups.
Estimating relative abundance is basically determining how many NGS
reads match a given gene within a genome. In this module you’ll use
Salmon
(<span class="citeproc-not-found" data-reference-id="Patro">**???**</span>)
to estimate relative abundance, tximport
(<span class="citeproc-not-found" data-reference-id="Soneson">**???**</span>)
to import the Salmon abundance estimates, and DESeq2
(<span class="citeproc-not-found" data-reference-id="Love">**???**</span>)
to perform statistical tests to identify differentially expressed genes.

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
