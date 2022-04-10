## 
# Module 11: Cancer Genomics

## Author: Rebecca Weiss

Overview:
The assignment was completed by using maftools vignette tutorial (https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#1_Introduction) in sections 6, 7, and 9, as detailed in the cancerGenomics.Rmd file. The cancerGenromics.sh script was used to render the markdown file to html format. 

---
title: "R Notebook"
author: "Rebecca Weiss"
date: "4/15/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(maftools)
```



# 6: Reading and summarizing maf files
The only required input is the MAF file, that can be in a zip/gz format. Other optional arguments are clinical data associated in MAF format, and the copy number data (GISTIC output or a customized data table). In this assignment, read.table is used since the metadata has already been obtained. read.maf is used to read and summarize the MAF files and stores that summary in MAF format. Below, we have read and summarized the data into the variable "laml", and can access the basic summary by inputting laml into the command line. 

Reading input files into R for parts 6, 7, and 9. 
```{r}
# TCGA laml with clinical data
clinical.data = read.table("tcga_laml_annot.tsv", sep="\t", header = TRUE) #read in metadata

laml.maf = read.maf(maf = "tcga_laml.maf.gz", clinicalData = clinical.data) #create a MAF object

# BRCA with no clinical data
brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)

#Primary and relapse APL MAFs with no clinical data
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
relapse.apl = read.maf(maf = relapse.apl)
```


Showing summary data from TCGA laml MAF file.
```{r}
#Typing laml shows basic summary of MAF file.
laml.maf

#Shows sample summry.
getSampleSummary(laml.maf)

#Shows gene summary.
getGeneSummary(laml.maf)

#shows clinical data associated with samples
getClinicalData(laml.maf)

#Shows all fields in MAF
getFields(laml.maf)

#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml.maf, basename = 'laml')
```


# 7: Visualization
'plotmafSummary' is useful to plot the summary of the maf file, which displays number of variants in each sample as a stacked barplot and variant types as a boxplot summarized by Variant_Classification. Average mean or median of variants within the cohort can be displayed on the stacked barplot. Oncoplots, aka waterfall plots, are an even better way to visualize maf files and can be customized in a variety of ways.

'titv' classifies SNPs into Transitions and Transversions, and returns a list of summarized tables in various ways such as a boxplot showing overall distribution of six different conversions, and as a stacked barplot showing fraction of conversions in each sample. Lollipop plots are a way to visualize amino acid changes. Rainfall plots are a way to visualize localized hyper-mutations that are commonly found in cancer genomes. 

'tcgaCompare'function uses mutation load from TCGA MC3 for comparing muttaion burden against 33 TCGA cohorts. 'plotVaf' plots Variant Allele Frequencies as a boxplot which quickly helps to estimate clonal status of top mutated genes.


```{r}
# Plotting MAF summary
plotmafSummary(maf = laml.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
```


```{r}
# oncoplot for top ten mutated genes. Note that Multi hit indicates genes are mutated more than once in a single sample
oncoplot(maf = laml.maf, top = 10)
```


```{r}
laml.titv = titv(maf = laml.maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
```


```{r}
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = laml.maf, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
```

```{r}
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
#NOTE: if labelPos is = 'all', all points will be highlighted
lollipopPlot(maf = laml.maf, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
```


```{r}
# Rainfall plot of brca MAF
# NOTE: If detectChangePoints = TRUE, rainfall plot also highlights regions where potential changes in inter-event distances are located
rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)
```


```{r}
# Compares mutation loads against TCGA cohorts
laml.mutload = tcgaCompare(maf = laml.maf, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)
```

```{r}
# Plot VAF
# NOTE: clonal genes usually have mean allele frequency around ~50% assuming pure sample
plotVaf(maf = laml.maf, vafCol = 'i_TumorVAF_WU')
```

# 9: Analysis (note: no 9.9)
There are several functions that provide different analyses within maftools. 'somaticInteractions' function performs pair-wise Fisher’s Exact test to detect such significant pair of genes, which detects mutually exclusive or co-occurring set of genes. 'oncodrive' identifies cancer genes (driver) from a given MAF, to identify cancer genes. 'pfamDomains' is a useful tool to add pfam domain information to AA changes.  Survival analysis is a crucial aspect of the analyses and can be done via the 'mafSurvive' function, which draws a kaplan meier curve by grouping samples based on mutation status. 

There are numerous ways to compare two cohorts and visualize the results. 'mafCompare' performs fisher test on all genes between two cohorts to detect differentially mutated genes -- in this case, primary and relapse APL.  The results of 'mafCompare' provides info on highly mutated in Relapse APL compared to Primary APL, which is then visualized as a forestplot. Another way to visualize the results is 
with 'coOncoplot', which takes two maf objects and plots them side by side for better comparison. Co-bar plots can be similarly generated using 'coBarplot'. 'lollipopPlot2' shows genewise differences.

Enrichment analyses are a useful too to identify enriched mutations for every category within a clinical feature. The 'clinicalEnrichment' function takes any clinical feature associated with the samples and performs various groupwise and pairwise comparisons. 'drugInteractions' function checks for gene druggability information and drug–gene interactions compiled from Drug Gene Interaction database, and the results can then be plotted. Enrichment of known Oncogenic Signaling Pathways in TCGA cohorts can be determined using the 'OncogenicPathways' function. Lastly, a number of tools demonstrated below, will demonstrate the use of mutational signatures / specific nucleotide substitutions that illustrate the progression of cancers. These analyses require a BSGenome object, and then take the surrounding bases of a mutated base to generate a mutation matrix.

```{r}
# Exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml.maf, top = 25, pvalue = c(0.05, 0.1))
```


```{r}
# Detecting cancer driver genes based on positional clustering
laml.sig = oncodrive(maf = laml.maf, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

head(laml.sig)

# Plot the results
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
```


```{r}
#  Adding and summarizing pfam domains
laml.pfam = pfamDomains(maf = laml.maf, AACol = 'Protein_Change', top = 10)

#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]

#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]
```

```{r}
# Mutation in any given genes: survival analysis based on grouping of DNMT3A mutation status
mafSurvival(maf = laml.maf, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)

# Predict genesets associated with survival
#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = laml.maf, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)

# Generate KM curve for group of 2 genes with associated with poor survival (p < 0.05) 
mafSurvGroup(maf = laml.maf, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")
```

```{r}
# Comparing two cohorts: primary and relapse APLs 
pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
print(pt.vs.rt)
```

```{r}
## Forest plot 
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, geneFontSize = 0.8)
```

```{r}
# Co-onco plot
genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)
```

```{r}
# Co-bar plots
coBarplot(m1 = primary.apl, m2 = relapse.apl, m1Name = "Primary", m2Name = "Relapse")
```

```{r}
# Lollipop plot-2
lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")
```

```{r}
# Clinical enrichment to ID mutations with FAB classification
fab.ce = clinicalEnrichment(maf = laml.maf, clinicalFeature = 'FAB_classification')

# Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]

# Plot the results
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
```

```{r}
# Drug-gene interactions
dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)

# Printing selected columns.
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
```


```{r}
# Oncogenic signaling pathways
OncogenicPathways(maf = laml.maf)

# Plot the complete pathway
PlotOncogenicPathways(maf = laml.maf, pathways = "RTK-RAS")
```



```{r}
#Mutational signatures -- requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)

# Matrix generation
laml.tnm = trinucleotideMatrix(maf = laml.maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")


#  Differences between APOBEC enriched and non-enriched samples
plotApobecDiff(tnm = laml.tnm, maf = laml.maf, pVal = 0.2)

# Signature analysis
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6, pConstant = 0.01)

#Draw elbow plot to visualize and decide optimal number of signatures from above results
plotCophenetic(res = laml.sign)
laml.sig = extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.01)

#Compare against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

#Compare against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

# Generate plot to compare similarities of detected signatures against validated signatures
library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

# Plot signatures
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")

# Visualize with barplot3d
library("barplot3d")
#Visualize first signature
sig1 = laml.sig$signatures[,1]
barplot3d::legoplot3d(contextdata = sig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)
```
