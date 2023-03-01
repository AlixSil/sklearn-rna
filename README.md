# sklearn-rna

Most RNA normalization method are aiming to correct for differences in library size (i.e. one sample has in total more read than others) and differences in library composition (i.e. gene G is highly transcribd in 20% of the samples, but not transcribed in the others). To that end, multiple samples are treated together.

However, this creates both replicability problems (Patil et. al.) as new data cannot be normalised into the same conditions, and data leakage problem as train and test are often normalised together, and hence part of test information is used to format the train dataset.

In this repository, I aim to extend sklearn with some classical RNA normalization method, modified for a data science use. More specifically I enforce the following demands.
- the fit method will learn what is needed for the normalization from the training dataset.
- the transform method will have the same results for samples that are processed separately or together.

I also aim to benchmark them in the case of survival study prediction, against themselves and against classical TPM and FPKM. I'm expecting the classical TPM and FPKM to be superior because of data leakage.


## Normalization methods
### Deseq(2) normalization
Normalization inspired by the scale factoring present in DESeq2 (Love MI et. al.)

This normalization had to have been severely changed from its origin to fit with the demands of this study.
In practice almost all the information on the size factor comes from the fit test.

### TMM normalization
Trimmed Mean of M values normalization method (Robinson MD et. al.).
This method applies itself very easily to our framework, as it normalises each sample against a reference sample. We hence select a reference sample from our fitted test and go from there.

### Rank based normalization
The value of gene g in sample s becomes the rank of gene g in sample s.
Each sample is treated independantly, the only fit is to guarantee that all samples have the same number of genes.


## Benchmarking
WIP


# References
Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8. PMID: 25516281; PMCID: PMC4302049.

Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 2010;11(3):R25. doi: 10.1186/gb-2010-11-3-r25. Epub 2010 Mar 2. PMID: 20196867; PMCID: PMC2864565.

Prasad Patil, Pierre-Olivier Bachant-Winner, Benjamin Haibe-Kains, Jeffrey T. Leek, Test set bias affects reproducibility of gene signatures, Bioinformatics, Volume 31, Issue 14, July 2015, Pages 2318–2323, https://doi.org/10.1093/bioinformatics/btv157