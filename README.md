# sklearn-rna

The aim of this project is to expand on sklearn classes to allow for better handling of RNA-seq data.

## Current features
- Normalisation
	- Deseq(2) normalisation
	- TTM normalisation

## Planned features
- Normalisation
	- Propose a multithreading of TTM transform

- Benchmarking
	- benchmarking of each normalisation method for the purpose of patient classification.

## Notes
Some modifications have been made to the some normalization (namely DESeq(2)) in order to make them compatible with a fit/transform paradigm where not every samples are available at the same time. Please read the code or the future documentation in detail if needed.



# References
Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8. PMID: 25516281; PMCID: PMC4302049.

Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 2010;11(3):R25. doi: 10.1186/gb-2010-11-3-r25. Epub 2010 Mar 2. PMID: 20196867; PMCID: PMC2864565.

