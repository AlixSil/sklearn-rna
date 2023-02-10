# sklearn-rna

The aim of this project is to expand on sklearn classes to allow for better handling of RNA-seq data.

## Planned features
- Nomalisation
	- Deseq2 normalisation
	- edgeR normalisation
	- TPM (?)
	- FPKM (?)
	- Comparison of normalisation efficiency for classification tasks

Notes, both TPM and FPKM are still used in pracice, but neither are made to compare samples between them, hence my current hesitancy to include them.

## Notes
Some modifications for fit and train will have to be made (in order to avoid any data leakage).

