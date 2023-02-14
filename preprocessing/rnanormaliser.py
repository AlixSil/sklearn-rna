import pandas as pd
import numpy as np
from sklearn.base import TransformerMixin, BaseEstimator
from scipy.stats.mstats import gmean
from statistics import median


class DESeqNormalizer(BaseEstimator, TransformerMixin):
    """This normalizer aims to replicate a DESeq (and DESeq 2) library normalisation in the context of a regression or a classification.
    Notably, genes used in the size factor computation will be defined only during the fit part.

    The parts of DESeq(2) model usefull only for Differential expression have been removed.

    This normalisation is supposed to be used on raw read numbers.

    For more information, please refer to the original article

    Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8. PMID: 25516281; PMCID: PMC4302049.
    """

    def fit(self, X, **fit_params):
        """Fit method will, in our case decide which are the genes that are used in the computation of size factor (by removing genes that have a zero expression in at least one sample),
        note that if a sample wasn't fitted on, he will not be taken into account in the computation of the geometric mean values."""
        if type(X) == pd.DataFrame:
            X_copy = X.copy().to_numpy()
        else:
            X_copy = X.copy()

        #check for negative values
        if (X_copy < 0).any() :
            raise ValueError("DESeqNormalizer cannot be fitted on negative values")

        self.non_zero_genes_indexes = np.where(X.all(axis=0))[0]

        #Check we have at least one gene with non zero values
        if len(self.non_zero_genes_indexes) == 0 :
            raise ValueError("DESeqNormalizer need at least one gene with no zero vvalues to be fitted")

        self.non_zero_genes_gmean_values = np.apply_along_axis(
            gmean, 0, X_copy[:, self.non_zero_genes_indexes]
        )

        return (self.non_zero_genes_indexes, self.non_zero_genes_gmean_values)

    def transform(self, X, **fit_params):
        """The transform method will compute the size factors of each sample, if the sample was has a gene with expression zero who was selected in the fit, the gene will be ignored only for this sample"""

        if type(X) == pd.DataFrame:
            X_copy = X.copy().to_numpy()
        else:
            X_copy = X.copy()
        print(X_copy)
        

        size_factor = []
        for i in range(X_copy.shape[0]):
            print("---------")
            print(i)
            sample_genes = X_copy[i, self.non_zero_genes_indexes]
            print("sample genes : {}".format(sample_genes))
            temp = sample_genes / self.non_zero_genes_gmean_values
            print("temp : {}".format(temp))
            temp = [n for n in temp if n != 0]
            size_factor = median(temp)
            print("size factor : {}".format(size_factor))

            X_copy = X_copy.astype(np.float32)
            X_copy[i] = X_copy[i] / size_factor

        print("final X_copy {}".format(X_copy))
        return(X_copy)