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

        # check for negative values
        if (X_copy < 0).any():
            raise ValueError("DESeqNormalizer cannot be fitted on negative values")

        # We store the number of gene to check during transform if the appropriate number of genes are passed
        self.number_of_genes = X_copy.shape[1]

        self.non_zero_genes_indexes = np.where(X.all(axis=0))[0]

        # Check we have at least one gene with non zero values
        if len(self.non_zero_genes_indexes) == 0:
            raise ValueError(
                "DESeqNormalizer need at least one gene with no zero vvalues to be fitted"
            )

        self.non_zero_genes_gmean_values = np.apply_along_axis(
            gmean, 0, X_copy[:, self.non_zero_genes_indexes]
        )

        return self

    def transform(self, X, **transform_params):
        """The transform method will compute the size factors of each sample, if the sample was has a gene with expression zero who was selected in the fit, the gene will be ignored only for this sample"""

        if type(X) == pd.DataFrame:
            X_copy = X.copy().to_numpy()
        else:
            X_copy = X.copy()
        print(X_copy)

        if X_copy.shape[1] != self.number_of_genes:
            raise ValueError(
                "Trying to transform a different number of gene that this transformer has been fitted on"
            )

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
        return X_copy


class TMM(BaseEstimator, TransformerMixin):
    """
    Trimmed Mean of M values normalisation method,
    The reference sample is picked out only during the fitting.
    This normalisation makes the assumption that most genes are not differentially expressed between classes.
    This normalisation is planned for raw reads.
    See the original article for more informations
    """

    def fit(self, X, library_sizes=None, reference_line=None, **fit_params):
        """TMM method uses a reference sample that will be picked among the fitted X, you can explicitly pick it out by submitting the index in reference_column.
        otherwise, the reference picked will be the sample whose upper quartile is closest to the mean upper quartile is used"""

        if type(X) == pd.DataFrame:
            X_copy = X.copy().to_numpy()
        else:
            X_copy = X.copy()

        if library_sizes is None:
            library_sizes = X.sum(axis=1)

        if reference_line is None:
            third_quartiles = np.quantile(X_copy / library_sizes[:, None], 0.75, axis=1)
            average_third_quartiles = third_quartiles.mean()
            reference_line = (
                np.abs(third_quartiles - average_third_quartiles)
            ).argmin()
            self.scaled_reference_sample = (
                X[reference_line] / library_sizes[reference_line]
            )

        else:
            if reference_line not in list(range(X_copy.shape[0])):
                raise ValueError("reference_line variable is out of bounds")
            else:
                self.scaled_reference_sample = (
                    X[reference_line] / library_sizes[reference_line]
                )
                
        return(self)

    def transform(self, X, library_sizes=None, **transform_params):
        if type(X) == pd.DataFrame:
            X_copy = X.copy().to_numpy()
        else:
            X_copy = X.copy()

        if library_sizes is None:
            library_sizes = X.sum(axis=1)


        reference_non_zero = np.nonzero(self.scaled_reference_sample)[0]


        for s in range(X_copy.shape[0]):

            # loop on sample
            current_sample = X_copy[s]
            scaled_current_sample = current_sample / library_sizes[s]

            current_sample_non_zero_index = np.nonzero(scaled_current_sample)[0]

            non_zero_indexes = np.intersect1d(reference_non_zero, current_sample_non_zero_index)


            M = np.log2(self.scaled_reference_sample[non_zero_indexes] / scaled_current_sample[non_zero_indexes])
            A = (
                np.log2(self.scaled_reference_sample[non_zero_indexes] * scaled_current_sample[non_zero_indexes])
            ) / 2

            m_lower, m_higher = np.percentile(M, [30, 70])
            a_lower, a_higher = np.percentile(A, [5, 95])


            references_genes_index = np.where(
                (M >= m_lower) & (M <= m_higher) & (A >= a_lower) & (A <= a_higher)
            )[0]

            print(references_genes_index)

            M_trimmed = M[references_genes_index]

            scaled_current_sample_trimmed = scaled_current_sample[non_zero_indexes][
                references_genes_index
            ]
            scaled_reference_sample_trimmed = self.scaled_reference_sample[non_zero_indexes][
                references_genes_index
            ]

            weights = (
                (1 - scaled_current_sample_trimmed) / scaled_current_sample_trimmed
            ) + (
                (1 - scaled_reference_sample_trimmed) / scaled_reference_sample_trimmed
            )


            tmm = 2 ** np.average(M_trimmed, axis = None, weights=weights)


            X_copy = X_copy.astype("float32")
            X_copy[s] = X_copy[s] / tmm

        return X_copy


class RankedExpression(BaseEstimator, TransformerMixin):
    """TODO"""

    def fit(self, X, **fit_params): 
        self.number_of_genes = X.shape[1]
        return(self)

    def transform(self, X, **transform_params):

        if X.shape[1] == self.number_of_genes:
            raise ValueError(
                "X has a different number of genes than the set this normalizer was fitted on"
            )

        X_copy = X.argsort(axis=1)/argsort(axis=1)
        return(X_copy)