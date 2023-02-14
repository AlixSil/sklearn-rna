from preprocessing.rnanormaliser import DESeqNormalizer
import pandas as pd
import numpy as np
import pytest


@pytest.fixture
def example_df_fitted_1():
	nrm = DESeqNormalizer()

	test_dataframe = pd.DataFrame({
		"A" : [1,1,1],
		"B" : [2,2,0],
		"C" : [8,1,1]
		}) 

	nrm.fit(test_dataframe)
	return(nrm)

@pytest.fixture
def example_np_fitted_1():
	nrm = DESeqNormalizer()
	test_mtx = np.array([
		[1,2,8],
		[1,2,1],
		[1,0,1]
		])	

	nrm.fit(test_mtx)
	return(nrm)

def test_np_negative_fit():
	nrm = DESeqNormalizer()
	test_mtx = np.array([
		[1,2,8],
		[1,2,-1],
		[1,0,1]
		])	

	with pytest.raises(ValueError):
		nrm.fit(test_mtx)

def test_np_zero_fit():
	nrm = DESeqNormalizer()
	test_mtx = np.array([
		[0,2,0],
		[1,2,-1],
		[1,0,1]
		])	

	with pytest.raises(ValueError):
		nrm.fit(test_mtx)

def test_removal_of_genes_df_1(example_df_fitted_1) :
	nrm = example_df_fitted_1
	assert np.array_equal(nrm.non_zero_genes_indexes, np.array([0,2]))

def test_removal_of_genes_np_1(example_np_fitted_1) :
	nrm = example_np_fitted_1
	assert np.array_equal(nrm.non_zero_genes_indexes, np.array([0,2]))

def test_computation_of_gmean_df_1(example_df_fitted_1) :
	nrm = example_df_fitted_1
	print(nrm.non_zero_genes_gmean_values)
	assert np.array_equal(nrm.non_zero_genes_gmean_values, np.array([1,2]))

def test_computation_of_gmean_np_1(example_np_fitted_1) :
	nrm = example_np_fitted_1
	assert np.array_equal(nrm.non_zero_genes_gmean_values, np.array([1,2]))

def test_transform_np_1(example_np_fitted_1) :
	nrm = DESeqNormalizer()
	test = np.array([
		[1,2,2],
		[4,2,8]
		])
	nrm.fit(test)
	after_transform = nrm.transform(test)
	assert np.array_equal(after_transform, np.array([
		[2.,4.,4.],
		[2.,1.,4.]
		]))

