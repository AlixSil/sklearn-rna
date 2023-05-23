from sklearnRNA.rnanormaliser import TMM
import pandas as pd
import numpy as np
import pytest
import random

@pytest.fixture
def example_np_fitted_1():
    test_np = np.array([
       [1, 9, 8, 8, 7, 7, 1, 2, 0, 8, 5, 8, 4, 0, 1, 7, 6, 5, 3, 6],
       [1, 2, 5, 1, 3, 4, 4, 2, 9, 0, 3, 7, 7, 8, 5, 1, 5, 4, 0, 5],
       [6, 7, 7, 7, 0, 8, 2, 7, 2, 9, 6, 9, 6, 6, 7, 2, 3, 0, 5, 7],
       [9, 9, 7, 3, 4, 6, 2, 5, 1, 2, 1, 4, 7, 5, 1, 1, 0, 3, 2, 9],
       [3, 7, 6, 5, 6, 7, 5, 8, 5, 8, 1, 6, 4, 5, 1, 7, 0, 2, 9, 6],
       [6, 6, 3, 8, 6, 9, 6, 2, 8, 8, 5, 5, 3, 4, 1, 5, 1, 0, 0, 2],
       [4, 9, 2, 1, 2, 8, 4, 3, 5, 5, 4, 7, 8, 6, 1, 6, 6, 5, 7, 8],
       [8, 2, 4, 3, 6, 1, 4, 0, 6, 9, 8, 9, 0, 3, 5, 7, 9, 0, 0, 0],
       [1, 1, 3, 2, 9, 5, 6, 4, 0, 2, 0, 8, 3, 4, 3, 1, 5, 1, 1, 6],
       [5, 5, 0, 5, 3, 4, 1, 9, 1, 6, 0, 9, 4, 0, 0, 2, 7, 2, 3, 5]
       ])
    #Reference sample should be the last by default
    return(test_np)

def test_normal_fitting(example_np_fitted_1):
    test_np = example_np_fitted_1
    nrm = TMM()
    nrm.fit(test_np)

    normaly_selected = np.array([5, 5, 0, 5, 3, 4, 1, 9, 1, 6, 0, 9, 4, 0, 0, 2, 7, 2, 3, 5])
    normaly_selected = normaly_selected/normaly_selected.sum()

    assert np.array_equal( nrm.scaled_reference_sample, normaly_selected)

def test_fitting_choose_line_and_library_size(example_np_fitted_1):
    test_np = example_np_fitted_1
    library_sizes = np.array([1 for k in range(10)])


    nrm = TMM()
    nrm.fit(test_np, library_sizes = library_sizes, reference_line = 1)

    assert np.array_equal(nrm.scaled_reference_sample, np.array([1, 2, 5, 1, 3, 4, 4, 2, 9, 0, 3, 7, 7, 8, 5, 1, 5, 4, 0, 5]))

def test_fit_transform(example_np_fitted_1):
    test_np = example_np_fitted_1
    nrm = TMM()
    nrm.fit_transform(test_np)