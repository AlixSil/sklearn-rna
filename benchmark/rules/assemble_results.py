import pandas as pd
import os
import glob
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

path_to_directory = "Temp"

all_files = glob.glob(os.path.join("Temp" "/**/*.csv"), recursive=True)

all_results = []

for f in all_files :
	all_results.append(pd.read_csv(f, index_col = None))

all_results = pd.concat(all_results, axis = 0)
all_results = all_results.drop(all_results.columns[0], axis = 1)

all_results = all_results.reset_index()

sns.catplot(data = all_results, y = "CI", x = "normalizerType", col = "project", kind = "box")