import sys

project = sys.argv[1]
normalizerType = sys.argv[2]
number_iteration = int(sys.argv[3])
outfile = sys.argv[4]


# Making the "package" accessible To be made properly later
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/home/alix/Programmation/sklearn-rna/')

import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sksurv.metrics import concordance_index_censored

if normalizerType == "TMM" :
	from preprocessing.rnanormaliser import TMM as nrm
elif normalizerType == "DESeqNormalizer" :
	from preprocessing.rnanormaliser import DESeqNormalizer as nrm
elif normalizerType == "RankedExpression" :
	from preprocessing.rnanormaliser import RankedExpression as nrm
else :
	raise ValueError("normalizerType is invalid")


# Load data
RNA_data = pd.read_csv("Formatted_data/{}/raw_counts.csv".format(project))
clinical_data = pd.read_csv("Formatted_data/{}/clinical.csv".format(project))

# Format data in numpy X and y

## Filter for mismatch (shouldn't happen)
patients = list(set(clinical_data.bcr_patient_barcode) & set(RNA_data.columns))
clinical_data = clinical_data.loc[ clinical_data.bcr_patient_barcode.isin(patients)]

## Check the order of patients is the correct one
RNA_data = RNA_data.loc[:,patients]
clinical_data = clinical_data.set_index("bcr_patient_barcode")
clinical_data = clinical_data.loc[patients, :]

##numpying and transposing RNA data
X = RNA_data.to_numpy().T

##Creating the correct y vector (time od death positive, censored time negative)
clinical_data["y"] = clinical_data.days_to_death.combine_first(clinical_data.days_to_last_follow_up)
clinical_data.loc[ clinical_data.vital_status =="Alive", "y"] = -clinical_data.loc[ clinical_data.vital_status =="Alive", "y"] 

y = clinical_data.y.to_numpy()

# separate train and test sets
final_results = []
for random_state in range(number_iteration):
	X_train, X_test, y_train, y_test = train_test_split(
			X,
			y,
			test_size = 0.3,
			random_state = random_state
		)

	# Run crossvalidation on the train test to choose the best set of parameters (TODO : run something smarter)

	## Define the general model with a pipeline
	pip = Pipeline( [
		("normalizer", nrm()), 
		("XGBoost", xgb.XGBRegressor(objective = "survival:cox", eval_metric = "cox-nloglik"))
		])

	## Define the space of exploration (Prototypal)

	params = {
		"XGBoost__learning_rate" : [0.01],
		"XGBoost__max_depth" : [2],
	}

	gscv = GridSearchCV(pip, params, verbose=2, n_jobs=1)
	gscv.fit(X_train,y_train)

	y_pred_risk = gscv.predict(X_test)
	y_test_events = y_test >= 0
	y_test_time = abs(y_test)


	CI = concordance_index_censored(y_test_events, y_test_time, y_pred_risk)[0]
	final_results.append({
		"CI" : CI,
		"random_state" : random_state,
		"normalizerType" : normalizerType,
		"project" : project
		})

final_results = pd.DataFrame.from_records(final_results)

final_results.to_csv(outfile)

