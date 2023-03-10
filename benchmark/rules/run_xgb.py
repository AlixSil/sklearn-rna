import sys

project = sys.argv[1]
normalizerType = sys.argv[2]
random_state = int(sys.argv[3])
outfile = sys.argv[4]


# Making the "package" accessible To be made properly later
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/home/alix/Programmation/sklearn-rna/')

import pandas as pd
import xgboost as xgb
import optuna
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
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


X_train, X_test, y_train, y_test = train_test_split(
			X,
			y,
			test_size = 0.5,
			random_state = random_state
		)


def objective(trial):
	pip = Pipeline( [
		("normalizer", nrm()), 
		("XGBoost", xgb.XGBRegressor(objective = "survival:cox",
			 eval_metric = "cox-nloglik",
			 learning_rate = trial.suggest_float("learning_rate", 0.01, 0.5),
			 max_depth = trial.suggest_int("max_depth", 1,8),
			 colsample_bytree = trial.suggest_float("colsample_bytree", 0.05, 1)
			 ))
		])
	scores = cross_val_score(pip, X_test, y_test)
	result = scores.mean()
	return(result)

study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=25)

print(study.best_params)

pip = Pipeline( [
		("normalizer", nrm()), 
		("XGBoost", xgb.XGBRegressor(objective = "survival:cox",
			 eval_metric = "cox-nloglik",
			 **study.best_params
			 ))
		])

pip.fit(X_train, y_train)

y_pred_risk = pip.predict(X_test)
y_test_events = y_test >= 0
y_test_time = abs(y_test)


try :
	CI = concordance_index_censored(y_test_events, y_test_time, y_pred_risk)[0]
except ValueError :
	CI = None

final_results.append({
	"CI" : CI,
	"random_state" : random_state,
	"normalizerType" : normalizerType,
	"project" : project,
	"learning_rate" : study.best_params["learning_rate"],
	"max_depth" : study.best_params["max_depth"],
	"colsample_bytree" : study.best_params["colsample_bytree"],
	"trial_score" : study.best_trial.value
	})

final_results = pd.DataFrame.from_records(final_results)

final_results.to_csv(outfile)