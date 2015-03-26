'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC, LinearSVC
from time import time
from utils import *

import cPickle as pickle
import FeatureBuilder as fb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random

random.seed(1)


def run_svm(X_train, X_test, Y_train, Y_test):
	# Set the parameters by cross-validation
	param_grid = {'kernel': ['rbf'], 
		'gamma': 10.0 ** np.arange(-5,4),
		'C': 10.0 ** np.arange(-2,9)}

	svm = SVC()

	grid_search = GridSearchCV(svm, param_grid=param_grid, cv=5)
	start = time()
	grid_search.fit(X_train, Y_train)

	print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
      % (time() - start, len(grid_search.grid_scores_)))
	report(grid_search.grid_scores_)

	with open("svm_gridsearch_95.pkl","wb") as f:
		pickle.dump(grid_search,f)

	labels = np.unique(Y_train)

	y_true, y_pred = Y_test, grid_search.predict(X_test)
	print classification_report(y_true, y_pred)
	cm = confusion_matrix(y_true, y_pred, labels)
	print cm

	accuracy_per_class(cm, labels)

	plt.figure()
	plot_confusion_matrix(cm, labels)

	cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

	plt.figure()
	plot_confusion_matrix(cm_normalized, labels, title='Normalized confusion matrix')

	plt.show()

	return 0


def run_random_forest(X_train, X_test, Y_train, Y_test):
	param_grid = {"max_depth": [3, None],
              "max_features": [10, "auto", None],
              "min_samples_split": [1, 3, 10],
              "min_samples_leaf": [1, 10, 15, 20],
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

	rf = RandomForestClassifier()
	# scores_rf = cross_validation.cross_val_score(rf, X_norm, Y, cv=skf)
	grid_search = GridSearchCV(rf, param_grid=param_grid, cv=5)
	start = time()
	grid_search.fit(X_train, Y_train)

	print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
      % (time() - start, len(grid_search.grid_scores_)))
	report(grid_search.grid_scores_)

	# print("Accuracy: %0.4f (+/- %0.4f)" % (scores_rf.mean(), scores_rf.std() * 2))
	# print scores_rf
	with open("rf_gridsearch_95.pkl","wb") as f:
		pickle.dump(grid_search,f)

	labels = np.unique(Y_train)

	y_true, y_pred = Y_test, grid_search.predict(X_test)
	print classification_report(y_true, y_pred)
	cm = confusion_matrix(y_true, y_pred, labels)
	print cm

	accuracy_per_class(cm, labels)

	plt.figure()
	plot_confusion_matrix(cm, labels)

	cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

	plt.figure()
	plot_confusion_matrix(cm_normalized, labels, title='Normalized confusion matrix')

	plt.show()

	return 0


def run_log_reg(X_train, X_test, Y_train, Y_test):
	param_grid = {"C": [0.1, 1, 10, 100, 1000],
              "penalty": ["l1", "l2"]}

	log_reg = LogisticRegression()
	grid_search = GridSearchCV(log_reg, param_grid=param_grid, cv=5)
	start = time()
	grid_search.fit(X_train, Y_train)

	print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
      % (time() - start, len(grid_search.grid_scores_)))
	report(grid_search.grid_scores_)

	with open("log_reg_gridsearch_95.pkl","wb") as f:
		pickle.dump(grid_search,f)

	labels = np.unique(Y_train)

	log_reg_best = grid_search.best_estimator_

	y_true, y_pred = Y_test, log_reg_best.predict(X_test)
	print classification_report(y_true, y_pred)
	cm = confusion_matrix(y_true, y_pred, labels)
	print cm

	accuracy_per_class(cm, labels)

	plt.figure()
	plot_confusion_matrix(cm, labels)

	cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

	plt.figure()
	plot_confusion_matrix(cm_normalized, labels, title='Normalized confusion matrix')

	plt.show()

	# log_reg_best = grid_search.best_estimator_
	# print log_reg_best.predict_proba(X_test)


def rf_with_svc_feature_selection(X_train, X_test, Y_train, Y_test):
	param_grid = {
		'clf__max_depth': [None],
		'clf__min_samples_split': [1, 3, 10],
		'clf__min_samples_leaf': [1, 10, 15, 20],
		'clf__max_features': [10, 'auto', None],
		'clf__criterion': ['gini', 'entropy'],
	}

	pipeline = Pipeline([
		('fs', LinearSVC(penalty="l1", dual=False)),
		('clf', RandomForestClassifier())
	])

	grid_search = GridSearchCV(pipeline, param_grid=param_grid, cv=5)
	start = time()
	grid_search.fit(X_train, Y_train)

	print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
      % (time() - start, len(grid_search.grid_scores_)))
	report(grid_search.grid_scores_)

	labels = np.unique(Y_train)

	y_true, y_pred = Y_test, grid_search.predict(X_test)
	print classification_report(y_true, y_pred)
	cm = confusion_matrix(y_true, y_pred, labels)
	print cm

	accuracy_per_class(cm, labels)

	plt.figure()
	plot_confusion_matrix(cm, labels)

	cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

	plt.figure()
	plot_confusion_matrix(cm_normalized, labels, title='Normalized confusion matrix')

	plt.show()


def main():
	pkl = True
	feature_model1 = load_data(pkl=pkl)

	feature_model1.compute_features()
	X,Y = feature_model1.get_dataset()

	labels = np.unique(Y)

	if not pkl:
		with open("dataframe.pkl","wb") as fp:
			pickle.dump(feature_model1,fp)

	# print X[:10]
	# print Y[:10]

	## Normalise the data
	# X_norm = preprocessing.normalize(X, axis=0, norm='l2')
	X_norm = X

	print X_norm[:10]

	## Split the dataset in train and test sets
	X_train, X_test, Y_train, Y_test = cross_validation.train_test_split(
		X_norm, Y, test_size=0.05, random_state=1)
	# skf = cross_validation.StratifiedKFold(y=Y, n_folds=5)

	# X_train, X_test, Y_train, Y_test = X_train[:1000], X_test[:100], Y_train[:1000], Y_test[:100]
	## Train SVM
	# run_svm(X_train, X_test, Y_train, Y_test)

	## Random Forest
	run_random_forest(X_train, X_test, Y_train, Y_test)

	## Logistic Regression
	run_log_reg(X_train, X_test, Y_train, Y_test)

	# rf_with_svc_feature_selection(X_train, X_test, Y_train, Y_test)


if __name__ == '__main__':
	main()