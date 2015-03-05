'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.svm import SVC

import cPickle as pickle
import FeatureBuilder as fb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random

random.seed(1)

def plot_confustion_matrix(cm):
	plt.matshow(cm)
	plt.title('Confusion matrix')
	plt.colorbar()
	plt.ylabel('True label')
	plt.xlabel('Predicted label')
	plt.show()
	return 0


def read_data(filename):
	'''
	Input:
		- filename: string representing a filename in FASTA format
	Output:
		- list of sequence records
	'''
	handle = open("data/" + filename, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	return records


def label_data(records, label):
	'''
	Input:
		- records: list of sequence records
		- label: string
	Output:
		- list of tuples (x,y) where x is the sequence record and y is the string denoting the label
	'''
	return zip(records, [label]*len(records))


def main():
	## Define the filenames
	files = ['cyto.fasta', 'mito.fasta', 'nucleus.fasta', 'secreted.fasta']

	## Read and label data ('C' for cytosolic, 'S' for secreted, 'M' for mitochondrial, 'N' for nuclear proteins respectively)
	labeled_data = []
	for datafile in files:
		r = read_data(datafile)
		labeled_data += label_data(r,datafile[0].upper())

	## Shuffle the data
	random.shuffle(labeled_data)

	## Model 1: sequence length
	feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len','sec_str',\
		'isoelec','gl_aac','gravy', 'mol_wght'])
	# feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len'])

	feature_model1.compute_features()
	X,Y = feature_model1.get_dataset()

	# with open("save.p","wb") as fp:
	# 	pickle.dump(X,fp)

	# print X[:10]
	# print Y[:10]

	## Normalise the data
	X_norm = preprocessing.normalize(X, axis=0, norm='l2')

	print X_norm[:10]

	## Split the dataset in train and test sets
	X_train, X_test, Y_train, Y_test = cross_validation.train_test_split(
		X_norm, Y, test_size=0.25, random_state=1)

	# ## Train SVM
	# clf = SVC(kernel='linear', cache_size=1000)

	skf = cross_validation.StratifiedKFold(y=Y, n_folds=5)
	# kf = cross_validation.KFold(n=len(Y), n_folds=5)
	# scores_svm = cross_validation.cross_val_score(clf, X_norm, Y, cv=skf)

	# print("Accuracy: %0.4f (+/- %0.4f)" % (scores_svm.mean(), scores_svm.std() * 2))
	# print scores_svm

	## Single train SVM
	# clf1 = SVC(kernel='rbf')
	# clf1.fit(X_train, Y_train)
	# print clf1.score(X_test, Y_test)

	# Set the parameters by cross-validation
	tuned_parameters = [{'kernel': ['rbf'], 'gamma': 10.0 ** np.arange(-5,4),\
		'C': 10.0 ** np.arange(-2,9)},
		{'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

	# scores = ['precision', 'recall']

	# # for score in scores:
	# print("# Tuning hyper-parameters for score")
	# print()

	# clf = GridSearchCV(SVC(C=1), tuned_parameters, cv=5)
	# clf.fit(X_train, Y_train)

	# print("Best parameters set found on development set:")
	# print()
	# print(clf.best_estimator_)
	# print()
	# print("Grid scores on development set:")
	# print()
	# for params, mean_score, scores in clf.grid_scores_:
	#     print("%0.3f (+/-%0.03f) for %r"
	#           % (mean_score, scores.std() / 2, params))
	# print()

	# print("Detailed classification report:")
	# print()
	# print("The model is trained on the full development set.")
	# print("The scores are computed on the full evaluation set.")
	# print()
	# y_true, y_pred = Y_test, clf.predict(X_test)
	# print(classification_report(y_true, y_pred))
	# print()

	## Random Forest
	rf = RandomForestClassifier()
	scores_rf = cross_validation.cross_val_score(rf, X_norm, Y, cv=skf)
	rf.fit(X_train, Y_train)

	print("Accuracy: %0.4f (+/- %0.4f)" % (scores_rf.mean(), scores_rf.std() * 2))
	print scores_rf

	y_true, y_pred = Y_test, rf.predict(X_test)
	print classification_report(y_true, y_pred)
	cm = confusion_matrix(y_true, y_pred)
	print cm
	plot_confustion_matrix(cm)

	## Single train RF
	# rf1 = RandomForestClassifier()
	# rf1.fit(X_train, Y_train)
	# print rf1.score(X_test, Y_test)

	## F! scores


if __name__ == '__main__':
	main()