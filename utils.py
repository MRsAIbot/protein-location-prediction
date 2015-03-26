'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from operator import itemgetter

import cPickle as pickle
import csv
import FeatureBuilder as fb
import matplotlib.pyplot as plt
import numpy as np
import random


def load_data(pkl, blind=False):
	if pkl:
		if not blind:
			with open("dataframe.pkl","rb") as fp:
				feature_model1 = pickle.load(fp)
		elif blind:
			with open("blind_dataframe.pkl","rb") as fp:
				feature_model1 = pickle.load(fp)
	else:
		## Define the filenames
		if blind:
			files = ['blind.fasta']
		else:
			files = ['cyto.fasta', 'mito.fasta', 'nucleus.fasta', 'secreted.fasta']

		## Read and label data ('C' for cytosolic, 'S' for secreted, 'M' for mitochondrial, 'N' for nuclear proteins respectively)
		labeled_data = []
		for datafile in files:
			r = read_data(datafile)
			# print "{0}: {1}".format(datafile,len(r))
			if blind:
				labeled_data += label_data(r)
			else:
				labeled_data += label_data(r,datafile[0].upper())

		## Shuffle the data
		if not blind:
			random.shuffle(labeled_data)

		## Model 1: sequence length
		feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len','sec_str',\
			'isoelec','gl_aac','gravy', 'mol_wght', 'loc_aac_first', 'loc_aac_last'])
		# feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len'])
	if blind:
		return feature_model1, r
	else:
		return feature_model1


def label_data(records, label='None'):
	'''
	Input:
		- records: list of sequence records
		- label: string
	Output:
		- list of tuples (x,y) where x is the sequence record and y is the string denoting the label
	'''
	return zip(records, [label]*len(records))


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


def blind_test_to_csv(clf):
	feature_model_blind, records = load_data(pkl=False, blind=True)
	seq_ids = [record.id for record in records]

	X = feature_model_blind.get_dataset()[0]

	class_probs = clf.predict_proba(X)
	predictions = clf.predict(X)
	class_probs_max = class_probs.max(axis=1)

	with open("blind_prediction.csv","wb") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(('Sequence ID', 'Predicted Label', 'Confidence'))
		for i in xrange(len(seq_ids)):
			writer.writerow((seq_ids[i], predictions[i], class_probs_max[i]))


def plot_confusion_matrix(cm, labels, title='Confusion matrix', cmap=plt.cm.Blues):
	full_names = {'C':'Cytosolic', 'M':'Mitochondrial', 'N':'Nuclear', 'S':'Secreted'}
	labels = [full_names[i] for i in labels]
	plt.imshow(cm, interpolation='nearest', cmap=cmap)
	plt.title(title)
	plt.colorbar()
	tick_marks = np.arange(len(labels))
	plt.xticks(tick_marks, labels, rotation=45)
	plt.yticks(tick_marks, labels)
	plt.tight_layout()
	plt.ylabel('True label')
	plt.xlabel('Predicted label')


def accuracy_per_class(cm, labels):
	accuracy_overall = float(np.trace(cm))/np.sum(np.sum(cm, axis=0))
	for l in xrange(len(labels)):
		accuracy = float(cm[l,l])/sum(cm[l,:])
		print "{0}: {1:.3f}".format(labels[l], accuracy)
	print "Acc. overall: {0:.3f}".format(accuracy_overall)
	return 0


def plot_feature_importance(forest):
	importances = forest.feature_importances_
	std = np.std([tree.feature_importances_ for tree in forest.estimators_],
	             axis=0)
	indices = np.argsort(importances)[::-1]

	# Print the feature ranking
	print("Feature ranking:")

	for f in range(len(importances)):
	    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

	# indices = indices[:15]
	# Plot the feature importances of the forest
	plt.figure()
	plt.title("Feature importances")
	plt.bar(range(len(importances)), importances[indices],
	       color="r", yerr=std[indices], align="center")
	plt.xticks(range(len(importances)), indices)
	plt.xlim([-1, len(importances)])
	plt.ylabel('Feature importances')
	plt.xlabel('Feature number')
	plt.show()


def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")