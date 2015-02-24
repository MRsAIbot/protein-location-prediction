'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from sklearn import cross_validation
from sklearn import preprocessing

import numpy as np
import FeatureBuilder as fb


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

	## Model 1: sequence length
	feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len',\
		'isoelec','gl_aac'])

	feature_model1.compute_features()
	X,Y = feature_model1.get_dataset()

	# print X[:10]
	# print Y[:10]

	# Normalise the data
	X_norm = preprocessing.normalize(X, axis=0, norm='l2')

	# print X_norm[:10]


	## Crossvalidation
	skf = cross_validation.StratifiedKFold(Y,n_folds=5)

	c_range = 10.0 ** np.arange(-2,9)
	gamma_range = 10.0 ** np.arange(-5,4)
	param_grid = dict(gamma=gamma_range, C=c_range)
	cv = cross_validation.StratifiedKFold(y=Y_train, n_folds=3)
	grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv)
	grid.fit(X_train,Y_train)
	print "Accuracy after grid search CV: {0}".format(grid.score(X_test, Y_test))

	# print len(labeled_data) # returns 9222


if __name__ == '__main__':
	main()