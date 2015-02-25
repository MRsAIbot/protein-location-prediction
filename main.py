'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from sklearn import cross_validation
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC

import cPickle as pickle
import numpy as np
import FeatureBuilder as fb
import matplotlib.pyplot as plt


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

	## Model 1: sequence length
	feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len','sec_str',\
		'isoelec','gl_aac','gravy'])
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

	## Train SVM
	clf = SVC(kernel='rbf', cache_size=1000)

	skf = cross_validation.StratifiedKFold(y=Y, n_folds=5)
	scores_svm = cross_validation.cross_val_score(clf, X_norm, Y, cv=skf)

	print("Accuracy: %0.4f (+/- %0.4f)" % (scores_svm.mean(), scores_svm.std() * 2))

	# ## Grid parameters
	# c_range = 10.0 ** np.arange(-2,9)
	# gamma_range = 10.0 ** np.arange(-5,4)
	# param_grid = dict(gamma=gamma_range, C=c_range)
	
	# ## Crossvalidation
	# skf = cross_validation.StratifiedKFold(y=Y_train, n_folds=5)
	# grid = GridSearchCV(SVC(), param_grid=param_grid, cv=skf)
	# grid.fit(X_train,Y_train)
	# print "Accuracy after grid search CV: {0}".format(grid.score(X_test, Y_test))

	## Random Forest
	rf = RandomForestClassifier()
	scores_rf = cross_validation.cross_val_score(rf, X_norm, Y, cv=skf)

	print("Accuracy: %0.4f (+/- %0.4f)" % (scores_rf.mean(), scores_rf.std() * 2))


if __name__ == '__main__':
	main()