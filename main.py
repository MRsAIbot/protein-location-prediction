'''
@author: Tobias Rijken
@date: 23 February 2015
'''

from Bio import SeqIO
from sklearn import cross_validation

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
	feature_model1 = fb.FeatureBuilder(labeled_data,['seq_len'])
	feature_model1.compute_features()
	X,Y = feature_model1.get_trainset()

	print X[:20]
	print Y[:10]

	## Crossvalidation
	skf = cross_validation.StratifiedKFold(Y,n_folds=5)

	# print len(labeled_data) # returns 9222
	print "Hello World"


if __name__ == '__main__':
	main()