'''
@author: Tobias Rijken
@date: 23 February 2015
@affiliation: University College London
'''

from Bio import SeqIO
import numpy as np
import FeatureBuilder as fb

files = ['cyto.fasta', 'mito.fasta', 'nucleus.fasta', 'secreted.fasta']

def read_data(filename):
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

def compute_features():
	pass

def construct_trainset():
	pass

def main():
	## Read and label data ('C' for cytosolic, 'S' for secreted, 'M' for mitochondrial, 'N' for nuclear proteins respectively)
	labeled_data = []
	for datafile in files:
		r = read_data(datafile)
		labeled_data += label_data(r,datafile[0].upper())

	# print len(labeled_data) # returns 9222
	print "Hello World"

if __name__ == '__main__':
	main()