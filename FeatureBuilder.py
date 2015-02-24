'''
@author: Tobias Rijken
@date: 23 February 2015
'''

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class FeatureBuilder(object):
	"""docstring for FeatureBuilder"""

	def __init__(self, data, feature_list):
		super(FeatureBuilder, self).__init__()
		self.raw_data = data
		self.X = np.empty([len(data),len(feature_list)])
		self.Y = np.array([tup[1] for tup in data])
		self.feature_list = feature_list
	
	def sequence_length(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- integer: that represents the length of the protein sequence
		'''
		return len(record.seq)

	def molecular_weight(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- float: representing the molecular weight of the protein
		'''
		PA = ProteinAnalysis(str(record.seq))
		return PA.molecular_weight

	def isoelectric_point(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- float: representing the isoelectric point
		'''
		PA = ProteinAnalysis(str(record.seq))
		return PA.isoelectric_point()

	def amino_acid_composition(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- dictionary: representing the distribution of amino acids in
			the sequence
		'''
		PA = ProteinAnalysis(str(record.seq))
		return PA.get_amino_acids_percent()

	def amino_acid_composition_first50(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- dictionary: representing the distribution of amino acids over the
			first 50 amino acids
		'''
		PA = ProteinAnalysis(str(record.seq)[:50])
		return PA.get_amino_acids_percent()

	def compute_features(self, feature_list = None):
		'''
		Input:
			- feature_list: list of strings
		Output:
			- computes the features as specified in feature_list for every
			for every data entry and outputs an numpy array
		'''

		'''
		TODO:
		- feature normalisation
		'''

		## Set default feature list to the class feature list
		feature_list = feature_list or self.feature_list

		if feature_list != self.feature_list:
			self.X = np.empty(len(self.raw_data),len(feature_list))

		for i,f in enumerate(feature_list):
			if f == "seq_len":
				# TODO: compute feature for every data point 
				self.X[:,i] = np.array([self.sequence_length(tup[0]) for tup in self.raw_data])
			elif f == "":
				# TODO: ...
				pass

	def get_trainset(self):
		return self.X, self.Y