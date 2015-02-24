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
		return PA.molecular_weight()

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
		amino_acids = ['G','P','A','V','L','I','M','C','F','Y','W','H','K',\
			'R','Q','N','E','D','S','T']

		## Set default feature list to the class feature list
		feature_list = feature_list or self.feature_list

		feature_count = len(feature_list)
		if "gl_aac" in feature_list:
			feature_count += 19
		self.X = np.empty([len(self.raw_data),feature_count])

		for f in feature_list:
			i = 0
			if f == "seq_len":
				self.X[:,i] = np.array([self.sequence_length(tup[0]) for tup in self.raw_data])
				i += 1
			elif f == "mol_wght":
				self.X[:,i] = np.array([self.molecular_weight(tup[0]) for tup in self.raw_data])
				i += 1
			elif f == "isoelec":
				self.X[:,i] = np.array([self.isoelectric_point(tup[0]) for tup in self.raw_data])
				i += 1
			elif f == "gl_aac":
				for aa in amino_acids:
					self.X[:,i] = np.array([self.amino_acid_composition(tup[0])[aa] for tup in self.raw_data])
					i += 1

	def get_dataset(self):
		return self.X, self.Y