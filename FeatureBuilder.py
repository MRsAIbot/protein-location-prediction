'''
@author: Tobias Rijken
@date: 23 February 2015
'''

import numpy as np

class FeatureBuilder(object):
	"""docstring for FeatureBuilder"""

	def __init__(self, data, feature_list):
		super(FeatureBuilder, self).__init__()
		self.raw_data = data
		self.X = np.empty([len(data),len(feature_list)])
		self.Y = np.array([tup[1] for tup in data])
		self.feature_list = feature_list
	
	def compute_sequence_length(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- integer that represents the length of 
		'''
		return len(record.seq)

	def compute_features(self, feature_list = None):
		'''
		Input:
			- feature_list: list of strings
		Output:
			- computes the features as specified in feature_list for every
			for every data entry and outputs an numpy array
		'''
		## Set default feature list to the class feature list
		feature_list = feature_list or self.feature_list

		if feature_list != self.feature_list:
			self.X = np.empty(len(self.raw_data),len(feature_list))

		for i,f in enumerate(feature_list):
			if f == "seq_len":
				# TODO: compute feature for every data point 
				self.X[:,i] = np.array([self.compute_sequence_length(tup[0]) for tup in self.raw_data])
			elif f == "":
				# TODO: ...
				pass

	def get_trainset(self):
		return self.X, self.Y