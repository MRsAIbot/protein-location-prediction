'''
@author: Tobias Rijken
@date: 23 February 2015
'''

class FeatureBuilder(object):
	"""docstring for FeatureBuilder"""

	tricks = []

	def __init__(self, data):
		super(FeatureBuilder, self).__init__()
		self.raw_data = data
	
	def compute_sequence_length(record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- integer that represents the length of 
		'''
		return len(record.seq)

	def compute_features(feature_list):
		pass

	def get_trainset():
		pass