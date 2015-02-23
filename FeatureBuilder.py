'''
@author: Tobias Rijken
@date: 23 February 2015
@affiliation: University College London
'''

class FeatureBuilder(object):
	"""docstring for FeatureBuilder"""

	tricks = []

	def __init__(self, data):
		super(FeatureBuilder, self).__init__()
		self.raw_data = data
	
	def compute_sequence_length(record):
		return len(record.seq)