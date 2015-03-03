'''
@author: Tobias Rijken
@date: 23 February 2015
'''

import numpy as np
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParamData, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter

class FeatureBuilder(object):
	"""docstring for FeatureBuilder"""

	amino_acids = ['G','P','A','V','L','I','M','C','F','Y','W','H','K',\
			'R','Q','N','E','D','S','T']

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
		# PA = ProteinAnalysis(str(record.seq))
		# return PA.length

	def molecular_weight(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- float: representing the molecular weight of the protein
		'''
		PA = ProteinAnalysis(str(record.seq))

		counter = Counter(str(record.seq))
		non_prot_count = sum([v for k,v in counter.items() if k not in self.amino_acids])
		cleaned_seq = Seq(''.join(c for c in str(record.seq) if c in self.amino_acids), IUPAC.protein)

		mol_weight = molecular_weight(seq=cleaned_seq, monoisotopic=PA.monoisotopic)
		avg_mol_weight = mol_weight/float(len(cleaned_seq))

		return mol_weight + non_prot_count * avg_mol_weight

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

	def gravy(self, record):
		'''
		Input:
			- record: a SeqRecord
		Ouput:
			- float: representing the gravy according to Kyte and Doolittle
		'''
		PA = ProteinAnalysis(str(record.seq))
		gravy_list = [ProtParamData.kd[aa] for aa in PA.sequence if aa in ProtParamData.kd]
		total_gravy = sum(gravy_list) 
		return total_gravy / float(len(gravy_list))

	def secondary_structure(self, record):
		'''
		Input:
			- record: a SeqRecord
		Output:
			- tuple of integers
		'''
		PA = ProteinAnalysis(str(record.seq))
		return PA.secondary_structure_fraction()

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
		if "sec_str" in feature_list:
			feature_count += 2
		self.X = np.empty([len(self.raw_data),feature_count])

		i = 0
		for f in feature_list:
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
			elif f == "sec_str":
				self.X[:,i:i+3] = a = np.array([np.asarray(self.secondary_structure(tup[0])) for tup in self.raw_data])
				i += 3
			elif f == "gravy":
				self.X[:,i] = np.array([self.gravy(tup[0]) for tup in self.raw_data])
				i += 1


	def get_dataset(self):
		return self.X, self.Y