# protein-location-prediction

## How to run

First, include data! We assumed the data is in the FASTA format.
In our dataset, every class has its own file with data entries. Thus, the label
is given by the name of the data file.

For example, mito.fasta contains a number of entries of mitochondrial 
protein sequences.

'''
|	README.md
|	main.py
|	FeatureBuilder.py
|---data
	|	mito.fasta
	|	nucleus.fasta
	|	secreted.fasta
	|	cyto.fasta
'''

'''
python main.py
'''