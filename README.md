# protein-location-prediction

## How to run

First, include data! We assumed the data is in the FASTA format.
In our dataset, every class has its own file with data entries. Thus, the label
is given by the name of the data file.

For example, mito.fasta contains a number of entries of mitochondrial 
protein sequences.

```
|	README.md
|	main.py
|	FeatureBuilder.py
|---data
	|	mito.fasta
	|	nucleus.fasta
	|	secreted.fasta
	|	cyto.fasta
```

## Software packages

I used the following software packages:

1. [Biopython](http://biopython.org/wiki/Main_Page) - tools for biological computation in Python
2. [scikit-learn](http://scikit-learn.org/stable/) - machine learning tools in Python

## Installation & running

### Installation

### Running

Now, you can run the experiment as follows:

```
python main.py
```