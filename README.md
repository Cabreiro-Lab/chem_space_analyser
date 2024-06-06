# chem_space_analyser
Generating multidimensional analyses of chemical compounds with ease

This script is designed to generate a multidimensional analysis of chemical compounds. It is based on the RDKit library and uses the functional groups from the Fragments in RDKit to generate a multidimensional representation of the chemical space. The script will read a csv file with 2 columns: the first column is the SMILES representation of the compound and the second column is metadata from the compound. The script will generate a multidimensional representation of the chemical space and will generate a t-SNE plot of the chemical space. 

The input must be a csv file with two columns where the first one is named as `smiles`. It may have the following format:

```
smiles, metadata
C1=CC=C(C=C1)C(=O)O, active
C1=CC=C(C=C1)C(=O)O, active
C1=CC=C(C=C1)C(=O)O, inactive
...
``` 

It will create a folder where it will save these output files:
- functional_groups.csv: a csv file with the functional groups of the compounds
- tsne_results.csv: a csv file with the t-SNE coordinates of the compounds
- tsne_plot.html: a html file with the t-SNE plot of the compounds

You need to install these dependencies:

```bash
pip install rdkit

conda install plotly pandas scikit-learn
```

Usage:

```bash 

python chem_spacer.py -i input.csv -o output_folder

```
