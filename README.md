# Chemical space analysis
Generating multidimensional analyses of chemical compounds with ease

This script is designed to generate a multidimensional analysis of chemical compounds. It is based on the RDKit library and uses the functional groups from the _Fragments_ library in RDKit to generate a multidimensional representation of the chemical space. The script will read a csv file with 2 columns: the first column is the SMILES representation of the compound and the second column is metadata from the compound. The script will generate a multidimensional representation of the chemical space and will generate a t-SNE plot of the chemical space. 

The input must be a csv file with three columns where the first **must** be named as `smiles` and the second **must** be named `name`. The third one can change or can be empty, **but the file should have 3 columns**. It should have the following format:

```
smiles, name, metadata
CC(C)(CN(CC1)CCC1NCc1nc([nH]cc2)c2cc1)O, name_1, active
CN(C)C1CCN(CCCNc2c3[s]c(N(C)C)nc3ncn2)CC1, name_2, active
CCCN1CC(CNC(NC(CN(C)C2)c3c2cccc3)=O)CC1, name_3, inactive
...
``` 

There is an example of a csv file with the right format.

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

Besides these required parameters, you can also set the following optional parameters:

```
- `--threads`: the number of threads to use. Default is 1.
- `--perplexity`: the perplexity of the t-SNE algorithm. Default is 30.
- `--iterations`: the number of iterations of the t-SNE algorithm. Default is 1000.
- `--perplexity_grid`: it will bypass the `--perplexity` parameter and will run the t-SNE algorithm with 
a grid of perplexities ([5, 10, 20, 30, 40, 50]). Default is False. This option takes a while to compute.
- `--umap`: it will calculate the UMAP coordinates of the compounds with default parameters (n_neighbors=15, min_dist=0.1). 
Default is False.
```

An example of usage where it computes a t-SNE with perplexity of 30, max_iter of 500, and a UMAP analysis with 8 threads is shown below:

```bash
python chem_spacer.py -i ./data/data_table.csv -o results -n 500 -t 8 -u 
```

TODO:
- [x] Add UMAP support
- [ ] Add the option to use the Morgan fingerprints
- [ ] Add the option to use the ECFP fingerprints
- [x] Support for Drug names


Any feedback is welcome!