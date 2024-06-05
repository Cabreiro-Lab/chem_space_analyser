#!/usr/bin/env python3
# description: script to extract functional groups from a list of molecules and analyse them with t-SNE and UMAP

from rdkit import Chem
import pandas as pd
from rdkit.Chem import Fragments
import re
from optparse import OptionParser 
from tqdm import tqdm
from multiprocessing import get_context
import glob
import os
from sklearn.manifold import TSNE
import plotly.express as px

# parameters
parser = OptionParser()

parser.add_option("-i", "--inp",
                    dest = "input",
                    help = "input file in form of a csv file with with at least one column named 'smiles', and other colums as metadata to be preserved",
                    metavar = "INPUT")

parser.add_option("-o", "--out",
                    dest = "output",
                    help = "name of the output file that will store the chemical groups in a csv",
                    metavar = "OUTPUT")

parser.add_option("-t", "--threads",
                    default=1,
                    dest = "threads",
                    help = "number of threads to use",
                    metavar = "THREADS")

parser.add_option("-p", "--perplexity",
                    default=30,
                    dest = "perplexity",
                    help = "perplexity for t-SNE",
                    metavar = "PERPLEXITY")

parser.add_option("-n", "--iterations",
                    default=1000,
                    dest = "iterations",
                    help = "number of iterations for t-SNE",
                    metavar = "ITERATIONS")

parser.add_option("-v", "--version",
                    action="store_true",
                    dest="version",
                    help="print the version of the script")

(options, args) = parser.parse_args()

threads = int(options.threads)

# version info from the script
def version():
    """
    Returns the version of the script.
    """
    return '0.0.1'

def create_output_folder(output):
    """
    Create the output folder if it does not exist.
    """
    if not os.path.exists(output):
        os.makedirs(output)

def read_input(input_file):
    """
    Read the input file and return a pandas dataframe.
    """

    # read the input file
    df = pd.read_csv(input_file)

    # check if the input file has the required columns
    if 'smiles' not in df.columns:
        raise ValueError("The input file does not have a column named 'smiles'.")
    
    # check for duplicates
    if df.duplicated(subset='smiles').any():
        print('The input file has duplicates. Removing them.')

    # remove duplicates and NA values
    df = df.drop_duplicates(subset='smiles')
    df = df.dropna(subset=['smiles'])

    return df


def get_func_gr(molecule):
    '''
    Parameters
    ----------
    molecule : RDKit molecule type
        Molecule object from the RDKit library in python.

    Returns
    -------
    Returns a vector of the functional groups present in the molecule.
    '''
    
    vector = []
    
    for func in frag_functions:
        str_to_eval = 'Fragments.'+func + '(molecule)'
        number = eval(str_to_eval)
        
        vector.append(number)
        
    return vector

# get the functions within Fragments and store them in a list
frag = dir(Fragments)
r = re.compile("fr_.*")
frag_functions = list(filter(r.match, frag)) 

def calculate_func_gr(input):
    print(f"\nReading input file: {options.input}\n")
    smiles = read_input(input)

    smiles['mols'] = smiles['smiles'].apply(Chem.MolFromSmiles) # convert smiles to mols
    smiles.dropna(inplace=True) # remove na in mols
    smiles.reset_index(drop=True, inplace=True) # reset index

    # calculate the functional groups
    print(f"\nCalculating functional groups with {threads} threads.\n")
    p = get_context("fork").Pool(int(threads))

    results = list(tqdm(p.imap(get_func_gr, smiles['mols']), total=smiles.shape[0]))

    p.close()

    # create a dataframe with the results
    func_groups_df = pd.DataFrame(results, columns=frag_functions)

    # concatenate the two dataframes
    final_df = pd.concat([smiles, func_groups_df], axis=1)
    # remove the mols column
    final_df.drop(columns=['mols'], inplace=True)

    # save the final dataframe
    final_df.to_csv(f"{options.output}/functional_groups.csv", index=False)

    return final_df

def tsne_calc(final_df):
    """
    Calculate the t-SNE for the final dataframe.
    
    Parameters
    
    final_df : pandas dataframe that comes from the calculate_func_gr function
    
    Returns
    
    tsne_df : pandas dataframe
    """
    print(f"\nCalculating t-SNE with perplexity {options.perplexity} and {options.iterations} iterations. This can take a while, go grab a coffee.\n")
    tsne = TSNE(n_components=3, random_state=123,
            perplexity=int(options.perplexity), 
            max_iter=int(options.iterations), 
            verbose=1,
            n_jobs=threads)

    # calculate the tSNE from the column fr_Al_COO to the last column
    X = final_df.iloc[:, 2:].values
    X_embedded = tsne.fit_transform(X)

    # create a dataframe with the tSNE results
    tsne_df = pd.DataFrame(X_embedded, columns=['tsne_1', 'tsne_2', 'tsne_3'])
    # add the two original columns from the final_df
    tsne_df = pd.concat([final_df.iloc[:, 0:2], tsne_df], axis=1)

    print("\nt-SNE calculation finished.\n")

    # save the t-SNE dataframe
    tsne_df.to_csv(f"{options.output}/tsne_results.csv", index=False)

    return tsne_df

def plotly_tsne(tsne_df):
    """
    Plot the t-SNE with plotly.
    
    Parameters
    
    tsne_df : pandas dataframe that comes from the tsne_calc function
    
    Returns
    
    None
    """
    columns = tsne_df.columns

    fig = px.scatter_3d(tsne_df, x='tsne_1', y='tsne_2', z='tsne_3', color=columns[1], hover_data=['smiles'])
    
    # save the plot in the output folder
    fig.write_html(f"{options.output}/tsne_plot.html")

    print(f"\nPlot saved in {options.output}/tsne_plot.html\n")



def main():
    # check the version
    if options.version:
        print(version())
        return

    # create the output folder
    create_output_folder(options.output)

    # calculate the functional groups
    final_df = calculate_func_gr(options.input)

    # calculate the t-SNE
    tsne_df = tsne_calc(final_df)

    # plot the t-SNE
    plotly_tsne(tsne_df)

if __name__ == '__main__':
    main()