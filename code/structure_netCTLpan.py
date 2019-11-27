# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:42:33 2019

@author: sebas

This script was made for the preparation of epitope prediction data based on the netCTLpan out.
"""

import pandas as pd
import numpy as np

pd.set_option('display.width', 200)

def structure_MHC(filename, separator=','):
    ''' 
    Takes a cleaned up netCTLpan file as input.
    The default separator is a ',' (comma).
    The function returns a cleaned up pandas dataframe.
    '''
    headers = ['Index in protein' ,'Protein', 'Allele', 'Peptide', 'MHC score', 'TAP score', 'Cle score', 'Comb score', '%Rank', 'sign']
    df = pd.read_csv(filename, sep=separator, 
                     header = 0, 
                     names=headers)
    df.reset_index(drop=True, inplace=True)
    new = df['Protein'] = df['Protein'].str.split('|', expand=True) # Splits the protein column in 3 new columns
    df['sp'] = new[0]
    df['Uniprot ID'] = new[1] # Uniprot ID of proteins from which the peptides originate
    df['Protein name'] = new[2] # Canonical name of proteins from which the peptides originate
    df.drop(columns=['Protein', 'sp', 'sign'], inplace=True)
    return df

test = structure_MHC('my_input_df.txt')
print(test.head())
