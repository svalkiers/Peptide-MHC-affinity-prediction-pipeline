# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 09:48:58 2019

@author: sebas

netCTLpan_pipeline v1.0

The netCTLpan pipeline was created for the automated affinity prediction of proteins/peptides to different HLA alleles.
The output of the epitope affinity prediction is processed and converted into a pandas DataFrame object.
"""

###############
# 1 - MODULES #
###############

import pandas as pd
import numpy as np
import re
import subprocess, shlex
from subprocess import Popen, PIPE


#################
# 2 - FUNCTIONS #
#################

# 2.1: execute netCTLpan loops over different HLA alleles

def execute_netCTLpan(allele, protein, length=9):
    '''
    Execute the netCTLpan algorithm from a python function.
    The function takes 3 arguments:
        - allele: path to file containing the names of different HLA alleles.
        - protein: filename of file containing different protein sequences in FASTA format.
        - length: peptide length (default = 9)
    The function loops over every HLA allele in the given file and returns one file containing all the outputs.
    '''
    with open(allele, 'r') as alleles:
        alleles = list(map(lambda s: s.strip(), alleles))    # Converts the allele file into a list
    
    path_to_output = '/home/svalkiers/data_folder/MHC_test_output.txt'
    with open(path_to_output, 'w+') as myoutput:     # output path
        pass

    c = 0    # initialize count at 0 -> corresponds to index of allele
    for i in alleles:
        program = "netCTLpan" + " -a " + str(alleles[c]) + " -f " + str(protein) + " -l " + str(length)
        
        p = Popen([program], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)     # enter commands in terminal
        output, stderr = p.communicate()    # return output
        
        with open(path_to_output, 'a') as myoutput:
            res = myoutput.write(str(output).replace("\\n", "\n"))  # append results to output file
        
        c += 1
        
    return


# 2.2: Clean the netCTLpan output data
    
def remove_subheaders(filename_in, filename_out):
    """ This function takes a raw netCTLpan output file as input
    and removes text inbetween the data columns. The function
    also removes the column headers. """
    with open(filename_in, 'r') as f_in:
        lines = f_in.readlines()
    with open(filename_out, 'w') as f_out:
        for line in lines:
            if line.strip('\n') != '----------------------------------------------------------------------':
                f_out.write(line) #writes all line w/o ----
    with open(filename_out, 'r') as f_out:
        lines = f_out.readlines()
    with open(filename_out, 'w') as f_out:
        for line in lines:
            if 'Number of MHC ligands' not in line.strip('\n') and 'NetCTLpan' not in line.strip('\n') and '#' not in line.strip('\n'):
                f_out.write(line)
                
def remove_empty_lines(filename_in, filename_out):
    """ Removes all whitespace inbetween data columns. """
    with open(filename_in, 'r') as f_in:
        lines = f_in.readlines()
    with open(filename_out, 'w') as f_out:
        for line in lines:
            if line.strip():
                f_out.write(line)
                
def replace_spaces(filename_in, filename_out):
    ''' Replaces every 1, 2 and 3 spaces by a comma. '''
    with open(filename_in, 'r') as fin:
        lines = fin.readlines()
    with open(filename_out, 'w') as fout:
        for line in lines:
            newline = re.sub('   ', ',', line)
            fout.write(newline)
    with open(filename_out, 'r') as fout:
        lines = fout.readlines()
    with open(filename_out, 'w') as fout:
        for line in lines:
            newline = re.sub('  ', ',', line)
            fout.write(newline)
    with open(filename_out, 'r') as fout:
        lines = fout.readlines()
    with open(filename_out, 'w') as fout:
        for line in lines:
            newline = re.sub(' ', ',', line)
            fout.write(str(newline).replace(',,', ','))    


# 2.3: Make a DataFrame from the cleaned netCTLpan output
            
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
    print(df)
    

###############
# 3 - EXECUTE #
###############
    
# 3.1: FASTA + HLA alleles ---> affinity score & rank
execute_netCTLpan('/home/svalkiers/data_folder/HLA_alleles.txt', 'FASTA_example.txt')

# 3.2: Clean output file
remove_subheaders('MHC_test_output.txt', 'MHC_subheader.txt')
remove_empty_lines('MHC_subheader.txt', 'MHC_emptylines.txt')
replace_spaces('MHC_emptylines.txt', 'MHC_cleaned.txt')

# 3.3: Make DataFrame object
structure_MHC('MHC_cleaned.txt')