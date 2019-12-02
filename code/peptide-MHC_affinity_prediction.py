# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:59:35 2019

@author: sebas
"""

# 1: Modules

from multiprocessing import cpu_count
from joblib import Parallel, delayed
from itertools import product
from time import sleep
import subprocess, shlex
from subprocess import Popen, PIPE
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# 2: Parallel exe version (~ number of cores) of prediction models

    
def predict_epitope_affinity(allele, protein, path, model='netMHCpan', length=9):
    
    commands = str(model) + " -a " + str(allele) + " -f " + str(protein) + " -l " + str(length)     # this will be entered in the command line

    p = Popen([commands], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)     # enter commands in terminal
    output, stderr = p.communicate()
    
    path_to_output = path + 'AFFINITY.txt'
    with open(path_to_output, 'a') as out:
        out.write(str(output).replace("\\n", "\n"))

    return


# 3: Prepare output for analysis

def remove_subheaders(filename_in, filename_out):
    """ This function takes a raw netCTLpan output file as input
    and removes text inbetween the data columns. The function
    also removes the column headers. """
    with open(filename_in, 'r') as f_in:
        lines = f_in.readlines()
    with open(filename_out, 'w') as f_out:
        for line in lines:
            if 'Number of MHC ligands' not in line.strip('\n') and 'NetCTLpan' not in line.strip('\n') and '#' not in line.strip('\n') and '---' not in line.strip('\n') and 'eptide' not in line.strip('\n') and 'Distance' not in line.strip('\n') and '<' not in line.strip('\n'):
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
            

#4: structure output and make immunogenic density plots
            
def structure_netMHCpan(filename, path, separator=','):
    ''' 
    Takes a cleaned up netCTLpan file as input.
    The default separator is a ',' (comma).
    The function returns a cleaned up pandas dataframe.
    
    Parameter description:
    
    - Pos: Residue number (starting from 0)
    - HLA: Molecule/allele name
    - Peptide: Amino acid sequence of the potential ligand
    - Core: The minimal 9 amino acid binding core directly in contact with the MHC
    - Of: The starting position of the Core within the Peptide (if > 0, the method predicts a N-terminal protrusion)
    - Gp: Position of the deletion, if any.
    - Gl: Length of the deletion.
    - Ip: Position of the insertions, if any.
    - Il: Length of the insertion.
    - Icore: Interaction core. This is the sequence of the binding core including eventual insertions of deletions.
    - Identity: Protein identifier, i.e. the name of the Fasta entry.
    - Score: The raw prediction score
    - Aff(nM): Predicted binding affinity in nanoMolar units (if binding affinity predictions is selected).
    - %Rank: Rank of the predicted affinity compared to a set of random natural peptides. 
    This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities. 
    Strong binders are defined as having %rank<0.5, and weak binders with %rank<2. We advise to select candidate binders based on %Rank rather than nM Affinity
    - BindLevel: (SB: strong binder, WB: weak binder). 
    The peptide will be identified as a strong binder if the % Rank is below the specified threshold for the strong binders, by default 0.5%. 
    The peptide will be identified as a weak binder if the % Rank is above the threshold of the strong binders but below the specified threshold for the weak binders, by default 2%.
    '''
    
    my_path = path
    
    header = ['A', 'Index in protein' ,'HLA', 'B', 'Peptide', 'Core', 'Of', 'Gp', 'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank']
    df = pd.read_csv(filename, sep=separator, 
                     names=header)
    new = df['Identity'] = df['Identity'].str.split('_', expand=True) # Splits the protein column in 3 new columns
    df['Uniprot ID'] = new[1] # Uniprot ID of proteins from which the peptides originate
    df['Protein name'] = new[2] # Canonical name of proteins from which the peptides originate
    df.drop(columns=['Identity', 'A', 'B'], inplace=True)
    
    proteins = df['Uniprot ID'].unique().tolist()
    proteins.pop()
    
    for protein in proteins:
        
        df_protein = df['Uniprot ID']==str(protein)
        new_df = df[df_protein]
            
        df_score = new_df.groupby('Index in protein', as_index=False)['Score'].mean()
        x = df_score['Index in protein']
        y = df_score['Score']
            
        df_score_stdev = new_df.groupby('Index in protein', as_index=False)['Score'].std()
            
        df_score['stdev'] = df_score_stdev['Score']
            
        df_score_error = df_score['stdev']
            
        immunogenicity_score = df_score.plot(kind='line', x='Index in protein', y='Score', color='red', figsize=(16,12))
        immunogenicity_score.fill_between(x, y-df_score_error, y+df_score_error, color='blue')
        immunogenicity_score.set_title('Immunogenic density plot ' + str(protein), fontsize=30)
        immunogenicity_score.set_ylabel('MHC score', fontsize=18)
        immunogenicity_score.set_xlabel('Index in protein', fontsize=18)
            
        plt.savefig(my_path + 'immunogenic_density_score-' + str(protein) + '.png')
            
        df_rank = new_df.groupby('Index in protein', as_index=False)['Rank'].mean()
        x = df_rank['Index in protein']
        y = df_rank['Rank']
            
        df_rank_stdev = new_df.groupby('Index in protein', as_index=False)['Rank'].std()
            
        df_rank['stdev'] = df_rank_stdev['Rank']
            
        df_rank_error = df_rank['stdev']
            
        immunogenicity_rank = df_rank.plot(kind='line', x='Index in protein', y='Rank', color='red', figsize=(16,12))
        immunogenicity_rank.fill_between(x, y-df_rank_error, y+df_rank_error, color='blue')
        immunogenicity_rank.set_title('Immunogenic density plot ' + str(protein), fontsize=30)
        immunogenicity_rank.set_ylabel('% Rank', fontsize=18)
        immunogenicity_rank.set_xlabel('Index in protein', fontsize=18)
            
        plt.axhline(y=2, color='lightgreen')
        plt.axhline(y=0.5, color='green')
            
        plt.savefig(my_path + 'immunogenic_density_rank-' + str(protein) + '.png')
    
    return


def structure_netCTLpan(filename, path, separator=','):
    ''' 
    Takes a cleaned up netCTLpan file as input.
    The default separator is a ',' (comma).
    The function returns a cleaned up pandas dataframe.
    '''
    
    my_path = path
    
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
            
 
#5: run pipeline()
           
def pipeline(model, path, allele_file, FASTA_file):
    
    my_path = path
    
    with open(path + allele_file, 'r') as alleles:
        alleles = list(map(lambda s: s.strip(), alleles))
    
    
    Parallel(n_jobs=cpu_count())(delayed(predict_epitope_affinity)(i, FASTA_file, path) for i in alleles)
    
    remove_subheaders('AFFINITY.txt', 'subheaders.txt')
    
    remove_empty_lines('subheaders.txt', 'emptylines.txt')
    
    replace_spaces('emptylines.txt', 'affinity_clean.txt')
    
    if model=='netMHCpan':
        structure_netMHCpan('affinity_clean.txt', path)
    else:
        structure_netCTLpan('affinity_clean.txt', path)
        
    return 
        
        
pipeline('netMHCpan', '/home/svalkiers/data_folder/', 'HLA_alleles.txt', 'FASTA_example.txt')
