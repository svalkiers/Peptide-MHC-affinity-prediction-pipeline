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


# 2: Parallel exe version (~ number of cores) of netCTLpan

my_path = '/home/svalkiers/data_folder/'

with open(my_path + 'HLA_alleles.txt', 'r') as alleles:
    alleles = list(map(lambda s: s.strip(), alleles))
    

with open(my_path + 'FASTA_example.txt', 'r') as f:
    FASTA = f.read()
    
FASTA_list = []
    
for protein in FASTA.split('>'):
    FASTA_list.append(protein)
    
prot_list = FASTA_list[1:]


with open('/home/svalkiers/data_folder/MHC_test_output.txt', 'w') as f:
    pass


def execute_netCTLpan(allele, protein, output_suffix, length=9):
    
    program = "netCTLpan" + " -a " + str(allele) + " -f " + str(protein) + " -l " + str(length)

    p = Popen([program], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)     # enter commands in terminal
    output, stderr = p.communicate()
    
    path_to_output = '/home/svalkiers/data_folder/MHC_test_output_' + output_suffix + '.txt'
    with open(path_to_output, 'a') as out:
        out.write(str(output).replace("\\n", "\n"))

    return

# Parallellize calculations over the number of cores and write output to one file
# EDIT: this code works

results = Parallel(n_jobs=cpu_count())(delayed(execute_netCTLpan)(i, 'FASTA_example.txt') for i in alleles)


# 3: Prepare output for analysis

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
            
remove_subheaders('MHC_test_output.txt', 'MHC_headers.txt')
remove_empty_lines('MHC_headers.txt', 'MHC_emptylines.txt')
replace_spaces('MHC_emptylines.txt', 'MHC_clean.txt')
