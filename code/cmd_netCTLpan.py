# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:59:11 2019

@author: sebas
"""

import subprocess, shlex
from subprocess import Popen, PIPE

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

execute_netCTLpan('/home/svalkiers/data_folder/HLA_alleles.txt', 'FASTA_example.txt')
