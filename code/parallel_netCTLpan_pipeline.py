# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:59:35 2019

@author: sebas
"""

from multiprocessing import cpu_count
from joblib import Parallel, delayed
from itertools import product
from time import sleep
import subprocess, shlex
from subprocess import Popen, PIPE


my_path = '/home/svalkiers/data_folder/'

with open(my_path + 'HLA_alleles.txt', 'r') as alleles:
    alleles = list(map(lambda s: s.strip(), alleles))


with open('/home/svalkiers/data_folder/MHC_test_output.txt', 'w') as f:
    pass
    

def execute_netCTLpan(allele, protein, length=9):
    
    program = "netCTLpan" + " -a " + str(allele) + " -f " + str(protein) + " -l " + str(length)

    p = Popen([program], stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)     # enter commands in terminal
    output, stderr = p.communicate()
    
    path_to_output = '/home/svalkiers/data_folder/MHC_test_output.txt'
    with open(path_to_output, 'a') as out:
        out.write(str(output).replace("\\n", "\n"))

    return


results = Parallel(n_jobs=cpu_count())(delayed(execute_netCTLpan)(i, 'FASTA_example.txt') for i in alleles)

print(results)
