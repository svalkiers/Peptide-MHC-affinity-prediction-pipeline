# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 09:27:42 2019

@author: sebas
"""

import re

##############
#INSTRUCTIONS#
##############

"""

Call functions by entering the desired input-file and the
desired name of your output-file. The output-file of the
first function can be used as an input for the next. 

"""
    
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



########################
#RUN THIS PIECE OF CODE#
###########|############
          #V#

remove_subheaders('MHC_test_output.txt', 'mytestoutput.txt')
remove_empty_lines('mytestoutput.txt', 'mytestresult.txt')
replace_spaces('mytestresult.txt', 'MHC_clean.txt')