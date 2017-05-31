#!/usr/bin/python

import sys
import pyfasta
import re
import os
import csv


def main(argv):
    
    # Path to mounted midway on Linux:
    #data_directory = "/run/user/1000/gvfs/smb-share:server=midwaysmb.rcc.uchicago.edu,share=project/cobey/mvieira/Ig_evolvability/Data/Doria-Rose_2014_VRC26/"

    #local directory:
    data_directory = "../../Data/Doria-Rose_2014_VRC26/"


    raw_file_name = data_directory + "Doria-Rose_2014_VRC26.fasta"
    raw_file = open(raw_file_name, 'r')

    processed_file_path = '../../results/alignments/VRC26_plus_GERMLINE.fasta'
    processed_file = open(processed_file_path, 'w')

    # Adding germline sequence:
    # Concatenate IGHV3-30*18 and IGHJ3 (germline according to authors)
    # I'll assume J gene is allele *01 for now.

    # J sequence from IMGT: 
    J_seq = 'TGATGCTTTTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG'

    processed_file.write('>GERMLINE' + '\n')
    # IGHV3-30*18

    processed_file.write('CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAG' + '\n')
    processed_file.write('CCTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTG' + '\n')
    processed_file.write('GGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATC' + '\n')
    processed_file.write('TCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGT' + '\n')
    processed_file.write('ATTACTGTGCGAAAGA')
 

    # IHGJ3*01
    processed_file.write(J_seq + '\n')
    processed_file.write('\n')

    # Adding VRC26 sequences, simplyfing the identifying line ('>')
    for line in raw_file:
        seq_name = re.search(r'KJ[^|]*', line)
        if seq_name is not None:
            seq_name = seq_name.group()
            seq_name = re.search('[^,\.]*', seq_name)
            seq_name = seq_name.group()
            replacement_line = '>' + seq_name + '\n'
            processed_file.write(replacement_line)
        else:
            processed_file.write(line)

    raw_file.close()
    processed_file.close()
    return 0

if(__name__== "__main__"):
    status = main(sys.argv)
    sys.exit(status)
