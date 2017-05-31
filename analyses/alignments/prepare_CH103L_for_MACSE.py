#!/usr/bin/python

import sys
import pyfasta
import re
import os
import csv


def main(argv):

    # Path to mounted midway on Linux:
    #data_directory = "/run/user/1000/gvfs/smb-share:server=midwaysmb.rcc.uchicago.edu,share=project/cobey/mvieira/Ig_evolvability/Data/Liao_2013_CH103L/"

    #local directory:
    data_directory = "../../Data/Liao_2013_CH103/"


    raw_file_name = data_directory + "Liao_2013_CH103L.fasta"
    raw_file = open(raw_file_name, 'r')

    processed_file_path = '../../results/alignments/CH103L_plus_GERMLINE.fasta'
    processed_file = open(processed_file_path, 'w')

    # Adding germline sequence:
    # Concatenate IGLV3-1 and IGLJ1 (germline according to authors)
    # I'll assume both are allele *01

    # J sequence from IMGT: 
    # J_seq = 'TTATGTCTTCGGAACTGGGACCAAGGTCACCGTCCTAG'

    processed_file.write('>GERMLINE' + '\n')
    
    # IGLV3-1*01     
    processed_file.write('TCCTATGAGCTGACTCAGCCACCCTCAGTGTCCGTGTCCCCAGGACAGACAGCCAGCATCACCTGCTCTG' + '\n')
    processed_file.write('GAGATAAATTGGGGGATAAATATGCTTGCTGGTATCAGCAGAAGCCAGGCCAGTCCCCTGTGCTGGTCAT' + '\n')
    processed_file.write('CTATCAAGATAGCAAGCGGCCCTCAGGGATCCCTGAGCGATTCTCTGGCTCCAACTCTGGGAACACAGCC' + '\n')
    processed_file.write('ACTCTGACCATCAGCGGGACCCAGGCTATGGATGAGGCTGACTATTACTGTCAGGCGTGGGACAGCAGCA' + '\n')
    processed_file.write('CTGCA')


    # IHGJ4*01
    #processed_file.write(J_seq + '\n')
    processed_file.write('\n\n')

    # Adding CH103L sequences, simplyfing the identifying line ('>')
    for line in raw_file:
        seq_name = re.search(r'KC[^|]*', line)
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

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
