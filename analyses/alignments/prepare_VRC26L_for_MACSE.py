#!/usr/bin/python

import sys
import pyfasta
import re
import os
import csv


def main(argv):
    
    # Path to mounted midway on Linux:
    #data_directory = "/run/user/1000/gvfs/smb-share:server=midwaysmb.rcc.uchicago.edu,share=project/cobey/mvieira/Ig_evolvability/Data/Doria-Rose_2014_VRC26L/"

    #local directory:
    data_directory = "../../Data/Doria-Rose_2014_VRC26/"


    raw_file_name = data_directory + "Doria-Rose_2014_VRC26L.fasta"
    raw_file = open(raw_file_name, 'r')

    processed_file_path = '../../results/alignments/VRC26L_plus_GERMLINE.fasta'
    processed_file = open(processed_file_path, 'w')

    # Adding germline sequence:
    # IGLV1-51*02 (germline according to authors -- see the paper, Fig. 1b)

    processed_file.write('>GERMLINE' + '\n')
    
    # IGLV1-51*02

    processed_file.write('CAGTCTGTGTTGACGCAGCCGCCCTCAGTGTCTGCGGCCCCAGGACAGAAGGTCACCATCTCCTGCTCTG' + '\n')
    processed_file.write('GAAGCAGCTCCAACATTGGGAATAATTATGTATCCTGGTACCAGCAGCTCCCAGGAACAGCCCCCAAACT' + '\n')
    processed_file.write('CCTCATCTATGAAAATAATAAGCGACCCTCAGGGATTCCTGACCGATTCTCTGGCTCCAAGTCTGGCACG' + '\n')
    processed_file.write('TCAGCCACCCTGGGCATCACCGGACTCCAGACTGGGGACGAGGCCGATTATTACTGCGGAACATGGGATA' + '\n')
    processed_file.write('GCAGCCTGAGTGCTGG')

    processed_file.write('\n\n')

    # Adding VRC26L sequences, simplyfing the identifying line ('>')
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
