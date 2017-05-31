#!/usr/bin/python

import sys
import pyfasta
import re
import os
import csv


def main(argv):
    #local directory:
    data_directory = "../../results/0_Clone_assignment/VRC01/"
    
    raw_file_name = data_directory + "VRC01_01.fa"
    raw_file = open(raw_file_name, 'r')

    processed_file_path = '../../results/alignments/VRC01_01_plus_GERMLINE.fasta'
    processed_file = open(processed_file_path, 'w')

    # Adding germline sequence:
    # Concatenate IGHV-ORF15-1*04 and IGHJ4*03 (germline according to partis)

    # J sequence from IMGT: 
    J_seq = 'GCTACTTTGACTACTGGGGCCAAGGGACCCTGGTCACCGTCTCCTCAG'

    processed_file.write('>GERMLINE' + '\n')
    
    # IGHV-ORF15-1*04
    processed_file.write('CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGG'+ '\n')
    processed_file.write('CTTCTGGATACATCTTCACCGACTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGAGCTTGGGTG'+ '\n')
    processed_file.write('GATGGGACGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGAGTCACCATG'+ '\n')
    processed_file.write('ACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCACGT'+ '\n')
    processed_file.write('ATTACTGTGCGAGAGA')
    
    
    # IHGJ4*03
    processed_file.write(J_seq + '\n')
    processed_file.write('\n')

    # Adding VRCO1 sequences, simplyfing the identifying line ('>')
    for line in raw_file:
        seq_name = re.search(r'KP[^|]*', line)
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
