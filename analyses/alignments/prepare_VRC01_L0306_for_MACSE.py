#!/usr/bin/python

import sys
import pyfasta
import re
import os
import csv


def main(argv):

    #sequence directory:
    data_directory = "../../Data/Wu_2015_VRC01/"
    
    raw_file_name = data_directory + "VRC01_L0306.fasta"
    raw_file = open(raw_file_name, 'r')

    processed_file_path = '../../results/alignments/VRC01_L0306_plus_GERMLINE.fasta'
    processed_file = open(processed_file_path, 'w')

    # Adding germline sequence:
    # IGKV3-20*01 (germline according to authors -- see the paper, Fig. 1B)

    processed_file.write('>GERMLINE' + '\n')
    
    # IGKV3-20*01
                          
    processed_file.write('GAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCA'+ '\n')
    processed_file.write('GGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCT'+ '\n')
    processed_file.write('CCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACA'+ '\n')
    processed_file.write('GACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAGTATGGTA'+ '\n')
    processed_file.write('GCTCACCTCC')
    
    processed_file.write('\n\n')
    
    # Adding VRCO1_L0306 sequences, simplyfing the identifying line ('>')
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
