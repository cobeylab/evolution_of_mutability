#!/usr/bin/python

import sys
import re
import os
import csv
import pyfasta

def main(argv):

    data_directory = "../../data/Wu_2015_VRC01/"

    original_fasta_filename = "Wu_2015_VRC01_heavy_curated.fasta"

    # Original fasta file with dates in sequence names:
    VDJ_VRC01 = pyfasta.Fasta(data_directory + original_fasta_filename)
    VDJ_VRC01_ids = VDJ_VRC01.keys()

    times = {'O': 60, 'A': 136, 'B': 145,
			 'C': 196, 'D': 202, 'E': 208,
		     'F': 214, 'G': 221, 'H': 231, 'I': 237}

    
    # Open original alignment file
    original_alignment_filename = "../../results/alignments/VRC01_19_MACSE_NT_EDITED.nex"
    original_alignment_file = open(original_alignment_filename, 'r')

    # Opening modified, dated combined nexus file:
    output_filename = "../../results/alignments/VRC01_19_final_alignment.nex"
    output_file = open(output_filename, 'w')

    # Adjust alignment file
    # For each line in the original alignment file:
    for line in original_alignment_file:
        # Does this line contain a sequence name?        
        seq_name = re.search(r'KP[0-9]+', line)
        if seq_name is not None:
            # If it does, go through the seq ids. in the original fasta
            seq_name = seq_name.group()
            for id in VDJ_VRC01_ids:
            	#If this id contains the focal sequence:
            	if id.find(seq_name) is not -1:
            		#Find time point:
            		time_point = re.search(r'\..-', id)
            		time_point = time_point.group()
            		time_point = time_point[1]
            		#Convert letter to number using "times" dictionary
            		time_point = times[time_point]
            	         
            	    #Add time point to modified alignment file    
            		dated_seq_name = seq_name + "_" + str(time_point)
            		dated_line = line.replace(seq_name, dated_seq_name)
            		output_file.write(dated_line) 
        else:
            line = line.replace("GERMLINE", "GERMLINE_00")
            output_file.write(line)

    # output_file.write(tree_string)
    original_alignment_file.close()
    output_file.close()
    return 0

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
