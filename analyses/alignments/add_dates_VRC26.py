#!/usr/bin/python

import sys
import re
import os
import csv
import pyfasta

def main(argv):

    data_directory = "../../data/Doria-Rose_2014_VRC26/"

    original_fasta_filename = "Doria-Rose_2014_VRC26.fasta"

    # Original fasta file with dates in sequence names:
    VDJ_VRC26 = pyfasta.Fasta(data_directory + original_fasta_filename)
    VRC26_ids = VDJ_VRC26.keys()

    # Open original and alingment file
    original_alignment_filename = "../../results/alignments/VRC26_MACSE_NT_EDITED.nex"
    original_alignment_file = open(original_alignment_filename, 'r')

    # Opening modified, dated combined nexus file:
    output_filename = "../../results/alignments/VRC26_final_alignment.nex"
    output_file = open(output_filename, 'w')

    # Adjust alignment file
    for line in original_alignment_file:
        # Does this line contain a sequence name?
        seq_name = re.search(r'KJ[0-9]+', line)
        if seq_name is not None:
            # If it does, go through the sequence names in the fasta
            seq_name = seq_name.group()
            for i in range(len(VRC26_ids)):
                # Find the corresponding sequence:
                if VRC26_ids[i].find(seq_name) is not -1:
                    # Find the sampling date
                    time_point = re.search(r'cap256-[^\s]*', VRC26_ids[i])
                    time_point = time_point.group()
                    time_point = re.search(r'-.*', time_point)
                    time_point = time_point.group()
                    time_point = time_point[1:len(time_point)]
                    time_point = re.search(r'[^\-]*', time_point)
                    time_point = time_point.group()

                    # Append it to the sequence name with _

                    dated_seq_name = seq_name + "_" + time_point 
                    dated_line = line.replace(seq_name, dated_seq_name)
                    output_file.write(dated_line)
        else:
            line = line.replace("GERMLINE", "GERMLINE_00")
            output_file.write(line)

    original_alignment_file.close()
    output_file.close()
    return 0

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
