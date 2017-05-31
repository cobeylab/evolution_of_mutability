#!/usr/bin/python

import sys
import re
import os
import csv
import pyfasta

def main(argv):

    data_directory = "../../data/Liao_2013_CH103/"

    original_fasta_filename = "Liao_2013_CH103.fasta"

    # Original fasta file with dates in sequence names:
    VDJ_CH103 = pyfasta.Fasta(data_directory + original_fasta_filename)
    CH103_ids = VDJ_CH103.keys()

    # Dictionary with dates for sequences with no sampling time in fasta sequence name (from Liao paper)

    # We couldn't find dates for KC575852 and KC575851 (coded 1CZVFE and 2GKPY3 by the authors). We'll assume that they were sampled at the median time point (53 weeks) 
    missing_times = {'KC575855': '14',
                     'KC575854': '14',
                     'KC575853': '140',
                     'KC575852': '53',
                     'KC575851': '53',
                     'KC575850': '66',
                     'KC575849': '66',
                     'KC575848': '136',
                     'KC575847': '136',
                     'KC575846': '136',
                     'KC575845': '136'}

    # Open original and alingment file
    original_alignment_filename = "../../results/alignments/CH103_MACSE_NT_EDITED.nex"
    original_alignment_file = open(original_alignment_filename, 'r')

    # Opening modified, dated combined nexus file:
    output_filename = "../../results/alignments/CH103_final_alignment.nex"
    output_file = open(output_filename, 'w')

    # Adjust alignment file
    for line in original_alignment_file:
        # Does this line contain a sequence name?
        seq_name = re.search(r'KC[0-9]+', line)
        if seq_name is not None:
            # If it does, go through the sequence names in the fasta
            seq_name = seq_name.group()
            for i in range(len(CH103_ids)):
                # Find the corresponding sequence:
                if CH103_ids[i].find(seq_name) is not -1:
                    time_point = re.search(r'_w[^_]*', CH103_ids[i])
                    # If it has a sampling date:
                    if time_point is not None:
                        # Find it and append to the sequence name with _
                        time_point = time_point.group()
                        time_point = re.search(r'[0-9]+', time_point)
                        time_point = time_point.group()
                    else:
                        time_point = missing_times[seq_name]
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
