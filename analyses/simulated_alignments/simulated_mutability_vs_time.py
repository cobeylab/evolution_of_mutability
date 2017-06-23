#!/usr/bin/python
'''Simple script to compute the mutability of sequences sampled at different time points during simulated B cell evolution.
    Used to check that simulations under scenarios 5 & 6 show changes in mutability for the specified sampling times and time breaks.
'''
import re
import sys
#import csv
import pyfasta
from Bio import AlignIO

# Import mutability functions from analyses/mutability folder
sys.path.insert(0, '../mutability/')
from mutability_function import compute_mean_S5F, count_WRCH, count_DGYW


def main(argv):
    alignment_file = str(argv[1])
    
    alignment = AlignIO.read(alignment_file, "nexus")
    
    time_point = []
    mean_S5F = []
    n_HS = []
    
    for sequence in alignment:
        # Read sequence name (which includes time point after _)
        name = sequence.name
        sampling_time = int(re.search('_[0-9]*', name).group().replace('_',''))

        # Append time point to list of sampling time points
        time_point.append(sampling_time)
        
        # append mean S5F score
        mean_S5F.append(compute_mean_S5F(str(sequence.seq)))
        
        # append number of hotspots
        n_HS.append(count_WRCH(str(sequence.seq)) + count_DGYW(str(sequence.seq)))
        
    output_filepath = alignment_file.replace('.nex', '_mutability_per_timepoint.csv')

    with open(output_filepath, 'w') as output_file:
        output_file.write('time_point,mean_S5F,n_HS\n')
        
        for i in range(len(time_point)):
            output_file.write(str(time_point[i]) + ',' + str(mean_S5F[i]) + ',' + str(n_HS[i]) + '\n')

if(__name__== "__main__"):
    status = main(sys.argv)
    sys.exit(status)    
        
            
    

    
    