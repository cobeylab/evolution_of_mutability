#!/usr/bin/python
"""Simulates sequence evolution based on S5F mutation rates and transitions for two sequences: a low-mutability sequence created by concatenating the lowest-mutability motif, and a high-mutability sequence created by concatenating the highest mutability motif

"""
import sys
import numpy as np
import csv

# Import mutability functions analyses/mutability folder
sys.path.insert(0, '../mutability/')
from mutability_function import S5F, compute_mean_S5F

# Import class BCR
sys.path.insert(0, '../simulated_alignments/modules')
from classBCR import BCR

# Import conditional transition probabilities for each five-mer (Yaari et al.):
S5F_transitions = {}
with open('../simulated_alignments/modules/Yaari_transition_probabilities.csv', 'r') as csvfile:
    csvread = csv.reader(csvfile, delimiter=' ', quotechar='|')
    csvread.next()
    for row in csvread:
        csvrow = row
        motif = csvrow[0].replace('"', '')
        S5F_transitions[motif] = {}
        S5F_transitions[motif]["A"] = float(csvrow[1].replace('"', ''))
        S5F_transitions[motif]["C"] = float(csvrow[2].replace('"', ''))
        S5F_transitions[motif]["G"] = float(csvrow[3].replace('"', ''))
        S5F_transitions[motif]["T"] = float(csvrow[4].replace('"', ''))

# Function for site specific transition probabilities
def site_specific_transitions(seq):
    '''Returns a list of transition probabilities for each site. Assigns transition probabilities based on five-mer centered at the site'''
    # Sites for which five-mers can't be determined (edges) are given uniform transition probs
    transitions_list = [{}] * len(seq)
    for site in range(len(seq)):
        # If site can be assigned a fivemer motif
        if site in range(2, (len(seq) - 2)):
            site_transitions = S5F_transitions[seq[site - 2:site + 3]]
        else:
            # Assign uniform transition probabilities to the first 2 and last 2 sites
            # (For the hybrid case w HS-based rates and S5F-based transitions)
            initial_nt = seq[site]
            site_transitions = {}
            for base in ['A', 'T', 'C', 'G']:
                if base is initial_nt:
                    site_transitions[base] = 0
                else:
                    site_transitions[base] = 1. / 3

        transitions_list[site] = site_transitions

    # As a test, check that the probability of a site mutating to the current nucleotide is zero for all sites:
    check = [transitions_list[i][seq[i]] == 0 for i in range(len(seq)) if len(transitions_list[i].keys()) > 0]
    assert check.count(False) == 0, 'Error in site specific transitions: non-zero probability of mutating to the current nucleotide'
    return transitions_list


# Function for site specific mutation rates
def site_specific_rates(seq, mean_reference_mutability, reference_mutation_rate):
    '''Finds the mutability of each site using the S5F mutability model. From site-specific mutabilities, the mean mutability of the reference sequence and the reference mutation rate, return a list of site-specific mutation rates (see equation for mu_i,j in model description). Note that h, the mean mutability of the reference sequence and the reference mutation rate are not input as arguments and are instead read from global names. h and the reference mutation rate are read from the control file, and the mean mutability of the reference sequence is calculated at the beginning of the simulation.
    '''
    # Sites at the edges (for which S5F can't be computed) will have mutability 0
    site_specific_mutabilities = [0] * len(seq)
    for site in range(2, (len(seq) - 2)):
        site_specific_mutabilities[site] = S5F[seq[site - 2:site + 3]]

    # Equation for mu_i,j in model description:
    rates = [(float(m) / mean_reference_mutability) * reference_mutation_rate for m in site_specific_mutabilities]

    return (rates)

# Number of times a motif is concatenated to create each sequence
n_repeats = 30

# Number of generations to simulate
nsteps = 5000

# Number of replicate trajectories starting from each sequence
n_trajectories = 10

# Find motifs with the lowest mutability and the highest mutability
lowest_mutability_index = np.array(S5F.values()).argsort()[0]
lowest_mutability_motif = [motif for motif in S5F.keys() if S5F[motif] == S5F.values()[lowest_mutability_index]][0]

highest_mutability_index = np.array(S5F.values()).argsort()[-1]
highest_mutability_motif = [motif for motif in S5F.keys() if S5F[motif] == S5F.values()[highest_mutability_index]][0]

# Initialize sequence consisting of concatenation of the lowest mutability motif:
initial_seq_low = lowest_mutability_motif * n_repeats

# Initialize sequence consisting of concatenation of the highest mutability motif:
initial_seq_high = highest_mutability_motif * n_repeats

# Reference mutation rate
reference_mutation_rate = float(1) / (4 * len(initial_seq_low))

# Reference mutability (i.e. mutability for which the mutation rate equals the reference mutation rate)
# Chosen to be the the mutability of the low-mutability seq
mean_reference_mutability = compute_mean_S5F(initial_seq_low)

# compute_mean_S5F excludes the 2 nucleotides at each edge from the mean
# adjusting value to take them into account with mutability 0 (in line with the site_specific_rates function)
mean_reference_mutability = mean_reference_mutability * float(len(initial_seq_low) - 4) / len(initial_seq_low)

#Initialize BCR objects
sequence_low = BCR(initial_seq_low)
sequence_high = BCR(initial_seq_high)


with open('../../results/equilibrium_mutability/equilibrium_mutability.csv', 'w') as output_file:
    output_file.write('generation,initial_mutability,replicate_trajectory,S5F_mutability\n')

    for trajectory in range(n_trajectories):

        mutability_low = [compute_mean_S5F(initial_seq_low)]
        mutability_high = [compute_mean_S5F(initial_seq_high)]

        for i in range(nsteps):
            # Simulate new sequence in the low mutability series
            site_rates_low = site_specific_rates(sequence_low.sequence, mean_reference_mutability = mean_reference_mutability,
                                                 reference_mutation_rate = reference_mutation_rate)
            site_transitions_low = site_specific_transitions(sequence_low.sequence)
            sequence_low = sequence_low.reproduce(site_mutation_probs = site_rates_low, site_transition_probs = site_transitions_low)

            mutability_low.append(compute_mean_S5F(sequence_low.sequence))

            # Simulate new sequence in the high mutability series
            site_rates_high = site_specific_rates(sequence_high.sequence, mean_reference_mutability=mean_reference_mutability,
                                                 reference_mutation_rate=reference_mutation_rate)
            site_transitions_high = site_specific_transitions(sequence_high.sequence)
            sequence_high = sequence_high.reproduce(site_mutation_probs=site_rates_high, site_transition_probs=site_transitions_high)

            mutability_high.append(compute_mean_S5F(sequence_high.sequence))

        for i in range(nsteps):
            output_file.write(str(i) + ',low,' + str(trajectory) + ',' + str(mutability_low[i]) + '\n')
        for i in range(nsteps):
            output_file.write(str(i) + ',high,' + str(trajectory) + ',' + str(mutability_high[i]) + '\n')
