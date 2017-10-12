#!/usr/bin/python

"""Implements a function for computing several mutability metrics for a DNA sequence based on a set of specified sequence partitions
"""

__author__ = "Marcos Vieira (mvieira@uchicago.edu)"


# imports
import csv
import os
import itertools
import numpy as np

# Get location of this script (allows files to be read by this script when calling this script from other scripts)
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


# constants
# Test sequence:

test_sequence = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCGTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
test_partition = [1, 34, 58, 109, 130, 244, 288]


test_WRCH = 'TGCTCCCCTGCT'
test_DGYW = 'GGTTGGGGAGCT'
test_WA = 'TAACCCTA'


test_sequence = test_sequence.upper()

# mutability scores from the S5F model by Yaari et al:
csvfile = open(__location__ + '/Yaari_fivemer_scores.csv', 'r')
csvread = csv.reader(csvfile, delimiter=' ', quotechar='|')
S5F = {}
firstline = True
for row in csvread:
    if firstline:
        firstline = False
    else:
        csvrow = row
        motif = csvrow[0].replace('"', '')
        mutability = csvrow[1].replace('"', '')
        S5F[motif] = float(mutability)

# constants:
W = {'A', 'T'}
Y = {'C', 'T'}
R = {'G', 'A'}
S = {'G', 'C'}
D = {'A', 'G', 'T'}
H = {'T', 'A', 'C'}
K = {'T','G'}
M = {'C','A'}
V = {'A','C','G'}
B = {'C','G','T'}
N = {'A','T','C','G'}
X = {'A','T','C','G'}

ambiguity_code = {'W' : W, 'Y' : Y, 'R' : R, 'S' : S,
                  'D' : D, 'H' : H, 'K' : K, 'M' : M,
                  'V' : V, 'B' : B, 'N' : N, 'X' : X}

# mutability scores from the background indep. component of the  7mer model of Elhanati et al:
csvfile = open(__location__ +'/Elhanati_files/7mer_freq_indep.csv', 'r')
csvread = csv.reader(csvfile, delimiter=',')
SevenMer = {}
for row in csvread:
    csvrow = row
    motif = csvrow[0]
    mutability = csvrow[1]
    SevenMer[motif] = float(mutability)


# Function for computing mean S5F score:

def compute_mean_S5F(seq):
    mean_S5F = 0
    mean_log_S5F = 0

    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}

    # For each nucleotide, except the first two and the last two:
    for i in range(2, (len(seq)-2)):
        motif = seq[i-2:i+3]
    
    # Check if motif contains ambiguous nucleotides
        # If it doesn't, find S5F from dictionary:
        if len(ambiguities.intersection(set(motif))) == 0:
            mean_S5F = mean_S5F + S5F[motif]
            mean_log_S5F = mean_log_S5F + np.log(S5F[motif])

        else:
            # If it does, find the ambiguous positions:
            positions = [p for p in range(len(motif)) if len(ambiguities.intersection(set(motif[p]))) > 0]

            # Find the exact ambiguity at each site:
            ambiguous_nts = [motif[p] for p in positions]

            # Find the sets of possible nucleotides for each site:
            alternative_nt_sets = [ambiguity_code[nt] for nt in ambiguous_nts]

            # Find all possible choices of alternative nucleotides for each position:
            replacements = list(itertools.product(*alternative_nt_sets))

            set_lengths = [len(nt_set) for nt_set in alternative_nt_sets]
            n_replacements = np.prod(np.array(set_lengths))

            assert len(replacements) == n_replacements

            S5F_values = []
            S5F_log_values = []
            for replacement in replacements:
                possible_motif = list(motif)
                for i in range(len(replacement)):
                    possible_motif[positions[i]] = replacement[i]
                possible_motif = ''.join(possible_motif)
                S5F_values.append(S5F[possible_motif])
                S5F_log_values.append(np.log(S5F[possible_motif]))

            mean_value = float(sum(S5F_values)) / len(S5F_values)
            mean_log_value = float(sum(S5F_log_values)) / len(S5F_log_values)

            mean_S5F = mean_S5F + mean_value
            mean_log_S5F = mean_log_S5F + mean_log_value

    mean_S5F = mean_S5F / float(len(seq) - 4)
    mean_log_S5F = mean_log_S5F / float(len(seq) - 4)

    # Compute geometric mean by exponentiating mean log
    geom_mean_S5F = np.exp(mean_log_S5F)
    return {'mean_S5F': mean_S5F, 'geom_mean_S5F': geom_mean_S5F}

# For a single motif, mean and geometric mean should be the same:
assert(compute_mean_S5F('CATAG')['mean_S5F'] == compute_mean_S5F('CATAG')['geom_mean_S5F'])



# Manual calculation for an example sequence
test_sum_log = np.log(S5F['AAGAA']) + np.log(S5F['AGAAA']) + np.log(S5F['GAAAA']) + np.log(S5F['AAAAG']) + \
               np.log(S5F['AAAGA']) + np.log(S5F['AAGAA'])
test_mean_log = test_sum_log / 6
test_geom_mean = np.exp(test_mean_log)

# Test that sequence computes the correct geometric mean:
assert(compute_mean_S5F('AAGAAAAGAA')['geom_mean_S5F'] == test_geom_mean)

# Function for computing mean 7mer score:

def compute_mean_7M(seq):
    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}
    
    mean_7M = 0

    # For each nucleotide, except the first three and the last three:
    for i in range(3, (len(seq)-3)):
        motif = seq[i-3:i+4]
    
        # Check if motif contains ambiguous nucleotides
        # If it doesn't, find 7-mer mutability from dictionary:
        if len(ambiguities.intersection(set(motif))) == 0:
            mean_7M = mean_7M + SevenMer[motif]
        else:
            # If it does, find the ambiguous positions:
            positions = [p for p in range(len(motif)) if len(ambiguities.intersection(set(motif[p]))) > 0]

            # Find the exact ambiguity at each site:
            ambiguous_nts = [motif[p] for p in positions]

            # Find the sets of possible nucleotides for each site:
            alternative_nt_sets = [ambiguity_code[nt] for nt in ambiguous_nts]

            # Find all possible choices of alternative nucleotides for each position:
            replacements = list(itertools.product(*alternative_nt_sets))

            set_lengths = [len(nt_set) for nt_set in alternative_nt_sets]
            n_replacements = np.prod(np.array(set_lengths))

            assert len(replacements) == n_replacements
            
            SevenMer_values = []
            for replacement in replacements:
                possible_motif = list(motif)
                for i in range(len(replacement)):
                    possible_motif[positions[i]] = replacement[i]
                possible_motif = ''.join(possible_motif)
                SevenMer_values.append(SevenMer[possible_motif])

            mean_value = float(sum(SevenMer_values)) / len(SevenMer_values)
            mean_7M = mean_7M + mean_value

    mean_7M = mean_7M / float(len(seq) - 6)
    return mean_7M

# Function for -counting- 7mers with non-zero mutability
def count_nonzero_7M(seq):
    # Set of ambiguities:
    ambiguities = {'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N'}

    nonzero_7M = 0
    # For each nucleotide, except the first three and the last three:
    for i in range(3, (len(seq)-3)):
        # If motif is non-ambiguous
        if len(ambiguities.intersection(set(seq[i-3:i+4]))) == 0:
            if SevenMer[seq[i-3:i+4]] > 0:
                nonzero_7M = nonzero_7M + 1
    return nonzero_7M


# Functions for counting AID hotspots:

def count_WRCH(seq):
    n_WRCH = 0
    # For each nucleotide, except the first two and the last one...
    for i in range(2, (len(seq)-1)):
        if (seq[i-2] in W) and (seq[i-1] in R) and (seq[i] == 'C') and (seq[i+1] in H):
            n_WRCH += 1
    return n_WRCH


def count_DGYW(seq):
    n_DGYW = 0
    # For each nucleotide, except the first and the last two...
    for i in range(1, (len(seq)-2)):
        if (seq[i-1] in D) and (seq[i] == 'G') and (seq[i+1] in Y) and (seq[i+2] in W):
            n_DGYW += 1
    return n_DGYW

# Functions for counting POL hotspots:


def count_WA(seq):
    n_WA = 0
    # For each nucleotide, except the first one:
    for i in range(1, len(seq)):
        if (seq[i-1] in W) and (seq[i] == 'A'):
            n_WA += 1
    return n_WA


def count_TW(seq):
    n_TA = 0
    # For each nucleotide, except the last one:
    for i in range(0, (len(seq) - 1)):
        if (seq[i] == 'T') and (seq[i + 1] in W):
            n_TA += 1
    return n_TA


# Functions for counting AID coldspots:
def count_SYC(seq):
    n_SYC = 0
    # For each nucleotide, except the first two
    for i in range(2, len(seq)):
        if (seq[i-2] in S) and (seq[i-1] in Y) and (seq[i] == "C"):
            n_SYC += 1
    return n_SYC


def count_GRS(seq):
    n_GRS = 0
    # For each nucleotide, except the last two:
    for i in range(0, (len(seq) - 2)):
        if (seq[i] == "G") and (seq[i+1] in R) and (seq[i+2] in S):
            n_GRS += 1
    return n_GRS


def count_trinucleotide_CS_coding(seq):
    n_trinucleotide_CS_coding = 0

    # Coding strand coldspot trinucleotides
    coding_CS = {"TTC", "CAC", "GGC", "GAC"}

    # For each nucleotide, except the first two:
    for i in range(2, len(seq)):
        tri = seq[i-2:i+1]
        if (tri in coding_CS):
            n_trinucleotide_CS_coding += 1
    return n_trinucleotide_CS_coding


def count_trinucleotide_CS_noncoding(seq):
    n_trinucleotide_CS_noncoding = 0
    # non-coding strand coldspot trinucleotides
    noncoding_CS = {"GAA", "GTG", "GCC", "GTC"}
    # For each nucleotide, except the last two:
    for i in range(0, (len(seq) - 2)):
        tri = seq[i:i+3]
        if (tri in noncoding_CS):
            n_trinucleotide_CS_noncoding += 1
    return n_trinucleotide_CS_noncoding
    

# Function for computing overlapping hotspots (Wei et al. 2015 PNAS):
def count_OHS(seq):
    n_OHS = 0
    # For each nucleotide, except the first two and the last one:
    for i in range(2, (len(seq)-1)):
        if (seq[i-2] + seq[i-1] + seq[i] + seq[i+1]) == 'AGCT':
            n_OHS += 1
    return n_OHS
    

# Main functions:
def seq_mutability(seq, partition_points=None):
    """Calculates several mutability metrics for a specified set of contiguous partitions in a DNA sequence (seq). The partitions are specified by argument partition_points. Each element of the list partition_points specifies the position in the sequence (NOT the Python index) where a partition begins. The last element indicates the last position of the last partition (in case the user does not want to cover the entire sequence). If partition_points is set to 'None', the whole sequence is treated as a single partition. Deals with gaps.
    """
    seq = seq.upper()
    if partition_points is None:
        partition_points = [1, len(seq)]

    # Convert partition_points into Python indices (i.e. subtract 1):
    partition_points = [(partition_points[i] - 1)
                        for i in range(len(partition_points))]

    # Output list of lists.First internal list gives partition lengths.
    # Each subsequent list gives mutability scores for a partition
    results_lists = [[]]

    # From partition points, store ungapped partitions in a list:
    ungapped_partitions = []

    for i in (range(len(partition_points) - 1)):
        if i < (len(partition_points) - 2):
            partition = seq[partition_points[i]: partition_points[i+1]]
        else:
            # For the last partition, the correct parsing expression is
            partition = seq[partition_points[i]: partition_points[i+1] + 1]
            # Since the last point is the last position of the las partition.
        ungapped_partitions.append(partition.replace("-", ""))
        results_lists[0].append(len(partition.replace("-", "")))

    # Find and ungap right and left "margins" - regions outside partitions
    # If 1st p. point not at the beggining of seq, then there's a left margin
    if partition_points[0] > 0:
        left_margin = seq[0:partition_points[0]]
        left_margin = left_margin.replace("-", "")
    else:
        left_margin = ""

    # If last p. point not the last position of seq, there's a right margin:
    if partition_points[len(partition_points) - 1] < (len(seq) - 1):
        right_margin = seq[(partition_points[len(partition_points) - 1]) + 1:(len(seq) - 1)]
        right_margin = right_margin.replace("-", "")
    else:
        right_margin = ""

    # For each (ungapped) partition:
    for i in range(len(ungapped_partitions)):
        partition = ungapped_partitions[i]

        # Find neighbor nucleotides
        # If this is the leftmost partition
        if i == 0:
            if len(left_margin) > 0:
                l_neighbors = ['', '',  left_margin[-1]]
            if len(left_margin) > 1:
                l_neighbors = ['', left_margin[-2], left_margin[-1]]
            if len(left_margin) > 2:
                l_neighbors = [left_margin[-3], left_margin[-2], left_margin[-1]]
            if len(left_margin) == 0:
                l_neighbors = ['','','']
            # If there are other partitions:
            if len(ungapped_partitions) > 1:
                r_neighbors = [ungapped_partitions[i+1][0], ungapped_partitions[i+1][1],
                               ungapped_partitions[i+1][2]]
            else:
                r_neighbors = ['','','']

        # If this is the rightmost partition:
        elif i == (len(ungapped_partitions) - 1):
            if len(right_margin) > 0:
                r_neighbors = [right_margin[0], '', '']
            if len(right_margin) > 1:
                r_neighbors = [right_margin[0], right_margin[1], '']
            if len(right_margin) > 2:
                r_neighbors = [right_margin[0], right_margin[1], right_margin[2]]
            if len(right_margin) == 0:
                r_neighbors = ['','','']
            # If there are other partitions:
            if len(ungapped_partitions) > 1:
                l_neighbors = [ungapped_partitions[i-1][-3], ungapped_partitions[i-1][-2],
                               ungapped_partitions[i-1][-1]]
            else:
                l_neighbors = ['','','']
        # If this partition has neighbor partitions on both sides:
        else:
            l_neighbors = [ungapped_partitions[i-1][-3], ungapped_partitions[i-1][-2],
                           ungapped_partitions[i-1][-1]]
            r_neighbors = [ungapped_partitions[i+1][0], ungapped_partitions[i+1][1],
                           ungapped_partitions[i+1][2]]

        # Calc. mutability metrics, appending neighbors where necessary:
        mean_S5F = compute_mean_S5F(l_neighbors[1] + l_neighbors[2] +
                                    partition +
                                    r_neighbors[0] + r_neighbors[1])['mean_S5F']

        geom_mean_S5F = compute_mean_S5F(l_neighbors[1] + l_neighbors[2] +
                                    partition +
                                    r_neighbors[0] + r_neighbors[1])['geom_mean_S5F']

        mean_7M = compute_mean_7M(l_neighbors[0] + l_neighbors[1] + l_neighbors[2] +
                                  partition +
                                  r_neighbors[0] + r_neighbors[1] + r_neighbors[2])

        n_nonzero_7M = count_nonzero_7M(l_neighbors[0] + l_neighbors[1] + l_neighbors[2] +
                                  partition +
                                  r_neighbors[0] + r_neighbors[1] + r_neighbors[2])

        n_WRCH = count_WRCH(l_neighbors[1] + l_neighbors[2] +
                            partition + r_neighbors[0])

        n_DGYW = count_DGYW(l_neighbors[2] + partition +
                            r_neighbors[0] + r_neighbors[1])

        n_AIDHS = n_WRCH + n_DGYW

        n_WA = count_WA(l_neighbors[2] + partition)

        n_TW = count_TW(partition + r_neighbors[0])

        n_POLHS = n_WA + n_TW

        C_fraction = partition.count('C')
        C_fraction = C_fraction / float(len(partition))

        n_SYC = count_SYC(l_neighbors[1] + l_neighbors[2] + partition)
        n_GRS = count_GRS(partition + r_neighbors[0] + r_neighbors[1])
        trinuc_CS_coding = count_trinucleotide_CS_coding(l_neighbors[1] +
                                                         l_neighbors[2] +
                                                         partition)

        trinuc_CS_noncoding = count_trinucleotide_CS_noncoding(partition +
                                                               r_neighbors[0] +
                                                               r_neighbors[1])

        n_AIDCS = n_SYC + n_GRS + trinuc_CS_coding + trinuc_CS_noncoding
        
        n_OHS = count_OHS(l_neighbors[1] + l_neighbors[2] +
                            partition + r_neighbors[0])

        partition_results = {'mean_S5F':mean_S5F, 'mean_7M':mean_7M, 'n_nonzero_7M':n_nonzero_7M,
                             'HS':n_AIDHS, 'POLHS':n_POLHS, 'C_fraction':C_fraction, 'CS':n_AIDCS,
                             'OHS':n_OHS, 'geom_mean_S5F' : geom_mean_S5F}

        results_lists.append(partition_results)
    return results_lists
    
    
def aggregated_mutability(seq, partition_points):
    '''Calls seq_mutability and aggregates results into FRs (partitions 1,3,5) vs CDRs (partitions 2,4,6)
       Hotspots, coldspots, polymerase hotspots, overlapping HS and non-zero 7Ms are summed, others are averaged weighted by the length of each partition"

        >>> aggregated_mutability(seq = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCGTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA', partition_points = [1, 34, 58, 109, 130, 244, 288])
        {'FR_mutability': {'C_fraction': 0.30808080808080807, 'POLHS': 23, 'HS': 19, 'geom_mean_S5F': 0.56619255059234308, 'n_nonzero_7M': 195, 'mean_S5F': 0.8343007913008704, 'CS': 65, 'mean_7M': 1.0339051547200506, 'OHS': 3}, 'CDR_mutability': {'C_fraction': 0.25555555555555554, 'POLHS': 23, 'HS': 14, 'geom_mean_S5F': 0.72673261557134272, 'n_nonzero_7M': 90, 'mean_S5F': 1.121063501136274, 'CS': 24, 'mean_7M': 1.2096232738333335, 'OHS': 0}}
    '''
    mutability = seq_mutability(seq, partition_points)

    # Dictionary with lengths of each region
    lengths = {'FR':[mutability[0][i] for i in [0,2,4]],
               'CDR':[mutability[0][i] for i in [1,3,5]]
               }

    # Get indices for FR or CDR results in mutability[i] (note mutability[0] gives region lengths)
    indices = {'FR': [1,3,5], 'CDR': [2,4,6]}

    results = {'FR_mutability' : {}, 'CDR_mutability' : {}}

    for metric in mutability[1].keys():
        for region in ['FR', 'CDR']:
            region_lengths = lengths[region]
            mutability_values = [mutability[i][metric] for i in indices[region]]

            # If metric is number of hotspots, coldspots, pol hotspots, OHS or non-zero 7ms, sum across regions
            if metric in ['HS','CS','n_nonzero_7M','POLHS','OHS']:
                aggregated_value = sum(mutability_values)

            # Else, do a weighted average
            # For all metrics except the geometric mean of S5F scores:
            elif metric != 'geom_mean_S5F':
                # Find weight of each FR (or each CDR) relative to all FRs (or all CDRs)
                region_weight = [float(length) / sum(region_lengths) for length in region_lengths]

                aggregated_value = [mutability_values[i] * region_weight[i] for i in range(len(mutability_values))]
                aggregated_value = sum(aggregated_value)

            # If metric is the geometric mean of S5F:
            elif metric == 'geom_mean_S5F':
                # Find weight of each FR (or each CDR) relative to all FRs (or all CDRs)
                region_weight = [float(length) / sum(region_lengths) for length in region_lengths]

                # Get aggregated geom. mean as the product of geom means for each region to the power of their weight

                aggregated_value = [mutability_values[i] ** region_weight[i] for i in range(len(mutability_values))]
                aggregated_value = np.prod(aggregated_value)

            results[region + '_mutability'][metric] = aggregated_value

    return {'FR_mutability': results['FR_mutability'], 'CDR_mutability': results['CDR_mutability']}