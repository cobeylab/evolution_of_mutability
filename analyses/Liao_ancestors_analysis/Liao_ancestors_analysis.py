#!/usr/bin/python

"""
Looks at ancestral CH103 sequences inferred by Liao et al. (2013, Nature), and partitions changes in mutability for ancestral-descendant pairs.
"""
import sys
from Bio import SeqIO

# Partition points obtained from IMGT
partition_points = [1,76,100,151,172,286,330]

# Import function for partitioning mutability changes into Syn. and Non-Syn. mutations
sys.path.insert(0, '../S_NS_mutability_changes/')
from mutation_functions import sequence_differences

# Read Liao ancestors sequences:
ancestors = {}
with open('Liao_experimentally_characterized_seqs.fasta', 'rU') as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        ancestors[record.id] = str(record.seq)

# Replace '-' with state at UCA
for anc in ancestors.keys():
    seq = ancestors[anc]
    seq = list(seq)
    seq = [ancestors['UCA_VH'][i] if seq[i] == '-' else seq[i] for i in range(len(seq))]
    seq = ''.join(seq)
    ancestors[anc] = seq



# Ancestral descendant pairs from the tree in the paper. [ancestor, descendant]
ancestral_descendant_pairs = [['UCA_VH', 'I8_VH'],
                              ['I8_VH', 'I7_VH'],
                              ['I8_VH', 'I4_VH'],
                              ['I4_VH', 'I3_VH'],
                              ['I3_VH', 'I2_VH'],
                              ['I2_VH', 'I1_VH'],
                              ['I7_VH', '1AZCETI5_VH'],
                              ['I7_VH', '1A102RI6_VH'],
                              ['I4_VH', '1AH92U_VH'],
                              ['I3_VH', 'CH105_VH'],
                              ['I2_VH', 'CH104_VH'],
                              ['I1_VH', 'CH103_VH'],
                              ['I1_VH', 'CH106_VH']
                              ]

with open('../../results/Liao_ancestors_analysis/Liao_ancestors_S_NS_changes.csv', 'w') as output_file:
    output_file.write('pair,substitution_class,region,S5F_mutability_change\n')
    for pair in ancestral_descendant_pairs:
        seq_diffs = sequence_differences(parent_sequence = ancestors[pair[0]], descendant_sequence = ancestors[pair[1]],
                         partition_points = partition_points)

        for region in ['WS', 'FR','CDR']:

            if region == 'WS':
                region_label = ''
            else:
                region_label = '_' + region

            for substitution_class in ['syn','nonsyn','total']:

                if substitution_class == 'total':
                    mutability_change = seq_diffs['mutability_change_syn' + region_label] + \
                                        seq_diffs['mutability_change_nonsyn' + region_label]
                else:
                    mutability_change = seq_diffs['mutability_change_' + substitution_class + region_label]

                mutability_change = str(mutability_change)

                output_file.write('-'.join(pair) + ',' + substitution_class + ',' + region + ',' + mutability_change)
                output_file.write('\n')





