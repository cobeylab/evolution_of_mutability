#!/usr/bin/python

"""
Looks at ancestral CH103 sequences inferred by Liao et al. (2013, Nature), and partitions changes in mutability for ancestral-descendant pairs.
"""
import sys
from Bio import SeqIO

# Partition points obtained from IMGT (starting from 1. Mutability functions will adjust to python indexing)
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

# Function for identifying region where codon change occurred (FR or CDR)
def region_of_difference(changed_codon_site, partition_points):

    # Partition points in codon positions (starting from 0)
    codon_partition_points = [(point - 1) / 3 for point in partition_points]

    if changed_codon_site > codon_partition_points[-1]:
        # Change is outside the FR1 to CDR3 region
        region_type = 'None'

    else:
        region_start = max([point for point in codon_partition_points[:-1] if changed_codon_site >= point])

        region_start_index = [i for i in range(len(codon_partition_points)) if codon_partition_points[i] == region_start][0]

        if region_start_index in [0, 2, 4]:
            region_type = 'FR'
        elif region_start_index in [1,3,5]:
            region_type = 'CDR'

    return region_type

assert(region_of_difference(25, partition_points) == 'CDR')
assert(region_of_difference(24, partition_points) == 'FR')
assert(region_of_difference(34, partition_points) == 'FR')


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
    output_file.write('pair,substitution_class,region,n_codon_changes,S5F_mutability_change,logS5F_mutability_change\n')

    for pair in ancestral_descendant_pairs:
        seq_diffs = sequence_differences(parent_sequence = ancestors[pair[0]], descendant_sequence = ancestors[pair[1]],
                         partition_points = partition_points)

        # Find codon sites that underwent syn. and non-syn. substitutions
        nonsyn_diff_codon_sites = seq_diffs['nonsyn_diff_sites']
        syn_diff_codon_sites = seq_diffs['syn_diff_sites']

        # Regions where syn. and non-syn. substitutions occurred
        nonsyn_diff_regions = [region_of_difference(site, partition_points) for site in nonsyn_diff_codon_sites]
        syn_diff_regions = [region_of_difference(site, partition_points) for site in syn_diff_codon_sites]

        # Count # of syn. and non-syn. substitutions in the pair, separately for FRs and CDRs
        n_nonsyn_diffs = {'WS':len(nonsyn_diff_codon_sites),
                          'CDR': nonsyn_diff_regions.count('CDR'),
                          'FR': nonsyn_diff_regions.count('FR')}

        n_syn_diffs = {'WS': len(syn_diff_codon_sites),
                       'CDR': syn_diff_regions.count('CDR'),
                       'FR': syn_diff_regions.count('FR')}

        n_total_diffs = {'WS': len(syn_diff_codon_sites) + len(nonsyn_diff_codon_sites),
                         'CDR': n_nonsyn_diffs['CDR'] + n_syn_diffs['CDR'],
                         'FR': n_nonsyn_diffs['FR'] + n_syn_diffs['FR']
                         }

        n_codon_diffs = {'total': n_total_diffs, 'syn': n_syn_diffs, 'nonsyn': n_nonsyn_diffs}


        for region in ['WS', 'FR','CDR']:

            if region == 'WS':
                region_label = ''
            else:
                region_label = '_' + region

            for substitution_class in ['syn','nonsyn','total']:

                if substitution_class == 'total':

                    # Get change in mean S5F
                    mutability_change = seq_diffs['mutability_change_syn' + region_label] + \
                                        seq_diffs['mutability_change_nonsyn' + region_label]

                    # Get change in mean log S5F
                    log_mutability_change = seq_diffs['log_mutability_change_syn' + region_label] + \
                                        seq_diffs['log_mutability_change_nonsyn' + region_label]

                else:
                    # Get change in mean S5F
                    mutability_change = seq_diffs['mutability_change_' + substitution_class + region_label]

                    # Get change in mean log S5F
                    log_mutability_change = seq_diffs['log_mutability_change_' + substitution_class + region_label]

                mutability_change = str(mutability_change)
                log_mutability_change = str(log_mutability_change)

                n_codon_changes = str(n_codon_diffs[substitution_class][region])

                output_file.write('-'.join(pair) + ',' + substitution_class + ',' + region + ',' + n_codon_changes)
                output_file.write(',' + mutability_change + ',' + log_mutability_change)
                output_file.write('\n')





