import sys

# Import mutability functions and partition points from analyses/mutability folder
sys.path.insert(0, '../mutability/')

from mutability_function import seq_mutability, aggregated_mutability, S5F
from partition_points import partition_points_dic
from gen_code_DNA import genetic_code
from dendropy import Tree

import re
from numpy import random
from numpy import log
from copy import deepcopy
import csv
from itertools import permutations

test_parent_sequence = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTCACTACTGGAGTTGGATCCGGCAGCCCCCAGGAAAGGGACTGGAGTGGATTGGCTATATCTATTACACTGGGAGCACCAACTACAATCCCTCTTTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAACTGACCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCCTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
test_descendant_sequence = 'TCAGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTGACTACTGGACTTGGATTCGGCAGGCCCCAGGAAAGGGACCGGAGTGGATTGGCTATATCTATTACACTGGGAGCACCAACTACAATCCCTCTTTCGCGAGTCGAATCACCATATCAGTAGACACGTCCAAGAAGCACTTCTCCCTGAAACTGACCTCTGTGACCGCTGCGGACACGGCCGTCTATTACTGTGCGAGTCTGCCCAGGGGGCAGTTAGTTAATGCCTACTTTGACTCCTGGGGCAAGGGAACCCTGGTCACCGTCGCCTCC'
test_partition_points = [1, 34, 58, 109, 130, 244, 288]

# Import conditional transition probabilities for each five-mer (Yaari et al.):
S5F_transitions = {}
with open('Yaari_transition_probabilities.csv', 'r') as csvfile:
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

nt_type = {'A':'purine','G':'purine','C':'pyrimidine', 'T':'pyrimidine'}

# Function for computing transition probabilities for a site based on the motif they're at the center of:
def transition_probs(motif,transition_model):
    if transition_model == 'S5F':
        tprobs = S5F_transitions[motif]
    elif transition_model == 'uniform':
        initial_nt = motif[2]
        tprobs = {}
        for base in ['A', 'T', 'C', 'G']:
            if base is initial_nt:
                tprobs[base] = 0
            else:
                tprobs[base] = 1. / 3

    return tprobs

# Function for computing relative mutability of sites based on its 5-nucleotide motif
def relative_mutabilities(seq, mutability_model):
    """
    :param seq:
    :param mutability_model: either a string ('uniform' or 'S5F') or a list with 3 floats, indicating the relative rates of codon positions 1,2 and 3
    :return:
    """
    if mutability_model == 'S5F':
        # Sites at the edges (for which S5F can't be computed) will have mutability 0
        mutabilities = [0] * len(seq)
        for site in range(2, (len(seq) - 2)):
            mutabilities[site] = S5F[seq[site - 2:site + 3]]

    elif mutability_model == 'uniform':
        # For consistency with the S5F based model, first 2 and last 2 sites will have mutability zero
        mutabilities = [0] * len(seq)
        for site in range(2, (len(seq) - 2)):
            mutabilities[site] = 1.0
    elif type(mutability_model) == list:
        mutabilities = [0] * len(seq)
        for site in range(2, (len(seq) - 2)):
            position = site %3
            mutabilities[site] = mutability_model[position]

    return mutabilities


def sequence_differences(parent_sequence, descendant_sequence, partition_points):
    '''Function for calculating differences in mean S5F and mean logS5F between a pair of sequences, partitioned into syn. and non-syn. contributions
        Also returns sequences with only syn. and only non-syn. differences, and a list of what 5-nucleotide motifs mutated
    '''

    assert (len(parent_sequence) == len(descendant_sequence))
    seq_length = len(parent_sequence)

    # Convert partition_points into Python indices (i.e. subtract 1):
    partition_points = [(partition_points[i] - 1) for i in range(len(partition_points))]

    # List of sites that belong to FRs (from partition_points)
    FR_sites = [range(partition_points[i], partition_points[i + 1]) for i in [0, 2, 4]]
    FR_sites = [j for i in FR_sites for j in i]

    # List of sites that belong to CDRs (from partition_points)
    CDR_sites = [range(partition_points[i], partition_points[i + 1]) for i in [1, 3, 5]]
    CDR_sites = [j for i in CDR_sites for j in i]
    CDR_sites += [partition_points[-1]]

    # Find number of nucleotide differences between parent and descendant
    diff_nts = [nt for nt in range(seq_length) if parent_sequence[nt] != descendant_sequence[nt]]
    diff_nts_FRs = [nt for nt in diff_nts if nt in FR_sites]
    diff_nts_CDRs = [nt for nt in diff_nts if nt in CDR_sites]

    n_nt_diffs = len(diff_nts)
    n_nt_diffs_FRs = len(diff_nts_FRs)
    n_nt_diffs_CDRs = len(diff_nts_CDRs)

    # Find number of AA differences in FRs and CDRs:
    n_aa_diffs_FRs = 0
    n_aa_diffs_CDRs = 0

    for nt in diff_nts_FRs:
        codon_site = nt / 3
        parent_codon = parent_sequence[codon_site * 3: codon_site * 3 + 3]
        descendant_codon = descendant_sequence[codon_site * 3: codon_site * 3 + 3]

        if set(parent_codon) <= {'A', 'G', 'C', 'T'} and set(descendant_codon) <= {'A', 'G', 'C', 'T'}:
            if genetic_code[parent_codon] != genetic_code[descendant_codon]:
                n_aa_diffs_FRs += 1

    for nt in diff_nts_CDRs:
        codon_site = nt / 3
        parent_codon = parent_sequence[codon_site * 3: codon_site * 3 + 3]
        descendant_codon = descendant_sequence[codon_site * 3: codon_site * 3 + 3]

        if set(parent_codon) <= {'A', 'G', 'C', 'T'} and set(descendant_codon) <= {'A', 'G', 'C', 'T'}:
            if genetic_code[parent_codon] != genetic_code[descendant_codon]:
                n_aa_diffs_CDRs += 1

                # % Divergence in FRs and CDRs:
    nt_divergence_FRs = float(n_nt_diffs_FRs) / len(FR_sites)
    nt_divergence_CDRs = float(n_nt_diffs_CDRs) / len(CDR_sites)

    aa_divergence_FRs = float(n_aa_diffs_FRs) / (len(FR_sites) / 3)
    aa_divergence_CDRs = float(n_aa_diffs_CDRs) / (len(CDR_sites) / 3)

    # Total number (in codon sites) of synonymous and non-synonymous differences
    n_syn_diffs = 0
    n_nonsyn_diffs = 0

    # Number of non-synonymous changes involving one, two and three nucleotide substitutions
    n_nonsyn_diffs_1nt = 0
    n_nonsyn_diffs_2nt = 0
    n_nonsyn_diffs_3nt = 0

    # Number of synonymous changes involving one, two and three nucleotide substitutions
    n_syn_diffs_1nt = 0
    n_syn_diffs_2nt = 0
    n_syn_diffs_3nt = 0

    syn_diff_sites = []
    nonsyn_diff_sites = []

    nts_in_syn_changes = []
    nts_in_nonsyn_changes = []

    # Change in mutability associated with synonymous changes
    mutability_change_syn = 0
    log_mutability_change_syn = 0

    # Change in mutability associated with nonsynonymous changes
    mutability_change_nonsyn = 0
    log_mutability_change_nonsyn = 0

    # Syn. and non-syn. changes partitioned by FRs and CDRs:
    mutability_change_syn_FR = 0
    mutability_change_nonsyn_FR = 0
    log_mutability_change_syn_FR = 0
    log_mutability_change_nonsyn_FR = 0

    mutability_change_syn_CDR = 0
    mutability_change_nonsyn_CDR = 0
    log_mutability_change_syn_CDR = 0
    log_mutability_change_nonsyn_CDR = 0

    # Number of motifs with more than one type of change (syn. and non-syn).
    motifs_with_syn_nonsyn = 0

    # Dictionary with number of transitions between each pair of amino acids
    aa_transitions = {}
    # aa_transitions['A']['B'] will give number of transitions from A to B

    # Dictionary of summed change in mutability associated with each amino acid transition
    # Summed across multiple occurrences of the same transition
    aa_trans_logmut_changes = {}

    # If at least one nucleotide changed, find all sites where syn. and non-syn. changes happened
    # Codon changes with >1 nt change are counted as only 1 syn. or non-syn. change.
    if n_nt_diffs > 0:
        diff_codon_sites = [i for i in range(seq_length / 3) if
                            descendant_sequence[i * 3: i * 3 + 3] != parent_sequence[i * 3: i * 3 + 3]]
        diff_codons = [[descendant_sequence[site * 3: site * 3 + 3], parent_sequence[site * 3: site * 3 + 3]] for site
                       in
                       diff_codon_sites]

        for i in range(len(diff_codon_sites)):
            site = diff_codon_sites[i]

            # Codons with non AGCT characters (i.e. gaps) will not be random. and don't contribute difs.
            if set(diff_codons[i][0]) <= {'A', 'G', 'C', 'T'} and set(diff_codons[i][1]) <= {'A', 'G', 'C', 'T'}:

                # Find positions in the codon (0,1,2) that have changed
                pos_changes = [k for k in range(3) if diff_codons[i][0][k] != diff_codons[i][1][k]]

                # get nt sites involved in change:
                nt_sites_changed = [range(site * 3, site * 3 + 3)[pos] for pos in pos_changes]

                if genetic_code[diff_codons[i][0]] == genetic_code[diff_codons[i][1]]:
                    n_syn_diffs += 1

                    nts_in_syn_changes += nt_sites_changed

                    pos_changes = len(pos_changes)

                    if pos_changes == 1:
                        n_syn_diffs_1nt += 1
                    if pos_changes == 2:
                        n_syn_diffs_2nt += 1
                    if pos_changes == 3:
                        n_syn_diffs_3nt += 1

                    syn_diff_sites.append(site)
                else:
                    pos_changes = len(pos_changes)

                    if pos_changes == 1:
                        n_nonsyn_diffs_1nt += 1
                    if pos_changes == 2:
                        n_nonsyn_diffs_2nt += 1
                    if pos_changes == 3:
                        n_nonsyn_diffs_3nt += 1

                    n_nonsyn_diffs += 1
                    nts_in_nonsyn_changes += nt_sites_changed
                    nonsyn_diff_sites.append(site)

                    descendant_aa = genetic_code[diff_codons[i][0]]
                    ancestor_aa = genetic_code[diff_codons[i][1]]

                    if ancestor_aa in aa_transitions.keys():
                        if descendant_aa in aa_transitions[ancestor_aa].keys():
                            aa_transitions[ancestor_aa][descendant_aa] += 1
                        else:
                            aa_transitions[ancestor_aa][descendant_aa] = 1
                    else:
                        aa_transitions[ancestor_aa] = {descendant_aa: 1}

    # For each five nucleotide motif position, compare parent and descendant
    for motif_pos in range(2, (seq_length - 2)):
        parent_motif = parent_sequence[motif_pos - 2:motif_pos + 3]
        descendant_motif = descendant_sequence[motif_pos - 2:motif_pos + 3]

        # Ignore gapped or ambiguous motifs
        if set(parent_motif) <= {'A', 'G', 'C', 'T'} and set(descendant_motif) <= {'A', 'G', 'C', 'T'}:

            # Do parent and descendant differ?
            if parent_motif != descendant_motif:
                delta_mutability = S5F[descendant_motif] - S5F[parent_motif]
                delta_log_mutability = log(S5F[descendant_motif]) - log(S5F[parent_motif])

                # Determine number of syn. and non-syn. substitutions associated with change in motif
                # List of nucleotide sites spanning motif
                nts_in_motif = range(motif_pos - 2, motif_pos + 3)

                # Number of nucleotides sites in motif that were involved in synonymous changes
                n_syn_nts = len([nt for nt in nts_in_motif if nt in nts_in_syn_changes])

                # Is the motif change associated with a synonymous codon change?
                syn_contribution = [1 if n_syn_nts > 0 else 0][0]

                # Number of nucleotides sites in motif that were involved in non-synonymous changes
                n_nonsyn_nts = len([nt for nt in nts_in_motif if nt in nts_in_nonsyn_changes])

                # Is the motif change associated with a non-synonymous codon change?
                nonsyn_contribution = [1 if n_nonsyn_nts > 0 else 0][0]

                if syn_contribution == 1 and nonsyn_contribution == 1:
                    motifs_with_syn_nonsyn += 1

                # In gapped sequences, some motif differences can't be determined to be syn. or non-syn.
                # Ignore those:
                if syn_contribution > 0 or nonsyn_contribution > 0:

                    # Change in S5F mutability attributable to syn changes:
                    delta_mutability_syn = delta_mutability * float(syn_contribution) / (
                    syn_contribution + nonsyn_contribution)
                    delta_log_mutability_syn = delta_log_mutability * float(syn_contribution) / (
                    syn_contribution + nonsyn_contribution)

                    # Change in S5F mutability attributable to nonsyn changes
                    delta_mutability_nonsyn = delta_mutability * float(nonsyn_contribution) / (
                    syn_contribution + nonsyn_contribution)
                    delta_log_mutability_nonsyn = delta_log_mutability * float(nonsyn_contribution) / (
                    syn_contribution + nonsyn_contribution)

                    # Add changes in mutability for this motif to sequence totals:
                    mutability_change_syn += delta_mutability_syn
                    log_mutability_change_syn += delta_log_mutability_syn

                    mutability_change_nonsyn += delta_mutability_nonsyn
                    log_mutability_change_nonsyn += delta_log_mutability_nonsyn

                    # Add changes in mutability to FR and CDR totals:

                    if motif_pos in FR_sites:
                        mutability_change_syn_FR += delta_mutability_syn
                        mutability_change_nonsyn_FR += delta_mutability_nonsyn

                        log_mutability_change_syn_FR += delta_log_mutability_syn
                        log_mutability_change_nonsyn_FR += delta_log_mutability_nonsyn

                    elif motif_pos in CDR_sites:
                        mutability_change_syn_CDR += delta_mutability_syn
                        mutability_change_nonsyn_CDR += delta_mutability_nonsyn

                        log_mutability_change_syn_CDR += delta_log_mutability_syn
                        log_mutability_change_nonsyn_CDR += delta_log_mutability_nonsyn

                    # If there is a non-synonymous contribution, record change into dictionary of changes per aa trans.
                    if nonsyn_contribution > 0:
                        # Determine which codon site spanned by the motif underwent NS changes
                        NS_diff_sites_in_motif = [nt / 3 for nt in nts_in_motif if nt in nts_in_nonsyn_changes]

                        # Number of AA changes that contributed to changes in this motif:
                        n_motif_aa_changes = len(NS_diff_sites_in_motif)

                        # For each AA change that contributed to change in motif:
                        for NS_diff_site in NS_diff_sites_in_motif:
                            # Find which aa transition this is
                            codon_pair = \
                            [diff_codons[i] for i in range(len(diff_codons)) if diff_codon_sites[i] == NS_diff_site][0]

                            # Divide nonsyn contribution equally among AA changes:
                            aa_logmut_change = float(delta_log_mutability_nonsyn) / n_motif_aa_changes

                            descendant_aa = genetic_code[codon_pair[0]]
                            ancestor_aa = genetic_code[codon_pair[1]]

                            # Add contribution of this AA trans to change in log mut. due to change in this motif
                            if ancestor_aa in aa_trans_logmut_changes.keys():
                                if descendant_aa in aa_trans_logmut_changes[ancestor_aa].keys():
                                    aa_trans_logmut_changes[ancestor_aa][descendant_aa] += aa_logmut_change
                                else:
                                    aa_trans_logmut_changes[ancestor_aa][descendant_aa] = aa_logmut_change
                            else:
                                aa_trans_logmut_changes[ancestor_aa] = {descendant_aa: aa_logmut_change}

    # Record motifs where substitutions happened, separately for syn. and non-syn. changes:
    mutated_motifs = {}

    for nt in diff_nts:
        if 1 < nt < (seq_length - 2):
            mut_motif = parent_sequence[nt - 2:nt + 3]

            # Ignore gapped or ambiguous motifs
            if set(mut_motif) <= {'A', 'G', 'C', 'T'}:

                # Was the motif in a codon site that underwent a synonymous or a non-synonymous substitution?
                if nt in nts_in_syn_changes:
                    sub_type = 'syn'
                elif nt in nts_in_nonsyn_changes:
                    sub_type = 'nonsyn'
                else:
                    # Some substitutions can'be classified into syn. or non-syn (e.g. trans. to gaps or ambiguous states)
                    sub_type = 'undefined'

                # Ignore those...
                if sub_type != 'undefined':
                    if mut_motif in mutated_motifs.keys():
                        if sub_type in mutated_motifs[mut_motif].keys():
                            mutated_motifs[mut_motif][sub_type] += 1
                        else:
                            mutated_motifs[mut_motif][sub_type] = 1
                    else:
                        mutated_motifs[mut_motif] = {sub_type: 1}

    # Averaging mutability differences (divide by seq_length - 4)
    mutability_change_syn = mutability_change_syn / (seq_length - 4)
    mutability_change_nonsyn = mutability_change_nonsyn / (seq_length - 4)
    log_mutability_change_syn = log_mutability_change_syn / (seq_length - 4)
    log_mutability_change_nonsyn = log_mutability_change_nonsyn / (seq_length - 4)

    mutability_change_syn_FR = mutability_change_syn_FR / (len(FR_sites) - 2)
    mutability_change_nonsyn_FR = mutability_change_nonsyn_FR / (len(FR_sites) - 2)
    log_mutability_change_syn_FR = log_mutability_change_syn_FR / (len(FR_sites) - 2)
    log_mutability_change_nonsyn_FR = log_mutability_change_nonsyn_FR / (len(FR_sites) - 2)

    mutability_change_syn_CDR = mutability_change_syn_CDR / (len(CDR_sites) - 2)
    mutability_change_nonsyn_CDR = mutability_change_nonsyn_CDR / (len(CDR_sites) - 2)
    log_mutability_change_syn_CDR = log_mutability_change_syn_CDR / (len(CDR_sites) - 2)
    log_mutability_change_nonsyn_CDR = log_mutability_change_nonsyn_CDR / (len(CDR_sites) - 2)

    # Converting contributions of AA transitions from changes in sum log mutability to changes in average log mutability
    for ancestor_aa in aa_trans_logmut_changes.keys():
        for descendant_aa in aa_trans_logmut_changes[ancestor_aa]:
            aa_trans_logmut_changes[ancestor_aa][descendant_aa] /= (seq_length - 4)

    output_dictionary = {}
    output_dictionary['n_nt_diffs'] = n_nt_diffs
    output_dictionary['n_syn_diffs'] = n_syn_diffs
    output_dictionary['n_syn_diffs_1nt'] = n_syn_diffs_1nt
    output_dictionary['n_syn_diffs_2nt'] = n_syn_diffs_2nt
    output_dictionary['n_syn_diffs_3nt'] = n_syn_diffs_3nt
    output_dictionary['n_nonsyn_diffs'] = n_nonsyn_diffs
    output_dictionary['n_nonsyn_diffs_1nt'] = n_nonsyn_diffs_1nt
    output_dictionary['n_nonsyn_diffs_2nt'] = n_nonsyn_diffs_2nt
    output_dictionary['n_nonsyn_diffs_3nt'] = n_nonsyn_diffs_3nt
    output_dictionary['nonsyn_diff_sites'] = nonsyn_diff_sites
    output_dictionary['syn_diff_sites'] = syn_diff_sites
    output_dictionary['mutated_motifs'] = mutated_motifs
    output_dictionary['mutability_change_syn'] = mutability_change_syn
    output_dictionary['mutability_change_nonsyn'] = mutability_change_nonsyn
    output_dictionary['mutability_change_syn_FR'] = mutability_change_syn_FR
    output_dictionary['mutability_change_nonsyn_FR'] = mutability_change_nonsyn_FR
    output_dictionary['mutability_change_syn_CDR'] = mutability_change_syn_CDR
    output_dictionary['mutability_change_nonsyn_CDR'] = mutability_change_nonsyn_CDR
    output_dictionary['log_mutability_change_syn'] = log_mutability_change_syn
    output_dictionary['log_mutability_change_nonsyn'] = log_mutability_change_nonsyn
    output_dictionary['log_mutability_change_syn_FR'] = log_mutability_change_syn_FR
    output_dictionary['log_mutability_change_nonsyn_FR'] = log_mutability_change_nonsyn_FR
    output_dictionary['log_mutability_change_syn_CDR'] = log_mutability_change_syn_CDR
    output_dictionary['log_mutability_change_nonsyn_CDR'] = log_mutability_change_nonsyn_CDR
    output_dictionary['motifs_with_syn_nonsyn'] = motifs_with_syn_nonsyn
    output_dictionary['nt_divergence_FRs'] = nt_divergence_FRs
    output_dictionary['nt_divergence_CDRs'] = nt_divergence_CDRs
    output_dictionary['aa_divergence_FRs'] = aa_divergence_FRs
    output_dictionary['aa_divergence_CDRs'] = aa_divergence_CDRs
    output_dictionary['aa_transitions'] = aa_transitions
    output_dictionary['aa_trans_logmut_changes'] = aa_trans_logmut_changes

    # Generate sequence with only syn. differences:
    syn_only_descendant = deepcopy(parent_sequence)
    syn_only_descendant = [syn_only_descendant[i:i + 3] for i in range(0, len(syn_only_descendant), 3)]

    for codon in output_dictionary['syn_diff_sites']:
        syn_only_descendant[codon] = [descendant_sequence[i:i + 3] for i in range(0, len(descendant_sequence), 3)][
            codon]
    syn_only_descendant = ''.join(syn_only_descendant)

    output_dictionary['syn_only_descendant'] = syn_only_descendant

    # Generate sequence with only non-syn. differences:
    nonsyn_only_descendant = deepcopy(parent_sequence)
    nonsyn_only_descendant = [nonsyn_only_descendant[i:i + 3] for i in range(0, len(nonsyn_only_descendant), 3)]

    for codon in output_dictionary['nonsyn_diff_sites']:
        nonsyn_only_descendant[codon] = [descendant_sequence[i:i + 3] for i in range(0, len(descendant_sequence), 3)][
            codon]
    nonsyn_only_descendant = ''.join(nonsyn_only_descendant)

    output_dictionary['nonsyn_only_descendant'] = nonsyn_only_descendant

    return output_dictionary

# Test that NS mean log S5F changes partitioned by amino acid transition type add to total nonsyn change
test_differences = sequence_differences(test_parent_sequence, test_descendant_sequence, test_partition_points)
change_nonsyn = test_differences['log_mutability_change_nonsyn']
aa_contribs = test_differences['aa_trans_logmut_changes']
contribs_list = []

for ancestor_aa in aa_contribs.keys():
    for descendant_aa in aa_contribs[ancestor_aa]:
        contribs_list.append(aa_contribs[ancestor_aa][descendant_aa])

assert (change_nonsyn - sum(contribs_list)) < 1e-10

# ================== SIMULATE DESCENDANT SEQUENCE KEEPING TRUE AA. DESCENDANT SEQ ====================

def randomize_sequence_constrained(parent_sequence, descendant_sequence, mutability_model, transition_model, partition_points):
    '''Function for generating random descendant sequence while constraining the amino acid sequence to be the same as the observed descendant
        Also returns sequences with only syn. and only non-syn. differences, and a list of what 5-nucleotide motifs mutated
    '''

    # Initialize randomized sequence:
    randomized_sequence = deepcopy(parent_sequence)

    seq_length = len(parent_sequence)
    
    seq_differences = sequence_differences(parent_sequence,descendant_sequence,partition_points)

    n_nt_diffs = seq_differences['n_nt_diffs']
    n_syn_diffs = seq_differences['n_syn_diffs']
    n_syn_diffs_1nt = seq_differences['n_syn_diffs_1nt']
    n_syn_diffs_2nt = seq_differences['n_syn_diffs_2nt'] 
    n_syn_diffs_3nt = seq_differences['n_syn_diffs_3nt']
    n_nonsyn_diffs = seq_differences['n_nonsyn_diffs']
    nonsyn_diff_sites = seq_differences['nonsyn_diff_sites']
    syn_diff_sites = seq_differences['syn_diff_sites']


    # 1. ------- Randomize non-synonymous changes exactly where they happened -----------

    # Count randomizable non-synonymous changes:
    n_nonsyn_randomizable = 0
    randomizable_nonsyn_sites = []

    # For each site where a non-synonymous codon change happened
    for ns_site in nonsyn_diff_sites:

        # Find ancestral and descendant codons
        ancestral_codon = parent_sequence[ns_site * 3: ns_site * 3 + 3]
        descendant_codon = descendant_sequence[ns_site * 3: ns_site * 3 + 3]

        # Ignoring gapped or ambiguous sites:
        if ancestral_codon in genetic_code.keys() and descendant_codon in genetic_code.keys():

            # Count nucleotide differences between ancestral and descendant codons
            codon_nt_diffs = len([j for j in range(3) if ancestral_codon[j] != descendant_codon[j]])

            # Find all codons that code for the same amino acid as the descendant and are at same distance from ancestor
            possible_codons = [codon for codon in genetic_code.keys() if
                               genetic_code[codon] == genetic_code[descendant_codon]]

            possible_codons = [codon for codon in possible_codons if
                               len([k for k in range(3) if codon[k] != ancestral_codon[k]]) == codon_nt_diffs]

            # A randomization is possible only if there's more than one possible codon:
            # If randomization not possible, change to DESCENDANT codon (skipping the loop to improve speed)
            if len(possible_codons) > 1:

                # Find the 2 closest nucleotides to the left and right of codon (from ANCESTRAL sequence)
                left_neighbors = parent_sequence[(ns_site * 3) - 2: ns_site * 3]
                right_neighbors = parent_sequence[(ns_site * 3) + 3: (ns_site * 3) + 5]

                # Because of S5F, randomization requires knowing the left and right neighbors
                # Otherwise skip site
                if len(left_neighbors) > 0 and len(right_neighbors) > 0 and \
                                set(left_neighbors) <= {'A', 'G', 'C', 'T'} and \
                                set(right_neighbors) <= {'A', 'G', 'C', 'T'}:

                    n_nonsyn_randomizable += 1
                    randomizable_nonsyn_sites.append(ns_site)

                    # List storing the relative probabilities of each candidate codon
                    candidate_probabilities = []

                    # For each candidate codon, find the relative probability of ancestral -> candidate
                    for candidate_codon in possible_codons:
                        diff_positions = [j for j in range(3) if ancestral_codon[j] != candidate_codon[j]]

                        # Multiple trajectories will be possible if codon change involves more than 1 nt change
                        # We consider all trajectories that don't have reversals or double hits
                        trajectories_list = []
                        trajectories_probabilities = []

                        for position_order in list(permutations(diff_positions)):
                            trajectory = [ancestral_codon]
                            current_codon = deepcopy(ancestral_codon)
                            for position in position_order:
                                next_codon = list(deepcopy(current_codon))
                                next_codon[position] = candidate_codon[position]
                                next_codon = ''.join(next_codon)

                                trajectory.append(next_codon)
                                current_codon = deepcopy(next_codon)

                            trajectories_list.append(trajectory)

                        for trajectory in trajectories_list:
                            trajectory_prob = 1
                            for step in range(len(trajectory) - 1):
                                initial_codon = trajectory[step]
                                next_codon = trajectory[step + 1]

                                mut_position = [k for k in range(3) if next_codon[k] != initial_codon[k]][0]

                                mutant_nt = next_codon[mut_position]

                                if mut_position == 0:
                                    motif = left_neighbors + initial_codon
                                if mut_position == 1:
                                    motif = left_neighbors[1] + initial_codon + right_neighbors[0]
                                if mut_position == 2:
                                    motif = initial_codon + right_neighbors

                                # Find relative mutation prob conditional on motif:
                                if mutability_model == 'S5F':
                                    motif_mutability = S5F[motif]
                                elif mutability_model == 'uniform':
                                    motif_mutability = 1
                                elif type(mutability_model) == list:
                                    motif_mutability = mutability_model[mut_position]

                                # Find transition prob conditional on motif:
                                trans_prob = transition_probs(motif, transition_model)
                                trans_prob = trans_prob[mutant_nt]

                                # Relative probability of trajectory:

                                trajectory_prob = trajectory_prob * trans_prob * motif_mutability

                            trajectories_probabilities.append(trajectory_prob)

                        # Sum prob. for all trajectories leading to candidate codon
                        candidate_probabilities.append(sum(trajectories_probabilities))

                    # Normalize codon probabilities
                    candidate_probabilities = [prob/sum(candidate_probabilities) for prob in candidate_probabilities]

                    # Sample codon based on its probability
                    new_codon = random.choice(possible_codons,1,p=candidate_probabilities)[0]

                    #print ancestral_codon + ' -> ' + new_codon + '; Site: ' + str(ns_site) + '(RANDOMIZATION)'
                else:
                    new_codon = descendant_codon
                    #print ancestral_codon + ' -> ' + new_codon + '; Site: ' + str(ns_site)
            else:
                new_codon = descendant_codon
                #print ancestral_codon + ' -> ' + new_codon + '; Site: ' + str(ns_site)

        randomized_sequence = list(randomized_sequence)
        randomized_sequence[ns_site * 3: ns_site * 3 + 3] = new_codon
        randomized_sequence = ''.join(randomized_sequence)

    # 2. ------- Randomize synonymous changes, avoiding double hits and codon sites where NS changes happened ----------

    # List of nucleotide sites to avoid, starting with sites in codons that underwent non-synonymous changes
    fixed_sites = []
    for ns_site in nonsyn_diff_sites:
        fixed_sites = fixed_sites + range(ns_site * 3, ns_site * 3 + 3)

    # We introduce the same number synonymous changes (in number of sites with a syn. difference)
    #total_syn_changes = n_syn_diffs_1nt + 2*n_syn_diffs_2nt + 3*n_syn_diffs_3nt
    total_syn_changes = n_syn_diffs

    # Syn mutations involving 2 or 3 nt changes are not allowed, but they are rarely seen in the observed trees


    # Counter for the number of syn changes with a single nt change
    syn_changes = 0
    while syn_changes < total_syn_changes:
        possible_sites = [site for site in range(seq_length) if site not in fixed_sites]

        site_mutabilities = relative_mutabilities(randomized_sequence, mutability_model)

        # Get only mutabilities for possible sites
        site_mutabilities = [site_mutabilities[i] for i in range(len(site_mutabilities)) if i in possible_sites]

        # Re-normalize mutability
        site_mutabilities = [mutability/sum(site_mutabilities) for mutability in site_mutabilities]

        # Choose a site to be mutated
        mutated_site = random.choice(possible_sites, 1, p=site_mutabilities)[0]

        # Codon position for checking mutation is synonymous
        codon_position = mutated_site / 3  # deliberately not a float division

        # Motif centered at mutated site
        motif = randomized_sequence[mutated_site - 2:mutated_site + 3]

        tprobs = transition_probs(motif, transition_model)
        nts = tprobs.keys()
        tprobs = [tprobs[key] for key in tprobs.keys()]

        #Choose a new nucleotide
        new_nt = random.choice(nts, 1, p=tprobs)[0]

        # Compare amino acids of simulated codon and ancestral sequence
        ancestral_codon = parent_sequence[codon_position * 3: codon_position * 3 + 3]

        sim_codon = list(randomized_sequence[codon_position * 3: codon_position * 3 + 3])
        sim_codon[mutated_site % 3] = new_nt
        sim_codon = ''.join(sim_codon)

        if genetic_code[sim_codon] == genetic_code[ancestral_codon]:
            syn_changes += 1
            randomized_sequence = list(randomized_sequence)
            randomized_sequence[mutated_site] = new_nt
            randomized_sequence = ''.join(randomized_sequence)

            fixed_sites = fixed_sites + range(codon_position * 3, codon_position * 3 + 3)

            #print ancestral_codon + ' -> ' + sim_codon + '; Site: ' + str(codon_position) + '(RANDOMIZATION)'

    return [randomized_sequence, n_nonsyn_randomizable, randomizable_nonsyn_sites]

def randomize_sequence_unconstrained(parent_sequence, descendant_sequence, mutability_model, transition_model, partition_points):
    '''Similar to randomize_sequence_constrained, but without constraining the amino acid sequence of the randomized descendant sequence to be the same as the observed descendant, and without constraining codon changes to involve a single nucleotide change.
    The number of non-syn. and syn. substitutions is still constrained to be the same as the observed numbers
    '''

    # Initialize randomized sequence:
    randomized_sequence = deepcopy(parent_sequence)

    seq_length = len(parent_sequence)

    seq_differences = sequence_differences(parent_sequence, descendant_sequence, partition_points)

    # Observed numbers of syn. and non-syn. differences
    obs_syn_diffs = seq_differences['n_syn_diffs']
    obs_nonsyn_diffs = seq_differences['n_nonsyn_diffs']

    # List of nucleotide sites to avoid. Sites in previously hit codons cannot be hit again.
    fixed_sites = []

    # Initialize counters for the number of syn and non-syn changes
    syn_changes = 0
    nonsyn_changes = 0

    # While numbers of syn. and non-syn changes are less than observed numbers
    while (syn_changes < obs_syn_diffs) or (nonsyn_changes < obs_nonsyn_diffs):

        # List of sites that can be mutated
        possible_sites = [site for site in range(seq_length) if site not in fixed_sites]

        # Check if ran out of sites to randomize
        assert len(possible_sites) > 0, 'Ran out of sites to randomize due to multiple hits not being allowed.'

        # Get only mutabilities for possible sites
        site_mutabilities = relative_mutabilities(randomized_sequence, mutability_model)
        site_mutabilities = [site_mutabilities[i] for i in range(len(site_mutabilities)) if i in possible_sites]

        # Re-normalize mutability
        site_mutabilities = [mutability/sum(site_mutabilities) for mutability in site_mutabilities]

        # Choose a site to be mutated
        mutated_site = random.choice(possible_sites, 1, p=site_mutabilities)[0]

        # Codon position of mutated site
        codon_position = mutated_site / 3  # deliberately not a float division

        # Motif centered at mutated site
        motif = randomized_sequence[mutated_site - 2:mutated_site + 3]

        tprobs = transition_probs(motif, transition_model)
        nts = tprobs.keys()
        tprobs = [tprobs[key] for key in tprobs.keys()]

        #Choose a new nucleotide
        new_nt = random.choice(nts, 1, p=tprobs)[0]

        # Compare amino acids of simulated codon and ancestral sequence
        ancestral_codon = parent_sequence[codon_position * 3: codon_position * 3 + 3]

        sim_codon = list(randomized_sequence[codon_position * 3: codon_position * 3 + 3])
        sim_codon[mutated_site % 3] = new_nt
        sim_codon = ''.join(sim_codon)

        # If simulated codon and ancestral code for the same amino acid...
        if genetic_code[sim_codon] == genetic_code[ancestral_codon] and syn_changes < obs_syn_diffs:
            #syn_changes += 1

            randomized_sequence = list(randomized_sequence)
            randomized_sequence[mutated_site] = new_nt
            randomized_sequence = ''.join(randomized_sequence)

            # Prevent all sites in the mutated codon to be mutated in the future
            #fixed_sites = fixed_sites + range(codon_position * 3, codon_position * 3 + 3)

        # else if simulated codon codes for an amino acid change
        elif genetic_code[sim_codon] != genetic_code[ancestral_codon] and nonsyn_changes < obs_nonsyn_diffs:
            #if sim. codon is not a stop codon (if it is, skip the attempt and randomize again)
            if genetic_code[sim_codon] != '*':
               # nonsyn_changes += 1

                randomized_sequence = list(randomized_sequence)
                randomized_sequence[mutated_site] = new_nt
                randomized_sequence = ''.join(randomized_sequence)

                # Prevent all sites in the mutated codon to be mutated in the future
                #fixed_sites = fixed_sites + range(codon_position * 3, codon_position * 3 + 3)

        #Update realized number of syn and nonsyn changes:
        randomized_vs_observed = sequence_differences(parent_sequence, randomized_sequence, partition_points)
        syn_changes = randomized_vs_observed['n_syn_diffs']
        nonsyn_changes = randomized_vs_observed['n_nonsyn_diffs']

        #print str(syn_changes) + '/' + str(nonsyn_changes)

    # Check that randomized sequence has same # of syn. and non-syn. diffs from ancestral seq. as observed descendant
    randomized_seq_differences = sequence_differences(parent_sequence, randomized_sequence, partition_points)
    assert randomized_seq_differences['n_syn_diffs'] == obs_syn_diffs, 'Randomized descendant has a different number of syn. differences from ancestor relative to observed descendant'
    assert randomized_seq_differences['n_nonsyn_diffs'] == obs_nonsyn_diffs, 'Randomized descendant has a different number of nonsyn. differences from ancestor relative to observed descendant'

    return randomized_sequence
