"""Takes each branch in an MCC tree, computes its change in mutability and compares it to changes in mutability given randomized descendant sequences under different models (with fixed amino acid sequence.
"""
import sys
# Import mutability functions from analyses/mutability folder
sys.path.insert(0, '../mutability/')
# Import partition_points dictionary from rate_correlations folder
sys.path.insert(0, '../rate_correlations/')

# Import genetic code from simulations folder:
sys.path.insert(0, '../simulations/modules/')

from mutation_functions import sequence_differences, randomize_sequence_constrained

from mutability_function import seq_mutability, aggregated_mutability, S5F
from partition_points import partition_points_dic
from gen_code_DNA import genetic_code
from dendropy import Tree

import re
from numpy import random, percentile, mean
from copy import deepcopy
import csv
from itertools import permutations

def main(argv):
    # Chain id, e.g. CH103_con_run1a, VRC01_L08_log_run1b, scenario2a_rep1
    chain_id = str(argv[1])

    # ================== RETRIEVE DIRECTORIES, OBS. SEQUENCES AND FR/CDR PARTITION POINTS FROM CHAIN ID ================
    # If chain is from an observed clone:
    if chain_id.find('scenario') is -1:
        if chain_id.find('CH103') > -1:
            if chain_id.find('CH103L') > -1:
                clone = 'CH103L'
            else:
                clone = 'CH103'
        elif chain_id.find('VRC26') > -1:
            if chain_id.find('VRC26L') > -1:
                clone = 'VRC26L'
            else:
                clone = 'VRC26'
        elif chain_id.find('VRC01') > -1:
            if chain_id.find('VRC01_01') > -1:
                clone = 'VRC01_01'
            elif chain_id.find('VRC01_13') > -1:
                clone = 'VRC01_13'
            elif chain_id.find('VRC01_19') > -1:
                clone = 'VRC01_19'
            elif chain_id.find('VRC01_H0306') > -1:
                clone = 'VRC01_H0306'
            elif chain_id.find('VRC01_L0306') > -1:
                clone = 'VRC01_L0306'
            elif chain_id.find('VRC01_H08') > -1:
                clone = 'VRC01_H08'
            elif chain_id.find('VRC01_L08') > -1:
                clone = 'VRC01_L08'

        if chain_id.find('_con_') > -1:
            prior = 'constant'
        elif chain_id.find('_exp_') > -1:
            prior = 'exponential'
        elif chain_id.find('_log_') > -1:
            prior = 'logistic'

        # (TEMPORARY) if chain is an interrupted chain, add 'int' to clone:
        if chain_id.find('int') > -1:
            clone = clone + 'int'

        partition_points = partition_points_dic[clone]
        chain_directory = '../../results/BEAST/observed_lineages/' + clone + '_' + prior + '/' + chain_id + '/'
        MCC_tree_file_path = chain_directory + chain_id + '_MCC_tree.tree'

        xml_file_path = '../../analyses/BEAST/observed_lineages/' + clone + '_' + prior + '/'
        xml_file_path += chain_id[0:len(chain_id) - 1] + '.xml'

        output_directory = '../../results/S_NS_mutability_changes/observed_lineages/' + clone + '_' + prior + '/'

        output_file_path_observed = output_directory + chain_id + '_observed_MCC.csv'
        output_file_path_simulated = output_directory + chain_id + '_simulated_MCC.csv'
        output_file_path_motifs = output_directory + chain_id + '_motif_changes_MCC.csv'



    # If chain is from a simulated scenario
    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        MCC_tree_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        MCC_tree_file_path += chain_id + '_MCC_tree.tree'

        xml_file_path = '../../analyses/BEAST/simulated_alignments/' + chain_id + '/' + chain_id + '.xml'

        output_directory = '../../results/S_NS_mutability_changes/simulated_alignments/' + chain_id + '/'
        output_file_path_observed = output_directory + chain_id + '_observed_MCC.csv'
        output_file_path_simulated = output_directory + chain_id + '_simulated_MCC.csv'
        output_file_path_motifs = output_directory + chain_id + '_motif_changes_MCC.csv'

        partition_points = partition_points_dic[chain_id]


    # ========================== READ TIP SEQUENCES FROM XML FILE AND COMPUTE THEIR MUTABILITY==========================
    # ------------------------------- (Tip sequences are not annotated on the BEAST trees) -----------------------------

    # Dictionaries with mutability for each tip sequence:
    mutability_WS_tips = {}
    mutability_aggregated_tips = {}

    # Dictionary with observed sequences
    obs_sequence = {}

    # Read XML as string
    with open(xml_file_path, 'r') as xml_file:
        xml = xml_file.readlines()
        xml = ''.join(xml)

    taxon_lines = re.findall(r'<sequence>.*</sequence>', xml, re.DOTALL)[0].split('</sequence>')
    for line in taxon_lines[0:len(taxon_lines) - 1]:
        taxon_id = re.search(r'taxon idref=".*"', line).group().replace("taxon idref=", '').replace('"', '')
        # Ignoring a weird VRC26 sequence whose CDR2 is entirely gaps:
        if taxon_id != 'KJ134124_119':
            sequence = re.search(r'/>\n\t\t\t.*\n\t\t', line).group().replace('/>\n\t\t\t', '').replace('\n\t\t', '')

            mutability_WS_tips[taxon_id] = seq_mutability(sequence)
            obs_sequence[taxon_id] = sequence

            if partition_points is not None:
                mutability_aggregated_tips[taxon_id] = aggregated_mutability(sequence, partition_points)
        else:
            print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'


    # Get length of sequences
    seq_length = len(obs_sequence[obs_sequence.keys()[0]])

    # ================ GET RELATIVE MUTATION RATES FOR EACH CODON POSITION FROM BEAST LOG SUMMARY FILE =================
    codon_position_rates = [0,0,0]

    log_summary_file_path = '/'.join(MCC_tree_file_path.split('/')[:-1]) + '/log_summary.csv'
    with open(log_summary_file_path,'r') as log_summary_file:
        for line in log_summary_file:
            if line.find('CP1.mu') > -1:
                codon_position_rates[0] = float(line.split('\t')[1])
            if line.find('CP2.mu') > -1:
                codon_position_rates[1] = float(line.split('\t')[1])
            if line.find('CP3.mu') > -1:
                codon_position_rates[2] = float(line.split('\t')[1])

    #global codon_relative_rates

    #randomize_sequence('AAACCCGGGTTT','AACCCGGGTTTA', codon_position_rates,'S5F')

    # ================================ DICTIONARY WITH N CHANGES FOR EACH MOTIF ========================================
    motif_changes = {}


    # ==================================== OPEN TREE FILE AND OUTPUT FILES =============================================

    with open(MCC_tree_file_path, 'r') as tree_file, open(output_file_path_observed, 'w') as output_file_obs, open(
            output_file_path_simulated, 'w') as output_file_sim:

        # Header shared by observed and simulated output files:
        output_header_base = 'parent,child,parent_distance_to_root,parent_time_to_root,branch_time,branch_exp_subs,'
        output_header_base += 'n_NT_diffs,n_NonSyn_diffs,n_NonSyn_diffs_1nt,n_NonSyn_diffs_2nt,n_NonSyn_diffs_3nt,'
        output_header_base += 'n_Syn_diffs,n_Syn_diffs_1nt,n_Syn_diffs_2nt,n_Syn_diffs_3nt,'
        output_header_base += 'n_nonsyn_changes_randomizable,n_motifs_with_syn_nonsyn_changes,branch_is_terminal,'
        output_header_base += 'branch_in_trunk,n_descendants_last_time'

        output_header_obs = deepcopy(output_header_base)
        # Header for output file with observed mutability
        # Ignoring other metrics and focusing on S5F:
        for metric in ['S5F']:
            output_header_obs += ',' + metric + '_parent'
            output_header_obs += ',' + metric + '_FR_parent'
            output_header_obs += ',' + metric + '_CDR_parent'

            output_header_obs += ',' + metric + '_child'
            output_header_obs += ',' + metric + '_FR_child'
            output_header_obs += ',' + metric + '_CDR_child'


            output_header_obs += ',' + metric + '_change_total'
            output_header_obs += ',' + metric + '_change_syn'
            output_header_obs += ',' + metric + '_change_nonsyn'

            output_header_obs += ',' + metric + '_change_syn_FR'
            output_header_obs += ',' + metric + '_change_nonsyn_FR'

            output_header_obs += ',' + metric + '_change_syn_CDR'
            output_header_obs += ',' + metric + '_change_nonsyn_CDR'


        output_header_obs += '\n'

        output_file_obs.write(output_header_obs)
        
        # Header for output file with simulation results:
        output_header_sim = 'replicate,' + output_header_base
        for region in ['WS']:
        # Looking at whole-sequence only for now
            #for metric in ['S5F','7M','HS','CS','OHS']:
            for metric in ['S5F']:

                # Mutability changes under S5F mutations/S5F transitions model
                output_header_sim += ',' + metric + '_' + region + '_change_S5FMut_S5FTrans_total'
                output_header_sim += ',' + metric + '_' + region + '_change_S5FMut_S5FTrans_syn'
                output_header_sim += ',' + metric + '_' + region + '_change_S5FMut_S5FTrans_nonsyn'

                # Mutability changes under uniform mutations/S5F transitions model
                output_header_sim += ',' + metric + '_' + region + '_change_uniformMut_S5FTrans_total'
                output_header_sim += ',' + metric + '_' + region + '_change_uniformMut_S5FTrans_syn'
                output_header_sim += ',' + metric + '_' + region + '_change_uniformMut_S5FTrans_nonsyn'

                # Mutability changes under codon-position-based mutations / S5F based transitions model:
                output_header_sim += ',' + metric + '_' + region + '_change_CPMut_S5FTrans_total'
                output_header_sim += ',' + metric + '_' + region + '_change_CPMut_S5FTrans_syn'
                output_header_sim += ',' + metric + '_' + region + '_change_CPMut_S5FTrans_nonsyn'

        output_header_sim += '\n'

        output_file_sim.write(output_header_sim)

        # ========================== DICTIONARY LINKING BEAST NUMBERING (FROM NEXUS) TO TAXON LABELS ===================
        number_to_id = {}

        # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
        n_semicolons = 0
        while n_semicolons < 6:
            line = tree_file.next()

            translation_search = re.search(r'[0-9]+ [A-Z|0-9|_]+', line)
            # If lines contain a taxon-id-to-number translation...
            if translation_search is not None:
                translation = translation_search.group()
                number = re.search(r'[0-9]*', translation).group()
                taxon_id = translation.replace(number + ' ', '')
                number_to_id[number] = taxon_id

            # Make sure no trees are being skipped
            assert line.find('STATE') is -1, 'Header-skipping loop skipped a tree.'
            if line.find(';') > -1:
                n_semicolons += 1

        tree_line = tree_file.next()
        tree_string = deepcopy(tree_line)

        # From sequence ids, find maximum sampling time
        sampling_times = [float(re.search('_[0-9]+',seq_id).group().replace('_','')) for seq_id in number_to_id.values()]
        max_sampling_time = max(sampling_times)

        # =================================== PROCESSING TREE FOR DENDROPY =============================================
        # Dendropy doesn't seem to read annotations of BEAST MCC trees correctly
        # Get annotation blocks for all nodes, simplify them and input edited tree string to dendropy

        # Basically find all instances of "[...]" in the newick string:
        annotation_blocks = re.findall(r'\[[^\]]*]',tree_string)
        annotation_blocks = annotation_blocks[1:len(annotation_blocks)]

        # Find alignment identifier (e.g. CH103_final_alignment in CH103_final_alignment.CP1)
        alignment_id = re.search(r'[^,^\[]*.CP1', annotation_blocks[0]).group()
        alignment_id = alignment_id.replace('&', '').replace('.CP1', '')

        # Annotations that will be kept in the edited tree string
        variables = [alignment_id + '.CP1', alignment_id + '.CP2', alignment_id + '.CP3','rate','length']

        for block in annotation_blocks:

            new_block = '[&'

            for variable in variables:
                value = re.search(variable + '=[^,^\]]*', block)
                if value is not None:
                    #assert value.group() is not 'height=64.72259232302139'
                    new_block += value.group() + ','

            #Remove extra comma at the end:
            new_block = new_block[0:(len(new_block) - 1)]

            new_block += ']'

            #assert(new_block.find('height=64.72259232302139') is -1)

            tree_string = tree_string.replace(block,new_block)

        # Search for a tree header in the line
        tree_header = re.search(r'.*\[\&R\]', tree_string)

        tree_header = tree_header.group()

        # Remove tree header from tree_string before reading it with dendropy
        tree_string = tree_string.replace(tree_header, '')
        tree_string = tree_string.replace('\n', '')

        # ====================================== READING TREE WITH DENDROPY ============================================

        # Read tree string as a tree object using dendropy:
        tree = Tree.get_from_string(tree_string, schema='newick')

        # Label all nodes according to their index in tree.nodes():
        for i in range(len(tree.nodes())):
            tree.nodes()[i].label = str(i)
        # Notice that this is different from the numbering BEAST does on the tips

        # Get number of sites from codon position 1 annotation for node 0
        node0_CP1 = tree.nodes()[0].annotations.get_value(alignment_id + '.CP1')

        # Check if node0_CP1 is ambiguous (i.e. has more than one sequence connected by +)
        if node0_CP1.find('+') > -1:
            # If it is, arbitrarily choose the first sequence
            node0_CP1 = node0_CP1.split('+')[0]

        n_sites = 3 * len(node0_CP1)

        # ================================== DO THE ANALYSIS FOR EACH NODE =============================================
        for node in tree.nodes():
            # Find node number (index in tree.nodes())
            node_number = node.label

            # If node is not the root, find changes in mutability (contrasts) from parent node
            if (node.parent_node is not None):
                print node_number

                # Ignoring a weird VRC26 sequence whose CDR2 is entirely gaps:
                if node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'KJ134124_119':
                    print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'

                # Ignoring germline sequence:
                elif node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'GERMLINE_00':
                    print 'Skipping GERMLINE sequence'
                else:

                    # Find node's parent node
                    parent_node = node.parent_node

                    # Find parent node number (index in tree.nodes())
                    parent_node_number = parent_node.label

                    assert node.parent_node == tree.nodes()[int(parent_node_number)]

                    # Find time, mol. clock distance and RC based distance between parent node and root:

                    parent_time = parent_node.distance_from_root()

                    # To find parent node's distance from root, first time-lengths and rate for all branches connecting it to the root
                    # (robust-counting per-site distances can be found directly from 'S' and 'N' without reference to time lengths)

                    parent_to_root_time_lengths = [0.0]
                    parent_to_root_rates = [0.0]
                    #parent_to_root_lengths_RC = [0.0]
                    focal_node = parent_node

                    # Move down the tree to the root node
                    while focal_node is not tree.nodes()[0]:
                        parent_to_root_time_lengths.append(focal_node.edge_length)
                        parent_to_root_rates.append(float(focal_node.annotations.get_value('rate')))
                        #RC_length = float(focal_node.annotations.get_value('S')) + float(
                            #focal_node.annotations.get_value('N'))
                        #RC_length = float(RC_length) / n_sites
                        #parent_to_root_lengths_RC.append(RC_length)

                        focal_node = focal_node.parent_node

                    assert sum(parent_to_root_time_lengths) == parent_time

                    # Find parent's distance to the root by multiplying molecular clock rates and time lengths, then summing
                    parent_to_root_distance = sum([parent_to_root_time_lengths[i] * parent_to_root_rates[i] for i in
                                                   range(len(parent_to_root_rates))])

                    # Find parent's robust counting distance by summing robust counting branch lengths
                    #parent_to_root_distance_RC = sum(parent_to_root_lengths_RC)

                    # Find node mutability, if it is a tip:
                    if node.is_leaf():
                        # Find taxon number from BEAST tree
                        taxon_number = str(node.taxon).replace("'", '')
                        # Find node taxon id (the original name for the tip sequence):
                        node_id = number_to_id[taxon_number]

                        # Get mutability from tip mutability dictionaries
                        node_mutability_WS = mutability_WS_tips[node_id]
                        node_mutability_aggregated = mutability_aggregated_tips[node_id]

                        node_sequence = obs_sequence[node_id]
                    else:
                        # If node is not a tip, find its sequence from the tree annotation:
                        node_CP1 = node.annotations.get_value(alignment_id + '.CP1')
                        node_CP2 = node.annotations.get_value(alignment_id + '.CP2')
                        node_CP3 = node.annotations.get_value(alignment_id + '.CP3')

                        # Arbitrarily choosing first listed reconstruction in case of ambiguity
                        if node_CP1.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP1 = node_CP1.split('+')[0]
                            print 'Ambiguous reconstruction in node ' + str(node_number) + '; arbitrarily choosing first listed reconstruction'
                        if node_CP2.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP2 = node_CP2.split('+')[0]
                            print 'Ambiguous reconstruction in node ' + str(
                                node_number) + '; arbitrarily choosing first listed reconstruction'
                        if node_CP3.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP3 = node_CP3.split('+')[0]
                            print 'Ambiguous reconstruction in node ' + str(
                                node_number) + '; arbitrarily choosing first listed reconstruction'

                        node_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
                        node_sequence = ''.join(node_sequence)

                        # Compute mutability
                        node_mutability_WS = seq_mutability(node_sequence)
                        node_mutability_aggregated = aggregated_mutability(node_sequence, partition_points)


                    # Find mutability of parent node:

                    parent_node_CP1 = parent_node.annotations.get_value(alignment_id + '.CP1')
                    parent_node_CP2 = parent_node.annotations.get_value(alignment_id + '.CP2')
                    parent_node_CP3 = parent_node.annotations.get_value(alignment_id + '.CP3')

                    # Arbitrarily choosing first listed reconstruction in case of ambiguity
                    if parent_node_CP1.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP1 = parent_node_CP1.split('+')[0]
                        print 'Ambiguous reconstruction in node ' + str(
                            parent_node_number) + '; arbitrarily choosing first listed reconstruction'
                    if parent_node_CP2.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP2 = parent_node_CP2.split('+')[0]
                        print 'Ambiguous reconstruction in node ' + str(
                            parent_node_number) + '; arbitrarily choosing first listed reconstruction'
                    if parent_node_CP3.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP3 = parent_node_CP3.split('+')[0]
                        print 'Ambiguous reconstruction in node ' + str(
                            parent_node_number) + '; arbitrarily choosing first listed reconstruction'

                    parent_node_sequence = [parent_node_CP1[i] + parent_node_CP2[i] + parent_node_CP3[i] for i in
                                            range(len(parent_node_CP1))]
                    parent_node_sequence = ''.join(parent_node_sequence)

                    # Compute mutability of parent node
                    parent_node_mutability_WS = seq_mutability(parent_node_sequence)
                    parent_node_mutability_aggregated = aggregated_mutability(parent_node_sequence, partition_points)

                    # Compare parent and descendant sequences:
                    differences = sequence_differences(parent_node_sequence, node_sequence, partition_points)
                    n_nt_diffs = differences['n_nt_diffs']

                    n_NonSyn_diffs = differences['n_nonsyn_diffs']
                    n_NonSyn_diffs_1nt = differences['n_nonsyn_diffs_1nt']
                    n_NonSyn_diffs_2nt = differences['n_nonsyn_diffs_2nt']
                    n_NonSyn_diffs_3nt = differences['n_nonsyn_diffs_3nt']

                    n_Syn_diffs = differences['n_syn_diffs']
                    n_Syn_diffs_1nt = differences['n_syn_diffs_1nt']
                    n_Syn_diffs_2nt = differences['n_syn_diffs_2nt']
                    n_Syn_diffs_3nt = differences['n_syn_diffs_3nt']

                    # Get S5F mutability changes associated with synonymous and non-synonymous changes
                    S5F_WS_change_syn = differences['mutability_change_syn']
                    S5F_WS_change_nonsyn = differences['mutability_change_nonsyn']

                    S5F_change_syn_FR = differences['mutability_change_syn_FR']
                    S5F_change_nonsyn_FR = differences['mutability_change_nonsyn_FR']

                    S5F_change_syn_CDR = differences['mutability_change_syn_CDR']
                    S5F_change_nonsyn_CDR = differences['mutability_change_nonsyn_CDR']

                    # Get list of codon sites that can be randomized from randomization function
                    rdm_nonsyn_sites = randomize_sequence_constrained(parent_node_sequence, node_sequence, 'S5F', 'S5F', partition_points)[2]

                    # Count number of motifs where both synonymous and nonsynonymous changes happened
                    n_motifs_with_syn_nonsyn_changes = differences['motifs_with_syn_nonsyn']

                    # Add motifs in the sequences to list of motifs present in the lineage
                    for i in range(2, (seq_length - 2)):
                        obs_motif = parent_node_sequence[i - 2:i + 3]
                        if obs_motif not in motif_changes.keys() and set(obs_motif) <= {'A', 'G', 'C', 'T'}:
                            motif_changes[obs_motif] = {'syn':0, 'nonsyn':0}

                    for i in range(2, (seq_length - 2)):
                        obs_motif = node_sequence[i - 2:i + 3]
                        if obs_motif not in motif_changes.keys() and set(obs_motif) <= {'A', 'G', 'C', 'T'}:
                            motif_changes[obs_motif] = {'syn': 0, 'nonsyn': 0}

                    # Record motif changes observed on branch:
                    branch_motif_changes = differences['mutated_motifs']
                    for key in branch_motif_changes:
                        for sub_type in branch_motif_changes[key].keys():
                            motif_changes[key][sub_type] += branch_motif_changes[key][sub_type]

                    # Determine if parent->node branch is part of the trunk, and count branch's n. descendants
                    # Branch is in trunk if it has desc. at last time but not if it only has desc. from last time
                    # Lemey et al. 2007 PlosCompBio

                    # (n. descendants at the last time point)

                    # branch_in_trunk = 0
                    if node.is_leaf():
                        n_descendants_last_time = 'NA'
                        branch_in_trunk = 0

                    else:
                        # Find all sequences that descend from the node:
                        node_descendants = node.inorder_iter()
                        n_descendants_last_time = 0
                        n_descendants_other_times = 0

                        # Count the number of its descendants at the last time point and other time points
                        for descendant in node_descendants:
                            if descendant.is_leaf():
                                seq_id = number_to_id[descendant.taxon.label]
                                seq_time = float(re.search('_[0-9]+', seq_id).group().replace('_', ''))

                                if seq_time == max_sampling_time:
                                    n_descendants_last_time += 1

                                else:
                                    n_descendants_other_times += 1

                        if n_descendants_last_time > 0 and n_descendants_other_times > 0:
                            branch_in_trunk = 1
                        else:
                            branch_in_trunk = 0

                    # Branch time (duration):
                    branch_time = node.annotations.get_value('length')

                    # Branch length (exp subs per site)
                    branch_exp_subs = float(node.annotations.get_value('rate')) * float(branch_time)



                    # ================================= SIMULATE 100 descendant sequences per model ============================
                    simulated_sequences = {}
                    simulated_sequences['S5FMut_S5FTrans'] = []
                    simulated_sequences['uniformMut_S5FTrans'] = []
                    simulated_sequences['CPMut_S5FTrans'] = []

                    n_reps = 100

                    # Dictionaries to store mutability (for whole sequence and aggregated) for replicate sequences under each model
                    # Subdictionaries for mutability of simulated descendant and for descendant with S or NS subs only:
                    mutability_change_WS_S5FMut_S5FTrans = {'total':{'mean_S5F': []},'syn_only': {'mean_S5F': []},
                                                     'nonsyn_only': {'mean_S5F': []}}

                    mutability_change_WS_uniformMut_S5FTrans = {'total': {'mean_S5F': []}, 'syn_only': {'mean_S5F': []},
                                                            'nonsyn_only': {'mean_S5F': []}}

                    mutability_change_WS_CPMut_S5FTrans = {'total': {'mean_S5F': []}, 'syn_only': {'mean_S5F': []},
                                                                'nonsyn_only': {'mean_S5F': []}}

                    # Get number of randomizable non-synonymous codon changes using the simulation function
                    n_nonsyn_randomizable = randomize_sequence(parent_node_sequence,node_sequence,'S5F','S5F', partition_points)[1]

                    for _ in range(n_reps):

                        # Simulate sequences
                        S5F_S5F_sequence = randomize_sequence(parent_node_sequence,node_sequence,'S5F','S5F',
                                                              partition_points)[0]
                        uniform_S5F_sequence = randomize_sequence(parent_node_sequence,node_sequence,'uniform','S5F',
                                                                  partition_points)[0]
                        CP_S5F_sequence = randomize_sequence(parent_node_sequence,node_sequence,codon_position_rates,
                                                             'S5F', partition_points)[0]

                        # Compute whole-sequence S5F mutability changes in simulated sequences (rel. to obs. parent):
                        S5F_S5F_mutability_change_syn = sequence_differences(parent_node_sequence, S5F_S5F_sequence,partition_points)[
                            'mutability_change_syn']
                        S5F_S5F_mutability_change_nonsyn = sequence_differences(parent_node_sequence, S5F_S5F_sequence,partition_points)[
                            'mutability_change_nonsyn']
                        S5F_S5F_mutability_change_total = S5F_S5F_mutability_change_syn + S5F_S5F_mutability_change_nonsyn

                        uniform_S5F_mutability_change_syn = sequence_differences(parent_node_sequence, uniform_S5F_sequence,partition_points)[
                            'mutability_change_syn']
                        uniform_S5F_mutability_change_nonsyn = sequence_differences(parent_node_sequence, uniform_S5F_sequence,partition_points)[
                            'mutability_change_nonsyn']
                        uniform_S5F_mutability_change_total = uniform_S5F_mutability_change_syn + uniform_S5F_mutability_change_nonsyn

                        CP_S5F_mutability_change_syn = sequence_differences(parent_node_sequence, CP_S5F_sequence,partition_points)[
                            'mutability_change_syn']
                        CP_S5F_mutability_change_nonsyn = sequence_differences(parent_node_sequence, CP_S5F_sequence,partition_points)[
                            'mutability_change_nonsyn']
                        CP_S5F_mutability_change_total = CP_S5F_mutability_change_syn + CP_S5F_mutability_change_nonsyn

                        # Store values in dictionaries containing results for all replicate sequences
                        #for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                        for metric in ['mean_S5F']:

                            mutability_change_WS_S5FMut_S5FTrans['total'][metric].append(S5F_S5F_mutability_change_total)

                            mutability_change_WS_uniformMut_S5FTrans['total'][metric].append(uniform_S5F_mutability_change_total)

                            mutability_change_WS_CPMut_S5FTrans['total'][metric].append(CP_S5F_mutability_change_total)


                            mutability_change_WS_S5FMut_S5FTrans['syn_only'][metric].append(S5F_S5F_mutability_change_syn)

                            mutability_change_WS_uniformMut_S5FTrans['syn_only'][metric].append(uniform_S5F_mutability_change_syn)

                            mutability_change_WS_CPMut_S5FTrans['syn_only'][metric].append(CP_S5F_mutability_change_syn)


                            mutability_change_WS_S5FMut_S5FTrans['nonsyn_only'][metric].append(S5F_S5F_mutability_change_nonsyn)

                            mutability_change_WS_uniformMut_S5FTrans['nonsyn_only'][metric].append(uniform_S5F_mutability_change_nonsyn)

                            mutability_change_WS_CPMut_S5FTrans['nonsyn_only'][metric].append(CP_S5F_mutability_change_nonsyn)

                    # ======================================= OUTPUT RESULTS ===================================================
                    # Follow exact same order as variable names in output_header!
                    new_line_base = str(parent_node_number) + ',' + str(node_number) + ','
                    new_line_base += str(parent_to_root_distance) + ',' + str(parent_time) + ',' + branch_time + ','
                    new_line_base += str(branch_exp_subs) + ',' + str(n_nt_diffs) + ','
                    new_line_base += str(n_NonSyn_diffs) + ',' + str(n_NonSyn_diffs_1nt) + ','
                    new_line_base += str(n_NonSyn_diffs_2nt) + ',' + str(n_NonSyn_diffs_3nt) + ','
                    new_line_base += str(n_Syn_diffs) + ',' + str(n_Syn_diffs_1nt) + ',' + str(n_Syn_diffs_2nt) + ','
                    new_line_base += str(n_Syn_diffs_3nt) + ',' + str(n_nonsyn_randomizable) + ','
                    new_line_base += str(n_motifs_with_syn_nonsyn_changes) + ','

                    new_line_base += str(node.is_leaf() * 1) + ',' + str(branch_in_trunk) + ','
                    new_line_base += str(n_descendants_last_time) + ','

                    new_line_obs = deepcopy(new_line_base)

                    # Observed mutability results
                    #for metric in ['mean_S5F','mean_7M','HS','CS','OHS']:
                    for metric in ['mean_S5F']:
                        # ------ 'observed mutabilities' ------
                        # Add parent mutability
                        new_line_obs += str(parent_node_mutability_WS[1][metric]) + ','
                        new_line_obs += str(parent_node_mutability_aggregated['FR_mutability'][metric]) + ','
                        new_line_obs += str(parent_node_mutability_aggregated['CDR_mutability'][metric]) + ','

                        # Add descendant mutability:
                        new_line_obs += str(node_mutability_WS[1][metric]) + ','
                        new_line_obs += str(node_mutability_aggregated['FR_mutability'][metric]) + ','
                        new_line_obs += str(node_mutability_aggregated['CDR_mutability'][metric]) + ','

                        # Add total mutability change:
                        new_line_obs += str(node_mutability_WS[1][metric] - parent_node_mutability_WS[1][metric]) + ','

                        # Add mutability change due to synonymous changes:
                        new_line_obs += str(S5F_WS_change_syn) + ','

                        # Add mutability change due to non-synonymous changes:
                        new_line_obs += str(S5F_WS_change_nonsyn) + ','
                        
                        # Mutability changes in FRs
                        new_line_obs += str(S5F_change_syn_FR) + ','
                        new_line_obs += str(S5F_change_nonsyn_FR) + ','
                        
                        # Mutability changes in CDRs
                        new_line_obs += str(S5F_change_syn_CDR) + ','
                        new_line_obs += str(S5F_change_nonsyn_CDR) + ','

                    new_line_obs = new_line_obs[:-1]
                    new_line_obs += '\n'

                    output_file_obs.write(new_line_obs)

                    # Results for simulations
                    for i in range(n_reps):
                        new_line_sim = str(i+1) + ',' + new_line_base
                        #for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                        for metric in ['mean_S5F']:
                            # ------ Results from the S5F-mutability/S5F-transitions model ------
                            # Add changes in mutability under S5F-S5F model:
                            new_line_sim += str(mutability_change_WS_S5FMut_S5FTrans['total'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_S5FMut_S5FTrans['syn_only'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_S5FMut_S5FTrans['nonsyn_only'][metric][i]) + ','

                            # Add changes in mutability under uniform-S5F model:
                            new_line_sim += str(mutability_change_WS_uniformMut_S5FTrans['total'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_uniformMut_S5FTrans['syn_only'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_uniformMut_S5FTrans['nonsyn_only'][metric][i]) + ','

                            # Add changes in mutability under codon-position-S5F model:
                            new_line_sim += str(mutability_change_WS_CPMut_S5FTrans['total'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_CPMut_S5FTrans['syn_only'][metric][i]) + ','
                            new_line_sim += str(mutability_change_WS_CPMut_S5FTrans['nonsyn_only'][metric][i]) + ','

                        new_line_sim = new_line_sim[:-1]
                        new_line_sim += '\n'

                        output_file_sim.write(new_line_sim)


if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)

