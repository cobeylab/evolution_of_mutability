"""Takes each branch in an MCC tree, computes its change in mutability and compares it to changes in mutability given randomized descendant sequences under different models (with fixed amino acid sequence.
"""
import sys
# Import mutability functions from analyses/mutability folder
sys.path.insert(0, '../mutability/')
# Import partition_points dictionary from rate_correlations folder
sys.path.insert(0, '../rate_correlations/')

# Import genetic code from simulations folder:
sys.path.insert(0, '../simulations/modules/')

from mutation_functions import sequence_differences, randomize_sequence_constrained, randomize_sequence_unconstrained

from mutability_function import seq_mutability, aggregated_mutability, S5F
from partition_points import partition_points_dic
from dendropy import Tree
from numpy import log

import re
from copy import deepcopy

# Import function to read observed sequences from XML file
sys.path.insert(0, '../contrasts/')
from contrasts_functions import get_mutability_from_XML


def main(argv):
    # Chain id, e.g. CH103_con_run1a, VRC01_L08_log_run1b, scenario2a_rep1
    chain_id = str(argv[1])

    # ================== RETRIEVE DIRECTORIES, OBS. SEQUENCES AND FR/CDR PARTITION POINTS FROM CHAIN ID ================
    # If chain is from an observed clone:
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
    output_file_path_simulated_constrained = output_directory + chain_id + '_simulated_MCC_constrained.csv'
    output_file_path_simulated_unconstrained = output_directory + chain_id + '_simulated_MCC_unconstrained.csv'

    output_file_path_aa_transitions_obs = output_directory + chain_id + '_aa_transitions_obs.csv'
    output_file_path_aa_transitions_unconstrained = output_directory + chain_id + '_aa_transitions_unconstrained.csv'


    # ========================== READ TIP SEQUENCES FROM XML FILE AND COMPUTE THEIR MUTABILITY==========================
    # ------------------------------- (Tip sequences are not annotated on the BEAST trees) -----------------------------
    mutability_of_obs_seqs = get_mutability_from_XML(xml_file_path, partition_points)

    # Dictionaries with mutability for each tip sequence:
    obs_sequence = mutability_of_obs_seqs['sequences']
    mutability_WS_tips = mutability_of_obs_seqs['whole_sequence']
    mutability_aggregated_tips = mutability_of_obs_seqs['aggregated_by_region']

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

    # ==================================== OPEN TREE FILE AND OUTPUT FILES =============================================
    with open(MCC_tree_file_path, 'r') as tree_file, open(output_file_path_observed, 'w') as output_file_obs, open(
            output_file_path_simulated_constrained, 'w') as output_file_sim_constrained, open(
        output_file_path_simulated_unconstrained, 'w') as output_file_sim_unconstrained:

        # Header shared by observed and simulated output files:
        output_header_base = 'parent,child,parent_distance_to_root,parent_time_to_root,branch_time,branch_exp_subs,'
        output_header_base += 'n_NT_diffs,n_NonSyn_diffs,n_NonSyn_diffs_1nt,n_NonSyn_diffs_2nt,n_NonSyn_diffs_3nt,'
        output_header_base += 'n_Syn_diffs,n_Syn_diffs_1nt,n_Syn_diffs_2nt,n_Syn_diffs_3nt,'
        output_header_base += 'n_nonsyn_changes_randomizable,n_motifs_with_syn_nonsyn_changes,branch_is_terminal,'
        output_header_base += 'branch_in_trunk,n_descendants_last_time'

        output_header_obs = deepcopy(output_header_base)
        # Header for output file with observed mutability
        # Ignoring other metrics and focusing on S5F:
        for metric in ['S5F','logS5F']:
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
        for region in ['WS', 'FR', 'CDR']:
            for metric in ['S5F', 'logS5F']:
                for mutability_model in ['S5F', 'uniform', 'CP']:
                    # Mutability changes under mutation model / S5F-transitions model
                    output_header_sim += ',' + metric + '_' + region + '_change_' + mutability_model + 'Mut_S5FTrans_total'
                    output_header_sim += ',' + metric + '_' + region + '_change_' + mutability_model + 'Mut_S5FTrans_syn'
                    output_header_sim += ',' + metric + '_' + region + '_change_' + mutability_model + 'Mut_S5FTrans_nonsyn'

        output_header_sim += '\n'

        output_file_sim_constrained.write(output_header_sim)
        output_file_sim_unconstrained.write(output_header_sim)

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
                    new_block += value.group() + ','

            #Remove extra comma at the end:
            new_block = new_block[0:(len(new_block) - 1)]

            new_block += ']'

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

        # AA transition results for entire observed tree (dictionary with one AA transition dic. per node)
        AA_transitions_obs = {}
        AA_meanLogS5F_changes_obs = {}

        # AA transition results for randomizations (Dictionary with a list of replicate AA trans. dicts. per model per node)
        AA_transitions_sim_unconstrained = {}
        AA_meanLogS5F_changes_sim_unconstrained = {}

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
                        focal_node = focal_node.parent_node

                    assert sum(parent_to_root_time_lengths) == parent_time

                    # Find parent's distance to the root by multiplying molecular clock rates and time lengths, then summing
                    parent_to_root_distance = sum([parent_to_root_time_lengths[i] * parent_to_root_rates[i] for i in
                                                   range(len(parent_to_root_rates))])

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
                    
                    # Get changes in the mean log S5F
                    logS5F_WS_change_syn = differences['log_mutability_change_syn']
                    logS5F_WS_change_nonsyn = differences['log_mutability_change_nonsyn']

                    logS5F_change_syn_FR = differences['log_mutability_change_syn_FR']
                    logS5F_change_nonsyn_FR = differences['log_mutability_change_nonsyn_FR']

                    logS5F_change_syn_CDR = differences['log_mutability_change_syn_CDR']
                    logS5F_change_nonsyn_CDR = differences['log_mutability_change_nonsyn_CDR']


                    # Count number of motifs where both synonymous and nonsynonymous changes happened
                    n_motifs_with_syn_nonsyn_changes = differences['motifs_with_syn_nonsyn']

                    # Count number of non-synonymous codon changes that can be randomized in constrained simulations
                    n_nonsyn_randomizable = randomize_sequence_constrained(parent_node_sequence, node_sequence,
                                                                           'S5F','S5F', partition_points)[1]

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

                    # Get branch amino acid transitions and their contributions to changes in mean log S5F mutability
                    AA_transitions_obs[node_number] =  differences['aa_transitions']
                    AA_meanLogS5F_changes_obs[node_number] = differences['aa_trans_logmut_changes']


                    # ============================== SIMULATE 100 descendant sequences per model =======================
                    AA_transitions_sim_unconstrained[node_number] = {'S5F':[],'uniform':[], 'CP': []}
                    AA_meanLogS5F_changes_sim_unconstrained[node_number] = {'S5F':[],'uniform':[], 'CP': []}

                    n_reps = 100

                    # Dictionary with S5F changes by simulation type (constrained / unconstrained) and mutability model
                    S5F_changes = {'constrained': {'S5F': {}, 'uniform': {}, 'CP': {}},
                                   'unconstrained': {'S5F': {}, 'uniform': {}, 'CP': {}}
                                   }
                    # Dictionary with mean log S5F changes by simulation type (constrained / unconstrained) and mutability model
                    logS5F_changes = {'constrained': {'S5F': {}, 'uniform': {}, 'CP': {}},
                                   'unconstrained': {'S5F': {}, 'uniform': {}, 'CP': {}}
                                   }

                    for key in S5F_changes.keys():
                        for subkey in S5F_changes[key].keys():
                            S5F_changes[key][subkey] = {'WS':{'total':[], 'syn': [], 'nonsyn': []},
                                                        'FR':{'total':[], 'syn': [], 'nonsyn': []},
                                                        'CDR':{'total':[], 'syn': [], 'nonsyn': []}
                                                        }
                            
                    for key in logS5F_changes.keys():
                        for subkey in logS5F_changes[key].keys():
                            logS5F_changes[key][subkey] = {'WS':{'total':[], 'syn': [], 'nonsyn': []},
                                                        'FR':{'total':[], 'syn': [], 'nonsyn': []},
                                                        'CDR':{'total':[], 'syn': [], 'nonsyn': []}
                                                        }

                    for _ in range(n_reps):
                        # Generate one simulated sequence per mutability model / constrained vs. unconstrained
                        for simulation_type in ['constrained','unconstrained']:

                            if simulation_type == 'constrained':
                                randomization_function = randomize_sequence_constrained
                            elif simulation_type == 'unconstrained':
                                randomization_function = randomize_sequence_unconstrained

                            for mutability_model in ['S5F','uniform', 'CP']:

                                mut_model = mutability_model
                                if mut_model == 'CP':
                                    mut_model = codon_position_rates

                                sim_seq = randomization_function(
                                    parent_node_sequence, node_sequence, mut_model, 'S5F', partition_points)

                                if simulation_type == 'constrained':
                                    sim_seq = sim_seq[0]

                                # Compute changes in S5F mutability between simulated sequences and parent sequence:
                                diffs_from_parent = sequence_differences(parent_node_sequence, sim_seq, partition_points)

                                # Record simulated AA trans. on this branch and their contribution to mean logS5F change
                                AA_transitions_sim_unconstrained[node_number][mutability_model].append(diffs_from_parent['aa_transitions'])
                                AA_meanLogS5F_changes_sim_unconstrained[node_number][mutability_model].append(diffs_from_parent['aa_trans_logmut_changes'])

                                # Store changes in dictionary
                                for region in ['WS','FR','CDR']:

                                    # WS is omitted in diffs_from_parent dictionary keys
                                    region_id = region
                                    if region_id == 'WS':
                                        region_id = ''
                                    else:
                                        region_id = '_' + region_id

                                    syn_change = diffs_from_parent['mutability_change_syn' + region_id]
                                    nonsyn_change = diffs_from_parent['mutability_change_nonsyn' + region_id]
                                    total_change = syn_change + nonsyn_change
                                    
                                    log_syn_change = diffs_from_parent['log_mutability_change_syn' + region_id]
                                    log_nonsyn_change = diffs_from_parent['log_mutability_change_nonsyn' + region_id]
                                    log_total_change = log_syn_change + log_nonsyn_change

                                    S5F_changes[simulation_type][mutability_model][region]['syn'].append(syn_change)
                                    S5F_changes[simulation_type][mutability_model][region]['nonsyn'].append(
                                        nonsyn_change)
                                    S5F_changes[simulation_type][mutability_model][region]['total'].append(
                                        total_change)

                                    logS5F_changes[simulation_type][mutability_model][region]['syn'].append(
                                        log_syn_change)
                                    logS5F_changes[simulation_type][mutability_model][region]['nonsyn'].append(
                                        log_nonsyn_change)
                                    logS5F_changes[simulation_type][mutability_model][region]['total'].append(
                                        log_total_change)
                                    

                    # ===================================== OUTPUT RESULTS =============================================
                    # Follow exact same order as variable names in output_header
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
                    
                        # ------ 'observed mutabilities' ------
                    # Add parent mean S5F
                    new_line_obs += str(parent_node_mutability_WS[1]['mean_S5F']) + ','
                    new_line_obs += str(parent_node_mutability_aggregated['FR_mutability']['mean_S5F']) + ','
                    new_line_obs += str(parent_node_mutability_aggregated['CDR_mutability']['mean_S5F']) + ','

                    # Add descendant mean S5F:
                    new_line_obs += str(node_mutability_WS[1]['mean_S5F']) + ','
                    new_line_obs += str(node_mutability_aggregated['FR_mutability']['mean_S5F']) + ','
                    new_line_obs += str(node_mutability_aggregated['CDR_mutability']['mean_S5F']) + ','

                    # Add total change in mean S5F:
                    new_line_obs += str(node_mutability_WS[1]['mean_S5F'] - parent_node_mutability_WS[1]['mean_S5F']) + ','

                    # Add change in mean S5F due to synonymous changes:
                    new_line_obs += str(S5F_WS_change_syn) + ','

                    # Add change in mean S5F due to non-synonymous changes:
                    new_line_obs += str(S5F_WS_change_nonsyn) + ','

                    # Mean S5F changes in FRs
                    new_line_obs += str(S5F_change_syn_FR) + ','
                    new_line_obs += str(S5F_change_nonsyn_FR) + ','

                    # Mean S5F changes in CDRs
                    new_line_obs += str(S5F_change_syn_CDR) + ','
                    new_line_obs += str(S5F_change_nonsyn_CDR) + ','

                    # ------- Add results for the mean log S5F, by taking the log of the geometric mean S5F
                    
                    # Add parent mean log S5F
                    new_line_obs += str(log(parent_node_mutability_WS[1]['geom_mean_S5F'])) + ','
                    new_line_obs += str(log(parent_node_mutability_aggregated['FR_mutability']['geom_mean_S5F'])) + ','
                    new_line_obs += str(log(parent_node_mutability_aggregated['CDR_mutability']['geom_mean_S5F'])) + ','

                    # Add descendant mean log S5F
                    new_line_obs += str(log(node_mutability_WS[1]['geom_mean_S5F'])) + ','
                    new_line_obs += str(log(node_mutability_aggregated['FR_mutability']['geom_mean_S5F'])) + ','
                    new_line_obs += str(log(node_mutability_aggregated['CDR_mutability']['geom_mean_S5F'])) + ','

                    # Add total change in mean log S5F:
                    new_line_obs += str(log(node_mutability_WS[1]['geom_mean_S5F']) - log(parent_node_mutability_WS[1]['geom_mean_S5F'])) + ','

                    # Add change in mean log S5F due to synonymous changes:
                    new_line_obs += str(logS5F_WS_change_syn) + ','
                    
                    # Add change in mean log S5F due to non-synonymous changes:
                    new_line_obs += str(logS5F_WS_change_nonsyn) + ','

                    # Mean log S5F changes in FRs
                    new_line_obs += str(logS5F_change_syn_FR) + ','
                    new_line_obs += str(logS5F_change_nonsyn_FR) + ','

                    # Mean log S5F changes in CDRs
                    new_line_obs += str(logS5F_change_syn_CDR) + ','
                    new_line_obs += str(logS5F_change_nonsyn_CDR)

                    new_line_obs += '\n'

                    output_file_obs.write(new_line_obs)

                    # Results for simulations, output to separate files for constrained / unconstrained

                    for simulation_type in ['constrained', 'unconstrained']:

                        if simulation_type == 'constrained':
                            results_file = output_file_sim_constrained
                        elif simulation_type == 'unconstrained':
                            results_file = output_file_sim_unconstrained

                        for i in range(n_reps):
                            new_line_sim = str(i + 1) + ',' + new_line_base[:-1]
                            for region in ['WS', 'FR', 'CDR']:
                                    for metric in ['S5F','logS5F']:
                                        if metric == 'S5F':
                                            change_dictionary = S5F_changes
                                        elif metric == 'logS5F':
                                            change_dictionary = logS5F_changes

                                        for mutability_model in ['S5F','uniform','CP']:

                                            total_change = change_dictionary[simulation_type][mutability_model][region]['total'][i]
                                            syn_change = change_dictionary[simulation_type][mutability_model][region]['syn'][i]
                                            nonsyn_change = change_dictionary[simulation_type][mutability_model][region]['nonsyn'][i]

                                            new_line_sim += ',' + str(total_change) + ',' + str(syn_change) + ',' + \
                                                            str(nonsyn_change)
                            new_line_sim += '\n'

                            results_file.write(new_line_sim)

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)

