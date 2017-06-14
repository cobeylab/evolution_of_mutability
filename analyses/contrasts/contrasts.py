#!/usr/bin/python

"""Takes each tree from a NEXUS file (output from BEAST). For each tree, computes mutability contrasts and rate contrasts. Mutability contrasts are defined as the difference in mutability between consecutive nodes (child - parent). Rate contrasts are defined as the difference in rate between two consecutive branches. For node N0 with child N1 and grandchild N1a, the rate contrast is defined as the rate on the branch connecting N1 to N1a minus the rate on the branch connecting N0 to N1.
"""
import sys
import re
from dendropy import Tree
from numpy import random

# Import mutability functions and partition points from analyses/mutability folder
sys.path.insert(0, '../mutability/')
from mutability_function import seq_mutability, aggregated_mutability
from partition_points import partition_points_dic

# Import function for partitioning S5F changes into syn. and non-syn. components:
sys.path.insert(0, '../S_NS_mutability_changes/')
from mutation_functions import sequence_differences


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
        trees_file_path = chain_directory + chain_id + '_tree_sample.trees'

        output_directory = '../../results/contrasts/observed_lineages/' + clone + '_' + prior + '/'
        mutability_output_file_path = output_directory + chain_id + '_mutability_contrasts.csv'
        rate_output_file_path = output_directory + chain_id + '_rate_contrasts.csv'

        xml_file_path = '../../analyses/BEAST/observed_lineages/' + clone + '_' + prior + '/'
        xml_file_path += chain_id[0:len(chain_id) - 1] + '.xml'

    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        trees_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        trees_file_path += chain_id + '_tree_sample.trees'

        output_directory = '../../results/contrasts/simulated_alignments/' + chain_id + '/'
        mutability_output_file_path = output_directory + chain_id + '_mutability_contrasts.csv'
        rate_output_file_path = output_directory + chain_id + '_rate_contrasts.csv'

        xml_file_path = '../../analyses/BEAST/simulated_alignments/' + chain_id + '/' + chain_id + '.xml'

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

    # ======================================== OPEN TREES FILE AND OUTPUT FILES ========================================

    with open(trees_file_path, 'r') as trees_file, open(rate_output_file_path, 'w') as rate_output_file, \
            open(mutability_output_file_path, 'w') as mutability_output_file:

        # Write header to mutability contrasts output file:
        mutability_output_header = 'tree,parent,child,branch_is_terminal,branch_in_trunk,n_descendants_last_time'
        for region in ['WS', 'FR', 'CDR']:
            for metric in ['S5F','7M','HS','CS','OHS']:
                mutability_output_header += ',' + metric + '_' + region + '_contrast'

        for region in ['WS', 'FR', 'CDR']:
            for metric in ['S5F','7M','HS','CS','OHS']:
                mutability_output_header += ',' + metric + '_' + region + '_parent'

        # For S5F only, changes partitioned into syn. and non-syn.
        mutability_output_header += ',S5F_WS_contrast_syn,S5F_WS_contrast_nonsyn,'
        mutability_output_header += 'S5F_FR_contrast_syn,S5F_FR_contrast_nonsyn,'
        mutability_output_header += 'S5F_CDR_contrast_syn,S5F_CDR_contrast_nonsyn,'

        mutability_output_header += 'parent_distance,parent_distance_RC,parent_time\n'

        mutability_output_file.write(mutability_output_header)

        # Write header to rate contrasts output file:
        rate_output_header = 'tree,parent,child,grandchild,grandchild_branch_is_terminal,total_rate_contrast,total_rate_RC_contrast,S_rate_contrast,N_rate_contrast,'
        rate_output_header += 'total_rate_to_child,total_rate_RC_to_child,S_rate_to_child,N_rate_to_child,'
        rate_output_header += 'child_distance,child_distance_RC,child_time\n'

        rate_output_file.write(rate_output_header)

        # Dictionary linking BEAST numbering (from NEXUS file) to taxon labels:
        number_to_id = {}

        # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
        n_semicolons = 0
        while n_semicolons < 6:
            line = trees_file.next()

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

        line_has_tree = True

        # From sequence ids, find maximum sampling time
        sampling_times = [float(re.search('_[0-9]+',seq_id).group().replace('_','')) for seq_id in number_to_id.values()]
        max_sampling_time = max(sampling_times)

        # ============================================ ANALYZING TREES =================================================

        # For each tree:
        while line_has_tree:
            line = trees_file.next()
            # Search for a tree header in the line
            tree_header = re.search(r'.*\[\&R\]', line)

            # If there's a tree in this line (i.e. if there is a header):
            if tree_header is not None:
                tree_string = line
                tree_header = tree_header.group()
                state = re.search(r'STATE_[0-9]*', tree_header).group().replace('STATE_', '')
                state = int(state)

                print 'STATE ' + str(state)

                # Remove tree header from tree_string before reading it with dendropy
                tree_string = tree_string.replace(tree_header, '')
                tree_string = tree_string.replace('\n', '')

                # Read tree string as a tree object using dendropy:
                tree = Tree.get_from_string(tree_string, schema='newick')

                # Label all nodes according to their index in tree.nodes():
                for i in range(len(tree.nodes())):
                    tree.nodes()[i].label = str(i)
                # Notice that this is different from the numbering BEAST does on the tips


                # Get number of sites from codon position 1 annotation for node 0
                assert tree.seed_node.annotations[0].name.find('CP1') > -1, 'Could not read internal node sequences from node 0'
                n_sites = 3 * len(tree.nodes()[0].annotations[0].value)

                # Alignment identifier in codon position annotations (e.g. CH103_final_alignment.CP1)
                alignment_id = [annotation.name for annotation in tree.seed_node.annotations if annotation.name.find('CP1') > 1]
                alignment_id = alignment_id[0].replace('.CP1', '')

                # For each node...
                for node in tree.nodes():
                    # Find node number (index in tree.nodes())
                    node_number = node.label

                    # ================================ COMPUTING MUTABILITY CONTRASTS ==================================

                    # If node is not the root, find changes in mutability (contrasts) from parent node
                    if node.parent_node is not None:
                        #print node_number

                        # Ignoring a weird VRC26 sequence whose CDR2 is entirely gaps:
                        if node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'KJ134124_119':
                            print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'
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
                            parent_to_root_lengths_RC = [0.0]
                            focal_node = parent_node

                            # Move down the tree to the root node
                            while focal_node is not tree.nodes()[0]:
                                parent_to_root_time_lengths.append(focal_node.edge_length)
                                parent_to_root_rates.append(float(focal_node.annotations.get_value('rate')))
                                RC_length = float(focal_node.annotations.get_value('S')) + float(focal_node.annotations.get_value('N'))
                                RC_length = float(RC_length) / n_sites
                                parent_to_root_lengths_RC.append(RC_length)

                                focal_node = focal_node.parent_node

                            assert sum(parent_to_root_time_lengths) == parent_time

                            # Find parent's distance to the root by multiplying molecular clock rates and time lengths, then summing
                            parent_to_root_distance = sum([parent_to_root_time_lengths[i] * parent_to_root_rates[i] for i in range(len(parent_to_root_rates))])

                            # Find parent's robust counting distance by summing robust counting branch lengths
                            parent_to_root_distance_RC = sum(parent_to_root_lengths_RC)

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

                                node_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
                                node_sequence = ''.join(node_sequence)

                                # Compute mutability
                                node_mutability_WS = seq_mutability(node_sequence)
                                node_mutability_aggregated = aggregated_mutability(node_sequence, partition_points)

                            # Find mutability of parent node:

                            parent_node_CP1 = parent_node.annotations.get_value(alignment_id + '.CP1')
                            parent_node_CP2 = parent_node.annotations.get_value(alignment_id + '.CP2')
                            parent_node_CP3 = parent_node.annotations.get_value(alignment_id + '.CP3')

                            parent_node_sequence = [parent_node_CP1[i] + parent_node_CP2[i] + parent_node_CP3[i] for i in range(len(parent_node_CP1))]
                            parent_node_sequence = ''.join(parent_node_sequence)

                            # Compute mutability of parent node
                            parent_node_mutability_WS = seq_mutability(parent_node_sequence)
                            parent_node_mutability_aggregated = aggregated_mutability(parent_node_sequence, partition_points)


                            # Determine if parent->node branch is part of the trunk, and count branch's n. descendants
                            # Branch is in trunk if it has desc. at last time but not if it only has desc. from last time
                            # Lemey et al. 2007 PlosCompBio

                            # (n. descendants at the last time point)

                            #branch_in_trunk = 0
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
                                        seq_time = float(re.search('_[0-9]+', seq_id).group().replace('_',''))

                                        if seq_time == max_sampling_time:
                                            n_descendants_last_time += 1

                                        else:
                                            n_descendants_other_times += 1

                                if n_descendants_last_time > 0 and n_descendants_other_times > 0:
                                    branch_in_trunk = 1
                                else:
                                    branch_in_trunk = 0


                            # --------------------------- Write wew line of mutability output file -------------------------
                            # Adding results to output line in the order specified in mutability_output_header

                            mutability_line = str(state) + ',' + str(parent_node_number) + ',' + str(node_number) + ','

                            mutability_line += str(node.is_leaf() * 1) + ',' + str(branch_in_trunk) + ','
                            mutability_line += str(n_descendants_last_time) + ','

                            # The following loops MUST BE in this exact order

                            # Add whole-sequence mutability contrasts:
                            for metric in ['mean_S5F','mean_7M','HS','CS','OHS']:
                                node_mutability = node_mutability_WS[1][metric]
                                parent_mutability = parent_node_mutability_WS[1][metric]
                                contrast = node_mutability - parent_mutability
                                mutability_line += str(contrast) + ','

                            # Add FR and CDR aggregated mutability contrasts:
                            for region in ['FR','CDR']:
                                for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                                    node_mutability = node_mutability_aggregated[region + '_mutability'][metric]
                                    parent_mutability = parent_node_mutability_aggregated[region + '_mutability'][metric]
                                    contrast = node_mutability - parent_mutability
                                    mutability_line += str(contrast) + ','

                            # Add parent whole-sequence mutability
                            for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                                parent_mutability = parent_node_mutability_WS[1][metric]
                                mutability_line += str(parent_mutability) + ','

                            # Add parent FR and CDR aggregated mutability:
                            for region in ['FR','CDR']:
                                for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                                    parent_mutability = parent_node_mutability_aggregated[region + '_mutability'][metric]
                                    mutability_line += str(parent_mutability) + ','


                            # For S5F only, add changes partitioned into syn. and non-syn. changes
                            S5F_contrast_WS_syn = sequence_differences(parent_node_sequence,node_sequence,
                                                                       partition_points)['mutability_change_syn']
                            S5F_contrast_WS_nonsyn = sequence_differences(parent_node_sequence,node_sequence,
                                                                          partition_points)['mutability_change_nonsyn']

                            S5F_contrast_FR_syn = sequence_differences(parent_node_sequence, node_sequence,
                                                                       partition_points)['mutability_change_syn_FR']
                            S5F_contrast_FR_nonsyn = sequence_differences(parent_node_sequence, node_sequence,
                                                                       partition_points)['mutability_change_nonsyn_FR']

                            S5F_contrast_CDR_syn = sequence_differences(parent_node_sequence, node_sequence,
                                                                       partition_points)['mutability_change_syn_CDR']
                            S5F_contrast_CDR_nonsyn = sequence_differences(parent_node_sequence, node_sequence,
                                                                      partition_points)['mutability_change_nonsyn_CDR']

                            mutability_line += str(S5F_contrast_WS_syn) + ',' + str(S5F_contrast_WS_nonsyn) + ','
                            mutability_line += str(S5F_contrast_FR_syn) + ',' + str(S5F_contrast_FR_nonsyn) + ','
                            mutability_line += str(S5F_contrast_CDR_syn) + ',' + str(S5F_contrast_CDR_nonsyn) + ','

                            mutability_line += ','.join([str(parent_to_root_distance), str(parent_to_root_distance_RC),
                                                         str(parent_time)])

                            mutability_line += '\n'

                            #print 'Mutability: ' + parent_node_number + ' ' + node_number #+ ':::::' + mutability_line
                            mutability_output_file.write(mutability_line)

                    # =================================== COMPUTING RATE CONTRASTS =====================================

                    # If node has at least one grandchild, compute rate(Child -> Grdchild) - rate (Node -> Child)
                    # If node has children...
                    if len(node.child_nodes()) > 0:
                        #print node_number
                        # For each child:
                        for child_node in node.child_nodes():

                            # If there's at least one grandchild:
                            if len(child_node.child_nodes()) > 0:

                                # Find child node dendropy number:
                                child_node_number = child_node.label

                                # ------------------ Find time and distance from root for the child node ---------------

                                # note tree is scaled by time so the command below gets time, not distance)
                                child_node_time = child_node.distance_from_root()

                                # To find child node's distance from root, first find time-lengths and rates for all branches connecting it to the root
                                # (robust-counting per-site distances can be found directly from 'S' and 'N' without reference to time lengths)
                                child_node_to_root_time_lengths = [0.0]
                                child_node_to_root_rates = [0.0]
                                child_node_to_root_lengths_RC = [0.0]

                                focal_node = child_node
                                while focal_node is not tree.nodes()[0]:
                                    child_node_to_root_time_lengths.append(focal_node.edge_length)
                                    child_node_to_root_rates.append(float(focal_node.annotations.get_value('rate')))
                                    RC_length = float(focal_node.annotations.get_value('S')) + float(focal_node.annotations.get_value('N'))
                                    RC_length = float(RC_length) / n_sites
                                    child_node_to_root_lengths_RC.append(RC_length)

                                    focal_node = focal_node.parent_node

                                # Check that summing times down to the root matches the .distance_from_root() dendropy method
                                assert sum(child_node_to_root_time_lengths) == child_node_time

                                # Find node's distance to the root by multiplying molecular clock rates and time lengths, then summing
                                child_node_to_root_distance = sum([child_node_to_root_time_lengths[i] * child_node_to_root_rates[i] for i in range(len(child_node_to_root_rates))])

                                # Find nodes's robust counting distance by summing robust counting branch lengths
                                child_node_to_root_distance_RC = sum(child_node_to_root_lengths_RC)

                                # ----------------------- Find rates on branch from node to child ----------------------

                                # From child node annotation, get rates on branch from node to child
                                # Get total rate from molecular clock:
                                total_rate_to_child = float(child_node.annotations.get_value('rate'))

                                # Synonymous and Non-synonymous counts from robust counting:
                                S_count_to_child = float(child_node.annotations.get_value('S'))

                                N_count_to_child = float(child_node.annotations.get_value('N'))

                                # Get branch length in time
                                time_length_to_child = child_node.edge_length

                                # Convert counts into rates
                                S_rate_to_child = S_count_to_child / time_length_to_child
                                N_rate_to_child = N_count_to_child / time_length_to_child

                                # Normalize rates by number of sites:
                                S_rate_to_child = S_rate_to_child / n_sites
                                N_rate_to_child = N_rate_to_child / n_sites

                                # From syn. and non-syn. rates, get total rate according to robust counting.
                                total_rate_RC_to_child = S_rate_to_child + N_rate_to_child

                                # =-------------------- For each grandchild, compute rate contrasts --------------------

                                for grandchild_node in child_node.child_nodes():

                                    grandchild_node_number = grandchild_node.label

                                    # From grandchild node annotation, get rates on branch from child to grandchild
                                    # Get total rate from molecular clock:
                                    total_rate_to_grandchild = float(grandchild_node.annotations.get_value('rate'))

                                    # Synonymous and Non-synonymous counts from robust counting:
                                    S_count_to_grandchild = float(grandchild_node.annotations.get_value('S'))
                                    N_count_to_grandchild = float(grandchild_node.annotations.get_value('N'))

                                    # Get branch length in time
                                    time_length_to_grandchild = grandchild_node.edge_length

                                    # Convert counts into rates
                                    S_rate_to_grandchild = S_count_to_grandchild / time_length_to_grandchild
                                    N_rate_to_grandchild = N_count_to_grandchild    / time_length_to_grandchild

                                    # Normalize rates by number of sites:
                                    S_rate_to_grandchild = S_rate_to_grandchild / n_sites
                                    N_rate_to_grandchild = N_rate_to_grandchild / n_sites

                                    # From syn. and non-syn. rates, get total rate according to robust counting.
                                    total_rate_RC_to_grandchild = S_rate_to_grandchild + N_rate_to_grandchild

                                    # Compute rate contrasts:
                                    total_rate_contrast = total_rate_to_grandchild - total_rate_to_child
                                    total_rate_RC_contrast = total_rate_RC_to_grandchild - total_rate_RC_to_child
                                    S_rate_contrast = S_rate_to_grandchild - S_rate_to_child
                                    N_rate_contrast = N_rate_to_grandchild - N_rate_to_child

                                    # ----------------------- Write new line of rate output file -----------------------
                                    # Concatenate results in the order specified by rate_output_header

                                    rate_line = ','.join([str(state), str(node_number),
                                                          str(child_node_number), str(grandchild_node_number)])
                                    rate_line += ',' + str(grandchild_node.is_leaf() * 1) + ','

                                    rate_line += ','.join([str(total_rate_contrast), str(total_rate_RC_contrast),
                                                           str(S_rate_contrast), str(N_rate_contrast)])

                                    rate_line += ','

                                    rate_line += ','.join([str(total_rate_to_child), str(total_rate_RC_to_child),
                                                           str(S_rate_to_child), str(N_rate_to_child),
                                                           str(child_node_to_root_distance),
                                                           str(child_node_to_root_distance_RC),
                                                           str(child_node_time)])
                                    rate_line += '\n'

                                    #print 'Rate: ' + node_number + ' ' + child_node_number + ' ' + grandchild_node_number #+ ':::::' + mutability_line
                                    rate_output_file.write(rate_line)

            else:
                line_has_tree = False

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)