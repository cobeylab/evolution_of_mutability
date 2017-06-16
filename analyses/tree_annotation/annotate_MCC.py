#!/usr/bin/python

"""Annotates the MCC tree for the chain specified by the chain_id in argv[1] with mutability, rates and contrasts
"""

import sys
import re
from dendropy import Tree
from copy import deepcopy

# Import mutability functions and partition points from analyses/mutability folder
sys.path.insert(0, '../mutability/')

from mutability_function import seq_mutability, aggregated_mutability
from partition_points import partition_points_dic

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
        MCC_file_path = chain_directory + chain_id + '_MCC_tree.tree'

        output_directory = '../../results/tree_annotation/'
        output_file_path = output_directory + chain_id + '_annotated_MCC.tree'


        xml_file_path = '../../analyses/BEAST/observed_lineages/' + clone + '_' + prior + '/'
        xml_file_path += chain_id[0:len(chain_id) - 1] + '.xml'

    # If chain is from a simulated scenario, use arbitrary partition points (based on VRC01_13)
    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        MCC_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        MCC_file_path += chain_id + '_MCC_tree.tree'

        output_directory = '../../results/tree_annotation/'
        output_file_path = output_directory + chain_id + '_annotated_MCC.tree'

        xml_file_path = '../../analyses/BEAST/simulated_alignments/' + chain_id + '/' + chain_id + '.xml'

        partition_points = [1, 76, 100, 151, 175, 289, 351]

    # ========================== READ TIP SEQUENCES FROM XML FILE AND COMPUTE THEIR MUTABILITY==========================
    # ------------------------------- (Tip sequences are not annotated on the BEAST trees) -----------------------------

    # Dictionaries with mutability for each tip sequence:
    mutability_WS_tips = {}
    mutability_aggregated_tips = {}

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

            if partition_points is not None:
                mutability_aggregated_tips[taxon_id] = aggregated_mutability(sequence, partition_points)
        else:
            print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'


    # ======================================== OPEN TREES FILE AND OUTPUT FILES ========================================

    with open(MCC_file_path, 'r') as tree_file, open(output_file_path, 'w') as output_file:

        # Dictionary linking BEAST numbering (from NEXUS file) to taxon labels:
        number_to_id = {}

        # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
        n_semicolons = 0
        while n_semicolons < 6:
            line = tree_file.next()

            # Write line to output file:
            output_file.write(line)

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

        # Dendropy doesn't seem to read annotations of BEAST MCC trees correctly
        # Get annotation blocks for all nodes, simplify them and input edited tree string to dendropy

        # Basically find all instances of "[...]" in the newick string:
        annotation_blocks = re.findall(r'\[[^\]]*]',tree_string)
        annotation_blocks = annotation_blocks[1:len(annotation_blocks)]

        # Find alignment identifier (e.g. CH103_final_alignment in CH103_final_alignment.CP1)
        alignment_id = re.search(r'[^,^\[]*.CP1', annotation_blocks[0]).group()
        alignment_id = alignment_id.replace('&', '').replace('.CP1', '')

        # Annotations that will be kept in the edited tree string
        variables = ['rate', 'S', 'N', 'height', 'length', 'b_u_S', 'b_u_N',
                     alignment_id + '.CP1', alignment_id + '.CP2',
                     alignment_id + '.CP3']

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

        # For each node...
        for node in tree.nodes():

            print 'Node number ' + str(node.label)

            # =========================================== ANNOTATIONS ==================================================

            # If node is not the root, find changes in mutability (contrasts) from parent node
            # (Root node will be annotated along with its immediate descendants)
            if node.parent_node is not None:

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
                        RC_length = float(focal_node.annotations.get_value('S')) + float(
                            focal_node.annotations.get_value('N'))
                        RC_length = float(RC_length) / n_sites
                        parent_to_root_lengths_RC.append(RC_length)

                        focal_node = focal_node.parent_node

                    assert sum(parent_to_root_time_lengths) == parent_time

                    # Find node mutability, if it is a tip:
                    if node.is_leaf():
                        # Find taxon number from BEAST tree
                        taxon_number = str(node.taxon).replace("'", '')
                        # Find node taxon id (the original name for the tip sequence):
                        node_id = number_to_id[taxon_number]

                        # Get mutability from tip mutability dictionaries
                        node_mutability_WS = mutability_WS_tips[node_id]
                        node_mutability_aggregated = mutability_aggregated_tips[node_id]
                    else:
                        # If node is not a tip, find its sequence from the tree annotation:
                        node_CP1 = node.annotations.get_value(alignment_id + '.CP1')
                        node_CP2 = node.annotations.get_value(alignment_id + '.CP2')
                        node_CP3 = node.annotations.get_value(alignment_id + '.CP3')

                        # Arbitrarily choosing first listed reconstruction in case of ambiguity
                        if node_CP1.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP1 = node_CP1.split('+')[0]
                        if node_CP2.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP2 = node_CP2.split('+')[0]
                        if node_CP3.find('+') > -1:
                            # If it is, arbitrarily choose the first sequence
                            node_CP3 = node_CP3.split('+')[0]

                        node_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
                        node_sequence = ''.join(node_sequence)

                        # Compute mutability
                        node_mutability_WS = seq_mutability(node_sequence)
                        node_mutability_aggregated = aggregated_mutability(node_sequence, partition_points)

                    # Find sequence of parent node:

                    parent_node_CP1 = parent_node.annotations.get_value(alignment_id + '.CP1')
                    parent_node_CP2 = parent_node.annotations.get_value(alignment_id + '.CP2')
                    parent_node_CP3 = parent_node.annotations.get_value(alignment_id + '.CP3')

                    # Arbitrarily choosing first listed reconstruction in case of ambiguity
                    if parent_node_CP1.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP1 = parent_node_CP1.split('+')[0]
                    if parent_node_CP2.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP2 = parent_node_CP2.split('+')[0]
                    if parent_node_CP3.find('+') > -1:
                        # If it is, arbitrarily choose the first sequence
                        parent_node_CP3 = parent_node_CP3.split('+')[0]

                    parent_node_sequence = [parent_node_CP1[i] + parent_node_CP2[i] + parent_node_CP3[i] for i in
                                            range(len(parent_node_CP1))]
                    parent_node_sequence = ''.join(parent_node_sequence)

                    # Compute mutability of parent node
                    parent_node_mutability_WS = seq_mutability(parent_node_sequence)
                    parent_node_mutability_aggregated = aggregated_mutability(parent_node_sequence, partition_points)

                    # --------------------------------- Annotate mutability on node ----------------------------------------

                    # Add whole-sequence mutability:
                    for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                        node_mutability = node_mutability_WS[1][metric]
                        parent_mutability = parent_node_mutability_WS[1][metric]
                        contrast = node_mutability - parent_mutability

                        if contrast == 0:
                            contrast_direction = ' '
                        else:
                            contrast_direction = '+' if contrast > 0 else '-'

                        node.annotations.add_new(metric + '_WS', str(node_mutability))
                        node.annotations.add_new(metric + '_WS_contrast', str(contrast))
                        node.annotations.add_new(metric + '_WS_contrast_direction', contrast_direction)


                        # If parent node is root node, annotate parent node:
                        if parent_node is tree.nodes()[0]:
                            parent_node.annotations.add_new(metric + '_WS', str(parent_mutability))

                    # Add FR and CDR aggregated mutability:
                    for region in ['FR', 'CDR']:
                        for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                            node_mutability = node_mutability_aggregated[region + '_mutability'][metric]
                            parent_mutability = parent_node_mutability_aggregated[region + '_mutability'][metric]
                            contrast = node_mutability - parent_mutability

                            if contrast == 0:
                                contrast_direction = ' '
                            else:
                                contrast_direction = '+' if contrast > 0 else '-'

                            node.annotations.add_new(metric + '_' + region, str(node_mutability))
                            node.annotations.add_new(metric + '_' + region + '_contrast', str(contrast))
                            node.annotations.add_new(metric + '_' + region + '_contrast_direction', contrast_direction)


                            # If parent node is root node, annotate parent node:
                            if parent_node is tree.nodes()[0]:
                                parent_node.annotations.add_new(metric + '_' + region, str(parent_mutability))


                    # ---------------------------------------- Rate annotation ---------------------------------------------
                    # Rate from molecular clock is already annotated, find robust counting rate

                    # Synonymous and Non-synonymous counts from robust counting:
                    S_count = float(node.annotations.get_value('S'))

                    N_count = float(node.annotations.get_value('N'))

                    # Get branch length in time
                    time_length = node.edge_length

                    # Convert counts into rates
                    S_rate = S_count / time_length
                    N_rate = N_count / time_length

                    # Normalize rates by number of sites:
                    S_rate = S_rate / n_sites
                    N_rate = N_rate / n_sites

                    # From syn. and non-syn. rates, get total rate according to robust counting.
                    total_rate_RC = S_rate + N_rate

                    # Annotate rates on node:
                    node.annotations.add_new('total_rate_RC', str(total_rate_RC))
                    node.annotations.add_new('S_rate', str(S_rate))
                    node.annotations.add_new('N_rate', str(N_rate))

        annotated_tree_string = tree.as_string('nexus')

        annotated_tree_string = annotated_tree_string.split('TREE 1 = ')[1]

        annotated_tree_string = 'tree TREE1 = [&R] ' + annotated_tree_string

        if annotated_tree_string.find('\nEND;\n\n') is -1:
            annotated_tree_string += '\nEND;\n\n'

        output_file.write(annotated_tree_string)

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)


