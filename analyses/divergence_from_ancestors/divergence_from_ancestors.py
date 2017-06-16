"""Takes sequences from MCC trees and computes their % A.A and NT divergence from ancestor in FRs and CDRs
"""
import sys
# Import partition_points dictionary from mutability folder
sys.path.insert(0, '../rate_correlations/')

# Import gen. code and function for identifying sequence differences from syn./non-syn analysis:
sys.path.insert(0, '../S_NS_mutability_changes/')

from mutation_functions import sequence_differences, randomize_sequence
from partition_points import partition_points_dic
from gen_code_DNA import genetic_code
from dendropy import Tree

import re
#from numpy import random, percentile, mean
from copy import deepcopy
import csv
#from itertools import permutations

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

        output_directory = '../../results/divergence_from_ancestors/'
        output_file_path = output_directory + chain_id + '_divergence_from_ancestors.csv'



    # ======================================= READ TIP SEQUENCES FROM XML FILE =========================================
    # ------------------------------- (Tip sequences are not annotated on the BEAST trees) -----------------------------

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

            obs_sequence[taxon_id] = sequence
        else:
            print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'


    # Get length of sequences
    seq_length = len(obs_sequence[obs_sequence.keys()[0]])


    # ==================================== OPEN TREE FILE AND OUTPUT FILES =============================================

    with open(MCC_tree_file_path, 'r') as tree_file, open(output_file_path, 'w') as output_file:

        output_header = 'node,time_to_root,nt_divergence_from_ancestor_FR,nt_divergence_from_ancestor_CDR,'
        output_header += 'aa_divergence_from_ancestor_FR,aa_divergence_from_ancestor_CDR\n'

        output_file.write(output_header)

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

        # =========== FIND ANCESTRAL SEQUENCE
        ancestral_node = tree.nodes()[0]

        # Find ancestral node sequence from the tree annotation:
        ancestral_node_CP1 = ancestral_node.annotations.get_value(alignment_id + '.CP1')
        ancestral_node_CP2 = ancestral_node.annotations.get_value(alignment_id + '.CP2')
        ancestral_node_CP3 = ancestral_node.annotations.get_value(alignment_id + '.CP3')

        # Arbitrarily choosing first listed reconstruction in case of ambiguity
        if ancestral_node_CP1.find('+') > -1:
            # If it is, arbitrarily choose the first sequence
            ancestral_node_CP1 = ancestral_node_CP1.split('+')[0]
            print 'Ambiguous reconstruction in ancestral_node ; arbitrarily choosing first listed reconstruction'
        if ancestral_node_CP2.find('+') > -1:
            # If it is, arbitrarily choose the first sequence
            ancestral_node_CP2 = ancestral_node_CP2.split('+')[0]
            print 'Ambiguous reconstruction in ancestral_node ; arbitrarily choosing first listed reconstruction'
        if ancestral_node_CP3.find('+') > -1:
            # If it is, arbitrarily choose the first sequence
            ancestral_node_CP3 = ancestral_node_CP3.split('+')[0]
            print 'Ambiguous reconstruction in ancestral_node ; arbitrarily choosing first listed reconstruction'

        ancestral_node_sequence = [ancestral_node_CP1[i] + ancestral_node_CP2[i] + ancestral_node_CP3[i] for i in range(len(ancestral_node_CP1))]
        ancestral_node_sequence = ''.join(ancestral_node_sequence)

        # ================================== DO THE ANALYSIS FOR EACH NODE =============================================
        for node in tree.nodes():
            # Find node number (index in tree.nodes())
            node_number = node.label

            # If node is not the root:
            if (node.parent_node is not None):
                print node_number

                # Ignoring a weird VRC26 sequence whose CDR2 is entirely gaps:
                if node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'KJ134124_119':
                    print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'

                # Ignoring germline sequence:
                elif node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'GERMLINE_00':
                    print 'Skipping GERMLINE sequence'
                else:

                    # Find node time to root:
                    node_time = node.distance_from_root()

                    # Find node sequence, if it is a tip:
                    if node.is_leaf():
                        # Find taxon number from BEAST tree
                        taxon_number = str(node.taxon).replace("'", '')
                        # Find node taxon id (the original name for the tip sequence):
                        node_id = number_to_id[taxon_number]
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


                    # Node vs. ancestor differences:
                    ancestor_differences = sequence_differences(ancestral_node_sequence, node_sequence, partition_points)

                    # New line of output file:
                    new_line = node_number + ',' + str(node_time) + ','
                    new_line += str(ancestor_differences['nt_divergence_FRs']) + ',' + str(ancestor_differences['nt_divergence_CDRs'])
                    new_line += ',' + str(ancestor_differences['aa_divergence_FRs']) + ',' + str(ancestor_differences['aa_divergence_CDRs'])
                    new_line += '\n'

                    output_file.write(new_line)

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)

