#!/usr/bin/python

"""Takes each node in an MCC tree (include obs. sequences), computes its mutability and compares it to the mutability expected by randomizing the sequence while keeping the amino acid sequence constant.
Employs two randomization procedures: one randomizes all sites ('allsites'), the other randomizes only sites that differ between the sequence and the UCA ('diffsites')
"""

import sys
import re
from dendropy import Tree
from numpy import random
from copy import deepcopy

# Import mutability functions and partition points from analyses/mutability folder
sys.path.insert(0, '../mutability/')

from mutability_function import seq_mutability, aggregated_mutability
from partition_points import partition_points_dic

from gen_code_DNA import genetic_code

#Dictionary of codon frequencies per amino acid:
codon_freqs = {}

# Read table of human codon frequencies (see reference in CSV file):
with open('human_codon_frequencies.csv','r') as codon_freq_file:
    for line in codon_freq_file:
        # Skip first 3 lines
        if line.find('#') == -1:
            aa = line.split(' ')[1]
            codon = line.split(' ')[0]
            frequency = line.split(' ')[2]

            if aa not in codon_freqs.keys():
                codon_freqs[aa] = {codon:frequency}
            else:
                codon_freqs[aa][codon] = frequency

    # Normalize so that freqs for each aa sum to 1 (some do not probably because of rounding error in the authors' report)
    for aa in codon_freqs.keys():
        sum = 0
        for codon in codon_freqs[aa].keys():
            sum = sum + float(codon_freqs[aa][codon])

        for codon in codon_freqs[aa].keys():
            codon_freqs[aa][codon] = float(codon_freqs[aa][codon])/sum


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

        output_directory = '../../results/relative_mutability/observed_lineages/' + clone + '_' + prior + '/'
        observed_mutability_file_path = output_directory + chain_id + '_observed_mutability_MCC.csv'
        randomized_mutability_file_path = output_directory + chain_id + '_randomized_mutability_MCC.csv'

    # If chain is from a simulated scenario
    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        MCC_tree_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        MCC_tree_file_path += chain_id + '_MCC_tree.tree'

        xml_file_path = '../../analyses/BEAST/simulated_alignments/' + chain_id + '/' + chain_id + '.xml'

        output_directory = '../../results/relative_mutability/simulated_alignments/' + chain_id + '/'
        observed_mutability_file_path = output_directory + chain_id + '_observed_mutability_MCC.csv'
        randomized_mutability_file_path = output_directory + chain_id + '_randomized_mutability_MCC.csv'

        partition_points = partition_points_dic[chain_id]


    # ========================== READ TIP SEQUENCES FROM XML FILE AND COMPUTE THEIR MUTABILITY==========================
    # ------------------------------- (Tip sequences are not annotated on the BEAST trees) -----------------------------

    # Dictionaries with mutability for each tip sequence:
    mutability_WS_tips = {}
    mutability_aggregated_tips = {}

    # Dictionary recording nt sequence of each observed sequence:
    observed_sequences = {}

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
            observed_sequences[taxon_id] = sequence

            if partition_points is not None:
                mutability_aggregated_tips[taxon_id] = aggregated_mutability(sequence, partition_points)
        else:
            print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'


    # Get length of sequences
    seq_length = len(observed_sequences[observed_sequences.keys()[0]])

    # ==================================== OPEN TREE FILE AND OUTPUT FILES =============================================
    with open(MCC_tree_file_path, 'r') as tree_file, open(observed_mutability_file_path, 'w') as output_file_obs, open(randomized_mutability_file_path, 'w') as output_file_random:

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
        variables = [alignment_id + '.CP1', alignment_id + '.CP2', alignment_id + '.CP3']

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

        # ====================================== WRITE OUTPUT FILES HEADERS ============================================

        output_file_obs.write('sequence_id,time_from_root,')

        output_file_obs.write('observed_S5F_WS,observed_7M_WS,observed_HS_WS,observed_CS_WS,observed_OHS_WS,')
        output_file_obs.write('observed_S5F_FR,observed_7M_FR,observed_HS_FR,observed_CS_FR,observed_OHS_FR,')
        output_file_obs.write('observed_S5F_CDR,observed_7M_CDR,observed_HS_CDR,observed_CS_CDR,observed_OHS_CDR\n')
        
        output_file_random.write('sequence_id,time_from_root,')

        output_file_random.write('randomized_S5F_WS_allsites,randomized_7M_WS_allsites,randomized_HS_WS_allsites,randomized_CS_WS_allsites,randomized_OHS_WS_allsites,')
        output_file_random.write('randomized_S5F_FR_allsites,randomized_7M_FR_allsites,randomized_HS_FR_allsites,randomized_CS_FR_allsites,randomized_OHS_FR_allsites,')
        output_file_random.write('randomized_S5F_CDR_allsites,randomized_7M_CDR_allsites,randomized_HS_CDR_allsites,randomized_CS_CDR_allsites,randomized_OHS_CDR_allsites,')
        
        output_file_random.write('randomized_S5F_WS_diffsites,randomized_7M_WS_diffsites,randomized_HS_WS_diffsites,randomized_CS_WS_diffsites,randomized_OHS_WS_diffsites,')
        output_file_random.write('randomized_S5F_FR_diffsites,randomized_7M_FR_diffsites,randomized_HS_FR_diffsites,randomized_CS_FR_diffsites,randomized_OHS_FR_diffsites,')
        output_file_random.write('randomized_S5F_CDR_diffsites,randomized_7M_CDR_diffsites,randomized_HS_CDR_diffsites,randomized_CS_CDR_diffsites,randomized_OHS_CDR_diffsites\n')


        # ================================== DO THE ANALYSIS FOR EACH NODE =============================================
        for node in tree.nodes():

            if node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'KJ134124_119':
                print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'
            else:
                print 'Node number ' + node.label
                # Get sequence time from root
                seq_time_from_root = node.distance_from_root()

                # IF NODE IS AN OBSERVED SEQUENCE
                if node.is_leaf():
                    # Find taxon number from BEAST tree
                    taxon_number = str(node.taxon).replace("'", '')
                    # Find node taxon id (the original name for the tip sequence):
                    node_id = number_to_id[taxon_number]

                    # Get original sequence (extracted from XML previously)
                    original_sequence = observed_sequences[node_id]

                    # Get mutability from tip mutability dictionaries
                    observed_mutability_WS = mutability_WS_tips[node_id]
                    observed_mutability_aggregated = mutability_aggregated_tips[node_id]

                # If node is an internal node:
                else:
                    node_id = 'Node_' + node.label
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

                    original_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
                    original_sequence = ''.join(original_sequence)

                    # Compute mutability
                    observed_mutability_WS = seq_mutability(original_sequence)
                    observed_mutability_aggregated = aggregated_mutability(original_sequence, partition_points)

                # Convert observed nt sequence into aa sequence
                original_aa_seq = ''
                for i in range(seq_length/3):
                    codon = original_sequence[i * 3:(i * 3 + 1) + 2]

                    # If codon contains gap or other non ACGT symbols, translate as X
                    if set(codon) <= {'A', 'G', 'C', 'T'}:
                        aa = genetic_code[codon]
                    else:
                        aa = 'X'
                    original_aa_seq = original_aa_seq + aa

                # If node is root node, record is AA sequence to be used as reference for subsequent nodes:
                if node == tree._get_seed_node():
                    ancestral_aa_seq = original_aa_seq


                # Write observed mutability results
                output_file_obs.write(node_id + ',' + str(seq_time_from_root) + ',')

                output_file_obs.write(str(observed_mutability_WS[1]['mean_S5F']) + ',')
                output_file_obs.write(str(observed_mutability_WS[1]['mean_7M']) + ',')
                output_file_obs.write(str(observed_mutability_WS[1]['HS']) + ',')
                output_file_obs.write(str(observed_mutability_WS[1]['CS']) + ',')
                output_file_obs.write(str(observed_mutability_WS[1]['OHS']) + ',')

                output_file_obs.write(str(observed_mutability_aggregated['FR_mutability']['mean_S5F']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['FR_mutability']['mean_7M']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['FR_mutability']['HS']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['FR_mutability']['CS']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['FR_mutability']['OHS']) + ',')

                output_file_obs.write(str(observed_mutability_aggregated['CDR_mutability']['mean_S5F']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['CDR_mutability']['mean_7M']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['CDR_mutability']['HS']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['CDR_mutability']['CS']) + ',')
                output_file_obs.write(str(observed_mutability_aggregated['CDR_mutability']['OHS']) + '\n')

                # Define lists of codon sites to randomize for each randomization procedure
                randomized_sites = {}

                # For the 'allsites' randomization, all sites are randomized
                randomized_sites['allsites'] = range(seq_length / 3)

                # For the 'diffsites' randomization, only codon sites that differ from UCA
                randomized_sites['diffsites'] = [j for j in range(seq_length / 3) if
                                                 original_aa_seq[j] != ancestral_aa_seq[j]]
                # For 1000 replicates:
                for i in range(1000):

                    output_file_random.write(node_id + ',' + str(seq_time_from_root) + ',')

                    # Generate empty randomized sequences:
                    randomized_sequences = {}
                    randomized_sequences['allsites'] = ''
                    randomized_sequences['diffsites'] = ''

                    # For each randomization type
                    for key in ['allsites','diffsites']:
                        sites_to_randomize = randomized_sites[key]

                        # For each codon site in the original sequence:
                        for k in range(seq_length/3):
                            original_codon = original_sequence[k * 3:(k * 3 + 1) + 2]

                            # If site is not to be randomized:
                            if k not in sites_to_randomize:
                                randomized_sequences[key] = randomized_sequences[key] + original_codon

                            # If site is to be randomized...
                            else:
                                #Observed sequences may contain gaps. If the codon is gapped, just repeat it
                                # Same if codon contains any symbol other than AGCT

                                # If codon is a good codon:
                                if set(original_codon) <= {'A','G','C','T'}:

                                    original_aa = genetic_code[original_codon]
                                    # List of all possible codons for the same AA, including the original codon
                                    possible_codons = [codon for codon in genetic_code.keys() if
                                                       genetic_code[codon] == original_aa]

                                    # Sample codons (inc. original) for original aa according to their relative frequencies
                                    codon_sampling_probs = []
                                    for codon in possible_codons:
                                        sampling_prob = codon_freqs[original_aa][codon]
                                        codon_sampling_probs.append(sampling_prob)

                                    randomized_codon = random.choice(possible_codons, p=codon_sampling_probs)
                                    #print randomized_codon

                                    assert(genetic_code[randomized_codon] == genetic_code[original_codon]), 'Randomized sequence has a different amino acid sequence'

                                else:
                                    randomized_codon = original_codon

                                randomized_sequences[key] += randomized_codon
                        #print randomized_sequence
                        #print original_sequence

                        # Calculate mutability of randomized sequence
                        randomized_mutability_WS = seq_mutability(randomized_sequences[key])
                        randomized_mutability_aggregated = aggregated_mutability(randomized_sequences[key], partition_points)

                        # Write mutability results for this randomization for the current randomization procedure

                        output_file_random.write(str(randomized_mutability_WS[1]['mean_S5F']) + ',')
                        output_file_random.write(str(randomized_mutability_WS[1]['mean_7M']) + ',')
                        output_file_random.write(str(randomized_mutability_WS[1]['HS']) + ',')
                        output_file_random.write(str(randomized_mutability_WS[1]['CS']) + ',')
                        output_file_random.write(str(randomized_mutability_WS[1]['OHS']) + ',')

                        output_file_random.write(str(randomized_mutability_aggregated['FR_mutability']['mean_S5F']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['FR_mutability']['mean_7M']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['FR_mutability']['HS']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['FR_mutability']['CS']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['FR_mutability']['OHS']) + ',')

                        output_file_random.write(str(randomized_mutability_aggregated['CDR_mutability']['mean_S5F']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['CDR_mutability']['mean_7M']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['CDR_mutability']['HS']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['CDR_mutability']['CS']) + ',')
                        output_file_random.write(str(randomized_mutability_aggregated['CDR_mutability']['OHS']))

                        if key == 'allsites':
                            output_file_random.write(',')
                        else:
                            output_file_random.write('\n')

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)