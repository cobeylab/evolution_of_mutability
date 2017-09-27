 #!/usr/bin/python

"""Takes each tree from a NEXUS file (output from BEAST). For each node on each tree, finds sequence mutability at that node based on the inferred ancestral sequence (or observed sequence, for the tips). Computes pair-wise linear regression statistics between mutability (as a response variable) and the distance and time between each the node and the root of the tree (predictors).
 This is a "by-node" analysis, as opposed to the "by-branch" analysis in "mutability_correlations"
"""

import sys
import re
from dendropy import Tree
from scipy.stats import linregress

# Import mutability functions and partition points from mutability folder
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
        trees_file_path = chain_directory + chain_id + '_tree_sample.trees'

        # Paths to the 2 output files: one with regression results and the other with the actual points

        correlations_file_path = '../../results/mutability_vs_time/observed_lineages/' + clone + '_' + prior + '/'
        correlations_file_path += chain_id + '_mutability_vs_time_correlations.csv'

        points_file_path = '../../results/mutability_vs_time/observed_lineages/' + clone + '_' + prior + '/'
        points_file_path += chain_id + '_mutability_vs_time_points.csv'

        xml_file_path = '../../analyses/BEAST/observed_lineages/' + clone + '_' + prior + '/'
        xml_file_path += chain_id[0:len(chain_id) - 1] + '.xml'

    # If chain is from a simulated scenario, use arbitrary partition points (based on VRC01_13)
    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        trees_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        trees_file_path += chain_id + '_tree_sample.trees'

        correlations_file_path = '../../results/mutability_vs_time/simulated_alignments/' + chain_id + '/'
        correlations_file_path += chain_id + '_mutability_vs_time_correlations.csv'

        points_file_path = '../../results/mutability_correlations/simulated_alignments/' + chain_id + '/'
        points_file_path += chain_id + '_mutability_vs_time_points.csv'

        xml_file_path = '../../analyses/BEAST/simulated_alignments/' + chain_id + '/'+ chain_id + '.xml'

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

    # ========================================= OPEN TREES FILE AND OUTPUT FILE ========================================

    with open(trees_file_path, 'r') as trees_file, open(correlations_file_path, 'w') as correlations_file, open(points_file_path, 'w') as points_file:
        # Dictionary linking BEAST numbering (from NEXUS file) to taxon labels:
        number_to_id = {}

        # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
        n_semicolons = 0
        while n_semicolons < 6:
            line = trees_file.next()

            translation_search =  re.search(r'[0-9]+ [A-Z|0-9|_]+',line)
            #If lines contain a taxon-id-to-number translation...
            if translation_search is not None:
                translation = translation_search.group()
                number = re.search(r'[0-9]*',translation).group()
                taxon_id = translation.replace(number + ' ', '')
                number_to_id[number] = taxon_id
                #print line

            # Make sure no trees are being skipped
            assert line.find('STATE') is -1, 'Header-skipping loop skipped a tree.'
            if line.find(';') > -1:
                n_semicolons += 1

        line_has_tree = True

        # List with results for each tree. Each element is a dictionary ("output", in the loop below) with stats results.
        tree_results_list = []

        # =============== WRITE HEADER OF POINTS CSV FILE (CORRELATIONS CSV FILE IS EXPORTED AT THE END) ===============
        points_file_header = 'tree,node,node_is_tip'

        for region in ['WS', 'FR', 'CDR']:
            for metric in ['S5F', '7M', 'HS', 'CS', 'OHS']:
                points_file_header += ',' + metric + '_' + region

        points_file_header += ',node_distance,node_distance_RC,node_time'
        points_file_header += '\n'
        points_file.write(points_file_header)

        # ============================================ ANALYZING TREES =================================================


        # For each tree, compute relationships
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

                # Each element of the lists below corresponds to a node

                # List with time at the node (since the MRCA in the tree)
                # Not to be confused with the branch length in time
                list_times = []

                # List with distance between each node and the root of the tree. For distances per molecular clock rates...
                list_distances = []
                # ...and per robust counting...
                list_distances_RC = []

                # List of node types (internal or external)
                list_node_type = []

                # Lists with node mutability in whole sequence...
                list_S5F_WS, list_7M_WS, list_HS_WS, list_CS_WS, list_OHS_WS = [], [], [], [], []

                # ...in frameworks...
                list_S5F_FR, list_7M_FR, list_HS_FR, list_CS_FR, list_OHS_FR = [], [], [], [], []

                # ...in CDRs...
                list_S5F_CDR, list_7M_CDR, list_HS_CDR, list_CS_CDR, list_OHS_CDR = [], [], [], [], []

                # Get values for each node
                for node in tree.nodes():

                    node_number = node.label
                    

                    if node.is_leaf() and number_to_id[str(node.taxon).replace("'", '')] == 'KJ134124_119':
                        print 'Skipping VRC26 sequence with missing CDR2 (KJ134124_119)'
                    else:

                        # Find time since root for the node
                        # (note tree is scaled by time so the command below gets time, not distance)
                        node_time = node.distance_from_root()

                        # To find node's distance from root, first find all branch time lengths and rates connecting it to the root
                        # (robust-counting per-site distances can be found directly from 'S' and 'N' without reference to time lengths)
                        node_to_root_time_lengths = [0.0]
                        node_to_root_rates = [0.0]
                        node_to_root_lengths_RC = [0.0]

                        # If node is not already the root node, move down the tree to the root node:
                        if node is not tree.nodes()[0]:
                            focal_node = node
                            while focal_node is not tree.nodes()[0]:
                                node_to_root_time_lengths.append(focal_node.edge_length)
                                node_to_root_rates.append(float(focal_node.annotations.get_value('rate')))
                                RC_length = float(focal_node.annotations.get_value('S')) + float(focal_node.annotations.get_value('N'))
                                RC_length = float(RC_length) / n_sites
                                node_to_root_lengths_RC.append(RC_length)

                                focal_node = focal_node.parent_node

                        # Check that summing times down to the root matches the .distance_from_root() dendropy method
                        assert sum(node_to_root_time_lengths) == node_time, "Distance/time to root computations are inconsistent with dendropy's .distance_from_root()"

                        # Find node's distance to the root by multiplying molecular clock rates and time lengths, then summing
                        node_to_root_distance = sum([node_to_root_time_lengths[i] * node_to_root_rates[i] for i in range(len(node_to_root_rates))])

                        # Find nodes's robust counting distance by summing robust counting branch lengths
                        node_to_root_distance_RC = sum(node_to_root_lengths_RC)

                        # Find node's sequence, if it is a tip:
                        if node.is_leaf():
                            # Find node number from BEAST tree
                            taxon_number = str(node.taxon).replace("'",'')
                            # Find node taxon id (the original name for the tip sequence):
                            node_id = number_to_id[taxon_number]

                            # Get mutability from tip mutability dictionaries
                            node_mutability_WS = mutability_WS_tips[node_id]
                            node_mutability_aggregated = mutability_aggregated_tips[node_id]

                            list_node_type.append('terminal')
                        else:
                            # If node is not a tip, find its sequence from the tree annotation:

                            node_CP1 = node.annotations.get_value(alignment_id + '.CP1')
                            node_CP2 = node.annotations.get_value(alignment_id + '.CP2')
                            node_CP3 = node.annotations.get_value(alignment_id + '.CP3')

                            node_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
                            node_sequence = ''.join(node_sequence)

                            node_mutability_WS = seq_mutability(node_sequence)
                            node_mutability_aggregated = aggregated_mutability(node_sequence, partition_points)

                            list_node_type.append('internal')

                        # Append values for this node to corresponding lists:
                        list_times.append(node_time)
                        list_distances.append(node_to_root_distance)
                        list_distances_RC.append(node_to_root_distance_RC)

                        list_S5F_WS.append(node_mutability_WS[1]['mean_S5F'])
                        list_7M_WS.append(node_mutability_WS[1]['mean_7M'])
                        list_HS_WS.append(node_mutability_WS[1]['HS'])
                        list_CS_WS.append(node_mutability_WS[1]['CS'])
                        list_OHS_WS.append(node_mutability_WS[1]['OHS'])

                        list_S5F_FR.append(node_mutability_aggregated['FR_mutability']['mean_S5F'])
                        list_7M_FR.append(node_mutability_aggregated['FR_mutability']['mean_7M'])
                        list_HS_FR.append(node_mutability_aggregated['FR_mutability']['HS'])
                        list_CS_FR.append(node_mutability_aggregated['FR_mutability']['CS'])
                        list_OHS_FR.append(node_mutability_aggregated['FR_mutability']['OHS'])

                        list_S5F_CDR.append(node_mutability_aggregated['CDR_mutability']['mean_S5F'])
                        list_7M_CDR.append(node_mutability_aggregated['CDR_mutability']['mean_7M'])
                        list_HS_CDR.append(node_mutability_aggregated['CDR_mutability']['HS'])
                        list_CS_CDR.append(node_mutability_aggregated['CDR_mutability']['CS'])
                        list_OHS_CDR.append(node_mutability_aggregated['CDR_mutability']['OHS'])

                        # ------------------------ Write new line of the points output file ----------------------------

                        # Concatenate results in the order specified by points_file_header

                        points_line = str(state) + ',' +node_number + ',' + str(node.is_leaf()*1)

                        # Add whole-sequence mutability
                        for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                            points_line += ',' + str(node_mutability_WS[1][metric])

                        # Add FR and CDR mutability
                        for region in ['FR', 'CDR']:
                            for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                                points_line += ',' + str(node_mutability_aggregated[region + '_mutability'][metric])

                        # Add distance, distance_RC and time
                        points_line += ',' + ','.join([str(node_to_root_distance), str(node_to_root_distance_RC),
                                                               str(node_time)])

                        points_line += '\n'

                        points_file.write(points_line)


                # Dictionary with all lists (response variables and predictors):
                lists = {"S5F_WS": list_S5F_WS, "7M_WS": list_7M_WS, "HS_WS": list_HS_WS, "CS_WS": list_CS_WS, "OHS_WS": list_OHS_WS,
                         "S5F_FR": list_S5F_FR, "7M_FR": list_7M_FR, "HS_FR": list_HS_FR, "CS_FR": list_CS_FR, "OHS_FR": list_OHS_FR,
                         "S5F_CDR": list_S5F_CDR, "7M_CDR": list_7M_CDR, "HS_CDR": list_HS_CDR, "CS_CDR": list_CS_CDR, "OHS_CDR": list_OHS_CDR,
                         "node_time": list_times, "node_distance": list_distances, "node_distance_RC": list_distances_RC,
                         "node_type": list_node_type
                        }

                # Output dictionary with statistical results for this tree (values for one line of the final output file)
                output = {}

                # Add maximum and minimum values times, distances and mutabilities to the dictionary
                for key in lists.keys():
                    output['max_' + key] = max(lists[key])
                    output['min_' + key] = min(lists[key])


                # Compute relationship between mutability and time (slope, intercept, r [correlation coefficient])
                # First argument of linregress is x, second is y. Slope, intercept and r are elements 0, 1 and 2 from output list
                for metric in ['S5F', '7M', 'HS', 'CS', 'OHS']:
                    for region in ['WS','FR','CDR']:

                        # Relationship for all nodes
                        predictor = lists['node_time']
                        predicted = lists[metric + '_' + region]

                        linear_regression = linregress(x = predictor, y = predicted)

                        output['slope_' + metric + '_' + region + '_vs_node_time'] = linear_regression[0]
                        output['intercept_' + metric + '_' + region + '_vs_node_time'] = linear_regression[1]
                        output['r_' + metric + '_' + region + '_vs_node_time'] = linear_regression[2]

                        # Relationship for observed values only
                        predictor_obs_only = [lists['node_time'][i] for i in range(len(lists['node_type'])) if lists['node_type'][i] == 'terminal']
                        predicted_obs_only = [lists[metric + '_' + region][i] for i in range(len(lists['node_type'])) if lists['node_type'][i] == 'terminal']

                        linear_regression_obs_only = linregress(x = predictor_obs_only, y = predicted_obs_only)

                        output['slope_' + metric + '_' + region + '_vs_node_time_obs_only'] = linear_regression_obs_only[0]
                        output['intercept_' + metric + '_' + region + '_vs_node_time_obs_only'] = linear_regression_obs_only[1]
                        output['r_' + metric + '_' + region + '_vs_node_time_obs_only'] = linear_regression_obs_only[2]

                # Append results for this tree
                tree_results_list.append(output)

            # If there's not a tree in this line:
            else:
                line_has_tree = False

        # ======================================== EXPORTING RESULTS TO CSV FILE =======================================

        # Write results for all trees to the output csv file:
        output_header = 'X'
        # Get column names from results for 1st tree
        for key in tree_results_list[0].keys():
            output_header += ',' + key
        output_header += '\n'
        output_header = output_header.replace('X,', '')
        correlations_file.write(output_header)

        # Write all lines
        for tree_results in tree_results_list:
            line = 'X'
            for key in tree_results.keys():
                line += ',' + str(tree_results[key])
            line += '\n'
            line = line.replace('X,', '')
            correlations_file.write(line)

if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
