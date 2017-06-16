#!/usr/bin/python

"""Takes each tree from a NEXUS file (output from BEAST). For each branch on each tree, finds the substitution rates along that branch, the time and distance between the branch's parent node and the root of the tree, and the mutability of the parent node.
    For each tree, computes pair-wise linear regression statistics for the relationship between substitution rates (response variable) and distance, time and mutability (predictors).
    Further analyzes the relationship between molecular clock rates/distances and robust counting rates/distances.
    This is a "by-branch" analysis, as opposed to the "by-node" analysis in "mutability_correlations"
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

    # ========================= RETRIEVE DIRECTORIES AND FR/CDR PARTITION POINTS FROM CHAIN ID =========================
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

        correlations_file_path = '../../results/rates_vs_time/observed_lineages/' + clone + '_' + prior + '/'
        correlations_file_path += chain_id + '_rates_vs_time_correlations.csv'

        points_file_path = '../../results/rates_vs_time/observed_lineages/' + clone + '_' + prior + '/'
        points_file_path += chain_id + '_rates_vs_time_points.csv'

        # If chain is from a simulated scenario, use arbitrary partition points (based on VRC01_13)
    else:
        scenario = re.search(r'scenario[1-9][abc]*', chain_id).group()
        trees_file_path = '../../results/BEAST/simulated_alignments/' + chain_id + '/'
        trees_file_path += chain_id + '_tree_sample.trees'

        correlations_file_path = '../../results/rates_vs_time/simulated_alignments/' + chain_id + '/'
        correlations_file_path += chain_id + '_rates_vs_time_correlations.csv'
        
        points_file_path = '../../results/rates_vs_time/simulated_alignments/' + chain_id + '/'
        points_file_path += chain_id + '_rates_vs_time_points.csv'

        partition_points = [1, 76, 100, 151, 175, 289, 351]

    # ========================================= OPEN TREES FILE AND OUTPUT FILE ========================================

    with open(trees_file_path, 'r') as trees_file, open(correlations_file_path, 'w') as correlations_file, open(points_file_path, 'w') as points_file:
        # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
        n_semicolons = 0
        while n_semicolons < 6:
            line = trees_file.next()
            # Make sure no trees are being skipped
            assert line.find('STATE') is -1, 'Header-skipping loop skipped a tree.'
            if line.find(';') > -1:
                n_semicolons += 1

        line_has_tree = True

        # List with results for each tree. Each element is a dictionary ("output", in the loop below) with stats results.
        tree_results_list = []

        # =============== WRITE HEADER OF POINTS CSV FILE (CORRELATIONS CSV FILE IS EXPORTED AT THE END) ===============
        points_file_header = 'tree,parent,child,branch_is_terminal'
        for rate in ['total_rate', 'total_rate_RC','S_rate', 'N_rate']:
            points_file_header += ',' + rate + '_to_child'

        for region in ['WS', 'FR', 'CDR']:
            for metric in ['S5F', '7M', 'HS', 'CS', 'OHS']:
                points_file_header += ',' + metric + '_' + region + '_parent'

        points_file_header += ',parent_distance,parent_distance_RC,parent_time'
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

                # Remove tree header from tree_string before passing reading it with dendropy
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
                n_sites = 3*len(tree.nodes()[0].annotations[0].value)

                # Alignment identifier in codon position annotations (e.g. CH103_final_alignment.CP1)
                alignment_id = [annotation.name for annotation in tree.seed_node.annotations if annotation.name.find('CP1') > 1]
                alignment_id = alignment_id[0].replace('.CP1','')

                # Each element of the lists below corresponds to a branch
                # Times/distances to root and mutability values are for the branch's parent node

                # List with total rate values from robust counting (i.e. sum of syn. and non.syn. rates)
                list_total_rate_RC = []

                # Total rate from the molecular clock alone (i.e. independent of robust counting):
                list_total_rate = []

                list_S_rate = []
                list_N_rate = []

                # List with time at the beginning of branch (i.e. parent node's distance in time from root)
                # Not to be confused with the branch length in time =
                list_times = []

                # List with distance between the parent node of each branch and the root of the tree. For distances per molecular clock rates...
                list_distances = []
                # ...and per robust counting...
                list_distances_RC = []

                # Lists with parent node mutability in whole sequence...
                list_S5F_WS, list_7M_WS, list_HS_WS, list_CS_WS, list_OHS_WS = [], [], [], [], []

                #...in frameworks...
                list_S5F_FR, list_7M_FR, list_HS_FR, list_CS_FR, list_OHS_FR = [], [], [], [], []

                #...in CDRs...
                list_S5F_CDR, list_7M_CDR, list_HS_CDR, list_CS_CDR, list_OHS_CDR = [], [], [], [], []

                # For each node (that has a parent) (branch info is stored in daughter node, not parent)
                for node in tree.nodes():
                    node_number = node.label

                    if node.parent_node is not None:

                        # Get total rate from molecular clock:
                        total_rate = float(node.annotations.get_value('rate'))

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

                        # Find parent node
                        parent_node = node.parent_node
                        parent_node_number = parent_node.label

                        # Find time since root for the parent node:
                        parent_time = parent_node.distance_from_root()

                        # To find parent node's distance from root, first find all branch time lengths and rates connecting it to the root
                        # (robust-counting per-site distances can be found directly from 'S' and 'N' without reference to time lengths)
                        parent_to_root_time_lengths = [0]
                        parent_to_root_rates = [0]
                        parent_to_root_lengths_RC = [0]
                        focal_node = parent_node

                        # Move down the tree to the root node
                        while focal_node is not tree.nodes()[0]:
                            parent_to_root_time_lengths.append(focal_node.edge_length)
                            parent_to_root_rates.append(float(focal_node.annotations.get_value('rate')))
                            RC_length = float(focal_node.annotations.get_value('S')) + float(focal_node.annotations.get_value('N'))
                            RC_length = float(RC_length) / n_sites
                            parent_to_root_lengths_RC.append(RC_length)

                            focal_node = focal_node.parent_node

                        # Find parent's distance to the root by multiplying molecular clock rates and time lengths, then summing
                        parent_to_root_distance = sum([parent_to_root_time_lengths[i] * parent_to_root_rates[i] for i in range(len(parent_to_root_rates))])

                        # Find parent's robust counting distance by summing robust counting branch lengths
                        parent_to_root_distance_RC = sum(parent_to_root_lengths_RC)

                        # To find parent node mutability, assemble sequence from annotations of the 3 codon positions:
                        parent_CP1 = parent_node.annotations.get_value(alignment_id + '.CP1')
                        parent_CP2 = parent_node.annotations.get_value(alignment_id + '.CP2')
                        parent_CP3 = parent_node.annotations.get_value(alignment_id + '.CP3')

                        parent_sequence = [parent_CP1[i] + parent_CP2[i] + parent_CP3[i] for i in range(len(parent_CP1))]
                        parent_sequence = ''.join(parent_sequence)

                        parent_mutability_WS = seq_mutability(parent_sequence)
                        parent_mutability_aggregated = aggregated_mutability(parent_sequence, partition_points)

                        # Append values to corresponding lists:
                        list_total_rate.append(total_rate)
                        list_total_rate_RC.append(total_rate_RC)
                        list_S_rate.append(S_rate)
                        list_N_rate.append(N_rate)
                        list_times.append(parent_time)
                        list_distances.append(parent_to_root_distance)
                        list_distances_RC.append(parent_to_root_distance_RC)

                        list_S5F_WS.append(parent_mutability_WS[1]['mean_S5F'])
                        list_7M_WS.append(parent_mutability_WS[1]['mean_7M'])
                        list_HS_WS.append(parent_mutability_WS[1]['HS'])
                        list_CS_WS.append(parent_mutability_WS[1]['CS'])
                        list_OHS_WS.append(parent_mutability_WS[1]['OHS'])

                        list_S5F_FR.append(parent_mutability_aggregated['FR_mutability']['mean_S5F'])
                        list_7M_FR.append(parent_mutability_aggregated['FR_mutability']['mean_7M'])
                        list_HS_FR.append(parent_mutability_aggregated['FR_mutability']['HS'])
                        list_CS_FR.append(parent_mutability_aggregated['FR_mutability']['CS'])
                        list_OHS_FR.append(parent_mutability_aggregated['FR_mutability']['OHS'])

                        list_S5F_CDR.append(parent_mutability_aggregated['CDR_mutability']['mean_S5F'])
                        list_7M_CDR.append(parent_mutability_aggregated['CDR_mutability']['mean_7M'])
                        list_HS_CDR.append(parent_mutability_aggregated['CDR_mutability']['HS'])
                        list_CS_CDR.append(parent_mutability_aggregated['CDR_mutability']['CS'])
                        list_OHS_CDR.append(parent_mutability_aggregated['CDR_mutability']['OHS'])

                        # ------------------------ Write new line of the points output file ----------------------------

                        # Concatenate results in the order specified by points_file_header

                        points_line = str(state) + ',' + parent_node_number + ',' + node_number + ','

                        points_line += str(node.is_leaf()*1) + ','

                        points_line += ','.join([str(total_rate), str(total_rate_RC), str(S_rate), str(N_rate)])

                        # Add parent whole-sequence mutability
                        for metric in ['mean_S5F','mean_7M','HS','CS','OHS']:
                            points_line += ',' + str(parent_mutability_WS[1][metric])

                        # Add parent FR and CDR mutability
                        for region in ['FR','CDR']:
                            for metric in ['mean_S5F', 'mean_7M', 'HS', 'CS', 'OHS']:
                                points_line += ',' + str(parent_mutability_aggregated[region + '_mutability'][metric])

                        # Add parent distance, distance_RC and time
                        points_line += ',' + ','.join([str(parent_to_root_distance), str(parent_to_root_distance_RC),
                                                           str(parent_time)])

                        points_line += '\n'

                        points_file.write(points_line)


                # Dictionary with all lists (rates, mutability, time and distance from root)
                lists = {"total_rate": list_total_rate, "total_rate_RC": list_total_rate_RC, "S_rate": list_S_rate,
                         "N_rate": list_N_rate, "parent_time": list_times, "parent_distance": list_distances, "parent_distance_RC": list_distances_RC,
                         "S5F_WS": list_S5F_WS, "7M_WS": list_7M_WS, "HS_WS": list_HS_WS, "CS_WS": list_CS_WS, "OHS_WS": list_OHS_WS,
                         "S5F_FR": list_S5F_FR, "7M_FR": list_7M_FR, "HS_FR": list_HS_FR, "CS_FR": list_CS_FR, "OHS_FR": list_OHS_FR,
                         "S5F_CDR": list_S5F_CDR, "7M_CDR": list_7M_CDR, "HS_CDR": list_HS_CDR, "CS_CDR": list_CS_CDR, "OHS_CDR": list_OHS_CDR,
                         }

                # Output dictionary with statistical results for this tree (values for one line of the final output file)
                output = {}
                
                # Add maximum and minimum values of rates, times, distances and mutabilities to the dictionary
                for key in lists.keys():
                    output['max_' + key] = max(lists[key])
                    output['min_' + key] = min(lists[key])

                # Compute relationship between each rate and time (slope, intercept, r [correlation coefficient])
                # First argument of linregress is x, second is y. Slope, intercept and r are elements 0, 1 and 2 from output list
                for rate in ['total_rate', 'total_rate_RC', 'S_rate', 'N_rate']:
                    output['slope_' + rate + '_vs_parent_time'] = linregress(lists['parent_time'], lists[rate])[0]
                    output['intercept_' + rate + '_vs_parent_time'] = linregress(lists['parent_time'], lists[rate])[1]
                    output['r_' + rate + '_vs_parent_time'] = linregress(lists['parent_time'], lists[rate])[2]

                # Compute relationship between each rate and parent distance from the root
                for rate in ['total_rate', 'total_rate_RC', 'S_rate', 'N_rate']:
                    # Mol clock distance
                    output['slope_' + rate + '_vs_parent_distance'] = linregress(lists['parent_distance'], lists[rate])[0]
                    output['intercept_' + rate + '_vs_parent_distance'] = linregress(lists['parent_distance'], lists[rate])[1]
                    output['r_' + rate + '_vs_parent_distance'] = linregress(lists['parent_distance'], lists[rate])[2]

                    # RC distance
                    output['slope_' + rate + '_vs_parent_distance_RC'] = linregress(lists['parent_distance_RC'], lists[rate])[0]
                    output['intercept_' + rate + '_vs_parent_distance_RC'] = linregress(lists['parent_distance_RC'], lists[rate])[1]
                    output['r_' + rate + '_vs_parent_distance_RC'] = linregress(lists['parent_distance_RC'], lists[rate])[2]

                # ------------------- Compute relationship between mol. clock rate and robust counting rate ------------------------
                output['slope_total_rate_RC_vs_total_rate'] = linregress(lists['total_rate'], lists['total_rate_RC'])[0]
                output['intercept_total_rate_RC_vs_total_rate'] = linregress(lists['total_rate'], lists['total_rate_RC'])[1]
                output['r_total_rate_RC_vs_total_rate'] = linregress(lists['total_rate'], lists['total_rate_RC'])[2]

                # ------------ Compute relationship btwn mol. clock distance from root and robust counting dist from root-----------
                output['slope_parent_distance_to_root_RC_vs_parent_distance_to_root'] = linregress(lists['parent_distance'],
                                                                                                   lists['parent_distance_RC'])[0]
                output['intercept_parent_distance_to_root_RC_vs_parent_distance_to_root'] = linregress(lists['parent_distance'],
                                                                                                   lists['parent_distance_RC'])[1]
                output['r_parent_distance_to_root_RC_vs_parent_distance_to_root'] = linregress(lists['parent_distance'],
                                                                                               lists['parent_distance_RC'])[2]

                # ------------------------------- Compute relationships between rates and mutability -------------------------------

                for rate in ['total_rate', 'total_rate_RC', 'S_rate', 'N_rate']:
                    for metric in ['S5F', '7M', 'HS', 'CS', 'OHS']:
                        for region in ['WS','FR','CDR']:
                            output['slope_'+rate+'_vs_'+metric+'_'+region] = linregress(lists[metric+'_'+region], lists[rate])[0]
                            output['intercept_'+rate+'_vs_'+metric+'_'+region] = linregress(lists[metric+'_'+region], lists[rate])[1]
                            output['r_'+rate+'_vs_'+metric+'_'+region] = linregress(lists[metric+'_'+region], lists[rate])[2]

                tree_results_list.append(output)

            #If there's not a tree in this line:
            else:
                line_has_tree = False

        # ======================================== EXPORTING RESULTS TO CSV FILE =======================================

        # Write results for all trees to the output csv file:
        output_header = 'X'
        # Get column names from results for 1st tree
        for key in tree_results_list[0].keys():
            output_header += ',' + key
        output_header += '\n'
        output_header = output_header.replace('X,','')
        correlations_file.write(output_header)

        # Write all lines
        for tree_results in tree_results_list:
            line = 'X'
            for key in tree_results.keys():
                line += ',' + str(tree_results[key])
            line += '\n'
            line = line.replace('X,','')
            correlations_file.write(line)

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
