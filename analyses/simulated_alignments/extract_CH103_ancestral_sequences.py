#!/usr/bin/python

"""Imports MCC trees for all CH103 chains that have them (i.e. excluding chains that are still running). Extracts the sequence from the root node of each tree, computes its mutability, and exports the sequences and the mutability results.
"""

import sys
import os
import re
from dendropy import Tree

# Import mutability functions from analyses/mutability folder
sys.path.insert(0, '../mutability/')
# Import partition_points dictionary from rate_correlations folder
sys.path.insert(0, '../rate_correlations/')

from mutability_function import seq_mutability, aggregated_mutability
from partition_points import partition_points_dic

clone = 'CH103'

output_file_name = 'ancestral_CH103_sequences.csv'

output_header = 'chain,sequence'
for region in ['WS','FR','CDR']:
    for metric in ['S5F','7M','HS','CS','OHS']:
        output_header+= ',' + metric + '_' + region
    
def main(argv):
    with open(output_file_name, 'w') as output_file:
        os.chdir('../../results/BEAST/observed_lineages/')
        output_file.write(output_header + '\n')
    
        partition_points = partition_points_dic[clone]

        priors = ['constant','logistic','exponential']

        for prior in priors:
            if prior is 'constant':
                chains = ['a','b','c','d']
            else:
                chains = ['a','b']
        
            prior_abbrev = prior[0:3]
    
            for chain in chains:               
                chain_directory = clone + '_' + prior + '/' + clone + '_' + prior_abbrev + '_run1' + chain + '/'    
                MCC_file = clone + '_' + prior_abbrev + '_run1' + chain + '_MCC_tree.tree'
        
                if MCC_file in os.listdir(chain_directory):
                    
                    print MCC_file                       
                   
                    MCC_tree = Tree.get_from_path(chain_directory + MCC_file,'nexus')

                    ancestor = MCC_tree.nodes()[0]
        
                    # Alignment identifier in codon position annotations (e.g. CH103_final_alignment.CP1)
                    alignment_id = [annotation.name for annotation in MCC_tree.seed_node.annotations if annotation.name.find('CP1') > 1]
                    alignment_id = alignment_id[0].replace('.CP1', '').replace('.set','').replace('.prob','')
        
                    # Assemble sequence from codon positions
                    ancestor_CP1 = ancestor.annotations.get_value(alignment_id + '.CP1')
                    ancestor_CP2 = ancestor.annotations.get_value(alignment_id + '.CP2')
                    ancestor_CP3 = ancestor.annotations.get_value(alignment_id + '.CP3')

                    ancestor_sequence = [ancestor_CP1[i] + ancestor_CP2[i] + ancestor_CP3[i] for i in range(len(ancestor_CP1))]
                    ancestor_sequence = ''.join(ancestor_sequence)
    
                    # Compute mutability of ancestral sequence
                    ancestor_mutability_WS = seq_mutability(ancestor_sequence)
                    ancestor_mutability_aggregated = aggregated_mutability(ancestor_sequence, partition_points)

                    ancestor_S5F_WS = str(ancestor_mutability_WS[1]['mean_S5F'])
                    ancestor_7M_WS = str(ancestor_mutability_WS[1]['mean_7M'])
                    ancestor_HS_WS = str(ancestor_mutability_WS[1]['HS'])
                    ancestor_CS_WS = str(ancestor_mutability_WS[1]['CS'])
                    ancestor_OHS_WS = str(ancestor_mutability_WS[1]['OHS'])

                    ancestor_S5F_FR = str(ancestor_mutability_aggregated['FR_mutability']['mean_S5F'])
                    ancestor_7M_FR = str(ancestor_mutability_aggregated['FR_mutability']['mean_7M'])
                    ancestor_HS_FR = str(ancestor_mutability_aggregated['FR_mutability']['HS'])
                    ancestor_CS_FR = str(ancestor_mutability_aggregated['FR_mutability']['CS'])
                    ancestor_OHS_FR = str(ancestor_mutability_aggregated['FR_mutability']['OHS'])
        
                    ancestor_S5F_CDR = str(ancestor_mutability_aggregated['CDR_mutability']['mean_S5F'])
                    ancestor_7M_CDR = str(ancestor_mutability_aggregated['CDR_mutability']['mean_7M'])
                    ancestor_HS_CDR = str(ancestor_mutability_aggregated['CDR_mutability']['HS'])
                    ancestor_CS_CDR = str(ancestor_mutability_aggregated['CDR_mutability']['CS'])
                    ancestor_OHS_CDR = str(ancestor_mutability_aggregated['CDR_mutability']['OHS'])
                
                    # Prepare new line of output file:
                    output_line = clone + '_' + prior_abbrev + '_run1' + chain + ','
                    
                    output_line += ancestor_sequence + ','
                
                    output_line += ','.join([ancestor_S5F_WS, ancestor_7M_WS, ancestor_HS_WS, ancestor_CS_WS, ancestor_OHS_WS,
                                             ancestor_S5F_FR, ancestor_7M_FR, ancestor_HS_FR, ancestor_CS_FR, ancestor_OHS_FR,
                                             ancestor_S5F_CDR, ancestor_7M_CDR, ancestor_HS_CDR, ancestor_CS_CDR, ancestor_OHS_CDR,
                                             ])
                                         
                    output_file.write(output_line + '\n')
                    
if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)                    
