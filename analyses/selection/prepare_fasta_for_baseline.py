"""Prepares fasta files for each lineage to be analyzed with Baseline (Yaari et al. 2012 Nuc. Acid. Res)
"""
import sys
import re
from copy import deepcopy
from dendropy import Tree

chain_ids = ['CH103_con_run1a','CH103L_con_run1a','VRC26int_con_run1a',
             'VRC26L_con_run1a','VRC01_01_log_run1a','VRC01_13_log_run1a',
             'VRC01_19_log_run1a']
def main(argv):
    for chain_id in chain_ids:

        # GET XML FILES AND MCC TREE
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

        chain_directory = '../../results/BEAST/observed_lineages/' + clone + '_' + prior + '/' + chain_id + '/'
        MCC_tree_file_path = chain_directory + chain_id + '_MCC_tree.tree'

        xml_file_path = '../../analyses/BEAST/observed_lineages/' + clone + '_' + prior + '/'
        xml_file_path += chain_id[0:len(chain_id) - 1] + '.xml'

        output_path = 'baseline_fasta_files/' + chain_id + '_baseline.fasta'

        fasta_string = ''


        # ============================= READ SEQUENCES FROM XML FILE AND WRITE THEM TO FASTA STRING ==========================
        # Read XML as string
        with open(xml_file_path, 'r') as xml_file:
            xml = xml_file.readlines()
            xml = ''.join(xml)

        taxon_lines = re.findall(r'<sequence>.*</sequence>', xml, re.DOTALL)[0].split('</sequence>')
        for line in taxon_lines[0:len(taxon_lines) - 1]:
            taxon_id = re.search(r'taxon idref=".*"', line).group().replace("taxon idref=", '').replace('"', '')
            # Ignoring a weird VRC26 sequence whose CDR2 is entirely gaps:
            if taxon_id != 'KJ134124_119':

                #Ignore concatenated V + J ("GERMLINE_00) sequence:
                if taxon_id != 'GERMLINE_00':
                    sequence = re.search(r'/>\n\t\t\t.*\n\t\t', line).group().replace('/>\n\t\t\t', '').replace('\n\t\t', '')

                    fasta_string += '>' + taxon_id + '\n' + sequence + '\n'

        # ==================== READ COMMON ANCESTOR SEQUENCE FROM MCC TREE AND WRITE IT TO FASTA STRING ====================
        with open(MCC_tree_file_path, 'r') as tree_file:

            # Skipping file header with NEXUS specifications. Header ends at the end of 'Translate' block, at the 6th ';'
            n_semicolons = 0
            while n_semicolons < 6:
                line = tree_file.next()

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

            # ====================================== READING TREE WITH DENDROPY ============================================
            # Read tree string as a tree object using dendropy:
            tree = Tree.get_from_string(tree_string, schema='newick')

            # Get root node sequence:
            node = tree.seed_node

            node_CP1 = node.annotations.get_value(alignment_id + '.CP1')
            node_CP2 = node.annotations.get_value(alignment_id + '.CP2')
            node_CP3 = node.annotations.get_value(alignment_id + '.CP3')

            # Arbitrarily choosing first listed reconstruction in case of ambiguity
            if node_CP1.find('+') > -1:
                # If it is, arbitrarily choose the first sequence
                node_CP1 = node_CP1.split('+')[0]
                print 'Ambiguous reconstruction in ancestral node; arbitrarily choosing first listed reconstruction'
            if node_CP2.find('+') > -1:
                # If it is, arbitrarily choose the first sequence
                node_CP2 = node_CP2.split('+')[0]
                print 'Ambiguous reconstruction in ancestral node; arbitrarily choosing first listed reconstruction'
            if node_CP3.find('+') > -1:
                # If it is, arbitrarily choose the first sequence
                node_CP3 = node_CP3.split('+')[0]
                print 'Ambiguous reconstruction in ancestral node; arbitrarily choosing first listed reconstruction'

            node_sequence = [node_CP1[i] + node_CP2[i] + node_CP3[i] for i in range(len(node_CP1))]
            node_sequence = ''.join(node_sequence)

        fasta_string = '>>UCA' + '\n' + node_sequence + '\n' + fasta_string

        with open(output_path,'w') as output_file:
            output_file.write(fasta_string)

if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)



