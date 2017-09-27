import re
import sys
sys.path.insert(0, '../mutability/')
from mutability_function import seq_mutability, aggregated_mutability


# Function for getting observed sequences from XML file, returning their mutability
def get_mutability_from_XML(xml_file_path, partition_points):

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

    return {'whole_sequence': mutability_WS_tips, 'aggregated_by_region' : mutability_aggregated_tips,
            'sequences':obs_sequence}

