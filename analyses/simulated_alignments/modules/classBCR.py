''' This script defines the object class "BCR", initialized with a string representing a nucleotide sequence. The "reproduce" method can be used to replicate the sequence while introducing mutations. This method takes as input a) a list of site-specific mutation rates and b) a list of dictionaries specifying transition probabilities for each site. This input structure is general enough to accomodate any of the three mutability models investigated ("uniform", "hotspots" and "S5F").
'''

import gen_code_DNA as GC
from numpy import random

class BCR(object):

    def __init__(self, sequence):
        self.sequence = sequence
    
    def reproduce(self, site_mutation_probs, site_transition_probs):
        '''Outputs a daughter sequence from a parent sequence, introducing mutations to each site according to site-specific mutation probabilities contained in 'site_mutation_probs'.
            Transition probabilities among bases are input as a list of dictionaries. Each element of the list corresponds to a site, and the dictionary specifies the probabilities that, given a mutation happened at that site, it will change to each base. 
        '''    
        mutated_sequence = list(self.sequence)
    
        for site in range(len(mutated_sequence)):
            if random.uniform() <= site_mutation_probs[site]:
                initial_nt = mutated_sequence[site]
                bases = [key for key in site_transition_probs[site]]
                transition_probs = [float(site_transition_probs[site][base]) for base in bases]
            
                new_nt = random.choice(bases,1, p = transition_probs)[0]
                assert new_nt is not initial_nt, 'Mutated base is equal to initial base; check site transition probabilities'
                
                mutated_sequence[site] = new_nt
          
        return BCR(''.join(mutated_sequence)) 
    
    def translate(self):
        protein_seq = ''
        for site in range(0, len(self.sequence), 3):
            codon = GC.genetic_code[self.sequence[site : site + 3]]
            protein_seq = protein_seq + codon
        return protein_seq   
