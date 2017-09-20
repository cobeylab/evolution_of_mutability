'''Simulates B cell evolution using the model described in simulations.tex.
    Does the following steps:
    1 - Imports mutability functions from analyses/mutability folder.
    2 - Imports class BCR from modules/classBCR.py
    3 - Reads a control file (from path in argv[1]) with simulation parameters.
    4 - Defines functions for computing site-specific mutation rates and transition probabilities. Does so based on the "mutability_model" argument read from the control file.
    5 - Simulates B cell evolution passing site specific mutation rates and transition probabilities to the "reproduce" method of BCR objects.
''' 
#!/usr/bin/python
import sys
import re
import csv
from Bio.Seq import Seq 
from Bio.Alphabet import generic_dna
import copy
import importlib
from numpy import random


# Add control_files folder to path
#sys.path.insert(0, 'control_files/')

# Import mutability functions from analyses/mutability folder
sys.path.insert(0, '../mutability/')

sys.path.insert(0,'control_files')

from mutability_function import compute_mean_S5F, count_WRCH, count_DGYW, S5F

# Import BCR class and genetic code from modules folder
import modules.classBCR
BCR = modules.classBCR.BCR

from modules.gen_code_DNA import genetic_code

# Reading control file from argv...
def main(argv):
    control_file =  str(argv[1])
    control_file = control_file.replace('.py','')
    control_file = control_file.replace('control_files/','')
    
    control_file = importlib.import_module(control_file)
        
        
    #=====IMPORTED PARAMETERS FROM CONTROL FILE (SEE README FOR DESCRIPTIONS)======
    # population size, sequence length, number of generations to be simulated 
    pop_size, seq_length, n_gen = control_file.pop_size, control_file.seq_length, control_file.n_gen
    
    # Reference sequence:
    initial_seq = control_file.initial_seq
    
    # Growth rate
    r = control_file.r
    
    # Carrying capacity
    K = control_file.K
     
    # sampling_times, sample size
    sampling_times, sample_size = control_file.sampling_times, control_file.sample_size
    
    # mutability_model, time brakes
    mutability_model, time_breaks = control_file.mutability_model, control_file.time_breaks

    # mutation rate list, reference mutation rate, h
    mutation_rate_list, reference_mutation_rate, h = control_file.mutation_rate_list, control_file.reference_mutation_rate, control_file.h
  
    # fitness cost list, output file path:
    fitness_cost_list, output_filepath = control_file.fitness_cost_list, control_file.output_filepath

    if mutation_rate_list is not None:
        assert(len(time_breaks) == len(mutation_rate_list)), "Number of mutation rate values does not match number of time intervals"
    
    assert(len(time_breaks) == len(fitness_cost_list)), "Number of fitness cost values does not match number of time intervals"
    
    assert(max(sampling_times) == n_gen), "Last sampling time does not match number of generations"
    
    if initial_seq is 'random':
        assert(seq_length is not None), 'Sequence length must be specified for a random reference sequence.'
    else:
        assert(seq_length is None), 'seq_length must be set to None if a reference sequence is directly specified.'
        
    if type(pop_size) is int:
        assert(r is None and K is None), 'Logistic growth parameters (r and K) must be set to None if population size is constant (i.e. pop_size is an integer in the control file).'
    elif pop_size is 'logistic':
        assert((type(r) is float or type(r) is int) and type(K) is int), 'Logistic growth parameters incorrectly specified.' 

    # If sequence was directly specified, find its length
    if initial_seq is not 'random':
        seq_length = len(initial_seq)
    

    # ===============FIND LENGTHS OF TIME INTERVALS FROM TIME_BREAKS======================
    interval_lengths = [time_breaks[0]]
    
    if len(time_breaks) > 1:
        interval_lengths = interval_lengths + [time_breaks[i] - time_breaks[i-1] for i in range(1,len(time_breaks))]
  
    
    #================= DEFINE FUNCTIONS FOR SITE-SPECIFIC MUTATION RATES =================
    # -------------- (only if using "hotspots" or "S5F" mutability models) ---------------

    if mutability_model == "hotspots" or mutability_model == "hotspots-S5F":
        def site_specific_rates(seq):
            '''Checks if each site is a WRCH or a DGYW hotspot given an input sequence. If it is one OR the other, assign it mutability h, else assign it mutability 1. From site-specific mutabilities, the mean mutability of the reference sequence and the reference mutation rate, return a list of site-specific mutation rates (see equation for mu_i,j in model description). Note that h, the mean mutability of the reference sequence and the reference mutation rate are not input as arguments and are instead read from global names. h and the reference mutation rate are read from the control file, and the mean mutability of the reference sequence is calculated at the beginning of the simulation.
            '''
            # Sites at the edges will have mutability 1
            is_WRCH = [0] * len(seq)
            is_DGYW = [0] * len(seq)
            
            # Check if each site (if applicable) is at the center of WRCH or DGYW
            for site in range(2, (len(seq)-1)):
                is_WRCH[site] = count_WRCH(seq[site-2:site+2])
            for site in range(1, (len(seq)-2)):
                is_DGYW[site] = count_DGYW(seq[site-1:site+3])
        
            site_specific_mutabilities = [h if (is_WRCH[i] is 1 or is_DGYW[i] is 1) else 1.0 for i in range(len(seq))]             
            
            # Equation for mu_i,j in model description:
            rates = [(float(m) / mean_reference_mutability) *
            reference_mutation_rate for 
            m in site_specific_mutabilities]

            return(rates)    

    elif mutability_model == "S5F":
        def site_specific_rates(seq):
            '''Finds the mutability of each site using the S5F mutability model. From site-specific mutabilities, the mean mutability of the reference sequence and the reference mutation rate, return a list of site-specific mutation rates (see equation for mu_i,j in model description). Note that h, the mean mutability of the reference sequence and the reference mutation rate are not input as arguments and are instead read from global names. h and the reference mutation rate are read from the control file, and the mean mutability of the reference sequence is calculated at the beginning of the simulation.
            '''    
            # Sites at the edges (for which S5F can't be computed) will have mutability 0
            site_specific_mutabilities = [0] * len(seq)
            for site in range(2, (len(seq)-2)):
                site_specific_mutabilities[site] = S5F[seq[site-2:site+3]]

            # Equation for mu_i,j in model description:
            rates = [(float(m) / mean_reference_mutability) * reference_mutation_rate for m in site_specific_mutabilities]

        
            return(rates)
    else:
        site_specific_rates = None


    #=========== DEFINE FUNCTIONS FOR SITE-SPECIFIC TRANSITION PROBABILITIES =============
    if mutability_model == "uniform" or mutability_model == "hotspots":
        def site_specific_transitions(seq):
            '''Returns a list of transition probabilities for each site. Assigns equal transition probabilities to all *alternative* bases
            '''
            transitions_list = []
            for nucleotide in seq:
                initial_nt = nucleotide
                site_transitions = {}
                for base in ['A','T','C','G']:
                    if base is initial_nt:
                        site_transitions[base] = 0
                    else:
                        site_transitions[base] = 1./3
                transitions_list.append(site_transitions)
                
            # As a test, check that the probability of a site mutating to the current nucleotide is zero for all sites:
            check = [transitions_list[i][seq[i]] == 0 for i in range(len(seq))]
            assert check.count(False) == 0, 'Error in site specific transitions: non-zero probability of mutating to the current nucleotide'
            return transitions_list

    elif mutability_model == "S5F" or mutability_model == 'hotspots-S5F':
       
        #Import conditional transition probabilities for each five-mer (Yaari et al.):
        S5F_transitions = {}
        with open('modules/Yaari_transition_probabilities.csv', 'r') as csvfile:
            csvread = csv.reader(csvfile, delimiter=' ', quotechar='|')
            csvread.next()
            for row in csvread:
                csvrow = row
                motif = csvrow[0].replace('"', '')
                S5F_transitions[motif] = {}
                S5F_transitions[motif]["A"] = float(csvrow[1].replace('"', ''))
                S5F_transitions[motif]["C"] = float(csvrow[2].replace('"', ''))
                S5F_transitions[motif]["G"] = float(csvrow[3].replace('"', ''))
                S5F_transitions[motif]["T"] = float(csvrow[4].replace('"', ''))  
       
        def site_specific_transitions(seq):
            '''Returns a list of transition probabilities for each site. Assigns transition probabilities based on five-mer centered at the site'''
            # Sites for which five-mers can't be determined (edges) are given uniform transition probs
            # (only relevant for the hybrid hotspots-S5F model, otherwise those sites have 0 mutation rate
            transitions_list = [{}] * len(seq)
            for site in range(len(seq)):
                # If site can be assigned a fivemer motif
                if site in range(2, (len(seq)-2)):
                    site_transitions = S5F_transitions[seq[site-2:site+3]]
                else:
                # Assign uniform transition probabilities to the first 2 and last 2 sites
                # (For the hybrid case w HS-based rates and S5F-based transitions)
                    initial_nt = seq[site]
                    site_transitions = {}
                    for base in ['A','T','C','G']:
                        if base is initial_nt:
                            site_transitions[base] = 0
                        else:
                            site_transitions[base] = 1./3
                
                transitions_list[site] = site_transitions
                    
             # As a test, check that the probability of a site mutating to the current nucleotide is zero for all sites:
            check = [transitions_list[i][seq[i]] == 0 for i in range(len(seq)) if len(transitions_list[i].keys()) > 0]
            assert check.count(False) == 0, 'Error in site specific transitions: non-zero probability of mutating to the current nucleotide' 
            return transitions_list
        

    #============ CONVERT PARAMETERS SPECIFIED AS STRINGS TO NUMERICAL VALUES ============    
     
    # If mutability model is "uniform"... 
    if mutability_model == "uniform":
    
        # Converting elements in mutation_rate_list to numerical values
        # mutation_rate_list is 'None' under the 'hotspot' and 'S5F' models
        for i in range(len(mutation_rate_list)):
            if type(mutation_rate_list[i]) is str:
               # Find multiplier
               assert mutation_rate_list[i].find('*1/4L') > -1, 'One or more mutation rates are incorrectly specified in "mutation_rate_list" in the control file'
               multiplier = re.search(r'[0-9|.]+', mutation_rate_list[i])
               
               assert multiplier.group() is not None, 'Multiplier for 1/4L in mutation rate must be a number'
               multiplier = float(multiplier.group())
           
               mutation_rate_list[i] = multiplier * 1./(4*seq_length)
               
    # If mutability model is "uniform"... 
    
    elif(mutability_model == "hotspots" or mutability_model == "S5F" or mutability_model == "hotspots-S5F"):
 
        # Converting reference mutation rate to numerical value (reference_mutation_rate is None under the 'uniform' mutability model)
        if type(reference_mutation_rate) is str:
            assert reference_mutation_rate.find('*1/4L') > -1, 'Reference mutation rate is incorrectly specified in control file'
            multiplier = re.search(r'[0-9|.]+', reference_mutation_rate)
            assert multiplier.group() is not None, 'Multiplier for 1/4L in reference mutation rate must be a number'
            multiplier = float(multiplier.group())
        
            print "rate multiplier: " + str(multiplier)
            reference_mutation_rate = multiplier * 1./(4*seq_length)
      
    #============ GENERATE LISTS WITH PARAMETER VALUES FOR EACH GENERATION ============
    #E.g. if 1st interval lasts 100 gen., 1st value of fitness_cost_list will be repeated 100 times in the list.
    fitness_cost_list = [[fitness_cost_list[i]] * interval_lengths[i] for i in
                         range(len(fitness_cost_list))]
    fitness_cost_list = [x for y in fitness_cost_list for x in y]
   
    # Same for mutation_rate_list (under the "uniform" model)
    if mutation_rate_list is not None:    
        mutation_rate_list = [[mutation_rate_list[i]] * interval_lengths[i] for i in 
                             range(len(mutation_rate_list))]
        mutation_rate_list = [x for y in mutation_rate_list for x in y]

      
    #================= INITIALIZE REFERENCE SEQUENCE AND POPULATION ======================
    if initial_seq == 'random':
        # Generate random reference sequence:
        # ----------- (no stop codons, all other codons sampled uniformly) ---------------- 

        nonstop_codons = [key for key in genetic_code.keys() if genetic_code[key] != '*']

        initial_seq = random.choice(nonstop_codons, seq_length/3)
        initial_seq = ''.join(initial_seq)
        assert(len(initial_seq) == seq_length), "Reference sequence length does not match seq_length"
    


    # Compute mutability of reference sequence, if under 'hotspots' or 'S5F' mut. model
    if mutability_model == "hotspots" or mutability_model == "hotspots-S5F":
        # each hotspot site has mutability h, each non-hotspot site has mutability 1
        ref_n_hotspots = count_WRCH(initial_seq) + count_DGYW(initial_seq)
        # Total summed mutability
        mean_reference_mutability = ref_n_hotspots * h + (seq_length - ref_n_hotspots)
        # Divided by number of sites
        mean_reference_mutability = float(mean_reference_mutability) / seq_length
        
    elif mutability_model == "S5F":
        #Sites at the edges for which S5F can't be computed are excluded
        mean_reference_mutability = compute_mean_S5F(initial_seq) 
        # compute_mean_S5F excludes the 2 nucleotides at each edge from the mean
        # adjusting value to take them into account with mutability 0 (in line with the site_specific_rates function)
        mean_reference_mutability = mean_reference_mutability * float(seq_length - 4)/seq_length

    # As a test of site_specific_rates, check that applying it to the reference sequence returns reference_mutation_rate
    if site_specific_rates is not None:
        error = (sum(site_specific_rates(initial_seq))/seq_length - reference_mutation_rate) / reference_mutation_rate
        assert abs(error) < 10 ** -10, "site_specific_rates function does not recover reference mutation rate when applied to reference sequence"

    
    # Initialize dictionary with output sequences
    # Each key is a time point, entries are lists of sequences sampled at time points
    output_sequences = {str(sampling_times[i]):[] for i in range(len(sampling_times))}


    #================================ SIMULATION =========================================

    # Initialize population 

    # If population is constant in size
    if type(pop_size) is int:
        r = 0
        K = 1
        
        # Initialize generation 1 from initial sequence, all with identical fitness
        parent_sequences = [BCR(initial_seq)] * pop_size
        parent_fitness = [1./pop_size] * pop_size

    # If pop_size is 'logistic' initialize population with 1 individual
    elif pop_size == 'logistic':
        pop_size = 1
        
        # Initialize one individual with reference sequence and fitness of 1
        parent_sequences = [BCR(initial_seq)]
        parent_fitness = [1] 


    # For each generation
    # Adjusting numbers so that initial population is labeled as generation 1
    for generation in range(1, n_gen + 1):
        
        if generation % 100 == 0:
            print "Generation " + str(generation)
    
        child_sequences = []
        child_fitness = []
        
        # Find population size of next generation
        next_pop_size = pop_size + pop_size * r * (1 - float(pop_size)/K)
        
        # Get parameter values for the time interval containing the generation
        # (Adjusting numbers so that values for generation 1 come from index 0, and so on)
        fitness_cost = fitness_cost_list[generation - 1]

        if mutability_model == 'uniform':
            mutation_rate = mutation_rate_list[generation - 1]

        # For each new sequence (rounding pop size to an integer)
        for i in range(int(round(next_pop_size))):
            # Randomly pick parent based on fitness of parental sequences
            # (Sample index from 0 to the nearest integer to the current pop. size)
            parent_index = random.choice(range(int(round(pop_size))),1, p = parent_fitness)[0]
            parent = parent_sequences[parent_index]
                
            # Find site-specific rates for parent sequence:
            if mutability_model == 'uniform':
                site_rates = [mutation_rate] * seq_length
            
            elif mutability_model == 'S5F' or mutability_model == 'hotspots' or mutability_model == 'hotspots-S5F':
                site_rates = site_specific_rates(parent.sequence)
 
            # Find site-specific transition probabilities for parent sequence:
            site_transitions = site_specific_transitions(parent.sequence)
        
            # Generate child from parent BCR object, introducing mutations
            child_sequences.append(parent.reproduce(site_rates, site_transitions))

            # Base fitness
            fitness = 1
        
            # Adjust child's fitness based on non-synonymous mutations relative to initial sequence:
        
            # Find sites that differ between new sequence and the INITIAL sequence
            dif_sites = [j for j in range(seq_length) if child_sequences[i].sequence[j] != initial_seq[j]]
        
            # For each site that differs:
            for site in dif_sites:

                # deliberately not a float division. Find codon position (first codon is position 0).
                codon_position = site/3
                initial_codon = initial_seq[codon_position*3 : codon_position*3 + 3]
                mutated_codon = child_sequences[i].sequence[codon_position*3 : codon_position*3 + 3]
            
                # Compare amino acids. If stop codon, make fitness 0, if non-syn., impose fitness cost
                if genetic_code[''.join(mutated_codon)]!=genetic_code[''.join(initial_codon)]:
                    # Assigning fitness 0 to sequences with stop codons
                    if genetic_code[''.join(mutated_codon)] == '*':
                        fitness = 0
                    else:
                        fitness += fitness_cost
                
                
            # If fitness is negative, set it to 0
            if fitness < 0:
                fitness = 0
            
            child_fitness.append(fitness)

        # Normalize fitness
        assert(sum(child_fitness) > 0), "All sequences have fitness zero. Mutation rate is probably too high for the chosen fitness cost of non-syn. mutations"
        
        child_fitness = [float(child_fitness[k]) / sum(child_fitness) for k in range(int(round(next_pop_size)))]
        
        # Set child sequences as parents for next generation
        parent_sequences = copy.copy(child_sequences)
        parent_fitness = copy.copy(child_fitness)    

        # Set population size at next generation as the current population size (for next cycle)
        pop_size = copy.copy(next_pop_size)
            
        # If generation is in sampling_times, record random sample of sequences
        if (generation) in sampling_times:
            #Sample only sequences without stop codons:
            viable_sequences = [seq for seq in child_sequences if seq.translate().find('*') == -1]   
            output_sequences[str(generation)] = random.choice(viable_sequences, sample_size)   
            print "Recording sequences (Generation " + str(generation) + ")"
         
    # Outputting sequences:

    # Export alignment of observed sequences (and germline) as a nexus file:
    with open(output_filepath,"w") as alignment_file:
        alignment_file.write("#NEXUS\n\n")
        alignment_file.write("BEGIN DATA;\n")
        alignment_file.write("\tDIMENSIONS NTAX=" + str(sample_size * len(sampling_times) + 1) + " NCHAR=" +
                                 str(seq_length) + ";\n")
        alignment_file.write("\tFORMAT DATATYPE=DNA\n")
        alignment_file.write("\tGAP=-\n")
        alignment_file.write("\t;\n")
        alignment_file.write("MATRIX\n")
        
        alignment_file.write('[1] GERMLINE_00\n')
        alignment_file.write(initial_seq + '\n')
    
        nexus_number = 2
        for time_point in output_sequences.keys():
            for seq_number in range(len(output_sequences[time_point])):
                 alignment_file.write("[" + str(nexus_number) + "] " + 'S' + str(seq_number) + '_' + time_point + "\n")
                 alignment_file.write(output_sequences[time_point][seq_number].sequence + '\n')
                 nexus_number += 1

        alignment_file.write(";\n")
        alignment_file.write("END;")   
    
if(__name__== "__main__"):
    status = main(sys.argv)
    sys.exit(status)