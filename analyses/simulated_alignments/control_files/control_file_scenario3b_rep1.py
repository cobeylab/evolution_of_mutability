pop_size = 'logistic'
initial_seq = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCCTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
r = 0.7
K = 1000
seq_length = None
n_gen = 2000
sampling_times = [500,750,1000,1250,1500,1750,2000]
sample_size = 25
mutability_model = 'hotspots'
time_breaks = [2000]
mutation_rate_list = None
reference_mutation_rate = '1*1/4L'
h = 3
fitness_cost_list = [-0.01]
output_filepath = '../../results/simulated_alignments/scenario3b/scenario3b_rep1.nex'