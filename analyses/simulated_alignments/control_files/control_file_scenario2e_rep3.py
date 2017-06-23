pop_size = 'logistic'
initial_seq = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCCTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
r = 0.7
K = 1000
seq_length = None
n_gen = 2000
sampling_times = [500,750,1000,1250,1500,1750,2000]
sample_size = 25
mutability_model = 'uniform'
time_breaks = [500,750,1000,1250,1500,1750,2000]
mutation_rate_list = ['1*1/4L', '1*1/4L', '0.9*1/4L', '0.8*1/4L', '0.7*1/4L', '0.6*1/4L', '0.5*1/4L']
reference_mutation_rate = None
h = None
fitness_cost_list = [-0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01]
output_filepath = '../../results/simulations/scenario2e/scenario2e_rep3.nex'