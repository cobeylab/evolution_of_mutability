pop_size = 'logistic'
initial_seq = 'TCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGTAGTTACTACTGGAGCTGGATCCGGCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGTGGGAGCACCAACTACAACCCCTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCTGCGGACACGGCCGTGTATTACTGTGCGAGCCTGCCCAGGGGGCAGTTAGTCAATGCCTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
r = 0.01
K = 2000
seq_length = None
n_gen = 2000
sampling_times = [500,750,1000,1250,1500,1750,2000]
sample_size = 25
mutability_model = 'uniform'
time_breaks = [500,750,1000,1250,1500,1750,2000]
mutation_rate_list = ['1*1/4L', '1*1/4L', '0.96*1/4L', '0.92*1/4L', '0.88*1/4L', '0.84*1/4L', '0.8*1/4L']
reference_mutation_rate = None
h = None
fitness_cost_list = [-0.005, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005]
output_filepath = '../../results/simulated_alignments/scenario4b/scenario4b_rep1.nex'