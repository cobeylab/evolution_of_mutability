background_and_model_all_evenodd.mat has the binding energies (in object h) and background 7mer frequencies (f_0_mer) necessary for computing the mutability of each 7mer.

running prob_Lmer.m:

prob_Lmer(h, f_0_mer)

outputs the mutability of all 16384 7mers.

To identify the 7mers, note that the order of 7mers in the f_0_mer vector (and in the output of prob_Lmer) is the same as in the rows of matrix ii generated internally by prob_Lmer. Each row of ii has seven elements, each element is a number. 1 = A, 2 = C, 3 = G, 4 = T.

To obtain mutability scores independent of background 7mer frequencies, modify the line multiplying f_0_mer in the prob_Lmer code.
