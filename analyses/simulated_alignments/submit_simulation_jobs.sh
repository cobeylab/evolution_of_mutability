#!/bin/bash

# Create list with file control file indices for all control files
n_control_files=$(ls control_files/*.py | grep -v 'pyc' | wc -l)

control_file_indices=$(for i in $(seq 1 $n_control_files); do echo $i; done | tr '\n' ',')

# Remove trailing comma
control_file_indices=${control_file_indices::-1}

# Run job array with list of ids. The array sbatch file will execute each file by its index
sbatch --array=$control_file_indices run_simulation_array.sbatch