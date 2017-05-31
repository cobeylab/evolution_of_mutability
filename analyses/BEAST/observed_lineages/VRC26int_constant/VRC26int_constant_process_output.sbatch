#!/bin/bash
#SBATCH --job-name=VRC26int_constant_process_output
#SBATCH --output=VRC26int_constant_process_output.out
#SBATCH --error=VRC26int_constant_process_output.err
#SBATCH --time=2:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16000

module load java
cd /project/cobey/mvieira/evolvability/analyses/BEAST/observed_lineages
./process_output.sh VRC26int constant
