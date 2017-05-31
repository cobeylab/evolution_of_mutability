#!/bin/bash
#SBATCH --job-name=VRC01_13_exponential_process_output
#SBATCH --output=VRC01_13_exponential_process_output.out
#SBATCH --error=VRC01_13_exponential_process_output.err
#SBATCH --time=2:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

module load java
cd /project/cobey/mvieira/evolvability/analyses/BEAST/observed_lineages
./process_output.sh VRC01_13 exponential
