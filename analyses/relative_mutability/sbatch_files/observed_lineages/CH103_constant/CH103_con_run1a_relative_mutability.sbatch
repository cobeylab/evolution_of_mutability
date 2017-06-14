#!/bin/bash
#SBATCH --job-name=CH103_con_run1a_relative_mutability
#SBATCH --output=CH103_con_run1a_relative_mutability.out
#SBATCH --error=CH103_con_run1a_relative_mutability.err
#SBATCH --time=600:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

cd ../../../
python relative_mutability.py CH103_con_run1a
