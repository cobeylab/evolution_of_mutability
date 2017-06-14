#!/bin/bash
#SBATCH --job-name=CH103L_con_run1d_S_NS_mutability_changes
#SBATCH --output=CH103L_con_run1d_S_NS_mutability_changes.out
#SBATCH --error=CH103L_con_run1d_S_NS_mutability_changes.err
#SBATCH --time=600:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

cd ../../../
python S_NS_mutability_changes.py CH103L_con_run1d
