#!/bin/bash
#SBATCH --job-name=VRC01_13_con_run1d_S_NS_mutability_changes
#SBATCH --output=VRC01_13_con_run1d_S_NS_mutability_changes.out
#SBATCH --error=VRC01_13_con_run1d_S_NS_mutability_changes.err
#SBATCH --time=600:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

cd ../../../
python S_NS_mutability_changes.py VRC01_13_con_run1d
