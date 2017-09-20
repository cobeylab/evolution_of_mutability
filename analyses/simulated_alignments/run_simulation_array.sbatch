#!/bin/bash
#SBATCH --job-name=simulate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o simulate_%A_%a.out        
#SBATCH -e simulate_%A_%a.err        
#SBATCH --time=200:00:00

module load python

# Get control file to execute from its index
index=$SLURM_ARRAY_TASK_ID

control_file=$(ls control_files/*.py | grep -v 'pyc' | head -$index | tail -1)

# Run simulation passing control file
python evolveSequences.py $control_file
