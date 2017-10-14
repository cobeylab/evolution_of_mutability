#!/bin/bash
#SBATCH --job-name=contrasts_plot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cobey
#SBATCH -o contrasts_plot.out        
#SBATCH -e contrasts_plot.err      
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=32000
module load R
Rscript contrasts_observed.R
