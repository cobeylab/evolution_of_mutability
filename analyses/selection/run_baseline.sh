#!/bin/bash

# Running baseline for each heavy/light chain:
# See Baseline_Main_Version1.3.R for description of arguments
# FR/CDR boundaries format: 1st AA of FR1: Last AA FR1: Last AA CDR1: ... : Last AA CDR3

# Partition points from 'partition_points.py' file

# CH103 heavy chain
echo Analyzing CH103 heavy chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 0 1:11:19:36:43:81:96 ./baseline_fasta_files/CH103_con_run1a_baseline.fasta ../../results/selection/CH103_con_run1a _baseline
mv ../../results/selection/CH103_con_run1a_baseline.txt ../../results/selection/CH103_con_run1a_baseline.tsv

# CH103 light chain
echo Analyzing CH103 light chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:25:31:48:51:87:95 ./baseline_fasta_files/CH103L_con_run1a_baseline.fasta ../../results/selection/CH103L_con_run1a _baseline
mv ../../results/selection/CH103L_con_run1a_baseline.txt ../../results/selection/CH103L_con_run1a_baseline.tsv

# VRC01-13 (heavy)
echo Analyzing VRC01-13 heavy chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:25:33:50:58:96:114 ./baseline_fasta_files/VRC01_13_log_run1a_baseline.fasta ../../results/selection/VRC01_13_log_run1a _baseline
mv ../../results/selection/VRC01_13_log_run1a_baseline.txt ../../results/selection/VRC01_13_log_run1a_baseline.tsv

# VRC01-01 (heavy)
echo Analyzing VRC01-01 heavy chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:25:33:50:58:96:121 ./baseline_fasta_files/VRC01_01_log_run1a_baseline.fasta ../../results/selection/VRC01_01_log_run1a _baseline
mv ../../results/selection/VRC01_01_log_run1a_baseline.txt ../../results/selection/VRC01_01_log_run1a_baseline.tsv

# VRC01_19 (heavy)
echo Analyzing VRC01-19 heavy chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:25:33:50:58:103:143 ./baseline_fasta_files/VRC01_19_log_run1a_baseline.fasta ../../results/selection/VRC01_19_log_run1a _baseline
mv ../../results/selection/VRC01_19_log_run1a_baseline.txt ../../results/selection/VRC01_19_log_run1a_baseline.tsv

# VRC26 heavy chain
echo Analyzing VRC26 heavy chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:27:35:52:60:98:151 ./baseline_fasta_files/VRC26int_con_run1a_baseline.fasta ../../results/selection/VRC26int_con_run1a _baseline
mv ../../results/selection/VRC26int_con_run1a_baseline.txt ../../results/selection/VRC26int_con_run1a_baseline.tsv

# VRC26 light chain
echo Analyzing VRC26 light chain
Rscript Baseline_Main_Version1.3.R 1 1 5 5 0 1 1:25:33:50:53:89:98 ./baseline_fasta_files/VRC26L_con_run1a_baseline.fasta ../../results/selection/VRC26L_con_run1a _baseline
mv ../../results/selection/VRC26L_con_run1a_baseline.txt ../../results/selection/VRC26L_con_run1a_baseline.tsv




