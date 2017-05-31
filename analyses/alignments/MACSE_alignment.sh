#!/bin/bash

# Concatenate germline sequences:
echo "- Preparing fasta files"
python prepare_CH103_for_MACSE.py
python prepare_CH103L_for_MACSE.py
python prepare_VRC01_01_for_MACSE.py
python prepare_VRC01_13_for_MACSE.py
python prepare_VRC01_19_for_MACSE.py
python prepare_VRC26_for_MACSE.py
python prepare_VRC26L_for_MACSE.py


# Align CH103:
echo "=========================Aligning CH103 sequences...================================="
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/CH103_plus_GERMLINE.fasta -../../results/alignments/CH103_MACSE_NT.fasta -../../results/alignments/CH103_MACSE_AA.fasta
mv ../../results/alignments/CH103_plus_GERMLINE_macse_AA.fasta ../../results/alignments/CH103_MACSE_AA.fasta
mv ../../results/alignments/CH103_plus_GERMLINE_macse_NT.fasta ../../results/alignments/CH103_MACSE_NT.fasta
rm ../../results/alignments/CH103_plus_GERMLINE.fasta

# Align CH103L:
echo "=========================Aligning CH103L sequences...================================"
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/CH103L_plus_GERMLINE.fasta -../../results/alignments/CH103L_MACSE_NT.fasta -../../results/alignments/CH103L_MACSE_AA.fasta
mv ../../results/alignments/CH103L_plus_GERMLINE_macse_AA.fasta ../../results/alignments/CH103L_MACSE_AA.fasta
mv ../../results/alignments/CH103L_plus_GERMLINE_macse_NT.fasta ../../results/alignments/CH103L_MACSE_NT.fasta
rm ../../results/alignments/CH103L_plus_GERMLINE.fasta

# Align VRC01_01:
echo "=========================Aligning VRC01_01 sequences...================================"
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/VRC01_01_plus_GERMLINE.fasta -../../results/alignments/VRC01_01_MACSE_NT.fasta -../../results/alignments/VRC01_01_MACSE_AA.fasta
mv ../../results/alignments/VRC01_01_plus_GERMLINE_macse_AA.fasta ../../results/alignments/VRC01_01_MACSE_AA.fasta
mv ../../results/alignments/VRC01_01_plus_GERMLINE_macse_NT.fasta ../../results/alignments/VRC01_01_MACSE_NT.fasta
rm ../../results/alignments/VRC01_01_plus_GERMLINE.fasta

# Align VRC01_13:
echo "=========================Aligning VRC01_13 sequences...================================"
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/VRC01_13_plus_GERMLINE.fasta -../../results/alignments/VRC01_13_MACSE_NT.fasta -../../results/alignments/VRC01_13_MACSE_AA.fasta
mv ../../results/alignments/VRC01_13_plus_GERMLINE_macse_AA.fasta ../../results/alignments/VRC01_13_MACSE_AA.fasta
mv ../../results/alignments/VRC01_13_plus_GERMLINE_macse_NT.fasta ../../results/alignments/VRC01_13_MACSE_NT.fasta
rm ../../results/alignments/VRC01_13_plus_GERMLINE.fasta

# Align VRC01_19:
echo "=========================Aligning VRC01_19 sequences...================================"
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/VRC01_19_plus_GERMLINE.fasta -../../results/alignments/VRC01_19_MACSE_NT.fasta -../../results/alignments/VRC01_19_MACSE_AA.fasta
mv ../../results/alignments/VRC01_19_plus_GERMLINE_macse_AA.fasta ../../results/alignments/VRC01_19_MACSE_AA.fasta
mv ../../results/alignments/VRC01_19_plus_GERMLINE_macse_NT.fasta ../../results/alignments/VRC01_19_MACSE_NT.fasta
rm ../../results/alignments/VRC01_19_plus_GERMLINE.fasta

# Align VRC26:
echo "=========================Aligning VRC26 sequences...================================="
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/VRC26_plus_GERMLINE.fasta -../../results/alignments/VRC26_MACSE_NT.fasta -../../results/alignments/VRC26_MACSE_AA.fasta
mv ../../results/alignments/VRC26_plus_GERMLINE_macse_AA.fasta ../../results/alignments/VRC26_MACSE_AA.fasta
mv ../../results/alignments/VRC26_plus_GERMLINE_macse_NT.fasta ../../results/alignments/VRC26_MACSE_NT.fasta
rm ../../results/alignments/VRC26_plus_GERMLINE.fasta

# Align VRC26L:
echo "=========================Aligning VRC26L sequences...================================"
java -jar ~/MACSE.jar -prog alignSequences -seq ../../results/alignments/VRC26L_plus_GERMLINE.fasta -../../results/alignments/VRC26L_MACSE_NT.fasta -../../results/alignments/VRC26L_MACSE_AA.fasta
mv ../../results/alignments/VRC26L_plus_GERMLINE_macse_AA.fasta ../../results/alignments/VRC26L_MACSE_AA.fasta
mv ../../results/alignments/VRC26L_plus_GERMLINE_macse_NT.fasta ../../results/alignments/VRC26L_MACSE_NT.fasta
rm ../../results/alignments/VRC26L_plus_GERMLINE.fasta

# Remove files with amino acid sequences output by MACSE
rm ../../results/alignments/*_AA.fasta
