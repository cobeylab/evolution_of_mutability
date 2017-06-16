#!/bin/bash

# Running divergence from ancestors script for each heavy and light chain:

# CH103 heavy chain
echo Analyzing CH103 heavy chain
python divergence_from_ancestors.py CH103_con_run1a

# CH103 light chain
echo Analyzing CH103 light chain
python divergence_from_ancestors.py CH103L_con_run1a

# VRC01-13 (heavy)
echo Analyzing VRC01-13 heavy chain
python divergence_from_ancestors.py VRC01_13_log_run1a

# VRC01-01 (heavy)
echo Analyzing VRC01-01 heavy chain
python divergence_from_ancestors.py VRC01_01_log_run1a

# VRC01_19 (heavy)
echo Analyzing VRC01-19 heavy chain
python divergence_from_ancestors.py VRC01_19_log_run1a

# VRC26 heavy chain
echo Analyzing VRC26 heavy chain
python divergence_from_ancestors.py VRC26int_con_run1a

# VRC26 light chain
echo Analyzing VRC26 light chain
python divergence_from_ancestors.py VRC26L_con_run1a



