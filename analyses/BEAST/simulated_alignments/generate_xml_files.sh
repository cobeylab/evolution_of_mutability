#!/bin/bash

# For each scenario:
for scenario in {1,2a,2b,2c,2d,2e,3a,3b,3c,3d}
do
    # For each replicate
    for replicate in {1..4}
    do
        echo ""  
        echo "Input:"
        echo ../../../results/simulated_alignments/scenario${scenario}/scenario${scenario}_rep${replicate}.nex
        echo "Output:"
        echo scenario${scenario}/scenario${scenario}_rep${replicate}/scenario${scenario}_rep${replicate}.xml
        echo ""
        ~/BEASTGen/bin/beastgen simulated_alignment_template.xml ../../../results/simulated_alignments/scenario${scenario}/scenario${scenario}_rep${replicate}.nex scenario${scenario}_rep${replicate}/scenario${scenario}_rep${replicate}.xml
    done    
done
