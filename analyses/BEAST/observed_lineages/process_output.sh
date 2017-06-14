#!/bin/bash
# First argument: clone name e.g. CH103, CH103L, VRC01_H08...
clone=$1

# Second argument: coalescent prior: 'constant', 'exponential' or 'logistic'
prior=$2

# Abbreviation of prior
prior_abbreviation=$(echo $prior | cut -c1-3)

# Go to results directory:
cd ../../../results/BEAST/observed_lineages/${clone}_${prior}/

#For each chain...
for chain in {a,b,c,d}
do
    chain_id=${clone}_${prior_abbreviation}_run1${chain}
    has_info=$(cat ../burnins_observed.csv | grep $chain_id | wc -l)

    # ...summarize results if the chain is listed in the burn-ins file:
    if [ $has_info -eq 1 ]
        then            
            # Go to chain directory:
            cd $chain_id

            # Make temporary copies of .log and .trees files:
            cp $chain_id.log temp_log.log
            cp $chain_id.trees temp_trees.trees
            cp $chain_id.${clone}_final_alignment.dNdS.log temp_dNdS.log
            
            # Create names for processed tree files
            tree_sample_file=${chain_id}_tree_sample.trees
            MCC_tree_file=${chain_id}_MCC_tree.tree
            
            # Get chain burn-in and length from burnins.csv file:
            length=$(cat ../../burnins_observed.csv | grep $chain_id | cut -d ',' -f 2)
            burnin=$(cat ../../burnins_observed.csv | grep $chain_id | cut -d ',' -f 3)
            
            # Get size of the NEXUS file header in (# lines)
            header_size=$(cat ../../burnins_observed.csv | grep $chain_id | cut -d ',' -f 4)
            
            # Check that length in burnins.csv is the same as obtained from the log file:
            length_from_file=$(tail -1 temp_log.log | grep -o '[0-9]*' | head -1)
            
            if [ $length -ne $length_from_file ]
            then
                echo "Error: Chain length in burnins.csv is incorrect"
                exit
            fi
            
            # Calculate number of steps post-burn-in
            post_burnin=$(($length - $burnin))
            
            # Summarize logs
            ../../../../../../beast1/bin/loganalyser -burnin $burnin temp_log.log log_summary.csv
            #../../../../../../beast1/bin/loganalyser -burnin $burnin temp_dNdS.log dNdS_summary.csv

            # Find sampling interval (in # states) required to get 1000 trees from post-burn-in
            sampling_interval=$(($post_burnin / 1000))

            # Find sampling interval in # lines. Note there's a tree every 10,000 states
            sampling_interval=$(($sampling_interval / 10000))

            # Find first line to be sampled. line # for last burn-in state + header + 1 (since there's a state 0) 
            first_line=$(($((burnin / 10000)) + $header_size + 1))
            
            # Take a sample of 1000 trees 
            head -$header_size temp_trees.trees > $tree_sample_file
            sed -n "$first_line~${sampling_interval}p" temp_trees.trees >> $tree_sample_file 
            echo 'End;' >> $tree_sample_file

            # Generate MCC trees from subsample of 1000 trees (no burn-in necessary since it was already skipped)
            ../../../../../../beast1/bin/treeannotator -burnin 0 -heights mean $tree_sample_file $MCC_tree_file

            # Remove temporary files
            rm temp_log.log
            rm temp_dNdS.log
            rm temp_trees.trees
            
            # Move up to the clone_prior directory
            cd ..
    fi
done
