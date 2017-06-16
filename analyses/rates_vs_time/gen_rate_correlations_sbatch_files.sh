#!/bin/bash

#Sbatch files for observed lineages
for clone in {CH103,CH103L,VRC26,VRC26int,VRC26L,VRC01_01,VRC01_19,VRC01_13}
do
   for prior in {constant,exponential,logistic}
   do
        prior_abbreviation=$(echo $prior | cut -c1-3)
        
        for chain in {a,b,c,d}
        do
            chain_id=${clone}_${prior_abbreviation}_run1${chain}

            echo "#!/bin/bash" > ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --job-name=${chain_id}_rates_vs_time" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --output=${chain_id}_rates_vs_time.out" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --error=${chain_id}_rates_vs_time.err" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --time=200:00:00" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --partition=cobey" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --nodes=1" >> ${chain_id}_rates_vs_time.sbatch
            echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_rates_vs_time.sbatch

            echo ""  >> ${chain_id}_rates_vs_time.sbatch

            echo "cd ../../../"  >> ${chain_id}_rates_vs_time.sbatch
            echo "python rates_vs_time.py $chain_id"  >> ${chain_id}_rates_vs_time.sbatch

            mv ${chain_id}_rates_vs_time.sbatch sbatch_files/observed_lineages/${clone}_${prior}
        
        done
   
   done
done

# Sbatch files for simulated alignments
# 
# for scenario in {1,2a,2b,3,4a,4b,5a,5b,5c,6a,6b,6c,7,8a,8b,8c,8d,8e,9,10a,10b,10c,10d,10e,11a,11b,11c,12a,12b,12c}
# do
#     
#     for replicate in {1,2,3,4}
#     do
#         chain_id=scenario${scenario}_rep${replicate}
# 
#         echo "#!/bin/bash" > ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --job-name=${chain_id}_rates_vs_time" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --output=${chain_id}_rates_vs_time.out" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --error=${chain_id}_rates_vs_time.err" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --time=200:00:00" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --partition=cobey" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --nodes=1" >> ${chain_id}_rates_vs_time.sbatch
#         echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_rates_vs_time.sbatch
# 
#         echo ""  >> ${chain_id}_rates_vs_time.sbatch
# 
#         echo "cd ../../../"  >> ${chain_id}_rates_vs_time.sbatch
#         echo "python rates_vs_time.py $chain_id"  >> ${chain_id}_rates_vs_time.sbatch
# 
#         mv ${chain_id}_rates_vs_time.sbatch sbatch_files/simulated_alignments/${chain_id}
#     
# 
#    done
# done



