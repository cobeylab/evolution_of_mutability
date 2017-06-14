#!/bin/bash

# Sbatch files for observed lineages
for clone in {CH103,CH103int,CH103L,VRC26,VRC26int,VRC26L,VRC26Lint,VRC01_01,VRC01_19,VRC01_13,VRC01_H0306,VRC01_H0306int,VRC01_L0306,VRC01_H08,VRC01_L08}
do
   for prior in {constant,exponential,logistic}
   do
        prior_abbreviation=$(echo $prior | cut -c1-3)
        
        for chain in {a,b,c,d}
        do
            chain_id=${clone}_${prior_abbreviation}_run1${chain}

            echo "#!/bin/bash" > ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --job-name=${chain_id}_S_NS_mutability_changes" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --output=${chain_id}_S_NS_mutability_changes.out" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --error=${chain_id}_S_NS_mutability_changes.err" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --time=600:00:00" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --partition=cobey" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --nodes=1" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "#SBATCH --mem-per-cpu=4000" >> ${chain_id}_S_NS_mutability_changes.sbatch
            
            echo ""  >> ${chain_id}_S_NS_mutability_changes.sbatch

            echo "cd ../../../"  >> ${chain_id}_S_NS_mutability_changes.sbatch
            echo "python S_NS_mutability_changes.py $chain_id"  >> ${chain_id}_S_NS_mutability_changes.sbatch

            mv ${chain_id}_S_NS_mutability_changes.sbatch sbatch_files/observed_lineages/${clone}_${prior}
        
        done
   
   done
done

# Sbatch files for simulated alignments

for scenario in {1,2a,2b,3,4a,4b,5a,5b,5c,6a,6b,6c,7,8a,8b,8c,8d,8e,9,10a,10b,10c,10d,10e,11a,11b,11c,11d,12a,12b,12c,12d}
do
    
    for replicate in {1,2,3,4}
    do
        chain_id=scenario${scenario}_rep${replicate}

        echo "#!/bin/bash" > ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --job-name=${chain_id}_S_NS_mutability_changes" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --output=${chain_id}_S_NS_mutability_changes.out" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --error=${chain_id}_S_NS_mutability_changes.err" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --time=30:00:00" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --partition=cobey" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --nodes=1" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "#SBATCH --mem-per-cpu=1000"  >> ${chain_id}_S_NS_mutability_changes.sbatch

        echo ""  >> ${chain_id}_S_NS_mutability_changes.sbatch

        echo "cd ../../../"  >> ${chain_id}_S_NS_mutability_changes.sbatch
        echo "python S_NS_mutability_changes.py $chain_id"  >> ${chain_id}_S_NS_mutability_changes.sbatch

        mv ${chain_id}_S_NS_mutability_changes.sbatch sbatch_files/simulated_alignments/${chain_id}
    

   done
done



