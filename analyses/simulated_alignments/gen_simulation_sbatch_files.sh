#!/bin/bash

# Sbatch files for simulated alignments

for scenario in {1,2a,2b,3,4a,4b,5a,5b,5c,6a,6b,6c,7,8a,8b,8c,8d,8e,9,10a,10b,10c,10d,10e,11a,11b,11c,12a,12b,12c}
do
    
    for replicate in {1,2,3,4}
    do
        chain_id=scenario${scenario}_rep${replicate}

        echo "#!/bin/bash" > ${chain_id}_simulation.sbatch
        echo "#SBATCH --job-name=${chain_id}_simulation" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --output=${chain_id}_simulation.out" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --error=${chain_id}_simulation.err" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --time=30:00:00" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --partition=cobey" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --nodes=1" >> ${chain_id}_simulation.sbatch
        echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_simulation.sbatch

        echo ""  >> ${chain_id}_simulation.sbatch

        echo "module load python/2.7.12+gcc-4.7" >> ${chain_id}_simulation.sbatch
        
        echo ""  >> ${chain_id}_simulation.sbatch

        echo "cd .."  >> ${chain_id}_simulation.sbatch
        echo "python evolveSequences.py control_files/control_file_$chain_id.py"  >> ${chain_id}_simulation.sbatch

        mv ${chain_id}_simulation.sbatch sbatch_files/
    

   done
done



