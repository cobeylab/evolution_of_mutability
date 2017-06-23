#!/bin/bash
for scenario in {1,2a,2b,2c,2d,2e,3a,3b,3c,3d}
do
   for replicate in {1..4}
   do
    chain_id=scenario${scenario}_rep${replicate}
    echo "#!/bin/bash" > ${chain_id}_process_output.sbatch
    echo "#SBATCH --job-name=${chain_id}_process_output" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --output=${chain_id}_process_output.out" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --error=${chain_id}_process_output.err" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --time=3:00:00" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --partition=cobey" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --nodes=1" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --ntasks-per-node=1" >> ${chain_id}_process_output.sbatch
    echo "#SBATCH --mem-per-cpu=10000" >> ${chain_id}_process_output.sbatch
    
    echo ""  >> ${chain_id}_process_output.sbatch
    echo "module load java" >> ${chain_id}_process_output.sbatch
    echo "cd /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments" >> ${chain_id}_process_output.sbatch
    echo "./process_output.sh $scenario $replicate" >> ${chain_id}_process_output.sbatch
    
    mv ${chain_id}_process_output.sbatch ${chain_id}/${chain_id}_process_output.sbatch
   
   done

done